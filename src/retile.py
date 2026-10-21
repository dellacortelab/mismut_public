"""Re-cut the Kather MSI tiles from the original TCGA slides at 0.5 um/px.

The Zenodo tiles are 256 um regions stored at 224 px (1.143 um/px; see check_resolution.py). For every tile we find
its position in its source slide (barcode from the filename, e.g. TCGA-AA-3858-01Z-00-DX1) and re-read the same
256 um region at 0.5 um/px (512 px, saved as JPEG q95, no colour normalization). Tile selection, per-patient counts
and therefore the patch-count confound stay exactly as in the dataset; only resolution (and stain normalization) differ.

Registration (grayscale, normalized cross-correlation):
  coarse: tile -> 56 px (4.571 um/px), matched over the whole slide at that scale; top K peaks kept
  fine:   full 224 px tile matched in a +-MARGIN window around each peak at 1.143 um/px; best NCC wins
Scale check: for SWEEP_N tiles per slide, the NCC at the found position is recomputed at several candidate pixel
sizes (local search), so the true scale is verified on every slide.

Slides (~239 GB) are kept in MSI_SLIDES (default data/external/tcga_crc_dx, second NVMe) for later
raw-slide work; with --delete-slides they are streamed instead (download, register, cut, delete).

  .venv-tia/bin/python src/retile.py [--workers 16] [--limit-slides N]
Outputs: data/retile/<split dir>/<tile name>.jpg, data/retile/registration.parquet, data/retile/scale_sweep.parquet
"""
import argparse
import json
import os
import re
import subprocess
import time
import traceback
import urllib.request
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import cv2
import numpy as np
import pandas as pd
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data/retile")
SLIDES = os.environ.get("MSI_SLIDES", "data/external/tcga_crc_dx")
TILE_UM, TILE_PX = 256.0, 224
T = TILE_UM / TILE_PX          # 1.143 um/px, tile scale
COARSE = 4 * T                 # 4.571 um/px
TOPK, MARGIN = 3, 32           # coarse candidates; fine search half-window (fine px)
OUT_MPP, OUT_PX = 0.5, 512
SWEEP = (0.9, 1.0, 1.07, 1.143, 1.21, 1.3, 1.4)
SWEEP_N = 3
BARCODE = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}Z-\d{2}-DX\d)")


def gdc_files(patients):
    filt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.submitter_id", "value": patients}},
        {"op": "=", "content": {"field": "data_format", "value": "SVS"}},
        {"op": "=", "content": {"field": "experimental_strategy", "value": "Diagnostic Slide"}}]}
    body = json.dumps({"filters": filt, "fields": "file_id,file_name,file_size", "size": 5000}).encode()
    req = urllib.request.Request("https://api.gdc.cancer.gov/files", data=body,
                                 headers={"Content-Type": "application/json"})
    hits = json.load(urllib.request.urlopen(req))["data"]["hits"]
    return pd.DataFrame([dict(file_id=h["file_id"], file_name=h["file_name"], size=h["file_size"],
                              barcode=h["file_name"].split(".")[0]) for h in hits])


def read_scaled(slide, x0, y0, w0, h0, mpp0, target_mpp, gray=True):
    """Read level-0 box (x0, y0, w0, h0) and return it resampled to target_mpp."""
    ds = target_mpp / mpp0
    lvl = slide.get_best_level_for_downsample(ds + 1e-6)
    lds = slide.level_downsamples[lvl]
    w, h = int(np.ceil(w0 / lds)), int(np.ceil(h0 / lds))
    im = slide.read_region((int(x0), int(y0)), lvl, (w, h)).convert("L" if gray else "RGB")
    out = (int(round(w0 / ds)), int(round(h0 / ds)))
    a = np.asarray(im)
    return cv2.resize(a, out, interpolation=cv2.INTER_AREA)


def peaks(r, k, sup):
    r = r.copy()
    out = []
    for _ in range(k):
        _, mx, _, (x, y) = cv2.minMaxLoc(r)
        out.append((mx, x, y))
        r[max(0, y - sup):y + sup + 1, max(0, x - sup):x + sup + 1] = -1
    return out


def process_slide(barcode, file_id, file_name, tiles, size, delete=False):
    """tiles: list of (row, path). Returns (registration rows, sweep rows, info)."""
    import openslide

    t0 = time.time()
    svs = os.path.join(SLIDES, file_name)
    for attempt in range(5):
        if os.path.exists(svs) and os.path.getsize(svs) == size:
            break
        part = svs + ".part"
        r = subprocess.run(["curl", "-sfL", "--retry", "5", "-o", part, f"https://api.gdc.cancer.gov/data/{file_id}"])
        if r.returncode == 0 and os.path.getsize(part) == size:
            os.replace(part, svs)
            break
        time.sleep(10)
    assert os.path.getsize(svs) == size, f"download failed {file_name}"
    t_dl = time.time() - t0
    try:
        s = openslide.OpenSlide(svs)
        mpp0 = float(s.properties.get("openslide.mpp-x") or s.properties["aperio.MPP"])
        W0, H0 = s.dimensions
        coarse = read_scaled(s, 0, 0, W0, H0, mpp0, COARSE).astype(np.float32)
        ct = int(round(TILE_PX * T / COARSE))  # 56
        rng = np.random.default_rng(zlib.crc32(barcode.encode()))
        sweep_rows = set(rng.choice(len(tiles), min(SWEEP_N, len(tiles)), replace=False).tolist())
        reg, sweep = [], []
        for j, (row, path) in enumerate(tiles):
            g = np.asarray(Image.open(os.path.join(ROOT, path)).convert("L"), dtype=np.float32)
            rc = cv2.matchTemplate(coarse, cv2.resize(g, (ct, ct), interpolation=cv2.INTER_AREA), cv2.TM_CCOEFF_NORMED)
            cands = peaks(rc, TOPK, ct // 2)
            best = (-2, 0, 0)
            for _, cx, cy in cands:
                # candidate top-left in level-0 px, window with margin at fine scale
                fx0 = cx * COARSE / mpp0 - MARGIN * T / mpp0
                fy0 = cy * COARSE / mpp0 - MARGIN * T / mpp0
                side0 = (TILE_PX + 2 * MARGIN) * T / mpp0
                fx0, fy0 = max(0, fx0), max(0, fy0)
                win = read_scaled(s, fx0, fy0, min(side0, W0 - fx0), min(side0, H0 - fy0), mpp0, T).astype(np.float32)
                if win.shape[0] < TILE_PX or win.shape[1] < TILE_PX:
                    continue
                rf = cv2.matchTemplate(win, g, cv2.TM_CCOEFF_NORMED)
                _, mx, _, (x, y) = cv2.minMaxLoc(rf)
                if mx > best[0]:
                    best = (mx, fx0 + x * T / mpp0, fy0 + y * T / mpp0)
            ncc, x0, y0 = best
            side0 = TILE_UM / mpp0
            region = read_scaled(s, x0, y0, min(side0, W0 - x0), min(side0, H0 - y0), mpp0, OUT_MPP, gray=False)
            if region.shape[:2] != (OUT_PX, OUT_PX):
                region = cv2.resize(region, (OUT_PX, OUT_PX), interpolation=cv2.INTER_AREA)
            dst = os.path.join(OUT, os.path.relpath(path, "data/raw"))[:-4] + ".jpg"
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            Image.fromarray(region).save(dst, quality=95)
            reg.append(dict(row=row, barcode=barcode, x0=x0, y0=y0, mpp0=mpp0, ncc=ncc,
                            ncc_coarse1=cands[0][0], ncc_coarse2=cands[1][0] if len(cands) > 1 else np.nan,
                            out=os.path.relpath(dst, ROOT)))
            if j in sweep_rows:
                cx, cy = x0 + side0 / 2, y0 + side0 / 2  # region centre, level 0
                for m in SWEEP:
                    half0 = (TILE_PX / 2 + 16) * m / mpp0
                    if cx - half0 < 0 or cy - half0 < 0 or cx + half0 > W0 or cy + half0 > H0:
                        continue
                    win = read_scaled(s, cx - half0, cy - half0, 2 * half0, 2 * half0, mpp0, m).astype(np.float32)
                    if win.shape[0] < TILE_PX or win.shape[1] < TILE_PX:
                        continue
                    _, mx, _, _ = cv2.minMaxLoc(cv2.matchTemplate(win, g, cv2.TM_CCOEFF_NORMED))
                    sweep.append(dict(row=row, barcode=barcode, mpp=m, ncc=mx))
        info = dict(barcode=barcode, mpp0=mpp0, dims=[W0, H0], n=len(tiles), t_download=t_dl,
                    t_total=time.time() - t0, ncc_median=float(np.median([r["ncc"] for r in reg])))
        return reg, sweep, info
    finally:
        if delete and os.path.exists(svs):
            os.remove(svs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--limit-slides", type=int, default=0)
    ap.add_argument("--delete-slides", action="store_true")
    ap.add_argument("--download-only", action="store_true", help="fetch slides whose tiles are already done")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(SLIDES, exist_ok=True)
    idx = pd.read_parquet(os.path.join(ROOT, "data/index.parquet"))
    idx["barcode"] = idx.path.str.extract(BARCODE)[0]
    assert idx.barcode.notna().all()
    files = gdc_files(sorted(idx.patient_id.unique()))
    missing = set(idx.barcode) - set(files.barcode)
    assert not missing, f"slides not on GDC: {missing}"
    files = files.drop_duplicates("barcode")  # one file per barcode (duplicates are re-uploads of one scan)
    done_dir = os.path.join(OUT, "_done")
    os.makedirs(done_dir, exist_ok=True)
    jobs = []
    for bc, g in idx.reset_index().groupby("barcode"):
        if os.path.exists(os.path.join(done_dir, bc + ".parquet")):
            continue
        f = files[files.barcode == bc].iloc[0]
        jobs.append((bc, f.file_id, f.file_name, list(zip(g["index"], g.path)), f["size"]))
    if a.download_only:
        miss = files[files.barcode.isin(idx.barcode)]
        miss = miss[[not os.path.exists(os.path.join(SLIDES, f)) for f in miss.file_name]]
        print(f"downloading {len(miss)} slides", flush=True)
        for f in miss.itertuples():
            subprocess.run(["curl", "-sfL", "--retry", "5", "-o", os.path.join(SLIDES, f.file_name + ".part"),
                            f"https://api.gdc.cancer.gov/data/{f.file_id}"])
            os.replace(os.path.join(SLIDES, f.file_name + ".part"), os.path.join(SLIDES, f.file_name))
        return
    jobs.sort(key=lambda j: -j[4])  # big slides first
    if a.limit_slides:
        jobs = jobs[:a.limit_slides]
    print(f"{len(jobs)} slides to process ({sum(j[4] for j in jobs) / 1e9:.0f} GB)", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        futs = {ex.submit(process_slide, *j[:5], a.delete_slides): j[0] for j in jobs}
        for k, fu in enumerate(as_completed(futs)):
            bc = futs[fu]
            try:
                reg, sweep, info = fu.result()
            except Exception:
                print("FAILED", bc, traceback.format_exc(), flush=True)
                continue
            pd.DataFrame(reg).to_parquet(os.path.join(done_dir, bc + ".parquet"), index=False)
            if sweep:
                pd.DataFrame(sweep).to_parquet(os.path.join(done_dir, bc + ".sweep.parquet"), index=False)
            sw = pd.DataFrame(sweep)
            peak = sw.loc[sw.groupby("row").ncc.idxmax(), "mpp"].tolist() if len(sw) else []
            print(f"[{k + 1}/{len(jobs)} {time.time() - t0:.0f}s] {bc} mpp0={info['mpp0']:.3f} n={info['n']} "
                  f"ncc_med={info['ncc_median']:.3f} sweep_peak={peak} dl={info['t_download']:.0f}s "
                  f"total={info['t_total']:.0f}s", flush=True)
    regs = [pd.read_parquet(f) for f in sorted(glob_done(done_dir, ".parquet"))]
    if len(regs) == idx.barcode.nunique():
        reg = pd.concat(regs).sort_values("row")
        reg.to_parquet(os.path.join(OUT, "registration.parquet"), index=False)
        sw = [pd.read_parquet(f) for f in glob_done(done_dir, ".sweep.parquet")]
        if sw:
            pd.concat(sw).to_parquet(os.path.join(OUT, "scale_sweep.parquet"), index=False)
        print(f"all {len(reg)} tiles registered; NCC median {reg.ncc.median():.3f}, "
              f"<0.5: {(reg.ncc < 0.5).sum()}", flush=True)


def glob_done(d, suffix):
    fs = [os.path.join(d, f) for f in os.listdir(d) if f.endswith(suffix)]
    return [f for f in fs if suffix == ".sweep.parquet" or not f.endswith(".sweep.parquet")]


if __name__ == "__main__":
    main()
