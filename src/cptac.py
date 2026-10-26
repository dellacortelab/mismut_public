"""External cohort: CPTAC-COAD (IDC, DICOM-WSM, 0.25 um/px), tiled and embedded exactly like the TCGA slides.

Slides: H&E "tumor_tissue" or unlabeled "HE" series of patients with a PCR MSI call (LinkedOmics clinical file,
Mutation_Phenotype; Vasaikar et al. 2019): 110 slides, 62 patients (14 MSI-H, 48 MSS); data/cptac/slides.csv.
Read with wsidicom (OpenSlide cannot decode their JPEG 2000 YBR_FULL_422 tiles).

  python src/cptac.py download                 # idc-index -> /mnt/slides/msimut/cptac_coad/<patient>/<series>/
  python src/cptac.py plan                     # tissue mask at 8 um/px, 224 px @ 0.5 um/px lattice -> data/cptac/plan.parquet
  python src/cptac.py embed --encoder flash    # -> /mnt/slides/msimut/wsi_cptac/<encoder>/<series>.npz
  .venv-tia/bin/python src/cptac.py cells      # KongNet on up to NMAX tumor-map tiles/patient -> data/cptac/cells_tumor.parquet
                                               # (needs data/cptac/tumor_map.parquet from external.py)
"""
import argparse
import json
import os
import sys
import time

import cv2
import numpy as np
import pandas as pd
import torch
from PIL import Image

sys.path.insert(0, os.path.dirname(__file__))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SLIDES = os.environ.get("MSI_CPTAC", "/mnt/slides/msimut/cptac_coad")
OUT = os.environ.get("MSI_WSI_CPTAC", "/mnt/slides/msimut/wsi_cptac")
DATA = os.path.join(ROOT, "data/cptac")
MASK_UM, TILE_PX, TILE_MPP, MIN_TISSUE, BLOCK, NMAX = 8.0, 224, 0.5, 0.5, 8, 1000
CELLS = int(TILE_PX * TILE_MPP / MASK_UM)


def series_dir(r):
    return os.path.join(SLIDES, r["PatientID"], r["SeriesInstanceUID"])


class Reader:
    """Minimal level-0-coordinate reader over wsidicom."""

    def __init__(self, path):
        from wsidicom import WsiDicom
        self.w = WsiDicom.open(path)
        self.mpp0 = float(self.w.mpp.width)
        self.W0, self.H0 = self.w.size.width, self.w.size.height
        self.levels = [(l.level, float(l.mpp.width)) for l in self.w.levels]

    def read(self, x0, y0, w0, h0, mpp, out_size):
        """Level-0 box -> RGB array resized to out_size, read from the coarsest level with mpp <= target."""
        lv, lm = max(((l, m) for l, m in self.levels if m <= mpp + 1e-6), key=lambda t: t[1])
        f = self.mpp0 / lm
        x, y = int(round(x0 * f)), int(round(y0 * f))
        w, h = max(1, int(np.ceil(w0 * f))), max(1, int(np.ceil(h0 * f)))
        L = next(l for l in self.w.levels if l.level == lv).size
        cw, ch = max(1, min(w, L.width - x)), max(1, min(h, L.height - y))
        a = np.asarray(self.w.read_region((x, y), lv, (cw, ch)).convert("RGB"))
        if (cw, ch) != (w, h):  # clipped at the slide edge: pad with background
            full = np.full((h, w, 3), 255, np.uint8)
            full[:ch, :cw] = a
            a = full
        return cv2.resize(a, out_size, interpolation=cv2.INTER_AREA)


def cmd_download(a):
    from idc_index import index
    s = pd.read_csv(os.path.join(DATA, "slides.csv"))
    index.IDCClient().download_from_selection(seriesInstanceUID=s.SeriesInstanceUID.tolist(), downloadDir=SLIDES,
                                              dirTemplate="%PatientID/%SeriesInstanceUID", show_progress_bar=False)


def plan_one(d):
    r = pd.Series(d)
    from wsi_regions import tissue_mask
    rd = Reader(series_dir(r))
    w, h = int(rd.W0 * rd.mpp0 / MASK_UM), int(rd.H0 * rd.mpp0 / MASK_UM)
    img = rd.read(0, 0, rd.W0, rd.H0, MASK_UM, (w, h))
    tis = tissue_mask(img)
    nty, ntx = h // CELLS, w // CELLS
    tf = tis[:nty * CELLS, :ntx * CELLS].reshape(nty, CELLS, ntx, CELLS).mean(axis=(1, 3))
    ty, tx = np.where(tf >= MIN_TISSUE)
    step0 = TILE_PX * TILE_MPP / rd.mpp0
    return pd.DataFrame(dict(series=r.SeriesInstanceUID, patient=r.PatientID, label=r.label, tx=tx, ty=ty,
                             x0=tx * step0, y0=ty * step0, mpp0=rd.mpp0, tissue=tf[ty, tx]))


def cmd_plan(a):
    from concurrent.futures import ProcessPoolExecutor
    s = pd.read_csv(os.path.join(DATA, "slides.csv"))
    with ProcessPoolExecutor(a.workers) as ex:
        plan = pd.concat(list(ex.map(plan_one, s.to_dict("records"))), ignore_index=True)
    plan.to_parquet(os.path.join(DATA, "plan.parquet"), index=False)
    print(f"{len(plan):,} tiles, {plan.series.nunique()} slides, {plan.patient.nunique()} patients; "
          f"median per slide {plan.groupby('series').size().median():.0f}")


class Blocks(torch.utils.data.IterableDataset):
    def __init__(self, plan, tf):
        self.plan, self.tf = plan, tf
        self.items = []
        for se, g in plan.groupby("series", sort=True):
            for _, gg in g.groupby([g.ty.values // BLOCK, g.tx.values // BLOCK]):
                self.items.append((se, gg.index.values))

    def __iter__(self):
        wi = torch.utils.data.get_worker_info()
        items = self.items if wi is None else self.items[wi.id::wi.num_workers]
        cur, rd = None, None
        for se, rows in items:
            g = self.plan.loc[rows]
            if se != cur:
                rd, cur = Reader(os.path.join(SLIDES, g.patient.iloc[0], se)), se
            step0 = TILE_PX * TILE_MPP / rd.mpp0
            tx0, ty0 = g.tx.min(), g.ty.min()
            ntx, nty = g.tx.max() - tx0 + 1, g.ty.max() - ty0 + 1
            blk = rd.read(tx0 * step0, ty0 * step0, ntx * step0, nty * step0, TILE_MPP,
                          (ntx * TILE_PX, nty * TILE_PX))
            tiles = [self.tf(Image.fromarray(blk[(ty - ty0) * TILE_PX:(ty - ty0 + 1) * TILE_PX,
                                                 (tx - tx0) * TILE_PX:(tx - tx0 + 1) * TILE_PX]))
                     for tx, ty in zip(g.tx.values, g.ty.values)]
            yield torch.from_numpy(rows.copy()), torch.stack(tiles)


def cmd_embed(a):
    import encoders
    plan = pd.read_parquet(os.path.join(DATA, "plan.parquet"))
    out = os.path.join(OUT, a.encoder)
    os.makedirs(out, exist_ok=True)
    model, tf, fwd, dim, info = encoders.load(a.encoder)
    model = model.cuda().eval()
    emb = np.zeros((len(plan), dim), np.float16)
    dl = torch.utils.data.DataLoader(Blocks(plan, tf), batch_size=None, num_workers=a.workers, prefetch_factor=4)
    t0 = time.time()
    with torch.inference_mode():
        for rows, x in dl:
            x = x.cuda(non_blocking=True)
            with torch.autocast("cuda", dtype=torch.float16, enabled=a.fp16):
                e = fwd(model, x)
            emb[rows.numpy()] = e.float().cpu().numpy().astype(np.float16)
    for se, g in plan.groupby("series"):
        np.savez(os.path.join(out, se + ".npz"), tx=g.tx.values, ty=g.ty.values, emb=emb[g.index.values])
    json.dump({"encoder": a.encoder, "tiles": len(plan), "seconds": time.time() - t0, "fp16_autocast": a.fp16},
              open(os.path.join(out, "meta.json"), "w"), indent=1)
    print(f"embedded {len(plan):,} tiles in {time.time() - t0:.0f}s")


class CellTiles(torch.utils.data.Dataset):
    def __init__(self, t):
        self.t, self.rd, self.cur = t, None, None

    def __len__(self):
        return len(self.t)

    def __getitem__(self, i):
        r = self.t.iloc[i]
        if r.series != self.cur:
            self.rd, self.cur = Reader(os.path.join(SLIDES, r.patient, r.series)), r.series
        side0 = TILE_PX * TILE_MPP / self.rd.mpp0
        return torch.from_numpy(self.rd.read(r.x0, r.y0, side0, side0, TILE_MPP, (TILE_PX, TILE_PX)))


def cmd_cells(a):
    import cells
    from wsi_cells import load_model
    tm = pd.read_parquet(os.path.join(DATA, "tumor_map.parquet"))
    tm = tm[tm.tumor_map]
    rng = np.random.default_rng(0)
    t = pd.concat([g.loc[np.sort(rng.choice(g.index.values, min(NMAX, len(g)), replace=False))]
                   for _, g in tm.groupby("patient", sort=True)]).sort_values(["series", "ty", "tx"]).reset_index(drop=True)
    model = load_model()
    dl = torch.utils.data.DataLoader(CellTiles(t), batch_size=64, num_workers=a.workers, prefetch_factor=4)
    feats = []
    with torch.inference_mode():
        for x in dl:
            feats += [cells.tile_features(*d) for d in cells.detect(model, x, torch.device("cuda"), True, 1)]
    f = pd.DataFrame(np.stack(feats), columns=cells.FEATURES)
    pd.concat([t[["series", "patient", "tx", "ty"]], f], axis=1).to_parquet(os.path.join(DATA, "cells_tumor.parquet"),
                                                                           index=False)
    print("cells", len(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["download", "plan", "embed", "cells"])
    ap.add_argument("--encoder", default="flash")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--fp16", action="store_true")
    a = ap.parse_args()
    {"download": cmd_download, "plan": cmd_plan, "embed": cmd_embed, "cells": cmd_cells}[a.cmd](a)


if __name__ == "__main__":
    main()
