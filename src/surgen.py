"""External cohort 2: SurGen SR386 (Myles et al. 2025; EBI BioImage Archive S-BIAD1285), Zeiss Axio Scan.Z1 CZI at
0.111 um/px, label = MMR loss by IHC (mmr_loss_binary). Subset downloaded via Globus: all 32 dMMR + 96 random pMMR
cases (data/surgen/subset.csv; one slide per case). Read with pylibCZIrw (pixels are BGR24 -> flipped to RGB).

  python src/surgen.py plan                      # tissue mask at 8 um/px, 224 px @ 0.5 um/px lattice
  python src/surgen.py embed --encoder flash     # all tissue tiles -> $MSI_WSI_SURGEN/<encoder>/<case>.npz
  python src/surgen.py tumormap                  # TCGA-trained tumor classifier (external.tumor_classifier)
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
SLIDES = os.environ.get("MSI_SURGEN", "/mnt/nas/msimut/surgen/SR386_WSIs")
OUT = os.environ.get("MSI_WSI_SURGEN", "/mnt/slides/msimut/wsi_surgen")
DATA = os.path.join(ROOT, "data/surgen")
MASK_UM, TILE_PX, TILE_MPP, MIN_TISSUE, BLOCK = 8.0, 224, 0.5, 0.5, 8
CELLS = int(TILE_PX * TILE_MPP / MASK_UM)


class Reader:
    """Level-0 (bounding-box-relative) coordinates, same interface as cptac.Reader."""

    def __init__(self, path):
        from pylibCZIrw import czi as pyczi
        self.d = pyczi.CziReader(path)
        self.bb = self.d.total_bounding_rectangle
        sc = self.d.metadata["ImageDocument"]["Metadata"]["Scaling"]["Items"]["Distance"]
        self.mpp0 = float(next(s for s in sc if s["@Id"] == "X")["Value"]) * 1e6
        self.W0, self.H0 = self.bb.w, self.bb.h

    def read(self, x0, y0, w0, h0, mpp, out_size):
        zoom = min(1.0, self.mpp0 / mpp)
        a = self.d.read(roi=(int(self.bb.x + x0), int(self.bb.y + y0), int(np.ceil(w0)), int(np.ceil(h0))), zoom=zoom,
                        background_pixel=(255, 255, 255))
        a = np.ascontiguousarray(a[..., ::-1])  # BGR -> RGB
        return cv2.resize(a, out_size, interpolation=cv2.INTER_AREA)


def slides():
    s = pd.read_csv(os.path.join(DATA, "subset.csv"), dtype={"case_id": str})
    s["path"] = [os.path.join(SLIDES, f) for f in s.file]
    return s[[os.path.exists(p) and os.path.getsize(p) == z for p, z in zip(s.path, s["size"])]]


def plan_one(d):
    from wsi_regions import tissue_mask
    rd = Reader(d["path"])
    w, h = int(rd.W0 * rd.mpp0 / MASK_UM), int(rd.H0 * rd.mpp0 / MASK_UM)
    img = rd.read(0, 0, rd.W0, rd.H0, MASK_UM, (w, h))
    tis = tissue_mask(img)
    nty, ntx = h // CELLS, w // CELLS
    tf = tis[:nty * CELLS, :ntx * CELLS].reshape(nty, CELLS, ntx, CELLS).mean(axis=(1, 3))
    ty, tx = np.where(tf >= MIN_TISSUE)
    step0 = TILE_PX * TILE_MPP / rd.mpp0
    return pd.DataFrame(dict(series=d["case_id"], patient=d["case_id"], label=d["mmr_loss_binary"], path=d["path"],
                             tx=tx, ty=ty, x0=tx * step0, y0=ty * step0, mpp0=rd.mpp0, tissue=tf[ty, tx]))


def cmd_plan(a):
    from concurrent.futures import ProcessPoolExecutor
    s = slides()
    with ProcessPoolExecutor(a.workers) as ex:
        plan = pd.concat(list(ex.map(plan_one, s.to_dict("records"))), ignore_index=True)
    plan.to_parquet(os.path.join(DATA, "plan.parquet"), index=False)
    print(f"{len(plan):,} tiles, {plan.series.nunique()} slides ({s.mmr_loss_binary.sum()} dMMR); "
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
                rd, cur = Reader(g.path.iloc[0]), se
            step0 = TILE_PX * TILE_MPP / rd.mpp0
            tx0, ty0 = g.tx.min(), g.ty.min()
            ntx, nty = g.tx.max() - tx0 + 1, g.ty.max() - ty0 + 1
            blk = rd.read(tx0 * step0, ty0 * step0, ntx * step0, nty * step0, TILE_MPP, (ntx * TILE_PX, nty * TILE_PX))
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
    t0, n = time.time(), 0
    with torch.inference_mode():
        for rows, x in dl:
            with torch.autocast("cuda", dtype=torch.float16, enabled=a.fp16):
                e = fwd(model, x.cuda(non_blocking=True))
            emb[rows.numpy()] = e.float().cpu().numpy().astype(np.float16)
            n += len(rows)
            if n % 200000 < len(rows):
                print(f"{n:,}/{len(plan):,} {n / (time.time() - t0):.0f} tiles/s", flush=True)
    for se, g in plan.groupby("series"):
        np.savez(os.path.join(out, f"{se}.npz"), tx=g.tx.values, ty=g.ty.values, emb=emb[g.index.values])
    json.dump({"encoder": a.encoder, "tiles": len(plan), "seconds": time.time() - t0},
              open(os.path.join(out, "meta.json"), "w"), indent=1)
    print(f"embedded {len(plan):,} in {time.time() - t0:.0f}s")


def cmd_tumormap(a):
    import external as X
    t, emb = X.load_tcga()
    clf = X.tumor_classifier(t, emb)
    plan = pd.read_parquet(os.path.join(DATA, "plan.parquet"))
    e = np.zeros((len(plan), emb.shape[1]), np.float16)
    for se, g in plan.groupby("series"):
        e[g.index.values] = np.load(os.path.join(OUT, "flash", f"{se}.npz"))["emb"]
    plan["p_tumor"] = clf.predict_proba(e.astype(np.float32))[:, 1]
    plan["tumor_map"] = X.smooth_map(plan, plan.p_tumor.values)
    plan.to_parquet(os.path.join(DATA, "tumor_map.parquet"), index=False)
    print("tumor fraction:", plan.groupby("series").tumor_map.mean().describe().round(2).to_dict())
    print("tumor tiles per case:", plan[plan.tumor_map].groupby("patient").size().describe().round(0).to_dict())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "embed", "tumormap"])
    ap.add_argument("--encoder", default="flash")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--fp16", action="store_true")
    a = ap.parse_args()
    {"plan": cmd_plan, "embed": cmd_embed, "tumormap": cmd_tumormap}[a.cmd](a)
