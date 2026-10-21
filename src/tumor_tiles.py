"""Lossless tile store of whole-tumor tiles for the large-encoder runs (E10b).

Up to CAP tumor-map tiles per patient (seeded random sample), 224 px at 0.5 um/px, for TCGA (data/wsi/regions_flash,
tumor_map) and CPTAC-COAD (data/cptac/tumor_map). Written once as PNG (no extra lossy compression on top of the
scanner's) so every encoder, on any machine, embeds identical pixels.

  python src/tumor_tiles.py build [--workers 24]   # -> $MSI_TILES/{tcga,cptac}/<slide>/<tx>_<ty>.png + index.parquet
  python src/tumor_tiles.py embed --encoder uni2_h [--fp16]   # -> $MSI_TILES/emb/<encoder>.f16.npy (row-aligned)
"""
import argparse
import glob
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
STORE = os.environ.get("MSI_TILES", "data/external/tumor_tiles")
CAP, SEED = 4000, 0


def make_index(cohorts=("tcga", "cptac")):
    rng = np.random.default_rng(SEED)
    dfs = []
    if "tcga" in cohorts:
        t = pd.read_parquet(os.path.join(ROOT, "data/wsi/regions_flash.parquet"))
        dfs.append(t[t.tumor_map].assign(cohort="tcga", patient=lambda d: d.barcode.str[:12], slide=lambda d: d.barcode))
    if "cptac" in cohorts:
        c = pd.read_parquet(os.path.join(ROOT, "data/cptac/tumor_map.parquet"))
        dfs.append(c[c.tumor_map].assign(cohort="cptac", slide=lambda d: d.series))
    if "surgen" in cohorts:
        g = pd.read_parquet(os.path.join(ROOT, "data/surgen/tumor_map.parquet"))
        dfs.append(g[g.tumor_map].assign(cohort="surgen", slide=lambda d: d.series))
    rows = []
    for df in dfs:
        for _, g in df.groupby("patient", sort=True):
            rows.append(g.loc[np.sort(rng.choice(g.index.values, min(CAP, len(g)), replace=False))])
    idx = pd.concat(rows)[["cohort", "patient", "slide", "tx", "ty", "x0", "y0", "mpp0"]].reset_index(drop=True)
    idx["path"] = [f"{co}/{s}/{x}_{y}.png" for co, s, x, y in zip(idx.cohort, idx.slide, idx.tx, idx.ty)]
    return idx


class Reader:
    def __init__(self, cohort, slide, patient):
        if cohort == "surgen":
            import surgen
            p = pd.read_parquet(os.path.join(ROOT, "data/surgen/plan.parquet"), columns=["series", "path"])
            self.r = surgen.Reader(p.path[p.series == slide].iloc[0])
            self.kind = "dcm"
        elif cohort == "tcga":
            import openslide
            self.s = openslide.OpenSlide(glob.glob(f"data/external/tcga_crc_dx/{slide}.*.svs")[0])
            self.kind = "os"
        else:
            import cptac
            self.r = cptac.Reader(os.path.join(cptac.SLIDES, patient, slide))
            self.kind = "dcm"

    def tile(self, x0, y0, mpp0):
        side0 = 224 * 0.5 / mpp0
        if self.kind == "dcm":
            return self.r.read(x0, y0, side0, side0, 0.5, (224, 224))
        lvl = self.s.get_best_level_for_downsample(0.5 / mpp0 + 1e-6)
        n = int(np.ceil(side0 / self.s.level_downsamples[lvl]))
        a = np.asarray(self.s.read_region((int(round(x0)), int(round(y0))), lvl, (n, n)).convert("RGB"))
        return cv2.resize(a, (224, 224), interpolation=cv2.INTER_AREA)


def build_slide(args):
    (cohort, slide, patient), g = args
    rd = Reader(cohort, slide, patient)
    os.makedirs(os.path.join(STORE, cohort, slide), exist_ok=True)
    for r in g.itertuples():
        out = os.path.join(STORE, r.path)
        if not os.path.exists(out):
            Image.fromarray(rd.tile(r.x0, r.y0, r.mpp0)).save(out, compress_level=1)
    return len(g)


def cmd_build(a):
    from concurrent.futures import ProcessPoolExecutor
    idx = make_index(tuple(a.cohorts.split(",")))
    os.makedirs(STORE, exist_ok=True)
    idx.to_parquet(os.path.join(STORE, "index.parquet"), index=False)
    print(f"{len(idx):,} tiles: {idx.cohort.value_counts().to_dict()}", flush=True)
    groups = list(idx.groupby(["cohort", "slide", "patient"], sort=False))
    t0, n = time.time(), 0
    with ProcessPoolExecutor(a.workers) as ex:
        for k in ex.map(build_slide, groups):
            n += k
    print(f"done {n:,} tiles in {time.time() - t0:.0f}s", flush=True)


class Pngs(torch.utils.data.Dataset):
    def __init__(self, paths, tf):
        self.paths, self.tf = paths, tf

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, i):
        with Image.open(os.path.join(STORE, self.paths[i])) as im:
            return self.tf(im.convert("RGB"))


def cmd_embed(a):
    import encoders
    idx = pd.read_parquet(os.path.join(STORE, "index.parquet"))
    model, tf, fwd, dim, info = encoders.load(a.encoder)
    model = model.cuda().eval()
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    os.makedirs(os.path.join(STORE, "emb"), exist_ok=True)
    out = os.path.join(STORE, "emb", f"{a.encoder}.f16.npy")
    emb = np.lib.format.open_memmap(out + ".partial", mode="w+", dtype=np.float16, shape=(len(idx), dim))
    dl = torch.utils.data.DataLoader(Pngs(idx.path.tolist(), tf), batch_size=a.bs, num_workers=a.workers,
                                     prefetch_factor=4)
    t0, o = time.time(), 0
    with torch.inference_mode():
        for x in dl:
            with torch.autocast("cuda", dtype=torch.float16, enabled=a.fp16):
                e = fwd(model, x.cuda(non_blocking=True))
            emb[o:o + len(x)] = e.float().cpu().numpy().astype(np.float16)
            o += len(x)
            if (o // a.bs) % 500 == 0:
                print(f"{o:,}/{len(idx):,} {o / (time.time() - t0):.0f} tiles/s", flush=True)
    emb.flush()
    del emb
    os.replace(out + ".partial", out)
    json.dump({**{k: v for k, v in info.items() if k != "processor"}, "encoder": a.encoder, "n": len(idx),
               "fp16_autocast": a.fp16, "seconds": time.time() - t0},
              open(os.path.join(STORE, "emb", f"{a.encoder}.json"), "w"), indent=1, default=str)
    print(f"done {o:,} in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "embed"])
    ap.add_argument("--encoder", default="flash")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--bs", type=int, default=256)
    ap.add_argument("--fp16", action="store_true")
    ap.add_argument("--cohorts", default="tcga,cptac")
    a = ap.parse_args()
    {"build": cmd_build, "embed": cmd_embed}[a.cmd](a)
