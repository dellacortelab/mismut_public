"""Measure the true pixel size of the Kather MSI tiles by registering them to their original TCGA slide.

Zenodo 2530835 documents the tiles as 224 px at 0.5 um/px. Kather et al. 2019 (Nat Med, Methods) state 256 um edge
length. This script downloads one diagnostic slide from GDC, rescales tiles over a range of candidate pixel sizes and
template-matches them (grayscale, normalized cross-correlation) against the slide at ~1 um/px. Runs in .venv-tia
(openslide + opencv come with tiatoolbox).

  .venv-tia/bin/python src/check_resolution.py [--patient TCGA-AD-6895] [--n 12]

Result (2026-10-16, TCGA-AD-6895, 12 tiles): every tile peaks at 1.15 um/px (NCC 0.80-0.97; mean NCC 0.88 at 1.15
vs 0.28 at 0.50), and match positions fall on a 256 um grid. Tiles are 256 um / 224 px = 1.143 um/px.
"""
import argparse
import json
import os
import urllib.request

import cv2
import numpy as np
import openslide
import pandas as pd
from PIL import Image

ROOT = os.environ.get("MSI_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def gdc_dx_file(patient):
    filt = {"op": "and", "content": [
        {"op": "=", "content": {"field": "cases.submitter_id", "value": patient}},
        {"op": "=", "content": {"field": "data_format", "value": "SVS"}},
        {"op": "=", "content": {"field": "experimental_strategy", "value": "Diagnostic Slide"}}]}
    body = json.dumps({"filters": filt, "fields": "file_id,file_name", "size": 20}).encode()
    req = urllib.request.Request("https://api.gdc.cancer.gov/files", data=body,
                                 headers={"Content-Type": "application/json"})
    hits = json.load(urllib.request.urlopen(req))["data"]["hits"]
    return next(h for h in hits if "DX1" in h["file_name"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", default="TCGA-AD-6895")
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--out", default="results/metrics/tile_resolution_check.json")
    a = ap.parse_args()
    out = os.environ.get("MSI_SLIDES", "data/external/tcga_crc_dx")
    os.makedirs(out, exist_ok=True)
    f = gdc_dx_file(a.patient)
    svs = os.path.join(out, f["file_name"])
    if not os.path.exists(svs):
        urllib.request.urlretrieve(f"https://api.gdc.cancer.gov/data/{f['file_id']}", svs)
    s = openslide.OpenSlide(svs)
    level = s.get_best_level_for_downsample(4)
    mpp = float(s.properties["openslide.mpp-x"]) * s.level_downsamples[level]
    big = np.asarray(s.read_region((0, 0), level, s.level_dimensions[level]).convert("L"), dtype=np.float32)

    idx = pd.read_parquet(os.path.join(ROOT, "data/index.parquet"))
    tiles = idx[idx.patient_id == a.patient].path.tolist()
    pick = np.random.default_rng(0).choice(len(tiles), min(a.n, len(tiles)), replace=False)
    scales = np.round(np.arange(0.45, 1.46, 0.05), 2)
    ncc, locs = [], []
    for k in pick:
        g = np.asarray(Image.open(os.path.join(ROOT, tiles[k])).convert("L"), dtype=np.float32)
        row, rl = [], []
        for m in scales:
            sz = int(round(224 * m / mpp))
            r = cv2.matchTemplate(big, cv2.resize(g, (sz, sz), interpolation=cv2.INTER_AREA), cv2.TM_CCOEFF_NORMED)
            _, mx, _, loc = cv2.minMaxLoc(r)
            row.append(mx)
            rl.append(loc)
        ncc.append(row)
        locs.append(rl[int(np.argmax(row))])
    ncc = np.array(ncc)
    best = scales[ncc.argmax(1)]
    res = {"slide": f["file_name"], "slide_level_mpp": mpp, "n_tiles": len(pick),
           "best_mpp_per_tile": best.tolist(), "best_ncc_per_tile": ncc.max(1).round(3).tolist(),
           "mean_ncc_by_mpp": dict(zip(scales.tolist(), ncc.mean(0).round(3).tolist())),
           "match_xy_um": [[round(x * mpp), round(y * mpp)] for x, y in locs]}
    print(json.dumps(res, indent=1))
    out = os.path.join(ROOT, a.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
