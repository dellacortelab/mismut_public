"""Crops for the overview figure (Fig. 1): one 224 um field (448 px at 0.5 um/px) from the top high-attention review
regions of the MSI review slides, with KongNet nuclei. Picks the field with the most intraepithelial lymphocytes.
Runs in .venv-tia:  .venv-tia/bin/python src/overview_crops.py   -> data/review/overview_crop.npz
"""
import glob
import json
import os
import sys

import numpy as np
import openslide
import pandas as pd
import torch
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(__file__))
import cells as C  # noqa: E402

ROOT = C.ROOT
OUT = os.path.join(ROOT, "data/review")
SLIDES = os.environ.get("MSI_SLIDES", "/mnt/slides/msimut") + "/tcga_crc_dx"


def main():
    from tiatoolbox.models.architecture import get_pretrained_model
    model, _ = get_pretrained_model("KongNet_CoNIC_1")
    model = model.cuda().eval()
    key = pd.read_csv(os.path.join(OUT, "regions_key.csv"))
    hk = key[(key.kind == "high") & (key.group == "MSI_correct")].sort_values("rel_attention", ascending=False)
    hk = hk.groupby("barcode").head(2)
    best = None
    for r in hk.itertuples():
        fc = {f["properties"]["region_id"]: f for f in json.load(open(os.path.join(OUT, "regions", r.barcode + ".geojson")))["features"]}
        cx, cy = np.array(fc[r.region_id]["geometry"]["coordinates"][0]).mean(0)
        sl = openslide.OpenSlide(glob.glob(os.path.join(SLIDES, r.barcode + ".*.svs"))[0])
        mpp0 = float(sl.properties["openslide.mpp-x"])
        size0 = int(round(448 * 0.5 / mpp0))
        im = sl.read_region((int(cx - size0 / 2), int(cy - size0 / 2)), 0, (size0, size0)).convert("RGB").resize((448, 448))
        a = np.asarray(im)
        y, x, c = C.detect(model, torch.from_numpy(a[None].copy()), torch.device("cuda"), False, 1)[0]
        f = C.tile_features(y, x, c)
        print(r.barcode, r.region_id, "IEL", f[6], "lym", f[2], "neu", f[0], "epi", f[1], flush=True)
        if best is None or f[6] > best[0]:
            best = (f[6], a, y, x, c, r.barcode, r.region_id)
    n, a, y, x, c, bc, rid = best
    lym, epi = c == 2, c == 1
    iel = np.zeros(len(c), bool)
    if lym.any() and epi.any():
        d, _ = cKDTree(np.c_[y[epi], x[epi]]).query(np.c_[y[lym], x[lym]], k=1)
        iel[np.where(lym)[0][d * C.MODEL_MPP <= C.IEL_UM]] = True
    np.savez_compressed(os.path.join(OUT, "overview_crop.npz"), img=a, y=y, x=x, c=c, iel=iel, barcode=bc, region=rid)
    print("picked", bc, rid, "IEL", n)


if __name__ == "__main__":
    main()
