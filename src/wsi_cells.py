"""E9: nucleus composition of the whole tumor at true 0.5 um/px (KongNet_CoNIC_1, native resolution).

Tiles: up to NMAX tumor-map tiles per patient (data/wsi/regions_flash.parquet, tumor_map == True), sampled with a
fixed seed, read from the original slides at 0.5 um/px (224 px). Densities are per-patient averages, so a capped
random sample estimates them without depending on how much tumor is on the slide. Per tile, cells.tile_features:
counts of neutrophil, epithelial, lymphocyte, plasma, eosinophil, connective nuclei, intraepithelial lymphocytes
(lymphocyte within 10 um of an epithelial nucleus) and all nuclei. Runs in .venv-tia.

  .venv-tia/bin/python src/wsi_cells.py qa          # detection overlay on a few tiles -> data/wsi/cells_qa.png
  .venv-tia/bin/python src/wsi_cells.py run         # -> data/wsi/cells_tumor.parquet
"""
import argparse
import glob
import os
import sys
import time

import cv2
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(__file__))
import cells  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SLIDES = os.environ.get("MSI_SLIDES", "/mnt/slides/msimut/tcga_crc_dx")
REGIONS = os.path.join(ROOT, "data/wsi/regions_flash.parquet")
OUT = os.path.join(ROOT, "data/wsi/cells_tumor.parquet")
NMAX, SEED = 1000, 0


def sample_tiles():
    r = pd.read_parquet(REGIONS)
    r = r[r.tumor_map].copy()
    r["patient"] = r.barcode.str[:12]
    rng = np.random.default_rng(SEED)
    parts = []
    for p, g in r.groupby("patient", sort=True):
        parts.append(g.loc[np.sort(rng.choice(g.index.values, min(NMAX, len(g)), replace=False))])
    return pd.concat(parts).sort_values(["barcode", "ty", "tx"]).reset_index(drop=True)


class Tiles(torch.utils.data.Dataset):
    def __init__(self, t):
        self.t = t
        self.cache = {}

    def __len__(self):
        return len(self.t)

    def __getitem__(self, i):
        import openslide
        r = self.t.iloc[i]
        if r.barcode not in self.cache:
            self.cache = {r.barcode: openslide.OpenSlide(glob.glob(os.path.join(SLIDES, r.barcode + ".*.svs"))[0])}
        s = self.cache[r.barcode]
        side0 = 224 * 0.5 / r.mpp0
        lvl = s.get_best_level_for_downsample(0.5 / r.mpp0 + 1e-6)
        lds = s.level_downsamples[lvl]
        n = int(np.ceil(side0 / lds))
        a = np.asarray(s.read_region((int(round(r.x0)), int(round(r.y0))), lvl, (n, n)).convert("RGB"))
        return torch.from_numpy(cv2.resize(a, (224, 224), interpolation=cv2.INTER_AREA))


def load_model():
    from tiatoolbox.models.architecture import get_pretrained_model
    m, _ = get_pretrained_model("KongNet_CoNIC_1")
    return m.cuda().eval()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["qa", "run"])
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--bs", type=int, default=64)
    a = ap.parse_args()
    t = sample_tiles()
    model = load_model()
    dev = torch.device("cuda")
    if a.cmd == "qa":
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sub = t.sample(8, random_state=1).reset_index(drop=True)
        ds = Tiles(sub)
        x = torch.stack([ds[i] for i in range(len(sub))])
        with torch.inference_mode():
            det = cells.detect(model, x, dev, True, 1)
        col = ["yellow", "lime", "blue", "magenta", "orange", "cyan"]
        fig, axs = plt.subplots(2, 8, figsize=(32, 8.5))
        for i in range(8):
            for r in (0, 1):
                axs[r, i].imshow(x[i].numpy())
                axs[r, i].axis("off")
            yy, xx, cc = det[i]
            for c in range(6):
                axs[1, i].scatter(xx[cc == c], yy[cc == c], s=14, c=col[c], label=cells.CLASSES[c])
            axs[0, i].set_title(f"n={len(cc)} " + " ".join(f"{cells.CLASSES[c][:3]}{(cc == c).sum()}" for c in range(6)),
                                fontsize=8)
        axs[1, 0].legend(fontsize=7, loc="lower left")
        plt.tight_layout()
        plt.savefig(os.path.join(ROOT, "data/wsi/cells_qa.png"), dpi=60)
        return
    print(f"{len(t):,} tumor tiles from {t.barcode.str[:12].nunique()} patients", flush=True)
    dl = torch.utils.data.DataLoader(Tiles(t), batch_size=a.bs, num_workers=a.workers, shuffle=False,
                                     prefetch_factor=4)
    feats, t0 = [], time.time()
    with torch.inference_mode():
        for bi, x in enumerate(dl):
            for y_, x_, c_ in cells.detect(model, x, dev, True, 1):
                feats.append(cells.tile_features(y_, x_, c_))
            if (bi + 1) % 200 == 0:
                n = (bi + 1) * a.bs
                print(f"{n:,}/{len(t):,} {n / (time.time() - t0):.0f} tiles/s", flush=True)
    f = pd.DataFrame(np.stack(feats), columns=cells.FEATURES)
    pd.concat([t[["barcode", "tx", "ty", "p_tumor"]], f], axis=1).to_parquet(OUT, index=False)
    print("done", len(f), f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
