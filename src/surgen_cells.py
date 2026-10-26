"""KongNet nucleus features (incl. intraepithelial lymphocytes) on SurGen tumor tiles, same protocol as wsi_cells.py:
up to 1,000 tumor tiles per case from the lossless tile store (0.5 um/px), seed 0. Runs in .venv-tia.

  .venv-tia/bin/python src/surgen_cells.py   # -> data/surgen/cells_tumor.parquet
"""
import os
import sys

import numpy as np
import pandas as pd
import torch
from PIL import Image

sys.path.insert(0, os.path.dirname(__file__))
import cells  # noqa: E402
from wsi_cells import load_model  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STORE = os.environ.get("MSI_TILES_SURGEN", "/mnt/slides/msimut/tumor_tiles_surgen")
NMAX = 1000


class Pngs(torch.utils.data.Dataset):
    def __init__(self, paths):
        self.paths = paths

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, i):
        return torch.from_numpy(np.asarray(Image.open(os.path.join(STORE, self.paths[i])).convert("RGB")))


def main():
    idx = pd.read_parquet(os.path.join(STORE, "index.parquet"))
    rng = np.random.default_rng(0)
    t = pd.concat([g.loc[np.sort(rng.choice(g.index.values, min(NMAX, len(g)), replace=False))]
                   for _, g in idx.groupby("patient", sort=True)]).reset_index(drop=True)
    model = load_model()
    dl = torch.utils.data.DataLoader(Pngs(t.path.tolist()), batch_size=64, num_workers=12, prefetch_factor=4)
    feats = []
    with torch.inference_mode():
        for x in dl:
            feats += [cells.tile_features(*d) for d in cells.detect(model, x, torch.device("cuda"), True, 1)]
    f = pd.DataFrame(np.stack(feats), columns=cells.FEATURES)
    pd.concat([t[["patient", "tx", "ty"]], f], axis=1).to_parquet(os.path.join(ROOT, "data/surgen/cells_tumor.parquet"),
                                                                  index=False)
    print("cells", len(f))


if __name__ == "__main__":
    main()
