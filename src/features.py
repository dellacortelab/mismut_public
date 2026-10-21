"""Interpretable per-tile features: concatenate tissue-type probabilities (tissue.py, 9-d) and nucleus counts
(cells.py, 8-d) into emb/feat_interp/emb.f16.npy, row-aligned with the index. These "encoders" then go through the
same pooling, LR protocols and confound probes as the foundation models. Mean pooling turns them into per-patient
tissue fractions and nucleus densities, which do not depend on patch count by construction.

  python src/features.py
"""
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    parts, names = [], []
    for enc in ("feat_tissue", "feat_cells"):
        d = os.path.join(ROOT, "emb", enc)
        parts.append(np.load(os.path.join(d, "emb.f16.npy")))
        names += [f"{enc[5:]}:{f}" for f in json.load(open(os.path.join(d, "meta.json")))["features"]]
    x = np.hstack(parts)
    assert np.isfinite(x).all()
    out = os.path.join(ROOT, "emb/feat_interp")
    os.makedirs(out, exist_ok=True)
    np.save(os.path.join(out, "emb.f16.npy"), x)
    json.dump({"encoder": "feat_interp", "dim": x.shape[1], "features": names,
               "sources": ["emb/feat_tissue", "emb/feat_cells"]}, open(os.path.join(out, "meta.json"), "w"), indent=1)
    print(x.shape, names)


if __name__ == "__main__":
    main()
