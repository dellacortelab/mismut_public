"""Place the archived Vardar & Vardar Phikon patch embeddings (ref/, one .npy per patient) into
emb/phikon_ref/emb.f16.npy, row-aligned with data/index.parquet.

Within-patient patch order of the archive is not documented. We assume sorted filename order (same as our index)
and test it against our own Phikon embeddings: per-patch cosine near 1 means the order matches; if it does not,
pooled features (mean/max) are still valid because they are order invariant, but fixed-N draws select different
patches than for our encoders.
"""
import glob
import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    idx = pd.read_parquet(os.path.join(ROOT, "data/index.parquet"))
    files = {os.path.basename(f)[:-4]: f for f in glob.glob(os.path.join(ROOT, "ref/emb_parts/**/patch_level_embeddings_phikon/*.npy"), recursive=True)
             if "__MACOSX" not in f}
    out_dir = os.path.join(ROOT, "emb/phikon_ref")
    os.makedirs(out_dir, exist_ok=True)
    mm = np.lib.format.open_memmap(os.path.join(out_dir, "emb.f16.npy"), mode="w+", dtype=np.float16,
                                   shape=(len(idx), 768))
    mism = []
    g = idx.reset_index().groupby("patient_id")["index"]
    for pid, rows in g:
        a = np.load(files[pid])
        if len(a) != len(rows):
            mism.append((pid, len(a), len(rows)))
            continue
        mm[rows.min():rows.max() + 1] = a.astype(np.float16)
    mm.flush()
    info = {"source": "Zenodo 10.5281/zenodo.22938170, patch_level_embeddings_phikon", "n_patients": len(files),
            "count_mismatches": mism}
    own = os.path.join(ROOT, "emb/phikon/emb.f16.npy")
    if os.path.exists(own) and not mism:
        o = np.load(own, mmap_mode="r")
        rng = np.random.default_rng(0)
        s = np.sort(rng.choice(len(idx), 5000, replace=False))
        x, y = np.asarray(mm[s], np.float32), np.asarray(o[s], np.float32)
        cos = (x * y).sum(1) / np.linalg.norm(x, axis=1) / np.linalg.norm(y, axis=1)
        info.update(patch_cos_vs_own_median=float(np.median(cos)), patch_cos_vs_own_p01=float(np.quantile(cos, 0.01)),
                    patch_cos_vs_own_min=float(cos.min()))
        # order-free check: per-patient mean embedding cosine
        pm = []
        for pid, rows in list(g)[:360]:
            a, b = np.asarray(mm[rows.values], np.float32).mean(0), np.asarray(o[rows.values], np.float32).mean(0)
            pm.append(a @ b / np.linalg.norm(a) / np.linalg.norm(b))
        info.update(patient_mean_cos_min=float(np.min(pm)), patient_mean_cos_median=float(np.median(pm)))
    json.dump(info, open(os.path.join(out_dir, "meta.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))
    if mism:
        raise SystemExit(f"{len(mism)} patients with patch-count mismatch")


if __name__ == "__main__":
    main()
