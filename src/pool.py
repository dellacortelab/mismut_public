"""Patient table, pooling (mean, max, mean+max), fixed-N and count-matched patch subsampling.

Patients are ordered by patient_id everywhere. Patch rows of a patient are contiguous in index.parquet.
Subsampling plans are deterministic in (scheme, N or repeat) and shared by every encoder and aggregation.

  python src/pool.py --encoder phikon    # caches emb/<encoder>/pooled.npz (all patches)
"""
import argparse
import os

import numpy as np
import pandas as pd

ROOT = os.environ.get("MSI_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIXED_N = (25, 50, 100, 200)
N_REPEATS = 10


def patients():
    idx = pd.read_parquet(os.path.join(ROOT, "data/index.parquet"), columns=["patient_id", "split", "label"])
    g = idx.reset_index().groupby("patient_id", sort=True)
    p = g.agg(split=("split", "first"), label=("label", "first"), start=("index", "min"),
              stop=("index", "max"), n=("index", "size")).reset_index()
    p["stop"] += 1
    assert (p.stop - p.start == p.n).all(), "patient rows not contiguous"
    p["label"] = p.label.astype(int)
    return p


def load_emb(encoder):
    return np.load(os.path.join(ROOT, "emb", encoder, "emb.f16.npy"), mmap_mode="r")


def plan_all(pats):
    """Selection plan: list of row-index arrays per patient (None = all rows)."""
    return [np.arange(s, e) for s, e in zip(pats.start, pats.stop)]


def plan_fixed_n(pats, n, repeat):
    """N patches per eligible patient (n_i >= N), without replacement. Ineligible patients -> None."""
    rng = np.random.default_rng([1, n, repeat])
    out = []
    for s, e in zip(pats.start, pats.stop):
        out.append(np.sort(rng.choice(np.arange(s, e), n, replace=False)) if e - s >= n else None)
    return out


def matched_targets(pats, repeat):
    """Target counts for TRAIN MSIMUT patients so their count distribution matches TRAIN MSS.

    Quantile matching with jitter: MSIMUT patient of count-rank r gets the MSS count quantile at (r + u) / n1,
    u ~ U(0,1) drawn per patient per repeat. Capped at the patient's own count.
    """
    rng = np.random.default_rng([2, repeat])
    tr = pats.split == "TRAIN"
    c0 = np.sort(pats.n[tr & (pats.label == 0)].values)
    pos = np.where(tr & (pats.label == 1))[0]
    order = pos[np.argsort(pats.n.values[pos], kind="stable")]
    q = (np.arange(len(order)) + rng.uniform(size=len(order))) / len(order)
    tgt = np.quantile(c0, q, method="inverted_cdf")
    targets = {}
    for i, t in zip(order, tgt):
        targets[i] = int(min(t, pats.n.values[i]))
    return targets


def plan_matched(pats, repeat):
    targets = matched_targets(pats, repeat)
    rng = np.random.default_rng([3, repeat])
    out = []
    for i, (s, e) in enumerate(zip(pats.start, pats.stop)):
        rows = np.arange(s, e)
        out.append(np.sort(rng.choice(rows, targets[i], replace=False)) if i in targets else rows)
    return out


def plan_counts(plan):
    return np.array([0 if r is None else len(r) for r in plan])


def pool(emb, plan):
    """Returns mean [P,d], max [P,d] (float32); NaN rows for patients with None."""
    d = emb.shape[1]
    mean = np.full((len(plan), d), np.nan, np.float32)
    mx = np.full((len(plan), d), np.nan, np.float32)
    for i, rows in enumerate(plan):
        if rows is None:
            continue
        x = np.asarray(emb[rows[0]:rows[-1] + 1] if _contig(rows) else emb[rows], dtype=np.float32)
        mean[i] = x.mean(0)
        mx[i] = x.max(0)
    return mean, mx


def _contig(rows):
    return len(rows) == rows[-1] - rows[0] + 1


def representations(mean, mx):
    return {"mean": mean, "max": mx, "mean+max": np.hstack([mean, mx])}


def pooled_all(encoder):
    """Cached all-patch pooling."""
    f = os.path.join(ROOT, "emb", encoder, "pooled.npz")
    if os.path.exists(f):
        z = np.load(f)
        return z["mean"], z["max"]
    mean, mx = pool(load_emb(encoder), plan_all(patients()))
    np.savez(f, mean=mean, max=mx)
    return mean, mx


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder", required=True)
    a = ap.parse_args()
    m, x = pooled_all(a.encoder)
    print(a.encoder, m.shape, x.shape)
