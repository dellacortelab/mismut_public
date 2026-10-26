"""Patch-level linear probe on the official split, as foundation-model benchmarks use the release (MSI-CRC task).

Each TRAIN tile inherits its patient's label; a logistic regression (standardised features, L2, C=1, lbfgs, fixed a
priori) is fit on all TRAIN tiles and evaluated on all TEST tiles. Arms per encoder:
  released  the Zenodo tiles as distributed (256 um at 224 px = 1.143 um/px, Macenko)        emb/<enc>
  raw114    the same regions re-read from the slides at 1.143 um/px, raw colour              emb/<enc>_rt114
  true05    one 224 px crop of the same region at 0.5 um/px (top-left crop, raw colour)      emb/<enc>_rt05/quads.f16.npy[:, 0]
Metrics on TEST: tile accuracy (threshold 0.5), balanced accuracy, tile ROC-AUC, patient ROC-AUC (mean tile probability).
Paired differences vs `released` with 95% CIs from 2,000 patient-cluster bootstrap resamples.

  python src/patch_probe.py   # -> results/metrics/patch_probe.csv, results/patch_probe.md
"""
import os
import sys

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(__file__))
import pool  # noqa: E402

ROOT = pool.ROOT
ENCS = ["phikon", "phikon_v2", "uni", "uni2_h", "virchow2", "gigapath", "flash"]
LABEL = {"phikon": "Phikon", "phikon_v2": "Phikon-v2", "uni": "UNI", "uni2_h": "UNI2-h", "virchow2": "Virchow2",
         "gigapath": "Prov-GigaPath", "flash": "Prov-GigaPath-flash"}
ARMS = ["released", "raw114", "true05"]
NBOOT = 2000


def load(enc, arm):
    if arm == "released":
        p = os.path.join(ROOT, "emb", enc, "emb.f16.npy")
        return np.load(p, mmap_mode="r") if os.path.exists(p) else None
    if arm == "raw114":
        p = os.path.join(ROOT, "emb", f"{enc}_rt114", "emb.f16.npy")
        return np.load(p, mmap_mode="r") if os.path.exists(p) else None
    p = os.path.join(ROOT, "emb", f"{enc}_rt05", "quads.f16.npy")
    return np.load(p, mmap_mode="r")[:, 0] if os.path.exists(p) else None


def fit_eval(enc, arm, tr, te, y_tr):
    X = load(enc, arm)
    if X is None:
        return None
    Xtr = np.asarray(X[tr], np.float32)
    sc = StandardScaler().fit(Xtr)
    clf = LogisticRegression(C=1.0, max_iter=2000, random_state=42).fit(sc.transform(Xtr), y_tr)
    p = clf.predict_proba(sc.transform(np.asarray(X[te], np.float32)))[:, 1]
    return enc, arm, p


def main():
    idx = pd.read_parquet(os.path.join(ROOT, "data/index.parquet"))
    tr = np.where(idx.split == "TRAIN")[0]
    te = np.where(idx.split == "TEST")[0]
    y = idx.label.values.astype(int)
    pats = idx.patient_id.values
    res = Parallel(n_jobs=6)(delayed(fit_eval)(e, a, tr, te, y[tr]) for e in ENCS for a in ARMS)
    P = {(e, a): p for r in res if r is not None for e, a, p in [r]}
    yt, pt = y[te], pats[te]
    up = np.unique(pt)
    pid = np.searchsorted(up, pt)
    rng = np.random.default_rng(0)
    boots = [rng.choice(len(up), len(up), replace=True) for _ in range(NBOOT)]
    pat_y = pd.Series(yt).groupby(pid).first().values

    def metrics(p):
        correct = (p >= 0.5) == yt
        pat_p = pd.Series(p).groupby(pid).mean().values
        bal = 0.5 * (correct[yt == 1].mean() + correct[yt == 0].mean())
        return dict(acc=correct.mean(), bal_acc=bal, tile_auc=roc_auc_score(yt, p), patient_auc=roc_auc_score(pat_y, pat_p))

    # per-patient sufficient statistics for the cluster bootstrap of accuracy
    def boot_acc(p):
        c = pd.Series(((p >= 0.5) == yt).astype(float)).groupby(pid).agg(["sum", "count"])
        s, n = c["sum"].values, c["count"].values
        return np.array([s[b].sum() / n[b].sum() for b in boots])

    rows = []
    for e in ENCS:
        base = P.get((e, "released"))
        bb = boot_acc(base) if base is not None else None
        for a in ARMS:
            if (e, a) not in P:
                continue
            m = metrics(P[(e, a)])
            row = dict(encoder=e, arm=a, **m)
            if a != "released" and base is not None:
                d = boot_acc(P[(e, a)]) - bb
                row.update(d_acc=m["acc"] - metrics(base)["acc"], d_acc_lo=np.percentile(d, 2.5),
                           d_acc_hi=np.percentile(d, 97.5))
            rows.append(row)
    r = pd.DataFrame(rows)
    os.makedirs(os.path.join(ROOT, "results/metrics"), exist_ok=True)
    r.to_csv(os.path.join(ROOT, "results/metrics/patch_probe.csv"), index=False)
    lines = ["# Patch-level linear probe on the official split (MSI-CRC task)", "",
             "Tile-level logistic regression, TRAIN tiles -> TEST tiles (93,408 / 99,904). Δ accuracy vs released tiles "
             "with 95% patient-cluster bootstrap CI.", "",
             "| Encoder | Arm | Tile acc | Balanced acc | Tile AUC | Patient AUC | Δ acc vs released [95% CI] |",
             "|---|---|---|---|---|---|---|"]
    for x in r.itertuples():
        d = "" if x.arm == "released" or pd.isna(getattr(x, "d_acc", np.nan)) else \
            f"{x.d_acc:+.3f} [{x.d_acc_lo:+.3f}, {x.d_acc_hi:+.3f}]"
        lines.append(f"| {LABEL[x.encoder]} | {x.arm} | {x.acc:.3f} | {x.bal_acc:.3f} | {x.tile_auc:.3f} | "
                     f"{x.patient_auc:.3f} | {d} |")
    # ranking stability
    from scipy.stats import spearmanr
    piv = r.pivot(index="encoder", columns="arm", values="acc").dropna()
    lines += ["", f"Encoder ranking by tile accuracy, Spearman between arms (n={len(piv)} encoders): "
              f"released vs true05 {spearmanr(piv['released'], piv['true05'])[0]:.2f}; "
              f"released vs raw114 {spearmanr(piv['released'], piv['raw114'])[0]:.2f}.",
              "", "Ranking (released): " + " > ".join(LABEL[e] for e in piv["released"].sort_values(ascending=False).index),
              "Ranking (true05):   " + " > ".join(LABEL[e] for e in piv["true05"].sort_values(ascending=False).index)]
    open(os.path.join(ROOT, "results/patch_probe.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
