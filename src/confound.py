"""E4 confound probes not covered by eval_lr.py (count-only LR and appended-count LR live there).

  python src/confound.py --encoder phikon     # ridge decoding of log1p(count) + score/count correlations
  python src/confound.py --counts             # encoder-free patch-count statistics (direct AUCs, matching check)

Ridge decoding: out-of-fold over the seed-42 stratified 5 folds (all 360 patients) and, to remove label mediation,
within MSS patients only (5-fold KFold, seed 42). In-fold StandardScaler + RidgeCV(alphas 1e-2..1e5, GCV on the
training fold only). Reports R^2, Spearman rho, MAE on log1p(count).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, r2_score, roc_auc_score
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(__file__))
import pool  # noqa: E402
import protocols as P  # noqa: E402

OUT = os.path.join(P.ROOT, "results/metrics")


def decode(X, t, folds):
    oof = np.full(len(t), np.nan)
    for tr, te in folds:
        m = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 5, 15)))
        oof[te] = m.fit(X[tr], t[tr]).predict(X[te])
    s = ~np.isnan(oof)
    return {"r2": r2_score(t[s], oof[s]), "spearman": spearmanr(t[s], oof[s])[0],
            "mae": mean_absolute_error(t[s], oof[s])}


def ridge_decoding(encoder):
    D = P.Design()
    t = np.log1p(D.n.astype(float))
    reps = pool.representations(*pool.pooled_all(encoder))
    mss = np.where(D.y == 0)[0]
    mss_folds = [(mss[a], mss[b]) for a, b in KFold(5, shuffle=True, random_state=P.SEED).split(mss)]
    rows = []
    for agg, X in reps.items():
        for cohort, folds in (("all", D.folds()), ("MSS_only", mss_folds)):
            for metric, v in decode(X, t, folds).items():
                rows.append(dict(encoder=encoder, aggregation=agg, protocol=f"E4_decode_{cohort}",
                                 metric=f"count_{metric}", value=v))
    return rows


def score_count_corr(encoder, learner):
    """Within-class Spearman between out-of-sample score and log1p(count). Uses E1rep (patient-averaged over
    repeats) and E2 predictions. A label-free model should show ~0 within class."""
    f = os.path.join(P.ROOT, f"results/preds/{encoder}_{learner}.parquet")
    if not os.path.exists(f):
        return []
    D = P.Design()
    pr = pd.read_parquet(f)
    pr = pr[(pr.scheme == "all") & pr.protocol.isin(["E1rep", "E2"])]
    rows = []
    for (agg, proto), g in pr.groupby(["aggregation", "protocol"]):
        s = g.groupby("pid").p.mean()
        for lab, name in ((0, "MSS"), (1, "MSIMUT")):
            ids = s.index[D.y[s.index] == lab]
            rows.append(dict(encoder=encoder, aggregation=agg, protocol=f"{proto}_score_count_corr",
                             metric=f"spearman_within_{name}",
                             value=spearmanr(s[ids], np.log1p(D.n[ids]))[0]))
    return rows


def count_stats():
    D = P.Design()
    pats = D.pats
    rows = []
    for split in ("TRAIN", "TEST", "ALL"):
        s = pats if split == "ALL" else pats[pats.split == split]
        rows.append(dict(encoder="count-only", aggregation="count", protocol=f"direct_auc_{split}", metric="roc_auc",
                         value=roc_auc_score(s.label, s.n)))
        for lab in (0, 1):
            c = s.n[s.label == lab]
            rows += [dict(encoder="count-only", aggregation="count", protocol=f"counts_{split}_label{lab}",
                          metric=m, value=v) for m, v in
                     (("n_patients", len(c)), ("mean", c.mean()), ("median", c.median()), ("total", c.sum()))]
    # count-matching check: AUC of counts within TRAIN after matching, per repeat
    tr = np.where(pats.split == "TRAIN")[0]
    aucs, ps = [], []
    for r in range(pool.N_REPEATS):
        n = pool.plan_counts(pool.plan_matched(pats, r))
        aucs.append(roc_auc_score(D.y[tr], n[tr]))
        ps.append(mannwhitneyu(n[tr][D.y[tr] == 1], n[tr][D.y[tr] == 0]).pvalue)
    rows.append(dict(encoder="count-only", aggregation="count", protocol="matched_direct_auc_TRAIN",
                     metric="roc_auc", value=np.mean(aucs), sd=np.std(aucs, ddof=1)))
    rows.append(dict(encoder="count-only", aggregation="count", protocol="matched_mannwhitney_TRAIN",
                     metric="p_min_over_repeats", value=np.min(ps)))
    # fixed-N eligibility
    for N in pool.FIXED_N:
        e = pats[pats.n >= N]
        for split in ("TRAIN", "TEST"):
            s = e[e.split == split]
            rows.append(dict(encoder="count-only", aggregation="count", protocol=f"fixedN_eligible_{split}", metric=f"N{N}",
                             value=len(s), sd=int(s.label.sum())))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder")
    ap.add_argument("--counts", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.counts:
        df = pd.DataFrame(count_stats())
        df.to_csv(os.path.join(OUT, "count_stats.csv"), index=False)
    else:
        df = pd.DataFrame(ridge_decoding(a.encoder) + score_count_corr(a.encoder, "lr")
                          + score_count_corr(a.encoder, "abmil"))
        df.to_csv(os.path.join(OUT, f"{a.encoder}_confound.csv"), index=False)
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
