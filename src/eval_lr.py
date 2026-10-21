"""Logistic-regression protocols (E1, E1rep, E2, E3a fixed-N, E3a all-patch comparator, E3b count-matched)
for mean / max / mean+max pooling. Config matches Vardar & Vardar: median impute + StandardScaler + LR(lbfgs, L2,
C=1, class_weight=balanced, max_iter=5000, random_state=42), all fit within training folds.

  python src/eval_lr.py --encoder phikon [--jobs 10]
Writes results/preds/<encoder>_lr.parquet and results/metrics/<encoder>_lr.csv.
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(__file__))
import pool  # noqa: E402
import protocols as P  # noqa: E402

AGGS = ("mean", "max", "mean+max")
# E4: log1p(patch count) alone and appended to each representation (schemes where count varies)
COUNT_AGGS = ("count", "mean+count", "max+count", "mean+max+count")


def lr():
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                         LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000, random_state=P.SEED))


def fit_predict(Xtr, ytr, Xte):
    return lr().fit(Xtr, ytr).predict_proba(Xte)[:, 1]


def run_group(encoder, plan_key, jobs, pooled=None):
    """All jobs sharing one patch-selection plan: pool once, fit every job."""
    D = P.Design()
    plan = P.get_plan(D, plan_key)
    if pooled is None:
        pooled = pool.pool(pool.load_emb(encoder), plan)
    reps = pool.representations(*pooled)
    c = np.log1p(pool.plan_counts(plan)).astype(np.float32)[:, None]
    reps["count"] = c
    for k in AGGS:
        reps[k + "+count"] = np.hstack([reps[k], c])
    out = []
    for j in jobs:
        X = reps[j["aggregation"]]
        out.append(P.job_rows(j, fit_predict(X[j["tr"]], D.y[j["tr"]], X[j["te"]])))
    return pd.concat(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder", required=True)
    ap.add_argument("--jobs", type=int, default=10)
    a = ap.parse_args()
    t0 = time.time()
    D = P.Design()
    seeds = {"single": (None,), "per_repeat": lambda r: (None,)}
    paper = a.encoder in ("phikon", "phikon_ref")
    jobs = P.make_jobs(D, AGGS, seeds, paper_folds=paper)
    jobs += P.make_jobs(D, COUNT_AGGS, seeds, paper_folds=paper, schemes=("all", "allN", "matched"))
    groups = {}
    for j in jobs:
        groups.setdefault((j["plan"], j["aggregation"]), []).append(j)
    all_pooled = pool.pooled_all(a.encoder)
    tasks = [delayed(run_group)(a.encoder, pk, js, all_pooled if pk == ("all",) else None)
             for (pk, _), js in sorted(groups.items(), key=lambda kv: -len(kv[1]))]
    preds = pd.concat(Parallel(n_jobs=a.jobs)(tasks), ignore_index=True)
    preds.insert(0, "encoder", a.encoder)
    os.makedirs(os.path.join(P.ROOT, "results/preds"), exist_ok=True)
    os.makedirs(os.path.join(P.ROOT, "results/metrics"), exist_ok=True)
    preds.to_parquet(os.path.join(P.ROOT, f"results/preds/{a.encoder}_lr.parquet"), index=False)
    m = P.summarize(preds, D.y)
    m.to_csv(os.path.join(P.ROOT, f"results/metrics/{a.encoder}_lr.csv"), index=False)
    show = m[(m.metric == "roc_auc") & (m.scheme == "all")][["aggregation", "protocol", "value", "sd", "CI_low", "CI_high"]]
    print(show.round(3).to_string(index=False))
    print(f"{len(jobs)} fits in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
