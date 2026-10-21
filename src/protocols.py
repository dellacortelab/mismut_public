"""Shared evaluation design: folds, protocol runner, metrics, bootstrap.

make_jobs() enumerates every (aggregation, scheme, N, repeat, protocol, cv_repeat, fold, seed) fit with its patient
train/test rows and patch-selection plan. Learners (eval_lr.py, abmil.py) run the same jobs and return out-of-sample
predictions as rows: aggregation, scheme, N, repeat, protocol, cv_repeat, fold, seed, pid (patient row), p.
summarize() turns prediction rows into metrics.
"""
import os

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold

import pool

ROOT = pool.ROOT
SEED = 42
N_BOOT = 5000


class Design:
    def __init__(self):
        self.pats = pool.patients()
        self.y = self.pats.label.values
        self.train = np.where(self.pats.split == "TRAIN")[0]
        self.test = np.where(self.pats.split == "TEST")[0]
        self.n = self.pats.n.values

    def folds(self, subset=None, seed=SEED):
        """Stratified 5-fold over patients (global row indices). subset: patient rows to restrict to."""
        sub = np.arange(len(self.y)) if subset is None else np.asarray(subset)
        skf = StratifiedKFold(5, shuffle=True, random_state=seed)
        return [(sub[a], sub[b]) for a, b in skf.split(sub, self.y[sub])]

    def repeated_folds(self, n_repeats=10):
        rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=n_repeats, random_state=SEED)
        sp = list(rskf.split(np.zeros(len(self.y)), self.y))
        return [sp[r * 5:(r + 1) * 5] for r in range(n_repeats)]

    def paper_folds(self):
        """Fold assignment archived by Vardar & Vardar (split_manifest.csv), for the sanity gate."""
        f = os.path.join(ROOT, "ref/main/Diagnostics_4542710_Zenodo_Revision_v2/repository_files/02_data_processed/"
                               "split_manifest.csv")
        sm = pd.read_csv(f)
        fold = sm[sm.split == "test"].set_index("patient_id").loc[self.pats.patient_id, "fold"].values
        return [(np.where(fold != k)[0], np.where(fold == k)[0]) for k in range(1, 6)]


# ---------------------------------------------------------------- metrics

def metrics(y, p):
    return {"roc_auc": roc_auc_score(y, p), "pr_auc": average_precision_score(y, p), "brier": brier_score_loss(y, p)}


def fast_auc(y, P):
    """AUC for each column of P (n x B) with labels y (n x B), via Mann-Whitney with average ranks."""
    from scipy.stats import rankdata
    r = rankdata(P, axis=0)
    n1 = y.sum(0)
    n0 = y.shape[0] - n1
    return ((r * y).sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)


def boot_indices(y, n_boot=N_BOOT, seed=SEED):
    """Stratified patient bootstrap: resample positives and negatives separately. Returns [n, B] index matrix."""
    rng = np.random.default_rng(seed)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    return np.vstack([rng.choice(pos, (len(pos), n_boot)), rng.choice(neg, (len(neg), n_boot))])


def boot_unit_aucs(df, y_all):
    """Holdout predictions with one or more units (seeds / draws) on the same test patients.
    Returns (test patient rows, bootstrap distribution of the mean-over-units AUC). Resamples depend only on the
    test labels, so they are identical across aggregations and encoders (paired comparisons)."""
    units = [g.set_index("pid").p for _, g in df.groupby(["repeat", "seed"])]
    pids = np.sort(units[0].index.values)
    yy = y_all[pids]
    bi = boot_indices(yy)
    aucs = np.mean([fast_auc(yy[bi], u.reindex(pids).values[bi]) for u in units], axis=0)
    return pids, aucs


def summarize_group(df, y_all):
    """One aggregation x scheme x protocol group of prediction rows -> metric value, spread, CI.

    E1: metrics per (repeat, cv_repeat, seed, fold) -> mean over folds -> mean/SD over (seed or cv_repeat or repeat).
    E2: metrics per (repeat, seed) -> mean/SD over them; 95% CI from the stratified bootstrap of the
        mean-over-units AUC.
    """
    proto = df.protocol.iloc[0]
    res = []
    if proto.startswith("E1"):
        per = []
        for _, g in df.groupby(["repeat", "cv_repeat", "seed"]):
            fm = pd.DataFrame([metrics(y_all[f.pid.values], f.p.values) for _, f in g.groupby("fold")])
            per.append(fm.mean())
        per = pd.DataFrame(per)
        for m in per.columns:
            # SD across units (repeats or seeds); with a single fold set, SD across its folds (paper style)
            sd = per[m].std(ddof=1) if len(per) > 1 else fm[m].std(ddof=1)
            res.append({"metric": m, "value": per[m].mean(), "sd": sd, "n_units": len(per)})
    else:
        per = pd.DataFrame([metrics(y_all[g.pid.values], g.p.values) for _, g in df.groupby(["repeat", "seed"])])
        aucs = boot_unit_aucs(df, y_all)[1]
        for m in per.columns:
            row = {"metric": m, "value": per[m].mean(), "sd": per[m].std(ddof=1) if len(per) > 1 else np.nan,
                   "n_units": len(per)}
            if m == "roc_auc":
                row.update(CI_low=np.percentile(aucs, 2.5), CI_high=np.percentile(aucs, 97.5))
            res.append(row)
    return res


def summarize(preds, y_all, keys=("encoder", "aggregation", "scheme", "N", "protocol")):
    out = []
    for k, g in preds.groupby(list(keys), dropna=False):
        for r in summarize_group(g, y_all):
            out.append({**dict(zip(keys, k)), **r, "n_patients": g.pid.nunique()})
    return pd.DataFrame(out)


# ---------------------------------------------------------------- job list shared by all learners

PLAN_CACHE = {}


def get_plan(D, key):
    """key: ('all',) | ('fixedN', N, repeat) | ('matched', repeat)"""
    if key not in PLAN_CACHE:
        if key[0] == "all":
            PLAN_CACHE[key] = pool.plan_all(D.pats)
        elif key[0] == "fixedN":
            PLAN_CACHE[key] = pool.plan_fixed_n(D.pats, key[1], key[2])
        elif key[0] == "matched":
            PLAN_CACHE[key] = pool.plan_matched(D.pats, key[1])
    return PLAN_CACHE[key]


def make_jobs(D, aggs, seeds, paper_folds=False, schemes=("all", "fixedN", "allN", "matched")):
    """seeds: dict with 'single' (tuple, for E1/E2/allN) and 'per_repeat' (fn repeat -> tuple, for E1rep, fixedN,
    matched). LR uses (None,) everywhere."""
    jobs = []

    def add(agg, scheme, N, rep, proto, cvr, folds, plan, seed_list):
        for seed in seed_list:
            for k, (tr, te) in enumerate(folds):
                jobs.append(dict(aggregation=agg, scheme=scheme, N=N, repeat=rep, protocol=proto, cv_repeat=cvr,
                                 fold=k if proto != "E2" else -1, seed=-1 if seed is None else seed,
                                 tr=tr, te=te, plan=plan))

    holdout = [(D.train, D.test)]
    for agg in aggs:
        if "all" in schemes:
            add(agg, "all", 0, 0, "E1", 0, D.folds(), ("all",), seeds["single"])
            for r, fs in enumerate(D.repeated_folds()):
                add(agg, "all", 0, 0, "E1rep", r, fs, ("all",), seeds["per_repeat"](r))
            add(agg, "all", 0, 0, "E2", 0, holdout, ("all",), seeds["single"])
            if paper_folds:
                add(agg, "all", 0, 0, "E1paper", 0, D.paper_folds(), ("all",), seeds["single"])
        for N in pool.FIXED_N:
            elig = np.where(D.n >= N)[0]
            ef = D.folds(elig)
            eh = [(np.intersect1d(D.train, elig), np.intersect1d(D.test, elig))]
            if "fixedN" in schemes:
                for r in range(pool.N_REPEATS):
                    add(agg, "fixedN", N, r, "E1", 0, ef, ("fixedN", N, r), seeds["per_repeat"](r))
                    add(agg, "fixedN", N, r, "E2", 0, eh, ("fixedN", N, r), seeds["per_repeat"](r))
            if "allN" in schemes:
                add(agg, "allN", N, 0, "E1", 0, ef, ("all",), seeds["single"])
                add(agg, "allN", N, 0, "E2", 0, eh, ("all",), seeds["single"])
        if "matched" in schemes:
            for r in range(pool.N_REPEATS):
                add(agg, "matched", 0, r, "E1", 0, D.folds(), ("matched", r), seeds["per_repeat"](r))
                add(agg, "matched", 0, r, "E2", 0, holdout, ("matched", r), seeds["per_repeat"](r))
    return jobs


def job_rows(job, p):
    meta = {k: job[k] for k in ("aggregation", "scheme", "N", "repeat", "protocol", "cv_repeat", "fold", "seed")}
    return pd.DataFrame({**{k: [v] * len(job["te"]) for k, v in meta.items()}, "pid": job["te"], "p": p})
