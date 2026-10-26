"""Clinical operating points and calibration for FM ABMIL and IEL density on identical tumor tiles.

Thresholds are fixed on TCGA TRAINING patients only (5-fold out-of-fold predictions) and applied unchanged to the
TCGA test split, CPTAC-COAD and SurGen:
  - rule-out: the highest threshold with >= 95% sensitivity out-of-fold (patients below it would skip MSI testing)
  - balanced: Youden's J out-of-fold
Model: Prov-GigaPath-flash gated ABMIL, 5 seeds, up to 1,000 tumor tiles per patient at 0.5 um/px (iel_fm.py setup).
Comparator: IEL density (mean IELs per tile), a fixed score; thresholds from the TRAIN distribution the same way.
Reported per cohort: sensitivity, specificity (Wilson 95% CI), PPV and NPV at the cohort prevalence and at an assumed
15% MSI prevalence, share of patients ruled out, Brier score, observed:expected ratio (mean label / mean predicted
probability) and calibration slope (logistic recalibration). Per-patient predictions -> results/metrics/op_preds.csv.

  python src/operating_points.py   # -> results/metrics/operating_points.csv, results/operating_points.md
"""
import os
import sys

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(__file__))
import iel_fm as I  # noqa: E402
import protocols as P  # noqa: E402

ROOT = P.ROOT
SEEDS = (0, 1, 2, 3, 4)
PREV = 0.15


def score(models, bags):
    out = []
    for k in sorted(bags):
        X = torch.from_numpy(bags[k][1].astype(np.float32)).cuda()
        out.append(np.mean([I.abmil_attention(m, (X - mu) / sd, "cuda")[0] for m, mu, sd in models]))
    return np.array(out)


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return p, c - h, c + h


def thresholds(y, s):
    order = np.sort(np.unique(s))[::-1]
    sens = np.array([(s[y == 1] >= t).mean() for t in order])
    spec = np.array([(s[y == 0] < t).mean() for t in order])
    rule_out = order[np.argmax(sens >= 0.95)]
    youden = order[np.argmax(sens + spec - 1)]
    return {"rule-out (sens >= 0.95 on TRAIN)": rule_out, "Youden (TRAIN)": youden}


def calib(y, p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    lg = np.log(p / (1 - p))
    m = LogisticRegression(C=1e6).fit(lg[:, None], y)
    return m.intercept_[0], m.coef_[0, 0]


def main():
    D = P.Design()
    pid = {p: i for i, p in enumerate(D.pats.patient_id)}
    T = I.merge(I.tcga())
    cl = pd.read_csv(os.path.join(ROOT, "data/cptac/labels.csv"), index_col=0)["label"]
    C = {k: v for k, v in I.merge(I.cptac()).items() if k in cl.index}
    sl = pd.read_csv(os.path.join(ROOT, "data/surgen/subset.csv"), dtype={"case_id": str}).set_index("case_id")
    S = {k: v for k, v in I.merge(I.surgen()).items() if k in sl.index}
    tr = sorted(k for k in T if pid[k] in set(D.train))
    te = sorted(k for k in T if pid[k] in set(D.test))
    ytr = np.array([D.y[pid[k]] for k in tr])
    # out-of-fold ABMIL scores on TRAIN
    oof = np.zeros(len(tr))
    for a, b in StratifiedKFold(5, shuffle=True, random_state=P.SEED).split(tr, ytr):
        ms = [I.train_abmil([T[tr[i]][1] for i in a], ytr[a], s) for s in SEEDS]
        oof[b] = score(ms, {tr[i]: T[tr[i]] for i in b})
    full = [I.train_abmil([T[k][1] for k in tr], ytr, s) for s in SEEDS]
    cohorts = {"tcga_test": ({k: T[k] for k in te}, np.array([D.y[pid[k]] for k in te])),
               "cptac": (C, np.array([int(cl[k]) for k in sorted(C)])),
               "surgen": (S, np.array([int(sl.loc[k, "mmr_loss_binary"]) for k in sorted(S)]))}
    iel = lambda bags: np.array([bags[k][0][:, 6].mean() for k in sorted(bags)])  # noqa: E731
    models = {"FM ABMIL": (oof, lambda b: score(full, b)), "IEL density": (iel({k: T[k] for k in tr}), iel)}
    rows, preds = [], []
    for mname, (s_tr, f) in models.items():
        th = thresholds(ytr, s_tr)
        for cname, (bags, y) in cohorts.items():
            s = f(bags)
            preds += [dict(model=mname, cohort=cname, patient=k, y=int(yy), score=float(v)) for k, yy, v in zip(sorted(bags), y, s)]
            prob = mname == "FM ABMIL"
            b = np.mean((s - y) ** 2) if prob else np.nan
            ci, cs = calib(y, s) if prob else (np.nan, np.nan)
            for tname, t in th.items():
                pos = s >= t
                se = wilson(int((pos & (y == 1)).sum()), int((y == 1).sum()))
                sp = wilson(int((~pos & (y == 0)).sum()), int((y == 0).sum()))
                prev = y.mean()
                ppv = lambda q: se[0] * q / (se[0] * q + (1 - sp[0]) * (1 - q))  # noqa: E731
                npv = lambda q: sp[0] * (1 - q) / (sp[0] * (1 - q) + (1 - se[0]) * q)  # noqa: E731
                rows.append(dict(model=mname, cohort=cname, n=len(y), pos=int(y.sum()), threshold_rule=tname,
                                 threshold=t, sens=se[0], sens_lo=se[1], sens_hi=se[2], spec=sp[0], spec_lo=sp[1],
                                 spec_hi=sp[2], ppv=ppv(prev), npv=npv(prev), ppv_15=ppv(PREV), npv_15=npv(PREV),
                                 ruled_out=(~pos).mean(), missed=int((~pos & (y == 1)).sum()), brier=b,
                                 oe=y.mean() / s.mean() if prob else np.nan,
                                 calib_intercept=ci, calib_slope=cs))
    r = pd.DataFrame(rows)
    os.makedirs(os.path.join(ROOT, "results/metrics"), exist_ok=True)
    r.to_csv(os.path.join(ROOT, "results/metrics/operating_points.csv"), index=False)
    pd.DataFrame(preds).to_csv(os.path.join(ROOT, "results/metrics/op_preds.csv"), index=False)
    with open(os.path.join(ROOT, "results/operating_points.md"), "w") as fh:
        fh.write("# Operating points (thresholds fixed on TCGA TRAIN out-of-fold)\n\n" + r.round(3).to_markdown(index=False) + "\n")
    print(r.round(3).to_string())


if __name__ == "__main__":
    main()
