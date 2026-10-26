"""Established biomarker vs foundation model: intraepithelial lymphocytes (IEL) and Prov-GigaPath-flash embeddings on the
SAME tumor tiles (KongNet-counted, up to 1,000 per patient, 0.5 um/px), TCGA TRAIN -> TCGA TEST, CPTAC-COAD, SurGen.

Models (all fit on TCGA TRAIN only):
  iel            IEL density (mean IEL per tile), no fitting
  iel_ratio      IEL per 100 epithelial nuclei, no fitting
  cells_lr       LR on log nucleus densities (8 features + ratio)
  fm_mean        LR on mean-pooled embeddings
  fm_abmil       gated ABMIL (5 seeds)
  late_fusion    LR on [logit(fm_mean out-of-fold score), log IEL density]; FM scores for TRAIN are 5-fold out-of-fold
  early_fusion   LR on [mean embedding, log nucleus densities]
  iel_weighted   LR on IEL-weighted mean embedding (tile weight ∝ 1 + IEL count): pathology-guided pooling
Analyses:
  - Spearman between FM score and IEL density per cohort (how much of the FM signal is IEL)
  - added value: AUC of FM score residualised on IEL (linear, fit on TRAIN)
  - ABMIL attention vs tile IEL count: within-patient Spearman, averaged (does the model look where pathologists look?)
  - site shift: IEL density in MSS patients per cohort

  python src/iel_fm.py   # -> results/iel_fm.md, results/metrics/iel_fm.csv
"""
import os
import sys

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(__file__))
import abmil as A  # noqa: E402
import eval_lr  # noqa: E402
import protocols as P  # noqa: E402
import wsi_eval as W  # noqa: E402
from external import auc_ci  # noqa: E402

ROOT = P.ROOT
CF = ["neutrophil", "epithelial", "lymphocyte", "plasma", "eosinophil", "connective", "iel", "nuclei"]
SEEDS = (0, 1, 2, 3, 4)


def tcga():
    plan = pd.read_parquet(W.PLAN)
    c = pd.read_parquet(os.path.join(ROOT, "data/wsi/cells_tumor.parquet"))
    out = []
    for bc, g in c.groupby("barcode"):
        pl = plan[plan.barcode == bc].reset_index(drop=True)
        e = np.load(os.path.join(W.WSI, "flash", bc + ".npz"))["emb"]
        pos = {(x, y): i for i, (x, y) in enumerate(zip(pl.tx.values, pl.ty.values))}
        rows = [pos[(x, y)] for x, y in zip(g.tx.values, g.ty.values)]
        out.append((bc[:12], g[CF].values, e[rows]))
    return out


def cptac():
    c = pd.read_parquet(os.path.join(ROOT, "data/cptac/cells_tumor.parquet"))
    out = []
    for se, g in c.groupby("series"):
        z = np.load(os.path.join(os.environ.get("MSI_WSI_CPTAC", "/mnt/slides/msimut/wsi_cptac"), f"flash/{se}.npz"))
        pos = {(x, y): i for i, (x, y) in enumerate(zip(z["tx"], z["ty"]))}
        rows = [pos[(x, y)] for x, y in zip(g.tx.values, g.ty.values)]
        out.append((g.patient.iloc[0], g[CF].values, z["emb"][rows]))
    return out


def surgen():
    c = pd.read_parquet(os.path.join(ROOT, "data/surgen/cells_tumor.parquet"))
    out = []
    for case, g in c.groupby("patient"):
        z = np.load(os.path.join(os.environ.get("MSI_WSI_SURGEN", "/mnt/slides/msimut/wsi_surgen"), f"flash/{case}.npz"))
        pos = {(x, y): i for i, (x, y) in enumerate(zip(z["tx"], z["ty"]))}
        rows = [pos[(x, y)] for x, y in zip(g.tx.values, g.ty.values)]
        out.append((case, g[CF].values, z["emb"][rows]))
    return out


def merge(items):
    """Several slides per patient -> one bag per patient."""
    d = {}
    for p, c, e in items:
        if p in d:
            d[p] = (np.vstack([d[p][0], c]), np.vstack([d[p][1], e]))
        else:
            d[p] = (c, e)
    return d


def patient_features(bags):
    keys = sorted(bags)
    cells = np.stack([bags[k][0].mean(0) for k in keys])
    ratio = np.array([100 * bags[k][0][:, 6].sum() / max(1, bags[k][0][:, 1].sum()) for k in keys])
    mean = np.stack([bags[k][1].astype(np.float32).mean(0) for k in keys])
    w = [1 + bags[k][0][:, 6] for k in keys]
    iel_w = np.stack([(bags[k][1].astype(np.float32) * (wk / wk.sum())[:, None]).sum(0) for k, wk in zip(keys, w)])
    cells_log = np.log1p(np.column_stack([cells, ratio]))
    return keys, cells, ratio, mean, iel_w, cells_log


def abmil_attention(model, X, device):
    model.eval()
    with torch.inference_mode():
        h = model.fc(X)
        a = torch.softmax(model.c(model.a(h) * model.b(h)).squeeze(-1), dim=0)
        return torch.sigmoid(model.cls(a @ h)).item(), a.cpu().numpy()


def train_abmil(bags_tr, y_tr, seed, device="cuda:0"):
    """abmil.py architecture and training; returns model and normalisation."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    allx = np.vstack([b for b in bags_tr])
    samp = allx[rng.choice(len(allx), min(50000, len(allx)), replace=False)].astype(np.float32)
    mu, sd = torch.tensor(samp.mean(0), device=device), torch.tensor(samp.std(0) + 1e-6, device=device)
    Xs = [torch.from_numpy(b.astype(np.float32)).to(device) for b in bags_tr]
    model = A.GatedABMIL(allx.shape[1]).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=A.LR, weight_decay=A.WD)
    lossf = torch.nn.BCEWithLogitsLoss()
    w = np.where(y_tr == 1, 1.0 / (y_tr == 1).sum(), 1.0 / (y_tr == 0).sum())
    w /= w.sum()
    model.train()
    for _ in range(A.EPOCHS):
        for k in rng.choice(len(Xs), len(Xs), replace=True, p=w):
            opt.zero_grad()
            lossf(model((Xs[k] - mu) / sd), torch.tensor(float(y_tr[k]), device=device)).backward()
            opt.step()
    return model, mu, sd


def main():
    D = P.Design()
    pid = {p: i for i, p in enumerate(D.pats.patient_id)}
    T = merge(tcga())
    cl = pd.read_csv(os.path.join(ROOT, "data/cptac/labels.csv"), index_col=0)["label"]
    C = {k: v for k, v in merge(cptac()).items() if k in cl.index}
    sl = pd.read_csv(os.path.join(ROOT, "data/surgen/subset.csv"), dtype={"case_id": str}).set_index("case_id")
    S = {k: v for k, v in merge(surgen()).items() if k in sl.index}
    tr_keys = sorted(k for k in T if pid[k] in set(D.train))
    te_keys = sorted(k for k in T if pid[k] in set(D.test))
    cohorts = {"tcga_train": ({k: T[k] for k in tr_keys}, np.array([D.y[pid[k]] for k in tr_keys])),
               "tcga_test": ({k: T[k] for k in te_keys}, np.array([D.y[pid[k]] for k in te_keys])),
               "cptac": (C, np.array([int(cl[k]) for k in sorted(C)])),
               "surgen": (S, np.array([int(sl.loc[k, "mmr_loss_binary"]) for k in sorted(S)]))}
    F = {name: patient_features(b) for name, (b, _) in cohorts.items()}
    Y = {name: y for name, (_, y) in cohorts.items()}
    ytr = Y["tcga_train"]
    _, ctr, rtr, mtr, wtr, cltr = F["tcga_train"]
    # fitted models
    lr_mean = eval_lr.lr().fit(mtr, ytr)
    lr_cells = eval_lr.lr().fit(cltr, ytr)
    lr_early = eval_lr.lr().fit(np.hstack([mtr, cltr]), ytr)
    lr_ielw = eval_lr.lr().fit(wtr, ytr)
    oof = np.zeros(len(ytr))
    for a, b in StratifiedKFold(5, shuffle=True, random_state=42).split(mtr, ytr):
        oof[b] = eval_lr.lr().fit(mtr[a], ytr[a]).predict_proba(mtr[b])[:, 1]
    logit = lambda p: np.log(np.clip(p, 1e-6, 1 - 1e-6) / np.clip(1 - p, 1e-6, 1))  # noqa: E731
    Ztr = np.column_stack([logit(oof), np.log1p(ctr[:, 6])])
    lr_late = LogisticRegression(class_weight="balanced").fit(Ztr, ytr)
    resid_model = LinearRegression().fit(np.log1p(ctr[:, [6]]), logit(oof))
    abm = [train_abmil([cohorts["tcga_train"][0][k][1] for k in tr_keys], ytr, s) for s in SEEDS]
    rows, extra = [], []
    for name in ("tcga_test", "cptac", "surgen"):
        keys, c, r, m, wv, clg = F[name]
        y = Y[name]
        bags = cohorts[name][0]
        s_mean = lr_mean.predict_proba(m)[:, 1]
        att_scores, att_rho = [], []
        for model, mu, sd in abm:
            ps, rhos = [], []
            for k in keys:
                p, a = abmil_attention(model, (torch.from_numpy(bags[k][1].astype(np.float32)).cuda() - mu) / sd, "cuda")
                ps.append(p)
                iel_t = bags[k][0][:, 6]
                if iel_t.std() > 0:
                    rhos.append(spearmanr(a, iel_t)[0])
            att_scores.append(ps)
            att_rho.append(np.nanmean(rhos))
        s_abmil = np.mean(att_scores, 0)
        preds = {"iel": c[:, 6], "iel_ratio": r, "cells_lr": lr_cells.predict_proba(clg)[:, 1], "fm_mean": s_mean,
                 "fm_abmil": s_abmil,
                 "late_fusion": lr_late.predict_proba(np.column_stack([logit(s_mean), np.log1p(c[:, 6])]))[:, 1],
                 "early_fusion": lr_early.predict_proba(np.hstack([m, clg]))[:, 1],
                 "iel_weighted": lr_ielw.predict_proba(wv)[:, 1]}
        for k, v in preds.items():
            a1, l1, h1, _ = auc_ci(y, np.asarray(v, float))
            rows.append(dict(cohort=name, model=k, auc=a1, ci_low=l1, ci_high=h1, n=len(y), pos=int(y.sum())))
        resid = logit(s_mean) - resid_model.predict(np.log1p(c[:, [6]]))
        extra.append(dict(cohort=name, rho_fm_iel=spearmanr(s_mean, c[:, 6])[0],
                          rho_fm_iel_MSS=spearmanr(s_mean[y == 0], c[y == 0, 6])[0],
                          auc_fm_residual_after_iel=roc_auc_score(y, resid),
                          attention_iel_rho=np.mean(att_rho),
                          iel_MSS_mean=c[y == 0, 6].mean(), iel_MSI_mean=c[y == 1, 6].mean(),
                          epi_MSS_mean=c[y == 0, 1].mean()))
    r_ = pd.DataFrame(rows)
    e_ = pd.DataFrame(extra)
    r_.to_csv(os.path.join(ROOT, "results/metrics/iel_fm.csv"), index=False)
    e_.to_csv(os.path.join(ROOT, "results/metrics/iel_fm_analysis.csv"), index=False)
    piv = r_.pivot(index="model", columns="cohort", values="auc")[["tcga_test", "cptac", "surgen"]]
    lines = ["# Established biomarker (IEL) vs foundation model on the same tumor tiles", "",
             "Prov-GigaPath-flash embeddings and KongNet nuclei on identical tiles (≤1,000 tumor tiles per patient, "
             "0.5 µm/px). Fit on TCGA TRAIN only.", "", "| Model | TCGA TEST | CPTAC | SurGen |", "|---|---|---|---|"]
    for mname in ["iel", "iel_ratio", "cells_lr", "fm_mean", "fm_abmil", "late_fusion", "early_fusion", "iel_weighted"]:
        cells_ = []
        for co in ("tcga_test", "cptac", "surgen"):
            x = r_[(r_.model == mname) & (r_.cohort == co)].iloc[0]
            cells_.append(f"{x.auc:.3f} [{x.ci_low:.2f}, {x.ci_high:.2f}]")
        lines.append(f"| {mname} | " + " | ".join(cells_) + " |")
    lines += ["", "## How much of the FM signal is IEL?", "", e_.round(3).to_markdown(index=False)]
    open(os.path.join(ROOT, "results/iel_fm.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
