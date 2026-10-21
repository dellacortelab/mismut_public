"""E9/E10: whole-tumor nucleus features and external validation on CPTAC-COAD.

All models are fit on the published TCGA TRAIN split only (260 patients) and evaluated, without any tuning, on
TCGA TEST (100) and on CPTAC-COAD (62 patients, 14 MSI-H by PCR). Encoder: Prov-GigaPath-flash.

Inputs per model family
  kather_*   Kather tiles re-cut at 0.5 um/px (emb/flash_rt05: per tile the mean of four 224 px crops of the
             256 um region). Applied to whole-slide tumor tiles as 2x2 tile means (224 um units, the closest
             analogue). LR mean / max / mean+max and ABMIL. The TRAIN split carries the patch-count confound.
  tumor_*    whole-tumor tiles from the tumor map (wsi_eval tumor classifier; for CPTAC the same classifier refit on
             all fully covered TCGA slides). LR mean pooling and ABMIL.
  iel        intraepithelial-lymphocyte density (KongNet, mean IEL per tumor tile), no fitting: score = density.
  cells_lr   LR on per-patient nucleus densities + IEL per 100 epithelial nuclei.
ABMIL: 5 seeds, scores averaged. 95% CIs: 5,000 stratified bootstrap resamples; paired comparisons on identical
resamples.

  python src/external.py tumormap      # CPTAC tumor map -> data/cptac/tumor_map.parquet
  python src/external.py eval          # -> results/external_flash.md, results/metrics/external_flash.csv
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter

sys.path.insert(0, os.path.dirname(__file__))
import eval_lr  # noqa: E402
import pool  # noqa: E402
import protocols as P  # noqa: E402
import wsi_eval as W  # noqa: E402

ROOT = P.ROOT
CPTAC = os.path.join(ROOT, "data/cptac")
WSI_C = os.environ.get("MSI_WSI_CPTAC", "data/external/wsi_cptac")
ENC = "flash"
SEEDS = (0, 1, 2, 3, 4)
MIN_TILES = 20
CELLF = ["neutrophil", "epithelial", "lymphocyte", "plasma", "eosinophil", "connective", "iel", "nuclei"]


def tumor_classifier(t, emb, seed=0):
    """Single tumor-vs-non-tumor LR fit on all fully covered TCGA slides (same sampling as wsi_eval)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    reg = pd.read_csv(os.path.join(ROOT, "results/tile_slide_registration.csv.gz")).groupby("slide").first()
    sparse = set(reg.index[(reg.split == "TRAIN") & (reg.label == 0)])
    rng = np.random.default_rng(seed)
    samp = []
    for bc, g in t[~t.barcode.isin(sparse)].groupby("barcode"):
        pos = g.index[g.region == 0].values
        neg = g.index[(g.region != 0) & (g.dist_um > W.NEG_UM)].values
        samp += [rng.choice(pos, min(1500, len(pos)), replace=False), rng.choice(neg, min(1500, len(neg)), replace=False)]
    samp = np.concatenate(samp)
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=2000, random_state=seed))
    return clf.fit(emb[samp].astype(np.float32), (t.region.values[samp] == 0).astype(int))


def load_tcga():
    plan = pd.read_parquet(W.PLAN)
    reg = pd.read_parquet(os.path.join(ROOT, f"data/wsi/regions_{ENC}.parquet"))
    assert len(reg) == len(plan)
    plan = plan.assign(tumor_map=reg.tumor_map.values, p_tumor=reg.p_tumor.values)
    embs = [np.load(os.path.join(W.WSI, ENC, bc + ".npz"))["emb"] for bc in sorted(plan.barcode.unique())]
    plan = plan.sort_values("barcode", kind="stable")
    assert (plan.index.values == np.arange(len(plan))).all()
    return plan, np.concatenate(embs)


def load_cptac():
    plan = pd.read_parquet(os.path.join(CPTAC, "plan.parquet"))
    embs = []
    for se, g in plan.groupby("series", sort=False):
        z = np.load(os.path.join(WSI_C, ENC, se + ".npz"))
        assert (z["tx"] == g.tx.values).all()
        embs.append((g.index.values, z["emb"]))
    emb = np.zeros((len(plan), embs[0][1].shape[1]), np.float16)
    for rows, e in embs:
        emb[rows] = e
    return plan, emb


def smooth_map(t, p):
    out = np.zeros(len(t), bool)
    for _, g in t.groupby("series"):
        H, Wd = g.ty.max() + 1, g.tx.max() + 1
        grid, have = np.zeros((H, Wd), np.float32), np.zeros((H, Wd), np.float32)
        grid[g.ty, g.tx], have[g.ty, g.tx] = p[g.index], 1
        sm = uniform_filter(grid, 3, mode="constant") / np.maximum(uniform_filter(have, 3, mode="constant"), 1e-6)
        out[g.index] = sm[g.ty, g.tx] >= 0.5
    return out


def cmd_tumormap():
    t, emb = load_tcga()
    clf = tumor_classifier(t, emb)
    c, cemb = load_cptac()
    c["p_tumor"] = clf.predict_proba(cemb.astype(np.float32))[:, 1]
    c["tumor_map"] = smooth_map(c, c.p_tumor.values)
    c.to_parquet(os.path.join(CPTAC, "tumor_map.parquet"), index=False)
    s = pd.read_csv(os.path.join(CPTAC, "slides.csv")).set_index("SeriesInstanceUID")
    frac = c.groupby("series").tumor_map.mean()
    print("tumor fraction of tissue by series type:")
    print(frac.groupby(s.SeriesDescription.reindex(frac.index).values).describe().round(2))
    print("tumor tiles per patient:", c[c.tumor_map].groupby("patient").size().describe().round(0).to_dict())


def blocks2x2(t, rows):
    """Group tumor tiles into 2x2 blocks (all four present) -> list of 4-row arrays."""
    g = t.loc[rows]
    key = (g.ty.values // 2) * 1_000_000 + g.tx.values // 2
    s = pd.Series(g.index.values).groupby(np.asarray(key))
    return [v.values for _, v in s if len(v) == 4]


def auc_ci(y, p):
    from sklearn.metrics import roc_auc_score
    bi = P.boot_indices(y)
    b = P.fast_auc(y[bi], p[bi])
    return roc_auc_score(y, p), np.percentile(b, 2.5), np.percentile(b, 97.5), b


def cmd_eval():
    D = P.Design()
    t, emb = load_tcga()
    t["patient"] = t.barcode.str[:12]
    pid = {p: i for i, p in enumerate(D.pats.patient_id)}
    c = pd.read_parquet(os.path.join(CPTAC, "tumor_map.parquet"))
    _, cemb = load_cptac()
    cl = pd.read_csv(os.path.join(CPTAC, "labels.csv"), index_col=0)["label"]
    # ---- bags: whole-tumor tiles
    tb = {pid[p]: g.index.values for p, g in t[t.tumor_map].groupby("patient") if len(g) >= MIN_TILES}
    cpats = sorted(p for p, g in c[c.tumor_map].groupby("patient") if len(g) >= MIN_TILES)
    cb = {p: c.index.values[(c.patient == p).values & c.tumor_map.values] for p in cpats}
    yc = cl.loc[cpats].values.astype(int)
    tr = np.array(sorted(i for i in tb if i in set(D.train)))
    te = np.array(sorted(i for i in tb if i in set(D.test)))
    print(f"TCGA TRAIN {len(tr)}, TEST {len(te)}; CPTAC {len(cpats)} ({yc.sum()} MSI-H)", flush=True)
    # one embedding table with TCGA rows then CPTAC rows; CPTAC bag keys are 10000 + k
    E = np.concatenate([emb, cemb])
    off = len(emb)
    bags = {**tb, **{10000 + k: cb[p] + off for k, p in enumerate(cpats)}}
    ckeys = np.array([10000 + k for k in range(len(cpats))])
    scores = {}  # model -> (tcga_test_scores, cptac_scores)

    def lr_pool(B, keys, kind):
        x = [E[B[k]].astype(np.float32) for k in keys]
        m = np.stack([a.mean(0) for a in x])
        mx = np.stack([a.max(0) for a in x])
        return {"mean": m, "max": mx, "mean+max": np.hstack([m, mx])}[kind]

    for kind in ("mean",):
        X = lr_pool(bags, np.concatenate([tr, te, ckeys]), kind)
        p = eval_lr.fit_predict(X[:len(tr)], D.y[tr], X[len(tr):])
        scores[f"tumor_lr_{kind}"] = (p[:len(te)], p[len(te):])
    pa = np.mean([W.abmil_fit(E, bags, tr, D.y[tr], np.concatenate([te, ckeys]), s, "cuda:0") for s in SEEDS], 0)
    scores["tumor_abmil"] = (pa[:len(te)], pa[len(te):])
    print("whole-tumor models done", flush=True)
    # ---- Kather-tile models (rt05), trained on Kather TRAIN tiles of all TRAIN patients
    rt = pool.load_emb(f"{ENC}_rt05")
    pats = pool.patients()
    ktr = np.array(sorted(D.train))
    kte = np.array(sorted(D.test))
    kb = {i: np.arange(s, e) for i, (s, e) in enumerate(zip(pats.start, pats.stop))}
    # whole-slide 2x2 units for CPTAC
    cunits, cE = {}, []
    o = 0
    for k, p in zip(ckeys, cpats):
        bl = blocks2x2(c, cb[p])
        if len(bl) < 5:  # fall back to single tiles if too few complete blocks
            v = cemb[cb[p]].astype(np.float32)
        else:
            v = np.stack([cemb[b].astype(np.float32).mean(0) for b in bl])
        cE.append(v)
        cunits[k] = np.arange(o, o + len(v))
        o += len(v)
    RT = np.concatenate([np.asarray(rt, np.float32), np.concatenate(cE)])
    roff = len(rt)
    kbags = {**kb, **{k: v + roff for k, v in cunits.items()}}
    for kind in ("mean", "max", "mean+max"):
        Xs = []
        for k in np.concatenate([ktr, kte, ckeys]):
            a = RT[kbags[k]]
            Xs.append({"mean": a.mean(0), "max": a.max(0), "mean+max": np.concatenate([a.mean(0), a.max(0)])}[kind])
        X = np.stack(Xs)
        p = eval_lr.fit_predict(X[:len(ktr)], D.y[ktr], X[len(ktr):])
        scores[f"kather_lr_{kind}"] = (p[:len(kte)], p[len(kte):])
    RT16 = RT.astype(np.float16)
    pk = np.mean([W.abmil_fit(RT16, kbags, ktr, D.y[ktr], np.concatenate([kte, ckeys]), s, "cuda:0") for s in SEEDS], 0)
    scores["kather_abmil"] = (pk[:len(kte)], pk[len(kte):])
    print("Kather-tile models done", flush=True)
    # ---- nucleus features
    tcf = pd.read_parquet(os.path.join(ROOT, "data/wsi/cells_tumor.parquet"))
    tcf["patient"] = tcf.barcode.str[:12]
    ccf = pd.read_parquet(os.path.join(CPTAC, "cells_tumor.parquet"))

    def cell_table(df):
        g = df.groupby("patient")
        m = g[CELLF].mean()
        m["iel_per_100_epi"] = 100 * g.iel.sum() / g.epithelial.sum().clip(lower=1)
        return m

    tct, cct = cell_table(tcf), cell_table(ccf)
    tct.index = [pid[p] for p in tct.index]
    tcells_te = tct.reindex(te)
    scores["iel_density"] = (tcells_te.iel.values, cct.reindex(cpats).iel.values)
    scores["iel_per_100_epi"] = (tcells_te.iel_per_100_epi.values, cct.reindex(cpats).iel_per_100_epi.values)
    Xtr = np.log1p(tct.reindex(tr).values)
    p = eval_lr.fit_predict(Xtr, D.y[tr], np.log1p(np.vstack([tcells_te.values, cct.reindex(cpats).values])))
    scores["cells_lr"] = (p[:len(te)], p[len(te):])
    # ---- metrics (scores cached so reporting can be rerun without refitting)
    import pickle
    pickle.dump(dict(scores=scores, yc=yc, te=te, kte=kte, cpats=cpats),
                open(os.path.join(ROOT, f"results/preds/external_{ENC}_scores.pkl"), "wb"))
    rows, boots = [], {}
    yte_t = D.y[te]
    for name, (st, sc) in scores.items():
        yt = D.y[kte] if name.startswith("kather") else yte_t
        a1, l1, h1, _ = auc_ci(yt, np.asarray(st, float))
        a2, l2, h2, b2 = auc_ci(yc, np.asarray(sc, float))
        boots[name] = b2
        rows.append(dict(model=name, tcga_test_auc=a1, tcga_ci_low=l1, tcga_ci_high=h1, cptac_auc=a2,
                         cptac_ci_low=l2, cptac_ci_high=h2, n_tcga_test=len(yt), n_cptac=len(yc),
                         n_cptac_pos=int(yc.sum())))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(ROOT, f"results/metrics/external_{ENC}.csv"), index=False)
    pairs = [("tumor_abmil", "kather_lr_max"), ("tumor_abmil", "kather_lr_mean+max"), ("tumor_abmil", "kather_abmil"),
             ("tumor_abmil", "tumor_lr_mean"), ("tumor_abmil", "iel_density"), ("tumor_abmil", "cells_lr"),
             ("kather_lr_mean", "kather_lr_max"), ("iel_density", "kather_lr_max")]
    pr = []
    for a_, b_ in pairs:
        d = boots[a_] - boots[b_]
        pr.append(dict(comparison=f"{a_} - {b_}", delta=res.set_index("model").cptac_auc[a_] -
                       res.set_index("model").cptac_auc[b_], CI_low=np.percentile(d, 2.5),
                       CI_high=np.percentile(d, 97.5)))
    pr = pd.DataFrame(pr)
    pr.to_csv(os.path.join(ROOT, f"results/metrics/external_{ENC}_paired.csv"), index=False)
    lines = ["# External validation: TCGA TRAIN -> TCGA TEST and CPTAC-COAD (Prov-GigaPath-flash)", "",
             f"CPTAC-COAD: {len(yc)} patients, {int(yc.sum())} MSI-H (PCR). Models fit on TCGA TRAIN only.", "",
             "| Model | TCGA TEST AUC [95% CI] | CPTAC AUC [95% CI] |", "|---|---|---|"]
    for r in res.itertuples():
        lines.append(f"| {r.model} | {r.tcga_test_auc:.3f} [{r.tcga_ci_low:.3f}, {r.tcga_ci_high:.3f}] | "
                     f"{r.cptac_auc:.3f} [{r.cptac_ci_low:.3f}, {r.cptac_ci_high:.3f}] |")
    lines += ["", "Paired CPTAC differences:", "", "| Comparison | Δ AUC [95% CI] |", "|---|---|"]
    lines += [f"| {r.comparison} | {r.delta:+.3f} [{r.CI_low:+.3f}, {r.CI_high:+.3f}] |" for r in pr.itertuples()]
    lines += ["", "TCGA whole-tumor nucleus features, TRAIN means (MSIMUT vs MSS):", ""]
    ytr_s = pd.Series(D.y[tr], index=tr)
    tt = tct.reindex(tr)
    for f in CELLF + ["iel_per_100_epi"]:
        lines.append(f"- {f}: {tt[f][ytr_s == 1].mean():.2f} vs {tt[f][ytr_s == 0].mean():.2f}")
    open(os.path.join(ROOT, f"results/external_{ENC}.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["tumormap", "eval"])
    a = ap.parse_args()
    {"tumormap": cmd_tumormap, "eval": cmd_eval}[a.cmd]()
