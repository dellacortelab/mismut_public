"""E10b: whole-tumor models with large encoders on the shared tile store, TCGA TRAIN -> TCGA TEST and CPTAC.

Same design as external.py, but every encoder sees identical pixels: the PNG tile store (tumor_tiles.py, up to
4,000 tumor-map tiles per patient at 0.5 um/px). Kather-tile models use emb/<encoder>_rt05 (Kather regions re-cut at
0.5 um/px) and are applied to single store tiles on CPTAC.

  python src/external_fm.py --encoder uni2_h        # -> results/metrics/external_store_<encoder>.csv (+ scores pkl)
  python src/external_fm.py --summary               # -> results/external_store.md (all encoders, paired tests)
"""
import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(__file__))
import eval_lr  # noqa: E402
import pool  # noqa: E402
import protocols as P  # noqa: E402
import wsi_eval as W  # noqa: E402
from external import auc_ci  # noqa: E402

ROOT = P.ROOT
STORE = os.environ.get("MSI_TILES", "/mnt/slides/msimut/tumor_tiles")
SEEDS = (0, 1, 2, 3, 4)
MIN_TILES = 20
LABEL = {"flash": "Prov-GigaPath-flash", "uni2_h": "UNI2-h", "virchow2": "Virchow2", "gigapath": "Prov-GigaPath"}


def pooled(E, bags, keys, kind):
    out = []
    for k in keys:
        a = E[bags[k]].astype(np.float32)
        out.append({"mean": a.mean(0), "max": a.max(0), "mean+max": np.concatenate([a.mean(0), a.max(0)])}[kind])
    return np.stack(out)


def run(enc):
    D = P.Design()
    idx = pd.read_parquet(os.path.join(STORE, "index.parquet"))
    E = np.load(os.path.join(STORE, "emb", f"{enc}.f16.npy"), mmap_mode="r")
    E = np.asarray(E)
    cl = pd.read_csv(os.path.join(ROOT, "data/cptac/labels.csv"), index_col=0)["label"]
    pid = {p: i for i, p in enumerate(D.pats.patient_id)}
    tc = idx[idx.cohort == "tcga"]
    tb = {pid[p]: g.index.values for p, g in tc.groupby("patient") if len(g) >= MIN_TILES}
    cc = idx[idx.cohort == "cptac"]
    cpats = sorted(p for p, g in cc.groupby("patient") if len(g) >= MIN_TILES and p in cl.index)
    yc = cl.loc[cpats].values.astype(int)
    ckeys = np.arange(10000, 10000 + len(cpats))
    bags = {**tb, **{k: cc.index.values[(cc.patient == p).values] for k, p in zip(ckeys, cpats)}}
    tr = np.array(sorted(i for i in tb if i in set(D.train)))
    te = np.array(sorted(i for i in tb if i in set(D.test)))
    print(f"{enc}: TCGA TRAIN {len(tr)} TEST {len(te)}; CPTAC {len(cpats)} ({yc.sum()} MSI-H)", flush=True)
    scores = {}
    X = pooled(E, bags, np.concatenate([tr, te, ckeys]), "mean")
    p = eval_lr.fit_predict(X[:len(tr)], D.y[tr], X[len(tr):])
    scores["tumor_lr_mean"] = (p[:len(te)], p[len(te):])
    pa = np.mean([W.abmil_fit(E, bags, tr, D.y[tr], np.concatenate([te, ckeys]), s, "cuda:0") for s in SEEDS], 0)
    scores["tumor_abmil"] = (pa[:len(te)], pa[len(te):])
    # Kather-tile models (published TRAIN tiles re-cut at 0.5 um/px), applied to single CPTAC store tiles
    rt = np.asarray(pool.load_emb(f"{enc}_rt05"))
    pats = pool.patients()
    kb = {i: np.arange(s, e) for i, (s, e) in enumerate(zip(pats.start, pats.stop))}
    RT = np.concatenate([rt, E[cc.index.values]])
    cpos = {r: j for j, r in enumerate(cc.index.values)}
    kbags = {**kb, **{k: len(rt) + np.array([cpos[r] for r in bags[k]]) for k in ckeys}}
    ktr, kte = np.array(sorted(D.train)), np.array(sorted(D.test))
    for kind in ("mean", "max", "mean+max"):
        X = pooled(RT, kbags, np.concatenate([ktr, kte, ckeys]), kind)
        p = eval_lr.fit_predict(X[:len(ktr)], D.y[ktr], X[len(ktr):])
        scores[f"kather_lr_{kind}"] = (p[:len(kte)], p[len(kte):])
    pk = np.mean([W.abmil_fit(RT, kbags, ktr, D.y[ktr], np.concatenate([kte, ckeys]), s, "cuda:0") for s in SEEDS], 0)
    scores["kather_abmil"] = (pk[:len(kte)], pk[len(kte):])
    pickle.dump(dict(scores=scores, yc=yc, y_te=D.y[te], y_kte=D.y[kte], cpats=cpats, te=te, kte=kte),
                open(os.path.join(ROOT, f"results/preds/external_store_{enc}.pkl"), "wb"))
    rows = []
    for name, (st, sc) in scores.items():
        yt = D.y[kte] if name.startswith("kather") else D.y[te]
        a1, l1, h1, _ = auc_ci(yt, np.asarray(st, float))
        a2, l2, h2, _ = auc_ci(yc, np.asarray(sc, float))
        rows.append(dict(encoder=enc, model=name, tcga_test_auc=a1, tcga_ci_low=l1, tcga_ci_high=h1, cptac_auc=a2,
                         cptac_ci_low=l2, cptac_ci_high=h2))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(ROOT, f"results/metrics/external_store_{enc}.csv"), index=False)
    print(res.round(3).to_string(index=False))


def run_surgen(enc):
    """Same TCGA TRAIN models (deterministic, refit), applied to SurGen SR386 subset (label: MMR loss by IHC)."""
    D = P.Design()
    idx = pd.read_parquet(os.path.join(STORE, "index.parquet"))
    E = np.asarray(np.load(os.path.join(STORE, "emb", f"{enc}.f16.npy"), mmap_mode="r"))
    s2 = os.environ.get("MSI_TILES_SURGEN", "/mnt/slides/msimut/tumor_tiles_surgen")
    gidx = pd.read_parquet(os.path.join(s2, "index.parquet"))
    G = np.asarray(np.load(os.path.join(s2, "emb", f"{enc}.f16.npy"), mmap_mode="r"))
    lab = pd.read_csv(os.path.join(ROOT, "data/surgen/subset.csv"), dtype={"case_id": str}).set_index("case_id")
    pid = {p: i for i, p in enumerate(D.pats.patient_id)}
    tc = idx[idx.cohort == "tcga"]
    tb = {pid[p]: g.index.values for p, g in tc.groupby("patient") if len(g) >= MIN_TILES}
    gpats = sorted(p for p, g in gidx.groupby("patient") if len(g) >= MIN_TILES and p in lab.index)
    yg = lab.loc[gpats, "mmr_loss_binary"].values.astype(int)
    gkeys = np.arange(20000, 20000 + len(gpats))
    EE = np.concatenate([E, G])
    off = len(E)
    bags = {**tb, **{k: gidx.index.values[(gidx.patient == p).values] + off for k, p in zip(gkeys, gpats)}}
    tr = np.array(sorted(i for i in tb if i in set(D.train)))
    print(f"{enc}: SurGen {len(gpats)} cases ({yg.sum()} dMMR)", flush=True)
    scores = {}
    X = pooled(EE, bags, np.concatenate([tr, gkeys]), "mean")
    scores["tumor_lr_mean"] = eval_lr.fit_predict(X[:len(tr)], D.y[tr], X[len(tr):])
    scores["tumor_abmil"] = np.mean([W.abmil_fit(EE, bags, tr, D.y[tr], gkeys, s, "cuda:0") for s in SEEDS], 0)
    rt = np.asarray(pool.load_emb(f"{enc}_rt05"))
    pats = pool.patients()
    kb = {i: np.arange(s_, e_) for i, (s_, e_) in enumerate(zip(pats.start, pats.stop))}
    RT = np.concatenate([rt, G])
    kbags = {**kb, **{k: len(rt) + (bags[k] - off) for k in gkeys}}
    ktr = np.array(sorted(D.train))
    for kind in ("mean", "max", "mean+max"):
        X = pooled(RT, kbags, np.concatenate([ktr, gkeys]), kind)
        scores[f"kather_lr_{kind}"] = eval_lr.fit_predict(X[:len(ktr)], D.y[ktr], X[len(ktr):])
    scores["kather_abmil"] = np.mean([W.abmil_fit(RT, kbags, ktr, D.y[ktr], gkeys, s, "cuda:0") for s in SEEDS], 0)
    pickle.dump(dict(scores=scores, y=yg, cases=gpats),
                open(os.path.join(ROOT, f"results/preds/external_surgen_{enc}.pkl"), "wb"))
    rows = []
    for name, sc in scores.items():
        a2, l2, h2, _ = auc_ci(yg, np.asarray(sc, float))
        rows.append(dict(encoder=enc, model=name, surgen_auc=a2, surgen_ci_low=l2, surgen_ci_high=h2, n=len(yg),
                         pos=int(yg.sum())))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(ROOT, f"results/metrics/external_surgen_{enc}.csv"), index=False)
    print(res.round(3).to_string(index=False))


def summary():
    encs = [e for e in LABEL if os.path.exists(os.path.join(ROOT, f"results/preds/external_store_{e}.pkl"))]
    S = {e: pickle.load(open(os.path.join(ROOT, f"results/preds/external_store_{e}.pkl"), "rb")) for e in encs}
    lines = ["# Whole-tumor vs Kather-tile models across encoders: TCGA TRAIN -> TCGA TEST / CPTAC-COAD", "",
             "Identical tile store for all encoders (<= 4,000 tumor tiles per patient, 0.5 µm/px, PNG). Models fit on "
             "TCGA TRAIN only. ROC-AUC [95% CI, stratified bootstrap].", "",
             "| Encoder | Model | TCGA TEST | CPTAC |", "|---|---|---|---|"]
    for e in encs:
        r = pd.read_csv(os.path.join(ROOT, f"results/metrics/external_store_{e}.csv"))
        for x in r.itertuples():
            lines.append(f"| {LABEL[e]} | {x.model} | {x.tcga_test_auc:.3f} [{x.tcga_ci_low:.3f}, {x.tcga_ci_high:.3f}] | "
                         f"{x.cptac_auc:.3f} [{x.cptac_ci_low:.3f}, {x.cptac_ci_high:.3f}] |")
    lines += ["", "Paired CPTAC differences (same resamples):", "", "| Encoder | Comparison | Δ AUC [95% CI] |",
              "|---|---|---|"]
    yc = S[encs[0]]["yc"]
    bi = P.boot_indices(yc)
    for e in encs:
        sc = {k: np.asarray(v[1], float) for k, v in S[e]["scores"].items()}
        b = {k: P.fast_auc(yc[bi], v[bi]) for k, v in sc.items()}
        for a_, b_ in (("tumor_abmil", "kather_lr_max"), ("tumor_abmil", "kather_abmil"),
                       ("kather_lr_mean", "kather_lr_max"), ("tumor_abmil", "tumor_lr_mean")):
            d = b[a_] - b[b_]
            dv = roc_auc_score(yc, sc[a_]) - roc_auc_score(yc, sc[b_])
            lines.append(f"| {LABEL[e]} | {a_} - {b_} | {dv:+.3f} [{np.percentile(d, 2.5):+.3f}, "
                         f"{np.percentile(d, 97.5):+.3f}] |")
    sg = [e for e in encs if os.path.exists(os.path.join(ROOT, f"results/metrics/external_surgen_{e}.csv"))]
    if sg:
        lines += ["", "## SurGen SR386 subset (MMR loss by IHC), models fit on TCGA TRAIN only", "",
                  "| Encoder | Model | SurGen AUC [95% CI] |", "|---|---|---|"]
        for e in sg:
            r = pd.read_csv(os.path.join(ROOT, f"results/metrics/external_surgen_{e}.csv"))
            lines += [f"| {LABEL[e]} | {x.model} | {x.surgen_auc:.3f} [{x.surgen_ci_low:.3f}, {x.surgen_ci_high:.3f}] "
                      f"(n={x.n}, {x.pos} dMMR) |" for x in r.itertuples()]
    open(os.path.join(ROOT, "results/external_store.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--surgen", action="store_true")
    a = ap.parse_args()
    summary() if a.summary else (run_surgen(a.encoder) if a.surgen else run(a.encoder))
