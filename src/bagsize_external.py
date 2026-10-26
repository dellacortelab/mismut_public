"""Is the external collapse of max pooling a learned count shortcut or a bag-size extrapolation artefact?

Units are consistent everywhere: single 224 px whole-tumor tiles at 0.5 um/px from the lossless tile store.
Training (TCGA TRAIN) reproduces the release's confound: each patient's bag is a random sample of exactly as many
store tiles as that patient has in the Kather release (MSIMUT ~1,198, MSS ~211). A deconfounded reference trains every
patient on a fixed 200 tiles. LR on mean and max pooling (eval_lr configuration).
Test bags (TCGA TEST, CPTAC-COAD, SurGen):
  full        all store tiles of the patient (<= 4,000)
  tcga_like   n drawn per patient from the TCGA TRAIN count distribution, independent of the label
  n200, n1000 fixed sizes
Ten draws per condition (bag sampling); AUC mean and range over draws.

  python src/bagsize_external.py --encoder uni2_h   # -> results/metrics/bagsize_<encoder>.csv
  python src/bagsize_external.py --summary          # -> results/bagsize.md
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(__file__))
import eval_lr  # noqa: E402
import protocols as P  # noqa: E402

ROOT = P.ROOT
STORE = os.environ.get("MSI_TILES", "/mnt/slides/msimut/tumor_tiles")
STORE_SG = os.environ.get("MSI_TILES_SURGEN", "/mnt/slides/msimut/tumor_tiles_surgen")
DRAWS = 10
LABEL = {"flash": "Prov-GigaPath-flash", "uni2_h": "UNI2-h", "virchow2": "Virchow2"}


def pool_bag(E, rows, kind):
    a = E[rows].astype(np.float32)
    return a.mean(0) if kind == "mean" else a.max(0)


def run(enc):
    D = P.Design()
    idx = pd.read_parquet(os.path.join(STORE, "index.parquet"))
    E = np.asarray(np.load(os.path.join(STORE, "emb", f"{enc}.f16.npy"), mmap_mode="r"))
    gidx = pd.read_parquet(os.path.join(STORE_SG, "index.parquet"))
    G = np.asarray(np.load(os.path.join(STORE_SG, "emb", f"{enc}.f16.npy"), mmap_mode="r"))
    pid = {p: i for i, p in enumerate(D.pats.patient_id)}
    tc = idx[idx.cohort == "tcga"]
    tcga_rows = {pid[p]: g.index.values for p, g in tc.groupby("patient")}
    tr = np.array(sorted(i for i in tcga_rows if i in set(D.train)))
    te = np.array(sorted(i for i in tcga_rows if i in set(D.test)))
    kather_n = D.n  # tiles per patient in the release
    train_dist = kather_n[D.train]
    cl = pd.read_csv(os.path.join(ROOT, "data/cptac/labels.csv"), index_col=0)["label"]
    cc = idx[idx.cohort == "cptac"]
    cptac = {p: g.index.values for p, g in cc.groupby("patient") if p in cl.index and len(g) >= 20}
    sl = pd.read_csv(os.path.join(ROOT, "data/surgen/subset.csv"), dtype={"case_id": str}).set_index("case_id")
    surgen = {p: g.index.values for p, g in gidx.groupby("patient") if p in sl.index and len(g) >= 20}
    cohorts = {"tcga_test": (E, {i: tcga_rows[i] for i in te}, {i: D.y[i] for i in te}),
               "cptac": (E, cptac, {p: int(cl[p]) for p in cptac}),
               "surgen": (G, surgen, {p: int(sl.loc[p, "mmr_loss_binary"]) for p in surgen})}
    rows = []
    for train_mode in ("kather_counts", "fixed200"):
        for kind in ("max", "mean"):
            rng = np.random.default_rng(0)
            Xtr = np.stack([pool_bag(E, rng.choice(tcga_rows[i], min(len(tcga_rows[i]),
                                                                       kather_n[i] if train_mode == "kather_counts" else 200),
                                                   replace=False), kind) for i in tr])
            model = eval_lr.lr().fit(Xtr, D.y[tr])
            for cname, (EE, bags, labs) in cohorts.items():
                keys = sorted(bags)
                y = np.array([labs[k] for k in keys])
                for cond in ("full", "tcga_like", "n200", "n1000"):
                    aucs, frac_hi = [], []
                    for d in range(DRAWS if cond != "full" else 1):
                        r2 = np.random.default_rng([1, d])
                        X = []
                        for k in keys:
                            b = bags[k]
                            if cond == "full":
                                n = len(b)
                            elif cond == "tcga_like":
                                n = int(r2.choice(train_dist))
                            else:
                                n = 200 if cond == "n200" else 1000
                            X.append(pool_bag(EE, r2.choice(b, min(n, len(b)), replace=False), kind))
                        p = model.predict_proba(np.stack(X))[:, 1]
                        aucs.append(roc_auc_score(y, p))
                        frac_hi.append((p > 0.99).mean())
                    rows.append(dict(encoder=enc, train=train_mode, pooling=kind, cohort=cname, test_bags=cond,
                                     auc=np.mean(aucs), auc_min=np.min(aucs), auc_max=np.max(aucs),
                                     frac_p_gt_099=np.mean(frac_hi), n=len(y), pos=int(y.sum())))
    r = pd.DataFrame(rows)
    r.to_csv(os.path.join(ROOT, f"results/metrics/bagsize_{enc}.csv"), index=False)
    print(r.round(3).to_string(index=False))


def summary():
    fs = [os.path.join(ROOT, f"results/metrics/bagsize_{e}.csv") for e in LABEL]
    r = pd.concat([pd.read_csv(f) for f in fs if os.path.exists(f)])
    lines = ["# External collapse of max pooling: shortcut or bag-size artefact?", "",
             "Training on TCGA TRAIN with each patient's Kather tile count (confounded) or 200 tiles (deconfounded); "
             "single whole-tumor tiles at 0.5 µm/px throughout. AUC (mean over 10 bag draws; full = one draw).", ""]
    for train in ("kather_counts", "fixed200"):
        lines += [f"## Training bags: {train}", "",
                  "| Encoder | Pooling | Cohort | full | TCGA-like n | n=200 | n=1000 | P>0.99 (full) |",
                  "|---|---|---|---|---|---|---|---|"]
        for (e, kind, c), g in r[r.train == train].groupby(["encoder", "pooling", "cohort"], sort=False):
            v = g.set_index("test_bags")
            lines.append(f"| {LABEL[e]} | {kind} | {c} | {v.auc['full']:.3f} | {v.auc['tcga_like']:.3f} | "
                         f"{v.auc['n200']:.3f} | {v.auc['n1000']:.3f} | {v.frac_p_gt_099['full']:.2f} |")
        lines.append("")
    open(os.path.join(ROOT, "results/bagsize.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder")
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()
    summary() if a.summary else run(a.encoder)
