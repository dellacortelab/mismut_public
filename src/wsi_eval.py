"""E8: where is the MSI signal? Region-wise evaluation on whole-slide tiles (wsi_regions.py embeddings).

Tumor map. Kather tile coverage cannot define tumor directly: MSS TRAIN slides were randomly undersampled, so
their Kather tiles are sparse dots inside the tumor and the rest of the tumor would be mislabelled "margin" for MSS
patients only (label leakage into the region definition). The NCT-CRC-HE tissue probe does not transfer to raw
(not Macenko-normalized) slide tiles (it misses most Kather tumor tiles on many slides). Instead we learn Kather's
tumor definition in-domain: on slides with complete Kather coverage (all except TRAIN MSS), tiles inside Kather
tiles are positives and tissue > NEG_UM away from any Kather tile negatives; an LR on the encoder embeddings is
cross-fitted (5 folds grouped by patient), so each complete slide is scored by a model that never saw it, and the
undersampled slides by the mean of the five fold models. The target is tumor vs non-tumor, never MSI. p(tumor) is
smoothed over the 3x3 tile neighbourhood; tumor = smoothed p >= 0.5.
  margin  = non-tumor tissue with centre within MARGIN_UM of a tumor tile
  distant = remaining tissue
NCT probe probabilities are kept only as secondary tissue-fraction features.

Region sets: tumor, margin, distant, tumor+margin, non-tumor, whole. Patients enter every analysis only if they
have >= MIN_TILES tiles in each of tumor, margin and distant (common cohort; paired comparisons need identical
test patients).

Learners: LR on mean pooling (eval_lr config) and gated ABMIL (abmil.py architecture and training; bags stored
fp16 on the GPU). Protocols: E1rep (10 x 5 stratified CV over the cohort, seed = repeat) and E2 (TRAIN -> TEST,
seeds 0,1,2 for ABMIL). Also log1p(tile count) per region alone. Paired stratified bootstrap on TEST.

  python src/wsi_eval.py --encoder flash [--jobs 3]
Outputs: results/metrics/wsi_<encoder>_*.csv, results/preds/wsi_<encoder>.parquet, results/wsi_<encoder>.md
"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from joblib import Parallel, delayed
from scipy.ndimage import distance_transform_edt, uniform_filter
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(__file__))
import abmil as A  # noqa: E402
import eval_lr  # noqa: E402
import protocols as P  # noqa: E402

ROOT = P.ROOT
WSI = os.environ.get("MSI_WSI", "data/external/wsi")
PLAN = os.path.join(ROOT, "data/wsi/plan.parquet")
TILE_UM, MARGIN_UM, MIN_TILES, NEG_UM = 112.0, 1000.0, 20, 300.0
SETS = {"tumor": (0,), "margin": (1,), "distant": (2,), "tumor+margin": (0, 1), "non-tumor": (1, 2),
        "whole": (0, 1, 2)}
TISSUE = ["ADI", "BACK", "DEB", "LYM", "MUC", "MUS", "NORM", "STR", "TUM"]


def tumor_probability(t, emb, seed=0):
    """Cross-fitted in-domain tumor classifier (see module docstring). Returns p(tumor) per tile and fold info."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    reg = pd.read_csv(os.path.join(ROOT, "results/tile_slide_registration.csv.gz")).groupby("slide").first()
    sparse = set(reg.index[(reg.split == "TRAIN") & (reg.label == 0)])
    dense_bc = np.array(sorted(set(t.barcode) - sparse))
    rng = np.random.default_rng(seed)
    samp = []
    for bc in dense_bc:
        g = t[t.barcode == bc]
        pos = g.index[g.region == 0].values
        neg = g.index[(g.region != 0) & (g.dist_um > NEG_UM)].values
        samp += [rng.choice(pos, min(1500, len(pos)), replace=False), rng.choice(neg, min(1500, len(neg)), replace=False)]
    samp = np.concatenate(samp)
    ys = (t.region.values[samp] == 0).astype(int)
    groups = t.barcode.values[samp].astype(str)
    pat = np.array([g[:12] for g in groups])
    p = np.zeros(len(t))
    models = []
    for k, (a, b) in enumerate(GroupKFold(5).split(samp, ys, pat)):
        clf = make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=2000, random_state=seed))
        clf.fit(emb[samp[a]].astype(np.float32), ys[a])
        models.append(clf)
        held = np.unique(pat[b])
        rows = np.where(np.isin(t.barcode.str[:12].values, held))[0]
        p[rows] = clf.predict_proba(emb[rows].astype(np.float32))[:, 1]
    rows = np.where(t.barcode.isin(sparse).values)[0]
    p[rows] = np.mean([m.predict_proba(emb[rows].astype(np.float32))[:, 1] for m in models], axis=0)
    return p


def load(encoder):
    """All tiles with plan coordinates, embeddings, tissue probabilities and data-driven regions."""
    import joblib
    plan = pd.read_parquet(PLAN)
    probe = joblib.load(os.path.join(ROOT, f"data/nct/probe_{encoder}.joblib"))
    parts, embs = [], []
    for bc, g in plan.groupby("barcode", sort=True):
        z = np.load(os.path.join(WSI, encoder, bc + ".npz"))
        assert len(z["emb"]) == len(g) and np.allclose(z["x0"], g.x0.values)
        g = g.copy()
        pr = probe.predict_proba(z["emb"].astype(np.float32))
        for k, c in enumerate(TISSUE):
            g["p_" + c] = pr[:, k]
        parts.append(g)
        embs.append(z["emb"])
    t = pd.concat(parts, ignore_index=True)
    emb = np.concatenate(embs)
    t["p_tumor"] = tumor_probability(t, emb)
    out = []
    for bc, g in t.groupby("barcode", sort=True):
        H, W = g.ty.max() + 1, g.tx.max() + 1
        grid = np.zeros((H, W), np.float32)
        have = np.zeros((H, W), np.float32)
        grid[g.ty, g.tx] = g.p_tumor.values
        have[g.ty, g.tx] = 1
        sm = uniform_filter(grid, 3, mode="constant") / np.maximum(uniform_filter(have, 3, mode="constant"), 1e-6)
        tum = np.zeros((H, W), bool)
        tum[g.ty, g.tx] = sm[g.ty, g.tx] >= 0.5
        dist = distance_transform_edt(~tum) * TILE_UM if tum.any() else np.full((H, W), np.inf)
        out.append(pd.DataFrame({"tumor_map": tum[g.ty, g.tx], "dist_tumor_um": dist[g.ty, g.tx]}, index=g.index))
    t = t.join(pd.concat(out))
    t["region_data"] = np.where(t.tumor_map, 0, np.where(t.dist_tumor_um <= MARGIN_UM, 1, 2))
    return t, emb


def validate_map(t):
    """Agreement of the data-driven tumor map with Kather's tiles (plan.region == 0 means inside a Kather tile)."""
    reg = pd.read_csv(os.path.join(ROOT, "results/tile_slide_registration.csv.gz")).groupby("slide").first()
    t = t.assign(kather=t.region == 0)
    rows = []
    for bc, g in t.groupby("barcode"):
        rows.append(dict(barcode=bc, split=reg.loc[bc, "split"], label=reg.loc[bc, "label"],
                         kather_tiles=int(g.kather.sum()), tumor_map_tiles=int(g.tumor_map.sum()),
                         kather_in_map=g.tumor_map[g.kather].mean() if g.kather.any() else np.nan,
                         map_covered_by_kather=g.kather[g.tumor_map].mean() if g.tumor_map.any() else np.nan))
    v = pd.DataFrame(rows)
    v["group"] = v.split + "_" + v.label.map({1: "MSIMUT", 0: "MSS"})
    return v


def patient_bags(t, D, cohort_rule="common"):
    pid_of = {p: i for i, p in enumerate(D.pats.patient_id)}
    t = t.assign(pid=t.barcode.str[:12].map(pid_of))
    assert t.pid.notna().all()
    cnt = t.groupby(["pid", "region_data"]).size().unstack(fill_value=0).reindex(range(len(D.y)), fill_value=0)
    need = [0, 1, 2] if cohort_rule == "common" else [0]
    cohort = np.where((cnt.reindex(columns=[0, 1, 2], fill_value=0)[need] >= MIN_TILES).all(axis=1))[0]
    bags = {}
    for name, regs in SETS.items():
        sel = t[t.region_data.isin(regs)]
        idx = sel.groupby("pid").apply(lambda g: g.index.values, include_groups=False)
        bags[name] = {int(p): idx[p] for p in cohort}
    return t, cohort, cnt, bags


def jobs_for(D, cohort):
    y = D.y[cohort]
    out = []
    for r in range(10):
        skf = StratifiedKFold(5, shuffle=True, random_state=P.SEED + r)
        for k, (a, b) in enumerate(skf.split(cohort, y)):
            out.append(dict(protocol="E1rep", cv_repeat=r, fold=k, tr=cohort[a], te=cohort[b], seeds=(r,)))
    tr = np.intersect1d(cohort, D.train)
    te = np.intersect1d(cohort, D.test)
    out.append(dict(protocol="E2", cv_repeat=0, fold=-1, tr=tr, te=te, seeds=(0, 1, 2)))
    return out


def rows_df(meta, te, p):
    return pd.DataFrame({**{k: [v] * len(te) for k, v in meta.items()}, "pid": te, "p": p})


def abmil_fit(emb, bags, tr, ytr, te, seed, device):
    """abmil.fit_predict with bags held in fp16 on the GPU (whole-slide bags are ~10k tiles per patient)."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    rows = np.concatenate([bags[i] for i in tr])
    samp = np.sort(rng.choice(rows, min(50_000, len(rows)), replace=False))
    xs = emb[samp].astype(np.float32)
    mu, sd = torch.tensor(xs.mean(0), device=device), torch.tensor(xs.std(0) + 1e-6, device=device)
    X = {i: torch.from_numpy(emb[bags[i]]).to(device) for i in np.concatenate([tr, te])}

    def bag(i):
        return (X[i].float() - mu) / sd

    model = A.GatedABMIL(emb.shape[1]).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=A.LR, weight_decay=A.WD)
    lossf = nn.BCEWithLogitsLoss()
    w = np.where(ytr == 1, 1.0 / (ytr == 1).sum(), 1.0 / (ytr == 0).sum())
    w /= w.sum()
    model.train()
    for _ in range(A.EPOCHS):
        for k in rng.choice(len(tr), len(tr), replace=True, p=w):
            opt.zero_grad()
            lossf(model(bag(tr[k])), torch.tensor(float(ytr[k]), device=device)).backward()
            opt.step()
    model.eval()
    with torch.inference_mode():
        return np.array([torch.sigmoid(model(bag(i))).item() for i in te])


def run_abmil(emb, bags, name, job, seed, y, slot):
    device = A.resolve_device("auto", slot)
    p = abmil_fit(emb, bags[name], job["tr"], y[job["tr"]], job["te"], seed, device)
    meta = dict(aggregation="abmil", scheme=name, N=0, repeat=0, protocol=job["protocol"],
                cv_repeat=job["cv_repeat"], fold=job["fold"], seed=seed)
    return rows_df(meta, job["te"], p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder", default="flash")
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--sets", default=",".join(SETS))
    ap.add_argument("--no-abmil", action="store_true")
    ap.add_argument("--save-regions", action="store_true", help="write data/wsi/regions_<encoder>.parquet and exit")
    ap.add_argument("--cohort", default="common", choices=["common", "tumor"],
                    help="common: >= MIN_TILES in tumor, margin and distant; tumor: >= MIN_TILES tumor tiles")
    a = ap.parse_args()
    t0 = time.time()
    D = P.Design()
    t, emb = load(a.encoder)
    if a.save_regions:
        t[["barcode", "tx", "ty", "x0", "y0", "mpp0", "region", "region_data", "p_tumor", "tumor_map",
           "dist_tumor_um"]].to_parquet(os.path.join(ROOT, f"data/wsi/regions_{a.encoder}.parquet"), index=False)
        print("saved regions", t.region_data.value_counts().sort_index().to_dict())
        return
    v = validate_map(t)
    v.to_csv(os.path.join(ROOT, f"results/metrics/wsi_{a.encoder}_tumormap_validation.csv"), index=False)
    print(v.groupby("group")[["kather_in_map", "map_covered_by_kather"]].median().round(3), flush=True)
    t, cohort, cnt, bags = patient_bags(t, D, a.cohort)
    print(f"cohort: {len(cohort)} patients with >= {MIN_TILES} tiles in each region "
          f"(TRAIN {np.isin(cohort, D.train).sum()}, TEST {np.isin(cohort, D.test).sum()})", flush=True)
    jobs = jobs_for(D, cohort)
    sets = a.sets.split(",")
    preds = []
    # LR: mean pooling, mean tissue probabilities, and log1p(count), per region set
    for name in sets:
        B = bags[name]
        mean = np.stack([emb[B[i]].astype(np.float32).mean(0) for i in cohort])
        tis = np.stack([t.loc[B[i], ["p_" + c for c in TISSUE]].values.mean(0) for i in cohort])
        cnt_f = np.log1p([len(B[i]) for i in cohort])[:, None]
        X = {"mean": mean, "tissue_frac": tis, "count": cnt_f}
        pos = {p: k for k, p in enumerate(cohort)}
        for agg, M in X.items():
            for j in jobs:
                tr, te = [pos[i] for i in j["tr"]], [pos[i] for i in j["te"]]
                p = eval_lr.fit_predict(M[tr], D.y[j["tr"]], M[te])
                preds.append(rows_df(dict(aggregation=agg, scheme=name, N=0, repeat=0, protocol=j["protocol"],
                                          cv_repeat=j["cv_repeat"], fold=j["fold"], seed=-1), j["te"], p))
    print(f"LR done {time.time() - t0:.0f}s", flush=True)
    if not a.no_abmil:
        tasks = [(name, j, s) for name in sets for j in jobs for s in j["seeds"]]
        out = Parallel(n_jobs=a.jobs, backend="threading")(
            delayed(run_abmil)(emb, bags, name, j, s, D.y, i) for i, (name, j, s) in enumerate(tasks))
        preds += out
        print(f"ABMIL done {time.time() - t0:.0f}s", flush=True)
    pr = pd.concat(preds, ignore_index=True)
    tag = a.encoder if a.cohort == "common" else f"{a.encoder}_tumorcohort"
    pr.insert(0, "encoder", f"wsi_{tag}")
    pr.to_parquet(os.path.join(ROOT, f"results/preds/wsi_{tag}.parquet"), index=False)
    m = P.summarize(pr, D.y)
    m.to_csv(os.path.join(ROOT, f"results/metrics/wsi_{tag}_metrics.csv"), index=False)
    # paired E2 bootstrap between region sets, same aggregation
    rows = []
    e2 = pr[pr.protocol == "E2"]
    for agg in ("mean", "abmil", "tissue_frac"):
        base = {s: e2[(e2.aggregation == agg) & (e2.scheme == s)] for s in sets}
        for x, yb in (("margin", "tumor"), ("tumor+margin", "tumor"), ("whole", "tumor"), ("whole", "tumor+margin"),
                      ("non-tumor", "tumor"), ("distant", "tumor"), ("margin", "distant")):
            if x not in base or yb not in base or not len(base[x]) or not len(base[yb]):
                continue
            ba, bb = P.boot_unit_aucs(base[x], D.y)[1], P.boot_unit_aucs(base[yb], D.y)[1]
            va, vb = (np.mean([P.metrics(D.y[u.pid.values], u.p.values)["roc_auc"] for _, u in
                               base[k].groupby(["repeat", "seed"])]) for k in (x, yb))
            d = ba - bb
            rows.append(dict(aggregation=agg, comparison=f"{x} - {yb}", auc_a=va, auc_b=vb, delta=va - vb,
                             CI_low=np.percentile(d, 2.5), CI_high=np.percentile(d, 97.5),
                             p_boot=min(1.0, max(2 * min((d <= 0).mean(), (d >= 0).mean()), 1 / len(d)))))
    pe = pd.DataFrame(rows)
    pe.to_csv(os.path.join(ROOT, f"results/metrics/wsi_{tag}_paired_e2.csv"), index=False)
    auc = m[m.metric == "roc_auc"].pivot_table(index=["aggregation", "scheme"], columns="protocol", values="value")
    ci = m[(m.metric == "roc_auc") & (m.protocol == "E2")].set_index(["aggregation", "scheme"])[["CI_low", "CI_high"]]
    print(auc.join(ci).round(3).to_string())
    print(pe.round(3).to_string(index=False))
    json.dump({"encoder": a.encoder, "cohort": len(cohort), "tiles": int(len(t)), "seconds": time.time() - t0,
               "region_tiles": t.region_data.value_counts().sort_index().tolist()},
              open(os.path.join(ROOT, f"results/metrics/wsi_{tag}_info.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
