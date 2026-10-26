"""Blinded pathologist review of attention regions: unblind the platform export and produce Fig. 4 and the text numbers.

  python src/review_analysis.py from-deposit <unzipped review.zip>   # Zenodo 10.5281/zenodo.23666864 -> data/review/
  python src/review_analysis.py simulate                      # fake export for layout work -> data/review/export_simulated.csv
  python src/review_analysis.py analyze data/review/export_simulated.csv --simulated
  python src/review_analysis.py analyze <survey_export.csv>   # real run

Input: the CSV from GET /gi/api/v1/annotate/survey/export/csv?project_id=<msi-attention-review>. Schema answers sit in
the platform's fixed columns: q1_growth_pattern = main tissue, q2_lumen_status = inflammation, q3_gland_separation =
tumor morphology, q8_cytology = tumor cytology (added on 22 Oct 2026 after 119 responses, so analysed on the
responses that have it), q4_gp_subtype = "looks MSI-associated?", q7_notes = free text. annotation_id = "<slug>:<GI-xx>:<Rk>".
Unblinding: data/review/slide_codes.csv (GI code -> barcode) and data/review/regions_key.csv (region kind, patient group).

Analyses (pre-specified):
  a) share of responses ticking each tissue / inflammation category, high-attention vs control regions (both raters
     pooled), Fisher's exact test per category, Benjamini-Hochberg across categories
  b) "MSI-associated?" yes / unsure / no by region type and patient group (MSI; MSS called MSI by the model; MSS)
  c) example high-attention regions with the longest free-text descriptions (dropped: the reviewers wrote almost no
     free text, 1 note in 354 responses)
  d) inter-rater agreement: Cohen's kappa for main tissue and MSI call, mean kappa over inflammation categories
Added after unblinding (exploratory; reported as such):
  e) high vs control restricted to responses calling the region tumor epithelium (controls were often not tumor)
  f) high-attention regions by patient group: per-slide share of responses naming each feature, MSI (14 slides) vs
     MSS (6 slides) by Mann-Whitney and slide-level AUC, both raters pooled and per rater
  g) attention-gated KongNet counts for every cell class vs all tiles, a label-free epithelial gate and a
     permuted-label attention gate (results/metrics/gate_controls.csv from src/gate_controls.py)
  h) molecular annotation of the review patients (data/review/molecular.csv: TMB, MANTIS score and POLE mutations of
     the TCGA PanCancer Atlas, cBioPortal study coadread_tcga_pan_can_atlas_2018)
Outputs: results/review/summary.md, results/review/responses_unblinded.csv,
results/fig10_review.{png,pdf} (watermarked when --simulated).
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REV = os.path.join(ROOT, "data/review")
OUT = os.path.join(ROOT, "results/review")
SLIDES = os.environ.get("MSI_SLIDES", "/mnt/slides/msimut") + "/tcga_crc_dx"
COLS = {"q1": "q1_growth_pattern", "q2": "q2_lumen_status", "q3": "q3_gland_separation", "q4": "q4_gp_subtype",
        "q8": "q8_cytology", "notes": "q7_notes"}
RATERS = ("rater_A", "rater_B")  # platform usernames are anonymised in sorted order, as in the Zenodo deposit
GROUP_LABEL = {"MSI_correct": "MSI", "MSS_false_positive": "MSS, called MSI", "MSS_correct": "MSS"}


def schema():
    p = os.path.join(REV, "project.json")
    s = json.load(open(p))["survey_schema"] if os.path.exists(p) else json.load(open(os.path.join(REV, "survey_schema.json")))
    return {q["id"]: q for q in s["questions"]}


def key():
    k = pd.read_csv(os.path.join(REV, "regions_key.csv"))
    codes = pd.read_csv(os.path.join(REV, "slide_codes.csv"))[["barcode", "code"]]
    return k.drop(columns="code", errors="ignore").merge(codes, on="barcode")


def cmd_from_deposit(a):
    """Rebuild the platform-format inputs from the Zenodo deposit (review.zip unpacked in a.dir)."""
    import shutil
    os.makedirs(REV, exist_ok=True)
    r = pd.read_csv(os.path.join(a.dir, "regions.csv")).rename(columns={"blind_code": "code", "region": "region_id"})
    r.to_csv(os.path.join(REV, "regions_key.csv"), index=False)
    r[["barcode", "code"]].drop_duplicates().to_csv(os.path.join(REV, "slide_codes.csv"), index=False)
    shutil.copy(os.path.join(a.dir, "survey_schema.json"), os.path.join(REV, "survey_schema.json"))
    d = pd.read_csv(os.path.join(a.dir, "responses.csv"), dtype=str, keep_default_na=False)
    inv = {"main_tissue": COLS["q1"], "inflammation": COLS["q2"], "morphology": COLS["q3"], "msi_impression": COLS["q4"],
           "cytology": COLS["q8"], "notes": COLS["notes"]}
    e = d.rename(columns=inv).rename(columns={"rater": "username"})
    e["annotation_id"] = "msi-attention-review:" + e.blind_code + ":" + e.region
    e["skipped"] = "False"
    out = os.path.join(REV, "export_deposit.csv")
    e.drop(columns=["blind_code", "region"]).to_csv(out, index=False)
    print(f"{len(e)} responses -> {out}; next: python src/fetch_molecular.py, then analyze {out}")


def cmd_simulate(a):
    """Random answers in the platform's export format. Placeholder only: no relation to the real review."""
    rng = np.random.default_rng(a.seed)
    Q = schema()
    k = key()
    rows = []
    for rater in RATERS:
        for r in k.itertuples():
            opts = {q: [o["value"] for o in Q[q]["options"]] for q in Q}
            q2 = [o for o in opts["q2"] if o != "none" and rng.random() < 0.3] or ["none"]
            q8 = ([o for o in opts["q8"][:4] if rng.random() < 0.3] or ["none"]) if rng.random() < 0.4 else []
            q3 = [o for o in opts["q3"] if rng.random() < 0.2]
            rows.append({"username": rater, "slide_filename": f"{r.code}.svs", "slide_id": "",
                         "annotation_id": f"msi-attention-review:{r.code}:{r.region_id}", "original_gp": "MSI review region",
                         "updated_gp": "", COLS["q1"]: rng.choice(opts["q1"], p=_p(len(opts["q1"]), rng)),
                         COLS["q2"]: "|".join(q2), COLS["q3"]: "|".join(q3), COLS["q8"]: "|".join(q8),
                         COLS["q4"]: rng.choice(opts["q4"]), COLS["notes"]: "[simulated description]",
                         "skipped": False})
    out = os.path.join(REV, "export_simulated.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print("wrote", out, len(rows), "rows (SIMULATED)")


def _p(n, rng):
    w = rng.dirichlet(np.ones(n) * 2)
    return w / w.sum()


def load(export):
    e = pd.read_csv(export, dtype=str, keep_default_na=False)
    e = e[e.skipped.str.lower() != "true"]
    if not set(e.username) <= set(RATERS):
        e["username"] = e.username.map(dict(zip(sorted(e.username.unique()), RATERS)))
    e["code"] = e.annotation_id.str.split(":").str[1]
    e["region_id"] = e.annotation_id.str.split(":").str[2]
    d = e.merge(key(), on=["code", "region_id"], how="left", validate="many_to_one")
    assert d.kind.notna().all(), "responses that do not match the key"
    if COLS["q8"] not in d:
        d[COLS["q8"]] = ""
    for q in ("q1", "q2", "q3", "q4", "q8"):
        d[q] = d[COLS[q]].apply(lambda s: [v for v in s.split("|") if v])
    d["notes"] = d[COLS["notes"]]
    return d


def bh(p):
    """Benjamini-Hochberg adjusted p-values."""
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def kappa(a, b):
    from sklearn.metrics import cohen_kappa_score
    return cohen_kappa_score(a, b) if len(set(a) | set(b)) > 1 else np.nan


def has(d, q, v):
    return d[q].apply(lambda xs: v in xs)


def fisher(a, b):
    from scipy.stats import fisher_exact
    return fisher_exact([[a.sum(), len(a) - a.sum()], [b.sum(), len(b) - b.sum()]])[1]


def analyze(d):
    from scipy.stats import mannwhitneyu
    from sklearn.metrics import roc_auc_score
    Q = schema()
    res = {}
    rows = []
    for q in ("q1", "q2", "q3", "q8"):
        dq = d[d[q].str.len() > 0]  # cytology exists only for responses given after it was added
        for o in Q[q]["options"]:
            v = o["value"]
            hit = has(dq, q, v)
            hi, lo = hit[dq.kind == "high"], hit[dq.kind == "low"]
            if hi.sum() + lo.sum() == 0:
                continue
            rows.append(dict(question=q, value=v, label=o["label"], high=hi.mean(), low=lo.mean(), n_high=len(hi),
                             n_low=len(lo), p=fisher(hi, lo)))
    cat = pd.DataFrame(rows)
    cat["q_bh"] = bh(cat.p.values)
    res["categories"] = cat
    d["msi_call"] = d.q4.str[0]
    d["cell"] = np.where(d.kind == "low", "Control regions", "High attention: " + d.group.map(GROUP_LABEL))
    res["calls"] = pd.crosstab(d.cell, d.msi_call, normalize="index").reindex(columns=["yes", "unsure", "no"], fill_value=0)
    res["calls_n"] = d.groupby("cell").size()
    w = d.pivot_table(index=["code", "region_id"], columns="username", values=["q1", "q2", "q4"], aggfunc="first")
    ag = {}
    if all(r in w["q1"].columns for r in RATERS):
        w = w.dropna()
        a, b = RATERS
        ag["Main tissue"] = kappa(w["q1"][a].str[0], w["q1"][b].str[0])
        ag["MSI-associated call"] = kappa(w["q4"][a].str[0], w["q4"][b].str[0])
        ks = [kappa(w["q2"][a].apply(lambda x: v in x), w["q2"][b].apply(lambda x: v in x))
              for v in [o["value"] for o in Q["q2"]["options"]]]
        ag["Inflammation (mean over categories)"] = np.nanmean(ks)
        ag["Neutrophils"] = kappa(w["q2"][a].apply(lambda x: "neutrophils" in x), w["q2"][b].apply(lambda x: "neutrophils" in x))
        ag["IEL"] = kappa(w["q2"][a].apply(lambda x: "iel" in x), w["q2"][b].apply(lambda x: "iel" in x))
        ag["n regions rated by both"] = len(w)
    res["agreement"] = ag
    # e) tumor epithelium only
    te = d[d.q1.str[0] == "tumor_epithelium"]
    rows = []
    for q, v in (("q2", "iel"), ("q2", "neutrophils"), ("q2", "stromal_lymphocytes"), ("q2", "plasma_cells"),
                 ("q8", "large_cells"), ("q8", "vesicular_chromatin"), ("q8", "prominent_nucleoli")):
        t = te[te[q].str.len() > 0]
        hi, lo = has(t[t.kind == "high"], q, v), has(t[t.kind == "low"], q, v)
        rows.append(dict(question=q, value=v, high=hi.mean(), low=lo.mean(), n_high=len(hi), n_low=len(lo), p=fisher(hi, lo)))
    res["tumor_epi"] = pd.DataFrame(rows)
    # f) high-attention regions by patient group, response level and slide level
    H = d[d.kind == "high"]
    feats = [("q2", "neutrophils"), ("q2", "iel"), ("q2", "stromal_lymphocytes"), ("q2", "none"), ("q3", "mucinous"),
             ("q3", "medullary"), ("q3", "poorly_diff_solid"), ("q3", "cribriform_comedo"), ("q8", "prominent_nucleoli"),
             ("q8", "vesicular_chromatin"), ("q8", "large_cells")]
    rows, slides = [], []
    for q, v in feats:
        t = H[H[q].str.len() > 0] if q == "q8" else H
        r = dict(question=q, value=v)
        for g, lab in GROUP_LABEL.items():
            x = has(t[t.group == g], q, v)
            r[lab] = x.mean()
            r["n " + lab] = len(x)
        ctl = d[d.kind == "low"]
        r["control"] = has(ctl[ctl[q].str.len() > 0] if q == "q8" else ctl, q, v).mean()
        r["p response (MSI vs all MSS)"] = fisher(has(t[t.group == "MSI_correct"], q, v), has(t[t.group != "MSI_correct"], q, v))
        s = t.assign(hit=has(t, q, v)).groupby(["code", "patient", "group"]).hit.mean().reset_index()
        y = (s.group == "MSI_correct").astype(int)
        r["slide AUC"] = roc_auc_score(y, s.hit)
        r["slide p (Mann-Whitney)"] = mannwhitneyu(s.hit[y == 1], s.hit[y == 0]).pvalue
        rows.append(r)
        slides.append(s.rename(columns={"hit": v}).set_index(["code", "patient", "group"]))
    res["groups"] = pd.DataFrame(rows)
    s = H.assign(yes=H.msi_call == "yes").groupby(["code", "patient", "group"]).yes.mean()
    res["slides"] = pd.concat(slides + [s], axis=1).reset_index()
    y = (res["slides"].group == "MSI_correct").astype(int)
    res["slide_auc_call"] = (roc_auc_score(y, res["slides"].yes), mannwhitneyu(res["slides"].yes[y == 1], res["slides"].yes[y == 0]).pvalue)
    rows = []
    for rater in RATERS:
        x = H[H.username == rater]
        r = dict(rater=rater)
        for v in ("neutrophils", "iel"):
            s = x.assign(hit=has(x, "q2", v)).groupby(["code", "group"]).hit.mean().reset_index()
            r[v + " slide AUC"] = roc_auc_score((s.group == "MSI_correct").astype(int), s.hit)
            for g, lab in GROUP_LABEL.items():
                r[f"{v} {lab}"] = has(x[x.group == g], "q2", v).mean()
        s = x.assign(hit=x.msi_call == "yes").groupby(["code", "group"]).hit.mean().reset_index()
        r["MSI call slide AUC"] = roc_auc_score((s.group == "MSI_correct").astype(int), s.hit)
        r["responses"] = len(d[d.username == rater])
        rows.append(r)
    res["raters"] = pd.DataFrame(rows)
    res["molecular"] = pd.read_csv(os.path.join(REV, "molecular.csv")) if os.path.exists(os.path.join(REV, "molecular.csv")) else None
    res["loop"] = loop()
    res["patient_stats"] = patient_stats(d)
    res["checks"] = checks(d)
    return res


def checks(d):
    """Numbers quoted in the text that are not in the tables above."""
    from sklearn.metrics import roc_auc_score
    H = d[d.kind == "high"]
    out = {}
    k = key()
    hk = k[k.kind == "high"]
    out["high-attention regions entirely inside tumor map"] = (hk.frac_in_tumor_map >= 0.999).mean()
    out["high-attention region median area mm2"] = hk.area_mm2.median()
    nm = H[(H.group == "MSI_correct") & has(H, "q2", "neutrophils")]
    out["MSI high neutrophil responses"] = len(nm)
    out["... of which dirty necrosis not ticked"] = int((~has(nm, "q3", "dirty_necrosis")).sum())
    for g, x in d[d.kind == "low"].groupby("group"):
        out[f"control neutrophils {GROUP_LABEL[g]}"] = has(x, "q2", "neutrophils").mean()
    for g, x in H.groupby("group"):
        out[f"high dirty necrosis {GROUP_LABEL[g]}"] = has(x, "q3", "dirty_necrosis").mean()

    def pauc(x, hit):
        s = x.assign(h=hit).groupby(["code", "group"]).h.mean().reset_index()
        return roc_auc_score((s.group == "MSI_correct").astype(int), s.h)
    out["neutrophil patient AUC, excluding responses with dirty necrosis"] = pauc(H, has(H, "q2", "neutrophils") & ~has(H, "q3", "dirty_necrosis"))
    te = H[H.q1.str[0] == "tumor_epithelium"]
    out["neutrophil patient AUC, tumor-epithelium responses"] = pauc(te, has(te, "q2", "neutrophils"))
    out["neutrophil patient AUC, without POLE patient"] = pauc(H[H.patient != "TCGA-AG-3892"], has(H[H.patient != "TCGA-AG-3892"], "q2", "neutrophils"))
    f = os.path.join(ROOT, "results/metrics/attention_enrichment.csv")
    if os.path.exists(f):
        e = pd.read_csv(f)
        e = e[e.patient.isin(d.patient.unique())]
        for c in ("iel", "neutrophil"):
            out[f"KongNet {c} density AUC on the {len(e)} review patients"] = roc_auc_score(e.y, e["all_" + c])
    return out


def loop():
    """Attention-gated KongNet counts and controls (src/gate_controls.py -> results/metrics/gate_controls.csv)."""
    f = os.path.join(ROOT, "results/metrics/gate_controls.csv")
    if not os.path.exists(f):
        return None
    o = pd.read_csv(f)
    fm = pd.read_csv(os.path.join(ROOT, "results/metrics/iel_fm.csv")).set_index(["cohort", "model"])
    o["fm_abmil"] = [fm.loc[(c, "fm_abmil"), "auc"] for c in o.cohort]
    return o


def patient_stats(d, n_boot=5000, seed=42):
    """Per-patient share of high-attention responses naming each category: AUC MSI vs MSS with stratified bootstrap
    CI, Mann-Whitney p and Benjamini-Hochberg q across categories."""
    from scipy.stats import mannwhitneyu
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    H = d[d.kind == "high"]
    feats = [("q2", v) for v in ("neutrophils", "iel", "stromal_lymphocytes", "plasma_cells", "eosinophils", "none")] + \
            [("q3", v) for v in ("mucinous", "medullary", "poorly_diff_solid", "cribriform_comedo", "dirty_necrosis", "glandular")] + \
            [("q8", v) for v in ("prominent_nucleoli", "vesicular_chromatin", "large_cells", "eosinophilic_amphophilic")]
    rows = []
    for q, v in feats:
        t = H[H[q].str.len() > 0] if q == "q8" else H
        s = t.assign(h=has(t, q, v)).groupby(["code", "group"]).h.mean().reset_index()
        y = (s.group == "MSI_correct").values.astype(int)
        x = s.h.values
        pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
        b = [roc_auc_score(y[i], x[i]) for i in (np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))] for _ in range(n_boot))]
        rows.append(dict(question=q, value=v, auc=roc_auc_score(y, x), lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5),
                         p=mannwhitneyu(x[y == 1], x[y == 0]).pvalue))
    r = pd.DataFrame(rows)
    r["q_bh"] = bh(r.p.values)
    return r


def figure(res, simulated):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import figs_extra as FX
    C_HI, C_LO = "#b8336a", "#8f8d87"
    C_G = {"MSI_correct": "#b8336a", "MSS_false_positive": "#e39a3b", "MSS_correct": "#4b6a88"}
    fig = plt.figure(figsize=(12, 9.2))
    # a: categories, high vs control
    cat = res["categories"].sort_values(["question", "high"], ascending=[True, False])
    ax = fig.add_axes([0.2, 0.45, 0.31, 0.53])
    y = np.arange(len(cat))
    ax.barh(y - 0.2, cat.high, 0.38, color=C_HI, label="High-attention regions")
    ax.barh(y + 0.2, cat.low, 0.38, color=C_LO, label="Control regions")
    for yi, (q_, l_) in enumerate(zip(cat.q_bh.values, cat[["high", "low"]].max(axis=1).values)):
        if q_ < 0.05:
            ax.text(l_ + 0.01, yi, "*" if q_ >= 0.01 else "**", va="center", fontsize=10)
    short = {"Tumor-stroma interface / invasive front": "Interface / invasive front",
             "Moderate eosinophilic or amphophilic cytoplasm": "Eosinophilic / amphophilic cytoplasm",
             "Not assessable (no tumor in region)": "Not assessable (no tumor)",
             "Neutrophils (epithelium or lumen)": "Neutrophils (epithelium/lumen)"}
    ax.set_yticks(y, [short.get(l, l) for l in cat.label], fontsize=7)
    ax.invert_yaxis()
    start = 0
    for qq, name in (("q1", "main tissue"), ("q2", "inflammation"), ("q3", "morphology"), ("q8", "cytology*")):
        n_ = (cat.question == qq).sum()
        if n_ == 0:
            continue
        if start:
            ax.axhline(start - 0.5, color=FX.INK2, lw=0.8)
        ax.text(1.0, start + n_ / 2 - 0.5, name, rotation=-90, va="center", ha="left", fontsize=8, color=FX.INK2,
                transform=ax.get_yaxis_transform())
        start += n_
    ax.set_xlim(0, 1)
    ax.set_xlabel("Share of responses (both pathologists)")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title("High-attention vs control regions (* q < 0.05, ** q < 0.01)", fontsize=9)
    fig.text(0.01, 0.985, "a", fontsize=14, fontweight="bold", va="top")
    # b: per-slide share naming neutrophils / IEL in high-attention regions, by patient group
    ax = fig.add_axes([0.62, 0.56, 0.36, 0.38])
    sl = res["slides"]
    mol = res["molecular"].set_index("patient") if res["molecular"] is not None else None
    rng = np.random.default_rng(0)
    gr = res["groups"].set_index("value")
    for i, (v, lab) in enumerate((("neutrophils", "Neutrophils"), ("iel", "Intraepithelial\nlymphocytes"),
                                  ("prominent_nucleoli", "Prominent\nnucleoli*"))):
        for j, g in enumerate(GROUP_LABEL):
            s = sl[sl.group == g]
            x = i * 4 + j + rng.uniform(-0.18, 0.18, len(s))
            ax.scatter(x, s[v], s=26, color=C_G[g], edgecolor="white", lw=0.5, zorder=3,
                       label=f"{GROUP_LABEL[g]} ({len(s)} patients)" if i == 0 else None)
            ax.plot([i * 4 + j - 0.3, i * 4 + j + 0.3], [s[v].median()] * 2, color=FX.INK, lw=1.4, zorder=4)
            if mol is not None:
                for xx, (_, r) in zip(x, s.iterrows()):
                    pole = str(mol.loc[r.patient, "pole"]) if r.patient in mol.index else ""
                    if "S459F" in pole or "P286R" in pole or "V411L" in pole:
                        ax.annotate("POLE", (xx, r[v]), xytext=(6, -2), textcoords="offset points", fontsize=7,
                                    color=C_G[g], fontweight="bold")
        ax.plot([i * 4 - 0.4, i * 4 + 2.4], [gr.loc[v, "control"]] * 2, color=C_LO, lw=1, ls="--", zorder=2)
        ax.text(i * 4 + 1, 1.04, f"patient AUC {gr.loc[v, 'slide AUC']:.2f}", ha="center", fontsize=7.5, color=FX.INK)
    ax.set_xticks([1, 5, 9], ["Neutrophils", "Intraepithelial\nlymphocytes", "Prominent\nnucleoli*"], fontsize=8)
    ax.set_ylim(-0.03, 1.1)
    ax.set_xlim(-0.8, 10.8)
    ax.set_ylabel("Share of responses per patient\n(high-attention regions)")
    ax.grid(axis="x", visible=False)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False)
    ax.set_title("High-attention regions by patient (dashed: control regions)", fontsize=9)
    fig.text(0.56, 0.985, "b", fontsize=14, fontweight="bold", va="top")
    # c: MSI-associated calls
    ax = fig.add_axes([0.13, 0.08, 0.26, 0.26])
    order = ["High attention: MSI", "High attention: MSS, called MSI", "High attention: MSS", "Control regions"]
    cl = res["calls"].reindex(order).fillna(0)
    left = np.zeros(len(cl))
    for col, c in (("yes", C_HI), ("unsure", "#e8c9d6"), ("no", "#d9d8d4")):
        ax.barh(range(len(cl)), cl[col], left=left, color=c, label=col, edgecolor="white")
        left += cl[col].values
    ax.set_yticks(range(len(cl)), [f"{o.replace('High attention: ', 'High attn. ').replace('Control regions', 'Control')} (n={res['calls_n'].get(o, 0)})"
                                   for o in order], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Share of responses")
    ax.grid(axis="y", visible=False)
    ax.legend(title="Looks MSI-associated?", ncol=3, fontsize=7.5, title_fontsize=8, loc="upper center",
              bbox_to_anchor=(0.45, -0.22), frameon=False)
    ag = res["agreement"]
    ax.set_title(f"Overall impression (patient AUC {res['slide_auc_call'][0]:.2f}; κ = {ag.get('MSI-associated call', np.nan):.2f})",
                 fontsize=9)
    fig.text(0.01, 0.4, "c", fontsize=14, fontweight="bold", va="top")
    # d: attention-gated counts and controls
    lp = res["loop"]
    if lp is not None:
        cols = {"tcga_test": FX.C_TCGA, "cptac": FX.C_CPTAC, "surgen": FX.C_SURGEN}
        names = {"tcga_test": "TCGA test", "cptac": "CPTAC", "surgen": "SurGen"}
        feats = (("neutrophil", "Neutrophils"), ("iel", "IEL"), ("lymphocyte", "Lymphocytes"), ("plasma", "Plasma cells"),
                 ("eosinophil", "Eosinophils"), ("connective", "Connective"), ("epithelial", "Epithelial"))
        marks = (("all", dict(marker="o", facecolor="white"), "All tumor tiles"),
                 ("epi5", dict(marker="^", facecolor="white"), "Most epithelial 5%"),
                 ("perm5", dict(marker="x"), "Permuted-label attention 5%"),
                 ("att5", dict(marker="o"), "Attention top 5%"))
        for j, co in enumerate(cols):
            ax = fig.add_axes([0.63 + j * 0.12, 0.1, 0.105, 0.26])
            g = lp[lp.cohort == co].set_index(["feature", "gate"]).auc
            for i, (f_, lab) in enumerate(feats):
                ax.plot([g[(f_, "all")], g[(f_, "att5")]], [i, i], color=cols[co], lw=1, zorder=2)
                for gate, kw, glab in marks:
                    kw = dict(kw)
                    fc = kw.pop("facecolor", cols[co])
                    ax.scatter([g[(f_, gate)]], [i], s=22, color=cols[co] if kw["marker"] == "x" else None,
                               facecolor=fc if kw["marker"] != "x" else None, edgecolor=cols[co], lw=1, zorder=3, **kw)
            ax.axvline(lp[lp.cohort == co].fm_abmil.iloc[0], color=cols[co], lw=0.8, ls=":")
            ax.axvline(0.5, color=FX.GRID, lw=0.8)
            ax.set_xlim(0.25, 1.0)
            ax.set_xticks([0.3, 0.5, 0.7, 0.9])
            ax.tick_params(labelsize=7)
            ax.set_yticks(range(len(feats)), [l for _, l in feats] if j == 0 else [], fontsize=7.5)
            ax.invert_yaxis()
            ax.grid(axis="y", visible=False)
            ax.set_title(names[co], fontsize=8.5, color=cols[co])
            if j == 1:
                ax.set_xlabel("ROC-AUC of mean cell count per tile (dotted: FM ABMIL)", fontsize=8)
            if j == 0:
                hs = [ax.scatter([], [], s=22, marker=kw["marker"], color=FX.INK2 if kw["marker"] == "x" else None,
                                 facecolor=("white" if kw.get("facecolor") == "white" else FX.INK2) if kw["marker"] != "x" else None,
                                 edgecolor=FX.INK2, lw=1) for _, kw, _ in marks]
                ax.legend(hs, [l for _, _, l in marks], fontsize=6.5, loc="upper left", bbox_to_anchor=(-0.05, -0.2),
                          ncol=4, frameon=False, handletextpad=0.2, columnspacing=0.8)
        fig.text(0.56, 0.4, "d", fontsize=14, fontweight="bold", va="top")
        fig.text(0.75, 0.395, "Attention-gated counts are not specific", fontsize=9, ha="center", va="top")
    if simulated:
        fig.text(0.5, 0.5, "SIMULATED DATA\nplaceholder layout", fontsize=60, color="red", alpha=0.18, ha="center",
                 va="center", rotation=25, fontweight="bold")
    fig.savefig(os.path.join(ROOT, "results/fig10_review.png"), dpi=200)
    fig.savefig(os.path.join(ROOT, "results/fig10_review.pdf"), dpi=300)


def cmd_analyze(a):
    os.makedirs(OUT, exist_ok=True)
    d = load(a.export)
    res = analyze(d)
    d.drop(columns=["q1", "q2", "q3", "q4", "q8"]).to_csv(os.path.join(OUT, "responses_unblinded.csv"), index=False)
    md = lambda t: t.round(3).to_markdown(index=False)  # noqa: E731
    parts = ["# SIMULATED DATA (placeholder)" if a.simulated else "# Pathologist review",
             f"{d.groupby(['code', 'region_id']).ngroups} regions, {len(d)} responses "
             f"({', '.join(f'{r}: {n}' for r, n in d.username.value_counts().items())}); "
             f"free-text notes: {(d.notes.str.strip() != '').sum()}",
             "## a) High-attention vs control regions\n\n" + md(res["categories"]),
             "## b) Looks MSI-associated?\n\n" + res["calls"].round(3).to_markdown() + "\n\n" + res["calls_n"].to_markdown()
             + f"\n\nslide AUC of the share of 'yes' (MSI vs MSS patients): {res['slide_auc_call'][0]:.3f} "
               f"(Mann-Whitney p = {res['slide_auc_call'][1]:.2g})",
             "## d) Agreement (Cohen's kappa)\n\n" + "\n".join(f"- {k}: {v:.3f}" if isinstance(v, float) else f"- {k}: {v}"
                                                            for k, v in res["agreement"].items()),
             "## e) Tumor epithelium only (exploratory)\n\n" + md(res["tumor_epi"]),
             "## f) High-attention regions by patient group (exploratory)\n\n" + md(res["groups"])
             + "\n\nPer rater:\n\n" + md(res["raters"]) + "\n\nPer slide:\n\n" + md(res["slides"].drop(columns="code"))]
    parts.append("## Checks (numbers quoted in the text)\n\n" + "\n".join(
        f"- {k_}: {v:.3f}" if isinstance(v, float) else f"- {k_}: {v}" for k_, v in res["checks"].items()))
    if res["molecular"] is not None:
        parts.append("## h) Molecular annotation of review patients\n\n" + md(res["molecular"]))
    parts.append("## f2) Patient-level AUC with bootstrap CI and BH q (high-attention regions)\n\n" + md(res["patient_stats"]))
    if res["loop"] is not None:
        parts.append("## g) Attention-gated counts and controls (see Supplementary Table 5)\n\n" + md(res["loop"]))
    with open(os.path.join(OUT, "summary.md"), "w") as fh:
        fh.write("\n\n".join(parts) + "\n")
    figure(res, a.simulated)
    print(open(os.path.join(OUT, "summary.md")).read())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("simulate")
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(fn=cmd_simulate)
    s = sp.add_parser("from-deposit")
    s.add_argument("dir")
    s.set_defaults(fn=cmd_from_deposit)
    s = sp.add_parser("analyze")
    s.add_argument("export")
    s.add_argument("--simulated", action="store_true")
    s.set_defaults(fn=cmd_analyze)
    a = ap.parse_args()
    a.fn(a)
