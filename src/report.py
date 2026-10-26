"""Collect metrics into results/summary.csv, paired E2 bootstrap comparisons, a markdown table, and Figures 1-4.

  python src/report.py
"""
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import protocols as P  # noqa: E402

R = os.path.join(P.ROOT, "results")
FM_ORDER = ["phikon_ref", "phikon", "phikon_v2", "uni", "uni2_h", "virchow2", "gigapath", "flash"]
FEAT_ORDER = ["feat_tissue", "feat_cells", "feat_interp"]
ENC_ORDER = FM_ORDER + FEAT_ORDER
ENC_LABEL = {"phikon_ref": "Phikon (archived)", "phikon": "Phikon", "phikon_v2": "Phikon-v2", "uni": "UNI",
             "uni2_h": "UNI2-h", "virchow2": "Virchow2", "gigapath": "Prov-GigaPath", "flash": "Prov-GigaPath-flash",
             "feat_tissue": "Tissue fractions", "feat_cells": "Nucleus densities", "feat_interp": "Tissue + nuclei"}
# image-source arms for the resolution comparison (extract.py SOURCES)
ARMS = [("", "Zenodo tiles"), ("_rt114", "Re-cut 1.14 µm/px"), ("_rt05", "Re-cut 0.5 µm/px")]
RT_ENCS = [e + suf for e in FM_ORDER[1:] for suf, _ in ARMS[1:]]
AGG_ORDER = ["mean", "max", "mean+max", "abmil"]
AGG_LABEL = {"mean": "Mean", "max": "Max", "mean+max": "Mean+max", "abmil": "ABMIL"}
# reference categorical palette, slots 1-4 in fixed order; marker shape is the secondary encoding
COLOR = {"mean": "#2a78d6", "max": "#eb6834", "mean+max": "#1baf7a", "abmil": "#eda100"}
MARKER = {"mean": "o", "max": "s", "mean+max": "D", "abmil": "^"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "legend.frameon": False, "figure.dpi": 150, "savefig.bbox": "tight"})


def proto_name(r):
    if r.scheme == "all":
        return {"E1": "E1_random_cv_seed42", "E1rep": "E1_random_cv_10x5", "E1paper": "E1_paper_folds",
                "E2": "E2_source_split"}[r.protocol]
    if r.scheme == "fixedN":
        return f"E3a_fixedN{int(r.N)}_{r.protocol}"
    if r.scheme == "allN":
        return f"E3a_allpatch_eligibleN{int(r.N)}_{r.protocol}"
    return f"E3b_count_matched_{r.protocol}"


def encoders_present():
    found = {os.path.basename(f).rsplit("_", 1)[0] for f in glob.glob(f"{R}/metrics/*_lr.csv")}
    # the CRC report covers the six foundation models; the fast Edge encoder is reported there only via the
    # whole-slide/external analyses. Other roots (e.g. STAD) include it.
    if "MSI_ROOT" not in os.environ:
        found.discard("flash")
    return [e for e in ENC_ORDER if e in found]


def load_metrics():
    fs = glob.glob(f"{R}/metrics/*_lr.csv") + glob.glob(f"{R}/metrics/*_abmil.csv")
    m = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    m["protocol_full"] = m.apply(proto_name, axis=1)
    # count-only is encoder independent: keep one copy
    cnt = m[m.aggregation == "count"]
    cnt = cnt[cnt.encoder == cnt.encoder.iloc[0]].assign(encoder="count-only") if len(cnt) else cnt
    return pd.concat([m[m.aggregation != "count"], cnt], ignore_index=True)


def paired_e2(encs):
    """Paired stratified bootstrap on identical TEST resamples: AUC differences between aggregations (all-patch
    training) and between count-matched and all-patch training (same aggregation). With several seeds or matching
    draws, each side is the mean AUC over units, on every resample."""
    D = P.Design()
    preds = {}
    for e in encs:
        for learner in ("lr", "abmil"):
            f = f"{R}/preds/{e}_{learner}.parquet"
            if not os.path.exists(f):
                continue
            pr = pd.read_parquet(f, filters=[("protocol", "==", "E2")])
            for (agg, sch), g in pr[pr.scheme.isin(["all", "matched"])].groupby(["aggregation", "scheme"]):
                preds[(e, agg, sch)] = g
    pairs = [(("max", "all"), ("mean", "all")), (("mean+max", "all"), ("mean", "all")),
             (("abmil", "all"), ("mean", "all")), (("max", "all"), ("abmil", "all")),
             (("mean+count", "all"), ("mean", "all")), (("max+count", "all"), ("max", "all")),
             (("mean+max+count", "all"), ("mean+max", "all"))]
    pairs += [((a, "matched"), (a, "all")) for a in AGG_ORDER]
    # interpretable features (mean pooling) vs each foundation model's mean pooling and ABMIL
    cross = [(("feat_interp", "mean", "all"), (fm, agg, "all")) for fm in encs if fm in FM_ORDER
             for agg in ("mean", "abmil")]
    within = [((e, *a), (e, *b)) for e in encs for a, b in pairs]
    rows = []
    for ka, kb in within + cross:
        if ka not in preds or kb not in preds:
            continue
        e = ka[0] if ka[0] == kb[0] else f"{ka[0]} vs {kb[0]}"
        a, b = ka[1:], kb[1:]
        ga, gb = preds[ka], preds[kb]
        ba, bb = P.boot_unit_aucs(ga, D.y)[1], P.boot_unit_aucs(gb, D.y)[1]
        va, vb = (np.mean([P.metrics(D.y[u.pid.values], u.p.values)["roc_auc"]
                           for _, u in g.groupby(["repeat", "seed"])]) for g in (ga, gb))
        da = ba - bb
        p = 2 * min((da <= 0).mean(), (da >= 0).mean())
        rows.append(dict(encoder=e, comparison=f"{a[0]} [{a[1]}] - {b[0]} [{b[1]}]", auc_a=va, auc_b=vb,
                         delta=va - vb, CI_low=np.percentile(da, 2.5), CI_high=np.percentile(da, 97.5),
                         p_boot=min(1.0, max(p, 1 / len(da)))))
    return pd.DataFrame(rows)


def build_summary(m, encs):
    auc = m[m.metric.isin(["roc_auc", "pr_auc", "brier"])]
    s = auc[["encoder", "aggregation", "protocol_full", "metric", "value", "CI_low", "CI_high", "sd", "n_units",
             "n_patients"]].rename(columns={"protocol_full": "protocol"})
    extra = []
    for e in encs:
        f = f"{R}/metrics/{e}_confound.csv"
        if os.path.exists(f):
            extra.append(pd.read_csv(f))
    f = f"{R}/metrics/count_stats.csv"
    if os.path.exists(f):
        extra.append(pd.read_csv(f))
    for e in RT_ENCS:
        f = f"{R}/metrics/{e}_confound.csv"
        if os.path.exists(f):
            extra.append(pd.read_csv(f))
    s = pd.concat([s] + extra, ignore_index=True)
    s["encoder"] = pd.Categorical(s.encoder, ENC_ORDER + RT_ENCS + ["count-only"], ordered=True)
    s = s.sort_values(["encoder", "aggregation", "protocol", "metric"])
    cols = ["encoder", "aggregation", "protocol", "metric", "value", "CI_low", "CI_high", "sd", "n_units", "n_patients"]
    s[cols].to_csv(f"{R}/summary.csv", index=False, float_format="%.4f")
    return s


def get(m, e, agg, proto, metric="roc_auc", col="value"):
    r = m[(m.encoder == e) & (m.aggregation == agg) & (m.protocol_full == proto) & (m.metric == metric)]
    return r[col].iloc[0] if len(r) else np.nan


def table(m, encs, conf):
    lines = ["| Encoder | Aggregation | E1 random CV 10x5 | E2 source split [95% CI] | Inflation (E1-E2) | "
             "E3b matched E1 | E3b matched E2 | Fixed-N 100 E1 | Fixed-N 100 E2 |",
             "|---|---|---|---|---|---|---|---|---|"]
    for e in encs:
        for agg in AGG_ORDER:
            e1, sd = get(m, e, agg, "E1_random_cv_10x5"), get(m, e, agg, "E1_random_cv_10x5", col="sd")
            e2, lo, hi = (get(m, e, agg, "E2_source_split", col=c) for c in ("value", "CI_low", "CI_high"))
            if np.isnan(e1):
                continue
            lines.append(f"| {ENC_LABEL[e]} | {AGG_LABEL[agg]} | {e1:.3f} ± {sd:.3f} | {e2:.3f} [{lo:.3f}, {hi:.3f}] | "
                         f"{e1 - e2:+.3f} | {get(m, e, agg, 'E3b_count_matched_E1'):.3f} | "
                         f"{get(m, e, agg, 'E3b_count_matched_E2'):.3f} | {get(m, e, agg, 'E3a_fixedN100_E1'):.3f} | "
                         f"{get(m, e, agg, 'E3a_fixedN100_E2'):.3f} |")
    c = "count-only"
    lines.append(f"| log1p(count) only | LR | {get(m, c, 'count', 'E1_random_cv_10x5'):.3f} ± "
                 f"{get(m, c, 'count', 'E1_random_cv_10x5', col='sd'):.3f} | {get(m, c, 'count', 'E2_source_split'):.3f} "
                 f"[{get(m, c, 'count', 'E2_source_split', col='CI_low'):.3f}, "
                 f"{get(m, c, 'count', 'E2_source_split', col='CI_high'):.3f}] | | "
                 f"{get(m, c, 'count', 'E3b_count_matched_E1'):.3f} | {get(m, c, 'count', 'E3b_count_matched_E2'):.3f} | | |")
    out = ["# Main results (ROC-AUC)", "", *lines, ""]
    if len(conf):
        out += ["# Count decoding (out-of-fold ridge, log1p(count))", "",
                "| Encoder | Cohort | Mean R² | Max R² | Mean+max R² | Mean ρ | Max ρ |", "|---|---|---|---|---|---|---|"]
        for e in encs:
            for coh in ("all", "MSS_only"):
                g = conf[(conf.encoder == e) & (conf.protocol == f"E4_decode_{coh}")]
                if not len(g):
                    continue
                v = {(r.aggregation, r.metric): r.value for r in g.itertuples()}
                out.append(f"| {ENC_LABEL[e]} | {coh} | {v[('mean', 'count_r2')]:.3f} | {v[('max', 'count_r2')]:.3f} | "
                           f"{v[('mean+max', 'count_r2')]:.3f} | {v[('mean', 'count_spearman')]:.3f} | "
                           f"{v[('max', 'count_spearman')]:.3f} |")
    open(f"{R}/tables.md", "w").write("\n".join(out) + "\n")


def feature_table():
    """Interpretable features: per-patient means (all patches), univariate ROC-AUC within TRAIN and TEST, and
    standardized coefficients of the E2 LR fit on TRAIN (same pipeline as eval_lr)."""
    import json

    from sklearn.metrics import roc_auc_score

    import eval_lr
    import pool
    f = os.path.join(P.ROOT, "emb/feat_interp/meta.json")
    if not os.path.exists(f):
        return ""
    names = json.load(open(f))["features"]
    X = pool.pooled_all("feat_interp")[0]
    D = P.Design()
    model = eval_lr.lr().fit(X[D.train], D.y[D.train])
    coef = model[-1].coef_[0]
    lines = ["# Interpretable features (patient mean over all patches)", "",
             "| Feature | Mean MSIMUT (TRAIN) | Mean MSS (TRAIN) | AUC TRAIN | AUC TEST | E2 LR coef (std.) |",
             "|---|---|---|---|---|---|"]
    rows = []
    for j, n in enumerate(names):
        tr, te = D.train, D.test
        rows.append(dict(feature=n, mean_msimut_train=X[tr][D.y[tr] == 1, j].mean(),
                         mean_mss_train=X[tr][D.y[tr] == 0, j].mean(), auc_train=roc_auc_score(D.y[tr], X[tr, j]),
                         auc_test=roc_auc_score(D.y[te], X[te, j]), coef=coef[j]))
        r = rows[-1]
        lines.append(f"| {n} | {r['mean_msimut_train']:.3f} | {r['mean_mss_train']:.3f} | {r['auc_train']:.3f} | "
                     f"{r['auc_test']:.3f} | {r['coef']:+.2f} |")
    pd.DataFrame(rows).to_csv(f"{R}/feature_univariate.csv", index=False, float_format="%.4f")
    return "\n".join(lines) + "\n"


def resolution(m):
    """Zenodo tiles vs the same regions re-cut from the slides at 1.14 and 0.5 um/px. Table, paired E2 bootstrap,
    Figure 5. Returns markdown (empty if no re-cut results yet)."""
    have = [e for e in FM_ORDER[1:] if any(os.path.exists(f"{R}/metrics/{e}{suf}_lr.csv") for suf, _ in ARMS[1:])
            and (e != "flash" or "MSI_ROOT" in os.environ)]
    if not have:
        return ""
    D = P.Design()
    lines = ["# Resolution: Zenodo tiles vs re-cut from the slides (ROC-AUC)", "",
             "Zenodo tiles: as published (256 µm at 224 px = 1.14 µm/px, Macenko). Re-cut: the same 256 µm regions read "
             "from the TCGA slides (raw colour) at 1.14 µm/px (one 224 px view) or 0.5 µm/px (mean of four 224 px crops).",
             "", "| Encoder | Aggregation | " + " | ".join(f"{lab} E1 | {lab} E2" for _, lab in ARMS) + " |",
             "|---|---|" + "---|---|" * len(ARMS)]
    for e in have:
        for agg in AGG_ORDER:
            cells = []
            for suf, _ in ARMS:
                e1, e2 = get(m, e + suf, agg, "E1_random_cv_10x5"), get(m, e + suf, agg, "E2_source_split")
                cells += [f"{e1:.3f}" if not np.isnan(e1) else "", f"{e2:.3f}" if not np.isnan(e2) else ""]
            lines.append(f"| {ENC_LABEL[e]} | {AGG_LABEL[agg]} | " + " | ".join(cells) + " |")
    rows = []
    for e in have:
        preds = {}
        for suf, _ in ARMS:
            for learner in ("lr", "abmil"):
                f = f"{R}/preds/{e}{suf}_{learner}.parquet"
                if os.path.exists(f):
                    pr = pd.read_parquet(f, filters=[("protocol", "==", "E2"), ("scheme", "==", "all")])
                    for agg, g in pr.groupby("aggregation"):
                        preds[(suf, agg)] = g
        for agg in AGG_ORDER:
            for a, b in (("_rt05", ""), ("_rt114", ""), ("_rt05", "_rt114")):
                if (a, agg) not in preds or (b, agg) not in preds:
                    continue
                ba = P.boot_unit_aucs(preds[(a, agg)], D.y)[1]
                bb = P.boot_unit_aucs(preds[(b, agg)], D.y)[1]
                va, vb = (np.mean([P.metrics(D.y[u.pid.values], u.p.values)["roc_auc"]
                                   for _, u in preds[(k, agg)].groupby(["repeat", "seed"])]) for k in (a, b))
                da = ba - bb
                p = 2 * min((da <= 0).mean(), (da >= 0).mean())
                name = {"": "zenodo", "_rt114": "recut_1.14", "_rt05": "recut_0.5"}
                rows.append(dict(encoder=e, aggregation=agg, comparison=f"{name[a]} - {name[b]}", auc_a=va, auc_b=vb,
                                 delta=va - vb, CI_low=np.percentile(da, 2.5), CI_high=np.percentile(da, 97.5),
                                 p_boot=min(1.0, max(p, 1 / len(da)))))
    pr = pd.DataFrame(rows)
    pr.to_csv(f"{R}/paired_e2_resolution.csv", index=False, float_format="%.4f")
    if len(pr):
        lines += ["", "Paired E2 differences (same TEST bootstrap resamples):", "",
                  "| Encoder | Aggregation | Comparison | Δ AUC [95% CI] |", "|---|---|---|---|"]
        for r in pr.itertuples():
            lines.append(f"| {ENC_LABEL[r.encoder]} | {AGG_LABEL[r.aggregation]} | {r.comparison} | "
                         f"{r.delta:+.3f} [{r.CI_low:+.3f}, {r.CI_high:+.3f}] |")
    fig5(m, have)
    return "\n".join(lines) + "\n"


def fig5(m, encs):
    """E2 (source split) AUC per encoder for each image source; mean pooling and ABMIL, plus max-pool E1 inflation."""
    arm_color = {"": "#9a9893", "_rt114": "#2a78d6", "_rt05": "#1baf7a"}
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4), sharey=False)
    for ax, (agg, proto, title) in zip(axs, (("mean", "E2_source_split", "Mean pooling, source split"),
                                            ("abmil", "E2_source_split", "ABMIL, source split"),
                                            ("max", "E1_random_cv_10x5", "Max pooling, random CV (inflated)"))):
        for j, (suf, lab) in enumerate(ARMS):
            x = np.arange(len(encs)) + (j - 1) * 0.25
            v = np.array([get(m, e + suf, agg, proto) for e in encs])
            if proto.startswith("E2"):
                lo = np.array([get(m, e + suf, agg, proto, col="CI_low") for e in encs])
                hi = np.array([get(m, e + suf, agg, proto, col="CI_high") for e in encs])
                ax.errorbar(x, v, yerr=[v - lo, hi - v], fmt="none", ecolor=arm_color[suf], elinewidth=0.8, alpha=0.5)
            ax.scatter(x, v, s=30, color=arm_color[suf], edgecolor="white", lw=0.6, zorder=3, label=lab)
        ax.set_xticks(range(len(encs)), [ENC_LABEL[e] for e in encs], rotation=30, ha="right")
        ax.set_title(title, fontsize=9)
        ax.grid(axis="x", visible=False)
    axs[0].set_ylabel("ROC-AUC")
    fig.legend(*axs[0].get_legend_handles_labels(), ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.06))
    fig.savefig(f"{R}/fig5_resolution.png")
    fig.savefig(f"{R}/fig5_resolution.pdf")
    plt.close(fig)


def fig1(m, encs, name="fig1_random_vs_source"):
    fig, ax = plt.subplots(figsize=(1.35 * len(encs) + 1.5, 3.6))
    off = {a: (i - 1.5) * 0.18 for i, a in enumerate(AGG_ORDER)}
    for i, e in enumerate(encs):
        for agg in AGG_ORDER:
            e1 = get(m, e, agg, "E1_random_cv_10x5")
            e2, lo, hi = (get(m, e, agg, "E2_source_split", col=c) for c in ("value", "CI_low", "CI_high"))
            if np.isnan(e1):
                continue
            x = i + off[agg]
            ax.plot([x, x], [e2, e1], color=COLOR[agg], lw=1.2, alpha=0.6, zorder=2)
            ax.errorbar(x, e2, yerr=[[e2 - lo], [hi - e2]], fmt="none", ecolor=COLOR[agg], elinewidth=0.8,
                        alpha=0.5, zorder=1)
            ax.scatter(x, e1, s=36, marker=MARKER[agg], color=COLOR[agg], edgecolor="white", lw=0.8, zorder=3)
            ax.scatter(x, e2, s=36, marker=MARKER[agg], facecolor="white", edgecolor=COLOR[agg], lw=1.4, zorder=3)
    c1 = get(m, "count-only", "count", "E1_random_cv_10x5")
    ax.axhline(c1, color=INK2, lw=0.8, ls="--")
    ax.text(-0.45, 0.55, f"dashed line: patch count alone, random CV ({c1:.2f})", ha="left", va="bottom",
            color=INK2, fontsize=7.5)
    ax.axhline(0.5, color=INK2, lw=0.6, ls=":")
    ax.set_xticks(range(len(encs)), [ENC_LABEL[e] for e in encs])
    ax.set_ylabel("ROC-AUC")
    ax.set_ylim(0.45, 1.0)
    ax.grid(axis="x", visible=False)
    h = [plt.Line2D([], [], color=COLOR[a], marker=MARKER[a], ls="", label=AGG_LABEL[a]) for a in AGG_ORDER]
    h += [plt.Line2D([], [], color=INK2, marker="o", ls="", label="Random CV (10x5)"),
          plt.Line2D([], [], color=INK2, marker="o", mfc="white", ls="", label="Source split, 95% CI")]
    ax.legend(handles=h, ncol=6, loc="upper center", bbox_to_anchor=(0.5, 1.13), fontsize=7.5, handletextpad=0.3)
    fig.savefig(f"{R}/{name}.png")
    fig.savefig(f"{R}/{name}.pdf")
    plt.close(fig)


def fig2(conf, encs):
    fig, axs = plt.subplots(1, 2, figsize=(1.0 * len(encs) + 3, 3.0), sharey=True)
    for ax, coh, title in zip(axs, ("all", "MSS_only"), ("All patients", "MSS patients only (no label mediation)")):
        for j, agg in enumerate(("mean", "max")):
            v = [conf[(conf.encoder == e) & (conf.aggregation == agg) & (conf.protocol == f"E4_decode_{coh}")
                      & (conf.metric == "count_r2")].value.squeeze() for e in encs]
            x = np.arange(len(encs)) + (j - 0.5) * 0.36
            ax.bar(x, v, width=0.34, color=COLOR[agg], label=AGG_LABEL[agg], edgecolor="white", lw=1)
            for xi, vi in zip(x, v):
                ax.text(xi, max(vi, 0) + 0.01, f"{vi:.2f}", ha="center", va="bottom", fontsize=6.5, color=INK2)
        ax.set_xticks(range(len(encs)), [ENC_LABEL[e] for e in encs], rotation=30, ha="right")
        ax.set_title(title, fontsize=9, color=INK)
        ax.text(-0.08, 1.06, "ab"[list(axs).index(ax)], transform=ax.transAxes, fontsize=11, fontweight="bold")
        ax.axhline(0, color=INK2, lw=0.6)
        ax.grid(axis="x", visible=False)
    axs[0].set_ylabel("Out-of-fold R², log1p(patch count)")
    axs[0].set_ylim(min(-0.1, axs[0].get_ylim()[0]), 1.05)
    fig.legend(*axs[0].get_legend_handles_labels(), ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.06))
    fig.savefig(f"{R}/fig2_count_decoding.png")
    fig.savefig(f"{R}/fig2_count_decoding.pdf")
    plt.close(fig)


def fig3(m, encs):
    """Fixed-N subsampling. Page-width layout: two encoders per row, each as a random-CV / test-split pair."""
    Ns = [25, 50, 100, 200]
    encs = [e for e in ("phikon", "phikon_v2", "uni", "uni2_h", "virchow2", "gigapath") if e in encs]  # the six FMs
    nrow = (len(encs) + 1) // 2
    fig, axs = plt.subplots(nrow, 4, figsize=(7.2, 1.85 * nrow + 0.6), sharey=True, sharex=True, squeeze=False)
    for i, e in enumerate(encs):
        for k, proto in enumerate(("E1", "E2")):
            ax = axs[i // 2, (i % 2) * 2 + k]
            for agg in AGG_ORDER:
                v = [get(m, e, agg, f"E3a_fixedN{N}_{proto}") for N in Ns]
                sd = [get(m, e, agg, f"E3a_fixedN{N}_{proto}", col="sd") for N in Ns]
                if np.all(np.isnan(v)):
                    continue
                ax.errorbar(Ns, v, yerr=sd, color=COLOR[agg], marker=MARKER[agg], ms=3, lw=1.2, capsize=0,
                            elinewidth=0.7)
                comp = [get(m, e, agg, f"E3a_allpatch_eligibleN{N}_{proto}") for N in Ns]
                ax.plot(Ns, comp, color=COLOR[agg], lw=0.9, ls="--", alpha=0.8)
            ax.set_xscale("log")
            ax.set_xticks(Ns, [str(n) for n in Ns], fontsize=7)
            ax.minorticks_off()
            ax.tick_params(labelsize=7)
            ax.set_title(f"{ENC_LABEL[e]}: {'random CV' if proto == 'E1' else 'test split'}", fontsize=7.5)
            if (i % 2) * 2 + k == 0:
                ax.set_ylabel("ROC-AUC", fontsize=8)
            if i // 2 == nrow - 1:
                ax.set_xlabel("Tiles per patient (N)", fontsize=7.5)
    for j in range(len(encs), nrow * 2):
        for k in range(2):
            axs[j // 2, (j % 2) * 2 + k].set_visible(False)
    h = [plt.Line2D([], [], color=COLOR[a], marker=MARKER[a], label=AGG_LABEL[a]) for a in AGG_ORDER]
    h += [plt.Line2D([], [], color=INK2, ls="-", label="Fixed N (mean ± SD, 10 draws)"),
          plt.Line2D([], [], color=INK2, ls="--", label="All tiles, same patients")]
    fig.legend(handles=h, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.07), fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{R}/fig3_fixed_n.png")
    fig.savefig(f"{R}/fig3_fixed_n.pdf")
    plt.close(fig)


def fig4():
    D = P.Design()
    p = D.pats
    groups = [("TRAIN", 1), ("TRAIN", 0), ("TEST", 1), ("TEST", 0)]
    lab_color = {1: "#eb6834", 0: "#2a78d6"}
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    rng = np.random.default_rng(0)
    for i, (s, l) in enumerate(groups):
        n = p.n[(p.split == s) & (p.label == l)].values
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(n)), n, s=9, color=lab_color[l], alpha=0.55, lw=0)
        ax.boxplot(n, positions=[i], widths=0.5, showfliers=False, medianprops=dict(color=INK, lw=1.4),
                   boxprops=dict(color=INK2), whiskerprops=dict(color=INK2), capprops=dict(color=INK2))
        ax.text(i, 2.6, f"n={len(n)}\nmedian {np.median(n):.0f}", ha="center", va="bottom", fontsize=7, color=INK2)
    ax.set_yscale("log")
    ax.set_ylim(2, 8000)
    ax.set_xticks(range(4), ["TRAIN\nMSIMUT", "TRAIN\nMSS", "TEST\nMSIMUT", "TEST\nMSS"])
    ax.set_ylabel("Patches per patient (log scale)")
    ax.grid(axis="x", visible=False)
    fig.savefig(f"{R}/fig4_patch_counts.png")
    fig.savefig(f"{R}/fig4_patch_counts.pdf")
    plt.close(fig)


def main():
    encs = encoders_present()
    m = load_metrics()
    s = build_summary(m, encs)
    conf = s[s.protocol.astype(str).str.startswith("E4_decode")]
    pe = paired_e2(encs)
    pe.to_csv(f"{R}/paired_e2_bootstrap.csv", index=False, float_format="%.4f")
    table(m, encs, conf)
    for extra in (feature_table(), resolution(m)):
        if extra:
            open(f"{R}/tables.md", "a").write("\n" + extra)
    fig1(m, encs)
    if len(conf):
        fig2(conf, encs)
    fig3(m, encs)
    fig4()
    print(open(f"{R}/tables.md").read())
    print(pe.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
