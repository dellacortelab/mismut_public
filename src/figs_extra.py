"""Figures 6-7 for the manuscript.

  fig6_resolution_proof: (a) a Zenodo tile next to the same region re-read from the slide; (b) NCC vs candidate pixel
                         size over all 1,092 per-slide scale sweeps (median and IQR).
  fig7_external:         per encoder and model, AUC on TCGA TEST vs CPTAC (and SurGen when present).

  python src/figs_extra.py
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from PIL import Image  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
C_TCGA, C_CPTAC, C_SURGEN = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def fig6():
    sw = pd.read_csv(os.path.join(R, "metrics/scale_sweep.csv"))
    piv = sw.pivot_table(index="row", columns="mpp", values="ncc")
    q = piv.quantile([0.25, 0.5, 0.75])
    idx = pd.read_parquet(os.path.join(ROOT, "data/index.parquet"))
    reg = pd.read_parquet(os.path.join(ROOT, "data/retile/registration.parquet"))
    row = 42280
    fig = plt.figure(figsize=(9, 3.2))
    for j, (img, title) in enumerate(((Image.open(os.path.join(ROOT, idx.path[row])).resize((512, 512), Image.BICUBIC),
                                       "Zenodo tile (224 px; labelled 0.5 µm/px)"),
                                      (Image.open(os.path.join(ROOT, reg.out[row])),
                                       "Same region re-read at 0.5 µm/px"))):
        ax = fig.add_axes([0.0 + j * 0.25, 0.08, 0.24, 0.8])
        ax.imshow(img)
        ax.set_axis_off()
        ax.set_title(title, fontsize=8, color=INK)
    fig.text(0.0, 0.95, "a", fontsize=11, fontweight="bold")
    fig.text(0.55, 0.95, "b", fontsize=11, fontweight="bold")
    ax = fig.add_axes([0.6, 0.18, 0.38, 0.7])
    x = q.columns.values
    ax.fill_between(x, q.loc[0.25], q.loc[0.75], color=C_TCGA, alpha=0.18, lw=0)
    ax.plot(x, q.loc[0.5], color=C_TCGA, lw=2, marker="o", ms=4)
    ax.axvline(256 / 224, color=INK2, lw=0.8, ls="--")
    ax.axvline(0.5, color=INK2, lw=0.8, ls=":")
    ax.text(256 / 224 + 0.01, 0.15, "256 µm / 224 px\n= 1.143 µm/px", fontsize=7.5, color=INK2)
    ax.text(0.51, 0.15, "documented\n0.5 µm/px", fontsize=7.5, color=INK2)
    ax.set_xlabel("Candidate tile pixel size (µm/px)")
    ax.set_ylabel("NCC with source slide")
    ax.set_ylim(0, 1.02)
    ax.set_title(f"{piv.shape[0]:,} tiles, {sw.barcode.nunique()} slides: median and IQR", fontsize=8)
    fig.savefig(os.path.join(R, "fig6_resolution_proof.png"))
    fig.savefig(os.path.join(R, "fig6_resolution_proof.pdf"))
    plt.close(fig)


def fig7():
    lab = {"flash": "Prov-GigaPath-flash", "uni2_h": "UNI2-h", "virchow2": "Virchow2"}
    models = [("tumor_abmil", "Whole tumor, ABMIL"), ("tumor_lr_mean", "Whole tumor, mean"),
              ("kather_abmil", "Kather tiles, ABMIL"), ("kather_lr_mean", "Kather tiles, mean"),
              ("kather_lr_max", "Kather tiles, max"), ("kather_lr_mean+max", "Kather tiles, mean+max")]
    encs = [e for e in lab if os.path.exists(os.path.join(R, f"metrics/external_store_{e}.csv"))]
    fig, axs = plt.subplots(1, len(encs), figsize=(3.3 * len(encs), 3.4), sharey=True, squeeze=False)
    for ax, e in zip(axs[0], encs):
        r = pd.read_csv(os.path.join(R, f"metrics/external_store_{e}.csv")).set_index("model")
        sg = os.path.join(R, f"metrics/external_surgen_{e}.csv")
        s = pd.read_csv(sg).set_index("model") if os.path.exists(sg) else None
        y = np.arange(len(models))[::-1]
        for yi, (m, _) in zip(y, models):
            a, b = r.loc[m, "tcga_test_auc"], r.loc[m, "cptac_auc"]
            ax.plot([a, b], [yi, yi], color="#b8b6b0", lw=1.5, zorder=1)
            ax.scatter(a, yi, s=34, color=C_TCGA, zorder=3, edgecolor="white", lw=0.6)
            ax.scatter(b, yi, s=34, color=C_CPTAC, zorder=3, edgecolor="white", lw=0.6)
            if s is not None and m in s.index:
                ax.scatter(s.loc[m, "surgen_auc"], yi, s=34, color=C_SURGEN, marker="D", zorder=3, edgecolor="white",
                           lw=0.6)
        ax.set_yticks(y, [t for _, t in models])
        ax.set_xlim(0.45, 1.0)
        ax.axvline(0.5, color=INK2, lw=0.6, ls=":")
        ax.set_title(lab[e], fontsize=9)
        ax.set_xlabel("ROC-AUC")
        ax.grid(axis="y", visible=False)
    h = [plt.Line2D([], [], color=C_TCGA, marker="o", ls="", label="TCGA TEST (held-out patients)"),
         plt.Line2D([], [], color=C_CPTAC, marker="o", ls="", label="CPTAC-COAD (external)")]
    if any(os.path.exists(os.path.join(R, f"metrics/external_surgen_{e}.csv")) for e in encs):
        h.append(plt.Line2D([], [], color=C_SURGEN, marker="D", ls="", label="SurGen (external)"))
    fig.legend(handles=h, ncol=len(h), loc="upper center", bbox_to_anchor=(0.55, 1.06))
    fig.savefig(os.path.join(R, "fig7_external.png"))
    fig.savefig(os.path.join(R, "fig7_external.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    fig6()
    fig7()


def fig2_paper():
    """Count decoding R² for the foundation models only (paper Fig. 2)."""
    s = pd.read_csv(os.path.join(R, "summary.csv"))
    encs = ["phikon", "phikon_v2", "uni", "uni2_h", "virchow2", "gigapath"]
    lab = {"phikon": "Phikon", "phikon_v2": "Phikon-v2", "uni": "UNI",
           "uni2_h": "UNI2-h", "virchow2": "Virchow2", "gigapath": "Prov-GigaPath"}
    col = {"mean": "#2a78d6", "max": "#eb6834"}
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.0), sharey=True)
    for k, (ax, coh, title) in enumerate(zip(axs, ("all", "MSS_only"), ("All patients", "MSS patients only"))):
        for j, agg in enumerate(("mean", "max")):
            v = [s[(s.encoder == e) & (s.aggregation == agg) & (s.protocol == f"E4_decode_{coh}") &
                   (s.metric == "count_r2")].value.squeeze() for e in encs]
            x = np.arange(len(encs)) + (j - 0.5) * 0.36
            ax.bar(x, v, width=0.34, color=col[agg], label=agg.capitalize(), edgecolor="white", lw=1)
            for xi, vi in zip(x, v):
                ax.text(xi, vi + 0.01, f"{vi:.2f}", ha="center", va="bottom", fontsize=6.5, color=INK2)
        ax.set_xticks(range(len(encs)), [lab[e] for e in encs], rotation=30, ha="right")
        ax.set_title(title, fontsize=9)
        ax.text(-0.08, 1.06, "ab"[k], transform=ax.transAxes, fontsize=11, fontweight="bold")
        ax.grid(axis="x", visible=False)
        ax.set_ylim(0, 1.05)
    axs[0].set_ylabel("Out-of-fold R², log patch count")
    fig.legend(*axs[0].get_legend_handles_labels(), ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.07))
    fig.savefig(os.path.join(R, "fig2_paper_count_decoding.png"))
    fig.savefig(os.path.join(R, "fig2_paper_count_decoding.pdf"))
    plt.close(fig)


def fig1_paper():
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    import report
    report.fig1(report.load_metrics(), ["phikon", "phikon_v2", "uni", "uni2_h", "virchow2", "gigapath"],
                name="fig1_paper_random_vs_source")


if __name__ == "__main__":
    fig2_paper()
    fig1_paper()


def fig_overview(target=None):
    """Paper Fig. 1: graphical overview. a) how the training split was balanced; b) max pooling tracks tile count;
    c) max vs mean pooling on identical tiles from random CV to two external hospitals."""
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    import protocols as P
    from scipy.stats import spearmanr
    D = P.Design()
    idx = pd.read_parquet(os.path.join(ROOT, "data/index.parquet"))
    C_MSI, C_MSS, C_MAX, C_MEAN = "#eb6834", "#2a78d6", "#eb6834", "#2a78d6"
    fig = target if target is not None else plt.figure(figsize=(12, 4.7))
    # ---- a: real tiles, one MSIMUT and one MSS TRAIN patient
    rng = np.random.default_rng(3)
    pats = D.pats
    pick = {1: pats[(pats.split == "TRAIN") & (pats.label == 1) & (pats.n.between(1100, 1300))].iloc[0],
            0: pats[(pats.split == "TRAIN") & (pats.label == 0) & (pats.n.between(180, 240))].iloc[0]}
    fig.text(0.01, 0.95, "a", fontsize=13, fontweight="bold")
    fig.text(0.035, 0.95, "Training data: MSS tiles discarded", fontsize=10, fontweight="bold", va="center")
    for r, lab in enumerate((1, 0)):
        p = pick[lab]
        rows = rng.choice(np.arange(p.start, p.stop), 12, replace=False)
        y0 = 0.53 - r * 0.42
        fig.text(0.035, y0 + 0.31, f"{'MSI' if lab else 'MSS'} patient", fontsize=9, color=C_MSI if lab else C_MSS,
                 fontweight="bold")
        for j, rw in enumerate(rows):
            ax = fig.add_axes([0.035 + (j % 6) * 0.042, y0 + 0.15 - (j // 6) * 0.13, 0.04, 0.12])
            ax.imshow(Image.open(os.path.join(ROOT, idx.path[rw])))
            ax.set_xticks([])
            ax.set_yticks([])
            kept = lab == 1 or j < 2
            for sp in ax.spines.values():
                sp.set_color(C_MSI if lab else C_MSS)
                sp.set_linewidth(1.2)
            if not kept:
                ax.imshow(np.ones((10, 10, 3)), alpha=0.75, extent=(0, 223, 223, 0))
                ax.plot([0, 223], [0, 223], color="#9a9893", lw=1)
                ax.plot([0, 223], [223, 0], color="#9a9893", lw=1)
                ax.set_xlim(0, 223)
                ax.set_ylim(223, 0)
        txt = "all tiles kept: ≈1,200 per patient" if lab else "most tiles discarded: ≈210 per patient"
        fig.text(0.035, y0 - 0.025, txt, fontsize=8.5, color=INK2)
    fig.text(0.035, 0.01, "Test split and new hospitals: no such imbalance.\nAlso: tiles are 1.14 µm/px, not the "
             "documented 0.5 µm/px.", fontsize=8, color=INK2, va="bottom")
    # ---- b: MSS patients only, mean predicted P(MSI) by tile-count bin (random CV, out-of-fold)
    fig.text(0.34, 0.95, "b", fontsize=13, fontweight="bold")
    fig.text(0.365, 0.95, "Max pooling learns tile count", fontsize=10, fontweight="bold", va="center")
    ax = fig.add_axes([0.385, 0.22, 0.235, 0.62])
    bins = [0, 100, 300, 1000, 1e6]
    blab = ["<100", "100–300", "300–1k", ">1k"]
    for e in ("uni2_h", "virchow2", "gigapath"):
        pr = pd.read_parquet(os.path.join(R, f"preds/{e}_lr.parquet"),
                             filters=[("protocol", "==", "E1rep"), ("scheme", "==", "all")])
        for agg, c in (("max", C_MAX), ("mean", C_MEAN)):
            s_ = pr[pr.aggregation == agg].groupby("pid").p.mean()
            n, y = D.n[s_.index], D.y[s_.index]
            v, nn = s_.values[y == 0], n[y == 0]
            b_ = np.digitize(nn, bins) - 1
            ax.plot(range(4), [v[b_ == k].mean() for k in range(4)], color=c, lw=1.6, marker="o", ms=4, alpha=0.85)
    ax.set_xticks(range(4), blab)
    ax.set_ylim(0, 0.32)
    ax.set_xlabel("Tiles per patient")
    ax.set_ylabel("Mean predicted P(MSI)\nin MSS patients (random CV)")
    ax.grid(axis="x", visible=False)
    ax.text(2.1, 0.27, "max pooling", color=C_MAX, fontsize=9, fontweight="bold", ha="right")
    ax.text(0.0, 0.14, "mean pooling", color=C_MEAN, fontsize=9, fontweight="bold")
    fig.text(0.385, 0.01, "None of these patients has MSI, yet more tiles raise the\nmax-pooling score "
             "(UNI2-h, Virchow2, Prov-GigaPath).", fontsize=8, color=INK2, va="bottom")
    # ---- c: bag-size control on identical whole-tumor tiles: confounded training bags vs equal (200-tile) bags
    fig.text(0.66, 0.95, "c", fontsize=13, fontweight="bold")
    fig.text(0.685, 0.95, "…which degrades it at new hospitals", fontsize=10, fontweight="bold", va="center")
    ax = fig.add_axes([0.69, 0.22, 0.29, 0.62])
    xs = np.arange(3)
    for e in ("flash", "uni2_h", "virchow2"):
        b_ = pd.read_csv(os.path.join(R, f"metrics/bagsize_{e}.csv"))
        b_ = b_[b_.test_bags == "full"].set_index(["train", "pooling", "cohort"]).auc
        for (tr, agg), c, ls in ((("kather_counts", "max"), C_MAX, "-"), (("fixed200", "max"), C_MAX, "--"),
                                 (("kather_counts", "mean"), C_MEAN, "-")):
            ax.plot(xs, [b_[(tr, agg, k)] for k in ("tcga_test", "cptac", "surgen")], color=c, ls=ls, lw=1.6,
                    alpha=0.85, marker="o", ms=4)
    ax.set_xticks(xs, ["Held-out\nTCGA", "CPTAC\n(USA)", "SurGen\n(UK)"])
    ax.set_ylim(0.5, 1.0)
    ax.set_ylabel("ROC-AUC")
    ax.grid(axis="x", visible=False)
    from matplotlib.lines import Line2D
    ax.legend([Line2D([], [], color=C_MAX, lw=1.6), Line2D([], [], color=C_MAX, lw=1.6, ls="--"),
               Line2D([], [], color=C_MEAN, lw=1.6)],
              ["max, released bag sizes", "max, equal bags (200 tiles)", "mean, released bag sizes"],
              loc="lower left", fontsize=7.5, handlelength=2.2)
    fig.text(0.69, 0.01, "Same whole-tumor tiles, fit on TCGA training patients;\nthree encoders (one line each)",
             fontsize=8, color=INK2, va="bottom")
    if target is None:
        fig.savefig(os.path.join(R, "fig0_overview.png"))
        fig.savefig(os.path.join(R, "fig0_overview.pdf"))
        plt.close(fig)


if __name__ == "__main__":
    fig_overview()


def fig_universality(target=None, letters="ab"):
    """Paper Fig. 2: a) random-CV minus held-out AUC per aggregation, one dot per encoder; b) count decoding R²."""
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    import report
    m = report.load_metrics()
    encs = ["phikon", "phikon_v2", "uni", "uni2_h", "virchow2", "gigapath"]
    aggs = [("mean", "Mean"), ("max", "Max"), ("mean+max", "Mean+max"), ("abmil", "ABMIL")]
    col = {"mean": "#2a78d6", "max": "#eb6834", "mean+max": "#1baf7a", "abmil": "#eda100"}
    rng = np.random.default_rng(0)
    fig = target if target is not None else plt.figure(figsize=(9, 3.4))
    axs = fig.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.15, 1]))
    ax = axs[0]
    for i, (agg, lab) in enumerate(aggs):
        v = np.array([report.get(m, e, agg, "E1_random_cv_10x5") - report.get(m, e, agg, "E2_source_split")
                      for e in encs])
        ax.bar(i, v.mean(), width=0.6, color=col[agg], alpha=0.25, edgecolor="none")
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=26, color=col[agg], edgecolor="white", lw=0.6, zorder=3)
        ax.text(i, max(v.max(), 0) + 0.015, f"{v.min():+.2f} to {v.max():+.2f}", ha="center", fontsize=7.5, color=INK2)
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xticks(range(len(aggs)), [l for _, l in aggs])
    ax.set_ylabel("Random-CV AUC minus held-out AUC")
    ax.set_ylim(-0.1, 0.3)
    ax.set_title("Inflation by aggregation (one dot per encoder)", fontsize=9)
    ax.grid(axis="x", visible=False)
    ax.text(-0.12, 1.04, letters[0], transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax = axs[1]
    s = pd.read_csv(os.path.join(R, "summary.csv"))
    groups = [("all", "All patients"), ("MSS_only", "MSS patients only")]
    for g, (coh, glab) in enumerate(groups):
        for j, agg in enumerate(("mean", "max")):
            v = np.array([s[(s.encoder == e) & (s.aggregation == agg) & (s.protocol == f"E4_decode_{coh}") &
                            (s.metric == "count_r2")].value.squeeze() for e in encs])
            x = g * 2.6 + j
            ax.bar(x, v.mean(), width=0.6, color=col[agg], alpha=0.25, edgecolor="none")
            ax.scatter(x + rng.uniform(-0.15, 0.15, len(v)), v, s=26, color=col[agg], edgecolor="white", lw=0.6,
                       zorder=3)
            ax.text(x, v.max() + 0.03, f"{v.min():.2f}–{v.max():.2f}", ha="center", fontsize=7.5, color=INK2)
    ax.set_xticks([0, 1, 2.6, 3.6], ["Mean", "Max", "Mean", "Max"])
    for g, (_, glab) in enumerate(groups):
        ax.text(g * 2.6 + 0.5, -0.2, glab, ha="center", fontsize=8.5, color=INK, transform=ax.transData)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Out-of-fold R², log patch count")
    ax.set_title("Patch count is decodable from max pooling", fontsize=9)
    ax.grid(axis="x", visible=False)
    ax.text(-0.12, 1.04, letters[1], transform=ax.transAxes, fontsize=13, fontweight="bold")
    if target is None:
        fig.tight_layout()
        fig.savefig(os.path.join(R, "fig2_universality.png"))
        fig.savefig(os.path.join(R, "fig2_universality.pdf"))
        plt.close(fig)


if __name__ == "__main__":
    fig_universality()


def fig_biomarker():
    """Paper Fig. 3: a) established biomarker (IEL density) vs FM models on identical tumor tiles, three cohorts;
    b) cell composition of the top-5% attention tiles relative to all tumor tiles (ABMIL, three cohorts)."""
    from matplotlib.lines import Line2D
    m = pd.read_csv(os.path.join(R, "metrics/iel_fm.csv"))
    coh = [("tcga_test", "Held-out TCGA", C_TCGA), ("cptac", "CPTAC (USA)", C_CPTAC), ("surgen", "SurGen (UK)", C_SURGEN)]
    models = [("iel", "IEL density\n(no fitting)"), ("cells_lr", "All nucleus\nfeatures (LR)"),
              ("fm_mean", "FM, mean\npooling"), ("fm_abmil", "FM, ABMIL"), ("late_fusion", "IEL + FM\n(late fusion)")]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), gridspec_kw=dict(width_ratios=[1, 1.05]))
    ax = axs[0]
    for j, (k, lab, c) in enumerate(coh):
        for i, (mk, _) in enumerate(models):
            r = m[(m.cohort == k) & (m.model == mk)].iloc[0]
            x = i + (j - 1) * 0.22
            ax.plot([x, x], [r.ci_low, r.ci_high], color=c, lw=1.2, alpha=0.7)
            ax.scatter(x, r.auc, s=30, color=c, edgecolor="white", lw=0.6, zorder=3, label=lab if i == 0 else None)
    ax.axvspan(-0.5, 1.5, color=GRID, alpha=0.35, lw=0)
    ax.text(0.5, 0.415, "interpretable", ha="center", fontsize=8, color=INK2)
    ax.text(3.0, 0.415, "foundation model", ha="center", fontsize=8, color=INK2)
    ax.set_xticks(range(len(models)), [l for _, l in models], fontsize=7.5)
    ax.set_xlim(-0.5, len(models) - 0.5)
    ax.set_ylim(0.4, 1.06)
    ax.set_ylabel("ROC-AUC (95% CI)")
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper left", fontsize=7.5, ncol=3, handletextpad=0.2, columnspacing=0.8)
    ax.set_title("Same tumor tiles, fit on TCGA training patients only", fontsize=9)
    ax.text(-0.1, 1.04, "a", transform=ax.transAxes, fontsize=12, fontweight="bold")
    enrichment_heatmap(fig, axs[1], "b")
    fig.tight_layout()
    fig.savefig(os.path.join(R, "fig8_biomarker.png"))
    fig.savefig(os.path.join(R, "fig8_biomarker.pdf"))
    plt.close(fig)


def fig_shortcut():
    """Paper Fig. 2: the benchmark audit in one figure. a-c) training-split imbalance, count learning, external
    effect of bag size; d-e) inflation and count decoding across six encoders."""
    fig = plt.figure(figsize=(12, 8.6))
    top, bot = fig.subfigures(2, 1, height_ratios=[4.7, 3.6], hspace=0.04)
    fig_overview(top)
    bot.subplots_adjust(left=0.07, right=0.97, bottom=0.17, top=0.86, wspace=0.28)
    fig_universality(bot, letters="de")
    fig.savefig(os.path.join(R, "fig2_shortcut.png"))
    fig.savefig(os.path.join(R, "fig2_shortcut.pdf"))
    plt.close(fig)


def enrichment_heatmap(fig, ax, letter):
    """Nucleus-class density in the top-5% attention tiles relative to all tumor tiles (attention_enrichment.csv)."""
    coh = [("tcga_test", "Held-out TCGA", C_TCGA), ("cptac", "CPTAC (USA)", C_CPTAC), ("surgen", "SurGen (UK)", C_SURGEN)]
    d = pd.read_csv(os.path.join(R, "metrics/attention_enrichment.csv"))
    cells = [("epithelial", "Epithelial"), ("iel", "Intraepithelial lymphocytes"), ("neutrophil", "Neutrophils"),
             ("eosinophil", "Eosinophils"), ("lymphocyte", "All lymphocytes"), ("plasma", "Plasma cells"),
             ("connective", "Connective tissue")]
    cols = [(y, k) for y in (1, 0) for k, _, _ in coh]
    M = np.array([[d[(d.cohort == k) & (d.y == y)][f"top_{f}"].mean() / d[(d.cohort == k) & (d.y == y)][f"all_{f}"].mean()
                   for y, k in cols] for f, _ in cells])
    L = np.log2(M)
    im = ax.imshow(L, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    for r_ in range(L.shape[0]):
        for c_ in range(L.shape[1]):
            ax.text(c_, r_, f"{M[r_, c_]:.2f}", ha="center", va="center", fontsize=7.5,
                    color="white" if abs(L[r_, c_]) > 1.3 else INK)
    short = {"tcga_test": "TCGA", "cptac": "CPTAC", "surgen": "SurGen"}
    ax.set_xticks(range(len(cols)), [short[k] for _, k in cols], fontsize=7.5)
    ax.set_yticks(range(len(cells)), [l for _, l in cells], fontsize=7.5)
    ax.tick_params(length=0)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.axvline(2.5, color="white", lw=3)
    ax.text(1, -0.85, "MSI patients", ha="center", fontsize=8.5, fontweight="bold")
    ax.text(4, -0.85, "MSS patients", ha="center", fontsize=8.5, fontweight="bold")
    ax.set_ylim(len(cells) - 0.5, -1.2)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, ticks=[-2, -1, 0, 1, 2])
    cb.ax.set_yticklabels(["¼×", "½×", "1×", "2×", "4×"], fontsize=7)
    cb.set_label("Top-5% attention tiles vs all tumor tiles", fontsize=7.5)
    cb.outline.set_visible(False)
    ax.set_title("Where ABMIL attends: nucleus density ratio", fontsize=9, pad=14)
    ax.text(-0.42, 1.04, letter, transform=ax.transAxes, fontsize=12, fontweight="bold")
