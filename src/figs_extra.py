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
    encs = ["phikon_ref", "phikon", "phikon_v2", "uni", "uni2_h", "virchow2", "gigapath"]
    lab = {"phikon_ref": "Phikon (archived)", "phikon": "Phikon", "phikon_v2": "Phikon-v2", "uni": "UNI",
           "uni2_h": "UNI2-h", "virchow2": "Virchow2", "gigapath": "Prov-GigaPath"}
    col = {"mean": "#2a78d6", "max": "#eb6834"}
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.0), sharey=True)
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
    report.fig1(report.load_metrics(), ["phikon_ref", "phikon", "phikon_v2", "uni", "uni2_h", "virchow2", "gigapath"],
                name="fig1_paper_random_vs_source")


if __name__ == "__main__":
    fig2_paper()
    fig1_paper()
