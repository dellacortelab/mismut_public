"""Paper Fig. 4: UNI2-h ABMIL attention over whole slides (all tissue, 0.5 um/px) for one MSI and one MSS held-out
TCGA patient, with the high-attention regions sent to pathologist review and views of the top regions.

  python src/fig_attention.py   # needs data/review/abmil_uni2h_seed*.pt and the whole-slide UNI2-h embeddings
"""
import glob
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
import numpy as np  # noqa: E402
import openslide  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from PIL import Image  # noqa: E402
from scipy import ndimage  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import figs_extra as FX  # noqa: E402,F401  (rcParams)
import iel_fm as I  # noqa: E402
import review_regions as RR  # noqa: E402
import wsi_eval as W  # noqa: E402

SLIDES = os.environ.get("MSI_SLIDES", "/mnt/slides/msimut") + "/tcga_crc_dx"
CASES = [("TCGA-CK-6746-01Z-00-DX1", "MSI patient"), ("TCGA-CL-5917-01Z-00-DX1", "MSS patient")]


def main():
    plan = pd.read_parquet(W.PLAN)
    key = pd.read_csv(os.path.join(RR.OUT, "regions_key.csv"))
    models = RR.load_models()
    whole = plt.figure(figsize=(12, 10.2))
    fig, bot = whole.subfigures(2, 1, height_ratios=[7.4, 2.8], hspace=0.0)
    for row, (bc, lab) in enumerate(CASES):
        g = plan[plan.barcode == bc].reset_index(drop=True)
        e = np.load(os.path.join(W.WSI, "uni2_h", bc + ".npz"))["emb"]
        X = torch.from_numpy(e.astype(np.float32)).cuda()
        res = [I.abmil_attention(m, (X - mu) / sd, "cuda") for m, mu, sd in models]
        a = np.mean([r[1] for r in res], 0) * len(g)
        p = float(np.mean([r[0] for r in res]))
        sl = openslide.OpenSlide(glob.glob(os.path.join(SLIDES, bc + ".*.svs"))[0])
        th = sl.get_thumbnail((1600, 1600))
        k = th.size[0] / sl.dimensions[0]
        step0 = 224 * 0.5 / g.mpp0.iloc[0]
        H, Wd = g.ty.max() + 1, g.tx.max() + 1
        grid = np.full((H, Wd), np.nan)
        grid[g.ty, g.tx] = np.log2(np.maximum(a, 1e-3))
        y0 = 0.53 - row * 0.5
        ax = fig.add_axes([0.0, y0, 0.27, 0.42])
        ax.imshow(th)
        ax.set_title(f"{lab}: H&E", fontsize=9)
        ax.axis("off")
        ax = fig.add_axes([0.27, y0, 0.27, 0.42])
        ax.imshow(np.asarray(th.convert("L")), cmap="gray", vmin=0, vmax=255, alpha=0.5)
        ext = (0, Wd * step0 * k, H * step0 * k, 0)
        im = ax.imshow(np.ma.masked_invalid(grid), cmap="magma", vmin=-3, vmax=3, extent=ext, interpolation="nearest")
        hk = key[(key.barcode == bc) & (key.kind == "high")].sort_values("rel_attention", ascending=False)
        fc = {f["properties"]["region_id"]: f for f in json.load(open(os.path.join(RR.OUT, "regions", bc + ".geojson")))["features"]}
        for n, r in enumerate(hk.itertuples(), 1):
            xy = np.array(fc[r.region_id]["geometry"]["coordinates"][0]) * k
            ax.plot(xy[:, 0], xy[:, 1], color="#00e5ff", lw=1.4,
                    path_effects=[pe.Stroke(linewidth=2.6, foreground="black"), pe.Normal()])
            if n <= 3:
                ax.text(xy[:, 0].max() + 10, xy[:, 1].mean(), str(n), color="white", fontsize=12, fontweight="bold",
                        va="center", path_effects=[pe.withStroke(linewidth=3, foreground="black")])
        ax.set_xlim(0, th.size[0])
        ax.set_ylim(th.size[1], 0)
        ax.axis("off")
        ax.set_title(f"UNI2-h attention, P(MSI) = {p:.2f}; cyan: high-attention regions", fontsize=9)
        if row == 0:
            cax = fig.add_axes([0.30, 0.50, 0.21, 0.012])
            cb = fig.colorbar(im, cax=cax, orientation="horizontal")
            cb.set_ticks([-3, 0, 3], labels=["⅛× mean", "mean", "8× mean"])
            cb.ax.tick_params(labelsize=7)
        mpp0 = g.mpp0.iloc[0]
        for n, r in enumerate(hk.head(3).itertuples()):
            poly = np.array(fc[r.region_id]["geometry"]["coordinates"][0])
            cx, cy = poly.mean(0)
            size0 = int(400 / mpp0)  # 400 um field
            lvl = sl.get_best_level_for_downsample(0.5 / mpp0)
            ds = sl.level_downsamples[lvl]
            im_ = sl.read_region((int(cx - size0 / 2), int(cy - size0 / 2)), lvl,
                                 (int(size0 / ds), int(size0 / ds))).convert("RGB").resize((800, 800))
            ax = fig.add_axes([0.56 + n * 0.148, y0 + 0.08, 0.142, 0.26])
            ax.imshow(im_)
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(True)
                sp.set_color("#00a6c0")
                sp.set_linewidth(1.5)
            ax.set_title(f"Region {n + 1}", fontsize=8.5)
            if n == 0:
                ax.plot([40, 40 + 200], [760, 760], color="white", lw=2.5)
                ax.text(140, 740, "100 µm", color="white", ha="center", fontsize=7)
    for t, x in (("a", 0.0), ("b", 0.27), ("c", 0.56)):
        fig.text(x, 0.98, t, fontsize=13, fontweight="bold")
    # d: what the attended tiles contain, three cohorts (Prov-GigaPath-flash ABMIL; attention_enrichment.py)
    FX.enrichment_heatmap(bot, bot.add_axes([0.30, 0.06, 0.42, 0.74]), "d")
    whole.savefig(os.path.join(FX.R, "fig9_attention.png"), dpi=200)
    whole.savefig(os.path.join(FX.R, "fig9_attention.pdf"), dpi=300)


def slide_attention(bc, models, plan):
    g = plan[plan.barcode == bc].reset_index(drop=True)
    e = np.load(os.path.join(W.WSI, "uni2_h", bc + ".npz"))["emb"]
    X = torch.from_numpy(e.astype(np.float32)).cuda()
    res = [I.abmil_attention(m, (X - mu) / sd, "cuda") for m, mu, sd in models]
    a = np.mean([r[1] for r in res], 0) * len(g)
    grid = np.full((g.ty.max() + 1, g.tx.max() + 1), np.nan)
    grid[g.ty, g.tx] = np.log2(np.maximum(a, 1e-3))
    return grid, 224 * 0.5 / g.mpp0.iloc[0]


def draw_released_tiles(fig, rect, barcode, color, title):
    """Slide thumbnail with the released (Kather) tiles of that slide marked."""
    from matplotlib.patches import Rectangle
    reg = pd.read_csv(os.path.join(FX.R, "tile_slide_registration.csv.gz"))
    t = reg[reg.slide == barcode]
    sl = openslide.OpenSlide(glob.glob(os.path.join(SLIDES, barcode + ".*.svs"))[0])
    th = sl.get_thumbnail((1400, 1400))
    k = th.size[0] / sl.dimensions[0]
    ax = fig.add_axes(rect)
    ax.imshow(th)
    w = t.region_um.iloc[0] / t.mpp0.iloc[0] * k
    for x0, y0 in zip(t.x0 * k, t.y0 * k):
        ax.add_patch(Rectangle((x0, y0), w, w, facecolor="#00c2ff", edgecolor="#003a66", lw=0.3, alpha=0.75))
    a = np.asarray(th)
    ys, xs = np.where((a.min(2) < 200) & (a.max(2) - a.min(2) > 25))  # tissue
    x_lo, x_hi, y_lo, y_hi = (np.percentile(xs, 0.5), np.percentile(xs, 99.5), np.percentile(ys, 0.5),
                              np.percentile(ys, 99.5))
    cx, cy, half = (x_lo + x_hi) / 2, (y_lo + y_hi) / 2, max(x_hi - x_lo, y_hi - y_lo) / 2 + 15
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy + half, cy - half)
    ax.set_facecolor("white")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    return len(t)


def fig_recipe():
    """Paper Fig. 1: three routes to MSI stratification from H&E, the benchmark pitfalls, and their performance."""
    from matplotlib.patches import FancyBboxPatch
    C_MSI, C_MSS = "#eb6834", "#2a78d6"
    fig = plt.figure(figsize=(12, 8.4))
    heads = [("a", "Established marker", "count intraepithelial lymphocytes (IEL)"),
             ("b", "Foundation-model embeddings", "pool tile embeddings per patient: mean | max | attention"),
             ("c", "Attention-guided discovery", "let attention point, pathologists name")]
    for (l, h, sub), x in zip(heads, (0.0, 0.335, 0.69)):
        fig.text(x, 0.985, l, fontsize=14, fontweight="bold", va="top")
        fig.text(x + 0.025, 0.985, h, fontsize=11, fontweight="bold", va="top")
        fig.text(x + 0.025, 0.957, sub, fontsize=8.5, color=FX.INK2, va="top")
    # a: IEL detections on a high-attention field (224 um)
    z = np.load(os.path.join(RR.OUT, "overview_crop.npz"))
    ax = fig.add_axes([0.02, 0.47, 0.29, 0.45])
    ax.imshow(z["img"])
    lym = (z["c"] == 2) & ~z["iel"]
    ax.scatter(z["x"][lym], z["y"][lym], s=30, facecolor="none", edgecolor="#00e5ff", lw=1.0)
    ax.scatter(z["x"][z["iel"]], z["y"][z["iel"]], s=55, facecolor="none", edgecolor="#ffe600", lw=1.6)
    ax.plot([20, 220], [425, 425], color="white", lw=3)
    ax.text(120, 412, "100 µm", color="white", ha="center", fontsize=8,
            path_effects=[pe.withStroke(linewidth=2, foreground="black")])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(f"yellow: IEL ({int(z['iel'].sum())}); cyan: other lymphocytes", fontsize=8, color=FX.INK2)
    # b: released tiles on two training slides of similar tumor size, then the two pitfalls
    n1 = draw_released_tiles(fig, [0.345, 0.64, 0.155, 0.235], "TCGA-F4-6570-01Z-00-DX1", C_MSI, "")
    n0 = draw_released_tiles(fig, [0.51, 0.64, 0.155, 0.235], "TCGA-F4-6569-01Z-00-DX1", C_MSS, "")
    fig.text(0.4225, 0.925, f"MSI training patient\n{n1:,} tiles released", ha="center", va="top", fontsize=8.5,
             color=C_MSI, fontweight="bold")
    fig.text(0.5875, 0.925, f"MSS training patient\n{n0:,} tiles released", ha="center", va="top", fontsize=8.5,
             color=C_MSS, fontweight="bold")
    fig.text(0.505, 0.668, "cyan: tiles in the released dataset (similar tumor size, same site)", ha="center", va="top",
             fontsize=8, color=FX.INK2)
    for k_, (t, body) in enumerate((("Pitfall 1: tile-count shortcut",
                                     "MSS training tiles were discarded to balance classes;\nmax pooling learns tile count."),
                                    ("Pitfall 2: wrong resolution",
                                     "Tiles documented and reused as 0.5 µm/px\nare 1.14 µm/px."))):
        y = 0.585 - k_ * 0.075
        fig.patches.append(FancyBboxPatch((0.345, y - 0.03), 0.32, 0.062, boxstyle="round,pad=0.004,rounding_size=0.01",
                                          transform=fig.transFigure, facecolor="#fdf0ea", edgecolor=C_MSI, lw=1))
        fig.text(0.355, y + 0.024, t, fontsize=8.5, fontweight="bold", color=C_MSI, va="top")
        fig.text(0.355, y + 0.003, body, fontsize=8, color=FX.INK, va="top")
    # c: attention map with review regions
    plan = pd.read_parquet(W.PLAN)
    models = RR.load_models()
    bc = "TCGA-CK-6746-01Z-00-DX1"
    grid, step0 = slide_attention(bc, models, plan)
    sl = openslide.OpenSlide(glob.glob(os.path.join(SLIDES, bc + ".*.svs"))[0])
    th = sl.get_thumbnail((1400, 1400))
    k = th.size[0] / sl.dimensions[0]
    ax = fig.add_axes([0.70, 0.47, 0.29, 0.45])
    ax.imshow(np.asarray(th.convert("L")), cmap="gray", vmin=0, vmax=255, alpha=0.5)
    H, Wd = grid.shape
    ax.imshow(np.ma.masked_invalid(grid), cmap="magma", vmin=-3, vmax=3, extent=(0, Wd * step0 * k, H * step0 * k, 0),
              interpolation="nearest")
    key = pd.read_csv(os.path.join(RR.OUT, "regions_key.csv"))
    fc = {f["properties"]["region_id"]: f for f in json.load(open(os.path.join(RR.OUT, "regions", bc + ".geojson")))["features"]}
    for r in key[(key.barcode == bc) & (key.kind == "high")].itertuples():
        xy = np.array(fc[r.region_id]["geometry"]["coordinates"][0]) * k
        ax.plot(xy[:, 0], xy[:, 1], color="#00e5ff", lw=1.3,
                path_effects=[pe.Stroke(linewidth=2.5, foreground="black"), pe.Normal()])
    ax.set_xlim(0.08 * th.size[0], 1.0 * th.size[0])
    ax.set_ylim(0.95 * th.size[1], 0.3 * th.size[1])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UNI2-h ABMIL attention (bright = high); cyan: regions\nsent to blinded review by two pathologists",
                  fontsize=8, color=FX.INK2)
    # d: all approaches on identical tumor tiles, three cohorts, 95% CIs
    m = pd.concat([pd.read_csv(os.path.join(FX.R, "metrics/iel_fm.csv")),
                   pd.read_csv(os.path.join(FX.R, "metrics/iel_fm_max.csv"))]).set_index(["cohort", "model"])
    coh = [("tcga_test", "Held-out TCGA", FX.C_TCGA), ("cptac", "CPTAC (USA)", FX.C_CPTAC),
           ("surgen", "SurGen (UK)", FX.C_SURGEN)]
    models_ = [("iel", "IEL density\n(no fitting)"), ("iel_ratio", "IEL per 100\nepithelial nuclei"),
               ("cells_lr", "All nucleus\nclasses (LR)"), ("fm_max_released", "Max pooling,\nreleased counts"),
               ("fm_max_equal", "Max pooling,\nequal bags"), ("fm_mean", "Mean\npooling"), ("fm_abmil", "Attention\nMIL"),
               ("late_fusion", "IEL + mean\npooling")]
    ax = fig.add_axes([0.06, 0.08, 0.92, 0.27])
    for j, (k_, lab, c) in enumerate(coh):
        for i, (mk, _) in enumerate(models_):
            r = m.loc[(k_, mk)]
            x = i + (j - 1) * 0.2
            ax.plot([x, x], [r.ci_low, r.ci_high], color=c, lw=1.2, alpha=0.7)
            ax.scatter(x, r.auc, s=32, color=c, edgecolor="white", lw=0.6, zorder=3, label=lab if i == 0 else None)
    for lo, hi, t in ((-0.5, 2.5, "interpretable (a)"), (2.5, 6.5, "foundation-model embeddings (b, c)"),
                      (6.5, 7.5, "combined")):
        ax.axvspan(lo, hi, color=FX.GRID, alpha=0.35 if lo < 0 or lo > 6 else 0, lw=0)
        ax.text((lo + hi) / 2, 1.02, t, ha="center", fontsize=8.5, color=FX.INK2)
    ax.axvspan(2.5, 3.5, color="#fdf0ea", lw=0)
    ax.text(3.0, 0.42, "pitfall 1", ha="center", fontsize=8, color=C_MSI, fontweight="bold")
    ax.set_xticks(range(len(models_)), [l for _, l in models_], fontsize=8)
    ax.set_xlim(-0.5, len(models_) - 0.5)
    ax.set_ylim(0.4, 1.05)
    ax.set_ylabel("ROC-AUC (95% CI)")
    ax.grid(axis="x", visible=False)
    ax.legend(loc="lower right", fontsize=8, ncol=3, handletextpad=0.2, columnspacing=0.8)
    fig.text(0.0, 0.405, "d", fontsize=14, fontweight="bold", va="top")
    fig.text(0.025, 0.405, "Which approach stratifies MSI across hospitals?", fontsize=11, fontweight="bold", va="top")
    fig.text(0.025, 0.378, "Identical tumor tiles (≤1,000 per patient; Prov-GigaPath-flash embeddings, KongNet nuclei); "
             "all models fit on TCGA training patients only", fontsize=8.5, color=FX.INK2, va="top")
    fig.savefig(os.path.join(FX.R, "fig1_recipe.png"), dpi=200)
    fig.savefig(os.path.join(FX.R, "fig1_recipe.pdf"), dpi=300)


if __name__ == "__main__":
    fig_recipe() if sys.argv[1:] == ["recipe"] else main()
