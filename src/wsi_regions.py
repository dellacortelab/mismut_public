"""Whole-slide tiling by region (tumor / margin / distant), for the "where is the MSI signal" test (E8).

The Kather tiles cover only tumor regions picked by an automatic detector; the invasive margin and the host
response around it (Crohn-like lymphoid reaction, peritumoral lymphocytes) are not in the dataset. We know where
every Kather tile sits in its slide (retile.py), so we can tile the whole tissue and label each tile:
  tumor   >= 50% of the tile inside the Kather tile mask
  margin  not tumor, centre within MARGIN_UM of the tumor mask
  distant everything else that is tissue

Lattice: 224 px tiles at 0.5 um/px (112 um), anchored at the slide origin (Prov-GigaPath-flash native lattice).
Only the slides the Kather tiles came from are used (364 slides, 360 patients).

  python src/wsi_regions.py plan            # tissue + region masks at 8 um/px -> data/wsi/plan.parquet, QA figure
  python src/wsi_regions.py embed --encoder flash [--workers 24]
                                            # -> /mnt/slides/msimut/wsi/<encoder>/<barcode>.npz (coords, region, emb)
"""
import argparse
import glob
import json
import os
import sys
import time

import cv2
import numpy as np
import pandas as pd
import torch
from PIL import Image

sys.path.insert(0, os.path.dirname(__file__))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SLIDES = os.environ.get("MSI_SLIDES", "/mnt/slides/msimut/tcga_crc_dx")
OUT = os.environ.get("MSI_WSI", "/mnt/slides/msimut/wsi")
PLAN = os.path.join(ROOT, "data/wsi/plan.parquet")
MASK_UM = 8.0                       # mask resolution
TILE_PX, TILE_MPP = 224, 0.5        # 112 um tiles
CELLS = int(TILE_PX * TILE_MPP / MASK_UM)   # 14 mask px per tile edge
MARGIN_UM = 1000.0
MIN_TISSUE = 0.5
REGIONS = ["tumor", "margin", "distant"]
BLOCK = 8                           # tiles per block edge when reading


def slide_path(barcode):
    fs = glob.glob(os.path.join(SLIDES, barcode + ".*.svs"))
    assert len(fs) == 1, (barcode, fs)
    return fs[0]


def tissue_mask(rgb):
    """rgb uint8 at MASK_UM. Tissue = saturated, not white, not pen ink (green to azure hues)."""
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)  # H in [0,180)
    h, s, v = hsv[..., 0].astype(int) * 2, hsv[..., 1], hsv[..., 2]
    tissue = (s > 20) & (v < 245) & (v > 30)
    ink = (h >= 60) & (h <= 215) & (s > 50)
    m = (tissue & ~ink).astype(np.uint8)
    m = cv2.morphologyEx(m, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))
    return m.astype(bool)


def plan_slide(barcode, g):
    import openslide

    s = openslide.OpenSlide(slide_path(barcode))
    mpp0 = float(s.properties.get("openslide.mpp-x") or s.properties["aperio.MPP"])
    W0, H0 = s.dimensions
    ds = MASK_UM / mpp0
    lvl = s.get_best_level_for_downsample(ds + 1e-6)
    w, h = int(W0 / ds), int(H0 / ds)
    img = np.asarray(s.read_region((0, 0), lvl, s.level_dimensions[lvl]).convert("RGB"))
    img = cv2.resize(img, (w, h), interpolation=cv2.INTER_AREA)
    tis = tissue_mask(img)
    tum = np.zeros((h, w), bool)
    for t in g.itertuples():  # Kather tiles: 256 um boxes, level-0 top-left
        x, y = int(round(t.x0 * mpp0 / MASK_UM)), int(round(t.y0 * mpp0 / MASK_UM))
        tum[y:y + 32, x:x + 32] = True
    dist = cv2.distanceTransform((~tum).astype(np.uint8), cv2.DIST_L2, 5) * MASK_UM  # um to nearest tumor px
    nty, ntx = h // CELLS, w // CELLS
    crop = (slice(0, nty * CELLS), slice(0, ntx * CELLS))

    def per_tile(a, f):
        return f(a[crop].reshape(nty, CELLS, ntx, CELLS), axis=(1, 3))

    tis_f = per_tile(tis.astype(np.float32), np.mean)
    tum_f = per_tile(tum.astype(np.float32), np.mean)
    dist_c = dist[CELLS // 2::CELLS, CELLS // 2::CELLS][:nty, :ntx]
    ty, tx = np.where(tis_f >= MIN_TISSUE)
    region = np.where(tum_f[ty, tx] >= 0.5, 0, np.where(dist_c[ty, tx] <= MARGIN_UM, 1, 2))
    step0 = TILE_PX * TILE_MPP / mpp0  # tile edge in level-0 px
    df = pd.DataFrame(dict(barcode=barcode, tx=tx, ty=ty, x0=tx * step0, y0=ty * step0, region=region,
                           dist_um=dist_c[ty, tx], tissue=tis_f[ty, tx], mpp0=mpp0))
    qa = dict(img=cv2.resize(img, (w // 4, h // 4)), tis=tis[::4, ::4], tum=tum[::4, ::4], h=h, w=w)
    return df, qa


def cmd_plan(a):
    reg = pd.read_csv(os.path.join(ROOT, "results/tile_slide_registration.csv.gz"))
    os.makedirs(os.path.dirname(PLAN), exist_ok=True)
    from concurrent.futures import ProcessPoolExecutor
    bcs = sorted(reg.slide.unique())
    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(plan_slide, bcs, [reg[reg.slide == b] for b in bcs]))
    plan = pd.concat([r[0] for r in res], ignore_index=True)
    plan.to_parquet(PLAN, index=False)
    c = plan.groupby(["barcode", "region"]).size().unstack(fill_value=0)
    kc = reg.groupby("slide").size()
    print(f"{len(plan):,} tiles over {plan.barcode.nunique()} slides")
    print("tiles per region:", plan.region.value_counts().sort_index().rename(dict(enumerate(REGIONS))).to_dict())
    print("per slide median:", c.median().rename(dict(enumerate(REGIONS))).to_dict())
    # tumor tiles vs Kather tiles: each Kather 256 um tile covers ~5.2 of our 112 um tiles
    print("tumor tiles / Kather tiles, median per slide:", float((c[0] / kc.reindex(c.index)).median()))
    qa_fig([r[1] for r in res[:12]], bcs[:12])


def qa_fig(qas, names):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(3, 4, figsize=(20, 12))
    for ax, q, n in zip(axs.flat, qas, names):
        img = q["img"].astype(np.float32)
        tis, tum = q["tis"][:img.shape[0], :img.shape[1]], q["tum"][:img.shape[0], :img.shape[1]]
        ov = img.copy()
        ov[tis & ~tum] = 0.5 * ov[tis & ~tum] + 0.5 * np.array([255, 160, 0])
        ov[tum] = 0.5 * ov[tum] + 0.5 * np.array([40, 120, 230])
        ax.imshow(ov.astype(np.uint8))
        ax.set_title(n, fontsize=8)
        ax.axis("off")
    fig.suptitle("blue = Kather tumor tiles, orange = other tissue in mask", fontsize=11)
    fig.savefig(os.path.join(ROOT, "data/wsi/plan_qa.png"), dpi=60, bbox_inches="tight")


class Blocks(torch.utils.data.IterableDataset):
    """Yields (slide index, row indices, tile tensor) for BLOCK x BLOCK groups of planned tiles."""

    def __init__(self, plan, tf):
        self.plan, self.tf = plan, tf
        items = []
        for si, (bc, g) in enumerate(plan.groupby("barcode", sort=True)):
            bx, by = g.tx.values // BLOCK, g.ty.values // BLOCK
            for (_, _), gg in g.groupby([by, bx]):
                items.append((si, bc, gg.index.values))
        self.items = items

    def __iter__(self):
        import openslide
        wi = torch.utils.data.get_worker_info()
        items = self.items if wi is None else self.items[wi.id::wi.num_workers]
        cur, s = None, None
        for si, bc, rows in items:
            if bc != cur:
                s, cur = openslide.OpenSlide(slide_path(bc)), bc
            g = self.plan.loc[rows]
            mpp0 = g.mpp0.iloc[0]
            step0 = TILE_PX * TILE_MPP / mpp0
            tx0, ty0 = g.tx.min(), g.ty.min()
            ntx, nty = g.tx.max() - tx0 + 1, g.ty.max() - ty0 + 1
            x0, y0 = int(round(tx0 * step0)), int(round(ty0 * step0))
            w0, h0 = int(round(ntx * step0)), int(round(nty * step0))
            ds = TILE_MPP / mpp0
            lvl = s.get_best_level_for_downsample(ds + 1e-6)
            lds = s.level_downsamples[lvl]
            blk = np.asarray(s.read_region((x0, y0), lvl, (int(np.ceil(w0 / lds)), int(np.ceil(h0 / lds))))
                             .convert("RGB"))
            blk = cv2.resize(blk, (ntx * TILE_PX, nty * TILE_PX), interpolation=cv2.INTER_AREA)
            tiles = [self.tf(Image.fromarray(blk[(ty - ty0) * TILE_PX:(ty - ty0 + 1) * TILE_PX,
                                                 (tx - tx0) * TILE_PX:(tx - tx0 + 1) * TILE_PX]))
                     for tx, ty in zip(g.tx.values, g.ty.values)]
            yield si, torch.from_numpy(rows), torch.stack(tiles)


def cmd_embed(a):
    import encoders

    plan = pd.read_parquet(PLAN)
    if a.limit_slides:
        plan = plan[plan.barcode.isin(sorted(plan.barcode.unique())[:a.limit_slides])].reset_index(drop=True)
    if a.patients:
        want = set(open(a.patients).read().split())
        plan = plan[plan.barcode.str[:12].isin(want)].reset_index(drop=True)
    out = os.path.join(OUT, a.encoder)
    os.makedirs(out, exist_ok=True)
    todo = [b for b in sorted(plan.barcode.unique()) if not os.path.exists(os.path.join(out, b + ".npz"))]
    plan = plan[plan.barcode.isin(todo)].reset_index(drop=True)
    print(f"{len(todo)} slides, {len(plan):,} tiles to embed with {a.encoder}", flush=True)
    model, tf, fwd, dim, info = encoders.load(a.encoder)
    device = torch.device("cuda")
    model = model.to(device).eval()
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    emb = np.zeros((len(plan), dim), np.float16)
    done = np.zeros(len(plan), bool)
    left = plan.groupby("barcode").size().to_dict()
    dl = torch.utils.data.DataLoader(Blocks(plan, tf), batch_size=None, num_workers=a.workers, prefetch_factor=4,
                                     persistent_workers=False)
    t0, n = time.time(), 0
    with torch.inference_mode():
        for si, rows, x in dl:
            x = x.to(device, non_blocking=True)
            if a.fp16:
                with torch.autocast("cuda", dtype=torch.float16):
                    e = fwd(model, x)
            else:
                e = fwd(model, x)
            rows = rows.numpy()
            emb[rows] = e.float().cpu().numpy().astype(np.float16)
            done[rows] = True
            n += len(rows)
            bc = plan.barcode.values[rows[0]]
            left[bc] -= len(rows)
            if left[bc] == 0:
                g = plan[plan.barcode == bc]
                assert done[g.index].all()
                np.savez(os.path.join(out, bc + ".npz"), x0=g.x0.values, y0=g.y0.values, region=g.region.values,
                         dist_um=g.dist_um.values, emb=emb[g.index.values])
            if n % 50000 < len(rows):
                print(f"{n:,}/{len(plan):,} {n / (time.time() - t0):.0f} tiles/s", flush=True)
    json.dump({**{k: v for k, v in info.items() if k != "processor"}, "encoder": a.encoder, "dim": dim,
               "fp16_autocast": a.fp16, "tiles": len(plan), "seconds": time.time() - t0},
              open(os.path.join(out, "meta.json"), "w"), indent=1, default=str)
    print("done", n, f"{time.time() - t0:.0f}s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "embed"])
    ap.add_argument("--encoder", default="flash")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--fp16", action="store_true", help="autocast (for the large encoders; flash runs its own fp16)")
    ap.add_argument("--limit-slides", type=int, default=0)
    ap.add_argument("--patients", help="file with TCGA patient ids to restrict to")
    a = ap.parse_args()
    {"plan": cmd_plan, "embed": cmd_embed}[a.cmd](a)


if __name__ == "__main__":
    main()
