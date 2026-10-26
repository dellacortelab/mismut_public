"""Regions for blinded pathologist review: where does a UNI2-h attention-MIL model look on whole slides?

  python src/review_regions.py train      # UNI2-h gated ABMIL on TCGA TRAIN tumor tiles (tile store), 5 seeds
  python src/review_regions.py regions    # attention over ALL tissue tiles of the selected slides -> polygons

Per slide (data/review/selection.csv: held-out TCGA TEST patients):
  - attention a_i over every tissue tile (224 px, 0.5 um/px), mean of 5 seeds; relative attention r_i = a_i * N
  - smoothed over the 3x3 tile neighbourhood
  - HIGH regions: tiles in the top 5% of smoothed attention, 8-connected components with >= MIN_TILES tiles,
    ranked by summed attention, up to N_HIGH per slide
  - LOW control regions: components of tumor-map tiles in the bottom 30% of attention, similar size, up to N_LOW per
    slide (blinded controls for the review)
  - polygon = union of tile squares, level-0 pixel coordinates (GeoJSON)
  - per region: area, relative attention, fraction inside the released (Kather) tumor tiles and inside our tumor map,
    nucleus composition is added later from the review tiles
Outputs: data/review/regions/<barcode>.geojson (blinded labels "R1..Rk"), data/review/regions_key.csv (unblinding key)
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import torch
from scipy import ndimage

sys.path.insert(0, os.path.dirname(__file__))
import iel_fm as I  # noqa: E402
import protocols as P  # noqa: E402
import wsi_eval as W  # noqa: E402

ROOT = P.ROOT
STORE = os.environ.get("MSI_TILES", "/mnt/slides/msimut/tumor_tiles")
OUT = os.path.join(ROOT, "data/review")
SEEDS = (0, 1, 2, 3, 4)
TOP_FRAC, LOW_FRAC, MIN_TILES, N_HIGH, N_LOW = 0.05, 0.30, 4, 6, 3


def cmd_train(a):
    D = P.Design()
    idx = pd.read_parquet(os.path.join(STORE, "index.parquet"))
    E = np.load(os.path.join(STORE, "emb", "uni2_h.f16.npy"), mmap_mode="r")
    pid = {p: i for i, p in enumerate(D.pats.patient_id)}
    tc = idx[idx.cohort == "tcga"]
    tr = sorted(p for p in tc.patient.unique() if pid[p] in set(D.train))
    bags = [np.asarray(E[tc.index.values[(tc.patient == p).values]]) for p in tr]
    y = np.array([D.y[pid[p]] for p in tr])
    os.makedirs(OUT, exist_ok=True)
    for s in SEEDS:
        m, mu, sd = I.train_abmil(bags, y, s)
        torch.save({"state": m.state_dict(), "mu": mu.cpu(), "sd": sd.cpu(), "dim": E.shape[1]},
                   os.path.join(OUT, f"abmil_uni2h_seed{s}.pt"))
        print("trained seed", s, flush=True)


def load_models():
    import abmil as A
    out = []
    for s in SEEDS:
        ck = torch.load(os.path.join(OUT, f"abmil_uni2h_seed{s}.pt"), map_location="cuda")
        m = A.GatedABMIL(ck["dim"]).cuda()
        m.load_state_dict(ck["state"])
        out.append((m.eval(), ck["mu"].cuda(), ck["sd"].cuda()))
    return out


def tile_polygon(g, step0):
    from shapely.geometry import box
    from shapely.ops import unary_union
    from shapely.geometry import Polygon
    u = unary_union([box(x * step0, y * step0, (x + 1) * step0, (y + 1) * step0) for x, y in zip(g.tx, g.ty)])
    u = u.buffer(step0 * 0.02).buffer(-step0 * 0.02)  # merge squares that only touch at corners
    if u.geom_type == "MultiPolygon":
        u = max(u.geoms, key=lambda q: q.area)
    return Polygon(u.exterior).simplify(step0 * 0.05)  # one simple outline, holes filled


def components(mask_grid, min_tiles):
    lab, n = ndimage.label(mask_grid, structure=np.ones((3, 3)))
    return [np.argwhere(lab == k) for k in range(1, n + 1) if (lab == k).sum() >= min_tiles]


def cmd_regions(a):
    from shapely.geometry import mapping
    sel = pd.read_csv(os.path.join(OUT, "selection.csv"))
    plan = pd.read_parquet(W.PLAN)
    reg = pd.read_parquet(os.path.join(ROOT, "data/wsi/regions_flash.parquet"))
    plan["tumor_map"] = reg.tumor_map.values
    models = load_models()
    os.makedirs(os.path.join(OUT, "regions"), exist_ok=True)
    key, rng = [], np.random.default_rng(0)
    for r in sel.itertuples():
        for bc, g in plan[plan.barcode.str[:12] == r.patient].groupby("barcode"):
            g = g.reset_index(drop=True)
            e = np.load(os.path.join(W.WSI, "uni2_h", bc + ".npz"))["emb"]
            X = torch.from_numpy(e.astype(np.float32)).cuda()
            att, ps = [], []
            for m, mu, sd in models:
                p, at = I.abmil_attention(m, (X - mu) / sd, "cuda")
                att.append(at)
                ps.append(p)
            a_ = np.mean(att, 0) * len(g)
            H, Wd = g.ty.max() + 1, g.tx.max() + 1
            grid, have = np.zeros((H, Wd)), np.zeros((H, Wd))
            grid[g.ty, g.tx], have[g.ty, g.tx] = a_, 1
            sm = ndimage.uniform_filter(grid, 3, mode="constant") / np.maximum(ndimage.uniform_filter(have, 3, mode="constant"), 1e-9)
            sm[have == 0] = np.nan
            g["att"] = a_
            g["att_s"] = sm[g.ty, g.tx]
            hi_thr = np.nanquantile(g.att_s, 1 - TOP_FRAC)
            lo_thr = np.nanquantile(g.att_s, LOW_FRAC)
            hi = components(np.nan_to_num(sm) >= hi_thr, MIN_TILES)
            tum_grid = np.zeros((H, Wd), bool)
            tum_grid[g.ty, g.tx] = g.tumor_map.values
            lo = components((np.nan_to_num(sm, nan=np.inf) <= lo_thr) & tum_grid, MIN_TILES)
            pos = {(x, y): i for i, (x, y) in enumerate(zip(g.tx, g.ty))}
            step0 = 224 * 0.5 / g.mpp0.iloc[0]

            def region_rows(comp):
                return g.loc[[pos[(x, y)] for y, x in comp if (x, y) in pos]]
            hi_r = sorted((region_rows(c) for c in hi), key=lambda d: -d.att.sum())[:N_HIGH]
            lo_sizes = [len(d) for d in hi_r] or [MIN_TILES * 2]
            lo_r = [region_rows(c) for c in lo]
            lo_r = [d for d in lo_r if len(d) <= 3 * max(lo_sizes)]
            lo_r = [lo_r[i] for i in rng.permutation(len(lo_r))[:N_LOW]]
            regions = [("high", d) for d in hi_r] + [("low", d) for d in lo_r]
            order = rng.permutation(len(regions))
            feats = []
            for k, j in enumerate(order, 1):
                kind, d = regions[j]
                poly = tile_polygon(d, step0)
                rid = f"R{k}"
                feats.append({"type": "Feature", "geometry": mapping(poly),
                              "properties": {"region_id": rid, "description": ""}})
                key.append(dict(patient=r.patient, barcode=bc, group=r.group, slide_p_msi=float(np.mean(ps)),
                                region_id=rid, kind=kind, n_tiles=len(d), area_mm2=len(d) * (0.112 ** 2),
                                rel_attention=float(d.att.mean()), frac_in_kather=float((d.region == 0).mean()),
                                frac_in_tumor_map=float(d.tumor_map.mean()),
                                x0=float(d.x0.min()), y0=float(d.y0.min())))
            json.dump({"type": "FeatureCollection", "slide": bc, "features": feats},
                      open(os.path.join(OUT, "regions", bc + ".geojson"), "w"))
            print(bc, r.group, f"p={np.mean(ps):.3f}", f"{len(hi_r)} high + {len(lo_r)} low regions", flush=True)
    k = pd.DataFrame(key)
    k.to_csv(os.path.join(OUT, "regions_key.csv"), index=False)
    hi = k[k.kind == "high"]
    print(f"\n{len(hi)} high-attention regions: median area {hi.area_mm2.median():.2f} mm2; "
          f"fraction outside Kather tumor tiles: median {1 - hi.frac_in_kather.median():.2f}; "
          f"outside our tumor map: median {1 - hi.frac_in_tumor_map.median():.2f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["train", "regions"])
    a = ap.parse_args()
    {"train": cmd_train, "regions": cmd_regions}[a.cmd](a)
