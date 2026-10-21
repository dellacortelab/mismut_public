"""Nucleus detection + classification with KongNet_CoNIC_1 (tiatoolbox; trained on Lizard/CoNIC colon tissue at
0.5 um/px). Runs in the separate .venv-tia environment.

Scale: Zenodo 2530835 documents the MSI tiles as 0.5 um/px, but they are coarser. KongNet works at native size on
NCT-CRC-HE tiles (documented 0.5 um/px; LYM tiles -> 87% lymphocytes) and badly under-detects on the MSI tiles;
the tiles are 256 um at 224 px = 1.143 um/px (Kather 2019 Methods; registration in check_resolution.py). MSI tiles
are upsampled SCALE=2x (bilinear, ~0.57 um/px) before detection. nct_val tiles run at native size.

  .venv-tia/bin/python src/cells.py --set kather      # -> emb/feat_cells/emb.f16.npy, row-aligned with the index
  .venv-tia/bin/python src/cells.py --set nct_val     # -> data/nct/cells_val.npz (sanity check vs tissue labels)
  .venv-tia/bin/python src/cells.py --set kather --check-fp16

Per tile: (upsampled) tile reflect-padded to a multiple of 32 (model input), peaks detected per class head (local max, window
2*min_distance+1, prob > threshold_abs, both from the model config), then cross-class suppression (keep a peak only
if it is the highest-probability peak of any class within the same window). Only peaks inside the central 224 px
are counted. Per-tile features (counts per tile; MSI tiles cover ~0.05 mm^2):
  0-5  Neutrophil, Epithelial, Lymphocyte, Plasma, Eosinophil, Connective
  6    intraepithelial lymphocyte proxy: lymphocytes with an epithelial nucleus within IEL_UM micrometres
  7    all nuclei
Detections (tile row, y, x, class) are saved alongside for QA.
"""
import argparse
import glob
import json
import os
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import torch
import torch.nn.functional as F
from PIL import Image
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLASSES = ["neutrophil", "epithelial", "lymphocyte", "plasma", "eosinophil", "connective"]
FEATURES = CLASSES + ["iel", "nuclei"]
MODEL_MPP, IEL_UM = 0.5, 10.0  # coordinates live on the model-input grid (~0.5 um/px)
SCALE = {"kather": 2, "nct_val": 1}
MEAN = torch.tensor([0.485, 0.456, 0.406]).view(1, 3, 1, 1)
STD = torch.tensor([0.229, 0.224, 0.225]).view(1, 3, 1, 1)


def paths_for(which):
    if which == "kather":
        import pandas as pd
        idx = pd.read_parquet(os.path.join(ROOT, "data/index.parquet"))
        return [os.path.join(ROOT, p) for p in idx.path], None
    fs = sorted(glob.glob(os.path.join(ROOT, "data/nct/CRC-VAL-HE-7K/**/*.tif"), recursive=True))
    return fs, [os.path.basename(os.path.dirname(f)) for f in fs]


def load_batch(paths):
    out = []
    for p in paths:
        with Image.open(p) as im:
            a = np.asarray(im.convert("RGB"))
        assert a.shape == (224, 224, 3), (p, a.shape)
        out.append(a)
    return torch.from_numpy(np.stack(out))


def detect(model, x_uint8, device, fp16, scale):
    """x_uint8 [B,224,224,3] -> list of (y, x, cls) int arrays per tile (coordinates on the upsampled grid)."""
    x = x_uint8.to(device).permute(0, 3, 1, 2).float() / 255.0
    if scale != 1:
        x = F.interpolate(x, scale_factor=scale, mode="bilinear", align_corners=False)
    x = (x - MEAN.to(device)) / STD.to(device)
    PAD = 16 * scale  # 224*s + 32*s is a multiple of 32
    x = F.pad(x, (PAD, PAD, PAD, PAD), mode="reflect")
    with torch.autocast(device.type, dtype=torch.float16, enabled=fp16):
        logits = model(x)
    p = torch.sigmoid(logits[:, model.target_channels].float())  # [B,6,256,256]
    k = 2 * model.min_distance + 1
    peak = (p == F.max_pool2d(p, k, 1, k // 2)) & (p > model.threshold_abs)
    pp = torch.where(peak, p, torch.zeros_like(p))
    best, cls = pp.max(1)  # strongest peak value per pixel over classes
    keep = (best > 0) & (best == F.max_pool2d(best[:, None], k, 1, k // 2)[:, 0])
    keep = keep[:, PAD:-PAD, PAD:-PAD]
    cls = cls[:, PAD:-PAD, PAD:-PAD]
    b, yy, xx = torch.nonzero(keep, as_tuple=True)
    c = cls[b, yy, xx]
    b, yy, xx, c = (t.cpu().numpy() for t in (b, yy, xx, c))
    splits = np.searchsorted(b, np.arange(1, x.shape[0]))
    return [(y_, x_, c_) for y_, x_, c_ in zip(np.split(yy, splits), np.split(xx, splits), np.split(c, splits))]


def tile_features(y, x, c):
    f = np.zeros(len(FEATURES), np.float32)
    f[:6] = np.bincount(c, minlength=6)
    f[7] = len(c)
    lym, epi = c == 2, c == 1
    if lym.any() and epi.any():
        d, _ = cKDTree(np.c_[y[epi], x[epi]]).query(np.c_[y[lym], x[lym]], k=1)
        f[6] = (d * MODEL_MPP <= IEL_UM).sum()
    return f


def run(model, paths, device, fp16, bs, workers, scale, log_every=0):
    feats, dets = [], []
    batches = [paths[i:i + bs] for i in range(0, len(paths), bs)]
    t0 = time.time()
    with ThreadPoolExecutor(workers) as ex, torch.inference_mode():
        futs = deque(ex.submit(load_batch, b) for b in batches[:2 * workers])
        nxt = len(futs)
        for bi in range(len(batches)):
            x = futs.popleft().result()
            if nxt < len(batches):
                futs.append(ex.submit(load_batch, batches[nxt]))
                nxt += 1
            for j, (y_, x_, c_) in enumerate(detect(model, x, device, fp16, scale)):
                feats.append(tile_features(y_, x_, c_))
                row = bi * bs + j
                dets.append(np.c_[np.full(len(c_), row), y_, x_, c_].astype(np.int32))
            if log_every and (bi + 1) % log_every == 0:
                done = (bi + 1) * bs
                print(f"{done}/{len(paths)} {done / (time.time() - t0):.0f} tiles/s", flush=True)
    return np.stack(feats), np.concatenate(dets)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["kather", "nct_val"], required=True)
    ap.add_argument("--bs", type=int, default=32)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--fp16", action="store_true")
    ap.add_argument("--check-fp16", action="store_true")
    a = ap.parse_args()
    from tiatoolbox.models.architecture import get_pretrained_model

    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, _ = get_pretrained_model("KongNet_CoNIC_1")
    model = model.to(device).eval()
    paths, labels = paths_for(a.set)

    if a.check_fp16:
        rng = np.random.default_rng(0)
        sub = [paths[i] for i in sorted(rng.choice(len(paths), 512, replace=False))]
        f32, _ = run(model, sub, device, False, a.bs, a.workers, SCALE[a.set])
        f16, _ = run(model, sub, device, True, a.bs, a.workers, SCALE[a.set])
        res = {"nuclei_fp32": float(f32[:, 7].sum()), "nuclei_fp16": float(f16[:, 7].sum()),
               "per_feature_corr": {n: float(np.corrcoef(f32[:, i], f16[:, i])[0, 1]) for i, n in enumerate(FEATURES)},
               "mean_abs_count_diff": float(np.abs(f32 - f16)[:, :6].sum(1).mean())}
        res["pass"] = bool(min(res["per_feature_corr"].values()) > 0.99)
        print(json.dumps(res, indent=1))
        os.makedirs(os.path.join(ROOT, "emb/feat_cells"), exist_ok=True)
        json.dump(res, open(os.path.join(ROOT, "emb/feat_cells/fp16_check.json"), "w"), indent=1)
        return

    t = time.time()
    feats, dets = run(model, paths, device, a.fp16, a.bs, a.workers, SCALE[a.set], log_every=200)
    dt = time.time() - t
    info = {"model": "tiatoolbox KongNet_CoNIC_1 (TIACentre/KongNet_pretrained_weights)", "features": FEATURES,
            "scale": SCALE[a.set], "model_mpp": MODEL_MPP, "iel_um": IEL_UM, "min_distance": model.min_distance,
            "threshold_abs": model.threshold_abs, "inference_dtype": "float16 autocast" if a.fp16 else "float32",
            "n": len(paths), "seconds": dt, "tiles_per_s": len(paths) / dt, "torch": torch.__version__}
    if a.set == "nct_val":
        np.savez(os.path.join(ROOT, "data/nct/cells_val.npz"), feats=feats, labels=np.array(labels),
                 paths=np.array([os.path.relpath(p, ROOT) for p in paths]))
    else:
        out = os.path.join(ROOT, "emb/feat_cells")
        os.makedirs(out, exist_ok=True)
        np.save(os.path.join(out, "emb.f16.npy"), feats.astype(np.float16))
        np.save(os.path.join(out, "detections.npy"), dets)
        json.dump({"encoder": "feat_cells", "dim": len(FEATURES), **info}, open(os.path.join(out, "meta.json"), "w"),
                  indent=1)
    print(json.dumps({k: info[k] for k in ("n", "seconds", "tiles_per_s")}))


if __name__ == "__main__":
    main()
