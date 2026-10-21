"""9-class tissue-type probe (Kather 2018 NCT-CRC-HE classes) on frozen UNI embeddings, applied to the MSI tiles.

  python src/tissue.py embed-nct      # UNI CLS for NCT-CRC-HE-100K (train) and CRC-VAL-HE-7K (val), native size
  python src/tissue.py embed-quads    # UNI CLS for the 4 quadrants of each MSI tile upsampled 2x
  python src/tissue.py fit            # LR probe on 100K, accuracy on 7K, then emb/feat_tissue/emb.f16.npy
  python src/tissue.py embed-nct --encoder flash   # NCT embeddings for another encoder (E8 tumor map)
  python src/tissue.py fit-probe --encoder flash   # probe only -> data/nct/probe_<encoder>.joblib

Scale: NCT tiles are 224 px at 0.5 um/px; the MSI tiles are 1.143 um/px (see check_resolution.py). Each MSI tile is upsampled 2x
(bilinear, 448 px) and cut into four 224 px quadrants at ~0.57 um/px, so the probe sees the scale it was trained on.
Per-tile feature = mean of the 4 quadrant class probabilities (9-d); mean pooling over a patient's tiles then gives
soft tissue-type fractions.
"""
import argparse
import glob
import json
import os
import sys
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import torch
from PIL import Image

sys.path.insert(0, os.path.dirname(__file__))
import encoders  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data/nct")
CLASSES = ["ADI", "BACK", "DEB", "LYM", "MUC", "MUS", "NORM", "STR", "TUM"]
ENCODER = "uni"
NO_AUTOCAST = {"flash"}  # encoders that run their own fixed fp16 numerics


def nct_paths(name):
    fs = sorted(glob.glob(os.path.join(ROOT, f"data/nct/{name}/**/*.tif"), recursive=True))
    return fs, np.array([CLASSES.index(os.path.basename(os.path.dirname(f))) for f in fs])


def quads(im):
    big = im.resize((448, 448), Image.BILINEAR)
    return [big.crop((x, y, x + 224, y + 224)) for y in (0, 224) for x in (0, 224)]


def embed(paths, model, tf, fwd, device, bs=256, workers=12, split_quads=False, autocast=True):
    def load(ps):
        out = []
        for p in ps:
            with Image.open(p) as im:
                im = im.convert("RGB")
                out += [tf(q) for q in quads(im)] if split_quads else [tf(im)]
        return torch.stack(out)

    per = 4 if split_quads else 1
    step = bs // per
    batches = [paths[i:i + step] for i in range(0, len(paths), step)]
    outs, t0 = [], time.time()
    with ThreadPoolExecutor(workers) as ex, torch.inference_mode():
        futs = deque(ex.submit(load, b) for b in batches[:2 * workers])
        nxt = len(futs)
        for bi in range(len(batches)):
            x = futs.popleft().result()
            if nxt < len(batches):
                futs.append(ex.submit(load, batches[nxt]))
                nxt += 1
            with torch.autocast(device.type, dtype=torch.float16, enabled=autocast):
                e = fwd(model, x.to(device, non_blocking=True))
            outs.append(e.float().cpu().numpy().astype(np.float16))
            if (bi + 1) % 100 == 0:
                print(f"{(bi + 1) * step}/{len(paths)} {(bi + 1) * step / (time.time() - t0):.0f} tiles/s", flush=True)
    e = np.concatenate(outs)
    return e.reshape(len(paths), per, -1) if split_quads else e


def fit_probe():
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    tr = np.load(os.path.join(OUT, f"{ENCODER}_NCT-CRC-HE-100K.npz"))
    va = np.load(os.path.join(OUT, f"{ENCODER}_CRC-VAL-HE-7K.npz"))
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=3000, random_state=42))
    clf.fit(tr["emb"].astype(np.float32), tr["y"])
    pv = clf.predict(va["emb"].astype(np.float32))
    res = {"encoder": ENCODER, "train": "NCT-CRC-HE-100K", "val": "CRC-VAL-HE-7K", "C": 0.1,
           "val_acc": accuracy_score(va["y"], pv), "val_bal_acc": balanced_accuracy_score(va["y"], pv),
           "classes": CLASSES, "confusion": confusion_matrix(va["y"], pv).tolist()}
    return clf, res


def load_model():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, tf, fwd, dim, info = encoders.load(ENCODER)
    model.eval()
    return model.to(device), tf, fwd, device


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["embed-nct", "embed-quads", "fit", "fit-probe"])
    ap.add_argument("--encoder", default="uni")
    a = ap.parse_args()
    global ENCODER
    ENCODER = a.encoder
    if a.step == "embed-nct":
        model, tf, fwd, device = load_model()
        for name in ("CRC-VAL-HE-7K", "NCT-CRC-HE-100K"):
            ps, y = nct_paths(name)
            e = embed(ps, model, tf, fwd, device, autocast=ENCODER not in NO_AUTOCAST)
            np.savez(os.path.join(OUT, f"{ENCODER}_{name}.npz"), emb=e, y=y)
            print(name, e.shape, flush=True)
    elif a.step == "embed-quads":
        model, tf, fwd, device = load_model()
        idx = pd.read_parquet(os.path.join(ROOT, "data/index.parquet"))
        e = embed([os.path.join(ROOT, p) for p in idx.path], model, tf, fwd, device, split_quads=True)
        np.save(os.path.join(ROOT, f"emb/{ENCODER}/quads2x.f16.npy"), e)
        print(e.shape)
    elif a.step == "fit-probe":
        import joblib
        clf, res = fit_probe()
        joblib.dump(clf, os.path.join(OUT, f"probe_{ENCODER}.joblib"))
        json.dump(res, open(os.path.join(OUT, f"probe_{ENCODER}.json"), "w"), indent=1)
        print(json.dumps({k: res[k] for k in ("val_acc", "val_bal_acc")}))
    else:
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        tr = np.load(os.path.join(OUT, f"{ENCODER}_NCT-CRC-HE-100K.npz"))
        va = np.load(os.path.join(OUT, f"{ENCODER}_CRC-VAL-HE-7K.npz"))
        clf = make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=3000, random_state=42))
        clf.fit(tr["emb"].astype(np.float32), tr["y"])
        pv = clf.predict(va["emb"].astype(np.float32))
        res = {"encoder": ENCODER, "train": "NCT-CRC-HE-100K", "val": "CRC-VAL-HE-7K", "C": 0.1,
               "val_acc": accuracy_score(va["y"], pv), "val_bal_acc": balanced_accuracy_score(va["y"], pv),
               "classes": CLASSES, "confusion": confusion_matrix(va["y"], pv).tolist()}
        print(json.dumps({k: res[k] for k in ("val_acc", "val_bal_acc")}))
        q = np.load(os.path.join(ROOT, f"emb/{ENCODER}/quads2x.f16.npy"), mmap_mode="r")
        n = q.shape[0]
        feats = np.zeros((n, len(CLASSES)), np.float32)
        for s in range(0, n, 20000):
            blk = np.asarray(q[s:s + 20000], np.float32)
            feats[s:s + len(blk)] = clf.predict_proba(blk.reshape(-1, blk.shape[-1])).reshape(len(blk), 4, -1).mean(1)
        out = os.path.join(ROOT, "emb/feat_tissue")
        os.makedirs(out, exist_ok=True)
        np.save(os.path.join(out, "emb.f16.npy"), feats.astype(np.float16))
        json.dump({"encoder": "feat_tissue", "dim": len(CLASSES), "features": CLASSES, "probe": res,
                   "input": f"emb/{ENCODER}/quads2x.f16.npy (2x upsampled quadrants), mean of 4 quadrant probs"},
                  open(os.path.join(out, "meta.json"), "w"), indent=1)
        print("tissue fractions (all tiles):", dict(zip(CLASSES, feats.mean(0).round(3))))


if __name__ == "__main__":
    main()
