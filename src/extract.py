"""Patch embedding extraction. Row-aligned with data/index.parquet. Chunked and resumable.

  python src/extract.py --encoder phikon --bench 1024     # throughput + ETA only
  python src/extract.py --encoder phikon --check-fp16     # fp16 autocast vs fp32 cosine on 256 patches
  python src/extract.py --encoder phikon [--fp16]         # full run
  python src/extract.py --encoder phikon --fp16 --source rt05    # re-tiled regions (retile.py), see SOURCES
"""
import argparse
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

ROOT = os.environ.get("MSI_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CHUNK = 10_000
# Image sources, all row-aligned with data/index.parquet:
#   kather: the Zenodo tiles as published (256 um at 224 px = 1.143 um/px, Macenko) -> emb/<enc>/
#   rt114:  retile.py regions (raw TCGA colour) downsampled to 224 px, i.e. the same 1.143 um/px -> emb/<enc>_rt114/
#   rt05:   retile.py regions at 0.5 um/px: four 224 px crops at offsets (16, 272) of the 512 px region (77% of the
#           area, at the encoders' training resolution); per-crop embeddings in quads.f16.npy, the per-tile embedding
#           is their mean -> emb/<enc>_rt05/
SOURCES = {"kather": "", "rt114": "_rt114", "rt05": "_rt05"}
CROPS = [(x, y) for y in (16, 272) for x in (16, 272)]


def auto_device():
    if torch.cuda.is_available():
        return "cuda"
    return "mps" if torch.backends.mps.is_available() else "cpu"


def _load_batch(paths, tf, source="kather"):
    out = []
    for p in paths:
        with Image.open(p) as im:
            im = im.convert("RGB")
            if source == "rt114":
                out.append(tf(im.resize((224, 224), Image.BOX)))
            elif source == "rt05":
                out += [tf(im.crop((x, y, x + 224, y + 224))) for x, y in CROPS]
            else:
                out.append(tf(im))
    return torch.stack(out)


def run(model, fwd, paths, tf, device, fp16, bs, workers, source="kather"):
    """Decode + transform in a bounded thread pool (PIL releases the GIL; spawned DataLoader workers were ~10x
    slower on macOS), model on device. Keeps 2*workers batches in flight."""
    torch.backends.cuda.matmul.allow_tf32 = False  # keep fp32 path true fp32 on Ampere (A6000)
    torch.backends.cudnn.allow_tf32 = False
    per = len(CROPS) if source == "rt05" else 1
    bs = max(1, bs // per)  # bs counts images on the GPU
    batches = [paths[i:i + bs] for i in range(0, len(paths), bs)]
    outs = []
    run.t_first = None
    with ThreadPoolExecutor(workers) as ex, torch.inference_mode():
        futs = deque(ex.submit(_load_batch, b, tf, source) for b in batches[:2 * workers])
        nxt = len(futs)
        for b in range(len(batches)):
            x = futs.popleft().result()
            if nxt < len(batches):
                futs.append(ex.submit(_load_batch, batches[nxt], tf, source))
                nxt += 1
            if b == 1:
                run.t_first = time.time()  # steady-state timer (excludes first batch)
            x = x.to(device, non_blocking=True)
            if fp16:
                with torch.autocast(device.type, dtype=torch.float16):
                    e = fwd(model, x)
            else:
                e = fwd(model, x)
            e = e.float().cpu().numpy()
            outs.append(e.reshape(-1, per, e.shape[-1]) if per > 1 else e)
    return np.concatenate(outs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder", required=True, choices=list(encoders.ENCODERS))
    ap.add_argument("--bench", type=int, default=0)
    ap.add_argument("--check-fp16", action="store_true")
    ap.add_argument("--fp16", action="store_true")
    ap.add_argument("--bs", type=int, default=64)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--device", default="auto", help="auto = cuda > mps > cpu")
    ap.add_argument("--source", default="kather", choices=list(SOURCES))
    a = ap.parse_args()

    torch.manual_seed(0)
    idx = pd.read_parquet(os.path.join(ROOT, "data/index.parquet"))
    if a.source == "kather":
        paths = [os.path.join(ROOT, p) for p in idx.path]
    else:
        reg = pd.read_parquet(os.path.join(ROOT, "data/retile/registration.parquet"))
        assert (reg.row.values == np.arange(len(idx))).all(), "registration not row-aligned with the index"
        paths = [os.path.join(ROOT, p) for p in reg.out]
    device = torch.device(auto_device() if a.device == "auto" else a.device)
    model, tf, fwd, dim, info = encoders.load(a.encoder)
    model = model.to(device)
    out_dir = os.path.join(ROOT, "emb", a.encoder + SOURCES[a.source])
    os.makedirs(out_dir, exist_ok=True)

    if a.check_fp16:
        rng = np.random.default_rng(0)
        sub = [paths[i] for i in sorted(rng.choice(len(paths), 256, replace=False))]
        e32 = run(model, fwd, sub, tf, device, False, a.bs, a.workers, a.source).reshape(-1, dim)
        e16 = run(model, fwd, sub, tf, device, True, a.bs, a.workers, a.source).reshape(-1, dim)
        cos = (e32 * e16).sum(1) / np.linalg.norm(e32, axis=1) / np.linalg.norm(e16, axis=1)
        res = {"encoder": a.encoder, "cos_min": float(cos.min()), "cos_mean": float(cos.mean()),
               "max_abs_diff": float(np.abs(e32 - e16).max()), "pass": bool(cos.min() > 0.999)}
        print(json.dumps(res))
        json.dump(res, open(os.path.join(out_dir, "fp16_check.json"), "w"), indent=1)
        return

    if a.bench:
        sub = paths[: a.bench]
        run(model, fwd, sub, tf, device, a.fp16, a.bs, a.workers, a.source)
        dt = time.time() - run.t_first
        rate = (a.bench - a.bs) / dt
        res = {"encoder": a.encoder, "fp16": a.fp16, "bs": a.bs, "patches_per_s": rate,
               "eta_min_192312": 192312 / rate / 60}
        print(json.dumps(res))
        json.dump(res, open(os.path.join(out_dir, f"bench_{'fp16' if a.fp16 else 'fp32'}.json"), "w"), indent=1)
        return

    emb_path = os.path.join(out_dir, "emb.f16.npy")
    prog_path = os.path.join(out_dir, "progress.json")
    prog = json.load(open(prog_path)) if os.path.exists(prog_path) else {"done": [], "seconds": 0.0}
    mode = "r+" if os.path.exists(emb_path) else "w+"
    mm = np.lib.format.open_memmap(emb_path, mode=mode, dtype=np.float16, shape=(len(paths), dim))
    qm = None
    if a.source == "rt05":
        qp = os.path.join(out_dir, "quads.f16.npy")
        qm = np.lib.format.open_memmap(qp, mode="r+" if os.path.exists(qp) else "w+", dtype=np.float16,
                                       shape=(len(paths), len(CROPS), dim))
    n_chunks = (len(paths) + CHUNK - 1) // CHUNK
    for c in range(n_chunks):
        if c in prog["done"]:
            continue
        s, e = c * CHUNK, min((c + 1) * CHUNK, len(paths))
        t = time.time()
        x = run(model, fwd, paths[s:e], tf, device, a.fp16, a.bs, a.workers, a.source)
        assert np.isfinite(x).all() and np.abs(x).max() < 65000, "non-finite or fp16 overflow"
        if qm is not None:
            qm[s:e] = x.astype(np.float16)
            qm.flush()
            x = x.mean(1)
        mm[s:e] = x.astype(np.float16)
        mm.flush()
        dt = time.time() - t
        prog["done"].append(c)
        prog["seconds"] += dt
        json.dump(prog, open(prog_path, "w"))
        print(f"[{a.encoder}] chunk {c + 1}/{n_chunks} {(e - s) / dt:.1f} patches/s", flush=True)

    meta = {**{k: v for k, v in info.items() if k != "processor"}, "encoder": a.encoder + SOURCES[a.source],
            "base_encoder": a.encoder, "source": a.source, "crops": CROPS if a.source == "rt05" else None,
            "processor": info.get("processor"), "dim": dim, "n": len(paths),
            "inference_dtype": "float16 autocast" if a.fp16 else "float32", "storage_dtype": "float16",
            "device": str(device), "batch_size": a.bs, "torch": torch.__version__,
            "extract_seconds": prog["seconds"], "patches_per_s": len(paths) / prog["seconds"],
            "index_sha": pd.util.hash_pandas_object(idx[["path"]]).sum().item()}
    try:
        import timm
        import transformers
        meta.update(timm=timm.__version__, transformers=transformers.__version__)
    except ImportError:
        pass
    json.dump(meta, open(os.path.join(out_dir, "meta.json"), "w"), indent=1, default=str)
    print("done", json.dumps({k: meta[k] for k in ("encoder", "revision", "patches_per_s")}))


if __name__ == "__main__":
    main()
