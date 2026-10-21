"""Gated-attention ABMIL (Ilse et al. 2018) on frozen patch embeddings, run on the same jobs as eval_lr.py.

Architecture: Linear(d,256)-ReLU-Dropout(0.25) -> gated attention (tanh/sigmoid branches, 256 hidden, dropout 0.25)
-> softmax over patches -> attention-weighted mean -> Linear(256,1).
Training: AdamW lr 1e-4, wd 1e-4, fixed 30 epochs, one bag per step, class-balanced sampling with replacement
(n_train draws per epoch), BCE loss. No early stopping, no use of test folds. Patch features are standardized with
mean/SD of training-patient patches (in-fold).
Seeds: E1, E2, all-patch comparator use seeds 0,1,2; E1rep, fixed-N and count-matched use seed = repeat index.

  python src/abmil.py --encoder phikon [--jobs 10] [--bench]
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from joblib import Parallel, delayed

sys.path.insert(0, os.path.dirname(__file__))
import pool  # noqa: E402
import protocols as P  # noqa: E402

HID, DROP, LR, WD, EPOCHS = 256, 0.25, 1e-4, 1e-4, 30


class GatedABMIL(nn.Module):
    def __init__(self, d):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(d, HID), nn.ReLU(), nn.Dropout(DROP))
        self.a = nn.Sequential(nn.Linear(HID, HID), nn.Tanh(), nn.Dropout(DROP))
        self.b = nn.Sequential(nn.Linear(HID, HID), nn.Sigmoid(), nn.Dropout(DROP))
        self.c = nn.Linear(HID, 1)
        self.cls = nn.Linear(HID, 1)

    def forward(self, x):  # x [n, d]
        h = self.fc(x)
        w = torch.softmax(self.c(self.a(h) * self.b(h)).squeeze(-1), dim=0)
        return self.cls(w @ h).squeeze(-1)


def fit_predict(emb, plan, tr, ytr, te, seed, device="cpu"):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    # in-fold feature standardization from (up to 50k) training patches
    rows = np.concatenate([plan[i] for i in tr])
    samp = np.sort(rng.choice(rows, min(50_000, len(rows)), replace=False))
    xs = np.asarray(emb[samp], dtype=np.float32)
    mu, sd = xs.mean(0), xs.std(0) + 1e-6
    # materialize this job's bags once: one normalized float32 tensor + per-patient slices
    need = np.concatenate([tr, te])
    slices, parts, o = {}, [], 0
    for i in need:
        r = plan[i]
        x = emb[r[0]:r[-1] + 1] if len(r) == r[-1] - r[0] + 1 else emb[r]
        parts.append((np.asarray(x, dtype=np.float32) - mu) / sd)
        slices[i] = (o, o + len(r))
        o += len(r)
    X = torch.from_numpy(np.concatenate(parts)).to(device)
    del parts

    def bag(i):
        a, b = slices[i]
        return X[a:b]

    model = GatedABMIL(emb.shape[1]).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    lossf = nn.BCEWithLogitsLoss()
    w = np.where(ytr == 1, 1.0 / (ytr == 1).sum(), 1.0 / (ytr == 0).sum())
    w /= w.sum()
    model.train()
    for _ in range(EPOCHS):
        for k in rng.choice(len(tr), len(tr), replace=True, p=w):
            opt.zero_grad()
            loss = lossf(model(bag(tr[k])), torch.tensor(float(ytr[k]), device=device))
            loss.backward()
            opt.step()
    model.eval()
    with torch.inference_mode():
        return np.array([torch.sigmoid(model(bag(i))).item() for i in te])


def resolve_device(device, slot):
    """'auto' -> cuda:<slot % n_gpus> if CUDA is available, else cpu (MPS per-step overhead made it slower)."""
    if device == "auto":
        return f"cuda:{slot % torch.cuda.device_count()}" if torch.cuda.is_available() else "cpu"
    return device


def run_job(encoder, j, device, threads, slot=0):
    torch.set_num_threads(threads)
    device = resolve_device(device, slot)
    D = P.Design()
    emb = pool.load_emb(encoder)
    plan = P.get_plan(D, j["plan"])
    p = fit_predict(emb, plan, j["tr"], D.y[j["tr"]], j["te"], j["seed"], device)
    return P.job_rows(j, p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder", required=True)
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--device", default="auto", help="auto = round-robin over CUDA GPUs, else cpu")
    ap.add_argument("--bench", action="store_true", help="time one full-bag and one N=100 fit, then exit")
    a = ap.parse_args()
    D = P.Design()
    jobs = P.make_jobs(D, ["abmil"], {"single": (0, 1, 2), "per_repeat": lambda r: (r,)},
                       paper_folds=a.encoder in ("phikon", "phikon_ref"))
    if a.bench:
        for want in [dict(scheme="all", protocol="E1"), dict(scheme="fixedN", N=100, protocol="E1")]:
            j = next(j for j in jobs if all(j[k] == v for k, v in want.items()))
            t = time.time()
            run_job(a.encoder, j, a.device, a.threads)
            print(want, f"{time.time() - t:.1f}s", flush=True)
        print(len(jobs), "jobs total")
        return
    out_f = os.path.join(P.ROOT, f"results/preds/{a.encoder}_abmil.parquet")
    # resumable: cache per-job predictions
    part_dir = os.path.join(P.ROOT, f"results/preds/_abmil_parts/{a.encoder}")
    os.makedirs(part_dir, exist_ok=True)

    def key(j):
        return f"{j['scheme']}_{j['N']}_{j['repeat']}_{j['protocol']}_{j['cv_repeat']}_{j['fold']}_{j['seed']}"

    todo = [j for j in jobs if not os.path.exists(os.path.join(part_dir, key(j) + ".parquet"))]
    # big bags first for load balance
    todo.sort(key=lambda j: -(j["plan"][0] == "all"))
    print(f"{len(todo)}/{len(jobs)} jobs to run", flush=True)
    t0 = time.time()

    def go(j, slot):
        df = run_job(a.encoder, j, a.device, a.threads, slot)
        df.to_parquet(os.path.join(part_dir, key(j) + ".parquet"), index=False)
        return 1

    done = 0
    for _ in Parallel(n_jobs=a.jobs, return_as="generator_unordered")(delayed(go)(j, i) for i, j in enumerate(todo)):
        done += 1
        if done % 20 == 0:
            print(f"{done}/{len(todo)} {time.time() - t0:.0f}s", flush=True)
    preds = pd.concat([pd.read_parquet(os.path.join(part_dir, key(j) + ".parquet")) for j in jobs], ignore_index=True)
    preds.insert(0, "encoder", a.encoder)
    preds.to_parquet(out_f, index=False)
    m = P.summarize(preds, D.y)
    m.to_csv(os.path.join(P.ROOT, f"results/metrics/{a.encoder}_abmil.csv"), index=False)
    print(m[(m.metric == "roc_auc") & (m.scheme == "all")][["protocol", "value", "sd", "CI_low", "CI_high"]]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
