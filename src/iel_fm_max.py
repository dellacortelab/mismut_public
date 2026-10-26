"""Max pooling on the SAME tumor tiles as iel_fm.py (up to 1,000 per patient, Prov-GigaPath-flash), for the overview
figure: does the tile-count shortcut appear when everything else is held identical?

  fm_max_released  LR on max-pooled embeddings; each TCGA TRAIN bag subsampled to the patient's released (Kather) tile
                   count (capped by the tiles available), i.e. the confounded bag sizes
  fm_max_equal     the same with 200 tiles per TRAIN patient
Test/external patients use all their tiles. Ten training draws; predicted probabilities averaged over draws.

  python src/iel_fm_max.py   # -> results/metrics/iel_fm_max.csv
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import eval_lr  # noqa: E402
import iel_fm as I  # noqa: E402
import protocols as P  # noqa: E402

DRAWS, EQUAL = 10, 200


def main():
    D = P.Design()
    pid = {p: i for i, p in enumerate(D.pats.patient_id)}
    T = I.merge(I.tcga())
    cl = pd.read_csv(os.path.join(I.ROOT, "data/cptac/labels.csv"), index_col=0)["label"]
    C = {k: v for k, v in I.merge(I.cptac()).items() if k in cl.index}
    sl = pd.read_csv(os.path.join(I.ROOT, "data/surgen/subset.csv"), dtype={"case_id": str}).set_index("case_id")
    S = {k: v for k, v in I.merge(I.surgen()).items() if k in sl.index}
    tr = sorted(k for k in T if pid[k] in set(D.train))
    te = sorted(k for k in T if pid[k] in set(D.test))
    ytr = np.array([D.y[pid[k]] for k in tr])
    ev = {"tcga_test": ([T[k][1] for k in te], np.array([D.y[pid[k]] for k in te])),
          "cptac": ([C[k][1] for k in sorted(C)], np.array([int(cl[k]) for k in sorted(C)])),
          "surgen": ([S[k][1] for k in sorted(S)], np.array([int(sl.loc[k, "mmr_loss_binary"]) for k in sorted(S)]))}
    Xev = {c: np.stack([b.astype(np.float32).max(0) for b in bags]) for c, (bags, _) in ev.items()}
    n_rel = np.array([D.n[pid[k]] for k in tr])
    n_av = np.array([len(T[k][1]) for k in tr])
    print(f"TRAIN bag sizes used (released, capped): MSI mean {np.minimum(n_rel, n_av)[ytr == 1].mean():.0f}, "
          f"MSS mean {np.minimum(n_rel, n_av)[ytr == 0].mean():.0f}")
    rows = []
    for name, sizes in (("fm_max_released", np.minimum(n_rel, n_av)), ("fm_max_equal", np.minimum(EQUAL, n_av))):
        P_ = {c: [] for c in ev}
        for d in range(DRAWS):
            rng = np.random.default_rng(d)
            Xtr = np.stack([T[k][1][rng.choice(len(T[k][1]), s, replace=False)].astype(np.float32).max(0)
                            for k, s in zip(tr, sizes)])
            m = eval_lr.lr().fit(Xtr, ytr)
            for c in ev:
                P_[c].append(m.predict_proba(Xev[c])[:, 1])
        for c, (_, y) in ev.items():
            p = np.mean(P_[c], 0)
            a, lo, hi, _ = I.auc_ci(y, p)
            rows.append(dict(cohort=c, model=name, auc=a, ci_low=lo, ci_high=hi, frac_p_gt_099=float((p > 0.99).mean()),
                             n=len(y), pos=int(y.sum())))
            print(name, c, f"{a:.3f} [{lo:.2f}, {hi:.2f}]", f"P>0.99: {(p > 0.99).mean():.2f}", flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(I.ROOT, "results/metrics/iel_fm_max.csv"), index=False)


if __name__ == "__main__":
    main()
