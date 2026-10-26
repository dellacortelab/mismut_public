"""TCGA label audit: released MSIMUT ("MSI or highly mutated") labels vs TCGA PanCancer Atlas molecular data.

  python src/label_audit.py   # needs data/review/tcga_molecular.csv (cBioPortal coadread_tcga_pan_can_atlas_2018:
                              # TMB_NONSYNONYMOUS, MSI_SCORE_MANTIS, MSI_SENSOR_SCORE, POLE mutations), and
                              # results/metrics/op_preds.csv from src/operating_points.py
Discordant labels: MSS-labelled patients that are hypermutated (TMB >= 10/Mb) or MSI by MSIsensor (>= 10), and
MSIMUT-labelled patients with TMB < 10, MSIsensor < 3.5 and MANTIS < 0.4.
Outputs: results/metrics/label_audit.csv and printed summary.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import protocols as P  # noqa: E402
from external import auc_ci  # noqa: E402

ROOT = P.ROOT
HOT = {"P286R", "V411L", "S459F", "A456P", "S297F", "F367S", "L424V", "M444K", "P436R", "D368Y"}  # POLE exonuclease hotspots


def main():
    D = P.Design()
    pats = D.pats.reset_index(drop=True)
    lab = pd.DataFrame({"patient": pats.patient_id, "y": D.y,
                        "split": np.where(np.isin(np.arange(len(pats)), D.test), "TEST", "TRAIN")}).set_index("patient")
    m = pd.read_csv(os.path.join(ROOT, "data/review/tcga_molecular.csv")).set_index("patient")
    m["pole_hot"] = m.pole.fillna("").apply(lambda s: ";".join(h for h in s.split(";") if h in HOT))
    a = lab.join(m[["tmb", "mantis", "sensor", "pole_hot"]])
    a.to_csv(os.path.join(ROOT, "results/metrics/label_audit.csv"))
    a["discordant"] = ""
    a.loc[(a.y == 0) & (a.tmb >= 10), "discordant"] = "MSS label, hypermutated"
    a.loc[(a.y == 0) & (a.tmb < 10) & (a.sensor >= 10), "discordant"] = "MSS label, MSIsensor >= 10"
    a.loc[(a.y == 1) & (a.tmb < 10) & (a.sensor < 3.5) & (a.mantis < 0.4), "discordant"] = "MSIMUT label, no MSI or hypermutation"
    a.to_csv(os.path.join(ROOT, "results/metrics/label_audit.csv"))
    dis = a[a.discordant != ""].sort_values(["discordant", "tmb"], ascending=[True, False])
    print(f"MSS median TMB {a[a.y == 0].tmb.median():.1f}; discordant labels: {len(dis)} "
          f"({(dis.split == 'TEST').sum()} in TEST)")
    print(dis.round(3).to_string())
    p = pd.read_csv(os.path.join(ROOT, "results/metrics/op_preds.csv"))
    t = p[(p.cohort == "tcga_test") & (p.model == "FM ABMIL")].set_index("patient")
    te = [k for k in dis.index if k in t.index]
    print("TCGA test FM ABMIL (5 seeds, <=1,000 tiles):", round(auc_ci(t.y.values, t.score.values)[0], 3))
    v = t.copy()
    v.loc[te, "y"] = 1 - v.loc[te, "y"]
    print(f"discordant TEST labels flipped ({len(te)}):", round(auc_ci(v.y.values, v.score.values)[0], 3))
    u = t.drop(te)
    print("discordant TEST patients excluded:", round(auc_ci(u.y.values, u.score.values)[0], 3), "n", len(u))
    print(t.loc[te, ["y", "score"]].join(a[["tmb", "mantis", "sensor", "pole_hot", "discordant"]]).round(4).to_string())


if __name__ == "__main__":
    main()
