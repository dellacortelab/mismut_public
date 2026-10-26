"""Build data/index.parquet (one row per patch) from the Kather CRC_DX folders and verify against the paper.

  python src/build_index.py                                        # CRC_DX (verified against the target paper)
  MSI_ROOT=/mnt/slides/msimut/stad python src/build_index.py --cohort STAD   # verified against the Zenodo tile counts
"""
import os
import re
import sys

import pandas as pd

ROOT = os.environ.get("MSI_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "data/raw")
PAT = re.compile(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})")

EXPECTED = {  # split, label -> (patients, patches)
    ("TRAIN", 1): (39, 46704), ("TRAIN", 0): (221, 46704),
    ("TEST", 1): (26, 28335), ("TEST", 0): (74, 70569),
}


STAD_TILES = {("TRAIN", 1): 50285, ("TRAIN", 0): 50285, ("TEST", 1): 27904, ("TEST", 0): 90104}  # Zenodo text


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="CRC_DX", choices=["CRC_DX", "STAD"])
    a = ap.parse_args()
    rows = []
    for split in ("TRAIN", "TEST"):
        for lab_name, lab in (("MSIMUT", 1), ("MSS", 0)):
            top = os.path.join(RAW, f"{a.cohort}_{split}_{lab_name}")
            for dp, _, fs in os.walk(top):
                for f in fs:
                    if f.startswith(".") or not f.lower().endswith((".png", ".jpg", ".jpeg", ".tif")):
                        continue
                    m = PAT.search(f)
                    if not m:
                        sys.exit(f"cannot parse patient id: {f}")
                    rows.append((os.path.relpath(os.path.join(dp, f), ROOT), m.group(1), split, lab))
    df = pd.DataFrame(rows, columns=["path", "patient_id", "split", "label"])
    df = df.sort_values(["split", "label", "patient_id", "path"], ascending=[False, False, True, True])
    df = df.reset_index(drop=True)
    df["label"] = df.label.astype("int8")

    ok = True
    print(f"patients {df.patient_id.nunique()}, patches {len(df)}")
    if a.cohort == "CRC_DX":
        ok &= df.patient_id.nunique() == 360 and len(df) == 192312
    if a.cohort == "STAD":
        for (s, l), nx in STAD_TILES.items():
            d = df[(df.split == s) & (df.label == l)]
            print(f"{s:5s} label={l}: patients {d.patient_id.nunique()}, patches {len(d)} (Zenodo {nx}), "
                  f"mean/patient {len(d) / max(1, d.patient_id.nunique()):.0f}")
            ok &= len(d) == nx
    for (s, l), (np_, nx) in (EXPECTED.items() if a.cohort == "CRC_DX" else []):
        d = df[(df.split == s) & (df.label == l)]
        good = d.patient_id.nunique() == np_ and len(d) == nx
        ok &= good
        print(f"{s:5s} label={l}: patients {d.patient_id.nunique()} (exp {np_}), patches {len(d)} (exp {nx}), "
              f"mean/patient {len(d) / max(1, d.patient_id.nunique()):.0f} {'OK' if good else 'MISMATCH'}")
    per = df.groupby("patient_id").agg(ns=("split", "nunique"), nl=("label", "nunique"))
    leak = int((per.ns > 1).sum() + (per.nl > 1).sum())
    print(f"patients in >1 split or >1 label: {leak}")
    ok &= leak == 0
    ok &= not df.path.duplicated().any()
    if not ok:
        sys.exit("VERIFICATION FAILED: stop and report")
    df.to_parquet(os.path.join(ROOT, "data/index.parquet"), index=False)
    print("wrote data/index.parquet")


if __name__ == "__main__":
    main()
