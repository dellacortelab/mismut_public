"""TCGA molecular annotation from cBioPortal (study coadread_tcga_pan_can_atlas_2018), public API, no login.

  python src/fetch_molecular.py   # -> data/review/tcga_molecular.csv (all patients), data/review/molecular.csv (review patients)

Per patient: TMB_NONSYNONYMOUS, MSI_SCORE_MANTIS, MSI_SENSOR_SCORE, SUBTYPE and POLE protein changes (';'-joined).
"""
import collections
import os

import pandas as pd
import requests

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
REV = os.path.join(ROOT, "data/review")
API = "https://www.cbioportal.org/api"
STUDY = "coadread_tcga_pan_can_atlas_2018"
ATTR = {"TMB_NONSYNONYMOUS": "tmb", "MSI_SCORE_MANTIS": "mantis", "MSI_SENSOR_SCORE": "sensor", "SUBTYPE": "subtype"}
POLE = 5426  # Entrez gene ID


def main():
    clin = requests.get(f"{API}/studies/{STUDY}/clinical-data",
                        params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"}, timeout=300).json()
    muts = requests.post(f"{API}/molecular-profiles/{STUDY}_mutations/mutations/fetch", params={"projection": "SUMMARY"},
                         json={"entrezGeneIds": [POLE], "sampleListId": f"{STUDY}_all"}, timeout=300).json()
    rows = collections.defaultdict(dict)
    for x in clin:
        if x["clinicalAttributeId"] in ATTR:
            rows[x["sampleId"][:12]][ATTR[x["clinicalAttributeId"]]] = x["value"]
    pole = collections.defaultdict(list)
    for x in muts:
        pole[x["sampleId"][:12]].append(x.get("proteinChange"))
    m = pd.DataFrame([{"patient": p, **{c: v.get(c) for c in ATTR.values()}, "pole": ";".join(pole.get(p, []))}
                      for p, v in rows.items()])
    os.makedirs(REV, exist_ok=True)
    m.to_csv(os.path.join(REV, "tcga_molecular.csv"), index=False)
    print(f"{len(m)} patients -> data/review/tcga_molecular.csv")

    sel = os.path.join(REV, "regions_key.csv")
    if os.path.exists(sel):
        g = pd.read_csv(sel)[["patient", "group"]].drop_duplicates()
        r = g.merge(m[["patient", "tmb", "mantis", "pole"]], on="patient", how="left")
        r.to_csv(os.path.join(REV, "molecular.csv"), index=False)
        print(f"{len(r)} review patients -> data/review/molecular.csv")


if __name__ == "__main__":
    main()
