"""Assemble the review part of the Zenodo deposit (paper/zenodo_deposit/review/) from data/review/.

Raters are anonymised (A/B), platform IDs and timestamps are dropped, question columns get their survey meaning.
"""
import json
import shutil
from pathlib import Path

import pandas as pd

SRC = Path("data/review")
OUT = Path("paper/zenodo_deposit/review")
COLS = {"q1_growth_pattern": "main_tissue", "q2_lumen_status": "inflammation", "q3_gland_separation": "morphology",
        "q8_cytology": "cytology", "q4_gp_subtype": "msi_impression", "q7_notes": "notes",
        "time_spent_seconds": "seconds"}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    key = pd.read_csv(SRC / "regions_key.csv")
    key.drop(columns=["slide_p_msi"]).rename(columns={"code": "blind_code", "region_id": "region"}).to_csv(
        OUT / "regions.csv", index=False)

    d = pd.read_csv(SRC / "survey_export_2026-10-25.csv")
    raters = {u: f"rater_{c}" for u, c in zip(sorted(d.username.unique()), "AB")}
    d["rater"] = d.username.map(raters)
    d["blind_code"] = d.annotation_id.str.split(":").str[1]
    d["region"] = d.annotation_id.str.split(":").str[2]
    d = d[~d.skipped.astype(str).eq("True")]
    d = d[["rater", "blind_code", "region", *COLS]].rename(columns=COLS)
    d.sort_values(["blind_code", "region", "rater"]).to_csv(OUT / "responses.csv", index=False)

    schema = json.load(open(SRC / "project.json"))["survey_schema"]
    json.dump(schema, open(OUT / "survey_schema.json", "w"), indent=1)

    reg = OUT / "regions"
    reg.mkdir(exist_ok=True)
    for f in sorted((SRC / "regions").glob("*.geojson")):
        shutil.copy(f, reg / f.name)
    print(f"{len(key)} regions, {len(d)} responses, {len(list(reg.glob('*.geojson')))} slides -> {OUT}")


if __name__ == "__main__":
    main()
