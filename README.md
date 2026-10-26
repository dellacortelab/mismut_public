# mismut: from prediction to morphology for MSI in colorectal cancer

Code and results for *"From prediction to morphology: a deconfounded recipe for H&E biomarker discovery, applied to
microsatellite instability"* (Della Corte, Hopson, Szymanski, Perry; manuscript in preparation). Data deposit: Zenodo
[10.5281/zenodo.23666864](https://doi.org/10.5281/zenodo.23666864). Code, all versions: Zenodo
[10.5281/zenodo.23667396](https://doi.org/10.5281/zenodo.23667396) (v1.0:
[10.5281/zenodo.23667397](https://doi.org/10.5281/zenodo.23667397)).

## What this repository shows

The tile dataset released with Kather et al. 2019 (Nat. Med.; Zenodo [2530835](https://doi.org/10.5281/zenodo.2530835))
is a common benchmark for predicting microsatellite instability (MSI) from H&E. We audited it across six pathology
foundation models and four aggregation methods.

1. **Patch-count shortcut.** The training split was balanced by discarding tiles of microsatellite-stable patients.
   Max-containing pooling learns the resulting tile count:
   - it inflates random cross-validation AUC by up to 0.26
   - count is linearly decodable from max-pooled embeddings (R² 0.87 to 0.93)
   - models trained this way collapse to chance on external cohorts (CPTAC-COAD, SurGen)

   Mean pooling and attention MIL are robust.
2. **Resolution error.** The tiles are 256 µm regions at 224 px, i.e. **1.143 µm/px, not the documented 0.5 µm/px**.
   We registered all 192,312 tiles to their 364 source TCGA slides. The coordinates are in
   `results/tile_slide_registration.csv.gz` and on Zenodo
   ([10.5281/zenodo.23666864](https://doi.org/10.5281/zenodo.23666864)).
3. **Where the signal is.** On whole slides the MSI signal is intratumoral. Whole-tumor models transfer to CPTAC-COAD
   (PCR MSI-H) and SurGen (IHC dMMR) with AUC 0.86 to 0.93.
4. **Label errors.** 14 released labels contradict TCGA PanCancer Atlas molecular data (hypermutated or MSIsensor-high
   patients labelled MSS, and MSI-labelled patients without MSI or hypermutation). Excluding the 5 in the test split
   raises attention-MIL test AUC from 0.919 to 0.973 (`results/metrics/label_audit.csv`).
5. **Interpretable comparator and clinical operating points.** Intraepithelial lymphocyte density on identical tumor
   tiles reaches AUC 0.79, 0.73 and 0.62 (TCGA, CPTAC, SurGen) against 0.92, 0.95 and 0.88 for attention MIL.
   Thresholds fixed on training out-of-fold predictions; external calibration under-predicts MSI
   (`results/operating_points.md`).
6. **Pathologist review of attended regions.** Two blinded pathologists described 180 high- and low-attention
   regions on 20 test patients. In this exploratory analysis neutrophil presence separated MSI from MSS patients
   (AUC 0.91); attention-gated cell counts could not validate it because the gate makes every cell type
   discriminative (`results/review/summary.md`, `results/metrics/gate_controls.csv`).

Findings with numbers: [`results/README.md`](results/README.md). Tables: [`results/tables.md`](results/tables.md),
[`results/external_store.md`](results/external_store.md).

## Setup

```bash
python3.12 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/huggingface-cli login          # gated encoders need approved access on Hugging Face
# nucleus detection (KongNet via TIAToolbox) uses a separate environment:
python3.12 -m venv .venv-tia && .venv-tia/bin/pip install tiatoolbox pyarrow wsidicom imagecodecs
```

Large data go under `data/` (git-ignored). Override storage locations with environment variables:
- `MSI_ROOT`: project root for a cohort, e.g. a separate root for the STAD tiles
- `MSI_SLIDES`: TCGA slide storage
- `MSI_WSI`, `MSI_WSI_CPTAC`, `MSI_WSI_SURGEN`: whole-slide embedding output
- `MSI_TILES`: tumor tile store
- `MSI_SURGEN`: SurGen CZI files
- `MSI_CPTAC`: CPTAC DICOM slides

## Pipeline

| Step | Script | What it does |
|---|---|---|
| Data | `src/download.py` | Zenodo 2530835 tiles (CRC, STAD), NCT-CRC-HE-100K; parallel range requests, md5-checked |
| Index | `src/build_index.py` | Patient-level index, verified against published counts |
| Embeddings | `src/extract.py`, `src/encoders.py` | Phikon, Phikon-v2, UNI, UNI2-h, Virchow2, Prov-GigaPath, Prov-GigaPath-flash; pinned revisions |
| Protocols | `src/protocols.py`, `src/pool.py` | Random CV, published split with bootstrap, fixed-N, count-matched |
| Learners | `src/eval_lr.py`, `src/abmil.py` | Logistic regression on mean/max/mean+max pooling; gated ABMIL |
| Count probes | `src/confound.py` | Count-only models, ridge decoding of patch count |
| Resolution | `src/check_resolution.py`, `src/retile.py` | Register tiles to TCGA slides; re-cut regions at 0.5 µm/px |
| Whole slides | `src/wsi_regions.py`, `src/wsi_eval.py` | Tile all tissue; tumor map; tumor / margin / distant regions |
| Nuclei | `src/cells.py`, `src/wsi_cells.py` | KongNet nucleus classes, intraepithelial lymphocytes |
| External | `src/cptac.py`, `src/surgen.py`, `src/tumor_tiles.py`, `src/external.py`, `src/external_fm.py` | CPTAC-COAD (IDC), SurGen (BioImage Archive); TCGA-trained models applied unchanged |
| Report | `src/report.py`, `src/figs_extra.py` | Tables and figures |

Later analyses (run after the core pipeline):

| Step | Script | What it does |
|---|---|---|
| Identical tumor tiles | `src/iel_fm.py`, `src/iel_fm_max.py`, `src/surgen_cells.py` | IEL density vs FM attention MIL, and max pooling with released bag sizes, on the same tiles in three cohorts |
| Attention | `src/attention_enrichment.py`, `src/fig_attention.py`, `src/bagsize_external.py` | Cell-class enrichment under attention; bag-size checks on external cohorts |
| Patch probe | `src/patch_probe.py` | Tile-level probe of the slide-level signal |
| Operating points | `src/operating_points.py` | Out-of-fold thresholds on TRAIN (rule-out, Youden), calibration, three cohorts |
| Label audit | `src/fetch_molecular.py`, `src/label_audit.py` | cBioPortal TMB / MANTIS / MSIsensor / POLE vs released labels |
| Review regions | `src/review_regions.py`, `src/overview_crops.py` | Select high- and low-attention regions for blinded review |
| Review analysis | `src/review_analysis.py`, `src/gate_controls.py` | Unblind and analyse responses (Fig. 4); attention-gate controls with permuted labels |
| Deposit | `src/zenodo_bundle.py` | Anonymised review files for Zenodo |

To reproduce the review analysis from the deposit: unzip `review.zip` from Zenodo, then
`python src/review_analysis.py from-deposit <dir>`, `python src/fetch_molecular.py` and
`python src/review_analysis.py analyze data/review/export_deposit.csv`.

`run_all.sh` runs the core TCGA pipeline (download, index, fp16 check, extraction, evaluation, report).

## Data sources

All public:
- TCGA tiles: Zenodo 2530835.
- TCGA diagnostic slides: [GDC](https://portal.gdc.cancer.gov/).
- CPTAC-COAD: [NCI Imaging Data Commons](https://portal.imaging.datacommons.cancer.gov/), with MSI calls from the
  CPTAC clinical data (Vasaikar et al., Cell 2019).
- SurGen SR386: BioImage Archive [S-BIAD1285](https://doi.org/10.6019/S-BIAD1285).
- NCT-CRC-HE-100K: Zenodo [1214456](https://doi.org/10.5281/zenodo.1214456).

## Not included

One secondary analysis in the manuscript (a zero-shot tumor map from a proprietary prostate cancer detector) uses a
model that is not publicly available. Its code is not part of this repository; its outputs are
(`results/p1_*.md`, `results/metrics/p1_*.csv`). All other analyses run from this code and public data.

## Citation

Please cite the manuscript (reference to follow) and the original dataset (Kather et al., Nat. Med. 25, 1054–1056,
2019).
