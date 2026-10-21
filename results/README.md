# Findings

Runs of 2026-10-16/17 on one RTX A6000. 360 patients (TRAIN 39 MSIMUT / 221 MSS, TEST 26 / 74), 6 foundation models
plus the authors' archived Phikon embeddings, 4 aggregations, and 3 interpretable feature sets. ROC-AUC unless
stated. E1 = random 5-fold CV, 10 repeats, mean ± SD over repeats. E2 = fit on TRAIN, test once on TEST, 95% CI
from 5,000 stratified patient bootstraps; paired comparisons reuse identical resamples. Full numbers:
`tables.md`, `summary.csv`, `paired_e2_bootstrap.csv`, `feature_univariate.csv`.

Main tables use the Zenodo tiles as published (as the target paper did); E7 repeats everything on regions re-cut
from the original slides at the correct resolution. The source TEST set has only 26 positives, so single E2 AUCs carry CIs about ±0.10 wide. Paired differences are much
tighter and are what the claims below rest on.

## Sanity gate

Our Phikon reproduces the paper (E1 paper folds mean+max 0.916 vs 0.915; E2 0.797 vs 0.788; count-only E2 0.517 vs
0.517; max-pool count decoding R² 0.872 vs 0.874). The authors' archived embeddings give the same numbers within
0.01 (patient-mean cosine to ours ≥ 0.988).

## H1. Max-containing pooling inflates random-CV AUC. Supported for every encoder.

Random-CV minus source-split AUC is +0.13 to +0.26 for max pooling and +0.10 to +0.20 for mean+max, against −0.01 to
+0.03 for mean pooling. Under random CV, max or mean+max looks best for every encoder (0.91 to 0.95). On the source
split max pooling is the worst aggregation for every encoder (0.69 to 0.79), significantly below mean pooling for UNI (−0.16,
95% CI −0.29 to −0.04), UNI2-h, Virchow2 and GigaPath (−0.10 to −0.12, CIs exclude 0). Count-matching TRAIN (count-only
AUC within TRAIN drops from 0.92 to 0.500) raises max-pooling source-split AUC by +0.08 to +0.19 (CI excludes 0 for
Phikon-v2, UNI, UNI2-h, Virchow2, GigaPath; Phikon +0.08, CI −0.01 to 0.17). Matching leaves mean pooling unchanged
(|Δ| ≤ 0.01) and ABMIL unchanged (|Δ| ≤ 0.02). Fixed-N subsampling (Fig. 3) closes the gap the same way: at equal
patch counts max and mean pooling give similar source-split AUCs. The mechanism is the one proposed in the target
paper: Kather undersampled MSS tiles in TRAIN only, so in TRAIN MSIMUT patients have ~5.7x more tiles than MSS
patients (1,198 vs 211), and max pooling encodes tile count. Patch count alone reaches 0.82 under random CV and 0.52 on
TEST. Appending log1p(count) to any representation changes source-split AUC by ≤ 0.01.

## H2. Patch count is decodable from max-pooled embeddings. Supported for every encoder.

Out-of-fold ridge R² for log1p(count): 0.87 to 0.93 from max pooling, 0.22 to 0.32 from mean pooling (Fig. 2).
Within MSS patients only, where count cannot be mediated by the label, the numbers are the same (max 0.86 to 0.92,
mean 0.26 to 0.30). UNI, UNI2-h and GigaPath encode count most strongly (R² ≈ 0.93); UNI also shows the largest max-pooling
inflation (+0.26). Out-of-sample max-pooling scores correlate with count within each class
(Spearman 0.57 to 0.69 within MSS on TEST); mean pooling and ABMIL scores do not (|ρ| ≤ 0.31, mostly ≤ 0.15).

## H3. ABMIL is less sensitive than max pooling. Supported, more strongly than expected.

For every encoder, ABMIL shows no inflation: its random-CV AUC is 0.02 to 0.07 below its source-split AUC.
Count-matching does not change it, and its scores do not track count within class. ABMIL gives the best
source-split AUC for every encoder (0.86 to 0.93), above mean pooling by +0.02 to +0.07. The paired CI excludes 0 for
both Phikon embeddings, Phikon-v2, UNI and UNI2-h. We did not test why ABMIL scores lower under random CV than on the source split.
One candidate is that the pooled cohort mixes TRAIN MSS patients with few tiles (median ~200) and TEST patients with
many; this is untested.

## H4. The source-robust AUC is encoder-dependent. Supported, but higher than the 0.75 to 0.85 we expected.

With count-robust aggregation, source-split AUC ranges from 0.79 (Phikon, mean) to 0.93 (UNI2-h, ABMIL). Mean pooling:
Phikon 0.79, Phikon-v2 0.82, UNI 0.85, Virchow2 0.88, UNI2-h 0.89, GigaPath 0.89. ABMIL: 0.86, 0.88, 0.89, 0.91,
0.93, 0.91 in the same order. The ranking roughly follows model generation, but the TEST CIs overlap heavily (±0.10).
So the main statements are relative: max-containing pooling under random CV overstates performance on this dataset,
and mean pooling or ABMIL on the source split is the honest estimate. The paper's headline 0.915 (random CV, Phikon
mean+max) corresponds to 0.79 on the source split for that model, while current encoders with ABMIL reach about 0.9
on the source split.

## E6. Interpretable count-invariant features (added)

Per tile: 9 tissue-type probabilities (UNI linear probe trained on NCT-CRC-HE-100K; 95.2% accuracy, 93.7% balanced,
on CRC-VAL-HE-7K) and nucleus counts by class from KongNet (Lizard/CoNIC colon), plus an intraepithelial-lymphocyte
proxy (lymphocytes within 10 µm of an epithelial nucleus). Mean pooling gives per-patient tissue fractions and
nucleus densities. Count decoding from these is R² ≈ 0.00 to 0.05 by construction, but max pooling over the same
17 features again decodes count (R² 0.61) and is inflated (+0.08).

- Tissue + nuclei, mean pooling: source split 0.840 [0.751, 0.917], random CV 0.808 ± 0.018, count-matched 0.854.
  That is +0.05 over Phikon mean pooling (CI −0.05 to 0.16), not different from UNI/UNI2-h/Virchow2/GigaPath mean
  pooling (Δ −0.01 to −0.05, CIs include 0), and below UNI2-h ABMIL (−0.09, CI −0.15 to −0.03).
- The intraepithelial-lymphocyte proxy alone has AUC 0.74 in TRAIN and 0.81 in TEST, the strongest single feature
  and the largest LR weight (MSIMUT 16.7 vs MSS 8.4 per tile in TRAIN). Neutrophils (0.69 TEST), lymphocytes (0.66),
  mucin fraction (0.62) and stroma fraction (0.62) follow. These match the known MSI morphology (tumor-infiltrating
  lymphocytes, mucinous differentiation, inflammatory infiltrate).
- Tissue fractions alone are weak (0.62). On tumor-only tiles there is little tissue variation left to exploit.
- So 17 count-invariant morphology numbers recover most of the source-split signal of mean-pooled foundation
  embeddings, but ABMIL on the strongest encoders remains better.

## E7. Tile resolution: documented 0.5 µm/px, actually 1.143 µm/px; impact is encoder-dependent

Zenodo 2530835 documents the tiles as 0.5 µm/px. They are 256 µm regions stored at 224 px (1.143 µm/px), as
Kather et al. 2019 state in their Methods ("256 μm edge length"). We registered all 192,312 tiles to their 364
source slides (34 TCGA tissue source sites, 7 scanner pixel sizes from 0.232 to 0.263 µm/px; `src/retile.py`).
Median NCC was 0.956; 120 tiles (0.06%) were below 0.5. In 1,092 per-slide scale sweeps, 99.45% peaked at exactly
1.143 µm/px. The foundation models above (and in the target paper) therefore saw tissue at ~0.44x their 0.5 µm/px
training resolution.

To measure the cost, we re-read the same 256 µm regions from the slides (identical tiles and per-patient counts) in
two arms. The first is 0.5 µm/px: four 224 px crops per region, embeddings averaged. The second is 1.14 µm/px with
raw colour, which separates the effect of Kather's Macenko normalization from resolution. Source-split AUC, 0.5
µm/px re-cut minus Zenodo tiles (paired bootstrap 95% CI, `paired_e2_resolution.csv`, Fig. 5):

| Encoder | Mean pooling | ABMIL |
|---|---|---|
| Phikon | 0.789 → 0.857, +0.068 [+0.007, +0.138] | 0.855 → 0.914, +0.059 [+0.009, +0.119] |
| Phikon-v2 | 0.821 → 0.822, +0.001 [−0.084, +0.091] | 0.879 → 0.832, −0.047 [−0.122, +0.018] |
| UNI | 0.851 → 0.871, +0.020 [−0.066, +0.095] | 0.886 → 0.906, +0.020 [−0.011, +0.056] |
| UNI2-h | 0.885 → 0.884, −0.001 [−0.042, +0.035] | 0.929 → 0.927, −0.002 [−0.018, +0.013] |
| Virchow2 | 0.877 → 0.879, +0.002 [−0.018, +0.023] | 0.906 → 0.917, +0.010 [−0.014, +0.042] |
| Prov-GigaPath | 0.891 → 0.889, −0.002 [−0.037, +0.035] | 0.908 → 0.932, +0.025 [−0.003, +0.062] |

- **Large current encoders are robust to the mismatch.** For UNI2-h, Virchow2 and GigaPath with mean pooling or
  ABMIL, |Δ| ≤ 0.03 and every CI includes 0.
- **Phikon, the smallest encoder and the one used in the target paper, was handicapped.** It gains +0.06 to +0.07
  at the correct resolution (CIs exclude 0). With ABMIL it reaches 0.914, level with the largest models. Most of
  that gain is resolution, not colour: the 1.14 µm/px raw-colour arm gives +0.03 (mean) and +0.02 (ABMIL), CIs
  including 0.
- **Max pooling improves somewhat on the source split at 0.5 µm/px** (UNI +0.09, Virchow2 +0.05, CIs exclude 0).
- **The confound conclusions do not change.** Max-pooling random-CV AUC stays at 0.90 to 0.95, and patch count
  remains decodable from max pooling (R² 0.82 to 0.92) versus mean pooling (0.25 to 0.36).
- **One unexplained drop:** Phikon-v2 ABMIL is lower on raw-colour re-cuts (1.14 µm/px: −0.086, CI −0.165 to
  −0.023). It is one of ~60 resolution comparisons with no multiplicity correction; we do not interpret it.

## E8. Where is the MSI signal on the whole slide? Mostly in the tumor (margin hypothesis not supported)

Encoder: Prov-GigaPath-flash (ViT-S/16), chosen for speed (4,266 tiles/s). Calibration on
the Kather regions re-cut at 0.5 µm/px: ABMIL source split 0.890 (UNI2-h 0.927, Virchow2 0.917, GigaPath 0.932,
UNI 0.906, Phikon 0.914); LR mean 0.850. Its max pooling shows the same shortcut (random CV 0.92, source split 0.79).

All 364 slides were tiled at 0.5 µm/px (3.9 M tissue tiles, `src/wsi_regions.py`). The tumor map comes from an
in-domain, cross-fitted tumor classifier trained on Kather's own tumor decisions on fully covered slides. Kather
coverage could not be used directly: it is sparse on the undersampled MSS TRAIN slides, which would leak the label.
The NCT probe also failed on raw-colour tiles. On held-out complete slides, 97 to 99% of Kather tumor tiles fall
inside the map. Regions: tumor; margin (non-tumor within 1 mm); distant. Cohort: 227 patients with ≥ 20 tiles in
every region (TRAIN 168 with 21 MSIMUT, TEST 59 with 19 MSIMUT); smaller slides drop out. `src/wsi_eval.py`.

Source-split ROC-AUC [95% CI] (random CV 10x5 in brackets after):

| Region | ABMIL | LR mean pooling |
|---|---|---|
| Tumor | 0.896 [0.774, 0.991] (0.847) | 0.838 [0.709, 0.943] (0.835) |
| Margin only | 0.650 [0.503, 0.787] (0.622) | 0.707 [0.542, 0.847] (0.631) |
| Distant only | 0.585 [0.437, 0.726] (0.452) | 0.553 (0.525) |
| Tumor + margin | 0.895 (0.852) | 0.782 (0.792) |
| Whole tissue | 0.883 (0.849) | 0.796 (0.804) |

- **The signal is in the tumor.** Adding the margin does not help ABMIL (tumor+margin − tumor −0.002, CI −0.017 to
  +0.013), and adding all tissue does not either (−0.013, CI −0.054 to +0.020). Both dilute mean pooling (−0.04 to
  −0.06).
- **Margin-only tissue carries a weak signal** (ABMIL 0.65, mean pooling 0.71; lower CI bounds at or near 0.5), well
  short of the pre-registered ≥ 0.75. Distant tissue is near chance.
- **Whole-slide tiling removes the patch-count shortcut.** Count-only AUC is 0.43 to 0.50 in random CV, and random
  CV no longer exceeds the source split.
- **Caveats:**
  - one (small) encoder
  - reduced cohort (19 TEST positives)
  - TCGA diagnostic sections do not necessarily include the deep invasive front, so "margin" here is
    peritumoral tissue in the section, not a curated invasive margin
  - the NCT tissue-fraction features are unreliable on raw-colour tiles (domain shift) and are not interpreted

Per the pre-registered decision rule, the large encoders were not run on this question.

## E9/E10. Whole-tumor nucleus features and external validation on CPTAC-COAD

All models were fit on TCGA TRAIN only (the published split) and applied unchanged to TCGA TEST and to CPTAC-COAD.
CPTAC-COAD: IDC DICOM, 62 patients with PCR MSI calls, 14 MSI-H, 110 tumor H&E slides; tiled and tumor-mapped
exactly like TCGA, with the tumor classifier trained on TCGA only. Encoder: Prov-GigaPath-flash. Nucleus features come
from KongNet at true 0.5 µm/px on up to 1,000 tumor tiles per patient. Code: `src/wsi_cells.py`, `src/cptac.py`,
`src/external.py`; full table in `external_flash.md`.

| Model | TCGA TEST | CPTAC (external) |
|---|---|---|
| Whole-tumor ABMIL | 0.915 [0.823, 0.988] | **0.935 [0.842, 0.994]** |
| Whole-tumor LR mean | 0.874 | 0.823 |
| Kather-tile ABMIL | 0.905 | 0.856 |
| Kather-tile LR mean | 0.850 | 0.844 |
| Kather-tile LR max | 0.787 | **0.607 [0.457, 0.753]** |
| Kather-tile LR mean+max | 0.807 | 0.616 |
| IEL density (no fitting) | 0.793 | 0.732 |
| IEL per 100 epithelial nuclei | 0.865 | 0.677 |
| Nucleus-feature LR | 0.790 | 0.653 |

- **The shortcut does not generalize.** Max-containing pooling trained on the count-confounded Kather tiles falls
  to near chance on CPTAC (0.61). Mean pooling on the same tiles holds (0.84; mean − max +0.24, paired CI +0.06 to
  +0.42).
- **Whole-tumor ABMIL generalizes best** (0.935 on CPTAC). It is above Kather-tile ABMIL (+0.08, CI −0.01 to +0.19)
  and well above every interpretable nucleus model (vs IEL density +0.20, CI +0.04 to +0.40; vs nucleus LR +0.28).
- **The biology is clear but the simple features are not robust.** At true resolution, whole-tumor IEL density is
  twice as high in MSIMUT (3.0 vs 1.5 per tile; 5.3 vs 2.8 per 100 epithelial nuclei in TRAIN). As a classifier it
  drops on the external site (IEL ratio 0.87 → 0.68), likely from detector calibration shifts across scanners and
  stains. On this evidence, the foundation-model route beats the handcrafted-feature route for MSI.
- **Caveats:** 14 external positives (wide CIs); one encoder; CPTAC labels are PCR MSI-H versus TCGA MSIMUT.


## E12. Second external cohort: SurGen SR386 (Scotland, Zeiss CZI, IHC-defined dMMR)

128 cases: all 32 dMMR plus 96 random pMMR from SR386, downloaded via Globus. They were processed with the
identical pipeline: tissue plan, Edge embedding, TCGA-trained tumor map, lossless tile store (≤ 4,000 tumor tiles per
case, 472k tiles), and three encoders. Models were fit on TCGA TRAIN only (`external_store.md`, Fig. 7).

| Model | Edge | UNI2-h | Virchow2 |
|---|---|---|---|
| Whole-tumor ABMIL | 0.862 | 0.899 | 0.895 |
| Whole-tumor mean | 0.809 | 0.919 | 0.810 |
| Kather-tile ABMIL | 0.861 | 0.902 | 0.867 |
| Kather-tile mean | 0.809 | 0.900 | 0.826 |
| **Kather-tile max** | **0.506** | **0.469** | **0.497** |
| Kather-tile mean+max | 0.573 | 0.573 | 0.647 |

- **The shortcut collapse replicates at a second independent site.** Max pooling falls to chance with every
  encoder and assigns P(dMMR) > 0.99 to 99% of cases. Mean minus max on the same Kather tiles is +0.30 to +0.43,
  with paired CIs excluding 0.
- **Honest aggregation transfers:** ABMIL 0.86 to 0.90 and mean pooling 0.81 to 0.92. Whole-tumor and Kather-tile
  ABMIL are equivalent here (|Δ| ≤ 0.03, CIs include 0).
- **Labels differ across cohorts** (TCGA MSIMUT, CPTAC PCR MSI-H, SurGen IHC dMMR). The consistent transfer
  suggests the morphology learned on TCGA is the MSI/dMMR phenotype rather than a TCGA-specific signal.

## Other limitations

- One dataset, one TEST set of 100 patients (26 positive). No external cohort.
- Tumor-only tiles: the invasive margin, including the Crohn-like lymphoid reaction, is not in the data.
- The label is MSIMUT (MSI or hypermutated, from the TCGA pan-cancer calls), not clinical MSI-PCR or MMR IHC.
- E6 nucleus classes and the IEL proxy come from a detector applied out of domain (Macenko-normalized tiles,
  upsampled 2x to ~0.57 µm/px). Their outputs were checked against NCT tissue classes and visually, not against nucleus annotations.
- ABMIL hyperparameters were fixed in advance (not tuned); a tuned ABMIL might score higher.
