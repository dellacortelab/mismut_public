# External validation: TCGA TRAIN -> TCGA TEST and CPTAC-COAD (Prov-GigaPath-flash)

CPTAC-COAD: 62 patients, 14 MSI-H (PCR). Models fit on TCGA TRAIN only.

| Model | TCGA TEST AUC [95% CI] | CPTAC AUC [95% CI] |
|---|---|---|
| tumor_lr_mean | 0.874 [0.775, 0.956] | 0.823 [0.671, 0.938] |
| tumor_abmil | 0.915 [0.823, 0.988] | 0.935 [0.842, 0.994] |
| kather_lr_mean | 0.850 [0.735, 0.945] | 0.844 [0.689, 0.963] |
| kather_lr_max | 0.787 [0.676, 0.886] | 0.607 [0.457, 0.753] |
| kather_lr_mean+max | 0.807 [0.694, 0.902] | 0.616 [0.455, 0.771] |
| kather_abmil | 0.905 [0.814, 0.978] | 0.856 [0.738, 0.951] |
| iel_density | 0.793 [0.663, 0.901] | 0.732 [0.548, 0.893] |
| iel_per_100_epi | 0.865 [0.763, 0.947] | 0.677 [0.500, 0.835] |
| cells_lr | 0.790 [0.675, 0.890] | 0.653 [0.475, 0.818] |

Paired CPTAC differences:

| Comparison | Δ AUC [95% CI] |
|---|---|
| tumor_abmil - kather_lr_max | +0.327 [+0.156, +0.500] |
| tumor_abmil - kather_lr_mean+max | +0.318 [+0.138, +0.500] |
| tumor_abmil - kather_abmil | +0.079 [-0.013, +0.189] |
| tumor_abmil - tumor_lr_mean | +0.112 [-0.004, +0.262] |
| tumor_abmil - iel_density | +0.202 [+0.036, +0.399] |
| tumor_abmil - cells_lr | +0.281 [+0.113, +0.467] |
| kather_lr_mean - kather_lr_max | +0.237 [+0.058, +0.420] |
| iel_density - kather_lr_max | +0.125 [-0.097, +0.332] |

TCGA whole-tumor nucleus features, TRAIN means (MSIMUT vs MSS):

- neutrophil: 2.15 vs 1.21
- epithelial: 59.26 vs 54.15
- lymphocyte: 16.86 vs 13.30
- plasma: 2.50 vs 2.90
- eosinophil: 0.45 vs 0.48
- connective: 17.46 vs 18.84
- iel: 3.02 vs 1.51
- nuclei: 98.68 vs 90.88
- iel_per_100_epi: 5.28 vs 2.77
