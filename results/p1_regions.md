# MSI with a Kather-independent tumor map (P1-flash, prostate-trained, zero-shot)

Tumor = 3x3-smoothed mean P1 >= 0.5 (fixed a priori). Models fit on TCGA TRAIN only.

Map agreement: tcga_p1_frac 0.429, tcga_kather_frac 0.629, tcga_overlap_jaccard 0.653, cptac_p1_frac 0.396, cptac_kather_frac 0.567, cptac_overlap_jaccard 0.532

| Map | Model | Cohort | n (pos) | AUC [95% CI] |
|---|---|---|---|---|
| kather_tumor | lr_mean | tcga_test | 100 (26) | 0.874 [0.775, 0.956] |
| kather_tumor | lr_mean | cptac | 62 (14) | 0.823 [0.671, 0.938] |
| kather_tumor | abmil | tcga_test | 100 (26) | 0.915 [0.823, 0.988] |
| kather_tumor | abmil | cptac | 62 (14) | 0.935 [0.842, 0.994] |
| p1_tumor | lr_mean | tcga_test | 100 (26) | 0.876 [0.774, 0.962] |
| p1_tumor | lr_mean | cptac | 62 (14) | 0.807 [0.676, 0.917] |
| p1_tumor | abmil | tcga_test | 100 (26) | 0.910 [0.819, 0.982] |
| p1_tumor | abmil | cptac | 62 (14) | 0.927 [0.847, 0.987] |

P1 map minus Kather-derived map (paired, common patients):

| Model | Cohort | n | Δ AUC [95% CI] |
|---|---|---|---|
| lr_mean | tcga_test | 100 | +0.002 [-0.025, +0.033] |
| lr_mean | cptac | 62 | -0.016 [-0.085, +0.054] |
| abmil | tcga_test | 100 | -0.006 [-0.029, +0.014] |
| abmil | cptac | 62 | -0.007 [-0.065, +0.051] |
