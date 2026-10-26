# Established biomarker (IEL) vs foundation model on the same tumor tiles

Prov-GigaPath-flash embeddings and KongNet nuclei on identical tiles (≤1,000 tumor tiles per patient, 0.5 µm/px). Fit on TCGA TRAIN only.

| Model | TCGA TEST | CPTAC | SurGen |
|---|---|---|---|
| iel | 0.793 [0.66, 0.90] | 0.732 [0.55, 0.89] | 0.624 [0.49, 0.75] |
| iel_ratio | 0.865 [0.76, 0.95] | 0.677 [0.50, 0.83] | 0.632 [0.51, 0.74] |
| cells_lr | 0.790 [0.68, 0.89] | 0.653 [0.47, 0.82] | 0.680 [0.56, 0.80] |
| fm_mean | 0.869 [0.77, 0.95] | 0.818 [0.66, 0.94] | 0.816 [0.72, 0.90] |
| fm_abmil | 0.919 [0.83, 0.98] | 0.951 [0.87, 1.00] | 0.880 [0.79, 0.96] |
| late_fusion | 0.884 [0.78, 0.97] | 0.844 [0.68, 0.97] | 0.825 [0.72, 0.91] |
| early_fusion | 0.869 [0.77, 0.95] | 0.832 [0.67, 0.95] | 0.822 [0.72, 0.91] |
| iel_weighted | 0.864 [0.77, 0.95] | 0.829 [0.67, 0.96] | 0.827 [0.74, 0.91] |

## How much of the FM signal is IEL?

| cohort    |   rho_fm_iel |   rho_fm_iel_MSS |   auc_fm_residual_after_iel |   attention_iel_rho |   iel_MSS_mean |   iel_MSI_mean |   epi_MSS_mean |
|:----------|-------------:|-----------------:|----------------------------:|--------------------:|---------------:|---------------:|---------------:|
| tcga_test |        0.341 |            0.09  |                       0.811 |               0.101 |          1.289 |          2.631 |         54.047 |
| cptac     |        0.505 |            0.388 |                       0.708 |               0.132 |          1.723 |          3.012 |         50.907 |
| surgen    |        0.164 |           -0.02  |                       0.769 |               0.181 |          1.546 |          2.238 |         44.393 |
