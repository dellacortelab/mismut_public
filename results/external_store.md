# Whole-tumor vs Kather-tile models across encoders: TCGA TRAIN -> TCGA TEST / CPTAC-COAD

Identical tile store for all encoders (<= 4,000 tumor tiles per patient, 0.5 µm/px, PNG). Models fit on TCGA TRAIN only. ROC-AUC [95% CI, stratified bootstrap].

| Encoder | Model | TCGA TEST | CPTAC |
|---|---|---|---|
| Prov-GigaPath-flash | tumor_lr_mean | 0.876 [0.779, 0.957] | 0.814 [0.661, 0.933] |
| Prov-GigaPath-flash | tumor_abmil | 0.903 [0.808, 0.979] | 0.929 [0.839, 0.990] |
| Prov-GigaPath-flash | kather_lr_mean | 0.850 [0.735, 0.945] | 0.815 [0.656, 0.941] |
| Prov-GigaPath-flash | kather_lr_max | 0.787 [0.676, 0.886] | 0.584 [0.451, 0.704] |
| Prov-GigaPath-flash | kather_lr_mean+max | 0.807 [0.694, 0.902] | 0.534 [0.380, 0.680] |
| Prov-GigaPath-flash | kather_abmil | 0.905 [0.814, 0.978] | 0.868 [0.746, 0.963] |
| UNI2-h | tumor_lr_mean | 0.899 [0.798, 0.980] | 0.885 [0.780, 0.964] |
| UNI2-h | tumor_abmil | 0.905 [0.799, 0.986] | 0.903 [0.815, 0.967] |
| UNI2-h | kather_lr_mean | 0.884 [0.772, 0.974] | 0.845 [0.708, 0.951] |
| UNI2-h | kather_lr_max | 0.800 [0.692, 0.893] | 0.531 [0.500, 0.562] |
| UNI2-h | kather_lr_mean+max | 0.850 [0.743, 0.936] | 0.531 [0.500, 0.562] |
| UNI2-h | kather_abmil | 0.921 [0.822, 0.994] | 0.903 [0.810, 0.973] |
| Virchow2 | tumor_lr_mean | 0.892 [0.797, 0.966] | 0.946 [0.866, 0.997] |
| Virchow2 | tumor_abmil | 0.915 [0.824, 0.983] | 0.914 [0.824, 0.976] |
| Virchow2 | kather_lr_mean | 0.878 [0.763, 0.972] | 0.940 [0.857, 0.994] |
| Virchow2 | kather_lr_max | 0.811 [0.703, 0.902] | 0.560 [0.443, 0.656] |
| Virchow2 | kather_lr_mean+max | 0.852 [0.748, 0.938] | 0.602 [0.471, 0.719] |
| Virchow2 | kather_abmil | 0.920 [0.831, 0.987] | 0.865 [0.753, 0.954] |

Paired CPTAC differences (same resamples):

| Encoder | Comparison | Δ AUC [95% CI] |
|---|---|---|
| Prov-GigaPath-flash | tumor_abmil - kather_lr_max | +0.344 [+0.193, +0.501] |
| Prov-GigaPath-flash | tumor_abmil - kather_abmil | +0.061 [-0.043, +0.177] |
| Prov-GigaPath-flash | kather_lr_mean - kather_lr_max | +0.231 [+0.082, +0.393] |
| Prov-GigaPath-flash | tumor_abmil - tumor_lr_mean | +0.115 [+0.010, +0.256] |
| UNI2-h | tumor_abmil - kather_lr_max | +0.372 [+0.283, +0.441] |
| UNI2-h | tumor_abmil - kather_abmil | +0.000 [-0.045, +0.049] |
| UNI2-h | kather_lr_mean - kather_lr_max | +0.314 [+0.177, +0.424] |
| UNI2-h | tumor_abmil - tumor_lr_mean | +0.018 [-0.049, +0.091] |
| Virchow2 | tumor_abmil - kather_lr_max | +0.353 [+0.214, +0.495] |
| Virchow2 | tumor_abmil - kather_abmil | +0.049 [+0.001, +0.115] |
| Virchow2 | kather_lr_mean - kather_lr_max | +0.380 [+0.248, +0.513] |
| Virchow2 | tumor_abmil - tumor_lr_mean | -0.033 [-0.082, +0.007] |

## SurGen SR386 subset (MMR loss by IHC), models fit on TCGA TRAIN only

| Encoder | Model | SurGen AUC [95% CI] |
|---|---|---|
| Prov-GigaPath-flash | tumor_lr_mean | 0.809 [0.708, 0.895] (n=128, 32 dMMR) |
| Prov-GigaPath-flash | tumor_abmil | 0.862 [0.765, 0.944] (n=128, 32 dMMR) |
| Prov-GigaPath-flash | kather_lr_mean | 0.809 [0.707, 0.895] (n=128, 32 dMMR) |
| Prov-GigaPath-flash | kather_lr_max | 0.506 [0.438, 0.565] (n=128, 32 dMMR) |
| Prov-GigaPath-flash | kather_lr_mean+max | 0.573 [0.477, 0.660] (n=128, 32 dMMR) |
| Prov-GigaPath-flash | kather_abmil | 0.861 [0.761, 0.944] (n=128, 32 dMMR) |
| UNI2-h | tumor_lr_mean | 0.919 [0.845, 0.972] (n=128, 32 dMMR) |
| UNI2-h | tumor_abmil | 0.899 [0.803, 0.976] (n=128, 32 dMMR) |
| UNI2-h | kather_lr_mean | 0.900 [0.823, 0.960] (n=128, 32 dMMR) |
| UNI2-h | kather_lr_max | 0.469 [0.406, 0.516] (n=128, 32 dMMR) |
| UNI2-h | kather_lr_mean+max | 0.573 [0.493, 0.645] (n=128, 32 dMMR) |
| UNI2-h | kather_abmil | 0.902 [0.810, 0.973] (n=128, 32 dMMR) |
| Virchow2 | tumor_lr_mean | 0.810 [0.727, 0.882] (n=128, 32 dMMR) |
| Virchow2 | tumor_abmil | 0.895 [0.824, 0.953] (n=128, 32 dMMR) |
| Virchow2 | kather_lr_mean | 0.826 [0.744, 0.897] (n=128, 32 dMMR) |
| Virchow2 | kather_lr_max | 0.497 [0.430, 0.554] (n=128, 32 dMMR) |
| Virchow2 | kather_lr_mean+max | 0.647 [0.553, 0.732] (n=128, 32 dMMR) |
| Virchow2 | kather_abmil | 0.867 [0.775, 0.943] (n=128, 32 dMMR) |
