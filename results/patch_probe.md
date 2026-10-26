# Patch-level linear probe on the official split (MSI-CRC task)

Tile-level logistic regression, TRAIN tiles -> TEST tiles (93,408 / 99,904). Δ accuracy vs released tiles with 95% patient-cluster bootstrap CI.

| Encoder | Arm | Tile acc | Balanced acc | Tile AUC | Patient AUC | Δ acc vs released [95% CI] |
|---|---|---|---|---|---|---|
| Phikon | released | 0.731 | 0.675 | 0.730 | 0.827 |  |
| Phikon | raw114 | 0.735 | 0.648 | 0.724 | 0.765 | +0.004 [-0.029, +0.034] |
| Phikon | true05 | 0.713 | 0.611 | 0.699 | 0.765 | -0.018 [-0.072, +0.036] |
| Phikon-v2 | released | 0.736 | 0.672 | 0.738 | 0.788 |  |
| Phikon-v2 | raw114 | 0.739 | 0.667 | 0.752 | 0.808 | +0.003 [-0.023, +0.030] |
| Phikon-v2 | true05 | 0.710 | 0.614 | 0.682 | 0.764 | -0.026 [-0.068, +0.014] |
| UNI | released | 0.741 | 0.689 | 0.749 | 0.857 |  |
| UNI | raw114 | 0.735 | 0.697 | 0.761 | 0.846 | -0.006 [-0.039, +0.031] |
| UNI | true05 | 0.728 | 0.697 | 0.757 | 0.833 | -0.013 [-0.057, +0.036] |
| UNI2-h | released | 0.772 | 0.716 | 0.778 | 0.898 |  |
| UNI2-h | raw114 | 0.772 | 0.718 | 0.779 | 0.889 | +0.001 [-0.019, +0.027] |
| UNI2-h | true05 | 0.754 | 0.701 | 0.766 | 0.836 | -0.017 [-0.055, +0.029] |
| Virchow2 | released | 0.786 | 0.719 | 0.782 | 0.874 |  |
| Virchow2 | raw114 | 0.801 | 0.736 | 0.806 | 0.913 | +0.015 [-0.002, +0.035] |
| Virchow2 | true05 | 0.757 | 0.705 | 0.770 | 0.852 | -0.029 [-0.075, +0.007] |
| Prov-GigaPath | released | 0.750 | 0.691 | 0.761 | 0.851 |  |
| Prov-GigaPath | raw114 | 0.761 | 0.713 | 0.774 | 0.851 | +0.011 [-0.020, +0.053] |
| Prov-GigaPath | true05 | 0.743 | 0.699 | 0.773 | 0.847 | -0.007 [-0.042, +0.041] |
| Prov-GigaPath-flash | released | 0.722 | 0.683 | 0.748 | 0.872 |  |
| Prov-GigaPath-flash | true05 | 0.720 | 0.685 | 0.750 | 0.845 | -0.002 [-0.039, +0.039] |

Encoder ranking by tile accuracy, Spearman between arms (n=6 encoders): released vs true05 0.94; released vs raw114 0.83.

Ranking (released): Virchow2 > UNI2-h > Prov-GigaPath > UNI > Phikon-v2 > Phikon
Ranking (true05):   Virchow2 > UNI2-h > Prov-GigaPath > UNI > Phikon > Phikon-v2
