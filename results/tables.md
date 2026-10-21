# Main results (ROC-AUC)

| Encoder | Aggregation | E1 random CV 10x5 | E2 source split [95% CI] | Inflation (E1-E2) | E3b matched E1 | E3b matched E2 | Fixed-N 100 E1 | Fixed-N 100 E2 |
|---|---|---|---|---|---|---|---|---|
| Phikon (archived) | Mean | 0.796 ± 0.015 | 0.793 [0.674, 0.899] | +0.004 | 0.805 | 0.784 | 0.783 | 0.797 |
| Phikon (archived) | Max | 0.912 ± 0.008 | 0.766 [0.654, 0.862] | +0.146 | 0.808 | 0.852 | 0.814 | 0.798 |
| Phikon (archived) | Mean+max | 0.910 ± 0.007 | 0.789 [0.681, 0.885] | +0.121 | 0.824 | 0.833 | 0.818 | 0.800 |
| Phikon (archived) | ABMIL | 0.816 ± 0.019 | 0.856 [0.766, 0.932] | -0.040 | 0.820 | 0.869 | 0.809 | 0.814 |
| Phikon | Mean | 0.795 ± 0.016 | 0.789 [0.666, 0.900] | +0.006 | 0.802 | 0.788 | 0.784 | 0.792 |
| Phikon | Max | 0.911 ± 0.010 | 0.777 [0.667, 0.872] | +0.134 | 0.813 | 0.854 | 0.814 | 0.783 |
| Phikon | Mean+max | 0.911 ± 0.006 | 0.797 [0.690, 0.890] | +0.114 | 0.828 | 0.839 | 0.818 | 0.797 |
| Phikon | ABMIL | 0.813 ± 0.009 | 0.855 [0.761, 0.932] | -0.041 | 0.831 | 0.872 | 0.810 | 0.818 |
| Phikon-v2 | Mean | 0.812 ± 0.016 | 0.821 [0.705, 0.921] | -0.009 | 0.810 | 0.819 | 0.798 | 0.826 |
| Phikon-v2 | Max | 0.929 ± 0.005 | 0.745 [0.630, 0.847] | +0.184 | 0.855 | 0.847 | 0.860 | 0.788 |
| Phikon-v2 | Mean+max | 0.923 ± 0.005 | 0.773 [0.663, 0.870] | +0.149 | 0.856 | 0.850 | 0.847 | 0.810 |
| Phikon-v2 | ABMIL | 0.810 ± 0.014 | 0.879 [0.784, 0.956] | -0.069 | 0.815 | 0.879 | 0.818 | 0.832 |
| UNI | Mean | 0.880 ± 0.017 | 0.851 [0.746, 0.937] | +0.029 | 0.891 | 0.851 | 0.869 | 0.821 |
| UNI | Max | 0.942 ± 0.006 | 0.687 [0.570, 0.798] | +0.255 | 0.877 | 0.876 | 0.868 | 0.825 |
| UNI | Mean+max | 0.945 ± 0.005 | 0.747 [0.635, 0.849] | +0.197 | 0.890 | 0.873 | 0.875 | 0.825 |
| UNI | ABMIL | 0.850 ± 0.016 | 0.886 [0.789, 0.965] | -0.037 | 0.852 | 0.880 | 0.847 | 0.838 |
| UNI2-h | Mean | 0.894 ± 0.010 | 0.885 [0.782, 0.966] | +0.009 | 0.900 | 0.887 | 0.878 | 0.874 |
| UNI2-h | Max | 0.952 ± 0.005 | 0.773 [0.659, 0.873] | +0.178 | 0.904 | 0.886 | 0.882 | 0.852 |
| UNI2-h | Mean+max | 0.945 ± 0.004 | 0.826 [0.719, 0.916] | +0.119 | 0.904 | 0.891 | 0.888 | 0.869 |
| UNI2-h | ABMIL | 0.886 ± 0.009 | 0.929 [0.848, 0.987] | -0.044 | 0.879 | 0.921 | 0.875 | 0.901 |
| Virchow2 | Mean | 0.887 ± 0.013 | 0.877 [0.767, 0.967] | +0.009 | 0.890 | 0.875 | 0.888 | 0.882 |
| Virchow2 | Max | 0.933 ± 0.004 | 0.756 [0.641, 0.861] | +0.177 | 0.893 | 0.873 | 0.890 | 0.844 |
| Virchow2 | Mean+max | 0.936 ± 0.005 | 0.833 [0.727, 0.920] | +0.104 | 0.902 | 0.882 | 0.896 | 0.857 |
| Virchow2 | ABMIL | 0.875 ± 0.014 | 0.906 [0.809, 0.982] | -0.031 | 0.884 | 0.906 | 0.875 | 0.856 |
| Prov-GigaPath | Mean | 0.901 ± 0.008 | 0.891 [0.800, 0.965] | +0.009 | 0.907 | 0.890 | 0.884 | 0.867 |
| Prov-GigaPath | Max | 0.951 ± 0.005 | 0.792 [0.686, 0.884] | +0.159 | 0.907 | 0.889 | 0.903 | 0.855 |
| Prov-GigaPath | Mean+max | 0.954 ± 0.004 | 0.838 [0.737, 0.921] | +0.116 | 0.917 | 0.896 | 0.902 | 0.868 |
| Prov-GigaPath | ABMIL | 0.884 ± 0.015 | 0.908 [0.819, 0.978] | -0.023 | 0.892 | 0.897 | 0.875 | 0.855 |
| Tissue fractions | Mean | 0.662 ± 0.017 | 0.624 [0.491, 0.750] | +0.038 | 0.673 | 0.616 | 0.669 | 0.570 |
| Tissue fractions | Max | 0.787 ± 0.010 | 0.695 [0.573, 0.807] | +0.092 | 0.726 | 0.700 | 0.706 | 0.622 |
| Tissue fractions | Mean+max | 0.805 ± 0.015 | 0.689 [0.568, 0.804] | +0.116 | 0.721 | 0.692 | 0.687 | 0.605 |
| Nucleus densities | Mean | 0.779 ± 0.009 | 0.840 [0.721, 0.937] | -0.061 | 0.778 | 0.848 | 0.758 | 0.828 |
| Nucleus densities | Max | 0.750 ± 0.012 | 0.668 [0.556, 0.780] | +0.081 | 0.655 | 0.620 | 0.683 | 0.666 |
| Nucleus densities | Mean+max | 0.800 ± 0.008 | 0.691 [0.557, 0.815] | +0.110 | 0.759 | 0.744 | 0.748 | 0.726 |
| Tissue + nuclei | Mean | 0.808 ± 0.018 | 0.840 [0.751, 0.917] | -0.032 | 0.817 | 0.854 | 0.812 | 0.818 |
| Tissue + nuclei | Max | 0.847 ± 0.010 | 0.767 [0.653, 0.870] | +0.081 | 0.786 | 0.788 | 0.782 | 0.766 |
| Tissue + nuclei | Mean+max | 0.893 ± 0.010 | 0.808 [0.717, 0.893] | +0.085 | 0.837 | 0.857 | 0.817 | 0.792 |
| log1p(count) only | LR | 0.819 ± 0.003 | 0.517 [0.382, 0.653] | | 0.558 | 0.500 | | |

# Count decoding (out-of-fold ridge, log1p(count))

| Encoder | Cohort | Mean R² | Max R² | Mean+max R² | Mean ρ | Max ρ |
|---|---|---|---|---|---|---|
| Phikon (archived) | all | 0.237 | 0.872 | 0.878 | 0.446 | 0.921 |
| Phikon (archived) | MSS_only | 0.262 | 0.861 | 0.848 | 0.453 | 0.909 |
| Phikon | all | 0.239 | 0.872 | 0.877 | 0.448 | 0.921 |
| Phikon | MSS_only | 0.270 | 0.859 | 0.848 | 0.463 | 0.909 |
| Phikon-v2 | all | 0.222 | 0.895 | 0.902 | 0.414 | 0.939 |
| Phikon-v2 | MSS_only | 0.269 | 0.887 | 0.888 | 0.465 | 0.927 |
| UNI | all | 0.278 | 0.933 | 0.931 | 0.491 | 0.965 |
| UNI | MSS_only | 0.289 | 0.921 | 0.917 | 0.475 | 0.955 |
| UNI2-h | all | 0.302 | 0.932 | 0.934 | 0.509 | 0.962 |
| UNI2-h | MSS_only | 0.296 | 0.920 | 0.917 | 0.471 | 0.953 |
| Virchow2 | all | 0.318 | 0.892 | 0.893 | 0.520 | 0.940 |
| Virchow2 | MSS_only | 0.296 | 0.891 | 0.875 | 0.467 | 0.939 |
| Prov-GigaPath | all | 0.290 | 0.928 | 0.930 | 0.510 | 0.955 |
| Prov-GigaPath | MSS_only | 0.291 | 0.924 | 0.921 | 0.495 | 0.952 |
| Tissue fractions | all | -0.000 | 0.498 | 0.606 | 0.162 | 0.693 |
| Tissue fractions | MSS_only | 0.015 | 0.501 | 0.604 | 0.181 | 0.704 |
| Nucleus densities | all | 0.013 | 0.487 | 0.635 | 0.113 | 0.677 |
| Nucleus densities | MSS_only | 0.037 | 0.474 | 0.645 | 0.212 | 0.662 |
| Tissue + nuclei | all | 0.001 | 0.615 | 0.676 | 0.180 | 0.768 |
| Tissue + nuclei | MSS_only | 0.053 | 0.616 | 0.686 | 0.244 | 0.770 |

# Interpretable features (patient mean over all patches)

| Feature | Mean MSIMUT (TRAIN) | Mean MSS (TRAIN) | AUC TRAIN | AUC TEST | E2 LR coef (std.) |
|---|---|---|---|---|---|
| tissue:ADI | 0.000 | 0.000 | 0.566 | 0.410 | +0.03 |
| tissue:BACK | 0.030 | 0.030 | 0.488 | 0.360 | -0.37 |
| tissue:DEB | 0.018 | 0.026 | 0.451 | 0.554 | -0.41 |
| tissue:LYM | 0.038 | 0.026 | 0.576 | 0.513 | +0.77 |
| tissue:MUC | 0.036 | 0.015 | 0.686 | 0.617 | +0.71 |
| tissue:MUS | 0.003 | 0.004 | 0.489 | 0.508 | -0.32 |
| tissue:NORM | 0.062 | 0.085 | 0.414 | 0.435 | -0.29 |
| tissue:STR | 0.162 | 0.133 | 0.595 | 0.615 | +0.51 |
| tissue:TUM | 0.650 | 0.681 | 0.443 | 0.466 | -0.24 |
| cells:neutrophil | 6.351 | 3.190 | 0.654 | 0.694 | +0.32 |
| cells:epithelial | 238.832 | 220.136 | 0.502 | 0.509 | +0.15 |
| cells:lymphocyte | 94.153 | 72.692 | 0.612 | 0.659 | -0.09 |
| cells:plasma | 7.599 | 8.605 | 0.464 | 0.600 | -0.38 |
| cells:eosinophil | 1.898 | 2.131 | 0.448 | 0.509 | +0.13 |
| cells:connective | 74.208 | 80.865 | 0.431 | 0.469 | -0.65 |
| cells:iel | 16.709 | 8.420 | 0.743 | 0.810 | +1.21 |
| cells:nuclei | 423.041 | 387.619 | 0.595 | 0.561 | -0.14 |

# Resolution: Zenodo tiles vs re-cut from the slides (ROC-AUC)

Zenodo tiles: as published (256 µm at 224 px = 1.14 µm/px, Macenko). Re-cut: the same 256 µm regions read from the TCGA slides (raw colour) at 1.14 µm/px (one 224 px view) or 0.5 µm/px (mean of four 224 px crops).

| Encoder | Aggregation | Zenodo tiles E1 | Zenodo tiles E2 | Re-cut 1.14 µm/px E1 | Re-cut 1.14 µm/px E2 | Re-cut 0.5 µm/px E1 | Re-cut 0.5 µm/px E2 |
|---|---|---|---|---|---|---|---|
| Phikon | Mean | 0.795 | 0.789 | 0.834 | 0.816 | 0.870 | 0.857 |
| Phikon | Max | 0.911 | 0.777 | 0.911 | 0.731 | 0.939 | 0.806 |
| Phikon | Mean+max | 0.911 | 0.797 | 0.916 | 0.764 | 0.936 | 0.830 |
| Phikon | ABMIL | 0.813 | 0.855 | 0.854 | 0.876 | 0.879 | 0.914 |
| Phikon-v2 | Mean | 0.812 | 0.821 | 0.816 | 0.786 | 0.828 | 0.822 |
| Phikon-v2 | Max | 0.929 | 0.745 | 0.914 | 0.760 | 0.904 | 0.766 |
| Phikon-v2 | Mean+max | 0.923 | 0.773 | 0.906 | 0.773 | 0.898 | 0.768 |
| Phikon-v2 | ABMIL | 0.810 | 0.879 | 0.816 | 0.793 | 0.817 | 0.832 |
| UNI | Mean | 0.880 | 0.851 | 0.855 | 0.854 | 0.875 | 0.871 |
| UNI | Max | 0.942 | 0.687 | 0.940 | 0.706 | 0.938 | 0.780 |
| UNI | Mean+max | 0.945 | 0.747 | 0.933 | 0.758 | 0.931 | 0.795 |
| UNI | ABMIL | 0.850 | 0.886 | 0.846 | 0.892 | 0.865 | 0.906 |
| UNI2-h | Mean | 0.894 | 0.885 | 0.889 | 0.885 | 0.892 | 0.884 |
| UNI2-h | Max | 0.952 | 0.773 | 0.955 | 0.774 | 0.946 | 0.800 |
| UNI2-h | Mean+max | 0.945 | 0.826 | 0.946 | 0.836 | 0.942 | 0.850 |
| UNI2-h | ABMIL | 0.886 | 0.929 | 0.883 | 0.931 | 0.913 | 0.927 |
| Virchow2 | Mean | 0.887 | 0.877 | 0.891 | 0.880 | 0.882 | 0.879 |
| Virchow2 | Max | 0.933 | 0.756 | 0.941 | 0.795 | 0.936 | 0.810 |
| Virchow2 | Mean+max | 0.936 | 0.833 | 0.940 | 0.857 | 0.936 | 0.853 |
| Virchow2 | ABMIL | 0.875 | 0.906 | 0.879 | 0.925 | 0.880 | 0.917 |
| Prov-GigaPath | Mean | 0.901 | 0.891 | 0.891 | 0.863 | 0.899 | 0.889 |
| Prov-GigaPath | Max | 0.951 | 0.792 | 0.951 | 0.782 | 0.946 | 0.800 |
| Prov-GigaPath | Mean+max | 0.954 | 0.838 | 0.947 | 0.819 | 0.943 | 0.849 |
| Prov-GigaPath | ABMIL | 0.884 | 0.908 | 0.875 | 0.898 | 0.911 | 0.932 |

Paired E2 differences (same TEST bootstrap resamples):

| Encoder | Aggregation | Comparison | Δ AUC [95% CI] |
|---|---|---|---|
| Phikon | Mean | recut_0.5 - zenodo | +0.068 [+0.007, +0.138] |
| Phikon | Mean | recut_1.14 - zenodo | +0.027 [-0.031, +0.093] |
| Phikon | Mean | recut_0.5 - recut_1.14 | +0.041 [+0.002, +0.084] |
| Phikon | Max | recut_0.5 - zenodo | +0.030 [-0.025, +0.083] |
| Phikon | Max | recut_1.14 - zenodo | -0.046 [-0.107, +0.009] |
| Phikon | Max | recut_0.5 - recut_1.14 | +0.075 [+0.023, +0.131] |
| Phikon | Mean+max | recut_0.5 - zenodo | +0.033 [-0.021, +0.084] |
| Phikon | Mean+max | recut_1.14 - zenodo | -0.033 [-0.086, +0.012] |
| Phikon | Mean+max | recut_0.5 - recut_1.14 | +0.067 [+0.018, +0.121] |
| Phikon | ABMIL | recut_0.5 - zenodo | +0.059 [+0.009, +0.119] |
| Phikon | ABMIL | recut_1.14 - zenodo | +0.021 [-0.015, +0.064] |
| Phikon | ABMIL | recut_0.5 - recut_1.14 | +0.038 [+0.005, +0.078] |
| Phikon-v2 | Mean | recut_0.5 - zenodo | +0.001 [-0.084, +0.091] |
| Phikon-v2 | Mean | recut_1.14 - zenodo | -0.035 [-0.112, +0.029] |
| Phikon-v2 | Mean | recut_0.5 - recut_1.14 | +0.036 [-0.038, +0.125] |
| Phikon-v2 | Max | recut_0.5 - zenodo | +0.021 [-0.047, +0.088] |
| Phikon-v2 | Max | recut_1.14 - zenodo | +0.016 [-0.032, +0.065] |
| Phikon-v2 | Max | recut_0.5 - recut_1.14 | +0.006 [-0.061, +0.072] |
| Phikon-v2 | Mean+max | recut_0.5 - zenodo | -0.006 [-0.072, +0.052] |
| Phikon-v2 | Mean+max | recut_1.14 - zenodo | -0.001 [-0.045, +0.041] |
| Phikon-v2 | Mean+max | recut_0.5 - recut_1.14 | -0.005 [-0.056, +0.043] |
| Phikon-v2 | ABMIL | recut_0.5 - zenodo | -0.047 [-0.122, +0.018] |
| Phikon-v2 | ABMIL | recut_1.14 - zenodo | -0.086 [-0.165, -0.023] |
| Phikon-v2 | ABMIL | recut_0.5 - recut_1.14 | +0.039 [-0.041, +0.129] |
| UNI | Mean | recut_0.5 - zenodo | +0.020 [-0.065, +0.095] |
| UNI | Mean | recut_1.14 - zenodo | +0.003 [-0.042, +0.046] |
| UNI | Mean | recut_0.5 - recut_1.14 | +0.017 [-0.035, +0.062] |
| UNI | Max | recut_0.5 - zenodo | +0.094 [+0.034, +0.157] |
| UNI | Max | recut_1.14 - zenodo | +0.020 [-0.020, +0.064] |
| UNI | Max | recut_0.5 - recut_1.14 | +0.074 [+0.023, +0.129] |
| UNI | Mean+max | recut_0.5 - zenodo | +0.048 [-0.015, +0.113] |
| UNI | Mean+max | recut_1.14 - zenodo | +0.011 [-0.033, +0.055] |
| UNI | Mean+max | recut_0.5 - recut_1.14 | +0.037 [-0.013, +0.087] |
| UNI | ABMIL | recut_0.5 - zenodo | +0.020 [-0.011, +0.056] |
| UNI | ABMIL | recut_1.14 - zenodo | +0.005 [-0.038, +0.041] |
| UNI | ABMIL | recut_0.5 - recut_1.14 | +0.014 [-0.005, +0.040] |
| UNI2-h | Mean | recut_0.5 - zenodo | -0.001 [-0.042, +0.035] |
| UNI2-h | Mean | recut_1.14 - zenodo | -0.000 [-0.026, +0.020] |
| UNI2-h | Mean | recut_0.5 - recut_1.14 | -0.001 [-0.038, +0.031] |
| UNI2-h | Max | recut_0.5 - zenodo | +0.027 [-0.015, +0.070] |
| UNI2-h | Max | recut_1.14 - zenodo | +0.001 [-0.030, +0.033] |
| UNI2-h | Max | recut_0.5 - recut_1.14 | +0.027 [-0.020, +0.073] |
| UNI2-h | Mean+max | recut_0.5 - zenodo | +0.024 [-0.009, +0.062] |
| UNI2-h | Mean+max | recut_1.14 - zenodo | +0.010 [-0.015, +0.037] |
| UNI2-h | Mean+max | recut_0.5 - recut_1.14 | +0.015 [-0.016, +0.046] |
| UNI2-h | ABMIL | recut_0.5 - zenodo | -0.002 [-0.018, +0.013] |
| UNI2-h | ABMIL | recut_1.14 - zenodo | +0.001 [-0.015, +0.019] |
| UNI2-h | ABMIL | recut_0.5 - recut_1.14 | -0.003 [-0.026, +0.015] |
| Virchow2 | Mean | recut_0.5 - zenodo | +0.002 [-0.018, +0.023] |
| Virchow2 | Mean | recut_1.14 - zenodo | +0.003 [-0.017, +0.022] |
| Virchow2 | Mean | recut_0.5 - recut_1.14 | -0.001 [-0.021, +0.021] |
| Virchow2 | Max | recut_0.5 - zenodo | +0.054 [+0.015, +0.097] |
| Virchow2 | Max | recut_1.14 - zenodo | +0.038 [+0.008, +0.070] |
| Virchow2 | Max | recut_0.5 - recut_1.14 | +0.016 [-0.024, +0.061] |
| Virchow2 | Mean+max | recut_0.5 - zenodo | +0.021 [-0.006, +0.051] |
| Virchow2 | Mean+max | recut_1.14 - zenodo | +0.024 [-0.003, +0.054] |
| Virchow2 | Mean+max | recut_0.5 - recut_1.14 | -0.003 [-0.027, +0.020] |
| Virchow2 | ABMIL | recut_0.5 - zenodo | +0.010 [-0.014, +0.042] |
| Virchow2 | ABMIL | recut_1.14 - zenodo | +0.019 [+0.001, +0.041] |
| Virchow2 | ABMIL | recut_0.5 - recut_1.14 | -0.008 [-0.035, +0.019] |
| Prov-GigaPath | Mean | recut_0.5 - zenodo | -0.002 [-0.037, +0.035] |
| Prov-GigaPath | Mean | recut_1.14 - zenodo | -0.028 [-0.062, +0.001] |
| Prov-GigaPath | Mean | recut_0.5 - recut_1.14 | +0.026 [-0.006, +0.064] |
| Prov-GigaPath | Max | recut_0.5 - zenodo | +0.008 [-0.035, +0.054] |
| Prov-GigaPath | Max | recut_1.14 - zenodo | -0.010 [-0.054, +0.030] |
| Prov-GigaPath | Max | recut_0.5 - recut_1.14 | +0.018 [-0.028, +0.068] |
| Prov-GigaPath | Mean+max | recut_0.5 - zenodo | +0.011 [-0.028, +0.051] |
| Prov-GigaPath | Mean+max | recut_1.14 - zenodo | -0.019 [-0.060, +0.017] |
| Prov-GigaPath | Mean+max | recut_0.5 - recut_1.14 | +0.030 [-0.016, +0.082] |
| Prov-GigaPath | ABMIL | recut_0.5 - zenodo | +0.025 [-0.003, +0.062] |
| Prov-GigaPath | ABMIL | recut_1.14 - zenodo | -0.010 [-0.036, +0.013] |
| Prov-GigaPath | ABMIL | recut_0.5 - recut_1.14 | +0.035 [+0.003, +0.076] |
