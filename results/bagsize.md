# External collapse of max pooling: shortcut or bag-size artefact?

Training on TCGA TRAIN with each patient's Kather tile count (confounded) or 200 tiles (deconfounded); single whole-tumor tiles at 0.5 µm/px throughout. AUC (mean over 10 bag draws; full = one draw).

## Training bags: kather_counts

| Encoder | Pooling | Cohort | full | TCGA-like n | n=200 | n=1000 | P>0.99 (full) |
|---|---|---|---|---|---|---|---|
| Prov-GigaPath-flash | max | tcga_test | 0.759 | 0.674 | 0.797 | 0.763 | 0.66 |
| Prov-GigaPath-flash | max | cptac | 0.601 | 0.562 | 0.656 | 0.622 | 0.69 |
| Prov-GigaPath-flash | max | surgen | 0.705 | 0.607 | 0.741 | 0.715 | 0.88 |
| Prov-GigaPath-flash | mean | tcga_test | 0.879 | 0.877 | 0.880 | 0.880 | 0.12 |
| Prov-GigaPath-flash | mean | cptac | 0.824 | 0.827 | 0.839 | 0.823 | 0.00 |
| Prov-GigaPath-flash | mean | surgen | 0.808 | 0.803 | 0.807 | 0.810 | 0.02 |
| UNI2-h | max | tcga_test | 0.846 | 0.730 | 0.889 | 0.881 | 0.54 |
| UNI2-h | max | cptac | 0.780 | 0.687 | 0.853 | 0.843 | 0.87 |
| UNI2-h | max | surgen | 0.768 | 0.648 | 0.817 | 0.785 | 0.70 |
| UNI2-h | mean | tcga_test | 0.902 | 0.902 | 0.903 | 0.903 | 0.11 |
| UNI2-h | mean | cptac | 0.884 | 0.889 | 0.887 | 0.888 | 0.00 |
| UNI2-h | mean | surgen | 0.919 | 0.916 | 0.918 | 0.918 | 0.01 |
| Virchow2 | max | tcga_test | 0.830 | 0.678 | 0.813 | 0.816 | 0.60 |
| Virchow2 | max | cptac | 0.610 | 0.551 | 0.592 | 0.601 | 0.77 |
| Virchow2 | max | surgen | 0.747 | 0.588 | 0.731 | 0.740 | 0.78 |
| Virchow2 | mean | tcga_test | 0.883 | 0.882 | 0.885 | 0.883 | 0.12 |
| Virchow2 | mean | cptac | 0.967 | 0.956 | 0.966 | 0.966 | 0.00 |
| Virchow2 | mean | surgen | 0.800 | 0.797 | 0.800 | 0.799 | 0.02 |

## Training bags: fixed200

| Encoder | Pooling | Cohort | full | TCGA-like n | n=200 | n=1000 | P>0.99 (full) |
|---|---|---|---|---|---|---|---|
| Prov-GigaPath-flash | max | tcga_test | 0.854 | 0.866 | 0.867 | 0.859 | 0.08 |
| Prov-GigaPath-flash | max | cptac | 0.884 | 0.885 | 0.898 | 0.882 | 0.00 |
| Prov-GigaPath-flash | max | surgen | 0.846 | 0.845 | 0.851 | 0.847 | 0.00 |
| Prov-GigaPath-flash | mean | tcga_test | 0.890 | 0.884 | 0.887 | 0.890 | 0.12 |
| Prov-GigaPath-flash | mean | cptac | 0.844 | 0.834 | 0.850 | 0.842 | 0.00 |
| Prov-GigaPath-flash | mean | surgen | 0.810 | 0.805 | 0.808 | 0.811 | 0.02 |
| UNI2-h | max | tcga_test | 0.899 | 0.905 | 0.907 | 0.903 | 0.01 |
| UNI2-h | max | cptac | 0.915 | 0.915 | 0.912 | 0.925 | 0.00 |
| UNI2-h | max | surgen | 0.887 | 0.892 | 0.895 | 0.892 | 0.00 |
| UNI2-h | mean | tcga_test | 0.898 | 0.898 | 0.899 | 0.899 | 0.11 |
| UNI2-h | mean | cptac | 0.876 | 0.879 | 0.877 | 0.875 | 0.00 |
| UNI2-h | mean | surgen | 0.912 | 0.909 | 0.911 | 0.911 | 0.00 |
| Virchow2 | max | tcga_test | 0.870 | 0.884 | 0.883 | 0.877 | 0.06 |
| Virchow2 | max | cptac | 0.885 | 0.867 | 0.890 | 0.878 | 0.02 |
| Virchow2 | max | surgen | 0.865 | 0.858 | 0.866 | 0.860 | 0.01 |
| Virchow2 | mean | tcga_test | 0.887 | 0.882 | 0.887 | 0.888 | 0.14 |
| Virchow2 | mean | cptac | 0.935 | 0.929 | 0.934 | 0.934 | 0.00 |
| Virchow2 | mean | surgen | 0.810 | 0.807 | 0.809 | 0.810 | 0.02 |

