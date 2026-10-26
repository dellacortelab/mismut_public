# Pathologist review

180 regions, 353 responses (rater_A: 177, rater_B: 176); free-text notes: 1

## a) High-attention vs control regions

| question   | value                    | label                                          |   high |   low |   n_high |   n_low |     p |   q_bh |
|:-----------|:-------------------------|:-----------------------------------------------|-------:|------:|---------:|--------:|------:|-------:|
| q1         | tumor_epithelium         | Tumor epithelium                               |  0.714 | 0.21  |      234 |     119 | 0     |  0     |
| q1         | tumor_stroma             | Tumor stroma                                   |  0.026 | 0.16  |      234 |     119 | 0     |  0     |
| q1         | interface_front          | Tumor-stroma interface / invasive front        |  0.214 | 0.176 |      234 |     119 | 0.483 |  0.539 |
| q1         | normal_mucosa            | Normal mucosa                                  |  0.009 | 0.235 |      234 |     119 | 0     |  0     |
| q1         | lymphoid_aggregate       | Lymphoid aggregate                             |  0     | 0.076 |      234 |     119 | 0     |  0     |
| q1         | mucin                    | Mucin                                          |  0.021 | 0     |      234 |     119 | 0.172 |  0.217 |
| q1         | necrosis_debris          | Necrosis / debris                              |  0.017 | 0.109 |      234 |     119 | 0     |  0.001 |
| q1         | other                    | Other                                          |  0     | 0.034 |      234 |     119 | 0.012 |  0.018 |
| q2         | iel                      | Intraepithelial lymphocytes                    |  0.427 | 0.118 |      234 |     119 | 0     |  0     |
| q2         | stromal_lymphocytes      | Stromal / peritumoral lymphocytes              |  0.205 | 0.235 |      234 |     119 | 0.584 |  0.627 |
| q2         | neutrophils              | Neutrophils (epithelium or lumen)              |  0.534 | 0.37  |      234 |     119 | 0.005 |  0.007 |
| q2         | plasma_cells             | Plasma cells                                   |  0.145 | 0.252 |      234 |     119 | 0.019 |  0.026 |
| q2         | eosinophils              | Eosinophils                                    |  0.154 | 0.235 |      234 |     119 | 0.079 |  0.104 |
| q2         | crohn_like               | Crohn-like lymphoid reaction                   |  0     | 0.059 |      234 |     119 | 0     |  0.001 |
| q2         | none                     | None                                           |  0.103 | 0.227 |      234 |     119 | 0.002 |  0.004 |
| q3         | glandular                | Gland-forming (well/moderately differentiated) |  0.491 | 0.294 |      234 |     119 | 0     |  0.001 |
| q3         | poorly_diff_solid        | Poorly differentiated / solid                  |  0.128 | 0.084 |      234 |     119 | 0.286 |  0.346 |
| q3         | medullary                | Medullary                                      |  0.167 | 0.042 |      234 |     119 | 0.001 |  0.001 |
| q3         | mucinous                 | Mucinous                                       |  0.235 | 0.084 |      234 |     119 | 0     |  0.001 |
| q3         | signet_ring              | Signet ring cells                              |  0.073 | 0.059 |      234 |     119 | 0.823 |  0.823 |
| q3         | cribriform_comedo        | Cribriform / comedo                            |  0.132 | 0.008 |      234 |     119 | 0     |  0     |
| q3         | dirty_necrosis           | Dirty (garland) necrosis                       |  0.137 | 0.092 |      234 |     119 | 0.302 |  0.35  |
| q3         | not_tumor                | Not tumor                                      |  0.03  | 0.538 |      234 |     119 | 0     |  0     |
| q8         | large_cells              | Large tumor cells                              |  0.713 | 0.269 |      164 |      78 | 0     |  0     |
| q8         | eosinophilic_amphophilic | Moderate eosinophilic or amphophilic cytoplasm |  0.701 | 0.333 |      164 |      78 | 0     |  0     |
| q8         | vesicular_chromatin      | Vesicular chromatin                            |  0.689 | 0.282 |      164 |      78 | 0     |  0     |
| q8         | prominent_nucleoli       | Prominent nucleoli                             |  0.537 | 0.205 |      164 |      78 | 0     |  0     |
| q8         | none                     | None of these                                  |  0.067 | 0.077 |      164 |      78 | 0.791 |  0.82  |
| q8         | not_assessable           | Not assessable (no tumor in region)            |  0.03  | 0.5   |      164 |      78 | 0     |  0     |

## b) Looks MSI-associated?

| cell                            |   yes |   unsure |    no |
|:--------------------------------|------:|---------:|------:|
| Control regions                 | 0.185 |    0.05  | 0.765 |
| High attention: MSI             | 0.457 |    0.043 | 0.5   |
| High attention: MSS             | 0.206 |    0.029 | 0.765 |
| High attention: MSS, called MSI | 0.361 |    0.111 | 0.528 |

| cell                            |   0 |
|:--------------------------------|----:|
| Control regions                 | 119 |
| High attention: MSI             | 164 |
| High attention: MSS             |  34 |
| High attention: MSS, called MSI |  36 |

slide AUC of the share of 'yes' (MSI vs MSS patients): 0.577 (Mann-Whitney p = 0.62)

## d) Agreement (Cohen's kappa)

- Main tissue: 0.516
- MSI-associated call: 0.342
- Inflammation (mean over categories): 0.303
- Neutrophils: 0.375
- IEL: 0.334
- n regions rated by both: 173

## e) Tumor epithelium only (exploratory)

| question   | value               |   high |   low |   n_high |   n_low |     p |
|:-----------|:--------------------|-------:|------:|---------:|--------:|------:|
| q2         | iel                 |  0.491 | 0.16  |      167 |      25 | 0.002 |
| q2         | neutrophils         |  0.503 | 0.48  |      167 |      25 | 1     |
| q2         | stromal_lymphocytes |  0.168 | 0.12  |      167 |      25 | 0.772 |
| q2         | plasma_cells        |  0.12  | 0.08  |      167 |      25 | 0.744 |
| q8         | large_cells         |  0.741 | 0.533 |      108 |      15 | 0.126 |
| q8         | vesicular_chromatin |  0.741 | 0.667 |      108 |      15 | 0.544 |
| q8         | prominent_nucleoli  |  0.611 | 0.4   |      108 |      15 | 0.163 |

## f) High-attention regions by patient group (exploratory)

| question   | value               |   MSI |   n MSI |   MSS, called MSI |   n MSS, called MSI |   MSS |   n MSS |   control |   p response (MSI vs all MSS) |   slide AUC |   slide p (Mann-Whitney) |
|:-----------|:--------------------|------:|--------:|------------------:|--------------------:|------:|--------:|----------:|------------------------------:|------------:|-------------------------:|
| q2         | neutrophils         | 0.671 |     164 |             0.306 |                  36 | 0.118 |      34 |     0.37  |                         0     |       0.911 |                    0.005 |
| q2         | iel                 | 0.445 |     164 |             0.472 |                  36 | 0.294 |      34 |     0.118 |                         0.471 |       0.601 |                    0.507 |
| q2         | stromal_lymphocytes | 0.171 |     164 |             0.361 |                  36 | 0.206 |      34 |     0.235 |                         0.053 |       0.351 |                    0.314 |
| q2         | none                | 0.037 |     164 |             0.194 |                  36 | 0.324 |      34 |     0.227 |                         0     |       0.131 |                    0.006 |
| q3         | mucinous            | 0.262 |     164 |             0.333 |                  36 | 0     |      34 |     0.084 |                         0.178 |       0.643 |                    0.319 |
| q3         | medullary           | 0.189 |     164 |             0.028 |                  36 | 0.206 |      34 |     0.042 |                         0.184 |       0.476 |                    0.889 |
| q3         | poorly_diff_solid   | 0.165 |     164 |             0     |                  36 | 0.088 |      34 |     0.084 |                         0.01  |       0.571 |                    0.594 |
| q3         | cribriform_comedo   | 0.11  |     164 |             0.083 |                  36 | 0.294 |      34 |     0.008 |                         0.141 |       0.333 |                    0.208 |
| q8         | prominent_nucleoli  | 0.617 |     120 |             0.321 |                  28 | 0.312 |      16 |     0.205 |                         0.001 |       0.75  |                    0.09  |
| q8         | vesicular_chromatin | 0.733 |     120 |             0.536 |                  28 | 0.625 |      16 |     0.282 |                         0.057 |       0.548 |                    0.771 |
| q8         | large_cells         | 0.708 |     120 |             0.643 |                  28 | 0.875 |      16 |     0.269 |                         0.848 |       0.518 |                    0.934 |

Per rater:

| rater   |   neutrophils slide AUC |   neutrophils MSI |   neutrophils MSS, called MSI |   neutrophils MSS |   iel slide AUC |   iel MSI |   iel MSS, called MSI |   iel MSS |   MSI call slide AUC |   responses |
|:--------|------------------------:|------------------:|------------------------------:|------------------:|----------------:|----------:|----------------------:|----------:|---------------------:|------------:|
| rater_A |                   0.869 |             0.807 |                         0.556 |             0.188 |           0.744 |     0.337 |                 0.333 |     0.062 |                0.69  |         177 |
| rater_B |                   0.893 |             0.531 |                         0.056 |             0.056 |           0.47  |     0.556 |                 0.611 |     0.5   |                0.512 |         176 |

Per slide:

| patient      | group              |   neutrophils |   iel |   stromal_lymphocytes |   none |   mucinous |   medullary |   poorly_diff_solid |   cribriform_comedo |   prominent_nucleoli |   vesicular_chromatin |   large_cells |   yes |
|:-------------|:-------------------|--------------:|------:|----------------------:|-------:|-----------:|------------:|--------------------:|--------------------:|---------------------:|----------------------:|--------------:|------:|
| TCGA-D5-6530 | MSI_correct        |         0.667 | 0.417 |                 0.25  |  0     |      0     |       0     |               0     |               0     |                0.6   |                 0.6   |         0.8   | 0     |
| TCGA-G4-6588 | MSI_correct        |         0.583 | 0.333 |                 0.083 |  0.167 |      0.25  |       0     |               0     |               0     |                0.778 |                 0.778 |         0.889 | 0.167 |
| TCGA-AA-3679 | MSS_correct        |         0.091 | 0.182 |                 0     |  0.636 |      0     |       0.182 |               0.091 |               0.455 |                0.2   |                 1     |         0.8   | 0.182 |
| TCGA-CI-6624 | MSS_correct        |         0     | 0.417 |                 0.333 |  0.167 |      0     |       0     |               0     |               0     |                0     |                 0     |         0.833 | 0     |
| TCGA-CM-4743 | MSI_correct        |         0.667 | 0.667 |                 0.167 |  0     |      0.083 |       0.5   |               0.5   |               0.083 |                1     |                 1     |         0.875 | 0.917 |
| TCGA-AA-3833 | MSI_correct        |         0.909 | 0.182 |                 0     |  0.091 |      0     |       0     |               0.273 |               0     |                1     |                 0.889 |         1     | 0.091 |
| TCGA-AY-A71X | MSS_false_positive |         0.25  | 0.25  |                 0.667 |  0.25  |      0.083 |       0.083 |               0     |               0.167 |                0     |                 0.111 |         0.667 | 0.167 |
| TCGA-D5-6540 | MSI_correct        |         0.909 | 0.364 |                 0.455 |  0     |      0.455 |       0     |               0     |               0.636 |                0.556 |                 0.778 |         0.667 | 0.455 |
| TCGA-AD-6964 | MSI_correct        |         0.583 | 0.583 |                 0.25  |  0     |      0     |       0.417 |               0.583 |               0     |                0.889 |                 1     |         0.778 | 1     |
| TCGA-CL-5917 | MSS_correct        |         0.273 | 0.273 |                 0.273 |  0.182 |      0     |       0.455 |               0.182 |               0.455 |                0.8   |                 1     |         1     | 0.455 |
| TCGA-CK-6746 | MSI_correct        |         0     | 0.909 |                 0.182 |  0     |      0     |       0.818 |               0.455 |               0     |                0.667 |                 0.889 |         0.778 | 1     |
| TCGA-AZ-6598 | MSI_correct        |         0.333 | 0.583 |                 0     |  0.083 |      0.083 |       0.833 |               0.5   |               0     |                0.5   |                 0.667 |         0.167 | 1     |
| TCGA-G4-6586 | MSI_correct        |         1     | 0.667 |                 0.25  |  0     |      0.75  |       0     |               0     |               0.167 |                0.375 |                 0.75  |         0.875 | 0.417 |
| TCGA-G4-6302 | MSI_correct        |         0.833 | 0.083 |                 0.083 |  0     |      0.417 |       0     |               0     |               0     |                0.444 |                 0.333 |         0.667 | 0.25  |
| TCGA-AG-3881 | MSS_false_positive |         0.5   | 0.167 |                 0     |  0.333 |      0.917 |       0     |               0     |               0     |                0.25  |                 0.5   |         0.5   | 0.417 |
| TCGA-AY-6197 | MSI_correct        |         0.818 | 0.182 |                 0     |  0.182 |      0.364 |       0     |               0     |               0     |                0.5   |                 0.7   |         0.9   | 0.091 |
| TCGA-WS-AB45 | MSI_correct        |         1     | 0     |                 0     |  0     |      0.833 |       0     |               0     |               0     |                0     |                 0.2   |         0     | 0.667 |
| TCGA-AZ-4615 | MSI_correct        |         0.333 | 0.667 |                 0.25  |  0     |      0     |       0.083 |               0     |               0.667 |                0.857 |                 1     |         1     | 0     |
| TCGA-AG-3892 | MSS_false_positive |         0.167 | 1     |                 0.417 |  0     |      0     |       0     |               0     |               0.083 |                0.636 |                 0.909 |         0.727 | 0.5   |
| TCGA-CK-4951 | MSI_correct        |         0.75  | 0.583 |                 0.417 |  0     |      0.417 |       0     |               0     |               0     |                0.571 |                 0.857 |         0.429 | 0.333 |

## Checks (numbers quoted in the text)

- high-attention regions entirely inside tumor map: 0.800
- high-attention region median area mm2: 0.383
- MSI high neutrophil responses: 110
- ... of which dirty necrosis not ticked: 86
- control neutrophils MSI: 0.381
- control neutrophils MSS: 0.353
- control neutrophils MSS, called MSI: 0.333
- high dirty necrosis MSI: 0.159
- high dirty necrosis MSS: 0.118
- high dirty necrosis MSS, called MSI: 0.056
- neutrophil patient AUC, excluding responses with dirty necrosis: 0.810
- neutrophil patient AUC, tumor-epithelium responses: 0.893
- neutrophil patient AUC, without POLE patient: 0.907
- KongNet iel density AUC on the 20 review patients: 0.905
- KongNet neutrophil density AUC on the 20 review patients: 0.833

## h) Molecular annotation of review patients

| patient      | group              |     tmb |   mantis | pole                            |
|:-------------|:-------------------|--------:|---------:|:--------------------------------|
| TCGA-AA-3833 | MSI_correct        |  19.133 |    0.955 | R1160H                          |
| TCGA-AD-6964 | MSI_correct        |  47.3   |    0.736 | Y473C                           |
| TCGA-AY-6197 | MSI_correct        |  35.267 |    0.965 | nan                             |
| TCGA-AZ-4615 | MSI_correct        |  27.833 |    1.033 | K325Sfs*19                      |
| TCGA-AZ-6598 | MSI_correct        |  77.833 |    1.378 | T1052M;X1517_splice             |
| TCGA-CK-4951 | MSI_correct        | 140.467 |    0.785 | E2137K                          |
| TCGA-CK-6746 | MSI_correct        |  40.133 |    1.168 | nan                             |
| TCGA-CM-4743 | MSI_correct        |  36.867 |    0.818 | nan                             |
| TCGA-D5-6530 | MSI_correct        |  35.167 |    1.086 | nan                             |
| TCGA-D5-6540 | MSI_correct        |  40.033 |    0.998 | A1224T                          |
| TCGA-G4-6302 | MSI_correct        |  59.733 |    0.459 | nan                             |
| TCGA-G4-6586 | MSI_correct        |  37.033 |    1.061 | nan                             |
| TCGA-G4-6588 | MSI_correct        |  65.567 |    1.007 | A2030Pfs*18                     |
| TCGA-WS-AB45 | MSI_correct        |  84.333 |    0.678 | V2152M;R924H;K1383R;A2030Pfs*18 |
| TCGA-AG-3881 | MSS_false_positive |   3.7   |    0.284 | nan                             |
| TCGA-AG-3892 | MSS_false_positive |  83.367 |    0.588 | S459F                           |
| TCGA-AY-A71X | MSS_false_positive |   3     |    0.306 | nan                             |
| TCGA-AA-3679 | MSS_correct        |   2.033 |    0.258 | nan                             |
| TCGA-CI-6624 | MSS_correct        |   8.5   |    0.34  | nan                             |
| TCGA-CL-5917 | MSS_correct        |   2.333 |    0.359 | nan                             |

## f2) Patient-level AUC with bootstrap CI and BH q (high-attention regions)

| question   | value                    |   auc |    lo |    hi |     p |   q_bh |
|:-----------|:-------------------------|------:|------:|------:|------:|-------:|
| q2         | neutrophils              | 0.911 | 0.75  | 1     | 0.005 |  0.047 |
| q2         | iel                      | 0.601 | 0.321 | 0.881 | 0.507 |  0.902 |
| q2         | stromal_lymphocytes      | 0.351 | 0.048 | 0.679 | 0.314 |  0.638 |
| q2         | plasma_cells             | 0.482 | 0.22  | 0.744 | 0.933 |  0.934 |
| q2         | eosinophils              | 0.435 | 0.143 | 0.738 | 0.675 |  0.934 |
| q2         | none                     | 0.131 | 0     | 0.345 | 0.006 |  0.047 |
| q3         | mucinous                 | 0.643 | 0.357 | 0.893 | 0.319 |  0.638 |
| q3         | medullary                | 0.476 | 0.226 | 0.696 | 0.889 |  0.934 |
| q3         | poorly_diff_solid        | 0.571 | 0.357 | 0.762 | 0.594 |  0.934 |
| q3         | cribriform_comedo        | 0.333 | 0.095 | 0.565 | 0.208 |  0.638 |
| q3         | dirty_necrosis           | 0.685 | 0.399 | 0.929 | 0.204 |  0.638 |
| q3         | glandular                | 0.333 | 0.089 | 0.619 | 0.257 |  0.638 |
| q8         | prominent_nucleoli       | 0.75  | 0.464 | 1     | 0.09  |  0.481 |
| q8         | vesicular_chromatin      | 0.548 | 0.202 | 0.869 | 0.771 |  0.934 |
| q8         | large_cells              | 0.518 | 0.262 | 0.786 | 0.934 |  0.934 |
| q8         | eosinophilic_amphophilic | 0.44  | 0.179 | 0.726 | 0.708 |  0.934 |

## g) Attention-gated counts and controls (see Supplementary Table 5)

| cohort    | gate   | feature    |   auc |    lo |    hi |   fm_abmil |
|:----------|:-------|:-----------|------:|------:|------:|-----------:|
| tcga_test | all    | neutrophil | 0.703 | 0.567 | 0.828 |      0.919 |
| tcga_test | all    | epithelial | 0.498 | 0.367 | 0.632 |      0.919 |
| tcga_test | all    | lymphocyte | 0.712 | 0.596 | 0.818 |      0.919 |
| tcga_test | all    | plasma     | 0.632 | 0.494 | 0.758 |      0.919 |
| tcga_test | all    | eosinophil | 0.458 | 0.334 | 0.58  |      0.919 |
| tcga_test | all    | connective | 0.541 | 0.416 | 0.667 |      0.919 |
| tcga_test | all    | iel        | 0.793 | 0.663 | 0.901 |      0.919 |
| tcga_test | all    | nuclei     | 0.584 | 0.443 | 0.716 |      0.919 |
| tcga_test | att1   | neutrophil | 0.794 | 0.684 | 0.89  |      0.919 |
| tcga_test | att1   | epithelial | 0.565 | 0.42  | 0.707 |      0.919 |
| tcga_test | att1   | lymphocyte | 0.879 | 0.784 | 0.952 |      0.919 |
| tcga_test | att1   | plasma     | 0.784 | 0.661 | 0.9   |      0.919 |
| tcga_test | att1   | eosinophil | 0.732 | 0.608 | 0.846 |      0.919 |
| tcga_test | att1   | connective | 0.798 | 0.687 | 0.894 |      0.919 |
| tcga_test | att1   | iel        | 0.868 | 0.77  | 0.943 |      0.919 |
| tcga_test | att1   | nuclei     | 0.696 | 0.555 | 0.828 |      0.919 |
| tcga_test | perm1  | neutrophil | 0.764 | 0.662 | 0.859 |      0.919 |
| tcga_test | perm1  | epithelial | 0.639 | 0.511 | 0.763 |      0.919 |
| tcga_test | perm1  | lymphocyte | 0.664 | 0.536 | 0.781 |      0.919 |
| tcga_test | perm1  | plasma     | 0.516 | 0.387 | 0.644 |      0.919 |
| tcga_test | perm1  | eosinophil | 0.49  | 0.37  | 0.605 |      0.919 |
| tcga_test | perm1  | connective | 0.477 | 0.36  | 0.598 |      0.919 |
| tcga_test | perm1  | iel        | 0.767 | 0.65  | 0.871 |      0.919 |
| tcga_test | perm1  | nuclei     | 0.695 | 0.569 | 0.808 |      0.919 |
| tcga_test | epi1   | neutrophil | 0.689 | 0.57  | 0.799 |      0.919 |
| tcga_test | epi1   | epithelial | 0.484 | 0.357 | 0.611 |      0.919 |
| tcga_test | epi1   | lymphocyte | 0.756 | 0.644 | 0.858 |      0.919 |
| tcga_test | epi1   | plasma     | 0.554 | 0.434 | 0.675 |      0.919 |
| tcga_test | epi1   | eosinophil | 0.524 | 0.419 | 0.633 |      0.919 |
| tcga_test | epi1   | connective | 0.541 | 0.41  | 0.669 |      0.919 |
| tcga_test | epi1   | iel        | 0.777 | 0.672 | 0.869 |      0.919 |
| tcga_test | epi1   | nuclei     | 0.507 | 0.378 | 0.638 |      0.919 |
| tcga_test | att5   | neutrophil | 0.801 | 0.682 | 0.904 |      0.919 |
| tcga_test | att5   | epithelial | 0.499 | 0.356 | 0.642 |      0.919 |
| tcga_test | att5   | lymphocyte | 0.881 | 0.799 | 0.949 |      0.919 |
| tcga_test | att5   | plasma     | 0.813 | 0.707 | 0.908 |      0.919 |
| tcga_test | att5   | eosinophil | 0.783 | 0.669 | 0.883 |      0.919 |
| tcga_test | att5   | connective | 0.826 | 0.727 | 0.91  |      0.919 |
| tcga_test | att5   | iel        | 0.844 | 0.732 | 0.934 |      0.919 |
| tcga_test | att5   | nuclei     | 0.654 | 0.514 | 0.786 |      0.919 |
| tcga_test | perm5  | neutrophil | 0.723 | 0.613 | 0.827 |      0.919 |
| tcga_test | perm5  | epithelial | 0.561 | 0.431 | 0.683 |      0.919 |
| tcga_test | perm5  | lymphocyte | 0.682 | 0.554 | 0.8   |      0.919 |
| tcga_test | perm5  | plasma     | 0.533 | 0.397 | 0.664 |      0.919 |
| tcga_test | perm5  | eosinophil | 0.433 | 0.311 | 0.561 |      0.919 |
| tcga_test | perm5  | connective | 0.506 | 0.38  | 0.635 |      0.919 |
| tcga_test | perm5  | iel        | 0.761 | 0.63  | 0.875 |      0.919 |
| tcga_test | perm5  | nuclei     | 0.626 | 0.492 | 0.752 |      0.919 |
| tcga_test | epi5   | neutrophil | 0.664 | 0.522 | 0.796 |      0.919 |
| tcga_test | epi5   | epithelial | 0.506 | 0.379 | 0.633 |      0.919 |
| tcga_test | epi5   | lymphocyte | 0.731 | 0.613 | 0.834 |      0.919 |
| tcga_test | epi5   | plasma     | 0.544 | 0.42  | 0.665 |      0.919 |
| tcga_test | epi5   | eosinophil | 0.437 | 0.318 | 0.558 |      0.919 |
| tcga_test | epi5   | connective | 0.511 | 0.39  | 0.633 |      0.919 |
| tcga_test | epi5   | iel        | 0.791 | 0.679 | 0.885 |      0.919 |
| tcga_test | epi5   | nuclei     | 0.516 | 0.385 | 0.644 |      0.919 |
| tcga_test | att10  | neutrophil | 0.787 | 0.667 | 0.892 |      0.919 |
| tcga_test | att10  | epithelial | 0.488 | 0.345 | 0.634 |      0.919 |
| tcga_test | att10  | lymphocyte | 0.877 | 0.786 | 0.948 |      0.919 |
| tcga_test | att10  | plasma     | 0.825 | 0.738 | 0.907 |      0.919 |
| tcga_test | att10  | eosinophil | 0.725 | 0.599 | 0.842 |      0.919 |
| tcga_test | att10  | connective | 0.854 | 0.773 | 0.92  |      0.919 |
| tcga_test | att10  | iel        | 0.848 | 0.729 | 0.944 |      0.919 |
| tcga_test | att10  | nuclei     | 0.64  | 0.498 | 0.772 |      0.919 |
| tcga_test | perm10 | neutrophil | 0.687 | 0.567 | 0.802 |      0.919 |
| tcga_test | perm10 | epithelial | 0.557 | 0.427 | 0.684 |      0.919 |
| tcga_test | perm10 | lymphocyte | 0.691 | 0.562 | 0.809 |      0.919 |
| tcga_test | perm10 | plasma     | 0.517 | 0.388 | 0.65  |      0.919 |
| tcga_test | perm10 | eosinophil | 0.402 | 0.28  | 0.53  |      0.919 |
| tcga_test | perm10 | connective | 0.496 | 0.367 | 0.63  |      0.919 |
| tcga_test | perm10 | iel        | 0.749 | 0.608 | 0.87  |      0.919 |
| tcga_test | perm10 | nuclei     | 0.61  | 0.469 | 0.744 |      0.919 |
| tcga_test | epi10  | neutrophil | 0.712 | 0.583 | 0.827 |      0.919 |
| tcga_test | epi10  | epithelial | 0.516 | 0.387 | 0.644 |      0.919 |
| tcga_test | epi10  | lymphocyte | 0.746 | 0.633 | 0.847 |      0.919 |
| tcga_test | epi10  | plasma     | 0.565 | 0.438 | 0.687 |      0.919 |
| tcga_test | epi10  | eosinophil | 0.439 | 0.316 | 0.562 |      0.919 |
| tcga_test | epi10  | connective | 0.526 | 0.4   | 0.657 |      0.919 |
| tcga_test | epi10  | iel        | 0.822 | 0.708 | 0.913 |      0.919 |
| tcga_test | epi10  | nuclei     | 0.534 | 0.4   | 0.663 |      0.919 |
| tcga_test | att20  | neutrophil | 0.782 | 0.667 | 0.89  |      0.919 |
| tcga_test | att20  | epithelial | 0.475 | 0.333 | 0.62  |      0.919 |
| tcga_test | att20  | lymphocyte | 0.865 | 0.774 | 0.939 |      0.919 |
| tcga_test | att20  | plasma     | 0.808 | 0.719 | 0.889 |      0.919 |
| tcga_test | att20  | eosinophil | 0.7   | 0.581 | 0.816 |      0.919 |
| tcga_test | att20  | connective | 0.83  | 0.741 | 0.907 |      0.919 |
| tcga_test | att20  | iel        | 0.854 | 0.739 | 0.946 |      0.919 |
| tcga_test | att20  | nuclei     | 0.626 | 0.489 | 0.755 |      0.919 |
| tcga_test | perm20 | neutrophil | 0.691 | 0.564 | 0.81  |      0.919 |
| tcga_test | perm20 | epithelial | 0.556 | 0.419 | 0.685 |      0.919 |
| tcga_test | perm20 | lymphocyte | 0.7   | 0.571 | 0.82  |      0.919 |
| tcga_test | perm20 | plasma     | 0.522 | 0.394 | 0.654 |      0.919 |
| tcga_test | perm20 | eosinophil | 0.382 | 0.266 | 0.503 |      0.919 |
| tcga_test | perm20 | connective | 0.483 | 0.356 | 0.618 |      0.919 |
| tcga_test | perm20 | iel        | 0.766 | 0.628 | 0.885 |      0.919 |
| tcga_test | perm20 | nuclei     | 0.608 | 0.466 | 0.743 |      0.919 |
| tcga_test | epi20  | neutrophil | 0.729 | 0.606 | 0.838 |      0.919 |
| tcga_test | epi20  | epithelial | 0.527 | 0.396 | 0.657 |      0.919 |
| tcga_test | epi20  | lymphocyte | 0.771 | 0.66  | 0.867 |      0.919 |
| tcga_test | epi20  | plasma     | 0.593 | 0.465 | 0.715 |      0.919 |
| tcga_test | epi20  | eosinophil | 0.463 | 0.34  | 0.584 |      0.919 |
| tcga_test | epi20  | connective | 0.554 | 0.423 | 0.689 |      0.919 |
| tcga_test | epi20  | iel        | 0.818 | 0.702 | 0.915 |      0.919 |
| tcga_test | epi20  | nuclei     | 0.556 | 0.422 | 0.682 |      0.919 |
| cptac     | all    | neutrophil | 0.613 | 0.448 | 0.768 |      0.951 |
| cptac     | all    | epithelial | 0.622 | 0.435 | 0.799 |      0.951 |
| cptac     | all    | lymphocyte | 0.574 | 0.405 | 0.74  |      0.951 |
| cptac     | all    | plasma     | 0.533 | 0.368 | 0.692 |      0.951 |
| cptac     | all    | eosinophil | 0.426 | 0.267 | 0.591 |      0.951 |
| cptac     | all    | connective | 0.481 | 0.29  | 0.667 |      0.951 |
| cptac     | all    | iel        | 0.732 | 0.549 | 0.893 |      0.951 |
| cptac     | all    | nuclei     | 0.625 | 0.436 | 0.804 |      0.951 |
| cptac     | att1   | neutrophil | 0.725 | 0.567 | 0.866 |      0.951 |
| cptac     | att1   | epithelial | 0.625 | 0.421 | 0.823 |      0.951 |
| cptac     | att1   | lymphocyte | 0.786 | 0.626 | 0.93  |      0.951 |
| cptac     | att1   | plasma     | 0.732 | 0.575 | 0.876 |      0.951 |
| cptac     | att1   | eosinophil | 0.576 | 0.473 | 0.693 |      0.951 |
| cptac     | att1   | connective | 0.514 | 0.321 | 0.7   |      0.951 |
| cptac     | att1   | iel        | 0.803 | 0.639 | 0.946 |      0.951 |
| cptac     | att1   | nuclei     | 0.677 | 0.47  | 0.863 |      0.951 |
| cptac     | perm1  | neutrophil | 0.624 | 0.478 | 0.753 |      0.951 |
| cptac     | perm1  | epithelial | 0.594 | 0.414 | 0.761 |      0.951 |
| cptac     | perm1  | lymphocyte | 0.485 | 0.324 | 0.651 |      0.951 |
| cptac     | perm1  | plasma     | 0.263 | 0.125 | 0.423 |      0.951 |
| cptac     | perm1  | eosinophil | 0.51  | 0.399 | 0.636 |      0.951 |
| cptac     | perm1  | connective | 0.411 | 0.247 | 0.588 |      0.951 |
| cptac     | perm1  | iel        | 0.638 | 0.485 | 0.783 |      0.951 |
| cptac     | perm1  | nuclei     | 0.536 | 0.354 | 0.71  |      0.951 |
| cptac     | epi1   | neutrophil | 0.691 | 0.516 | 0.858 |      0.951 |
| cptac     | epi1   | epithelial | 0.561 | 0.378 | 0.74  |      0.951 |
| cptac     | epi1   | lymphocyte | 0.734 | 0.583 | 0.871 |      0.951 |
| cptac     | epi1   | plasma     | 0.548 | 0.367 | 0.726 |      0.951 |
| cptac     | epi1   | eosinophil | 0.515 | 0.406 | 0.638 |      0.951 |
| cptac     | epi1   | connective | 0.4   | 0.246 | 0.559 |      0.951 |
| cptac     | epi1   | iel        | 0.804 | 0.673 | 0.914 |      0.951 |
| cptac     | epi1   | nuclei     | 0.579 | 0.397 | 0.754 |      0.951 |
| cptac     | att5   | neutrophil | 0.746 | 0.599 | 0.873 |      0.951 |
| cptac     | att5   | epithelial | 0.592 | 0.388 | 0.795 |      0.951 |
| cptac     | att5   | lymphocyte | 0.81  | 0.65  | 0.949 |      0.951 |
| cptac     | att5   | plasma     | 0.798 | 0.669 | 0.906 |      0.951 |
| cptac     | att5   | eosinophil | 0.64  | 0.477 | 0.807 |      0.951 |
| cptac     | att5   | connective | 0.6   | 0.417 | 0.778 |      0.951 |
| cptac     | att5   | iel        | 0.768 | 0.585 | 0.933 |      0.951 |
| cptac     | att5   | nuclei     | 0.67  | 0.466 | 0.857 |      0.951 |
| cptac     | perm5  | neutrophil | 0.534 | 0.368 | 0.691 |      0.951 |
| cptac     | perm5  | epithelial | 0.661 | 0.485 | 0.821 |      0.951 |
| cptac     | perm5  | lymphocyte | 0.493 | 0.326 | 0.661 |      0.951 |
| cptac     | perm5  | plasma     | 0.387 | 0.218 | 0.568 |      0.951 |
| cptac     | perm5  | eosinophil | 0.384 | 0.216 | 0.564 |      0.951 |
| cptac     | perm5  | connective | 0.378 | 0.219 | 0.551 |      0.951 |
| cptac     | perm5  | iel        | 0.754 | 0.588 | 0.898 |      0.951 |
| cptac     | perm5  | nuclei     | 0.612 | 0.424 | 0.79  |      0.951 |
| cptac     | epi5   | neutrophil | 0.709 | 0.527 | 0.874 |      0.951 |
| cptac     | epi5   | epithelial | 0.583 | 0.39  | 0.769 |      0.951 |
| cptac     | epi5   | lymphocyte | 0.684 | 0.528 | 0.821 |      0.951 |
| cptac     | epi5   | plasma     | 0.555 | 0.396 | 0.708 |      0.951 |
| cptac     | epi5   | eosinophil | 0.45  | 0.307 | 0.594 |      0.951 |
| cptac     | epi5   | connective | 0.519 | 0.35  | 0.682 |      0.951 |
| cptac     | epi5   | iel        | 0.747 | 0.58  | 0.891 |      0.951 |
| cptac     | epi5   | nuclei     | 0.598 | 0.415 | 0.777 |      0.951 |
| cptac     | att10  | neutrophil | 0.742 | 0.585 | 0.88  |      0.951 |
| cptac     | att10  | epithelial | 0.586 | 0.381 | 0.792 |      0.951 |
| cptac     | att10  | lymphocyte | 0.815 | 0.679 | 0.939 |      0.951 |
| cptac     | att10  | plasma     | 0.836 | 0.713 | 0.937 |      0.951 |
| cptac     | att10  | eosinophil | 0.589 | 0.411 | 0.769 |      0.951 |
| cptac     | att10  | connective | 0.645 | 0.47  | 0.815 |      0.951 |
| cptac     | att10  | iel        | 0.775 | 0.6   | 0.93  |      0.951 |
| cptac     | att10  | nuclei     | 0.67  | 0.466 | 0.857 |      0.951 |
| cptac     | perm10 | neutrophil | 0.56  | 0.394 | 0.722 |      0.951 |
| cptac     | perm10 | epithelial | 0.67  | 0.493 | 0.833 |      0.951 |
| cptac     | perm10 | lymphocyte | 0.475 | 0.308 | 0.647 |      0.951 |
| cptac     | perm10 | plasma     | 0.406 | 0.238 | 0.585 |      0.951 |
| cptac     | perm10 | eosinophil | 0.404 | 0.254 | 0.567 |      0.951 |
| cptac     | perm10 | connective | 0.403 | 0.243 | 0.573 |      0.951 |
| cptac     | perm10 | iel        | 0.754 | 0.586 | 0.903 |      0.951 |
| cptac     | perm10 | nuclei     | 0.61  | 0.424 | 0.787 |      0.951 |
| cptac     | epi10  | neutrophil | 0.743 | 0.587 | 0.881 |      0.951 |
| cptac     | epi10  | epithelial | 0.592 | 0.396 | 0.781 |      0.951 |
| cptac     | epi10  | lymphocyte | 0.676 | 0.521 | 0.816 |      0.951 |
| cptac     | epi10  | plasma     | 0.546 | 0.385 | 0.696 |      0.951 |
| cptac     | epi10  | eosinophil | 0.474 | 0.331 | 0.618 |      0.951 |
| cptac     | epi10  | connective | 0.517 | 0.353 | 0.682 |      0.951 |
| cptac     | epi10  | iel        | 0.763 | 0.592 | 0.909 |      0.951 |
| cptac     | epi10  | nuclei     | 0.612 | 0.42  | 0.795 |      0.951 |
| cptac     | att20  | neutrophil | 0.699 | 0.531 | 0.851 |      0.951 |
| cptac     | att20  | epithelial | 0.598 | 0.393 | 0.799 |      0.951 |
| cptac     | att20  | lymphocyte | 0.786 | 0.635 | 0.923 |      0.951 |
| cptac     | att20  | plasma     | 0.766 | 0.628 | 0.89  |      0.951 |
| cptac     | att20  | eosinophil | 0.603 | 0.447 | 0.754 |      0.951 |
| cptac     | att20  | connective | 0.6   | 0.424 | 0.778 |      0.951 |
| cptac     | att20  | iel        | 0.742 | 0.567 | 0.903 |      0.951 |
| cptac     | att20  | nuclei     | 0.676 | 0.478 | 0.862 |      0.951 |
| cptac     | perm20 | neutrophil | 0.567 | 0.399 | 0.732 |      0.951 |
| cptac     | perm20 | epithelial | 0.688 | 0.509 | 0.85  |      0.951 |
| cptac     | perm20 | lymphocyte | 0.442 | 0.277 | 0.613 |      0.951 |
| cptac     | perm20 | plasma     | 0.456 | 0.277 | 0.64  |      0.951 |
| cptac     | perm20 | eosinophil | 0.321 | 0.172 | 0.481 |      0.951 |
| cptac     | perm20 | connective | 0.35  | 0.199 | 0.512 |      0.951 |
| cptac     | perm20 | iel        | 0.721 | 0.539 | 0.888 |      0.951 |
| cptac     | perm20 | nuclei     | 0.603 | 0.418 | 0.78  |      0.951 |
| cptac     | epi20  | neutrophil | 0.705 | 0.531 | 0.859 |      0.951 |
| cptac     | epi20  | epithelial | 0.609 | 0.415 | 0.796 |      0.951 |
| cptac     | epi20  | lymphocyte | 0.641 | 0.482 | 0.789 |      0.951 |
| cptac     | epi20  | plasma     | 0.56  | 0.405 | 0.709 |      0.951 |
| cptac     | epi20  | eosinophil | 0.46  | 0.31  | 0.606 |      0.951 |
| cptac     | epi20  | connective | 0.539 | 0.363 | 0.713 |      0.951 |
| cptac     | epi20  | iel        | 0.736 | 0.554 | 0.896 |      0.951 |
| cptac     | epi20  | nuclei     | 0.626 | 0.434 | 0.81  |      0.951 |
| surgen    | all    | neutrophil | 0.52  | 0.394 | 0.639 |      0.88  |
| surgen    | all    | epithelial | 0.513 | 0.391 | 0.635 |      0.88  |
| surgen    | all    | lymphocyte | 0.562 | 0.437 | 0.683 |      0.88  |
| surgen    | all    | plasma     | 0.485 | 0.374 | 0.596 |      0.88  |
| surgen    | all    | eosinophil | 0.426 | 0.308 | 0.544 |      0.88  |
| surgen    | all    | connective | 0.366 | 0.25  | 0.49  |      0.88  |
| surgen    | all    | iel        | 0.624 | 0.493 | 0.748 |      0.88  |
| surgen    | all    | nuclei     | 0.497 | 0.369 | 0.621 |      0.88  |
| surgen    | att1   | neutrophil | 0.682 | 0.573 | 0.784 |      0.88  |
| surgen    | att1   | epithelial | 0.447 | 0.328 | 0.573 |      0.88  |
| surgen    | att1   | lymphocyte | 0.744 | 0.63  | 0.845 |      0.88  |
| surgen    | att1   | plasma     | 0.706 | 0.598 | 0.807 |      0.88  |
| surgen    | att1   | eosinophil | 0.618 | 0.515 | 0.719 |      0.88  |
| surgen    | att1   | connective | 0.757 | 0.647 | 0.856 |      0.88  |
| surgen    | att1   | iel        | 0.643 | 0.511 | 0.766 |      0.88  |
| surgen    | att1   | nuclei     | 0.603 | 0.481 | 0.72  |      0.88  |
| surgen    | perm1  | neutrophil | 0.497 | 0.384 | 0.609 |      0.88  |
| surgen    | perm1  | epithelial | 0.532 | 0.418 | 0.646 |      0.88  |
| surgen    | perm1  | lymphocyte | 0.376 | 0.255 | 0.499 |      0.88  |
| surgen    | perm1  | plasma     | 0.34  | 0.228 | 0.459 |      0.88  |
| surgen    | perm1  | eosinophil | 0.341 | 0.239 | 0.449 |      0.88  |
| surgen    | perm1  | connective | 0.331 | 0.225 | 0.439 |      0.88  |
| surgen    | perm1  | iel        | 0.583 | 0.467 | 0.699 |      0.88  |
| surgen    | perm1  | nuclei     | 0.413 | 0.284 | 0.549 |      0.88  |
| surgen    | epi1   | neutrophil | 0.544 | 0.446 | 0.642 |      0.88  |
| surgen    | epi1   | epithelial | 0.51  | 0.387 | 0.628 |      0.88  |
| surgen    | epi1   | lymphocyte | 0.594 | 0.477 | 0.708 |      0.88  |
| surgen    | epi1   | plasma     | 0.512 | 0.407 | 0.617 |      0.88  |
| surgen    | epi1   | eosinophil | 0.54  | 0.455 | 0.628 |      0.88  |
| surgen    | epi1   | connective | 0.56  | 0.45  | 0.674 |      0.88  |
| surgen    | epi1   | iel        | 0.589 | 0.464 | 0.712 |      0.88  |
| surgen    | epi1   | nuclei     | 0.515 | 0.388 | 0.636 |      0.88  |
| surgen    | att5   | neutrophil | 0.604 | 0.469 | 0.727 |      0.88  |
| surgen    | att5   | epithelial | 0.455 | 0.333 | 0.58  |      0.88  |
| surgen    | att5   | lymphocyte | 0.708 | 0.585 | 0.824 |      0.88  |
| surgen    | att5   | plasma     | 0.692 | 0.582 | 0.792 |      0.88  |
| surgen    | att5   | eosinophil | 0.62  | 0.507 | 0.724 |      0.88  |
| surgen    | att5   | connective | 0.71  | 0.59  | 0.817 |      0.88  |
| surgen    | att5   | iel        | 0.685 | 0.555 | 0.805 |      0.88  |
| surgen    | att5   | nuclei     | 0.619 | 0.495 | 0.736 |      0.88  |
| surgen    | perm5  | neutrophil | 0.486 | 0.366 | 0.606 |      0.88  |
| surgen    | perm5  | epithelial | 0.613 | 0.503 | 0.721 |      0.88  |
| surgen    | perm5  | lymphocyte | 0.462 | 0.341 | 0.583 |      0.88  |
| surgen    | perm5  | plasma     | 0.412 | 0.305 | 0.521 |      0.88  |
| surgen    | perm5  | eosinophil | 0.354 | 0.243 | 0.468 |      0.88  |
| surgen    | perm5  | connective | 0.295 | 0.19  | 0.405 |      0.88  |
| surgen    | perm5  | iel        | 0.641 | 0.518 | 0.756 |      0.88  |
| surgen    | perm5  | nuclei     | 0.477 | 0.351 | 0.607 |      0.88  |
| surgen    | epi5   | neutrophil | 0.587 | 0.473 | 0.695 |      0.88  |
| surgen    | epi5   | epithelial | 0.515 | 0.394 | 0.634 |      0.88  |
| surgen    | epi5   | lymphocyte | 0.547 | 0.434 | 0.663 |      0.88  |
| surgen    | epi5   | plasma     | 0.42  | 0.312 | 0.538 |      0.88  |
| surgen    | epi5   | eosinophil | 0.424 | 0.309 | 0.543 |      0.88  |
| surgen    | epi5   | connective | 0.472 | 0.353 | 0.595 |      0.88  |
| surgen    | epi5   | iel        | 0.6   | 0.476 | 0.721 |      0.88  |
| surgen    | epi5   | nuclei     | 0.506 | 0.38  | 0.63  |      0.88  |
| surgen    | att10  | neutrophil | 0.608 | 0.479 | 0.729 |      0.88  |
| surgen    | att10  | epithelial | 0.425 | 0.302 | 0.553 |      0.88  |
| surgen    | att10  | lymphocyte | 0.711 | 0.586 | 0.826 |      0.88  |
| surgen    | att10  | plasma     | 0.683 | 0.566 | 0.788 |      0.88  |
| surgen    | att10  | eosinophil | 0.579 | 0.456 | 0.695 |      0.88  |
| surgen    | att10  | connective | 0.69  | 0.575 | 0.794 |      0.88  |
| surgen    | att10  | iel        | 0.675 | 0.542 | 0.797 |      0.88  |
| surgen    | att10  | nuclei     | 0.573 | 0.449 | 0.699 |      0.88  |
| surgen    | perm10 | neutrophil | 0.508 | 0.388 | 0.623 |      0.88  |
| surgen    | perm10 | epithelial | 0.614 | 0.501 | 0.724 |      0.88  |
| surgen    | perm10 | lymphocyte | 0.485 | 0.361 | 0.609 |      0.88  |
| surgen    | perm10 | plasma     | 0.424 | 0.316 | 0.534 |      0.88  |
| surgen    | perm10 | eosinophil | 0.318 | 0.21  | 0.431 |      0.88  |
| surgen    | perm10 | connective | 0.282 | 0.181 | 0.389 |      0.88  |
| surgen    | perm10 | iel        | 0.626 | 0.503 | 0.741 |      0.88  |
| surgen    | perm10 | nuclei     | 0.504 | 0.377 | 0.632 |      0.88  |
| surgen    | epi10  | neutrophil | 0.542 | 0.421 | 0.653 |      0.88  |
| surgen    | epi10  | epithelial | 0.516 | 0.393 | 0.635 |      0.88  |
| surgen    | epi10  | lymphocyte | 0.544 | 0.428 | 0.662 |      0.88  |
| surgen    | epi10  | plasma     | 0.437 | 0.328 | 0.553 |      0.88  |
| surgen    | epi10  | eosinophil | 0.406 | 0.292 | 0.524 |      0.88  |
| surgen    | epi10  | connective | 0.453 | 0.338 | 0.578 |      0.88  |
| surgen    | epi10  | iel        | 0.61  | 0.486 | 0.729 |      0.88  |
| surgen    | epi10  | nuclei     | 0.504 | 0.378 | 0.629 |      0.88  |
| surgen    | att20  | neutrophil | 0.626 | 0.503 | 0.738 |      0.88  |
| surgen    | att20  | epithelial | 0.429 | 0.308 | 0.554 |      0.88  |
| surgen    | att20  | lymphocyte | 0.704 | 0.583 | 0.815 |      0.88  |
| surgen    | att20  | plasma     | 0.652 | 0.535 | 0.757 |      0.88  |
| surgen    | att20  | eosinophil | 0.547 | 0.421 | 0.666 |      0.88  |
| surgen    | att20  | connective | 0.646 | 0.532 | 0.751 |      0.88  |
| surgen    | att20  | iel        | 0.665 | 0.534 | 0.786 |      0.88  |
| surgen    | att20  | nuclei     | 0.534 | 0.41  | 0.658 |      0.88  |
| surgen    | perm20 | neutrophil | 0.535 | 0.416 | 0.652 |      0.88  |
| surgen    | perm20 | epithelial | 0.593 | 0.476 | 0.708 |      0.88  |
| surgen    | perm20 | lymphocyte | 0.485 | 0.356 | 0.613 |      0.88  |
| surgen    | perm20 | plasma     | 0.44  | 0.33  | 0.551 |      0.88  |
| surgen    | perm20 | eosinophil | 0.351 | 0.239 | 0.466 |      0.88  |
| surgen    | perm20 | connective | 0.289 | 0.182 | 0.406 |      0.88  |
| surgen    | perm20 | iel        | 0.602 | 0.476 | 0.723 |      0.88  |
| surgen    | perm20 | nuclei     | 0.502 | 0.377 | 0.631 |      0.88  |
| surgen    | epi20  | neutrophil | 0.523 | 0.401 | 0.637 |      0.88  |
| surgen    | epi20  | epithelial | 0.519 | 0.397 | 0.638 |      0.88  |
| surgen    | epi20  | lymphocyte | 0.549 | 0.434 | 0.667 |      0.88  |
| surgen    | epi20  | plasma     | 0.418 | 0.31  | 0.531 |      0.88  |
| surgen    | epi20  | eosinophil | 0.41  | 0.3   | 0.526 |      0.88  |
| surgen    | epi20  | connective | 0.427 | 0.314 | 0.546 |      0.88  |
| surgen    | epi20  | iel        | 0.619 | 0.49  | 0.743 |      0.88  |
| surgen    | epi20  | nuclei     | 0.504 | 0.378 | 0.631 |      0.88  |
