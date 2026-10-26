# Zero-shot: prostate patch-cancer head P1 on colorectal tissue

## NCT CRC-VAL-HE-7K (labelled colon tissue classes)

| Variant | Class | n | mean P1 | tokens > 0.5 | tiles with >10% cancer tokens | AUC TUM vs class |
|---|---|---|---|---|---|---|
| v1 | ADI | 1338 | 0.072 | 0.009 | 0.000 | 0.998 |
| v1 | BACK | 847 | 0.034 | 0.000 | 0.000 | 0.999 |
| v1 | DEB | 339 | 0.144 | 0.096 | 0.274 | 0.994 |
| v1 | LYM | 634 | 0.109 | 0.011 | 0.014 | 0.998 |
| v1 | MUC | 1035 | 0.028 | 0.007 | 0.015 | 0.999 |
| v1 | MUS | 592 | 0.066 | 0.002 | 0.007 | 0.998 |
| v1 | NORM | 741 | 0.497 | 0.550 | 0.989 | 0.958 |
| v1 | STR | 421 | 0.079 | 0.006 | 0.014 | 0.998 |
| v1 | TUM | 1233 | 0.811 | 0.919 | 0.996 |  |
| flash | ADI | 1338 | 0.034 | 0.000 | 0.000 | 1.000 |
| flash | BACK | 847 | 0.038 | 0.000 | 0.000 | 1.000 |
| flash | DEB | 339 | 0.246 | 0.195 | 0.431 | 0.972 |
| flash | LYM | 634 | 0.308 | 0.163 | 0.522 | 0.979 |
| flash | MUC | 1035 | 0.032 | 0.011 | 0.033 | 0.998 |
| flash | MUS | 592 | 0.115 | 0.008 | 0.020 | 0.999 |
| flash | NORM | 741 | 0.470 | 0.484 | 0.962 | 0.920 |
| flash | STR | 421 | 0.130 | 0.019 | 0.055 | 0.997 |
| flash | TUM | 1233 | 0.800 | 0.882 | 0.991 |  |

## TCGA: Kather tumor tiles vs tissue > 300 µm away (20944 tiles, 40 slides)

| Scorer | AUC | mean tumor | mean non-tumor |
|---|---|---|---|
| v1 | 0.879 | 0.632 | 0.187 |
| flash | 0.876 | 0.616 | 0.197 |
| in-domain tumor classifier (cross-fitted) | 0.893 | nan | nan |
