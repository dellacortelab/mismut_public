# Literature check: pixel size of the Kather MSI tiles (Zenodo 2530835)

Searched 2026-10-16. Sources searched:
- the web
- Europe PMC
- about 25 arXiv full texts
- GitHub: the jnkather/MSIfromHE issues, plus the THUNDER, kaiko eva and Patho-Bench repos

Not checked: Kaggle discussions, MDPI pages (blocked) and X.

## Prior reports that 0.5 um/px is wrong

None found. The record still says "224 px x 224 px at a resolution of 0.5 µm/px".

## Sources consistent with 256 um at 224 px (1.143 um/px)

- **Kather et al. 2019, Nat Med, Methods:** "unique image tiles of 256 μm edge length"; Fig. 1 scale bar 256 μm.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7423299/
- **Kather lab STAMP protocol** (El Nahhas et al., arXiv 2312.10944): "The default MPP is 256/224 ≈ 1.14 which
  corresponds to ~9x magnification." https://arxiv.org/abs/2312.10944
- **DeepSMILE** (Schirris et al., arXiv 2107.09405): own slides "tiled at a spacing of 1.14 microns per pixel (mpp)
  with a tile size of 224 × 224px, resulting in an edge length of 256 µm".
- **Kather's later dataset, Zenodo 3832231:** "tiles of 256 µm edge length, saved as 512 px images (effective
  magnification 0.5 µm/px)". This is a different dataset (51,918 tiles).
- **This repo:** registration of the 2530835 tiles to the original TCGA slides (`src/check_resolution.py`,
  `src/retile.py`).

## Works using 2530835 that state 0.5 um/px or 112 um

| Work | Statement | URL |
|---|---|---|
| EXAONEPath 1.0 (foundation model, 2024) | "All images are 224² pixels at 0.5 MPP"; MSI-CRC benchmark | https://arxiv.org/abs/2408.00380 |
| Hibou (foundation model, 2024) | "MSI-CRC ... (224×224 pixels, 0.5 µm/px)"; also MSI-STAD | https://arxiv.org/abs/2406.05074 |
| Wang et al., TEMI, PLoS Comput Biol 2026 | "224 pixel × 224 pixel at a resolution of 0.5 μm/pixel" | https://doi.org/10.1371/journal.pcbi.1013950 |
| Swanson et al., bioRxiv 2025 | "224 × 224-pixel tiles (i.e. tiles that are 112μm x 112μm), ... 0.5 μm/px" | https://www.biorxiv.org/content/10.1101/2025.03.30.646246v1.full |
| Vardar & Vardar, Diagnostics 2026 (target paper) | "Source spatial resolution 0.5 µm/pixel" | doi:10.3390/diagnostics16203271 |

## Other uses of 2530835 that give no pixel size

- DeepSMILE ("TCGA-CRCk")
- Saillard et al. 2021
- Schwen et al. 2021
- Ke et al. 2025, Sci Rep: Kaggle mirror, 99.44% accuracy, no patient-level split stated

## Large foundation-model benchmarks

UNI, THUNDER and GPFM use the 512 px dataset (3832231, about 0.5 um/px) rather than 2530835. kaiko eva and
Patho-Bench contain neither dataset.

## Patch-count confound

Apart from Vardar & Vardar, no other report found. Kather documents the undersampling ("The number of image tiles
per class was equalized by undersampling") but does not discuss it as a confound.

## Contact

J. N. Kather:
- **Since 1 Oct 2026:** Hopp Foundation Professor and founding director of the Institute for Medical AI (IMAI),
  Heidelberg University Hospital. jakob.kather@med.uni-heidelberg.de (idw press release, 16 Sep 2026).
- **Before that:** TU Dresden / EKFZ.
- **GitHub:** https://github.com/KatherLab and https://github.com/jnkather/MSIfromHE (the latter has open issues
  with replies from him).

## Side note

The Zenodo description lists CRC_DX_TEST_MSIMUT as 29,335 patches; the archive contains 28,335 (our index; also
noted in AI-HHMI/mia-agentic-search PR #1399).
