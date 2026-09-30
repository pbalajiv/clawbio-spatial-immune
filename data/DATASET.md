# Dataset

## Source

Visium HD Spatial Gene Expression Library, Human Breast Cancer (Fresh Frozen), Ultima Sequencing.

- URL: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-human-breast-cancer-ff-ultima-4
- Date accessed: 2026-09-25
- Licence: Creative Commons Attribution 4.0 International (CC BY 4.0). Cite per 10x Genomics' citation guidelines linked from the dataset page.
- Product / software: HD Spatial Gene Expression, processed with Space Ranger v4.0.1 (`spaceranger count`).

## Sample description

- One donor, one sample (slide H1-649ZJJX, area A1).
- Fresh frozen human female breast tissue, ductal carcinoma in situ (DCIS), RIN 7.7, obtained from BioIVT.
- 10 µm cryostat section, H&E stained and imaged (Olympus VS200 Slide Scanner, 20x, NA 0.8).
- Probe set: Visium Human Transcriptome Probe Set v2.0. Sequenced on Ultima UG 100 (2.02 billion reads, downsampled), converted for Ultima compatibility.
- 10x notes local freezing artifacts in parts of the H&E image that impair cell segmentation performance in those areas; expect to see this in QC (T3).

## HER2 status

**Unknown.** The dataset page documents disease state as "cancer" / DCIS but does not report ER/PR/HER2 receptor status. Per `DATA.md`, `sample_subtype` in the final answer is set to `"unknown"`. Any literature cited for HER2+ phenotypes must carry an explicit mismatch/unknown acknowledgment (see `CLASSIFICATION_SPEC.md` section 11 and `DATA.md` recommendation 2).

## DCIS caveat

This is an in-situ (not invasive) lesion: ducts are bounded by a myoepithelial layer and basement membrane, so T cells sitting outside ducts may look like an "excluded" spatial pattern for purely anatomical reasons rather than a biological immune-exclusion mechanism. Most immune-phenotype literature (`curated_references.json`, T7) describes invasive tumours. State this caveat in the final answer whenever the immune phenotype call is `excluded`.

## Files downloaded (via `scripts/download_data.sh`)

Only the files this project's spec (`CLASSIFICATION_SPEC.md` section 1) needs were fetched — not the full batch-download bundle (which also includes ~20 GB of binned outputs and multi-GB Loupe/molecule-info files unused here).

| File | Type | Size (per dataset page) | Role |
|---|---|---|---|
| `Visium_HD_FF_Human_Breast_Cancer_segmented_outputs.tar.gz` | GZ | 2.44 GB | Contains `filtered_feature_cell_matrix.h5` and `cell_segmentations.geojson` (loader input, T2) |
| `Visium_HD_FF_Human_Breast_Cancer_spatial.tar.gz` | GZ | 100 MB | Contains `scalefactors_json.json` and the hires H&E image (loader input + T2 overlay check) |
| `Visium_HD_FF_Human_Breast_Cancer_metrics_summary.csv` | CSV | 2.05 kB | Space Ranger run metrics, for reference |

Total download: ~2.5 GB (under the ~10 GB ask-first threshold in `CLAUDE.md`).

md5sums as published on the dataset page (for the downloader's own cross-check; this project's checksums of record are the sha256 hashes written by the download script to `data/CHECKSUMS.txt`):

- `segmented_outputs.tar.gz`: `3f633d8f008ff90a7ed9f40bf0b09ef4`
- `spatial.tar.gz`: `8fea53dbe80ace55f0fda8878de2c0b5`
- `metrics_summary.csv`: `4200fa254f3360d0d005ec63b718d3a3`

## T3 pre-flight (before any phenotype metric)

Run via `scripts/qc_preflight.py --h5ad data/tissue.h5ad --out data/qc_preflight.png`. Both findings below are dataset-specific, not phenotype-metric-driven, and are applied as a **config override** (`data/real_data_config.json`, deep-merged over `classification_config.json`), so the frozen v1.0.0 config used for synthetic calibration is untouched (`CLASSIFICATION_SPEC.md` section 0). `classification.json`'s `config` field always logs the fully resolved config actually used for a given run, so this override is fully traceable per run without needing to bump `spec_version`.

### Signature gene presence

All tumor (`EPCAM, KRT8, KRT18, KRT19`), stroma (`COL1A1, COL1A2, DCN, PDGFRB`) and biomarker (`ERBB2`) genes are present.

**`CD8B` is absent from this dataset's probe panel** (Visium Human Transcriptome Probe Set v2.1.0, per the matrix h5's `target_sets`). Verified by Ensembl ID (`ENSG00000172116`), not a symbol/alias mismatch — the panel does contain a different gene, `CD8B2` (`ENSG00000254126`), a distinct, poorly-characterized paralog with no established co-expression with cytotoxic T cells, so it is not used as a substitute.

**Decision:** the `t_cell` signature for this dataset is `{PTPRC, CD3E, CD8A, GZMB}` (4 genes, `CD8B` dropped). Reason: forced by probe-panel gene availability, not tuned to any phenotype outcome — decided before computing any spatial/phenotype metric.

### UMI / genes-per-cell distribution

n = 38,993 cells. Median 6,637 UMI / 2,728 genes per cell; unimodal distribution, no distinct bimodal low-quality population, but the frozen `qc.min_umi=20` / `qc.min_genes=5` (calibrated for the 23-gene synthetic panel) passes essentially 100% of cells on this ~18,000-gene panel, i.e. no real filtering.

| Threshold | Fraction of cells passing |
|---|---|
| UMI ≥ 20 | 1.000 |
| UMI ≥ 500 | 0.963 |
| genes ≥ 5 | 1.000 |
| genes ≥ 200 | 0.985 |
| genes ≥ 500 | 0.947 |

**Decision:** `qc.min_umi = 500`, `qc.min_genes = 200` for this dataset. Reason: conventional scRNA-seq-style cutoffs for whole-transcriptome data, drop the bottom ~4-6% tail of low-complexity cells, consistent with 10x's own note of local freezing-artifact segmentation issues in parts of this section. Chosen from the QC distribution alone, before looking at any immune-phenotype or ERBB2 metric.

Both decisions are recorded here per `CLASSIFICATION_SPEC.md` section 0's pre-registration rule, before classification was run on this data.

### Additional QC diagnostics (reporting only — not yet applied as gates)

Added after an initial review flagged these as gaps: mitochondrial fraction, segmentation cell-area, per-gene detection rate (a gene can be present but barely detected), and a lightweight doublet-candidate heuristic. None of these currently drop cells or change `qc.n_cells_passed` — they are computed and plotted by `scripts/qc_preflight.py` for review, exactly like `min_umi`/`min_genes` were before becoming a gate. Whether/where to turn any of these into a blocking filter is a separate decision, to be made from these distributions alone, not from any phenotype outcome.

- **Per-gene detection rate**: every signature/biomarker gene has healthy detection (24–90%among passing-QC cells; lowest is `EPCAM` at 24.4%). No additional weak-but-present gene beyond `CD8B`'s outright absence.
- **Mitochondrial fraction** (12 `MT-` genes found in this panel): median 1.29%, p95 7.2%, p99 9.2%, max 22.6%. Mostly low and unimodal, with a small secondary bump around 5–8% — a candidate stressed/dying subpopulation worth a closer look before any gate is set, but not an alarming overall picture.
- **Segmented cell area**: median 318.8 µm², p1 51.8 µm², p99 722.7 µm², max 1299.6 µm². Unimodal, no obvious bimodal "clearly-merged-cell" population, but the upper tail (largest ~1%) is a plausible place to look for segmentation-merged cells.
- **Doublet-candidate heuristic** (cells exceeding 2× the cohort median in *both* total UMI and genes detected — see `scripts/qc_preflight.py`'s docstring for why this crude proxy was chosen over a trained classifier like Scrublet): **5,457 of 38,993 cells (14.0%) flagged.** This is almost certainly an overestimate of true doublets — large, highly transcriptionally active real tumour cells will also clear this bar — but 14% is high enough that it should not be waved away either; it's a signal that the segmentation-merge concern (gap #4 in `qa_notes/gaps_and_independent_checklist.md`) is worth taking seriously before this pipeline's outputs are treated as final.

Figure: `data/qc_preflight.png` (2×3 panel: UMI, genes, mitochondrial %, cell area, and the doublet scatter).

### Spatially-local outlier QC (gate, now active)

Per literature review (OSTA book's QC chapter; `SpotSweeper`, Totty/Hicks/Guo 2025), a single global cutoff on library size/gene count/mitochondrial fraction risks conflating genuine tissue biology (a naturally lower-RNA-content region) with technical failure — spatial data specifically needs outlier detection relative to each cell's local neighbourhood, not the whole section. Implemented as `local_outlier_qc()` in `classify_reference.py` (see `CLASSIFICATION_SPEC.md` section 2.1) and enabled for this dataset via `data/real_data_config.json`'s `spatial_qc: {"k": 30, "z_threshold": 3.0}`.

**Result on this dataset** (layered on top of the existing `min_umi=500`/`min_genes=200` global floor, applied only to cells that already passed it): an additional **1,962 cells (5.0%)** dropped — 786 for anomalously low local library size, 346 for anomalously low local gene count, 830 for anomalously high local mitochondrial fraction. Combined with the global floor: **35,576 of 38,993 cells pass (91.2%)**. The classification result is unchanged (`excluded` / `homogeneous`), i.e. this additional QC layer removed more low-quality cells without flipping the phenotype call — a mild robustness check in itself.

Doublet-candidate heuristic and cell-area remain reporting-only; `SpotSweeper`'s method targets library size, gene count, and mitochondrial fraction specifically, not segmentation-merge detection.
