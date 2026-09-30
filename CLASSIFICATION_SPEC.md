# Classification spec: spatial immune phenotype and biomarker pattern

**Spec version 1.0.0.** Machine-readable parts: `classification_config.json` (thresholds), `classification.schema.json` (output). Reference implementation: `classify_reference.py`. Conformance tests: `test_classify_reference.py`.

This document is the contract between the spatial statistics skill(s) and the orchestrator. The skill computes; the orchestrator and agents only read `classification.json` and never recompute or reinterpret a label.

## 0. Pre-registration rule

The thresholds in `classification_config.json` v1.0.0 are frozen **before** the real dataset is analysed. They were calibrated only on the synthetic suite (section 10). If a threshold is changed after seeing real data, bump `spec_version`, rerun the synthetic suite, and report both the original and the new result. This is what lets you say to judges that the rule was not tuned to produce the answer.

## 1. Inputs

The skill accepts one of:

1. **Space Ranger v4 segmented output** (the target for real data). Read from `outs/segmented_outputs/`:
   - `filtered_feature_cell_matrix.h5` (cells × genes counts)
   - `cell_segmentations.geojson` (cell polygons; centroid = polygon centroid)
   - `spatial/scalefactors_json.json` (pixel to micrometre conversion)

   GeoJSON `cell_id: 1` corresponds to matrix barcode `cellid_000000001-1`; the loader must map these explicitly and fail if any matrix cell has no polygon. Polygon coordinates are in full-resolution image pixels; convert to micrometres using the microns-per-pixel value in `scalefactors_json.json` (**verify the key name on the downloaded file before relying on it**). The loader writes an `.h5ad` with `obsm["spatial"]` in micrometres and `uns["spatial_units"] = "um"`, plus the input checksums.
2. **`.h5ad`** with raw integer counts in `X`, gene symbols as `var_names`, `obsm["spatial"]` in micrometres, and `uns["spatial_units"] == "um"`.
3. **Synthetic tissue directory** from `synthetic_tissue.py` (`counts.npz`, `genes.txt`, `cells.csv`). `truth_cells.csv` and `ground_truth.json` must never be read by the skill.

Also required: `seed` (integer; logged) and optionally a config override file (the fully resolved config is written into the output).

All distances in this spec are micrometres.

## 2. QC

| Parameter | v1.0.0 |
|---|---|
| `qc.min_umi` | 20 |
| `qc.min_genes` | 5 |
| `qc.min_pass_fraction` | 0.5 |

A cell is dropped if total UMI < `min_umi` (reason `low_umi`), else if detected genes < `min_genes` (reason `low_genes`). Every drop is counted by reason; nothing is filtered silently. If the passing fraction is below `min_pass_fraction`, **both** calls are `unclassifiable` with reason `qc_fail`.

`min_genes` = 5 suits the 23-gene synthetic panel. For whole-transcriptome Visium HD data, inspect the distribution first and, if needed, change it **before** looking at any phenotype metric, recording the change as in section 0.

### 2.1 Optional spatially-local outlier QC (`spatial_qc`)

Absent from the frozen v1.0.0 config, so behaviour above is unchanged unless a config override adds a `spatial_qc: {"k": int, "z_threshold": float}` section. When present, `qc()` additionally flags a cell as an outlier relative to its `k` spatially-nearest neighbours (not a single global threshold across the section) for three metrics: log-transformed library size, log-transformed gene count, and mitochondrial fraction (untransformed) — computed only if the panel has `MT-`-prefixed genes. A cell is flagged if its value is `z_threshold` local standard deviations below (library size, genes) or above (mitochondrial fraction) its neighbourhood's mean, using `graph.max_edge_um` for the neighbour search. Reasons: `spatial_outlier_low_libsize`, `spatial_outlier_low_genes`, `spatial_outlier_high_mito` (mutually exclusive with each other and with `low_umi`/`low_genes` — a cell already dropped globally is not double-counted). A cell with fewer than one valid neighbour cannot be locally evaluated and is not flagged by this mechanism.

Rationale: in spatial data, these metrics can shift for genuine biological reasons across a tissue section (e.g. a naturally lower-RNA-content region), which a single global cutoff cannot distinguish from technical failure. This follows the local-outlier-detection approach used by `SpotSweeper` (Totty, Hicks & Guo, 2025) and described in the *Orchestrating Spatial Transcriptomics Analysis with Bioconductor* (OSTA) book's QC chapter, adapted here to numpy/scipy without external spatial-stats dependencies.

## 3. Normalisation and cell typing

1. Normalise each cell to 10,000 counts and take `log1p`.
2. For every signature gene, scale by its `cell_typing.scale_quantile` (0.99) quantile across passing cells and clip to [0, 1].
3. Signature score = mean of its scaled genes:

| Signature | Genes |
|---|---|
| `tumor` | EPCAM, KRT8, KRT18, KRT19 |
| `t_cell` | PTPRC, CD3E, CD8A, CD8B, GZMB |
| `stroma` | COL1A1, COL1A2, DCN, PDGFRB |

4. Cell type = argmax signature if its score ≥ `min_score` (0.2) and it beats the runner-up by ≥ `min_margin` (0.1); otherwise `ambiguous`.

**The biomarker gene must not be in the tumor signature.** The loader rejects such a config. If ERBB2 defined tumour cells, ERBB2-low tumour cells would be missed and the heterogeneity test would be circular.

`t_cell` targets cytotoxic T cells because the inflamed/excluded/desert framework is defined by T cell localisation. `ambiguous` cells stay in the graph but count as neither tumour nor T cell.

## 4. Signature concordance

For each signature and each marker, compute the Spearman correlation between the marker and the leave-one-out signature score across passing cells. A signature is **trusted** if the median of these correlations ≥ `concordance.min_median_loo_spearman` (0.2). Values are always reported.

The tumour signature must be trusted for any immune call. The T cell signature must be trusted for `inflamed`/`excluded`, but **not** for `desert`: with almost no T cells, T markers cannot correlate (the synthetic desert tissue gives a median of 0.04), so requiring it would make desert impossible to call.

## 5. Spatial graph

kNN graph on cell centroids, `graph.k` = 6 neighbours, edges longer than `graph.max_edge_um` = 50 µm removed (prevents linking across tissue gaps or folds). Directed edges; Moran's I uses row-standardised weights over the remaining neighbours.

## 6. Immune phenotype

### Metrics

- `t_cell_fraction` = T cells / passing cells.
- **Tumour region**: a cell (of any type) whose ≥ `region_min_tumor_neighbour_fraction` (0.5) of its `region_k` (15) nearest neighbours are tumour cells.
- `infiltration_fraction` = T cells in tumour regions / all T cells (reported, descriptive).
- **`infiltration_ratio`** (primary) = (T cells / all cells, within tumour regions) ÷ (T cells / all cells, whole section). 1 means T cells are as dense inside tumour as elsewhere; 0 means none inside.
- `nhood_z`, `nhood_p` (supporting): tumour-T edge count vs 1000 permutations of cell-type labels across passing cells. z = (observed − null mean) / null SD; two-sided empirical p = (1 + #|null − mean| ≥ |obs − mean|) / (n + 1). Minimum achievable p with 1000 permutations is 1/1001.

**Why the ratio is primary and z is supporting:** z depends on how many T cells sit in stroma far from any tumour, so even an inflamed tumour can give negative z. In the synthetic spillover tissue, a truly inflamed section gave z = −7.4 while the ratio stayed at 0.91. The ratio compares T cell density inside vs outside tumour directly.

### Decision rule (evaluated in order; first match wins)

| # | Condition | Label | reason |
|---|---|---|---|
| 1 | pass fraction < 0.5 | unclassifiable | `qc_fail` |
| 2 | tumour cells < `min_tumor_cells` (200) | unclassifiable | `insufficient_tumor_cells` |
| 3 | tumour signature not trusted | unclassifiable | `tumor_signature_discordant` |
| 4 | `t_cell_fraction` < `desert_max_t_fraction` (0.01) | **desert** | |
| 5 | T cells < `min_t_cells` (50) | unclassifiable | `insufficient_t_cells` |
| 6 | T cell signature not trusted | unclassifiable | `t_cell_signature_discordant` |
| 7 | no tumour-region cells | unclassifiable | `no_tumor_region` |
| 8 | `infiltration_ratio` ≥ `inflamed_min_ratio` (0.5) | **inflamed** | |
| 9 | `infiltration_ratio` ≤ `excluded_max_ratio` (0.2) **and** `nhood_z` ≤ −2 **and** `nhood_p` ≤ 0.01 | **excluded** | |
| 10 | otherwise | unclassifiable | `borderline` |

Each step taken is appended to `rule_path` as a human-readable string with the actual values, so the agent can explain the call by quoting the path rather than reasoning about it.

The gap between 0.2 and 0.5 is a deliberate abstention zone. The thresholds are heuristic: clinical immunophenotyping is usually done on CD8 IHC in epithelium vs stroma, and there is no single agreed cutoff. The output must describe these as rule-based categories, not clinical immunophenotypes.

## 7. Biomarker pattern (default ERBB2)

Computed **among tumour cells only**.

- `detection_fraction` = tumour cells with ≥ 1 ERBB2 count.
- `morans_i` on log-normalised ERBB2 over a kNN graph built on tumour cells (same k and max edge); one-sided permutation p for positive autocorrelation (1000 permutations of values).
- `morans_i_depth_control` = Moran's I of log total counts among tumour cells (technical spatial structure such as depth gradients).
- **`morans_i_excess`** = `morans_i` − max(`morans_i_depth_control`, 0).

| # | Condition | Label | reason |
|---|---|---|---|
| 1 | pass fraction < 0.5 | unclassifiable | `qc_fail` |
| 2 | tumour cells < 200 | unclassifiable | `insufficient_tumor_cells` |
| 3 | `detection_fraction` < 0.05 | unclassifiable | `biomarker_not_expressed` |
| 4 | `morans_i_excess` ≥ 0.10 **and** `morans_p` ≤ 0.01 | **heterogeneous** | |
| 5 | `morans_i_excess` < 0.05 | **homogeneous** | |
| 6 | otherwise | unclassifiable | `borderline` |

The decision uses effect size plus significance, never p alone: with thousands of cells, tiny autocorrelations reach significance.

This is **ERBB2 transcript** heterogeneity. Clinical HER2 status and heterogeneity are defined by IHC/ISH (protein and gene amplification). The final answer must say "ERBB2 transcript" and must not state or imply a clinical HER2 score.

## 8. Warnings (do not change labels; must be reported)

- `spillover_suspected`: among T-like cells (T score ≥ 0.5), the rate of also scoring as tumour (tumour score ≥ 0.2) in cells with a tumour neighbour, divided by the rate in cells without one, exceeds `max_spillover_ratio` (2.0), given ≥ 30 cells in each group. Spillover at segmentation boundaries only affects cells touching tumour, so the ratio sits near 1 without it (0.5 to 1.1 on clean synthetic tissues) and far above with it (12.5).
- `nhood_z_conflict`: label is inflamed but `nhood_z` ≤ −2.
- `depth_autocorrelation`: `morans_i_depth_control` > 0.10.

## 9. Output

`outputs/classification.json`, validated against `classification.schema.json`. Invariants the schema enforces:

- `label == "unclassifiable"` ⇔ `reason` is a string; otherwise `reason` is null.
- `rule_path` is non-empty for both calls.
- The fully resolved `config`, `seed`, `spec_version`, and input `sha256` are always present.

Field names used by `faithfulness_check.py` (do not rename without updating it): `immune_phenotype.label`, `.nhood_z`, `.nhood_p`; `biomarker.label`, `.gene`, `.morans_i`, `.morans_p`; `qc.n_cells_input`, `.n_cells_passed`, `.n_dropped`.

**The orchestrator's obligations:** copy both labels verbatim into `final_answer.json`; if a label is `unclassifiable`, set `abstention_reason` from `reason` and do not offer a phenotype anyway; report every warning; report effect sizes alongside p-values.

## 10. Calibration evidence (synthetic suite, 20,000 cells, seed 0, 1000 permutations)

| Tissue | t_cell_fraction | infiltration_ratio | nhood_z | Immune label | ERBB2 I_excess | ERBB2 label | spillover ratio |
|---|---|---|---|---|---|---|---|
| desert_homogeneous | 0.004 | 1.03 | −0.1 | desert | 0.00 | homogeneous | 1.06 |
| excluded_heterogeneous | 0.124 | 0.11 | −76.9 | excluded | 0.49 | heterogeneous | 0.85 |
| excluded_homogeneous | 0.125 | 0.10 | −78.1 | excluded | −0.01 | homogeneous | 0.93 |
| inflamed_heterogeneous | 0.141 | 1.43 | 35.8 | inflamed | 0.49 | heterogeneous | 1.00 |
| inflamed_homogeneous | 0.141 | 1.43 | 35.8 | inflamed | 0.02 | homogeneous | 0.86 |
| inflamed_homogeneous_lowqc | 0.142 | n/a | n/a | unclassifiable (qc_fail) | n/a | unclassifiable (qc_fail) | 0.52 |
| inflamed + spillover 0.3 | 0.100 | 0.91 | −7.4 | inflamed (+2 warnings) | 0.02 | homogeneous | 12.47 |

All seven match ground truth. These tissues are deliberately clear-cut; they prove the plumbing and the rule logic, not that the thresholds are right for real tissue. Section 11 lists what that means for claims.

## 11. Known limitations (state these in the demo)

- Thresholds are heuristic and calibrated on synthetic data only.
- Cell typing uses a handful of markers; no reference-based annotation. Sparse markers (CD8A) at cell level can under-call T cells, pushing borderline tissue toward desert.
- `desert` is only as good as T marker detection: a technical failure of T markers would also look like desert. The T signature concordance is reported so a reader can judge this.
- One tissue section, one timepoint; no claim about the rest of the tumour or about outcome.
- ERBB2 RNA is a proxy for, not a measurement of, clinical HER2.

## 12. Change control

Any change to thresholds, signatures, graph parameters, or rule order is a spec change: bump `spec_version` in this file, the config, the schema, and `classify_reference.py`, rerun `pytest`, and update section 10.
