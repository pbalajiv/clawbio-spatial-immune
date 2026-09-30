# CLAUDE.md

Project context for Claude Code. Read this, then `CLASSIFICATION_SPEC.md`, before changing anything.

## Goal

ClawBio hackathon project. Answer, for one breast cancer tissue section: is it immune-inflamed, immune-excluded or immune-desert; is ERBB2 transcript expression spatially homogeneous or heterogeneous among tumour cells; and what does curated literature associate with each category. Every number must be traceable to a deterministic skill output, and the system must abstain when the data cannot support a call. Output is pattern classification against known categories, **not** a per-patient outcome prediction.

## Architecture

- **ClawBio skills (deterministic, no LLM):** load + QC, cell scoring, spatial statistics + classification (per `CLASSIFICATION_SPEC.md`), plots/report.
- **Agents:** orchestrator (calls skills, writes `final_answer.json`), literature agent (uses only `curated_references.json`), critic (= `faithfulness_check.py`, plain code).
- **Orchestration wrapper** (plain code, not the LLM) writes `run_log.jsonl` with skill name, version, seed, and input/output SHA-256.

## Files

| File | Role |
|---|---|
| `CLASSIFICATION_SPEC.md` | The rule. Source of truth. |
| `classification_config.json` | Frozen v1.0.0 thresholds |
| `classification.schema.json` | Output contract for `outputs/classification.json` |
| `classify_reference.py` | Reference implementation of the spec (numpy/scipy) |
| `synthetic_tissue.py` | Ground-truth tissues (inflamed/excluded/desert × ERBB2 pattern, low-QC, spillover) |
| `faithfulness_check.py` | Deterministic checks on the agent's final answer + stability across runs |
| `make_examples.py` | Example good/bad runs, curated-reference template |
| `test_harness.py`, `test_classify_reference.py` | 32 tests; must stay green |
| `DATA.md` | Dataset options and the default choice |

## Commands

```bash
pip install -r requirements.txt
pytest -q                                              # all tests
python synthetic_tissue.py --suite --out synthetic/    # ground-truth tissues
python classify_reference.py --input synthetic/excluded_heterogeneous --out runs/demo/outputs/classification.json
python faithfulness_check.py check --run-dir runs/demo --refs curated_references.json --truth synthetic/excluded_heterogeneous/ground_truth.json
```

## Invariants (do not violate without asking the user)

1. Thresholds, signatures, graph parameters and rule order live in the config and spec. Changing any of them is a spec change: bump `spec_version` everywhere, rerun `pytest`, update spec section 10. Never tune thresholds after looking at real-data phenotype metrics.
2. The biomarker gene is never part of the tumour signature (the loader enforces this).
3. Skills never read `truth_cells.csv` or `ground_truth.json`.
4. The LLM never computes, rounds creatively, or recalls a statistic; it copies from skill JSON. Labels are copied verbatim; `unclassifiable` stays unclassifiable.
5. Every dropped cell is counted by reason. No silent filtering.
6. Any ClawBio skill implementation must reproduce `classify_reference.py` labels on the synthetic suite (reuse `test_classify_reference.py` against it).
7. Cite references by `REF-###` only; entries stay `verified: false` until a human has read the paper.

## Task list (in order; each has acceptance criteria)

**T1. Dataset.** Follow `DATA.md`. Default: Visium HD Human Breast Cancer (Fresh Frozen), Ultima, Space Ranger 4.0.1. Confirm on the dataset page that `segmented_outputs` is downloadable and note file sizes; ask the user before any download over ~10 GB. Write `scripts/download_data.sh` (URLs from the page's batch-download section, plus `sha256sum` of every file into `data/CHECKSUMS.txt`) and `data/DATASET.md` (URL, date accessed, licence, sample description, HER2 status = unknown unless documented).
*Done when:* files present, checksums recorded, metadata written.

**T2. Space Ranger v4 loader** (`load_spaceranger_segmented.py`). Per spec section 1: read `filtered_feature_cell_matrix.h5` and `cell_segmentations.geojson`, map GeoJSON `cell_id` N to barcode `cellid_{N:09d}-1` (verify against the actual files), polygon centroids, convert pixels to µm via `scalefactors_json.json` (verify the key name), write `.h5ad` with `obsm["spatial"]` in µm and `uns["spatial_units"]="um"`, and log input checksums.
*Done when:* unit tests with a tiny hand-made GeoJSON + matrix pass; loader fails loudly on unmatched cells; a centroid scatter overlaid on the hires H&E image lines up (save the figure, and show it to the user).

**T3. Real-data pre-flight (before any phenotype metric).** Confirm all signature genes and ERBB2 are present; plot UMI and genes-per-cell distributions; decide `qc.min_genes`/`min_umi` for whole-transcriptome data and record the decision + reason in `data/DATASET.md` (spec section 0). Only then run classification.

**T4. ClawBio skills.** Clone ClawBio, read its `AGENTS.md`, `skills/catalog.json`, the SKILL.md format and `docs/reproducibility.md`. Check whether the spatial-transcriptomics skill PR (#467, Moran's I + neighbourhood enrichment) is merged and whether it can be reused or extended. Implement the spec as skill(s) with SKILL.md (including the abstention behaviour), demo, and tests; declare chaining partners in the catalog.
*Done when:* skill passes the conformance tests on the synthetic suite and emits schema-valid `classification.json` plus a reproducibility bundle.

**T5. Plots.** Cell-type map, tumour-region overlay with T cells, ERBB2 map over tumour cells, and a stats panel (labels, rule_path, warnings). Plots are skill outputs, logged with checksums.

**T6. Orchestrator.** Calls skills in order, writes `run_log.jsonl`, produces `final_answer.json` (schema in README), renders prose from it, and runs `faithfulness_check.py` as a gate: a failing answer is not shown.

**T7. Literature.** Build `curated_references.json` from the template; a human verifies each entry. Tag subtype and phenotype for each.

**T8. Validation for the demo.** Synthetic suite (all labels correct), 5 to 10 repeated agent runs on real data (`faithfulness_check.py stability`), and one deliberately broken run showing the checker catching it.

## Style

- Python 3.10+, numpy/scipy/pandas; type hints welcome; no notebooks as the source of truth.
- In reports, READMEs and generated prose for this user: avoid em-dashes.
- Keep functions small and testable; every new behaviour gets a test.
