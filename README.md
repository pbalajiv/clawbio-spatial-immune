# Spatial immune phenotype workflow: spec, reference implementation, verification harness

Deterministic classification of one breast cancer tissue section's spatial immune phenotype
(inflamed / excluded / desert) and ERBB2 transcript spatial pattern (homogeneous / heterogeneous)
from segmented spatial transcriptomics data (Visium HD, Space Ranger v4), with a literature layer
and an automated faithfulness gate on top.

## What this is

- A frozen, pre-registered decision rule (`CLASSIFICATION_SPEC.md`), calibrated only on a
  synthetic ground-truth suite before any real dataset was analysed, with every threshold and
  dataset-specific override disclosed and logged per run.
- A deterministic, numpy/scipy pipeline: no ML model is trained, and the LLM never computes,
  rounds, or recomputes a statistic — it only copies numbers out of skill JSON into prose, which
  is then checked, not trusted, by `faithfulness_check.py` before being shown to anyone.
- Every run is reproducible: checksummed inputs/outputs, a logged config (including overrides),
  and a replay command.
- A real contribution back to the open-source [ClawBio](https://github.com/ClawBio/ClawBio) skill
  library (`skills/spatial-immune-phenotype/` in that repo), not just a one-off script.

## What this is NOT

- **Not a clinical tool and not a per-patient outcome prediction.** It is pattern classification
  against known categories on a single tissue section, one timepoint, no replicate.
- **Not a validated HER2 clinical score.** ERBB2 spatial pattern here is a transcript-level
  statistic; clinical HER2 status is determined by IHC/ISH, a different assay entirely.
- **Not literature-verified yet.** `curated_references.json` entries are drafted with real
  citations but every one is `verified: false` until a human actually reads the paper — that flag
  cannot be flipped by an LLM, by design.
- **Not independently validated against a second cell-typing method.** Cell typing here is
  4-marker-signature scoring only; there is no unsupervised-clustering cross-check yet.

## Results on the real dataset

Visium HD Human Breast Cancer (Fresh Frozen), Ultima sequencing (see `data/DATASET.md` for full
provenance, licence, and QC decisions). This sample is DCIS (in situ), and its HER2 receptor
status is undocumented (`sample_subtype: "unknown"`) — both are load-bearing caveats, not footnotes.

| | |
|---|---|
| Immune phenotype | **excluded** (infiltration_ratio 0.049, nhood_z −62.99, p = 0.001) |
| ERBB2 pattern | **homogeneous** (Moran's I excess −0.017) |
| QC | 35,576 / 38,993 cells passed (91.2%) |
| Warnings | `spillover_suspected` (29.65×), `depth_autocorrelation` (0.420) |

The DCIS caveat matters here specifically: an "excluded" T-cell pattern in an in-situ lesion may
partly reflect an intact myoepithelial/basement-membrane barrier rather than only a biological
immune-evasion mechanism (`REF-005`). See `examples/real_data_demo/final_answer.json` for the
full, faithfulness-gated answer (and its figures/report), `examples/t8_validation/` for the
synthetic ground-truth, stability, and deliberately-broken-run validation artifacts.

## Quick start

```bash
pip install -r requirements.txt
pytest -q                                              # 68 tests
python synthetic_tissue.py --suite --out synthetic/    # ground-truth tissues
python classify_reference.py --input synthetic/excluded_heterogeneous --out runs/demo/outputs/classification.json
python orchestrator.py run --input data/tissue.h5ad --config data/real_data_config.json --out runs/demo --seed 0
# ... draft runs/demo/final_answer.json from outputs/classification.json ...
python orchestrator.py gate --run-dir runs/demo --refs curated_references.json
```

Start with `CLAUDE.md` (project context and task list) and `CLASSIFICATION_SPEC.md` (the rule).

## Files

| File | What it does |
|---|---|
| `synthetic_tissue.py` | Generates tissues whose immune phenotype and ERBB2 pattern are known by construction |
| `faithfulness_check.py` | Deterministic checks that the agent's final answer is faithful to skill outputs; also run-to-run stability |
| `make_examples.py` | Writes a faithful run, a run with planted errors, and a curated-reference template |
| `test_harness.py` | 21 tests for the generator and checker |
| `CLASSIFICATION_SPEC.md` | The classification rule (inputs, thresholds, decision order, output contract) |
| `classification_config.json` | Frozen v1.0.0 thresholds |
| `classification.schema.json` | JSON Schema for `outputs/classification.json` |
| `classify_reference.py` | Reference implementation of the spec |
| `test_classify_reference.py` | 11 conformance tests incl. end-to-end with the checker |
| `load_spaceranger_segmented.py` | Space Ranger v4 segmented-output loader (T2) |
| `plots.py` | Cell-type map, tumour-region/T-cell overlay, biomarker map, stats panel (T5) |
| `orchestrator.py` | Runs the skill, writes `run_log.jsonl`, gates `final_answer.json` via `faithfulness_check.py` (T6) |
| `curated_references.json` | Curated literature, tagged by subtype/phenotype, `verified: false` until human-checked (T7) |
| `scripts/validate_synthetic_suite.py` | T8: confirms every synthetic tissue's label matches ground truth |
| `scripts/run_stability_demo.py` | T8: N repeated real-data runs + `faithfulness_check.py stability` |
| `scripts/run_broken_demo.py` | T8: plants known errors in a real run and shows the gate catching them |
| `DATA.md` | Dataset options for the real-data run |
| `data/DATASET.md` | Chosen dataset, provenance, QC decisions, DCIS/HER2 caveats |
| `CLAUDE.md` | Handoff for Claude Code |

Requirements: `pip install -r requirements.txt`, then `pytest -q` (68 tests).

## 1. Ground-truth tissues

```bash
python synthetic_tissue.py --suite --out synthetic/
```

Produces six tissues: inflamed and excluded (each with homogeneous and heterogeneous ERBB2), desert, and a low-quality tissue whose expected answer is `unclassifiable` for both calls. Add `--spillover 0.3` to simulate keratin reads bleeding into immune cells that touch tumour cells, which is how you test whether your marker-concordance check catches segmentation artefacts.

Give the pipeline only `tissue.h5ad` (or `counts.npz` + `genes.txt` + `cells.csv`). `truth_cells.csv` and `ground_truth.json` are for scoring and must never be an input.

`ground_truth.json` includes a self-check computed from the true labels, so you can confirm the scenarios are separable before blaming your pipeline. On the default suite, immune-in-nest fraction is about 0.49 (inflamed) vs 0.02 (excluded), immune fraction is about 0.004 for desert, and within-tumour ERBB2 Moran's I is about 0.57 (heterogeneous) vs about 0 (homogeneous).

Also test borderline tissues (e.g. edit `SCENARIOS` to put the immune fraction near your desert threshold). Clear-cut cases prove the plumbing; borderline cases show how your thresholds behave.

## 2. Faithfulness check

```bash
python make_examples.py
python faithfulness_check.py check --run-dir examples/run_good --refs examples/curated_references.template.json
python faithfulness_check.py check --run-dir examples/run_bad  --refs examples/curated_references.template.json
# on a synthetic tissue, also score against ground truth:
python faithfulness_check.py check --run-dir my_run --refs refs.json --truth synthetic/excluded_heterogeneous/ground_truth.json
```

Exit code is 0 on pass and 1 on any failure, so it can gate CI or the orchestrator itself (refuse to show an answer that fails).

Checks performed:

| Check | Fails when |
|---|---|
| label_faithfulness | final answer's phenotype, pattern or gene differs from `classification.json` |
| abstention | skill returned `unclassifiable` but the answer gives no `abstention_reason` |
| numeric_faithfulness | any number in `narrative`/`conclusion` matches no skill output (rounding and % allowed; `p < x` checked against p-value fields) |
| paired_statistics | a p-value is reported without its effect size, or vice versa |
| citation_grounding | a cited or mentioned `REF-###` is missing from the curated file or from the citations list |
| subtype_match | a cited reference is from another subtype and the mismatch is not acknowledged |
| caveat | caveat lacks the required terms |
| qc_reporting | dropped-cell count differs from the skill, or skill QC counts are internally inconsistent |
| trajectory_order | a required skill was skipped or called out of order |
| provenance | a skill has no version, or a stochastic skill has no seed |
| chaining / output_checksums | a skill consumed a file with a different checksum than was produced, or an output was altered after being written |
| ground_truth_* | (with `--truth`) the skill's label differs from the known answer |

Warnings (do not fail the run): unverified references, a cited reference that does not cover the reported label, other phenotype names in the conclusion, missing agent metadata.

## 3. Stability across repeated agent runs

Run the full agent pipeline N times (5 to 10) on the same input into separate run directories, then:

```bash
python faithfulness_check.py stability --runs run_1 run_2 run_3 run_4 run_5 --refs refs.json
```

Reports label agreement (should be 1.0; labels come from deterministic skills), faithfulness pass rate (should be 1.0), and mean/min citation Jaccard (should be high, need not be 1.0).

## Contracts your skills and agent must follow

**`outputs/classification.json`** (written by the classify skill): full contract in `classification.schema.json` and `CLASSIFICATION_SPEC.md` section 9. The checker relies on `immune_phenotype.{label,nhood_z,nhood_p}`, `biomarker.{label,gene,morans_i,morans_p}` and `qc.{n_cells_input,n_cells_passed,n_dropped}`.
Any other `outputs/*.json` is also indexed for number matching, keyed as `<file stem>.<json path>`.

**`run_log.jsonl`** (one record per line; written by your orchestration wrapper, not by the LLM)
```json
{"type": "agent", "model": "...", "prompt_version": "orchestrator-v1", "temperature": 0}
{"type": "skill", "skill": "spatial_stats", "version": "0.1.0", "seed": 42,
 "inputs": {"outputs/scores.json": "<sha256>"}, "outputs": {"outputs/spatial_stats.json": "<sha256>"}}
```

**`final_answer.json`** (the agent's only output; prose is rendered from it)
```json
{"immune_phenotype": "excluded", "biomarker_pattern": "heterogeneous", "biomarker_gene": "ERBB2",
 "sample_subtype": "HER2+", "subtype_mismatch_acknowledged": false,
 "citations": ["REF-001", "REF-004"], "dropped_cells_reported": 688,
 "caveat": "...", "abstention_reason": null, "conclusion": "...", "narrative": "..."}
```

**Curated references**: list of `{id, verified, subtypes, phenotypes, citation, claim}`. The template entries are starting points and are all `verified: false`; read each paper and confirm the claim before flipping the flag.

Skill names in `expected_order` and `stochastic_skills`, caveat terms, and allowed constants are set in `DEFAULT_CONFIG` and can be overridden with `--config my_config.json`.

## Known limitations

- Number matching works by value. If two statistics share a value (as both p-values do in the example), the checker cannot tell which one a sentence refers to. Keep `p` next to its statistic in prose, and rely on the paired-statistics check.
- Integers in prose ("3 phenotypes") must match a skill value or be listed in `allowed_constants`. Cite papers by `REF-###` in the narrative rather than author-year, or years will be flagged.
- Caveat terms are substring matches; tune `required_caveat_terms` to your actual caveat wording.
- These checks verify faithfulness and provenance, not scientific correctness. Correctness comes from the ground-truth tissues and the skills' own tests.
