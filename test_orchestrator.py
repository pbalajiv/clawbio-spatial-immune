"""Tests for orchestrator.py (T6): run() output contract, and the faithfulness gate."""
import json
from pathlib import Path

import pytest

import orchestrator as orch
import synthetic_tissue as st

REFS = []  # no citations exercised by these tests; T7 supplies real content


@pytest.fixture(scope="module")
def tissue(tmp_path_factory):
    out = tmp_path_factory.mktemp("orch_tissue")
    st.generate("excluded", "heterogeneous", out / "tissue", n_cells=3000, seed=0)
    return out / "tissue"


@pytest.fixture(scope="module")
def run_dir(tissue, tmp_path_factory):
    out = tmp_path_factory.mktemp("orch_run")
    result = orch.run(tissue, out / "run1", config_path=None, seed=0)
    return out / "run1", result


def _minimal_valid_final_answer(result):
    imm, bio, qc = result["immune_phenotype"], result["biomarker"], result["qc"]
    return {
        "immune_phenotype": imm["label"], "biomarker_pattern": bio["label"], "biomarker_gene": bio["gene"],
        "sample_subtype": "unknown", "subtype_mismatch_acknowledged": False,
        "citations": [], "dropped_cells_reported": qc["n_dropped"],
        "caveat": "This is pattern classification, not a validated prediction of anything clinical.",
        "abstention_reason": imm["reason"] if imm["label"] == "unclassifiable" else None,
        "conclusion": f"{bio['gene']} pattern reported.",
        "narrative": f"{qc['n_dropped']} cells were dropped at QC out of {qc['n_cells_input']}.",
    }


def test_run_writes_expected_output_tree(run_dir):
    out, result = run_dir
    assert (out / "outputs" / "classification.json").exists()
    assert (out / "outputs" / "qc.json").exists()
    assert (out / "outputs" / "scores.json").exists()
    assert (out / "outputs" / "spatial_stats.json").exists()
    assert (out / "outputs" / "report.md").exists()
    assert (out / "run_log.jsonl").exists()
    for name in ("cell_type_map", "tumor_region_t_cells", "biomarker_map", "stats_panel"):
        assert (out / "figures" / f"{name}.png").exists()


def test_run_log_has_agent_record_and_expected_skill_order(run_dir):
    out, _ = run_dir
    records = [json.loads(l) for l in (out / "run_log.jsonl").read_text().splitlines() if l.strip()]
    assert records[0]["type"] == "agent"
    assert records[0]["model"] and records[0]["prompt_version"] is not None
    skill_names = [r["skill"] for r in records if r["type"] == "skill"]
    assert skill_names == ["load_qc", "score_cells", "spatial_stats", "classify", "report"]
    spatial_stats_record = next(r for r in records if r.get("skill") == "spatial_stats")
    assert spatial_stats_record["seed"] == 0  # stochastic skill must log its seed


def test_gate_passes_on_a_faithful_final_answer(run_dir, tmp_path):
    out, result = run_dir
    (out / "final_answer.json").write_text(json.dumps(_minimal_valid_final_answer(result)))
    refs_path = tmp_path / "refs.json"
    refs_path.write_text(json.dumps(REFS))
    report = orch.gate(out, refs_path)
    assert report["status"] == "pass", report["results"]


def test_gate_fails_on_a_wrong_label(run_dir, tmp_path):
    out, result = run_dir
    final = _minimal_valid_final_answer(result)
    final["immune_phenotype"] = "inflamed" if result["immune_phenotype"]["label"] != "inflamed" else "excluded"
    (out / "final_answer_bad.json").write_text(json.dumps(final))
    # swap in the bad answer, run the gate, then restore the good one for other tests in this module
    good = (out / "final_answer.json").read_text()
    (out / "final_answer.json").write_text(json.dumps(final))
    try:
        refs_path = tmp_path / "refs.json"
        refs_path.write_text(json.dumps(REFS))
        report = orch.gate(out, refs_path)
        assert report["status"] == "fail"
        assert any(r["check"] == "label_faithfulness" and r["status"] == "fail" for r in report["results"])
    finally:
        (out / "final_answer.json").write_text(good)


def test_gate_fails_on_wrong_dropped_cell_count(run_dir, tmp_path):
    out, result = run_dir
    final = _minimal_valid_final_answer(result)
    final["dropped_cells_reported"] = result["qc"]["n_dropped"] + 1
    good = (out / "final_answer.json").read_text()
    (out / "final_answer.json").write_text(json.dumps(final))
    try:
        refs_path = tmp_path / "refs.json"
        refs_path.write_text(json.dumps(REFS))
        report = orch.gate(out, refs_path)
        assert report["status"] == "fail"
        assert any(r["check"] == "qc_reporting" and r["status"] == "fail" for r in report["results"])
    finally:
        (out / "final_answer.json").write_text(good)
