"""Tests for the verification harness itself. Run: pytest -q test_harness.py"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import faithfulness_check as fc
import make_examples
import synthetic_tissue as st

HERE = Path(__file__).parent
REFS = make_examples.REFS


# ---------------------------------------------------------------- synthetic tissues

@pytest.fixture(scope="module")
def suite(tmp_path_factory):
    out = tmp_path_factory.mktemp("suite")
    return {f"{s}_{h}" + ("_lowqc" if x.get("low_quality") else ""):
            st.generate(s, h, out / f"{s}_{h}", n_cells=8000, seed=1, **x)
            for s, h, x in st.SUITE}


def test_labels_fixed_by_construction(suite):
    for name, t in suite.items():
        if name.endswith("_lowqc"):
            assert t["expected_immune_phenotype"] == "unclassifiable"
        else:
            assert t["expected_immune_phenotype"] == name.split("_")[0]
            assert t["expected_biomarker_pattern"] == name.split("_")[1]


def test_immune_scenarios_separable(suite):
    sc = {k: v["self_check_from_true_labels"] for k, v in suite.items()}
    assert sc["desert_homogeneous"]["immune_fraction"] < 0.01
    assert sc["excluded_homogeneous"]["immune_fraction"] > 0.05
    assert sc["inflamed_homogeneous"]["immune_in_nest_fraction"] > 5 * sc["excluded_homogeneous"]["immune_in_nest_fraction"]


def test_erbb2_heterogeneity_separable(suite):
    sc = {k: v["self_check_from_true_labels"] for k, v in suite.items()}
    for phen in ("inflamed", "excluded"):
        assert sc[f"{phen}_heterogeneous"]["erbb2_morans_i_within_tumor"] > 0.3
        assert abs(sc[f"{phen}_homogeneous"]["erbb2_morans_i_within_tumor"]) < 0.1


def test_generator_is_deterministic(tmp_path):
    a = st.generate("excluded", "heterogeneous", tmp_path / "a", n_cells=3000, seed=7)
    b = st.generate("excluded", "heterogeneous", tmp_path / "b", n_cells=3000, seed=7)
    assert a["self_check_from_true_labels"] == b["self_check_from_true_labels"]
    assert fc.sha256(tmp_path / "a" / "counts.npz") == fc.sha256(tmp_path / "b" / "counts.npz")


def test_truth_not_in_pipeline_inputs(tmp_path):
    st.generate("inflamed", "homogeneous", tmp_path, n_cells=2000, seed=0)
    header = (tmp_path / "cells.csv").read_text().splitlines()[0]
    assert header == "cell_id,x,y"


# ---------------------------------------------------------------- faithfulness checker

@pytest.fixture()
def runs(tmp_path, monkeypatch):
    monkeypatch.setattr(make_examples, "HERE", tmp_path)
    make_examples.write_run("run_good", make_examples.GOOD)
    make_examples.write_run("run_bad", make_examples.BAD, drop_seed=True)
    return tmp_path


def _failed_checks(rep):
    return {r["check"] for r in rep["results"] if r["status"] == "fail"}


def test_good_run_passes(runs):
    rep = fc.run_checks(runs / "run_good", REFS, fc.DEFAULT_CONFIG)
    assert rep["status"] == "pass", rep


def test_bad_run_catches_each_planted_error(runs):
    rep = fc.run_checks(runs / "run_bad", REFS, fc.DEFAULT_CONFIG)
    assert _failed_checks(rep) >= {"label_faithfulness", "numeric_faithfulness", "paired_statistics",
                                   "subtype_match", "citation_grounding", "caveat", "qc_reporting",
                                   "provenance"}


def test_tampered_output_detected(runs):
    p = runs / "run_good" / "outputs" / "classification.json"
    obj = json.loads(p.read_text())
    obj["immune_phenotype"]["nhood_z"] = -9.99
    p.write_text(json.dumps(obj))
    rep = fc.run_checks(runs / "run_good", REFS, fc.DEFAULT_CONFIG)
    assert "output_checksums" in _failed_checks(rep)


def test_abstention_requires_reason(runs):
    d = runs / "run_good"
    cls = json.loads((d / "outputs" / "classification.json").read_text())
    cls["immune_phenotype"]["label"] = "unclassifiable"
    (d / "outputs" / "classification.json").write_text(json.dumps(cls))
    final = json.loads((d / "final_answer.json").read_text())
    final["immune_phenotype"] = "unclassifiable"
    (d / "final_answer.json").write_text(json.dumps(final))
    rep = fc.run_checks(d, REFS, fc.DEFAULT_CONFIG)
    assert "abstention" in _failed_checks(rep)


def test_ground_truth_check(runs):
    truth = {"expected_immune_phenotype": "inflamed", "expected_biomarker_pattern": "heterogeneous"}
    rep = fc.run_checks(runs / "run_good", REFS, fc.DEFAULT_CONFIG, truth)
    assert "ground_truth_immune" in _failed_checks(rep)
    assert "ground_truth_pattern" not in _failed_checks(rep)


@pytest.mark.parametrize("token,value,ok", [
    ("-4.21", -4.2137, True), ("-4.2", -4.21, True), ("-3.9", -4.21, False),
    ("8.3%", 0.083, True), ("1.0e-03", 0.001, True), ("688", 688, True), ("0.41", 0.412, True),
])
def test_number_matching(token, value, ok):
    assert fc.number_matches(token, value) is ok


def test_gene_names_not_parsed_as_numbers():
    toks = [t for t, _ in fc.extract_numbers("ERBB2, CD8A, MKI67, PD-1, Space Ranger v4, REF-004", r"\bREF-\d{3}\b")]
    assert toks == []


def test_threshold_pvalue(runs):
    d = runs / "run_good"
    final = json.loads((d / "final_answer.json").read_text())
    final["narrative"] = final["narrative"].replace("(p = 0.001, 1000 permutations)", "(p < 0.01)")
    (d / "final_answer.json").write_text(json.dumps(final))
    assert fc.run_checks(d, REFS, fc.DEFAULT_CONFIG)["status"] == "pass"


def test_stability(runs):
    copies = []
    for i in range(3):
        dst = runs / f"rep{i}"
        shutil.copytree(runs / "run_good", dst)
        copies.append(dst)
    rep = fc.stability(copies, REFS, fc.DEFAULT_CONFIG)
    assert rep["immune_phenotype"]["agreement"] == 1.0
    assert rep["faithfulness_pass_rate"] == 1.0
    assert rep["citation_jaccard_mean"] == 1.0


def test_cli_exit_codes(runs):
    refs = runs / "refs.json"
    refs.write_text(json.dumps(REFS))
    good = subprocess.run([sys.executable, str(HERE / "faithfulness_check.py"), "check",
                           "--run-dir", str(runs / "run_good"), "--refs", str(refs)], capture_output=True)
    bad = subprocess.run([sys.executable, str(HERE / "faithfulness_check.py"), "check",
                          "--run-dir", str(runs / "run_bad"), "--refs", str(refs)], capture_output=True)
    assert good.returncode == 0 and bad.returncode == 1
