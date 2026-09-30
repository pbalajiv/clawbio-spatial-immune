"""Tests that classify_reference.py implements CLASSIFICATION_SPEC.md v1.0.0.

Run: pytest -q test_classify_reference.py
Any ClawBio skill implementing the spec must pass the same ground-truth tests.
"""
import copy
import json
from pathlib import Path

import jsonschema
import numpy as np
import pytest
from scipy import sparse

import classify_reference as cr
import faithfulness_check as fc
import make_examples
import synthetic_tissue as st

HERE = Path(__file__).parent
SCHEMA = json.loads((HERE / "classification.schema.json").read_text())
N_CELLS = 8000


@pytest.fixture(scope="module")
def cfg():
    c = cr.load_config()
    c["permutations"]["n"] = 200          # speed; min p = 1/201 < alpha 0.01
    return c


@pytest.fixture(scope="module")
def suite(tmp_path_factory):
    out = tmp_path_factory.mktemp("suite")
    tissues = {}
    for s, h, x in st.SUITE:
        name = f"{s}_{h}" + ("_lowqc" if x.get("low_quality") else "")
        tissues[name] = (out / name, st.generate(s, h, out / name, n_cells=N_CELLS, seed=3, **x))
    spill = out / "spill"
    tissues["spill"] = (spill, st.generate("inflamed", "homogeneous", spill, n_cells=N_CELLS, seed=3, spillover=0.3))
    return tissues


@pytest.fixture(scope="module")
def results(suite, cfg):
    res = {}
    for name, (d, _) in suite.items():
        X, genes, xy, prov = cr.load_input(d)
        res[name] = cr.classify(X, genes, xy, cfg, seed=0, provenance=prov)
    return res


def test_ground_truth_labels(suite, results):
    for name, (_, truth) in suite.items():
        r = results[name]
        assert r["immune_phenotype"]["label"] == truth["expected_immune_phenotype"], (name, r["immune_phenotype"])
        assert r["biomarker"]["label"] == truth["expected_biomarker_pattern"], (name, r["biomarker"])


def test_outputs_match_schema(results):
    for name, r in results.items():
        jsonschema.validate(json.loads(json.dumps(r)), SCHEMA)


def test_lowqc_abstains_with_reason(results):
    r = results["inflamed_homogeneous_lowqc"]
    assert r["immune_phenotype"]["reason"] == "qc_fail" and r["biomarker"]["reason"] == "qc_fail"


def test_zero_cells_pass_qc_abstains_without_crashing(suite, cfg):
    """A pathological config (or input) that drops every cell must abstain, not raise.

    Regression test: score_and_type() used to run unconditionally before the
    pass_fraction check, crashing on an empty post-QC array instead of hitting
    the qc_fail branch.
    """
    d, _ = suite["inflamed_homogeneous"]
    X, genes, xy, prov = cr.load_input(d)
    c = copy.deepcopy(cfg)
    c["qc"]["min_umi"] = 10**9
    r = cr.classify(X, genes, xy, c, seed=0, provenance=prov)
    assert r["qc"]["n_cells_passed"] == 0
    assert r["immune_phenotype"]["label"] == "unclassifiable" and r["immune_phenotype"]["reason"] == "qc_fail"
    assert r["biomarker"]["label"] == "unclassifiable" and r["biomarker"]["reason"] == "qc_fail"
    jsonschema.validate(json.loads(json.dumps(r)), SCHEMA)


def test_spillover_flagged_only_when_present(results):
    assert any(w.startswith("spillover_suspected") for w in results["spill"]["warnings"])
    for name, r in results.items():
        if name != "spill":
            assert not any(w.startswith("spillover_suspected") for w in r["warnings"]), name


def test_spillover_does_not_flip_label(results):
    assert results["spill"]["immune_phenotype"]["label"] == "inflamed"


def test_every_decision_has_rule_path(results):
    for r in results.values():
        assert r["immune_phenotype"]["rule_path"] and r["biomarker"]["rule_path"]


def test_deterministic(suite, cfg):
    d, _ = suite["excluded_heterogeneous"]
    X, genes, xy, prov = cr.load_input(d)
    a = cr.classify(X, genes, xy, cfg, seed=5, provenance=prov)
    b = cr.classify(X, genes, xy, cfg, seed=5, provenance=prov)
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def test_biomarker_in_tumor_signature_rejected(tmp_path):
    bad = {"signatures": {"tumor": ["EPCAM", "KRT8", "ERBB2"]}}
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(bad))
    with pytest.raises(ValueError, match="circularity"):
        cr.load_config(p)


def test_biomarker_not_expressed(suite, cfg):
    d, _ = suite["inflamed_homogeneous"]
    X, genes, xy, prov = cr.load_input(d)
    X = X.tolil()
    X[:, genes.index("ERBB2")] = 0
    r = cr.classify(sparse.csr_matrix(X), genes, xy, cfg, seed=0, provenance=prov)
    assert r["biomarker"]["label"] == "unclassifiable"
    assert r["biomarker"]["reason"] == "biomarker_not_expressed"


def test_borderline_zone_abstains(suite, cfg):
    """Force the observed ratio into the borderline zone by moving thresholds around it."""
    d, _ = suite["inflamed_homogeneous"]
    X, genes, xy, prov = cr.load_input(d)
    c = copy.deepcopy(cfg)
    c["immune_rule"]["inflamed_min_ratio"] = 5.0
    c["immune_rule"]["excluded_max_ratio"] = 0.01
    r = cr.classify(X, genes, xy, c, seed=0, provenance=prov)
    assert r["immune_phenotype"]["label"] == "unclassifiable"
    assert r["immune_phenotype"]["reason"] == "borderline"


def test_end_to_end_with_faithfulness_checker(suite, results, tmp_path):
    """Reference output -> faithful final answer -> checker passes, including ground truth."""
    name = "excluded_heterogeneous"
    r, truth = results[name], suite[name][1]
    run = tmp_path / "run"
    (run / "outputs").mkdir(parents=True)
    (run / "outputs" / "classification.json").write_text(json.dumps(r))
    sha = fc.sha256(run / "outputs" / "classification.json")
    log = [{"type": "agent", "model": "test", "prompt_version": "v1", "temperature": 0}]
    for skill in ["load_qc", "score_cells", "spatial_stats", "classify", "report"]:
        rec = {"type": "skill", "skill": skill, "version": "1.0.0"}
        if skill == "spatial_stats":
            rec["seed"] = 0
        if skill == "classify":
            rec["outputs"] = {"outputs/classification.json": sha}
        log.append(rec)
    (run / "run_log.jsonl").write_text("\n".join(json.dumps(x) for x in log))
    imm, bio, qc = r["immune_phenotype"], r["biomarker"], r["qc"]
    final = {
        "immune_phenotype": imm["label"], "biomarker_pattern": bio["label"], "biomarker_gene": bio["gene"],
        "sample_subtype": "HER2+", "subtype_mismatch_acknowledged": False,
        "citations": ["REF-001", "REF-004"], "dropped_cells_reported": qc["n_dropped"],
        "caveat": make_examples.CAVEAT, "abstention_reason": None,
        "conclusion": "Immune-excluded tumour with spatially heterogeneous ERBB2.",
        "narrative": (f"Infiltration ratio {imm['infiltration_ratio']:.2f} with tumour-T neighbourhood "
                      f"enrichment z = {imm['nhood_z']:.1f} (p = {imm['nhood_p']:.3f}). ERBB2 Moran's I = "
                      f"{bio['morans_i']:.3f} (p = {bio['morans_p']:.3f}) across {bio['n_tumor_cells']} "
                      f"tumour cells. {qc['n_dropped']} cells were removed at QC (REF-001, REF-004)."),
    }
    (run / "final_answer.json").write_text(json.dumps(final))
    rep = fc.run_checks(run, make_examples.REFS, fc.DEFAULT_CONFIG, truth)
    assert rep["status"] == "pass", [x for x in rep["results"] if x["status"] == "fail"]
