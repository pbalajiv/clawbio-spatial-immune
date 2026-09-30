"""Tests for plots.py (T5): figures must be consistent with classify_reference.py's own labels."""
import json
from pathlib import Path

import numpy as np
import pytest

import classify_reference as cr
import plots as pl
import synthetic_tissue as st


@pytest.fixture(scope="module")
def tissue(tmp_path_factory):
    out = tmp_path_factory.mktemp("plots_tissue")
    st.generate("excluded", "heterogeneous", out / "tissue", n_cells=3000, seed=0)
    return out / "tissue"


@pytest.fixture(scope="module")
def cfg():
    c = cr.load_config()
    c["permutations"]["n"] = 200
    return c


@pytest.fixture(scope="module")
def result(tissue, cfg):
    X, genes, xy, prov = cr.load_input(tissue)
    return cr.classify(X, genes, xy, cfg, seed=0, provenance=prov)


def test_compute_plot_data_matches_classify_cell_counts(tissue, cfg, result):
    data = pl.compute_plot_data(tissue, cfg)
    counts = result["cell_types"]["counts"]
    for cell_type in ("tumor", "t_cell", "stroma", "ambiguous"):
        assert int((data["types"] == cell_type).sum()) == counts[cell_type]
    assert int(data["region"].sum()) == result["immune_phenotype"]["n_tumor_region_cells"]
    assert len(data["xy"]) == result["qc"]["n_cells_passed"]


def test_render_plots_writes_all_expected_figures(tissue, cfg, result, tmp_path):
    paths, checksums = pl.render_plots(tissue, result, cfg, tmp_path / "figures")
    assert set(paths) == {"cell_type_map", "tumor_region_t_cells", "biomarker_map", "stats_panel"}
    for p in paths.values():
        assert Path(p).exists() and Path(p).stat().st_size > 0
    assert set(checksums) == {f"{name}.png" for name in paths}
    for name, digest in checksums.items():
        assert len(digest) == 64


def test_render_plots_skips_biomarker_map_when_gene_absent(tissue, cfg, result, tmp_path):
    c = json.loads(json.dumps(cfg))
    c["biomarker"]["gene"] = "NOT_A_REAL_GENE"
    paths, checksums = pl.render_plots(tissue, result, c, tmp_path / "figures_no_gene")
    assert "biomarker_map" not in paths
    assert "biomarker_map.png" not in checksums


def test_render_plots_degrades_to_stats_panel_only_on_qc_fail(tissue, cfg, tmp_path):
    """Regression test: compute_plot_data used to crash on a zero-cell qc_fail input."""
    X, genes, xy, prov = cr.load_input(tissue)
    c = json.loads(json.dumps(cfg))
    c["qc"]["min_umi"] = 10**9
    r = cr.classify(X, genes, xy, c, seed=0, provenance=prov)
    assert r["qc"]["n_cells_passed"] == 0

    paths, checksums = pl.render_plots(tissue, r, c, tmp_path / "figures_qc_fail")
    assert set(paths) == {"stats_panel"}
    assert Path(paths["stats_panel"]).exists()


def test_tumor_region_mask_matches_classify_internal_definition(tissue, cfg):
    """plots.py must reuse classify_reference.tumor_region_mask, not a re-derived definition."""
    X, genes, xy, prov = cr.load_input(tissue)
    keep, qcd = cr.qc(X, genes, xy, cfg)
    X, xy = X[keep], xy[keep]
    L, total = cr.lognorm(X)
    types, scores, double_pos, conc = cr.score_and_type(L, genes, cfg)
    expected = cr.tumor_region_mask(xy, types, cfg)

    data = pl.compute_plot_data(tissue, cfg)
    np.testing.assert_array_equal(data["region"], expected)
