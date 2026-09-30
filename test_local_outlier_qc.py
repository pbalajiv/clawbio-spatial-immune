"""Tests for the spatially-local outlier QC (SpotSweeper-style), classify_reference.local_outlier_qc / qc()."""
import numpy as np
import pytest
from scipy import sparse

import classify_reference as cr

BASE_CFG = {
    "qc": {"min_umi": 20, "min_genes": 2, "min_pass_fraction": 0.5},
    "graph": {"max_edge_um": 100.0},
    "spatial_qc": {"k": 8, "z_threshold": 3.0},
}


def _line_coords(n):
    return np.column_stack([np.arange(n, dtype=float) * 10.0, np.zeros(n)])


def test_spatial_qc_absent_from_config_leaves_qc_unchanged():
    """No 'spatial_qc' key -> qc() behaves exactly as before (backward compatible with frozen v1.0.0 config)."""
    n_genes = 3
    X = sparse.csr_matrix(np.full((10, n_genes), 100, dtype=np.int32))
    xy = _line_coords(10)
    cfg = {"qc": BASE_CFG["qc"]}  # no spatial_qc key
    keep, qcd = cr.qc(X, ["EPCAM", "PTPRC", "MT-ND1"], xy, cfg)
    assert set(qcd["drop_reasons"]) == {"low_umi", "low_genes"}
    assert keep.all()


def test_local_outlier_flags_cell_with_anomalously_low_libsize():
    genes = ["EPCAM", "PTPRC", "COL1A1"]
    counts = np.full((9, 3), 1000 // 3, dtype=np.int32)  # ~333 per gene -> ~1000 total, 3 genes detected
    counts = np.vstack([counts, [[5, 5, 5]]])  # cell 9: total=15... too close to global min_umi; use bigger baseline
    # Use a bigger baseline so the outlier cell is well above min_umi=20 but still a local outlier.
    counts = np.full((9, 3), 2000 // 3, dtype=np.int32)
    counts = np.vstack([counts, [[30, 30, 30]]])  # cell 9: total=90 >> min_umi=20, but far below neighbours' ~2000
    X = sparse.csr_matrix(counts)
    xy = _line_coords(10)
    total = np.asarray(X.sum(1)).ravel()
    ngenes = np.asarray((X > 0).sum(1)).ravel()

    low_libsize, low_genes, high_mito = cr.local_outlier_qc(xy, total, ngenes, None, BASE_CFG)
    assert low_libsize[9] and not low_libsize[:9].any()
    assert not low_genes.any()  # every cell detects all 3 genes; no gene-count anomaly
    assert not high_mito.any()  # no mito_frac provided


def test_local_outlier_flags_cell_with_anomalously_high_mito():
    genes = ["EPCAM", "PTPRC", "MT-ND1"]
    counts = np.tile([300, 300, 20], (10, 1)).astype(np.int32)  # ~3% mito baseline
    counts[9] = [150, 150, 500]  # cell 9: same total-ish, but mostly mitochondrial
    X = sparse.csr_matrix(counts)
    xy = _line_coords(10)
    total = np.asarray(X.sum(1)).ravel()
    ngenes = np.asarray((X > 0).sum(1)).ravel()
    mito_frac = cr._mitochondrial_fraction(X, genes, total)

    low_libsize, low_genes, high_mito = cr.local_outlier_qc(xy, total, ngenes, mito_frac, BASE_CFG)
    assert high_mito[9] and not high_mito[:9].any()


def test_mitochondrial_fraction_none_when_no_mt_genes():
    genes = ["EPCAM", "PTPRC"]
    X = sparse.csr_matrix(np.full((5, 2), 100, dtype=np.int32))
    total = np.asarray(X.sum(1)).ravel()
    assert cr._mitochondrial_fraction(X, genes, total) is None


def test_qc_drop_reasons_are_mutually_exclusive_and_sum_to_n_dropped():
    """A cell already dropped by the global min_umi floor must not also be counted as a spatial outlier."""
    genes = ["EPCAM", "PTPRC", "MT-ND1"]
    counts = np.tile([300, 300, 10], (9, 1)).astype(np.int32)
    counts = np.vstack([counts, [[2, 2, 1]]])  # cell 9: total=5 < min_umi=20 -> caught globally, not locally
    X = sparse.csr_matrix(counts)
    xy = _line_coords(10)
    cfg = {**BASE_CFG}

    keep, qcd = cr.qc(X, genes, xy, cfg)
    assert qcd["drop_reasons"]["low_umi"] == 1
    assert qcd["drop_reasons"]["spatial_outlier_low_libsize"] == 0  # not double-counted
    assert sum(qcd["drop_reasons"].values()) == qcd["n_dropped"]
    assert not keep[9]


def test_qc_output_is_schema_shaped_with_spatial_qc_enabled():
    import json
    from pathlib import Path
    genes = ["EPCAM", "PTPRC", "MT-ND1"]
    X = sparse.csr_matrix(np.full((10, 3), 300, dtype=np.int32))
    xy = _line_coords(10)
    keep, qcd = cr.qc(X, genes, xy, BASE_CFG)
    assert set(qcd) == {"n_cells_input", "n_cells_passed", "n_dropped", "drop_reasons", "pass_fraction"}
    assert all(isinstance(v, int) for v in qcd["drop_reasons"].values())
