"""Tests for scripts/qc_preflight.py's diagnostic functions (reporting-only QC extensions)."""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from qc_preflight import check_genes, doublet_candidates, mitochondrial_fraction, umi_and_gene_counts


def _make_adata(gene_names, counts):
    X = sparse.csr_matrix(np.asarray(counts, dtype=np.float64))
    a = ad.AnnData(X=X)
    a.var_names = gene_names
    return a


def test_check_genes_reports_detection_rate_not_just_presence():
    # EPCAM detected in 2/4 cells, PTPRC in 0/4 (present but never detected)
    genes = ["EPCAM", "PTPRC", "ERBB2"]
    counts = [
        [5, 0, 3],
        [0, 0, 0],
        [7, 0, 0],
        [0, 0, 1],
    ]
    a = _make_adata(genes, counts)
    cfg = {"signatures": {"tumor": ["EPCAM"], "t_cell": ["PTPRC"]}, "biomarker": {"gene": "ERBB2"}}
    report = check_genes(a, cfg)
    assert report["tumor"]["detection_rate"]["EPCAM"] == pytest.approx(0.5)
    assert report["t_cell"]["detection_rate"]["PTPRC"] == pytest.approx(0.0)
    assert report["biomarker"]["detection_rate"]["ERBB2"] == pytest.approx(0.5)


def test_check_genes_missing_gene_has_none_detection_rate():
    a = _make_adata(["EPCAM"], [[1], [2]])
    cfg = {"signatures": {"t_cell": ["CD8B"]}, "biomarker": {"gene": "EPCAM"}}
    report = check_genes(a, cfg)
    assert report["t_cell"]["missing"] == ["CD8B"]
    assert report["t_cell"]["detection_rate"]["CD8B"] is None


def test_mitochondrial_fraction_none_when_no_mt_genes():
    a = _make_adata(["EPCAM", "PTPRC"], [[1, 2], [3, 4]])
    total_umi, _ = umi_and_gene_counts(a)
    mt_frac, mt_genes = mitochondrial_fraction(a, total_umi)
    assert mt_frac is None
    assert mt_genes == []


def test_mitochondrial_fraction_computed_correctly():
    # cell 0: 8 non-MT + 2 MT = 20% mito; cell 1: 0 non-MT + 5 MT = 100% mito
    genes = ["EPCAM", "MT-ND1", "MT-CO1"]
    counts = [[8, 1, 1], [0, 2, 3]]
    a = _make_adata(genes, counts)
    total_umi, _ = umi_and_gene_counts(a)
    mt_frac, mt_genes = mitochondrial_fraction(a, total_umi)
    assert set(mt_genes) == {"MT-ND1", "MT-CO1"}
    np.testing.assert_allclose(mt_frac, [0.2, 1.0])


def test_doublet_candidates_flags_only_cells_high_in_both_dimensions():
    total_umi = np.array([100, 100, 1000, 1000])
    n_genes = np.array([50, 500, 50, 500])  # median umi=550, median genes=275
    flag, median_umi, median_genes = doublet_candidates(total_umi, n_genes, umi_multiplier=1.5, genes_multiplier=1.5)
    # only the last cell (1000 umi, 500 genes) exceeds 1.5x median in BOTH dimensions
    np.testing.assert_array_equal(flag, [False, False, False, True])
    assert median_umi == 550
    assert median_genes == 275


def test_doublet_candidates_high_in_only_one_dimension_not_flagged():
    """A naturally large, highly-transcribing cell (high UMI, modest gene count) should not be flagged alone."""
    total_umi = np.array([100, 100, 100, 10000])
    n_genes = np.array([50, 50, 50, 55])  # huge UMI outlier but unremarkable gene count
    flag, _, _ = doublet_candidates(total_umi, n_genes, umi_multiplier=2.0, genes_multiplier=2.0)
    assert not flag[-1]
