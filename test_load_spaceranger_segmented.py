"""Unit tests for load_spaceranger_segmented.py, against a tiny hand-made fixture."""
import json

import h5py
import numpy as np
import pytest
import scipy.sparse as sparse

from load_spaceranger_segmented import load_spaceranger_segmented

GENES = ["EPCAM", "PTPRC"]
MICRONS_PER_PIXEL = 0.5


def _write_matrix_h5(path, barcodes, genes, counts):
    """counts: dense array [n_cells x n_genes]."""
    X = sparse.csr_matrix(np.asarray(counts, dtype=np.int32))
    X_gc = X.T.tocsc()  # genes x cells, as Space Ranger writes it
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("data", data=X_gc.data)
        g.create_dataset("indices", data=X_gc.indices)
        g.create_dataset("indptr", data=X_gc.indptr)
        g.create_dataset("shape", data=np.array(X_gc.shape, dtype=np.int32))
        g.create_dataset("barcodes", data=np.array(barcodes, dtype="S18"))
        feat = g.create_group("features")
        feat.create_dataset("name", data=np.array(genes, dtype="S15"))


def _square_ring(cx, cy, half=2.0):
    return [
        [cx - half, cy - half],
        [cx - half, cy + half],
        [cx + half, cy + half],
        [cx + half, cy - half],
        [cx - half, cy - half],
    ]


def _write_geojson(path, cell_ids_and_centers):
    features = [
        {
            "type": "Feature",
            "geometry": {"type": "Polygon", "coordinates": [_square_ring(cx, cy)]},
            "properties": {"cell_id": cell_id},
        }
        for cell_id, (cx, cy) in cell_ids_and_centers
    ]
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))


def _write_scalefactors(path, microns_per_pixel=MICRONS_PER_PIXEL, **extra):
    d = {"microns_per_pixel": microns_per_pixel, **extra}
    path.write_text(json.dumps(d))


def _make_fixture(tmp_path, n_cells=3, drop_polygon_for=None, extra_polygon=False,
                   scalefactors_key="microns_per_pixel"):
    seg_dir = tmp_path / "segmented_outputs"
    (seg_dir / "spatial").mkdir(parents=True)

    barcodes = [f"cellid_{i:09d}-1" for i in range(1, n_cells + 1)]
    counts = np.arange(1, n_cells * len(GENES) + 1).reshape(n_cells, len(GENES))
    _write_matrix_h5(seg_dir / "filtered_feature_cell_matrix.h5", barcodes, GENES, counts)

    centers = [(cell_id, (10.0 * cell_id, 20.0 * cell_id)) for cell_id in range(1, n_cells + 1)]
    if drop_polygon_for is not None:
        centers = [c for c in centers if c[0] != drop_polygon_for]
    if extra_polygon:
        centers.append((n_cells + 1, (999.0, 999.0)))
    _write_geojson(seg_dir / "cell_segmentations.geojson", centers)

    scalefactors_path = seg_dir / "spatial" / "scalefactors_json.json"
    if scalefactors_key == "microns_per_pixel":
        _write_scalefactors(scalefactors_path)
    else:
        _write_scalefactors(scalefactors_path, **{scalefactors_key: MICRONS_PER_PIXEL})
        scalefactors_path.write_text(json.dumps({scalefactors_key: MICRONS_PER_PIXEL}))

    return seg_dir


def test_happy_path_centroids_and_units(tmp_path):
    seg_dir = _make_fixture(tmp_path, n_cells=3)
    adata = load_spaceranger_segmented(seg_dir, seed=7)

    assert adata.n_obs == 3
    assert adata.n_vars == len(GENES)
    assert list(adata.var_names) == GENES
    assert list(adata.obs_names) == [f"cellid_{i:09d}-1" for i in (1, 2, 3)]
    assert adata.uns["spatial_units"] == "um"
    assert adata.uns["seed"] == 7

    expected_px = np.array([[10.0, 20.0], [20.0, 40.0], [30.0, 60.0]])
    np.testing.assert_allclose(adata.obsm["spatial"], expected_px * MICRONS_PER_PIXEL)


def test_cell_area_computed_in_square_micrometres(tmp_path):
    seg_dir = _make_fixture(tmp_path, n_cells=2)  # squares are 4x4 px (half=2.0) -> 16 px^2 each
    adata = load_spaceranger_segmented(seg_dir)
    expected_area_um2 = 16.0 * MICRONS_PER_PIXEL ** 2
    np.testing.assert_allclose(adata.obs["cell_area_um2"].to_numpy(), [expected_area_um2] * 2)


def test_counts_transposed_correctly(tmp_path):
    seg_dir = _make_fixture(tmp_path, n_cells=3)
    adata = load_spaceranger_segmented(seg_dir)
    expected = np.arange(1, 3 * len(GENES) + 1).reshape(3, len(GENES))
    np.testing.assert_array_equal(np.asarray(adata.X.todense()), expected)


def test_checksums_logged(tmp_path):
    seg_dir = _make_fixture(tmp_path, n_cells=2)
    adata = load_spaceranger_segmented(seg_dir)
    checksums = adata.uns["source"]["input_checksums"]
    assert set(checksums) == {
        "filtered_feature_cell_matrix.h5",
        "cell_segmentations.geojson",
        "scalefactors_json.json",
    }
    assert all(len(v) == 64 for v in checksums.values())  # sha256 hex digest length


def test_missing_polygon_for_matrix_cell_fails_loudly(tmp_path):
    seg_dir = _make_fixture(tmp_path, n_cells=3, drop_polygon_for=2)
    with pytest.raises(ValueError, match="no polygon"):
        load_spaceranger_segmented(seg_dir)


def test_extra_polygon_with_no_matrix_barcode_fails_loudly(tmp_path):
    seg_dir = _make_fixture(tmp_path, n_cells=3, extra_polygon=True)
    with pytest.raises(ValueError, match="match no matrix barcode"):
        load_spaceranger_segmented(seg_dir)


def test_missing_scalefactors_key_fails_loudly(tmp_path):
    seg_dir = _make_fixture(tmp_path, n_cells=2, scalefactors_key="mpp_wrong_name")
    with pytest.raises(KeyError, match="microns_per_pixel"):
        load_spaceranger_segmented(seg_dir)


def test_missing_input_file_fails_loudly(tmp_path):
    seg_dir = _make_fixture(tmp_path, n_cells=2)
    (seg_dir / "cell_segmentations.geojson").unlink()
    with pytest.raises(FileNotFoundError):
        load_spaceranger_segmented(seg_dir)
