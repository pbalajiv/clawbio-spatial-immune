#!/usr/bin/env python3
"""
Space Ranger v4 segmented-output loader (CLASSIFICATION_SPEC.md section 1, input type 1).

Reads `filtered_feature_cell_matrix.h5`, `cell_segmentations.geojson` and
`spatial/scalefactors_json.json` from a Space Ranger v4 `segmented_outputs/`
directory and writes an `.h5ad` with `obsm["spatial"]` in micrometres,
`uns["spatial_units"] = "um"`, and the input file checksums, ready to be
consumed by `classify_reference.py` (input type 2).

Usage:
    python load_spaceranger_segmented.py --input data/raw/segmented_outputs --out data/tissue.h5ad
"""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
from scipy import sparse


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_matrix_h5(path):
    """Return (X csr [cells x genes], gene symbols, barcodes) from a 10x feature-cell matrix h5."""
    with h5py.File(path, "r") as f:
        grp = f["matrix"]
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
        n_genes, n_barcodes = grp["shape"][:]
        genes = [g.decode() for g in grp["features"]["name"][:]]
        barcodes = [b.decode() for b in grp["barcodes"][:]]
    if len(genes) != n_genes or len(barcodes) != n_barcodes:
        raise ValueError(
            f"{path}: shape {(n_genes, n_barcodes)} does not match "
            f"{len(genes)} feature names / {len(barcodes)} barcodes"
        )
    # 10x h5 stores a genes x barcodes CSC matrix; transpose to cells x genes.
    X_gc = sparse.csc_matrix((data, indices, indptr), shape=(n_genes, n_barcodes))
    return X_gc.T.tocsr(), genes, barcodes


def _polygon_centroid_and_area(ring):
    """Area-weighted centroid and unsigned area of a closed polygon ring [[x, y], ...]."""
    coords = np.asarray(ring, dtype=float)
    if not np.allclose(coords[0], coords[-1]):
        coords = np.vstack([coords, coords[0]])
    x, y = coords[:, 0], coords[:, 1]
    x1, y1, x2, y2 = x[:-1], y[:-1], x[1:], y[1:]
    cross = x1 * y2 - x2 * y1
    area2 = cross.sum()
    if abs(area2) < 1e-9:  # degenerate polygon: fall back to the vertex mean, zero area
        return float(x[:-1].mean()), float(y[:-1].mean()), 0.0
    cx = ((x1 + x2) * cross).sum() / (3 * area2)
    cy = ((y1 + y2) * cross).sum() / (3 * area2)
    return float(cx), float(cy), float(abs(area2) / 2)


def _load_geojson_centroids(path):
    """Return {barcode: (x_px, y_px, area_px2)} from cell_segmentations.geojson, keyed by mapped barcode."""
    doc = json.loads(Path(path).read_text())
    out = {}
    for feat in doc["features"]:
        cell_id = feat["properties"]["cell_id"]
        barcode = f"cellid_{cell_id:09d}-1"
        geom = feat["geometry"]
        if geom["type"] == "Polygon":
            ring = geom["coordinates"][0]
        elif geom["type"] == "MultiPolygon":
            ring = max((poly[0] for poly in geom["coordinates"]), key=len)
        else:
            raise ValueError(f"cell_id {cell_id}: unsupported geometry type {geom['type']!r}")
        if barcode in out:
            raise ValueError(f"duplicate cell_id {cell_id} maps to barcode {barcode} already seen")
        out[barcode] = _polygon_centroid_and_area(ring)
    return out


def _microns_per_pixel(scalefactors_path):
    scalefactors = json.loads(Path(scalefactors_path).read_text())
    if "microns_per_pixel" not in scalefactors:
        raise KeyError(
            f"{scalefactors_path}: expected key 'microns_per_pixel' not found; "
            f"available keys: {sorted(scalefactors)}"
        )
    return scalefactors["microns_per_pixel"], scalefactors


def load_spaceranger_segmented(segmented_outputs_dir, seed=0):
    """Load a Space Ranger v4 segmented_outputs/ directory into an AnnData object."""
    seg_dir = Path(segmented_outputs_dir)
    matrix_path = seg_dir / "filtered_feature_cell_matrix.h5"
    geojson_path = seg_dir / "cell_segmentations.geojson"
    scalefactors_path = seg_dir / "spatial" / "scalefactors_json.json"
    for p in (matrix_path, geojson_path, scalefactors_path):
        if not p.exists():
            raise FileNotFoundError(f"required Space Ranger v4 file missing: {p}")

    X, genes, barcodes = _read_matrix_h5(matrix_path)
    centroids_px = _load_geojson_centroids(geojson_path)
    um_per_px, scalefactors = _microns_per_pixel(scalefactors_path)

    missing = [b for b in barcodes if b not in centroids_px]
    if missing:
        raise ValueError(
            f"{len(missing)} of {len(barcodes)} matrix cells have no polygon in "
            f"{geojson_path} (e.g. {missing[:5]})"
        )
    extra = sorted(set(centroids_px) - set(barcodes))
    if extra:
        raise ValueError(
            f"{len(extra)} polygons in {geojson_path} match no matrix barcode "
            f"(e.g. {extra[:5]})"
        )

    xy_area_px = np.array([centroids_px[b] for b in barcodes], dtype=float)
    xy_px, area_px2 = xy_area_px[:, :2], xy_area_px[:, 2]
    xy_um = xy_px * um_per_px
    area_um2 = area_px2 * (um_per_px ** 2)

    adata = ad.AnnData(X=X)
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.obs["cell_area_um2"] = area_um2
    adata.obsm["spatial"] = xy_um
    adata.uns["spatial_units"] = "um"
    adata.uns["seed"] = int(seed)
    adata.uns["source"] = {
        "loader": "load_spaceranger_segmented.py",
        "segmented_outputs_dir": str(seg_dir),
        "microns_per_pixel": um_per_px,
        "scalefactors": scalefactors,
        "input_checksums": {
            "filtered_feature_cell_matrix.h5": _sha256(matrix_path),
            "cell_segmentations.geojson": _sha256(geojson_path),
            "scalefactors_json.json": _sha256(scalefactors_path),
        },
    }
    return adata


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="path to a Space Ranger v4 segmented_outputs/ directory")
    ap.add_argument("--out", required=True, help="output .h5ad path")
    ap.add_argument("--seed", type=int, default=0, help="seed to log with the output (no randomness in this loader)")
    args = ap.parse_args()

    adata = load_spaceranger_segmented(args.input, seed=args.seed)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)
    print(f"wrote {out_path}: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
