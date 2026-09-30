#!/usr/bin/env python3
"""
Synthetic breast-tumour tissues with KNOWN immune phenotype and ERBB2 pattern.

Purpose: ground-truth test cases for the ClawBio classification skills.
Each generated tissue has a label fixed by construction, so the pipeline's
output can be scored as right or wrong.

Layout: tumour nests (discs) on a jittered grid, each surrounded by a stromal
ring, embedded in background stroma. Immune cells are placed with a
scenario-specific probability in each region:

    inflamed : immune cells inside nests and in stroma
    excluded : immune cells concentrated in the stromal ring, nearly absent in nests
    desert   : almost no immune cells anywhere

ERBB2 pattern (tumour cells only):
    homogeneous   : all nests ERBB2-high
    heterogeneous : half the nests ERBB2-high, half ERBB2-low

Extras:
    --low-quality : most cells get very low depth -> expected "unclassifiable"
    --spillover   : immune cells touching tumour cells pick up keratin reads,
                    mimicking segmentation-boundary transcript bleed

Outputs per tissue (the pipeline must only read the first two):
    tissue.h5ad          counts (cells x genes), obsm["spatial"]  [if anndata installed]
    counts.npz, genes.txt, cells.csv   same data without anndata
    truth_cells.csv      true cell types / regions  (NEVER give this to the pipeline)
    ground_truth.json    expected labels + generation parameters
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

GENE_GROUPS = {
    "tumor": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "biomarker": ["ERBB2"],
    "prolif": ["MKI67"],
    "immune": ["PTPRC", "CD3E", "CD8A", "CD8B", "GZMB"],
    "stroma": ["COL1A1", "COL1A2", "DCN", "PDGFRB"],
    "background": ["ACTB", "GAPDH", "RPLP0", "B2M", "MALAT1", "TMSB4X", "EEF1A1", "FTL"],
}
GENES = [g for grp in GENE_GROUPS.values() for g in grp]
GROUP_OF = {g: grp for grp, gs in GENE_GROUPS.items() for g in gs}

# Mean UMIs per gene per cell, by cell type and gene group. Deliberately sparse,
# roughly in the range of cell-level Visium HD output.
BASE_MEANS = {
    "tumor": {"tumor": 6.0, "prolif": 0.8, "immune": 0.05, "stroma": 0.1, "background": 8.0},
    "immune": {"tumor": 0.1, "prolif": 0.4, "immune": 3.0, "stroma": 0.05, "background": 6.0},
    "stroma": {"tumor": 0.1, "prolif": 0.2, "immune": 0.05, "stroma": 5.0, "background": 7.0},
}
ERBB2_HIGH, ERBB2_LOW, ERBB2_NONTUMOR = 6.0, 0.4, 0.1
NB_DISPERSION = 2.0

# Probability that a cell is immune, by region: (inside nest, stromal ring, distant stroma)
SCENARIOS = {
    "inflamed": (0.20, 0.15, 0.08),
    "excluded": (0.005, 0.35, 0.06),
    "desert": (0.003, 0.004, 0.004),
}

SUITE = [
    ("inflamed", "homogeneous", {}),
    ("inflamed", "heterogeneous", {}),
    ("excluded", "homogeneous", {}),
    ("excluded", "heterogeneous", {}),
    ("desert", "homogeneous", {}),
    ("inflamed", "homogeneous", {"low_quality": True}),
]


def _layout(rng, n_cells, size, n_nests, nest_r, ring_w):
    grid = int(np.ceil(np.sqrt(n_nests)))
    step = size / grid
    if nest_r + ring_w >= step / 2:
        raise ValueError("nest_r + ring_w must be < half the grid step so nests don't overlap")
    centers = []
    for i in range(grid):
        for j in range(grid):
            if len(centers) < n_nests:
                jitter = rng.uniform(-0.1, 0.1, 2) * step
                centers.append(((i + 0.5) * step + jitter[0], (j + 0.5) * step + jitter[1]))
    centers = np.array(centers)
    xy = rng.uniform(0, size, (n_cells, 2))
    dist, nest_id = cKDTree(centers).query(xy)
    region = np.where(dist <= nest_r, "nest", np.where(dist <= nest_r + ring_w, "ring", "stroma"))
    return xy, nest_id, region


def _assign_types(rng, region, scenario):
    p_in, p_ring, p_out = SCENARIOS[scenario]
    p_immune = np.select([region == "nest", region == "ring"], [p_in, p_ring], p_out)
    is_immune = rng.random(len(region)) < p_immune
    base = np.where(region == "nest", "tumor", "stroma")
    return np.where(is_immune, "immune", base)


def _nb(rng, mu):
    mu = np.maximum(mu, 0.0)
    return rng.negative_binomial(NB_DISPERSION, NB_DISPERSION / (NB_DISPERSION + mu))


def _expression(rng, types, erbb2_high, depth):
    n = len(types)
    X = np.zeros((n, len(GENES)), dtype=np.int32)
    for j, gene in enumerate(GENES):
        grp = GROUP_OF[gene]
        if grp == "biomarker":
            mu = np.where(types == "tumor", np.where(erbb2_high, ERBB2_HIGH, ERBB2_LOW), ERBB2_NONTUMOR)
        else:
            mu = np.select([types == "tumor", types == "immune"],
                           [BASE_MEANS["tumor"][grp], BASE_MEANS["immune"][grp]],
                           BASE_MEANS["stroma"][grp])
        X[:, j] = _nb(rng, mu * depth)
    return X


def _add_spillover(rng, X, xy, types, frac, radius=10.0):
    """Immune cells within `radius` um of a tumour cell gain a fraction of tumour keratin reads."""
    tumor_idx = np.where(types == "tumor")[0]
    immune_idx = np.where(types == "immune")[0]
    if len(tumor_idx) == 0 or len(immune_idx) == 0:
        return X, np.zeros(len(types), bool)
    d, _ = cKDTree(xy[tumor_idx]).query(xy[immune_idx])
    touching = immune_idx[d <= radius]
    for gene in GENE_GROUPS["tumor"]:
        j = GENES.index(gene)
        X[touching, j] += rng.poisson(BASE_MEANS["tumor"]["tumor"] * frac, len(touching))
    flag = np.zeros(len(types), bool)
    flag[touching] = True
    return X, flag


def morans_i(values, xy, k=6):
    """Moran's I with row-standardised kNN weights (no significance test; for self-check only)."""
    n = len(values)
    if n <= k:
        return float("nan")
    _, nbrs = cKDTree(xy).query(xy, k=k + 1)
    nbrs = nbrs[:, 1:]
    z = values - values.mean()
    denom = (z ** 2).sum()
    if denom == 0:
        return float("nan")
    return float((z * z[nbrs].mean(axis=1)).sum() / denom)


def self_check(X, xy, types, region):
    """Summary metrics computed from TRUE labels, to confirm scenarios are separable."""
    immune = types == "immune"
    tumor = types == "tumor"
    out = {
        "immune_fraction": float(immune.mean()),
        "immune_in_nest_fraction": float((immune & (region == "nest")).sum() / max(immune.sum(), 1)),
    }
    if immune.sum() and tumor.sum():
        _, nbrs = cKDTree(xy).query(xy[immune], k=7)
        out["immune_with_tumor_neighbour"] = float((types[nbrs[:, 1:]] == "tumor").any(axis=1).mean())
    erbb2 = np.log1p(X[:, GENES.index("ERBB2")].astype(float))
    out["erbb2_morans_i_within_tumor"] = morans_i(erbb2[tumor], xy[tumor])
    return out


def generate(scenario, her2, out_dir, n_cells=20000, size=1000.0, n_nests=9, nest_r=110.0,
             ring_w=40.0, seed=0, low_quality=False, spillover=0.0):
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {list(SCENARIOS)}")
    if her2 not in ("homogeneous", "heterogeneous"):
        raise ValueError("her2 must be 'homogeneous' or 'heterogeneous'")
    rng = np.random.default_rng(seed)
    xy, nest_id, region = _layout(rng, n_cells, size, n_nests, nest_r, ring_w)
    types = _assign_types(rng, region, scenario)

    if her2 == "homogeneous":
        high_nests = np.arange(n_nests)
    else:
        high_nests = rng.permutation(n_nests)[: n_nests // 2]
    erbb2_high = np.isin(nest_id, high_nests)

    depth = rng.lognormal(0, 0.35, n_cells)
    low_depth = np.zeros(n_cells, bool)
    if low_quality:
        low_depth = rng.random(n_cells) < 0.7
        depth[low_depth] *= 0.05
    X = _expression(rng, types, erbb2_high, depth)
    spill_flag = np.zeros(n_cells, bool)
    if spillover > 0:
        X, spill_flag = _add_spillover(rng, X, xy, types, spillover)

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cell_ids = [f"cell_{i:06d}" for i in range(n_cells)]
    Xs = sparse.csr_matrix(X)

    sparse.save_npz(out / "counts.npz", Xs)
    (out / "genes.txt").write_text("\n".join(GENES) + "\n")
    pd.DataFrame({"cell_id": cell_ids, "x": xy[:, 0], "y": xy[:, 1]}).to_csv(out / "cells.csv", index=False)
    pd.DataFrame({
        "cell_id": cell_ids, "true_type": types, "region": region, "nest_id": nest_id,
        "erbb2_high_nest": erbb2_high & (types == "tumor"), "low_depth": low_depth,
        "spillover_affected": spill_flag,
    }).to_csv(out / "truth_cells.csv", index=False)

    try:
        import anndata as ad
        adata = ad.AnnData(X=Xs.astype(np.float32),
                           obs=pd.DataFrame(index=cell_ids),
                           var=pd.DataFrame(index=GENES))
        adata.obsm["spatial"] = xy
        adata.write_h5ad(out / "tissue.h5ad")
        wrote_h5ad = True
    except ImportError:
        wrote_h5ad = False

    expected_immune = "unclassifiable" if low_quality else scenario
    expected_pattern = "unclassifiable" if low_quality else her2
    truth = {
        "expected_immune_phenotype": expected_immune,
        "expected_biomarker_pattern": expected_pattern,
        "biomarker_gene": "ERBB2",
        "params": {
            "scenario": scenario, "her2": her2, "n_cells": n_cells, "size_um": size,
            "n_nests": n_nests, "nest_r_um": nest_r, "ring_w_um": ring_w, "seed": seed,
            "low_quality": low_quality, "spillover": spillover,
            "immune_probabilities_nest_ring_stroma": SCENARIOS[scenario],
        },
        "self_check_from_true_labels": self_check(X, xy, types, region),
        "files": {"h5ad": wrote_h5ad},
    }
    (out / "ground_truth.json").write_text(json.dumps(truth, indent=2))
    return truth


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--suite", action="store_true", help="generate the standard 6-tissue test suite")
    ap.add_argument("--scenario", choices=list(SCENARIOS), default="inflamed")
    ap.add_argument("--her2", choices=["homogeneous", "heterogeneous"], default="homogeneous")
    ap.add_argument("--n-cells", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--low-quality", action="store_true")
    ap.add_argument("--spillover", type=float, default=0.0,
                    help="fraction of tumour keratin mean added to immune cells touching tumour (e.g. 0.3)")
    args = ap.parse_args()

    runs = SUITE if args.suite else [(args.scenario, args.her2, {"low_quality": args.low_quality})]
    for scenario, her2, extra in runs:
        name = f"{scenario}_{her2}" + ("_lowqc" if extra.get("low_quality") else "") + \
               (f"_spill{args.spillover}" if args.spillover else "")
        out = Path(args.out) / name if args.suite else Path(args.out)
        t = generate(scenario, her2, out, n_cells=args.n_cells, seed=args.seed,
                     spillover=args.spillover, **extra)
        sc = t["self_check_from_true_labels"]
        print(f"{name:40s} immune={t['expected_immune_phenotype']:15s} "
              f"erbb2={t['expected_biomarker_pattern']:15s} "
              f"immune_frac={sc['immune_fraction']:.3f} "
              f"in_nest={sc['immune_in_nest_fraction']:.3f} "
              f"moranI={sc['erbb2_morans_i_within_tumor']:.3f}")


if __name__ == "__main__":
    main()
