#!/usr/bin/env python3
"""
Reference implementation of CLASSIFICATION_SPEC.md (spec v1.0.0).

Deterministic given (input, config, seed). Numpy/scipy only, so it can be checked
line by line against the spec. The ClawBio skill may use Scanpy/Squidpy internally,
but must reproduce this implementation's labels on the synthetic suite
(see test_classify_reference.py).

Usage:
    python classify_reference.py --input synthetic/excluded_heterogeneous --out run/outputs/classification.json
    python classify_reference.py --input tissue.h5ad --config classification_config.json --out ...
"""
import argparse
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

SPEC_VERSION = "1.0.0"
HERE = Path(__file__).parent


# ----------------------------------------------------------------------------- io

def load_config(path=None):
    cfg = json.loads((HERE / "classification_config.json").read_text())
    if path:
        cfg = _deep_update(cfg, json.loads(Path(path).read_text()))
    bm = cfg["biomarker"]["gene"]
    if bm in cfg["signatures"]["tumor"]:
        raise ValueError(f"{bm} must not be in the tumor signature (circularity; see spec section 3)")
    return cfg


def _deep_update(base, upd):
    for k, v in upd.items():
        base[k] = _deep_update(base.get(k, {}), v) if isinstance(v, dict) else v
    return base


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_input(path):
    """Return (counts csr [cells x genes], genes list, xy in micrometres, provenance dict)."""
    p = Path(path)
    if p.is_dir() and (p / "counts.npz").exists():           # synthetic_tissue.py output
        import pandas as pd
        X = sparse.load_npz(p / "counts.npz").tocsr()
        genes = (p / "genes.txt").read_text().split()
        xy = pd.read_csv(p / "cells.csv")[["x", "y"]].to_numpy(float)
        return X, genes, xy, {"path": str(p / "counts.npz"), "sha256": _sha256(p / "counts.npz"),
                              "format": "synthetic_npz"}
    if p.suffix == ".h5ad":                                   # e.g. output of the Space Ranger loader
        import anndata as ad
        a = ad.read_h5ad(p)
        units = a.uns.get("spatial_units", "um")
        if units != "um":
            raise ValueError(f"obsm['spatial'] must be in micrometres, got uns['spatial_units']={units!r}")
        return (sparse.csr_matrix(a.X), list(a.var_names), np.asarray(a.obsm["spatial"], float),
                {"path": str(p), "sha256": _sha256(p), "format": "h5ad"})
    raise ValueError("input must be a synthetic tissue directory or an .h5ad with obsm['spatial'] in um")


# ----------------------------------------------------------------------------- steps

def _mitochondrial_fraction(X, genes, total):
    """Fraction of each cell's counts from MT-prefixed genes, or None if the panel has none."""
    mt_idx = [i for i, g in enumerate(genes) if g.upper().startswith("MT-")]
    if not mt_idx:
        return None
    mt_counts = np.asarray(X[:, mt_idx].sum(1)).ravel()
    return mt_counts / np.maximum(total, 1)


def local_outlier_qc(xy, total, ngenes, mito_frac, cfg):
    """Spatially-local outlier QC (SpotSweeper-style; Totty, Hicks & Guo 2025).

    Flags a cell as an outlier only relative to its k spatially-nearest
    neighbours, not against one global threshold across the whole section.
    In spatial data, library size / gene count / mitochondrial fraction can
    shift for real biological reasons across a tissue section (e.g. a
    naturally lower-RNA-content region), which a single global cutoff cannot
    distinguish from genuine technical failure (see OSTA book, QC chapter).

    log1p-transforms library size and gene count (assumed log-normal);
    leaves mitochondrial fraction untransformed (assumed already ~normal),
    matching SpotSweeper's convention. Returns three boolean masks:
    (low_libsize, low_genes, high_mito). high_mito is all-False if mito_frac
    is None (no MT- genes in the panel).
    """
    sq = cfg["spatial_qc"]
    nb = knn(xy, sq["k"], cfg["graph"]["max_edge_um"])
    valid = nb >= 0
    z_thresh = sq["z_threshold"]

    def local_z(values, log_transform):
        v = np.log1p(values) if log_transform else np.asarray(values, dtype=float)
        safe_nb = np.where(valid, nb, 0)
        neighbour_vals = np.where(valid, v[safe_nb], np.nan)
        # A cell with zero valid neighbours (isolated beyond max_edge_um) gives an
        # all-NaN row here -- nanmean/nanstd warn about that expected, handled case
        # (local_mean/local_sd come back NaN, z becomes NaN below, mapped to 0 = not
        # flagged, since there's no local context to judge it against).
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            local_mean = np.nanmean(neighbour_vals, axis=1)
            local_sd = np.nanstd(neighbour_vals, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            z = (v - local_mean) / local_sd
        # NaN means "can't evaluate" (no valid neighbours, or v == a zero-variance
        # neighbourhood) -> don't flag. +-inf means "differs from a *zero-variance*
        # neighbourhood" -- a genuine, extreme local outlier -- so clip to a large
        # finite value rather than silently zeroing it out.
        z = np.where(np.isnan(z), 0.0, z)
        return np.clip(z, -1e6, 1e6)

    low_libsize = local_z(total, log_transform=True) <= -z_thresh
    low_genes = local_z(ngenes, log_transform=True) <= -z_thresh
    if mito_frac is None:
        high_mito = np.zeros(len(total), dtype=bool)
    else:
        high_mito = local_z(mito_frac, log_transform=False) >= z_thresh
    return low_libsize, low_genes, high_mito


def qc(X, genes, xy, cfg):
    total = np.asarray(X.sum(1)).ravel()
    ngenes = np.asarray((X > 0).sum(1)).ravel()
    low_umi = total < cfg["qc"]["min_umi"]
    low_genes = (~low_umi) & (ngenes < cfg["qc"]["min_genes"])
    dropped = low_umi | low_genes
    drop_reasons = {"low_umi": int(low_umi.sum()), "low_genes": int(low_genes.sum())}

    # Spatially-local outlier QC is opt-in: only runs when the config carries a
    # "spatial_qc" section (absent from the frozen, synthetic-calibrated v1.0.0
    # config, so synthetic-suite behaviour is unchanged unless explicitly enabled
    # via a config override -- same pattern as the real-data qc.min_umi override).
    if "spatial_qc" in cfg:
        mito_frac = _mitochondrial_fraction(X, genes, total)
        low_libsize_sp, low_genes_sp, high_mito_sp = local_outlier_qc(xy, total, ngenes, mito_frac, cfg)
        low_libsize_sp = low_libsize_sp & ~dropped
        low_genes_sp = low_genes_sp & ~dropped & ~low_libsize_sp
        high_mito_sp = high_mito_sp & ~dropped & ~low_libsize_sp & ~low_genes_sp
        dropped = dropped | low_libsize_sp | low_genes_sp | high_mito_sp
        drop_reasons.update({
            "spatial_outlier_low_libsize": int(low_libsize_sp.sum()),
            "spatial_outlier_low_genes": int(low_genes_sp.sum()),
            "spatial_outlier_high_mito": int(high_mito_sp.sum()),
        })

    keep = ~dropped
    return keep, {"n_cells_input": int(len(keep)), "n_cells_passed": int(keep.sum()),
                  "n_dropped": int((~keep).sum()),
                  "drop_reasons": drop_reasons,
                  "pass_fraction": float(keep.mean())}


def lognorm(X):
    total = np.asarray(X.sum(1)).ravel()
    Xn = sparse.diags(1e4 / np.maximum(total, 1)) @ X
    return Xn.log1p().tocsc(), total


def _gene(L, genes, g):
    if g not in genes:
        raise KeyError(f"gene {g} not in matrix")
    return np.asarray(L[:, genes.index(g)].todense()).ravel()


def _scaled(v, q):
    top = np.quantile(v, q)
    return np.clip(v / top, 0, 1) if top > 0 else np.zeros_like(v)


def score_and_type(L, genes, cfg):
    q = cfg["cell_typing"]["scale_quantile"]
    sigs = cfg["signatures"]
    scaled = {s: {g: _scaled(_gene(L, genes, g), q) for g in gs} for s, gs in sigs.items()}
    scores = {s: np.mean(list(d.values()), axis=0) for s, d in scaled.items()}
    names = list(sigs)
    M = np.vstack([scores[s] for s in names]).T
    srt = np.sort(M, axis=1)
    best = np.array(names)[M.argmax(1)]
    ok = (srt[:, -1] >= cfg["cell_typing"]["min_score"]) & \
         (srt[:, -1] - srt[:, -2] >= cfg["cell_typing"]["min_margin"])
    types = np.where(ok, best, "ambiguous")
    mn = cfg["cell_typing"]["min_score"]
    double_pos = (scores["tumor"] >= mn) & (scores["t_cell"] >= mn)

    conc = {}
    for s, d in scaled.items():
        rho = {}
        for g in d:
            loo = np.mean([v for h, v in d.items() if h != g], axis=0)
            r = spearmanr(_gene(L, genes, g), loo)[0]
            rho[g] = float(0.0 if np.isnan(r) else r)
        med = float(np.median(list(rho.values())))
        conc[s] = {"marker_loo_spearman": rho, "median": med,
                   "trusted": med >= cfg["concordance"]["min_median_loo_spearman"]}
    return types, scores, double_pos, conc


def knn(xy, k, max_edge_um):
    """Neighbour index array (n, k); entries of -1 mark edges longer than max_edge_um or missing."""
    kk = min(k + 1, len(xy))
    d, nb = cKDTree(xy).query(xy, k=kk, distance_upper_bound=max_edge_um)
    d, nb = d[:, 1:], nb[:, 1:]
    nb = np.where(np.isfinite(d), nb, -1)
    if nb.shape[1] < k:
        nb = np.hstack([nb, -np.ones((len(xy), k - nb.shape[1]), int)])
    return nb


def tumor_region_mask(xy, types, cfg):
    """Boolean mask: cell has >= region_min_tumor_neighbour_fraction tumour cells among its region_k neighbours.

    Exposed (not just inlined in classify()) so plotting code can reproduce the exact
    tumour-region definition used by the immune rule, rather than re-deriving it.
    """
    ic = cfg["immune_rule"]
    nb_region = knn(xy, ic["region_k"], cfg["graph"]["max_edge_um"])
    valid = nb_region >= 0
    tum_mask = types == "tumor"
    frac_tum_nbrs = (tum_mask[np.where(valid, nb_region, 0)] & valid).sum(1) / np.maximum(valid.sum(1), 1)
    return frac_tum_nbrs >= ic["region_min_tumor_neighbour_fraction"]


def _edges(a, b, nb, valid):
    """Directed kNN edges between label sets a and b, counted in both directions."""
    safe = np.where(valid, nb, 0)
    return int(((a[:, None] & b[safe]) & valid).sum() + ((b[:, None] & a[safe]) & valid).sum())


def nhood_enrichment(types, nb, n_perm, rng):
    valid = nb >= 0
    tum, t = types == "tumor", types == "t_cell"
    obs = _edges(tum, t, nb, valid)
    null = np.empty(n_perm)
    for i in range(n_perm):
        p = rng.permutation(len(types))
        null[i] = _edges(tum[p], t[p], nb, valid)
    sd = null.std()
    z = float((obs - null.mean()) / sd) if sd > 0 else 0.0
    p = float((1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum()) / (n_perm + 1))
    return z, p


def morans_i(v, nb):
    valid = nb >= 0
    safe = np.where(valid, nb, 0)
    cnt = valid.sum(1)
    z = v - v.mean()
    denom = (z ** 2).sum()
    if denom == 0:
        return float("nan")
    lag = np.where(cnt > 0, (z[safe] * valid).sum(1) / np.maximum(cnt, 1), 0.0)
    return float((z * lag).sum() / denom)


def morans_perm(v, nb, n_perm, rng):
    obs = morans_i(v, nb)
    null = np.array([morans_i(v[rng.permutation(len(v))], nb) for _ in range(n_perm)])
    p = float((1 + (null >= obs).sum()) / (n_perm + 1))   # one-sided: positive autocorrelation
    return obs, p


def spillover_check(scores, types, nb, cfg):
    """Double-positive rate among T-like cells adjacent to tumour vs distant from tumour.

    Transcript spillover at segmentation boundaries only affects cells touching tumour, so a
    ratio well above 1 flags it; a uniform background of double positives does not.
    """
    mn = cfg["cell_typing"]["min_score"]
    t_like = scores["t_cell"] >= cfg["warnings"]["spillover_t_like_min_score"]
    dp = t_like & (scores["tumor"] >= mn)
    valid = nb >= 0
    adj = ((types == "tumor")[np.where(valid, nb, 0)] & valid).any(1)
    a, d = t_like & adj, t_like & ~adj
    min_n = cfg["warnings"]["spillover_min_cells"]
    ra = float(dp[a].mean()) if a.sum() else None
    rd = float(dp[d].mean()) if d.sum() else None
    ratio = float(ra / rd) if (a.sum() >= min_n and d.sum() >= min_n and rd) else None
    return {"n_t_like_adjacent": int(a.sum()), "n_t_like_distant": int(d.sum()),
            "double_positive_rate_adjacent": ra, "double_positive_rate_distant": rd, "ratio": ratio}


# ----------------------------------------------------------------------------- rule

def classify(X, genes, xy, cfg, seed, provenance=None):
    rng = np.random.default_rng(seed)
    out = {"schema_version": "1.0.0", "spec_version": SPEC_VERSION, "seed": seed, "config": cfg,
           "input": provenance or {}, "warnings": []}
    keep, qcd = qc(X, genes, xy, cfg)
    out["qc"] = qcd
    qc_failed = qcd["pass_fraction"] < cfg["qc"]["min_pass_fraction"]

    # Below min_pass_fraction, both calls abstain immediately (spec section 6/7, rule 1).
    # Skip cell typing and spatial statistics entirely: with too few (possibly zero)
    # passing cells, per-signature quantile scaling has nothing to operate on.
    if qc_failed:
        n_t = n_tum = n = 0
        out["cell_types"] = {"counts": {k: 0 for k in list(cfg["signatures"]) + ["ambiguous"]},
                              "n_double_positive_tumor_t": 0}
        out["signature_concordance"] = {s: {"marker_loo_spearman": {}, "median": 0.0, "trusted": False}
                                        for s in cfg["signatures"]}
        out["spillover_check"] = {"n_t_like_adjacent": 0, "n_t_like_distant": 0,
                                  "double_positive_rate_adjacent": None,
                                  "double_positive_rate_distant": None, "ratio": None}
    else:
        X, xy = X[keep], xy[keep]
        L, total = lognorm(X)
        types, scores, double_pos, conc = score_and_type(L, genes, cfg)
        uniq, cnt = np.unique(types, return_counts=True)
        counts = {k: 0 for k in list(cfg["signatures"]) + ["ambiguous"]}
        counts.update({str(k): int(v) for k, v in zip(uniq, cnt)})
        n_t, n_tum, n = counts["t_cell"], counts["tumor"], len(types)
        out["cell_types"] = {"counts": counts, "n_double_positive_tumor_t": int(double_pos.sum())}
        out["signature_concordance"] = conc
        nb = knn(xy, cfg["graph"]["k"], cfg["graph"]["max_edge_um"])
        out["spillover_check"] = spillover_check(scores, types, nb, cfg)
        sp = out["spillover_check"]
        if sp["ratio"] is not None and sp["ratio"] > cfg["warnings"]["max_spillover_ratio"]:
            out["warnings"].append(f"spillover_suspected: T-like cells touching tumour are double-positive "
                                   f"{sp['ratio']:.2f}x as often as distant T-like cells")

    conc = out["signature_concordance"]
    imm = {"label": None, "reason": None, "rule_path": [], "n_cells_passed": n, "n_tumor_cells": n_tum,
           "n_t_cells": n_t, "t_cell_fraction": float(n_t / max(n, 1)), "n_tumor_region_cells": None,
           "infiltration_fraction": None, "infiltration_ratio": None, "nhood_z": None, "nhood_p": None,
           "n_permutations": cfg["permutations"]["n"]}
    ic = cfg["immune_rule"]

    def done(label, reason=None):
        imm["label"], imm["reason"] = label, reason

    if qc_failed:
        imm["rule_path"].append(f"pass_fraction {qcd['pass_fraction']:.3f} < {cfg['qc']['min_pass_fraction']}")
        done("unclassifiable", "qc_fail")
    elif n_tum < ic["min_tumor_cells"]:
        imm["rule_path"].append(f"n_tumor_cells {n_tum} < {ic['min_tumor_cells']}")
        done("unclassifiable", "insufficient_tumor_cells")
    elif not conc["tumor"]["trusted"]:
        imm["rule_path"].append("tumor signature concordance below threshold")
        done("unclassifiable", "tumor_signature_discordant")
    else:
        imm["rule_path"].append("qc_ok, tumor signature trusted")
        region = tumor_region_mask(xy, types, cfg)
        t_mask = types == "t_cell"
        imm["n_tumor_region_cells"] = int(region.sum())
        z, p = nhood_enrichment(types, nb, cfg["permutations"]["n"], rng)
        imm["nhood_z"], imm["nhood_p"] = z, p
        if n_t > 0:
            imm["infiltration_fraction"] = float((t_mask & region).sum() / n_t)
        if region.sum() > 0 and n_t > 0:
            imm["infiltration_ratio"] = float(((t_mask & region).sum() / region.sum()) / (n_t / n))

        tf = imm["t_cell_fraction"]
        if tf < ic["desert_max_t_fraction"]:
            imm["rule_path"].append(f"t_cell_fraction {tf:.4f} < {ic['desert_max_t_fraction']}")
            done("desert")
        elif n_t < ic["min_t_cells"]:
            imm["rule_path"].append(f"n_t_cells {n_t} < {ic['min_t_cells']}")
            done("unclassifiable", "insufficient_t_cells")
        elif not conc["t_cell"]["trusted"]:
            imm["rule_path"].append("t_cell signature concordance below threshold")
            done("unclassifiable", "t_cell_signature_discordant")
        elif imm["infiltration_ratio"] is None:
            imm["rule_path"].append("no tumour-region cells")
            done("unclassifiable", "no_tumor_region")
        else:
            r = imm["infiltration_ratio"]
            imm["rule_path"].append(f"t_cell_fraction {tf:.4f} >= {ic['desert_max_t_fraction']}, t_cell signature trusted")
            if r >= ic["inflamed_min_ratio"]:
                imm["rule_path"].append(f"infiltration_ratio {r:.3f} >= {ic['inflamed_min_ratio']}")
                done("inflamed")
                if z <= -ic["excluded_max_z"]:
                    out["warnings"].append(f"nhood_z_conflict: inflamed by infiltration ratio but nhood_z={z:.2f}")
            elif r <= ic["excluded_max_ratio"] and z <= ic["excluded_max_z"] and p <= ic["alpha"]:
                imm["rule_path"].append(f"infiltration_ratio {r:.3f} <= {ic['excluded_max_ratio']} and "
                                        f"nhood_z {z:.2f} <= {ic['excluded_max_z']} (p={p:.4f})")
                done("excluded")
            else:
                imm["rule_path"].append(f"infiltration_ratio {r:.3f} / nhood_z {z:.2f} in the borderline zone")
                done("unclassifiable", "borderline")
    out["immune_phenotype"] = imm

    # ---- biomarker
    bc = cfg["biomarker"]
    bio = {"gene": bc["gene"], "label": None, "reason": None, "rule_path": [], "n_tumor_cells": n_tum,
           "detection_fraction": None, "morans_i": None, "morans_p": None, "morans_i_depth_control": None,
           "morans_i_excess": None, "n_permutations": cfg["permutations"]["n"]}

    def bdone(label, reason=None):
        bio["label"], bio["reason"] = label, reason

    if qcd["pass_fraction"] < cfg["qc"]["min_pass_fraction"]:
        bio["rule_path"].append("qc_fail")
        bdone("unclassifiable", "qc_fail")
    elif n_tum < bc["min_tumor_cells"]:
        bio["rule_path"].append(f"n_tumor_cells {n_tum} < {bc['min_tumor_cells']}")
        bdone("unclassifiable", "insufficient_tumor_cells")
    else:
        tm = types == "tumor"
        v = _gene(L, genes, bc["gene"])[tm]
        bio["detection_fraction"] = float((v > 0).mean())
        if bio["detection_fraction"] < bc["min_detection_fraction"]:
            bio["rule_path"].append(f"detection_fraction {bio['detection_fraction']:.3f} < {bc['min_detection_fraction']}")
            bdone("unclassifiable", "biomarker_not_expressed")
        else:
            nbt = knn(xy[tm], cfg["graph"]["k"], cfg["graph"]["max_edge_um"])
            i_obs, p = morans_perm(v, nbt, cfg["permutations"]["n"], rng)
            i_depth = morans_i(np.log(total[tm]), nbt)
            excess = i_obs - max(i_depth, 0.0)
            bio.update(morans_i=i_obs, morans_p=p, morans_i_depth_control=i_depth, morans_i_excess=excess)
            if i_depth > cfg["warnings"]["max_depth_morans_i"]:
                out["warnings"].append(f"depth_autocorrelation: Moran's I of log total counts = {i_depth:.3f}")
            if excess >= bc["heterogeneous_min_excess"] and p <= bc["alpha"]:
                bio["rule_path"].append(f"morans_i_excess {excess:.3f} >= {bc['heterogeneous_min_excess']} (p={p:.4f})")
                bdone("heterogeneous")
            elif excess < bc["homogeneous_max_excess"]:
                bio["rule_path"].append(f"morans_i_excess {excess:.3f} < {bc['homogeneous_max_excess']}")
                bdone("homogeneous")
            else:
                bio["rule_path"].append(f"morans_i_excess {excess:.3f} in the borderline zone")
                bdone("unclassifiable", "borderline")
    out["biomarker"] = bio
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--config")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    X, genes, xy, prov = load_input(args.input)
    res = classify(X, genes, xy, cfg, args.seed, prov)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(res, indent=2))
    print(f"immune: {res['immune_phenotype']['label']} ({res['immune_phenotype']['reason']})  "
          f"{res['biomarker']['gene']}: {res['biomarker']['label']} ({res['biomarker']['reason']})")
    for w in res["warnings"]:
        print("WARNING", w)


if __name__ == "__main__":
    main()
