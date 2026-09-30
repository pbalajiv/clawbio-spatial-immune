#!/usr/bin/env python3
"""
Plots for CLASSIFICATION_SPEC.md runs (T5): cell-type map, tumour-region
overlay with T cells, biomarker map over tumour cells, and a stats panel.

These are skill outputs, not decoration: every plot is built from exactly the
same intermediate arrays classify_reference.py's classify() computes (same QC
mask, same cell typing, same tumour_region_mask()), so a figure can never show
a cell type or region that disagrees with classification.json. No new
statistics are computed here; this module only visualises what classify()
already decided.

Usage:
    python plots.py --input synthetic/excluded_heterogeneous \
        --classification runs/demo/outputs/classification.json --out runs/demo/figures
"""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import classify_reference as cr

TYPE_COLORS = {
    "tumor": "#c53030",
    "t_cell": "#2b6cb0",
    "stroma": "#2f855a",
    "ambiguous": "#a0aec0",
}


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compute_plot_data(input_path, cfg):
    """Recompute the same QC / typing / region intermediates classify() uses internally.

    Returns None for the cell-typing fields when QC fails outright (matches
    classify_reference.classify()'s qc_fail short-circuit): with too few, possibly
    zero, passing cells there is nothing to type or plot spatially.
    """
    X, genes, xy, provenance = cr.load_input(input_path)
    keep, qcd = cr.qc(X, genes, xy, cfg)
    if qcd["pass_fraction"] < cfg["qc"]["min_pass_fraction"]:
        return {"xy": None, "types": None, "region": None,
                "biomarker_gene": cfg["biomarker"]["gene"], "biomarker_expr": None,
                "qc": qcd, "provenance": provenance}
    X, xy = X[keep], xy[keep]
    L, total = cr.lognorm(X)
    types, scores, double_pos, conc = cr.score_and_type(L, genes, cfg)
    region = cr.tumor_region_mask(xy, types, cfg)
    gene = cfg["biomarker"]["gene"]
    biomarker_expr = cr._gene(L, genes, gene) if gene in genes else None
    return {
        "xy": xy, "types": types, "region": region,
        "biomarker_gene": gene, "biomarker_expr": biomarker_expr,
        "qc": qcd, "provenance": provenance,
    }


def plot_cell_type_map(xy, types, out_path, max_points=40000):
    fig, ax = plt.subplots(figsize=(7, 7))
    rng = np.random.default_rng(0)
    idx = np.arange(len(xy)) if len(xy) <= max_points else rng.choice(len(xy), max_points, replace=False)
    for t, color in TYPE_COLORS.items():
        m = (types == t)[idx] if len(xy) > max_points else (types == t)
        pts = xy[idx][m] if len(xy) > max_points else xy[m]
        ax.scatter(pts[:, 0], pts[:, 1], s=2, c=color, label=f"{t} (n={int((types == t).sum())})", linewidths=0)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("Cell type map")
    ax.legend(markerscale=6, loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_tumor_region_with_tcells(xy, types, region, out_path, max_points=40000):
    fig, ax = plt.subplots(figsize=(7, 7))
    rng = np.random.default_rng(0)
    idx = np.arange(len(xy)) if len(xy) <= max_points else rng.choice(len(xy), max_points, replace=False)
    xy_s, region_s = xy[idx], region[idx]
    ax.scatter(xy_s[~region_s, 0], xy_s[~region_s, 1], s=2, c="#e2e8f0", label="outside tumour region", linewidths=0)
    ax.scatter(xy_s[region_s, 0], xy_s[region_s, 1], s=2, c="#fbb6b6", label="tumour region", linewidths=0)
    t_mask = (types == "t_cell")[idx] if len(xy) > max_points else (types == "t_cell")
    xy_t = xy[idx][t_mask] if len(xy) > max_points else xy[t_mask]
    ax.scatter(xy_t[:, 0], xy_t[:, 1], s=4, c="#2b6cb0", label=f"T cells (n={int((types == 't_cell').sum())})", linewidths=0)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(f"Tumour region ({int(region.sum())} cells) with T cells overlaid")
    ax.legend(markerscale=6, loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_biomarker_map(xy, types, expr, gene, out_path, max_points=40000):
    fig, ax = plt.subplots(figsize=(7, 7))
    tum = types == "tumor"
    xy_t, expr_t = xy[tum], expr[tum]
    if len(xy_t) > max_points:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(xy_t), max_points, replace=False)
        xy_t, expr_t = xy_t[idx], expr_t[idx]
    ax.scatter(xy[~tum, 0], xy[~tum, 1], s=1, c="#edf2f7", linewidths=0, label="non-tumour")
    sc = ax.scatter(xy_t[:, 0], xy_t[:, 1], s=3, c=expr_t, cmap="viridis", linewidths=0)
    fig.colorbar(sc, ax=ax, label=f"{gene} (log-normalised)", shrink=0.7)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(f"{gene} among tumour cells (n={int(tum.sum())})")
    ax.legend(markerscale=6, loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_stats_panel(result, out_path):
    imm, bio, qc = result["immune_phenotype"], result["biomarker"], result["qc"]
    fig, ax = plt.subplots(figsize=(8, 8.5))
    ax.axis("off")
    lines = [
        f"spec_version {result['spec_version']}   seed {result['seed']}",
        "",
        f"QC: {qc['n_cells_passed']}/{qc['n_cells_input']} passed "
        f"({qc['pass_fraction']:.1%}), dropped {qc['n_dropped']} {qc['drop_reasons']}",
        "",
        f"Immune phenotype: {imm['label']}" + (f"  (reason: {imm['reason']})" if imm["label"] == "unclassifiable" else ""),
    ]
    lines += [f"  {s}" for s in imm["rule_path"]]
    lines += [
        "",
        f"Biomarker ({bio['gene']}): {bio['label']}" + (f"  (reason: {bio['reason']})" if bio["label"] == "unclassifiable" else ""),
    ]
    lines += [f"  {s}" for s in bio["rule_path"]]
    warnings = result.get("warnings") or []
    lines += ["", "Warnings:"] + ([f"  {w}" for w in warnings] if warnings else ["  none"])
    ax.text(0.02, 0.98, "\n".join(lines), va="top", ha="left", fontsize=9, family="monospace", wrap=True,
            transform=ax.transAxes)
    ax.set_title("Classification stats panel")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def render_plots(input_path, classification_result, cfg, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data = compute_plot_data(input_path, cfg)

    paths = {}
    if data["xy"] is not None:
        paths["cell_type_map"] = out_dir / "cell_type_map.png"
        plot_cell_type_map(data["xy"], data["types"], paths["cell_type_map"])

        paths["tumor_region_t_cells"] = out_dir / "tumor_region_t_cells.png"
        plot_tumor_region_with_tcells(data["xy"], data["types"], data["region"], paths["tumor_region_t_cells"])

        if data["biomarker_expr"] is not None:
            paths["biomarker_map"] = out_dir / "biomarker_map.png"
            plot_biomarker_map(data["xy"], data["types"], data["biomarker_expr"], data["biomarker_gene"], paths["biomarker_map"])

    paths["stats_panel"] = out_dir / "stats_panel.png"
    plot_stats_panel(classification_result, paths["stats_panel"])

    checksums = {str(p.relative_to(out_dir)): _sha256(p) for p in paths.values()}
    return {str(k): str(v) for k, v in paths.items()}, checksums


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="same input classify_reference.py was run against")
    ap.add_argument("--classification", required=True, help="path to the classification.json this run produced")
    ap.add_argument("--config", default=None, help="config override, if one was used for --classification")
    ap.add_argument("--out", required=True, help="output directory for figures")
    args = ap.parse_args()

    result = json.loads(Path(args.classification).read_text())
    cfg = cr.load_config(args.config)
    paths, checksums = render_plots(args.input, result, cfg, args.out)

    checksum_path = Path(args.out) / "checksums.sha256"
    checksum_path.write_text("\n".join(f"{v}  {k}" for k, v in checksums.items()) + "\n")
    print(f"wrote {len(paths)} figures to {args.out}")
    for name, p in paths.items():
        print(f"  {name}: {p}")


if __name__ == "__main__":
    main()
