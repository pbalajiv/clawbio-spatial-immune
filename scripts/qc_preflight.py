#!/usr/bin/env python3
"""
T3 real-data pre-flight (CLASSIFICATION_SPEC.md section 0 and section 2).

Run before any phenotype metric on real data:
  1. confirm all signature genes and the biomarker gene are present, and report
     their per-gene detection rate (a gene can be "present" but barely detected),
  2. plot UMI, genes-per-cell, mitochondrial fraction, cell-area, and a
     lightweight doublet-candidate heuristic,
  3. report candidate qc.min_umi / qc.min_genes for whole-transcriptome data.

Mitochondrial fraction, the doublet heuristic, and cell area are REPORTING ONLY:
they are not applied as QC gates here. Any of them becoming a blocking filter is
a spec change (CLASSIFICATION_SPEC.md section 12) applied separately, with the
user's sign-off, after this report is reviewed -- the same pattern already used
for qc.min_umi/min_genes and the t_cell signature on this dataset.

The doublet heuristic is deliberately simple (no scikit-learn / Scrublet): a
cell is flagged as a doublet candidate if its total UMI AND its detected-gene
count both exceed a multiple of the cohort median. This is cruder than a
trained classifier but consistent with this project's numpy/scipy-only,
audit-by-reading design, and is meant to say "here's roughly how many cells
look suspicious," not to make a validated per-cell doublet call.

Usage:
    python scripts/qc_preflight.py --h5ad data/tissue.h5ad --out data/qc_preflight.png
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent.parent


def check_genes(adata, cfg):
    """Presence AND detection rate per signature/biomarker gene (a gene can be present but barely detected)."""
    present = set(adata.var_names)
    X = adata.X
    n = adata.n_obs

    def detection_rate(gene):
        if gene not in present:
            return None
        col = X[:, list(adata.var_names).index(gene)]
        return float((np.asarray(col.todense()).ravel() > 0).sum() / n) if n else 0.0

    report = {}
    for sig, genes in cfg["signatures"].items():
        report[sig] = {"genes": genes,
                       "missing": [g for g in genes if g not in present],
                       "detection_rate": {g: detection_rate(g) for g in genes}}
    bm = cfg["biomarker"]["gene"]
    report["biomarker"] = {"genes": [bm], "missing": [] if bm in present else [bm],
                           "detection_rate": {bm: detection_rate(bm)}}
    return report


def umi_and_gene_counts(adata):
    X = adata.X
    total_umi = np.asarray(X.sum(1)).ravel()
    n_genes = np.asarray((X > 0).sum(1)).ravel()
    return total_umi, n_genes


def mitochondrial_fraction(adata, total_umi):
    """Fraction of a cell's counts coming from MT-prefixed genes. Returns (mt_frac, mt_gene_list)."""
    mt_genes = [g for g in adata.var_names if g.upper().startswith("MT-")]
    if not mt_genes:
        return None, []
    idx = [list(adata.var_names).index(g) for g in mt_genes]
    mt_counts = np.asarray(adata.X[:, idx].sum(1)).ravel()
    mt_frac = mt_counts / np.maximum(total_umi, 1)
    return mt_frac, mt_genes


def doublet_candidates(total_umi, n_genes, umi_multiplier=2.0, genes_multiplier=2.0):
    """Lightweight heuristic: flag cells anomalously high in BOTH total UMI and gene count.

    Not a trained classifier (see module docstring) -- a coarse, auditable proxy for
    "candidate segmentation-merged cells," reported so the *scale* of the concern
    can be judged, not to make a validated per-cell call.
    """
    median_umi, median_genes = np.median(total_umi), np.median(n_genes)
    is_candidate = (total_umi > umi_multiplier * median_umi) & (n_genes > genes_multiplier * median_genes)
    return is_candidate, median_umi, median_genes


def plot_distributions(total_umi, n_genes, mt_frac, area_um2, doublet_flag, out_path,
                        candidate_min_umi, candidate_min_genes):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    axes[0, 0].hist(np.log10(total_umi + 1), bins=100, color="#4c72b0")
    axes[0, 0].axvline(np.log10(candidate_min_umi + 1), color="crimson", linestyle="--",
                        label=f"candidate min_umi = {candidate_min_umi}")
    axes[0, 0].set_xlabel("log10(total UMI + 1) per cell")
    axes[0, 0].set_ylabel("cells")
    axes[0, 0].set_title("UMI per cell")
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].hist(n_genes, bins=100, range=(0, np.percentile(n_genes, 99.5)), color="#55a868")
    axes[0, 1].axvline(candidate_min_genes, color="crimson", linestyle="--",
                        label=f"candidate min_genes = {candidate_min_genes}")
    axes[0, 1].set_xlabel("genes detected per cell")
    axes[0, 1].set_ylabel("cells")
    axes[0, 1].set_title("Genes per cell")
    axes[0, 1].legend(fontsize=8)

    if mt_frac is not None:
        axes[0, 2].hist(mt_frac * 100, bins=100, range=(0, min(100, np.percentile(mt_frac * 100, 99.5))),
                         color="#c44e52")
        axes[0, 2].set_xlabel("% mitochondrial counts")
        axes[0, 2].set_ylabel("cells")
        axes[0, 2].set_title("Mitochondrial fraction (reporting only)")
    else:
        axes[0, 2].text(0.5, 0.5, "no MT- genes\nin this panel", ha="center", va="center")
        axes[0, 2].set_axis_off()

    if area_um2 is not None:
        axes[1, 0].hist(area_um2, bins=100, range=(0, np.percentile(area_um2, 99.5)), color="#8172b2")
        axes[1, 0].set_xlabel("segmented cell area (um^2)")
        axes[1, 0].set_ylabel("cells")
        axes[1, 0].set_title("Cell area (reporting only)")
    else:
        axes[1, 0].text(0.5, 0.5, "no cell_area_um2\nin obs", ha="center", va="center")
        axes[1, 0].set_axis_off()

    n = len(total_umi)
    idx = np.arange(n) if n <= 40000 else np.random.default_rng(0).choice(n, 40000, replace=False)
    axes[1, 1].scatter(total_umi[idx], n_genes[idx], s=2, c=np.where(doublet_flag[idx], "#c44e52", "#4c72b0"),
                        linewidths=0, alpha=0.5)
    axes[1, 1].set_xlabel("total UMI")
    axes[1, 1].set_ylabel("genes detected")
    axes[1, 1].set_title(f"Doublet-candidate heuristic (reporting only): {int(doublet_flag.sum())} cells")

    axes[1, 2].set_axis_off()

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--config", default=str(HERE / "classification_config.json"))
    ap.add_argument("--out", required=True, help="output PNG for the distribution plots")
    ap.add_argument("--candidate-min-umi", type=int, default=None,
                     help="override the reported candidate min_umi (default: current config value)")
    ap.add_argument("--candidate-min-genes", type=int, default=None,
                     help="override the reported candidate min_genes (default: current config value)")
    ap.add_argument("--doublet-umi-multiplier", type=float, default=2.0)
    ap.add_argument("--doublet-genes-multiplier", type=float, default=2.0)
    args = ap.parse_args()

    adata = ad.read_h5ad(args.h5ad)
    cfg = json.loads(Path(args.config).read_text())

    gene_report = check_genes(adata, cfg)
    total_umi, n_genes = umi_and_gene_counts(adata)
    mt_frac, mt_genes = mitochondrial_fraction(adata, total_umi)
    area_um2 = adata.obs["cell_area_um2"].to_numpy() if "cell_area_um2" in adata.obs else None
    doublet_flag, median_umi, median_genes = doublet_candidates(
        total_umi, n_genes, args.doublet_umi_multiplier, args.doublet_genes_multiplier)

    candidate_min_umi = args.candidate_min_umi if args.candidate_min_umi is not None else cfg["qc"]["min_umi"]
    candidate_min_genes = args.candidate_min_genes if args.candidate_min_genes is not None else cfg["qc"]["min_genes"]

    print(f"n_cells = {adata.n_obs}, n_genes_in_panel = {adata.n_vars}")
    print()
    print("Signature / biomarker gene presence and detection rate:")
    all_present = True
    for sig, d in gene_report.items():
        if d["missing"]:
            all_present = False
            print(f"  {sig}: {d['genes']} -> MISSING: {d['missing']}")
        else:
            rates = ", ".join(f"{g}={d['detection_rate'][g]:.1%}" for g in d["genes"])
            weak = [g for g, r in d["detection_rate"].items() if r is not None and r < 0.05]
            flag = f"  <-- WEAK (<5% detection): {weak}" if weak else ""
            print(f"  {sig}: {rates}{flag}")
    print()
    if not all_present:
        print("ABSTAIN: at least one required gene is missing from this panel. "
              "Do not proceed to phenotype metrics until this is resolved (see CLAUDE.md invariants).")

    print("UMI per cell:      min=%d  p1=%.0f  median=%.0f  p99=%.0f  max=%d" % (
        total_umi.min(), np.percentile(total_umi, 1), np.median(total_umi),
        np.percentile(total_umi, 99), total_umi.max()))
    print("Genes per cell:    min=%d  p1=%.0f  median=%.0f  p99=%.0f  max=%d" % (
        n_genes.min(), np.percentile(n_genes, 1), np.median(n_genes),
        np.percentile(n_genes, 99), n_genes.max()))
    print()
    for min_umi in (20, 50, 100, 200, 500):
        frac = float((total_umi >= min_umi).mean())
        print(f"  fraction of cells with UMI >= {min_umi}: {frac:.3f}")
    print()
    for min_genes in (5, 50, 100, 200, 300, 500):
        frac = float((n_genes >= min_genes).mean())
        print(f"  fraction of cells with genes >= {min_genes}: {frac:.3f}")

    print()
    if mt_frac is not None:
        print(f"Mitochondrial fraction ({len(mt_genes)} MT- genes: {mt_genes}) [reporting only, no gate applied]:")
        print("  median=%.2f%%  p95=%.2f%%  p99=%.2f%%  max=%.2f%%" % (
            np.median(mt_frac) * 100, np.percentile(mt_frac, 95) * 100,
            np.percentile(mt_frac, 99) * 100, mt_frac.max() * 100))
    else:
        print("Mitochondrial fraction: no MT- prefixed genes found in this panel; skipped.")

    print()
    if area_um2 is not None:
        print("Cell area (um^2) [reporting only, no gate applied]:")
        print("  p1=%.1f  median=%.1f  p99=%.1f  max=%.1f" % (
            np.percentile(area_um2, 1), np.median(area_um2), np.percentile(area_um2, 99), area_um2.max()))
    else:
        print("Cell area: obs['cell_area_um2'] not present (re-run the loader to populate it).")

    print()
    print(f"Doublet-candidate heuristic [reporting only, no gate applied]: "
          f"UMI > {args.doublet_umi_multiplier}x median ({median_umi:.0f}) "
          f"AND genes > {args.doublet_genes_multiplier}x median ({median_genes:.0f})")
    print(f"  {int(doublet_flag.sum())} of {len(doublet_flag)} cells flagged ({doublet_flag.mean():.2%})")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    plot_distributions(total_umi, n_genes, mt_frac, area_um2, doublet_flag, args.out,
                        candidate_min_umi, candidate_min_genes)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
