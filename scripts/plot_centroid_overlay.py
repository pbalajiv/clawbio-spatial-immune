#!/usr/bin/env python3
"""
T2 acceptance check: overlay loader centroids on the hires H&E image to confirm alignment.

Usage:
    python scripts/plot_centroid_overlay.py --h5ad data/tissue.h5ad \
        --segmented-outputs data/raw/segmented_outputs --out data/centroid_overlay.png
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
from PIL import Image, ImageDraw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--segmented-outputs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-points", type=int, default=20000, help="subsample for a readable plot")
    args = ap.parse_args()

    a = ad.read_h5ad(args.h5ad)
    if a.uns.get("spatial_units") != "um":
        raise ValueError("expected obsm['spatial'] in micrometres")

    seg_dir = Path(args.segmented_outputs)
    scalefactors = json.loads((seg_dir / "spatial" / "scalefactors_json.json").read_text())
    um_per_px = scalefactors["microns_per_pixel"]
    hires_scalef = scalefactors["tissue_hires_scalef"]

    xy_um = a.obsm["spatial"]
    xy_fullres_px = xy_um / um_per_px
    xy_hires_px = xy_fullres_px * hires_scalef

    img = Image.open(seg_dir / "spatial" / "tissue_hires_image.png").convert("RGB")
    draw = ImageDraw.Draw(img)

    n = xy_hires_px.shape[0]
    idx = np.arange(n) if n <= args.max_points else np.random.default_rng(0).choice(n, args.max_points, replace=False)
    r = 1.2
    for x, y in xy_hires_px[idx]:
        draw.ellipse([x - r, y - r, x + r, y + r], fill=(255, 0, 0))

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    img.save(args.out)
    print(f"wrote {args.out}: {len(idx)} of {n} centroids overlaid on {img.size} image")


if __name__ == "__main__":
    main()
