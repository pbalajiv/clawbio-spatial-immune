#!/usr/bin/env bash
# Download the T1 dataset: Visium HD Human Breast Cancer (Fresh Frozen), Ultima sequencing,
# Space Ranger 4.0.1. Source page:
# https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-human-breast-cancer-ff-ultima-4
#
# Only the files this project needs are fetched (not the full batch-download bundle,
# which also includes ~20 GB of binned outputs and Loupe files this spec does not use).
set -euo pipefail

OUT_DIR="${1:-data/raw}"
mkdir -p "$OUT_DIR"

BASE_URL="https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_FF_Human_Breast_Cancer"

# file, expected size (informational only, not enforced)
# Visium_HD_FF_Human_Breast_Cancer_segmented_outputs.tar.gz   2.44 GB
# Visium_HD_FF_Human_Breast_Cancer_spatial.tar.gz             100 MB
# Visium_HD_FF_Human_Breast_Cancer_metrics_summary.csv        2.05 kB

FILES=(
  "Visium_HD_FF_Human_Breast_Cancer_segmented_outputs.tar.gz"
  "Visium_HD_FF_Human_Breast_Cancer_spatial.tar.gz"
  "Visium_HD_FF_Human_Breast_Cancer_metrics_summary.csv"
)

for f in "${FILES[@]}"; do
  echo "Downloading $f ..."
  curl -fL --retry 3 -o "$OUT_DIR/$f" "$BASE_URL/$f"
done

echo "Extracting archives ..."
tar -xzf "$OUT_DIR/Visium_HD_FF_Human_Breast_Cancer_segmented_outputs.tar.gz" -C "$OUT_DIR"
tar -xzf "$OUT_DIR/Visium_HD_FF_Human_Breast_Cancer_spatial.tar.gz" -C "$OUT_DIR"

echo "Writing checksums to data/CHECKSUMS.txt ..."
(
  cd "$OUT_DIR"
  sha256sum "${FILES[@]}" > ../CHECKSUMS.txt
)

echo "Done. See data/DATASET.md for provenance and data/CHECKSUMS.txt for hashes."
