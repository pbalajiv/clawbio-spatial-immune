# Dataset selection

Goal: a public 10x breast cancer dataset that is (a) cell-level, ideally Space Ranger v4 native segmentation, and (b) HER2+ with documented status. As of 2026-09-25 no candidate below was confirmed to satisfy both. Choose deliberately and record the choice in `data/DATASET.md`.

## Candidates

| Dataset | Cell-level? | HER2 status documented? | Notes |
|---|---|---|---|
| Visium HD, Human Breast Cancer (Fresh Frozen), Ultima sequencing. https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-human-breast-cancer-ff-ultima-4 | Yes, processed with Space Ranger v4.0.1 | No (page says ductal carcinoma in situ, BioIVT) | 10x notes freezing artifacts in parts of the H&E that impair cell segmentation. **Default choice for the build.** |
| Visium HD, Human Breast Cancer (Fresh Frozen). https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-breast-cancer-fresh-frozen | Only after rerunning Space Ranger v4 from FASTQs | Not seen on page | Space Ranger v3.1.1 outputs (binned only). |
| Visium HD, Human Breast Cancer (Fixed Frozen). https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-breast-cancer-fixed-frozen | Only after rerun | Not seen on page | Invasive ductal carcinoma (Avaden). v3.1.1 outputs. |
| Visium HD, Human Breast Cancer, IF (FFPE). https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-human-breast-cancer-ffpe-if | Unlikely: v4 segmentation is documented for H&E images | Not seen on page | v3.1.2 outputs. |
| Visium (standard, 55 µm spots), Human Breast Cancer Block A Section 1. https://www.10xgenomics.com/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0 | **No** (spots span many cells) | **Yes**: IDC, stage IIA, ER+, PR−, HER2+ | Not usable with this spec, which requires cells. |
| Xenium FFPE human breast (Janesick et al.) | Yes (native) | Reported as HER2+ in the literature; verify on the dataset page | Targeted panel (~300 genes), so check every signature gene is present; needs a Xenium loader instead of the Space Ranger one. |

## Recommendation

1. **Build and demo on the Visium HD FF Ultima dataset** (native v4 segmentation, no rerun). Set `sample_subtype` to `"unknown"` in the final answer: its clinical HER2 status is not documented, so HER2+-specific references must be cited with an explicit mismatch/unknown acknowledgment. The pipeline still measures ERBB2 transcript level and pattern, which is exactly what it claims to do.
2. **Scientific caveat specific to DCIS:** in situ lesions are ducts bounded by a myoepithelial layer and basement membrane, so T cells sitting outside ducts may produce an "excluded"-like pattern for anatomical reasons. The immune phenotype literature is largely about invasive tumours. If the section contains invasive regions, consider restricting the analysis to them (annotated manually or by a documented rule) and say so; otherwise state this caveat prominently.
3. **If a documented HER2+ sample is essential**, the options are: rerun Space Ranger v4 with segmentation on the FF or fixed-frozen IDC Visium HD dataset from FASTQs (hours of compute and large downloads; do it before the hackathon, not during), or use the Xenium HER2+ dataset with a Xenium loader.

## Before downloading, confirm on the dataset page ("Output and supplemental files" tab)

- A `segmented_outputs` archive exists (for Space Ranger v4 datasets it is provided alongside `binned_outputs` and `spatial`).
- File sizes, so the download fits your disk and the hackathon timeline.
- The licence (10x datasets are typically CC BY 4.0) and the citation guideline.
