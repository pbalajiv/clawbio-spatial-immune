# Spatial Immune Phenotype Report

**Input**: `data\tissue.h5ad`
**spec_version**: 1.0.0  **seed**: 0

## QC
- Cells in: 38993, passed: 35576, dropped: 3417 (91.2% pass fraction)
- Drop reasons: {'low_umi': 1455, 'low_genes': 0, 'spatial_outlier_low_libsize': 786, 'spatial_outlier_low_genes': 346, 'spatial_outlier_high_mito': 830}

## Immune phenotype
**Label: excluded**
- qc_ok, tumor signature trusted
- t_cell_fraction 0.0970 >= 0.01, t_cell signature trusted
- infiltration_ratio 0.049 <= 0.2 and nhood_z -62.99 <= -2.0 (p=0.0010)

## Biomarker pattern
**Gene: ERBB2. Label: homogeneous**
- morans_i_excess -0.017 < 0.05

## Warnings
- spillover_suspected: T-like cells touching tumour are double-positive 29.65x as often as distant T-like cells
- depth_autocorrelation: Moran's I of log total counts = 0.420

## Figures
- `figures/cell_type_map.png`
- `figures/tumor_region_t_cells.png`
- `figures/biomarker_map.png`
- `figures/stats_panel.png`
