import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set up scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Read the single-cell data
adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad')

print(f'Initial data shape: {adata.shape}')
print(f'Number of cells: {adata.n_obs}')
print(f'Number of genes: {adata.n_vars}')

# Make variable names unique (in case of duplicates)
adata.var_names_unique()

# Calculate QC metrics
# Calculate mitochondrial gene percentage
adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Calculate mitochondrial gene percentage manually if needed
adata.obs['pct_counts_mt'] = adata[:, adata.var['mt']].X.sum(1).A1 / adata.obs['total_counts'] * 100

# Save QC metrics before filtering
print(f'\nQC metrics before filtering:')
print(f'Mean genes per cell: {adata.obs["n_genes_by_counts"].mean():.2f}')
print(f'Mean UMIs per cell: {adata.obs["total_counts"].mean():.2f}')
print(f'Mean mitochondrial percentage: {adata.obs["pct_counts_mt"].mean():.2f}%')

# Filter cells with low gene counts (< 200 genes)
sc.pp.filter_cells(adata, min_genes=200)

# Filter genes expressed in very few cells (< 3 cells)
sc.pp.filter_genes(adata, min_cells=3)

# Filter cells with high mitochondrial gene percentage (> 20%)
adata = adata[adata.obs['pct_counts_mt'] < 20, :]

# Filter cells with low total UMI counts (< 500)
adata = adata[adata.obs['total_counts'] >= 500, :]

# Save QC metrics after filtering
print(f'\nQC metrics after filtering:')
print(f'Final data shape: {adata.shape}')
print(f'Number of cells after filtering: {adata.n_obs}')
print(f'Number of genes after filtering: {adata.n_vars}')
print(f'Mean genes per cell: {adata.obs["n_genes_by_counts"].mean():.2f}')
print(f'Mean UMIs per cell: {adata.obs["total_counts"].mean():.2f}')
print(f'Mean mitochondrial percentage: {adata.obs["pct_counts_mt"].mean():.2f}%')

# Save raw counts for later use
adata.raw = adata

# Normalize to 10,000 reads per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

print(f'\nNumber of highly variable genes: {sum(adata.var["highly_variable"])}')

# Save the preprocessed data
adata.write('./output/cell_annotation_sonnet4/preprocessed_data.h5ad')

print('\nData preprocessing completed successfully!')
print('Preprocessed data saved to: ./output/cell_annotation_sonnet4/preprocessed_data.h5ad')
