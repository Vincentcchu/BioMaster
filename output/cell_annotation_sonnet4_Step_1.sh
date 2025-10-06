#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_sonnet4
conda install -c conda-forge -c bioconda scanpy pandas numpy -y
python -c "
import scanpy as sc
import pandas as pd
import numpy as np

# Load the dataset
adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad')

# Basic dataset information
print('Dataset Overview:')
print(f'Number of cells: {adata.n_obs}')
print(f'Number of genes: {adata.n_vars}')

# Check data structure
print(f'\nData matrix shape: {adata.X.shape}')
print(f'Data type: {type(adata.X)}')

# Cell metadata
print(f'\nCell metadata columns: {list(adata.obs.columns)}')
if len(adata.obs.columns) > 0:
    print('Cell metadata summary:')
    print(adata.obs.describe())

# Gene metadata
print(f'\nGene metadata columns: {list(adata.var.columns)}')
if len(adata.var.columns) > 0:
    print('Gene metadata summary:')
    print(adata.var.describe())

# Calculate basic QC metrics
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Check for mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Quality control metrics summary
print('\nQuality Control Metrics:')
print(f'Total UMIs per cell - Mean: {adata.obs.total_counts.mean():.2f}, Median: {adata.obs.total_counts.median():.2f}')
print(f'Genes per cell - Mean: {adata.obs.n_genes_by_counts.mean():.2f}, Median: {adata.obs.n_genes_by_counts.median():.2f}')
print(f'Mitochondrial gene percentage - Mean: {adata.obs.pct_counts_mt.mean():.2f}%, Median: {adata.obs.pct_counts_mt.median():.2f}%')

# Save overview to file
with open('./output/cell_annotation_sonnet4/data_overview.txt', 'w') as f:
    f.write('Single-Cell RNA-seq Dataset Overview\n')
    f.write('=====================================\n\n')
    f.write(f'Number of cells: {adata.n_obs}\n')
    f.write(f'Number of genes: {adata.n_vars}\n')
    f.write(f'Data matrix shape: {adata.X.shape}\n')
    f.write(f'Data type: {type(adata.X)}\n\n')
    f.write(f'Cell metadata columns: {list(adata.obs.columns)}\n')
    f.write(f'Gene metadata columns: {list(adata.var.columns)}\n\n')
    f.write('Quality Control Metrics:\n')
    f.write(f'Total UMIs per cell - Mean: {adata.obs.total_counts.mean():.2f}, Median: {adata.obs.total_counts.median():.2f}\n')
    f.write(f'Genes per cell - Mean: {adata.obs.n_genes_by_counts.mean():.2f}, Median: {adata.obs.n_genes_by_counts.median():.2f}\n')
    f.write(f'Mitochondrial gene percentage - Mean: {adata.obs.pct_counts_mt.mean():.2f}%, Median: {adata.obs.pct_counts_mt.median():.2f}%\n')

print('\nData overview saved to ./output/cell_annotation_sonnet4/data_overview.txt')
"
