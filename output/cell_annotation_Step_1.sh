#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation
conda install -y scanpy
mkdir -p ./output/cell_annotation
cat << 'EOF' > ./output/cell_annotation/qc_filter.py
import scanpy as sc

# Read the single-cell RNA-seq dataset
adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_cleaned.h5ad')

# Identify mitochondrial genes assuming they start with 'MT-'
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter out low quality cells: remove cells with total counts less than 500 and with mitochondrial percentage 20% or higher
adata = adata[adata.obs.total_counts >= 500, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

# Filter out low expressed genes: retain genes expressed in at least 3 cells
sc.pp.filter_genes(adata, min_cells=3)

# Save the filtered AnnData object
adata.write('./output/cell_annotation/filtered_dataset.h5ad')
EOF
python ./output/cell_annotation/qc_filter.py
