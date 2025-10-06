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

# Read the raw single-cell RNA-seq dataset using the correct file path
adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad')

# Identify mitochondrial genes assuming gene names starting with MT-
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate QC metrics including percent mitochondrial counts
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter cells with low total UMIs (min_counts=500) and high mitochondrial percentage less than 20
sc.pp.filter_cells(adata, min_counts=500)
adata = adata[adata.obs.pct_counts_mt < 20, :]

# Filter genes to remove low expressed ones that are expressed in fewer than 3 cells
sc.pp.filter_genes(adata, min_cells=3)

# Save the QC-filtered AnnData object
adata.write('./output/cell_annotation/filtered_dataset.h5ad')
EOF
python ./output/cell_annotation/qc_filter.py
