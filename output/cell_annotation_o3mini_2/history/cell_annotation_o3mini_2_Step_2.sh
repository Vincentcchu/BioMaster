#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_2
conda install -y scanpy
cat << "EOF" > ./output/cell_annotation_2/run_qc_filtering.py
import scanpy as sc

# Read the initialized AnnData object containing raw gene expression data
adata = sc.read_h5ad('./output/cell_annotation_2/initialized_dataset.h5ad')

# Identify mitochondrial genes (assuming gene names starting with MT-)
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate QC metrics including percentage of mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter cells: mitochondrial percentage < 20 percent, total UMIs >= 500, and number of detected genes >= 200
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata = adata[adata.obs.total_counts >= 500, :]
adata = adata[adata.obs.n_genes_by_counts >= 200, :]

# Save the filtered AnnData object
adata.write('./output/cell_annotation_2/filtered_dataset.h5ad')
EOF
python ./output/cell_annotation_2/run_qc_filtering.py
