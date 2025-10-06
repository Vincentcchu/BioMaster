#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_with_cluster2
conda install -y scanpy
cat << 'EOF' > ./output/o3_with_cluster2/normalize_hvg.py
import scanpy as sc

# Load the QC-filtered AnnData object
adata = sc.read_h5ad('./output/o3_with_cluster2/filtered_adata.h5ad')

# Normalize the data so that total counts per cell equal 1e4
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data
sc.pp.log1p(adata)

# Identify highly variable genes (HVGs)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Save the normalized AnnData object including HVG annotations
adata.write('./output/o3_with_cluster2/normalized_adata.h5ad')
EOF
python ./output/o3_with_cluster2/normalize_hvg.py
