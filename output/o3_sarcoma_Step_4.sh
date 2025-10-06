#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_sarcoma
conda install -y scanpy
cat << 'EOF' > ./output/o3_sarcoma/dimension_reduction_leiden_clustering.py
import scanpy as sc

# Read the initialized AnnData object
adata = sc.read_h5ad('./output/o3_sarcoma/initialized_data.h5ad')

# Basic preprocessing: normalization, log transformation and HVG selection
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
adata = adata[:, adata.var.highly_variable]

# Perform PCA for dimensionality reduction
sc.tl.pca(adata, svd_solver='arpack')

# Compute the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Perform Leiden clustering
sc.tl.leiden(adata)

# Compute UMAP for visualization
sc.tl.umap(adata)

# Save the AnnData object with clustering labels and dimensionality reduction results
adata.write('./output/o3_sarcoma/clustering_results.h5ad')
EOF
python ./output/o3_sarcoma/dimension_reduction_leiden_clustering.py
