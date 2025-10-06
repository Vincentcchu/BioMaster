#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation
conda install -y scanpy
mkdir -p ./output/cell_annotation
cat << 'EOF' > ./output/cell_annotation/clustering_analysis.py
import scanpy as sc

# Read the normalized single-cell RNA-seq dataset
adata = sc.read_h5ad("./output/cell_annotation/normalized_dataset.h5ad")

# Perform PCA for dimensionality reduction
sc.tl.pca(adata, svd_solver='arpack')

# Construct the nearest-neighbors graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Apply Leiden clustering algorithm
sc.tl.leiden(adata)

# Generate UMAP visualization
sc.tl.umap(adata)

# Save the annotated AnnData object with clustering results
adata.write("./output/cell_annotation/clusters_annotated.h5ad")
EOF
