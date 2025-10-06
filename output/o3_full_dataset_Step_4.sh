#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_full_dataset
conda install -y scanpy matplotlib
cat << 'EOF' > ./output/o3_full_dataset/scanpy_clustering.py
import scanpy as sc
import matplotlib.pyplot as plt

# Read the initialized AnnData object
adata = sc.read_h5ad("./output/o3_full_dataset/initialized_dataset.h5ad")

# Perform PCA for dimensionality reduction
sc.tl.pca(adata, svd_solver='arpack')

# Construct the KNN graph using the first 40 principal components
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Perform Leiden clustering
sc.tl.leiden(adata)

# Compute UMAP for visualization
sc.tl.umap(adata)

# Plot UMAP with Leiden cluster labels
sc.pl.umap(adata, color='leiden', show=False)
plt.savefig("./output/o3_full_dataset/umap_clusters.png")

# Save the AnnData object with computed PCA, neighborhood graph, Leiden clusters, and UMAP coordinates
adata.write("./output/o3_full_dataset/clusters_dataset.h5ad")
EOF
python ./output/o3_full_dataset/scanpy_clustering.py
