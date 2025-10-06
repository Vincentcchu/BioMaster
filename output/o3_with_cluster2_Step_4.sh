#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_with_cluster2
conda install -y scanpy matplotlib
cat << 'EOF' > ./output/o3_with_cluster2/cluster_analysis.py
import scanpy as sc
import matplotlib.pyplot as plt

# Load the normalized AnnData object with highly variable genes
adata = sc.read_h5ad("./output/o3_with_cluster2/normalized_adata.h5ad")

# Perform PCA for dimensionality reduction
sc.tl.pca(adata, svd_solver="arpack")

# Compute the neighborhood graph
sc.pp.neighbors(adata)

# Apply the Leiden algorithm for clustering
sc.tl.leiden(adata)

# Generate UMAP visualization
sc.tl.umap(adata)

# Save the AnnData object with clustering annotations
adata.write("./output/o3_with_cluster2/clustered_adata.h5ad")

# Plot UMAP with clusters and save the figure
sc.pl.umap(adata, color="leiden", show=False)
plt.savefig("./output/o3_with_cluster2/umap.png")
EOF
python ./output/o3_with_cluster2/cluster_analysis.py
