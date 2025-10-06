#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_with_cluster
conda install -y scanpy
cat << 'EOF' > ./output/o3_with_cluster/cluster_cells.py
import scanpy as sc

# Read the preprocessed and normalized single-cell dataset
adata = sc.read_h5ad("./output/o3_with_cluster/preprocessed_dataset.h5ad")

# Compute the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Perform clustering using the Leiden algorithm
sc.tl.leiden(adata, resolution=1.0)

# Compute UMAP for visualization
sc.tl.umap(adata)

# Save the clustered dataset with cell clusters annotated
adata.write("./output/o3_with_cluster/clustered_dataset.h5ad")
EOF
python ./output/o3_with_cluster/cluster_cells.py
