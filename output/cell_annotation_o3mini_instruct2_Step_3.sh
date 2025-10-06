#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct2
conda install -y scanpy
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct2/run_scanpy_clustering.py
import scanpy as sc

# Read the filtered single-cell dataset
adata = sc.read_h5ad("./output/cell_annotation_o3mini_instruct2/filtered_dataset.h5ad")

# Perform PCA on the dataset
sc.tl.pca(adata, svd_solver="arpack")

# Compute the neighborhood graph of cells
sc.pp.neighbors(adata)

# Compute UMAP embedding for visualization
sc.tl.umap(adata)

# Perform clustering using the Leiden algorithm
sc.tl.leiden(adata)

# Save the dimension-reduced dataset with clustering labels
adata.write("./output/cell_annotation_o3mini_instruct2/dim_reduced_dataset.h5ad")
EOF
python ./output/cell_annotation_o3mini_instruct2/run_scanpy_clustering.py
