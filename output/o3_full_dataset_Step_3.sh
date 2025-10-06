#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_full_dataset
conda install -y scanpy
mkdir -p ./output/o3_full_dataset
cat << 'EOF' > ./output/o3_full_dataset/normalize_hvg.py
import scanpy as sc

# Read the filtered AnnData object
adata = sc.read_h5ad("./output/o3_full_dataset/filtered_dataset.h5ad")

# Normalize total counts per cell to a target sum (default: 1e4)
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data
sc.pp.log1p(adata)

# Identify highly variable genes (HVGs) using the Seurat method and selecting top 2000 genes
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# Save the normalized dataset annotated with HVGs
adata.write("./output/o3_full_dataset/normalized_dataset.h5ad")
EOF
python ./output/o3_full_dataset/normalize_hvg.py
