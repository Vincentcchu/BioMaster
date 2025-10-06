#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct
conda install -y scanpy
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct/normalize_hvg.py
import scanpy as sc

# Read the QC-filtered AnnData object
adata = sc.read("./output/cell_annotation_o3mini_instruct/qc_filtered_dataset.h5ad")

# Normalize the total counts per cell to 1e4
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data
sc.pp.log1p(adata)

# Identify highly variable genes (HVGs) using default parameters
sc.pp.highly_variable_genes(adata)

# Subset the AnnData object to keep only highly variable genes
adata = adata[:, adata.var.highly_variable]

# Save the normalized AnnData object with HVGs
adata.write("./output/cell_annotation_o3mini_instruct/normalized_hvg_dataset.h5ad")
EOF
python ./output/cell_annotation_o3mini_instruct/normalize_hvg.py
