#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_full_dataset
conda install -y scanpy
mkdir -p ./output/o3_full_dataset
cat << 'EOF' > ./output/o3_full_dataset/scanpy_init.py
import scanpy as sc

# Load the two single-cell RNA-seq datasets
adata1 = sc.read_h5ad('/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/dataset_restricted.h5ad')
adata2 = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad')

# Concatenate the two AnnData objects along the observations (cells)
adata = adata1.concatenate(adata2)

# Write the initialized AnnData object to output
adata.write('./output/o3_full_dataset/initialized_dataset.h5ad')
EOF
python ./output/o3_full_dataset/scanpy_init.py
