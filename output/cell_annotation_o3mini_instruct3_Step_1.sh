#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct3
conda install -y scanpy
mkdir -p ./output/cell_annotation_o3mini_instruct3
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct3/initialize_data.py
import scanpy as sc

# Read the single-cell RNA-seq dataset and initialize the AnnData object
adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad')

# Save the initialized AnnData object for downstream analyses
adata.write('./output/cell_annotation_o3mini_instruct3/initialized_data.h5ad')
EOF
python ./output/cell_annotation_o3mini_instruct3/initialize_data.py
