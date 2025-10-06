#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct
conda install -y scanpy
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct/initialize_adat.py
import scanpy as sc

# Read the raw single-cell RNA-seq dataset
adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad')

# Write the initialized AnnData object to the output file
adata.write('./output/cell_annotation_o3mini_instruct/initialized_dataset.h5ad')
EOF
python ./output/cell_annotation_o3mini_instruct/initialize_adat.py
