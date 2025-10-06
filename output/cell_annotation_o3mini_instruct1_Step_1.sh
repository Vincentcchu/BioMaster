#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct1
mkdir -p ./output/cell_annotation_o3mini_instruct1
conda install -y scanpy
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct1/read_adata.py
import scanpy as sc

# Read the single-cell RNA-seq dataset using Scanpy
adata = sc.read_h5ad("./data/cell_annotation/dataset_debug_restricted.h5ad")

# Save the AnnData object to the specified output file
adata.write("./output/cell_annotation_o3mini_instruct1/adata_raw.h5ad")
EOF
python ./output/cell_annotation_o3mini_instruct1/read_adata.py
