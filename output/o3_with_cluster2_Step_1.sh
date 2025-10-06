#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_with_cluster2
mkdir -p ./output/o3_with_cluster2
conda install -y scanpy
cat << 'EOF' > ./output/o3_with_cluster2/load_adata.py
import scanpy as sc

# Read the single-cell RNA-seq dataset from the provided .h5ad file
adata = sc.read_h5ad("/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/dataset_restricted_with_cluster.h5ad")

# Save the AnnData object to the output file
adata.write("./output/o3_with_cluster2/initial_adata.h5ad")
EOF
python ./output/o3_with_cluster2/load_adata.py
