#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_sarcoma
conda install -y scanpy
mkdir -p ./output/o3_sarcoma
cat << 'EOF' > ./output/o3_sarcoma/read_initial_data.py
import scanpy as sc

# Read the sarcoma single-cell RNA-seq dataset in h5ad format
adata = sc.read_h5ad('/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/sarcoma_dataset_restricted.h5ad')

# Write the initialized AnnData object to output
adata.write('./output/o3_sarcoma/initial_data.h5ad')
EOF
python ./output/o3_sarcoma/read_initial_data.py
