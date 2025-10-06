#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_dormancy_with_malignancy
mkdir -p ./output/o3_dormancy_with_malignancy
conda install -y scanpy
echo -e "import scanpy as sc\nadata = sc.read_h5ad('/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/dataset_restricted_with_labels.h5ad')\nadata.write('./output/o3_dormancy_with_malignancy/adata_initial.h5ad')" > ./output/o3_dormancy_with_malignancy/read_adata.py
python ./output/o3_dormancy_with_malignancy/read_adata.py
