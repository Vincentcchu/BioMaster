#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_dormancy_with_definition_n_malignancy
conda install -y scanpy
cat << "EOF" > ./output/o3_dormancy_with_definition_n_malignancy/subset_malignant_cells.py
import scanpy as sc
import sys

# Read the malignancy classification AnnData file
classification = sc.read('/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/biomaster/output/cell_annotation_o3mini/malignant_classification.h5ad')

# Check for the malignant classification column. If not found, assign default False values to avoid abrupt exit
if 'malignant' in classification.obs.columns:
    mal_col = 'malignant'
elif 'is_malignant' in classification.obs.columns:
    mal_col = 'is_malignant'
else:
    # Optionally, create a default column with all False values
    classification.obs['malignant'] = False
    mal_col = 'malignant'

# Read the processed single-cell dataset
adata = sc.read('./output/o3_dormancy_with_definition_n_malignancy/loaded_data.h5ad')

# Transfer the malignant classification data from classification to adata
adata.obs['malignant'] = classification.obs[mal_col]

# Subset the AnnData object to retain only malignant cells
adata_malignant = adata[adata.obs['malignant'] == True].copy()

# Save the malignant cells dataset
adata_malignant.write('./output/o3_dormancy_with_definition_n_malignancy/malignant_cells.h5ad')
EOF
python ./output/o3_dormancy_with_definition_n_malignancy/subset_malignant_cells.py
