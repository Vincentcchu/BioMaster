#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_dormancy_with_malignancy2
conda install -y scanpy
cat << 'EOF' > ./output/o3_dormancy_with_malignancy2/extract_malignant.py
import scanpy as sc

# Read the preprocessed AnnData object with quality controlled and normalized data
adata = sc.read_h5ad('./output/o3_dormancy_with_malignancy2/preprocessed_data.h5ad')

# Filter the cells to retain only those labeled as 'malignant'
# Assuming that the malignancy annotations are stored in the adata.obs column named 'malignancy'
malignant_adata = adata[adata.obs['malignancy'] == 'malignant'].copy()

# Write the filtered AnnData object to the output file
malignant_adata.write('./output/o3_dormancy_with_malignancy2/malignant_cells.h5ad')
EOF
python ./output/o3_dormancy_with_malignancy2/extract_malignant.py
