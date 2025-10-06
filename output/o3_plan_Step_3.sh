#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_plan
conda install -y scanpy
cat << 'EOF' > ./output/o3_plan/save_malignancy_classification.py
import scanpy as sc

# Load the AnnData object with malignancy scores and binary classification labels
adata = sc.read_h5ad('./output/o3_plan/scored_classified_dataset.h5ad')

# Save the final AnnData object with malignancy annotations
adata.write('./output/o3_plan/malignancy_classification.h5ad')
EOF
python ./output/o3_plan/save_malignancy_classification.py
