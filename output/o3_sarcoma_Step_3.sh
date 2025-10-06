#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_sarcoma
conda install -y scanpy
mkdir -p ./output/o3_sarcoma
cat << 'EOF' > ./output/o3_sarcoma/normalize_hvg_script.py
import scanpy as sc

# Read the QC filtered AnnData object
adata = sc.read_h5ad("./output/o3_sarcoma/qc_filtered.h5ad")

# Normalize total counts per cell to a target sum (default 1e4) and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes (default parameters)
sc.pp.highly_variable_genes(adata)

# Save the normalized data with HVG information
adata.write("./output/o3_sarcoma/normalized_data.h5ad")
EOF
python ./output/o3_sarcoma/normalize_hvg_script.py
