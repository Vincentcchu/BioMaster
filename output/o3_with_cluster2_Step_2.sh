#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_with_cluster2
conda install -y scanpy
cat << 'EOF' > ./output/o3_with_cluster2/qc_filtering.py
import scanpy as sc

# Read the initial AnnData object
adata = sc.read_h5ad("./output/o3_with_cluster2/initial_adata.h5ad")

# Identify mitochondrial genes assuming they start with "MT-"
adata.var['mt'] = adata.var_names.str.startswith("MT-")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter cells: keep cells with less than 20%% mitochondrial counts, at least 550 total UMIs, and at least 200 detected genes
adata = adata[adata.obs.pct_counts_mt < 20, :]
sc.pp.filter_cells(adata, min_counts=550)
sc.pp.filter_cells(adata, min_genes=200)

# Save the filtered AnnData object
adata.write("./output/o3_with_cluster2/filtered_adata.h5ad")
EOF
python ./output/o3_with_cluster2/qc_filtering.py
