#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct
conda install -y scanpy
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct/qc_filter.py
import scanpy as sc

# Read the initialized AnnData object containing raw gene expression matrix and metadata
adata = sc.read_h5ad('./output/cell_annotation_o3mini_instruct/initialized_dataset.h5ad')

# Identify mitochondrial genes by checking if gene names start with 'MT-' and save as a new column in adata.var
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate QC metrics using the newly created 'mt' indicator
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# Filter cells: mitochondrial gene percentage < 20%, total UMIs >= 500, and number of detected genes >= 200
adata = adata[(adata.obs.pct_counts_mt < 20) & (adata.obs.total_counts >= 500) & (adata.obs.n_genes_by_counts >= 200)]

# Save the quality-controlled AnnData object
adata.write('./output/cell_annotation_o3mini_instruct/qc_filtered_dataset.h5ad')
EOF
python ./output/cell_annotation_o3mini_instruct/qc_filter.py
