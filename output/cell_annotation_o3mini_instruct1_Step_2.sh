#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/cell_annotation_o3mini_instruct1
conda install -y scanpy
cat << 'EOF' > ./output/cell_annotation_o3mini_instruct1/run_qc.py
import scanpy as sc

# Read the raw AnnData object; using the provided raw file
adata = sc.read_h5ad("./output/cell_annotation_o3mini_instruct1/adata_raw.h5ad")

# Annotate mitochondrial genes (assuming gene names starting with MT-)
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calculate quality control metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# Filter cells with high mitochondrial gene percentage (less than 20%), low total UMIs (at least 500), and low number of detected genes (at least 200)
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata = adata[adata.obs.total_counts >= 500, :]
adata = adata[adata.obs.n_genes_by_counts >= 200, :]

# Save the filtered AnnData object
adata.write("./output/cell_annotation_o3mini_instruct1/adata_qc.h5ad")
EOF
python ./output/cell_annotation_o3mini_instruct1/run_qc.py
