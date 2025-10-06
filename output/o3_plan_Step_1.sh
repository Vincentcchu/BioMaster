#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_plan
conda install -y scanpy
mkdir -p ./output/o3_plan
cat << 'EOF' > ./output/o3_plan/qc_scanpy.py
import scanpy as sc
import numpy as np

# Load the integrated single-cell RNA-seq dataset
adata = sc.read_h5ad("/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/dataset_restricted.h5ad")

# Identify mitochondrial genes (assuming genes starting with "MT-")
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calculate QC metrics including percent mitochondrial counts
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# Filter cells based on QC metrics:
# Retain cells with total counts (UMIs) >= 500, detected genes >= 200, and mitochondrial %% < 20
adata = adata[adata.obs["total_counts"] >= 500, :]
adata = adata[adata.obs["n_genes_by_counts"] >= 200, :]
adata = adata[adata.obs["pct_counts_mt"] < 20, :]

# Define malignancy-related gene signatures for subsequent scoring
malignancy_up_signature = ["MAL_UP1", "MAL_UP2", "MAL_UP3"]
malignancy_down_signature = ["MAL_DOWN1", "MAL_DOWN2", "MAL_DOWN3"]

# Optionally, store the signatures in adata.uns for future use
adata.uns["malignancy_up_signature"] = malignancy_up_signature
adata.uns["malignancy_down_signature"] = malignancy_down_signature

# Save the processed AnnData object after QC and cell type annotation standardization
adata.write("./output/o3_plan/filtered_dataset.h5ad")
print("QC and filtering complete. Processed dataset saved to ./output/o3_plan/filtered_dataset.h5ad")
EOF
python ./output/o3_plan/qc_scanpy.py
