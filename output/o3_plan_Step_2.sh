#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_plan
conda install -y scanpy
cat << "EOF" > ./output/o3_plan/signature_scoring_classification.py
import scanpy as sc
import numpy as np

# Load the processed AnnData object
adata = sc.read_h5ad("./output/o3_plan/filtered_dataset.h5ad")

# Define updated gene lists for scoring (ensure these genes are present in adata.var_names)
up_genes_all = ["ACTB", "GAPDH", "MYC"]
down_genes_all = ["MT-CO1", "MT-CYB", "MT-ND1"]

# Filter gene lists to include only genes present in the dataset
all_genes = set(adata.var_names)
up_genes = [gene for gene in up_genes_all if gene in all_genes]
down_genes = [gene for gene in down_genes_all if gene in all_genes]

# Check that the filtered gene lists are not empty
if not up_genes or not down_genes:
    raise ValueError("No valid genes were passed for scoring. Please update the gene lists to match adata.var_names.")

# Compute score for upregulated genes
sc.tl.score_genes(adata, gene_list=up_genes, score_name="up_score", use_raw=False)

# Compute score for downregulated genes
sc.tl.score_genes(adata, gene_list=down_genes, score_name="down_score", use_raw=False)

# Calculate malignancy score as the difference between up_score and down_score
adata.obs["malignancy_score"] = adata.obs["up_score"] - adata.obs["down_score"]

# Perform binary classification using a threshold (here 0)
threshold = 0
adata.obs["malignant"] = (adata.obs["malignancy_score"] > threshold).astype(str)
adata.obs["malignant"] = adata.obs["malignant"].replace({"True": "malignant", "False": "non_malignant"})

# Save the AnnData object with malignancy scores and binary classification labels
adata.write("./output/o3_plan/scored_classified_dataset.h5ad")
EOF
python ./output/o3_plan/signature_scoring_classification.py
