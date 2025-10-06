#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_dormancy_with_malignancy
conda install -y scanpy
cat <<EOF > ./output/o3_dormancy_with_malignancy/cell_cycle_scoring.py
import scanpy as sc
import numpy as np

# Read the malignant cells AnnData object
adata = sc.read("./output/o3_dormancy_with_malignancy/adata_malignant.h5ad")

# Convert gene names in adata to uppercase for consistency
data = adata.copy()
data.var_names = data.var_names.str.upper()
adata = data

# Define cell cycle marker genes (using uppercase for consistency)
s_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA2", "MCM7", "RFC2", "RPA3", "RPA1", "RFC4", "RFC5", "POLD3", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "PHF6"]
g2m_genes = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "MKI67"]

# Filter marker genes ensuring they exist in adata.var_names
s_genes_filtered = [gene for gene in s_genes if gene in adata.var_names]
g2m_genes_filtered = [gene for gene in g2m_genes if gene in adata.var_names]

if not s_genes_filtered or not g2m_genes_filtered:
    raise ValueError("Filtered cell cycle marker genes list is empty. Check gene names in adata.var_names")

# Perform cell cycle scoring using filtered marker genes
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_filtered, g2m_genes=g2m_genes_filtered)

# Classify malignant cells based on cell cycle scores

def classify_cell(row):
    if row["S_score"] >= 0.5 or row["G2M_score"] >= 0.5:
        return "proliferative"
    elif row["S_score"] <= 0.2 and row["G2M_score"] <= 0.2:
        return "quiescent"
    else:
        return "slow-cycling"

adata.obs["cell_cycle_class"] = adata.obs.apply(lambda row: classify_cell(row), axis=1)

# Save the updated AnnData object with cell cycle classification
adata.write("./output/o3_dormancy_with_malignancy/adata_classified.h5ad")
EOF
python ./output/o3_dormancy_with_malignancy/cell_cycle_scoring.py
