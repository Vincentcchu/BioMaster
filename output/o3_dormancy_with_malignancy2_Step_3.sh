#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_dormancy_with_malignancy2
conda install -y scanpy
mkdir -p ./output/o3_dormancy_with_malignancy2
cat << 'EOF' > ./output/o3_dormancy_with_malignancy2/compute_cell_cycle_scores.py
import scanpy as sc

# Define marker gene lists for cell cycle phases
s_genes = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "RFC2", "NASP", "GMNN",
    "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
    "RAD51AP1", "RRM2", "CDC45", "CDC6", "EXO1"
]

g2m_genes = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CDC20", "TTK", "CDC25C", "KIF11",
    "RRM2", "PBK", "NCAPG", "DLGAP5", "KIF20B", "HMMR", "CCNB2",
    "CKS1B", "MKI67", "TMPO", "CENPF"
]

# Read the malignant cells AnnData object
adata = sc.read_h5ad('./output/o3_dormancy_with_malignancy2/malignant_cells.h5ad')

# Compute cell cycle scores
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Save the AnnData object with the computed scores
adata.write('./output/o3_dormancy_with_malignancy2/malignant_cells_with_scores.h5ad')
EOF
python ./output/o3_dormancy_with_malignancy2/compute_cell_cycle_scores.py
