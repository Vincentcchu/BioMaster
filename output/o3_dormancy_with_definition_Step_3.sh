#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/o3_dormancy_with_definition
conda install -y scanpy
cat << 'EOF' > ./output/o3_dormancy_with_definition/compute_cell_cycle_scores.py
import scanpy as sc

# Load the malignant cells dataset
adata = sc.read_h5ad("./output/o3_dormancy_with_definition/malignant_cells.h5ad")

# Define cell cycle gene sets for scoring (example gene lists)
s_genes = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2",
    "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1"
]
g2m_genes = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
    "TACC3", "CDC20", "TTK", "CDC25C", "KIF11", "ANP32E", "LMNB1",
    "KIF20B", "DLGAP5", "KIF2C", "RRM2", "HMGB1", "CCNA2", "CENPE", "CENPA"
]

# Compute cell cycle scores using Scanpy's built-in function
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Save the updated AnnData object with cell cycle scores
adata.write("./output/o3_dormancy_with_definition/malignant_cells_cell_cycle_scores.h5ad")
EOF
python ./output/o3_dormancy_with_definition/compute_cell_cycle_scores.py
