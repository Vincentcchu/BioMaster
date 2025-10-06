import scanpy as sc
import sys

# Read the malignant cells dataset
adata = sc.read_h5ad("./output/o3_dormancy_with_definition_n_malignancy/malignant_cells.h5ad")

# Convert gene names in the dataset to uppercase to ensure matching
adata.var_names = adata.var_names.str.upper()

# Define cell cycle gene sets for S phase and G2/M phase and convert them to uppercase
s_genes = [gene.upper() for gene in ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CHK1", "RFC2", "RPA3", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51AP1"]]

g2m_genes = [gene.upper() for gene in ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "CCNB1", "CDC20", "TTK", "CDC2"]]

# Filter gene lists to include only genes present in the dataset
s_genes = [gene for gene in s_genes if gene in adata.var_names]
g2m_genes = [gene for gene in g2m_genes if gene in adata.var_names]

# Exit if no overlapping gene is found to avoid KeyError
if not s_genes or not g2m_genes:
    sys.exit("No overlapping cell cycle genes found in dataset")

# Compute cell cycle scores with use_raw set to False to avoid gene lookup errors in raw data
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, use_raw=False)

# Save the AnnData object with computed cell cycle scores
adata.write("./output/o3_dormancy_with_definition_n_malignancy/cell_cycle_scored_malignant.h5ad")
