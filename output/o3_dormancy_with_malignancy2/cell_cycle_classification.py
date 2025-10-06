import scanpy as sc
import pandas as pd

# Load the AnnData object containing malignant cells with computed cell cycle scores
adata = sc.read_h5ad("./output/o3_dormancy_with_malignancy2/malignant_cells_with_scores.h5ad")

# Calculate a combined cell cycle score assuming scores for S phase and G2M phase are present
# Adjust the threshold values based on your data distribution
if 'S_score' in adata.obs and 'G2M_score' in adata.obs:
    adata.obs['cell_cycle_score'] = adata.obs['S_score'] + adata.obs['G2M_score']
else:
    # Fallback: assume a pre-computed 'cell_cycle_score' exists
    adata.obs['cell_cycle_score'] = adata.obs['cell_cycle_score']

# Define classification based on thresholds (these thresholds are examples and may need tuning)
# For example, if cell_cycle_score >= 1.0, classify as 'proliferative'; if <= 0.5, classify as 'quiescent'; otherwise 'slow-cycling'

def classify_cell(score):
    if score >= 1.0:
        return 'proliferative'
    elif score <= 0.5:
        return 'quiescent'
    else:
        return 'slow-cycling'

adata.obs['cell_cycle_state'] = adata.obs['cell_cycle_score'].apply(classify_cell)

# Extract the cell barcodes and their classification
classification_df = adata.obs[['cell_cycle_state']].copy()
classification_df.index.name = 'cell_barcode'
classification_df.reset_index(inplace=True)

# Save the classification results and cell barcodes to a CSV file
classification_df.to_csv("./output/o3_dormancy_with_malignancy2/cell_cycle_classification.csv", index=False)
