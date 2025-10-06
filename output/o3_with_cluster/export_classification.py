import scanpy as sc
import pandas as pd

# Read the annotated single-cell dataset with malignant/non-malignant labels
adata = sc.read_h5ad('./output/o3_with_cluster/cell_classification_results.h5ad')

# Extract cell barcodes (index) and their classification (assuming the classification is stored in the 'classification' column)
df = adata.obs[['classification']].reset_index().rename(columns={'index': 'cell_barcode'})

# Export the results to a CSV report
df.to_csv('./output/o3_with_cluster/final_cell_classification_report.csv', index=False)
