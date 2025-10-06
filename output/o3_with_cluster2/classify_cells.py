import scanpy as sc
import pandas as pd
import numpy as np

# Load the AnnData object with clustering annotations
adata = sc.read_h5ad("./output/o3_with_cluster2/clustered_adata.h5ad")

# Load the marker gene table
markers = pd.read_csv("./output/o3_with_cluster2/marker_genes.csv")

# Assuming the marker gene table has a column named gene Define a threshold expression value
threshold = 1.0

# Filter marker genes that are present in the AnnData object
available_markers = [gene for gene in markers["gene"] if gene in adata.var_names]
if len(available_markers) == 0:
    raise ValueError("No marker genes found in the AnnData object!")

# Extract the expression of the available marker genes for all cells
marker_expression = adata[:, available_markers].X

# Compute the average expression of marker genes for each cell without converting the whole sparse matrix to dense
if hasattr(marker_expression, "mean"):
    avg_marker_expression = np.array(marker_expression.mean(axis=1)).flatten()
else:
    avg_marker_expression = marker_expression.mean(axis=1)

# Classify cells as malignant if average marker expression is above or equal to the threshold
cell_classification = ["malignant" if expr >= threshold else "non-malignant" for expr in avg_marker_expression]

# Create a DataFrame with cell barcodes and their classification
df = pd.DataFrame({
    "cell_barcode": adata.obs_names,
    "classification": cell_classification
})

# Save the classification results to a CSV file
df.to_csv("./output/o3_with_cluster2/cell_malignancy_annotation.csv", index=False)
