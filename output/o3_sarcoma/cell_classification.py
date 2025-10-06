import scanpy as sc
import numpy as np
import os

# Load the AnnData object that contains clustering labels and marker gene information
adata = sc.read_h5ad("./output/o3_sarcoma/clustering_results.h5ad")

# Load marker gene list from file (assuming one gene per line)
marker_file = "./output/o3_sarcoma/marker_genes.txt"
if not os.path.exists(marker_file):
    raise FileNotFoundError("Marker gene file not found: " + marker_file)
with open(marker_file, "r") as f:
    marker_genes = [line.strip() for line in f if line.strip() != ""]

# Filter marker genes that exist in the dataset
existing_markers = [gene for gene in marker_genes if gene in adata.var_names]

if len(existing_markers) == 0:
    print("Warning: None of the marker genes were found in the dataset, classifying all cells as non-malignant by default.")
    cell_labels = np.array(["non-malignant"] * adata.n_obs)
else:
    # Compute the average expression of the marker genes for each cell
    marker_expr = adata[:, existing_markers].X
    if hasattr(marker_expr, 'toarray'):
        marker_expr = marker_expr.toarray()
    avg_marker_expr = np.mean(marker_expr, axis=1)
    threshold = np.median(avg_marker_expr)
    cell_labels = np.where(avg_marker_expr > threshold, "malignant", "non-malignant")

# Annotate the AnnData object with the new classification
adata.obs['cell_classification'] = cell_labels

# Save the fully annotated AnnData object
adata.write("./output/o3_sarcoma/cell_classification.h5ad")

print("Cell classification completed and saved to ./output/o3_sarcoma/cell_classification.h5ad")
