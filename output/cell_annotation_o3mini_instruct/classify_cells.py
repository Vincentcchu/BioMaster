import scanpy as sc
import pandas as pd
import numpy as np

# Load the normalized high variable gene dataset
adata = sc.read_h5ad("./output/cell_annotation_o3mini_instruct/normalized_hvg_dataset.h5ad")

# Load the clustering results and add clustering information if available
clustering = sc.read_h5ad("./output/cell_annotation_o3mini_instruct/clustering_results.h5ad")
if 'leiden' in clustering.obs.columns:
    adata.obs['leiden'] = clustering.obs['leiden']

# Load the marker genes table (assuming it has a column named 'gene')
markers = pd.read_csv("./output/cell_annotation_o3mini_instruct/marker_genes_table.txt", sep="\t")
if 'gene' in markers.columns:
    marker_genes = markers['gene'].tolist()
else:
    marker_genes = []

# Filter marker genes that are present in the dataset
valid_markers = [gene for gene in marker_genes if gene in adata.var_names]

# Compute average expression of the valid marker genes per cell
if valid_markers:
    marker_expr = adata[:, valid_markers].X
    if hasattr(marker_expr, "toarray"):
        marker_expr = marker_expr.toarray()
    avg_expr = np.mean(marker_expr, axis=1)
else:
    avg_expr = np.zeros(adata.n_obs)

# Classify cells: use a simple threshold on the average expression of marker genes (e.g., threshold = 1.0)
threshold = 1.0
adata.obs['malignant'] = ['malignant' if expr > threshold else 'non-malignant' for expr in avg_expr]

# Prepare the final classification result with cell identifiers, classification and average marker expression
classification = adata.obs.copy()
classification['avg_marker_expr'] = avg_expr

# Save the classification results to the output file
classification.to_csv("./output/cell_annotation_o3mini_instruct/malignant_cells_classification.txt", sep="\t")
