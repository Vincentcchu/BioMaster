import scanpy as sc
import numpy as np

# Read the normalized single-cell RNA-seq dataset
adata = sc.read_h5ad("./output/cell_annotation_o3mini_instruct2/normalized_dataset.h5ad")

# Load known malignant marker genes
with open("./output/cell_annotation_o3mini_instruct2/marker_genes.txt", "r") as f:
    markers = [line.strip() for line in f if line.strip()]

# Filter the marker genes to those present in the dataset
markers_present = [gene for gene in markers if gene in adata.var_names]

if not markers_present:
    print("Warning: None of the marker genes were found in the dataset. All cells will be labeled as non-malignant.")
    adata.obs["malignant"] = "non-malignant"
else:
    # For each cell, classify as malignant if any of the marker genes are expressed above 0
    expr = adata[:, markers_present].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    malignant_labels = np.where(expr.sum(axis=1) > 0, "malignant", "non-malignant")
    adata.obs["malignant"] = malignant_labels

# Save the annotated AnnData object
adata.write("./output/cell_annotation_o3mini_instruct2/cell_annotation.h5ad")
