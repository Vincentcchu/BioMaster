import scanpy as sc
import pandas as pd
import numpy as np

# Load the preprocessed single-cell dataset with clustering information
adata = sc.read_h5ad('./output/o3_with_cluster/preprocessed_dataset.h5ad')

# Load the marker genes table assuming the file has a header and a column named gene containing marker gene names
markers_df = pd.read_csv('./output/o3_with_cluster/marker_genes.txt', sep='\t', header=0)
if 'gene' in markers_df.columns:
    marker_genes = markers_df['gene'].tolist()
else:
    marker_genes = markers_df.iloc[:, 0].tolist()

# Filter marker genes to only those that exist in the dataset
available_markers = [gene for gene in marker_genes if gene in adata.var_names]
if len(available_markers) == 0:
    raise ValueError('None of the marker genes are present in the dataset.')

# Compute the mean expression of the available marker genes per cell
expr_data = adata[:, available_markers].X
if hasattr(expr_data, "toarray"):
    expr_data = expr_data.toarray()
mean_expr = np.mean(expr_data, axis=1)

# Classify cells using a threshold value; cells with average expression above the threshold are labeled as malignant
threshold = 1.0
adata.obs['malignancy'] = ['malignant' if val >= threshold else 'non-malignant' for val in mean_expr]

# Save the annotated dataset with cell classifications
adata.write('./output/o3_with_cluster/classified_cells.h5ad')

print('Classification complete. Annotated dataset saved to ./output/o3_with_cluster/classified_cells.h5ad')
