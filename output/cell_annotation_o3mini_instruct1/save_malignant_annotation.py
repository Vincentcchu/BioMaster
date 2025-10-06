import scanpy as sc

# Load the clustered AnnData object with dimension reduction coordinates and cluster labels
adata_clustered = sc.read_h5ad("./output/cell_annotation_o3mini_instruct1/adata_clustered.h5ad")

# Load the AnnData object with malignant classification appended to cell metadata
adata_annotated = sc.read_h5ad("./output/cell_annotation_o3mini_instruct1/adata_annotated.h5ad")

# Merge malignant classification into the clustered AnnData object
# Assuming that the cells are in the same order and have matching indices
adata_clustered.obs['malignant'] = adata_annotated.obs['malignant']

# Save the final annotated AnnData object
adata_clustered.write("./output/cell_annotation_o3mini_instruct1/malignant_annotation.h5ad")
