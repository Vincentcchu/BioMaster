import scanpy as sc

# Read the dimensionally reduced AnnData object with PCA and UMAP embeddings
adata = sc.read_h5ad('./output/cell_annotation_2/dimensionality_reduced_adata.h5ad')

# Compute the neighborhood graph using the PCA representation
sc.pp.neighbors(adata, use_rep='X_pca')

# Perform Leiden clustering with a default resolution of 0.5
sc.tl.leiden(adata, resolution=0.5)

# Write the updated AnnData object with Leiden cluster labels to the output file
adata.write('./output/cell_annotation_2/clustered_adata.h5ad')
