import scanpy as sc

# Read the normalized AnnData object with highly variable genes computed
adata = sc.read_h5ad('./output/cell_annotation_2/normalized_adata.h5ad')

# Perform PCA
sc.tl.pca(adata, svd_solver='arpack')

# Compute the neighborhood graph using the PCA representation
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

# Compute UMAP for visualization
sc.tl.umap(adata)

# Save the AnnData object with PCA and UMAP embeddings
adata.write('./output/cell_annotation_2/dimensionality_reduced_adata.h5ad')
