import scanpy as sc

# Read the normalized AnnData object (with highly variable genes marked) for clustering
adata = sc.read_h5ad('./output/cell_annotation_o3mini_instruct1/adata_norm.h5ad')

# Perform PCA for dimensionality reduction
sc.tl.pca(adata, svd_solver='arpack')

# Construct the K-nearest neighbor graph using the PCA results
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Apply the Leiden clustering algorithm
sc.tl.leiden(adata, resolution=1.0)

# Compute UMAP for visualization
sc.tl.umap(adata)

# Save the AnnData object with clustering and dimension reduction results
adata.write('./output/cell_annotation_o3mini_instruct1/adata_clustered.h5ad')
