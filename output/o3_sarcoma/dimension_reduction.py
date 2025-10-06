import scanpy as sc

# Read the normalized AnnData object
adata = sc.read_h5ad('./output/o3_sarcoma/normalized_data.h5ad')

# Perform PCA on the normalized data
sc.pp.pca(adata)

# Compute the K-nearest neighbor graph using the PCA representation
sc.pp.neighbors(adata, use_rep='X_pca')

# Save the AnnData object with PCA results and KNN graph
adata.write('./output/o3_sarcoma/dimension_reduced_data.h5ad')
