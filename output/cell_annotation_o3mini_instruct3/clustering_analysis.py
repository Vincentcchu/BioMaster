import scanpy as sc

# Load the initialized AnnData object (assumed to be preprocessed and QC filtered)
adata = sc.read_h5ad('./output/cell_annotation_o3mini_instruct3/initialized_data.h5ad')

# Perform PCA for dimensionality reduction
sc.tl.pca(adata, svd_solver='arpack')

# Construct a KNN neighborhood graph using the PCA results
d = 40  # default number of principal components
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=d)

# Apply the Leiden clustering algorithm
sc.tl.leiden(adata, resolution=1.0)

# Compute UMAP for visualization
sc.tl.umap(adata)

# Save the AnnData object with clustering results
adata.write('./output/cell_annotation_o3mini_instruct3/clustering_results.h5ad')
