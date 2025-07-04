import scanpy as sc

# Read the normalized and highly variable gene selected dataset
adata = sc.read_h5ad("./output/cell_annotation/normalized_dataset.h5ad")

# Perform principal component analysis using scanpy pca
sc.tl.pca(adata, svd_solver='arpack')

# Construct the neighborhood graph of cells
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Compute Leiden clustering
sc.tl.leiden(adata, resolution=1.0)

# Compute UMAP for visualization
sc.tl.umap(adata)

# Save the AnnData object with clustering results and UMAP coordinates
adata.write("./output/cell_annotation/clusters_annotated.h5ad")
