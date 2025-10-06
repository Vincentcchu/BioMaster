import scanpy as sc

# Read the filtered single-cell dataset
adata = sc.read_h5ad("./output/cell_annotation_o3mini_instruct2/filtered_dataset.h5ad")

# Perform PCA on the dataset
sc.tl.pca(adata, svd_solver="arpack")

# Compute the neighborhood graph of cells
sc.pp.neighbors(adata)

# Compute UMAP embedding for visualization
sc.tl.umap(adata)

# Perform clustering using the Leiden algorithm
sc.tl.leiden(adata)

# Save the dimension-reduced dataset with clustering labels
adata.write("./output/cell_annotation_o3mini_instruct2/dim_reduced_dataset.h5ad")
