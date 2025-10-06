import scanpy as sc

# Read the normalized AnnData object with selected highly variable genes
adata = sc.read('./output/cell_annotation_o3mini_instruct/normalized_hvg_dataset.h5ad')

# Perform PCA for dimensionality reduction
sc.tl.pca(adata)

# Compute the neighborhood graph
sc.pp.neighbors(adata)

# Perform clustering using the Leiden algorithm
sc.tl.leiden(adata)

# Generate UMAP visualization
sc.tl.umap(adata)

# Save the results to a new AnnData object containing the PCA, neighborhood graph, cluster labels, and UMAP coordinates
adata.write('./output/cell_annotation_o3mini_instruct/clustering_results.h5ad')
