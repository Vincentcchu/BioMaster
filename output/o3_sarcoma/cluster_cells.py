import scanpy as sc

# Read the dimension-reduced AnnData object (assumed to contain PCA results and KNN graph)
adata = sc.read_h5ad('./output/o3_sarcoma/dimension_reduced_data.h5ad')

# Perform Leiden clustering (using default resolution)
sc.tl.leiden(adata, resolution=1.0)

# Compute UMAP for visualization
sc.tl.umap(adata)

# Generate and save the UMAP plot with clustering labels
sc.pl.umap(adata, color=['leiden'], show=False, save='_cluster_plots.png')

# Write the AnnData object with clustering results and UMAP coordinates
adata.write('./output/o3_sarcoma/clustering_results.h5ad')
