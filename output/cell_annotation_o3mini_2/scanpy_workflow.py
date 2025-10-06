import scanpy as sc
import matplotlib.pyplot as plt

# Read the normalized AnnData object annotated with highly variable genes
adata = sc.read_h5ad('./output/cell_annotation_2/normalized_dataset.h5ad')

# Perform PCA for dimensionality reduction
sc.tl.pca(adata, svd_solver='arpack')

# Construct the KNN graph for neighborhood estimation
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Apply Leiden clustering algorithm
sc.tl.leiden(adata, resolution=1.0)

# Compute UMAP for visualization
sc.tl.umap(adata)

# Generate UMAP plot colored by Leiden clusters
ax = sc.pl.umap(adata, color='leiden', show=False)
plt.savefig('./output/cell_annotation_2/umap_clusters.png')

# Save the AnnData object with computed PCA, neighborhood graph, Leiden clusters, and UMAP coordinates
adata.write('./output/cell_annotation_2/clusters_dataset.h5ad')
