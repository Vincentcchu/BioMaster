import scanpy as sc
import nichepca as npc

# Load the normalized AnnData object with HVGs identified
adata = sc.read_h5ad('./output/cell_annotation/adata_normalized.h5ad')

# Run nichePCA with default normalization pipeline (norm -> log1p -> agg -> pca)
npc.wf.nichepca(adata, knn=25)

# Construct neighborhood graph based on nichePCA results
sc.pp.neighbors(adata, use_rep='X_npca')

# Save the output AnnData object
adata.write('./output/cell_annotation/adata_pca_neighbors.h5ad')
