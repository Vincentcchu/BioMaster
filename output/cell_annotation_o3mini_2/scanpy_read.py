import scanpy as sc

# Read the single-cell RNA-seq dataset from the provided .h5ad file
adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad')

# Save the AnnData object to the specified output file
adata.write('./output/cell_annotation_2/initial_adata.h5ad')

print('AnnData object successfully saved to ./output/cell_annotation_2/initial_adata.h5ad')
