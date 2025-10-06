import scanpy as sc

# Read the single-cell RNA-seq dataset and initialize the AnnData object
adata = sc.read_h5ad('./data/cell_annotation/dataset_debug_restricted.h5ad')

# Save the initialized AnnData object for downstream analyses
adata.write('./output/cell_annotation_o3mini_instruct3/initialized_data.h5ad')
