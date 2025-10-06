import scanpy as sc

# Read the QC-filtered single-cell RNA-seq dataset
adata = sc.read_h5ad('./output/cell_annotation/filtered_dataset.h5ad')

# Total count normalization (default target sum 1e4)
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data
sc.pp.log1p(adata)

# Identify highly variable genes using default parameters
sc.pp.highly_variable_genes(adata)

# Save the normalized dataset with HVG information
adata.write('./output/cell_annotation/normalized_dataset.h5ad')
