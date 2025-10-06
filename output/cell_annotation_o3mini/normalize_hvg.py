import scanpy as sc

# Read the QC-filtered single-cell RNA-seq dataset
adata = sc.read_h5ad('./output/cell_annotation/filtered_dataset.h5ad')

# Total count normalization with a target sum (default: 1e4)
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data
sc.pp.log1p(adata)

# Identify highly variable genes (HVGs) using default parameters
sc.pp.highly_variable_genes(adata)

# Optionally, you can store the HVG information in the AnnData object
# For example, adata.var['highly_variable'] is a boolean array indicating HVGs

# Write the normalized and HVG-selected dataset to an output file
adata.write('./output/cell_annotation/normalized_dataset.h5ad')
