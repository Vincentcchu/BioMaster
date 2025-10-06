import scanpy as sc

# Read the filtered AnnData object
adata = sc.read_h5ad('./output/o3_sarcoma/filtered_data.h5ad')

# Normalize the total counts per cell to 1e4 and log-transform the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes (HVGs) using default parameters
sc.pp.highly_variable_genes(adata)

# Save the normalized AnnData object with HVGs
adata.write('./output/o3_sarcoma/normalized_data.h5ad')
