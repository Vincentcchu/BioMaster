import scanpy as sc

# Read the QC filtered AnnData object
adata = sc.read_h5ad("./output/o3_sarcoma/qc_filtered.h5ad")

# Normalize total counts per cell to a target sum (default 1e4) and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes (default parameters)
sc.pp.highly_variable_genes(adata)

# Save the normalized data with HVG information
adata.write("./output/o3_sarcoma/normalized_data.h5ad")
