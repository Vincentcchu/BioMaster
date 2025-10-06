import scanpy as sc

# Load the raw AnnData object
adata = sc.read_h5ad('./output/o3_full_dataset/initialized_dataset.h5ad')

# Identify mitochondrial genes (assuming they start with 'MT-')
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate QC metrics including the percentage of mitochondrial counts
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)

# Filter cells: keep cells with total counts (UMIs) >= 500 and mitochondrial percentage < 20%
adata = adata[adata.obs['total_counts'] >= 500, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

# Save the filtered AnnData object
adata.write('./output/o3_full_dataset/filtered_dataset.h5ad')
