import scanpy as sc

# Read the initial AnnData object
adata = sc.read_h5ad('./output/o3_sarcoma/initial_data.h5ad')

# Annotate mitochondrial genes (assuming gene names starting with 'MT-')
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate QC metrics including mitochondrial gene percentage
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter cells based on total UMIs (counts), number of detected genes, and mitochondrial percentage
# Retain cells with total counts >= 500, detected genes >= 200, and mitochondrial percentage < 20%
adata = adata[adata.obs.total_counts >= 500, :]
adata = adata[adata.obs.n_genes_by_counts >= 200, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

# Save the QC-filtered AnnData object
adata.write('./output/o3_sarcoma/qc_filtered.h5ad')
