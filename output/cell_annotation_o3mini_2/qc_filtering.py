import scanpy as sc

# Read the initial AnnData object (loaded dataset)
adata = sc.read_h5ad("./output/cell_annotation_2/initial_adata.h5ad")

# Identify mitochondrial genes assuming gene names that start with "MT-"
adata.var['mt'] = adata.var_names.str.upper().str.startswith("MT-")

# Calculate QC metrics including the percentage of mitochondrial gene counts
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter out low-quality cells
# Keep cells with <20% mitochondrial counts, >=500 total counts, and >=200 detected genes
adata = adata[adata.obs['pct_counts_mt'] < 20, :]
adata = adata[adata.obs['total_counts'] >= 500, :]
adata = adata[adata.obs['n_genes_by_counts'] >= 200, :]

# Save the filtered AnnData object
adata.write("./output/cell_annotation_2/qc_filtered_adata.h5ad")
