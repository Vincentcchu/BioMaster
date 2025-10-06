import scanpy as sc

# Read the single-cell RNA-seq dataset
adata = sc.read_h5ad("./data/cell_annotation/dataset_debug_restricted.h5ad")

# Identify mitochondrial genes (assuming gene names starting with 'MT-')
adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')

# Compute QC metrics, including the percentage of counts from mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter out low-quality cells based on total UMIs and mitochondrial gene percentage
# Keeping cells with at least 500 total counts and less than 20% mitochondrial counts
adata = adata[adata.obs['total_counts'] >= 500, :]
adata = adata[adata.obs['pct_counts_mt'] < 20, :]

# Optionally, additional filtering on detected genes can be applied if needed (e.g., adata.obs['n_genes_by_counts'] >= 200)

# Save the filtered AnnData object
adata.write("./output/cell_annotation_o3mini_instruct2/filtered_dataset.h5ad")
