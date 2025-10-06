import scanpy as sc

# Read the single-cell RNA-seq dataset
adata = sc.read_h5ad('/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/dataset_restricted_with_labels.h5ad')

# Basic Quality Control (QC) metrics
# Annotate mitochondrial genes assuming they start with 'MT-'
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter cells: keep cells with <20% mitochondrial counts and >=500 genes expressed
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata = adata[adata.obs.n_genes_by_counts >= 500, :]

# Normalization and log transformation
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# Save the preprocessed AnnData object
adata.write('./output/o3_dormancy_with_malignancy2/preprocessed_data.h5ad')
