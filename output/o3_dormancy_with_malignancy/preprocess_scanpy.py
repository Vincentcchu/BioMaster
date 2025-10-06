import scanpy as sc

# Read the initial AnnData objects
adata_initial = sc.read_h5ad('./output/o3_dormancy_with_malignancy/adata_initial.h5ad')
adata_restricted = sc.read_h5ad('/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/dataset_restricted_with_labels.h5ad')

# Mark mitochondrial genes by creating a boolean column in .var
adata_initial.var["mt"] = adata_initial.var_names.str.startswith("MT-")
adata_restricted.var["mt"] = adata_restricted.var_names.str.startswith("MT-")

# Merge the two datasets (adding a batch column to distinguish them if needed)
ad = adata_initial.concatenate(adata_restricted, batch_key="source", batch_categories=["initial", "restricted"])

# Calculate quality control metrics using the correct qc variable name
sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# Filter cells based on QC metrics
ad = ad[ad.obs.total_counts >= 500, :]
ad = ad[ad.obs.pct_counts_mt < 20, :]

# Normalize total counts per cell
sc.pp.normalize_total(ad, target_sum=1e4)

# Log-transform the data
sc.pp.log1p(ad)

# Select highly variable genes (HVGs)
sc.pp.highly_variable_genes(ad, flavor="seurat", n_top_genes=2000)

# Save the preprocessed AnnData object
ad.write('./output/o3_dormancy_with_malignancy/adata_normalized.h5ad')
