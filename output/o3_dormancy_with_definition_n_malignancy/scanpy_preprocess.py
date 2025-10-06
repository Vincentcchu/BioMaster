import scanpy as sc

# Read the single-cell RNA-seq dataset
adata = sc.read_h5ad("/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/biomaster/output/cell_annotation_o3mini/malignant_classification.h5ad")

# Annotate mitochondrial genes assuming gene names starting with MT-
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calculate quality control metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# Filter out cells with more than 20 percent mitochondrial counts and with less than 500 total counts
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata = adata[adata.obs.total_counts >= 500, :]

# Normalize: scale total counts per cell to 10000 and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Scale the data to unit variance and zero mean with clipping at 10
sc.pp.scale(adata, max_value=10)

# Write the processed dataset to a temporary directory with more available space
adata.write("/tmp/o3_dormancy_with_definition_n_malignancy/processed_dataset.h5ad")
