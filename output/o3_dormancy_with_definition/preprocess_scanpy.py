import scanpy as sc
import os

# Read the single-cell RNA-seq dataset
adata = sc.read_h5ad("/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/biomaster/output/o3_full_dataset/annotated_cells_fixed_obs_names.h5ad")

# Basic Quality Control: annotate mitochondrial genes (assuming human gene names start with 'MT-')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter cells: keep cells with at least 500 total counts
sc.pp.filter_cells(adata, min_counts=500)

# Filter cells: remove cells with >20% mitochondrial counts
adata = adata[adata.obs.pct_counts_mt < 20, :]

# Normalization: scale total counts per cell to 10,000 and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save the preprocessed dataset
output_path = "./output/o3_dormancy_with_definition/preprocessed_annotated_cells.h5ad"
adata.write(output_path)
