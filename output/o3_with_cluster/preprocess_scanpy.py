import scanpy as sc
import os

# Step 1: Data Reading and Initialization
adata = sc.read_h5ad("/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/dataset_restricted_with_cluster.h5ad")

# Step 2: Quality Control (QC)
sc.pp.calculate_qc_metrics(adata, inplace=True)
# Filter cells based on total counts and mitochondrial gene percentage
adata = adata[adata.obs["total_counts"] >= 550, :]
if "pct_counts_mt" in adata.obs.columns:
    adata = adata[adata.obs["pct_counts_mt"] < 20, :]

# Step 3: Normalization and Highly Variable Gene Selection
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# Step 4: Dimensionality Reduction and Visualization
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# Step 5: Clustering Analysis
sc.tl.leiden(adata)

# Save the preprocessed AnnData object
adata.write("./output/o3_with_cluster/preprocessed_dataset.h5ad")
