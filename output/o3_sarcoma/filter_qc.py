import scanpy as sc
import anndata as ad

# Read the two AnnData objects
adata1 = sc.read_h5ad("./output/o3_sarcoma/initialized_data.h5ad")
adata2 = sc.read_h5ad("/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/sarcoma_dataset_restricted.h5ad")

# Concatenate the two datasets
adata = ad.concat([adata1, adata2], join='outer', merge='same')

# Identify mitochondrial genes (assumes gene names starting with 'MT-')
adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter cells:
#   - Mitochondrial gene percentage less than 20%
#   - Total UMIs at least 500
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata = adata[adata.obs.total_counts >= 500, :]

# Save the filtered AnnData object
adata.write("./output/o3_sarcoma/filtered_data.h5ad")
