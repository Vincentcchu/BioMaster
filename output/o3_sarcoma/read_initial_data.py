import scanpy as sc

# Read the sarcoma single-cell RNA-seq dataset in h5ad format
adata = sc.read_h5ad('/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/sarcoma_dataset_restricted.h5ad')

# Write the initialized AnnData object to output
adata.write('./output/o3_sarcoma/initial_data.h5ad')
