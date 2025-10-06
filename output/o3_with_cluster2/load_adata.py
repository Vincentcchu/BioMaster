import scanpy as sc

# Read the single-cell RNA-seq dataset from the provided .h5ad file
adata = sc.read_h5ad("/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/dataset_restricted_with_cluster.h5ad")

# Save the AnnData object to the output file
adata.write("./output/o3_with_cluster2/initial_adata.h5ad")
