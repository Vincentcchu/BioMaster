import scanpy as sc
adata = sc.read_h5ad('/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/remove_label/dataset_restricted_with_labels.h5ad')
adata.write('./output/o3_dormancy_with_malignancy/adata_initial.h5ad')
