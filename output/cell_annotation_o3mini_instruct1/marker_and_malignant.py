import scanpy as sc
import numpy as np

# Read the clustered AnnData object
adata = sc.read_h5ad('./output/cell_annotation_o3mini_instruct1/adata_clustered.h5ad')

# Perform differential expression analysis to identify marker genes per cluster
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')

def get_top_markers(adata, group, top_n=10):
    # Extract the DE results for each group
    gene_names = adata.uns['rank_genes_groups']['names'][group][:top_n]
    return set(gene_names)

# Define a set of known malignant cell markers (this is a dummy list, adjust as needed)
known_malignant_markers = {'MAL1', 'MAL2', 'MAL3'}

# Determine malignant status per cluster
cluster_labels = adata.obs['leiden'].unique()
malignant_classification = {}
for clust in cluster_labels:
    top_markers = get_top_markers(adata, clust, top_n=10)
    # If there's any intersection with known malignant markers, classify cluster as malignant
    if len(top_markers.intersection(known_malignant_markers)) > 0:
        malignant_classification[clust] = 1
    else:
        malignant_classification[clust] = 0

# Append malignant status to each cell in the AnnData object's metadata
adata.obs['malignant'] = adata.obs['leiden'].map(malignant_classification).astype('int')

# Save the annotated AnnData object with malignant classification
adata.write('./output/cell_annotation_o3mini_instruct1/adata_annotated.h5ad')

