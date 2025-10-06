import scanpy as sc
import pandas as pd

# Read the AnnData object with cluster labels and dimensionality reduction results
adata = sc.read_h5ad("./output/cell_annotation/clusters_annotated.h5ad")

# Perform differential expression analysis using t-test between clusters (assuming cluster labels are stored in 'leiden')
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test', n_genes=adata.shape[1])

# Extract marker genes for each cluster and concatenate the results
df_list = []
for group in adata.obs['leiden'].cat.categories:
    df = sc.get.rank_genes_groups_df(adata, group=group)
    df['cluster'] = group
    df_list.append(df)

markers_df = pd.concat(df_list, ignore_index=True)

# Rename columns to match desired output names
markers_df.rename(columns={'logfoldchanges': 'logFC', 'pvals_adj': 'adj_pvals'}, inplace=True)

# Filter marker genes: adjusted p-value < 0.05 and logFC > 1
filtered_markers = markers_df[(markers_df['adj_pvals'] < 0.05) & (markers_df['logFC'] > 1)]

# Save the differential expression marker gene list to CSV
filtered_markers.to_csv("./output/cell_annotation/cluster_markers.csv", index=False)
