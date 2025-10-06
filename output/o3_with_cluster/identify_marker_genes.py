import scanpy as sc
import pandas as pd

# Read the clustered single-cell dataset with cell clusters annotated
adata = sc.read_h5ad('./output/o3_with_cluster/clustered_dataset.h5ad')

# Perform differential expression analysis to identify marker genes for each cluster
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test', n_genes=adata.shape[1])

# Initialize an empty DataFrame to store marker genes per cluster
marker_list = []

# Get unique clusters from the AnnData object
clusters = adata.obs['leiden'].unique()

# For each cluster, extract marker genes and filter based on criteria (adj p-value < 0.05 and logFC > 1)
for cluster in clusters:
    df = sc.get.rank_genes_groups_df(adata, group=cluster)
    # Rename columns to match desired output
    df.rename(columns={'pvals_adj': 'adj_pvals', 'logfoldchanges': 'logFC'}, inplace=True)
    # Filter markers
    df_filtered = df[(df['adj_pvals'] < 0.05) & (df['logFC'] > 1)]
    df_filtered['cluster'] = cluster
    marker_list.append(df_filtered)

# Concatenate marker genes from all clusters into a single DataFrame
if marker_list:
    markers = pd.concat(marker_list, ignore_index=True)
else:
    markers = pd.DataFrame()

# Save the final list of marker genes to the output file
markers.to_csv('./output/o3_with_cluster/marker_genes.txt', sep='\t', index=False)
