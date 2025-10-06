import scanpy as sc
import pandas as pd

# Read the AnnData object with clustering annotations
adata = sc.read_h5ad('./output/o3_with_cluster2/clustered_adata.h5ad')

# Perform differential expression analysis to identify marker genes for each cluster using t-test
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')

# Extract the results into a combined DataFrame
groups = adata.obs['leiden'].cat.categories
result = adata.uns['rank_genes_groups']
df_list = []
for group in groups:
    df = pd.DataFrame({
        'gene': result['names'][group],
        'logFC': result['logfoldchanges'][group],
        'pvals': result['pvals'][group],
        'adj_pvals': result['pvals_adj'][group]
    })
    df['cluster'] = group
    df_list.append(df)
markers = pd.concat(df_list)

# Filter marker genes with adjusted p-value < 0.05 and logFC > 1
markers_filtered = markers[(markers['adj_pvals'] < 0.05) & (markers['logFC'] > 1)]

# Save the results as a CSV file
markers_filtered.to_csv('./output/o3_with_cluster2/marker_genes.csv', index=False)
