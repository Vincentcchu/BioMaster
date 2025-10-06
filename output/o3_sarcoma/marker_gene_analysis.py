import scanpy as sc
import pandas as pd

# Read the AnnData object with clustering labels
adata = sc.read_h5ad('./output/o3_sarcoma/clustering_results.h5ad')

# Perform differential expression analysis using Scanpy's rank_genes_groups
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')

# Extract the results from the analysis
result = adata.uns['rank_genes_groups']
groups = list(result['names'].dtype.names) if hasattr(result['names'], 'dtype') else list(result['names'].keys())

all_markers = []
for group in groups:
    df = pd.DataFrame({
        'gene': result['names'][group],
        'logFC': result['logfoldchanges'][group],
        'pvals': result['pvals'][group],
        'adj_pvals': result['pvals_adj'][group]
    })
    df['cluster'] = group
    # Filter marker genes: adjusted p-value < 0.05 and logFC > 1
    filtered_df = df[(df['adj_pvals'] < 0.05) & (df['logFC'] > 1)]
    all_markers.append(filtered_df)

# Concatenate results from all clusters and save to a text file
marker_genes = pd.concat(all_markers, ignore_index=True)
marker_genes.to_csv('./output/o3_sarcoma/marker_genes.txt', sep='\t', index=False)
