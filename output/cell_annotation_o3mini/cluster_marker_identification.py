import scanpy as sc
import pandas as pd

# Read the AnnData object with cluster labels and dimensionality reduction results
adata = sc.read_h5ad("./output/cell_annotation/clusters_annotated.h5ad")

# Perform differential expression analysis among clusters using a t-test
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')

# Extract results from the DE analysis
de_results = adata.uns['rank_genes_groups']
groups = de_results['names'].dtype.names
records = []
for grp in groups:
    names = de_results['names'][grp]
    logfoldchanges = de_results['logfoldchanges'][grp]
    pvals = de_results['pvals'][grp]
    pvals_adj = de_results['pvals_adj'][grp]
    for gene, lfc, pval, padj in zip(names, logfoldchanges, pvals, pvals_adj):
        records.append({"cluster": grp, "gene": gene, "logFC": lfc, "pvals": pval, "adj_pvals": padj})

# Create a DataFrame from the collected records
df = pd.DataFrame(records)

# Filter the results based on significance thresholds: adjusted p-value < 0.05 and logFC > 1
df_filtered = df[(df['adj_pvals'] < 0.05) & (df['logFC'] > 1)]

# Save the filtered marker gene list to a CSV file
df_filtered.to_csv("./output/cell_annotation/cluster_markers.csv", index=False)
