import scanpy as sc
import pandas as pd

# Read the clustered AnnData object with Leiden cluster labels
adata = sc.read_h5ad("./output/cell_annotation_2/clustered_adata.h5ad")

# Perform differential expression analysis between clusters using t-test
sc.tl.rank_genes_groups(adata, groupby="leiden", method="t-test")

df = sc.get.rank_genes_groups_df(adata, group=None)

# Rename columns to match specifications: logFC, pvals, and adj_pvals
df.rename(columns={"logfoldchanges": "logFC", "pvals_adj": "adj_pvals"}, inplace=True)

# Save the differential expression results to CSV
df.to_csv("./output/cell_annotation_2/deg_results.csv", index=False)
