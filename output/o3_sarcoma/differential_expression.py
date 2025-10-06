import scanpy as sc
import pandas as pd

# Load the AnnData object with clustering results (assumed to contain the cluster labels in adata.obs['leiden'])
adata = sc.read_h5ad("./output/o3_sarcoma/clustering_results.h5ad")

# Perform differential expression analysis using t-test on the clusters
sc.tl.rank_genes_groups(adata, groupby="leiden", method="t-test")

# Extract the results
result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names

table_list = []
for group in groups:
    genes = result["names"][group]
    logfc = result["logfoldchanges"][group]
    pvals = result["pvals"][group]
    padj = result["pvals_adj"][group]
    for gene, lfc, pval, adj in zip(genes, logfc, pvals, padj):
        table_list.append({"group": group, "gene": gene, "logFC": lfc, "pvals": pval, "adj_pvals": adj})

# Create a DataFrame and save as text file with tab delimiter
df = pd.DataFrame(table_list)
df.to_csv("./output/o3_sarcoma/differential_expression_results.txt", sep="\t", index=False)

print("Differential expression analysis completed. Results saved to ./output/o3_sarcoma/differential_expression_results.txt")
