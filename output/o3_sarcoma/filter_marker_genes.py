import pandas as pd

# Read the differential expression results
df = pd.read_csv('./output/o3_sarcoma/differential_expression_results.txt', sep='\t', index_col=0)

# Filter high-confidence marker genes: adjusted p-value < 0.05 and logFC > 1
filtered_df = df[(df['adj_pvals'] < 0.05) & (df['logFC'] > 1)]

# Save the filtered marker genes
filtered_df.to_csv('./output/o3_sarcoma/filtered_marker_genes.txt', sep='\t')
