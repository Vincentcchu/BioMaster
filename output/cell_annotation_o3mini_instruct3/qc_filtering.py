import scanpy as sc

# Read the initialized AnnData object
adata = sc.read_h5ad("./output/cell_annotation_o3mini_instruct3/initialized_data.h5ad")

# Create a boolean column in adata.var to mark mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calculate QC metrics using the boolean marker column
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# Filter cells: retain cells with mitochondrial percentage less than 20 and total counts at least 500
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata = adata[adata.obs.total_counts >= 500, :]

# Write the filtered AnnData object to output file
adata.write("./output/cell_annotation_o3mini_instruct3/qc_filtered_data.h5ad")
