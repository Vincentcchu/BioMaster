import scanpy as sc
import pandas as pd

# Load the AnnData object with clustering annotations
adata = sc.read_h5ad("./output/o3_with_cluster2/clustered_adata.h5ad")

# Load the cell malignancy classification file; assuming the CSV file has cell barcodes/index as the first column
malignancy = pd.read_csv("./output/o3_with_cluster2/cell_malignancy_annotation.csv", index_col=0)

# Add the malignant status annotation to the cell metadata
adata.obs = adata.obs.join(malignancy, how="left")

# Save the final annotated AnnData object with malignant cell annotations
adata.write("./output/o3_with_cluster2/final_annotated_adata.h5ad")
