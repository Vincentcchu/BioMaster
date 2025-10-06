import scanpy as sc
import sys

# Read the preprocessed AnnData object with normalized expression and HVG annotation
adata = sc.read_h5ad("./output/o3_dormancy_with_malignancy/adata_normalized.h5ad")

# Check for the malignant classification column in a flexible manner
if "malignant" in adata.obs.columns:
    col_name = "malignant"
elif "malignancy" in adata.obs.columns:
    col_name = "malignancy"
else:
    # Fallback to alternative possible column names
    for alt in ["tumor", "is_malignant"]:
        if alt in adata.obs.columns:
            col_name = alt
            break
    else:
        raise KeyError("No valid malignancy classification column found. Available columns: " + ", ".join(adata.obs.columns))

# Filter the AnnData object to retain only malignant cells
adata_malignant = adata[adata.obs[col_name] == True].copy()

# Save the filtered AnnData object containing only malignant cells
adata_malignant.write("./output/o3_dormancy_with_malignancy/adata_malignant.h5ad")
