import scanpy as sc
import matplotlib.pyplot as plt

# Read the annotated AnnData object which contains cell type annotations (malignant vs non-malignant)
adata = sc.read_h5ad('./output/cell_annotation_2/annotated_cells.h5ad')

# If UMAP coordinates have not been computed, compute neighbors and UMAP
if 'X_umap' not in adata.obsm:
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

# Generate UMAP plot colored by the cell type
sc.pl.umap(adata, color='cell_type', show=False)
plt.savefig('./output/cell_annotation_2/umap_annotation.png')
