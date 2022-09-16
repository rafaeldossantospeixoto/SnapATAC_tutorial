# module load python
# cd hpc/SnapATAC_tutorial/

### Introduction
# In addition to SnapATAC2, we will utilize scanpy and scglue to perform the integration.
import scanpy as sc
import snapatac2 as snap
import anndata as ad
import scglue
import itertools

### Analyze gene expression data
rna = sc.read_h5ad("data/10x-Multiome-Pbmc10k-RNA.h5ad")
rna

# First we backup the raw counts in anndata and perform the preprocessing.
rna.layers["counts"] = rna.X.copy()
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna)
rna = rna[:, rna.var.highly_variable]
sc.pp.scale(rna, max_value=10)

# We then perform the dimension reduction, followed by UMAP embedding.
sc.tl.pca(rna)
sc.pp.neighbors(rna)
sc.tl.umap(rna)

sc.pl.umap(rna, color="cell_type")

### Analyze chromatin accessibility data
atac = snap.read("data/10x-Multiome-Pbmc10k-ATAC.h5ad")
atac

snap.tl.spectral(atac, features=None, n_comps=30)
snap.tl.umap(atac)

snap.pl.umap(atac, color="cell_type", interactive = False)

atac.close()

### Perform data integration
# First we need to assign genomic coordinates to the genes in scRNA-seq experiment.
