# module load python
# module load cellranger
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

snap.pl.umap(atac, color="cell_type", interactive=False, out_file='results/umap_atac.png')

atac.close()

### Perform data integration
# First we need to assign genomic coordinates to the genes in scRNA-seq experiment.

scglue.data.get_gene_annotation(
    rna, gtf="data/gencode.v40.chr_patch_hapl_scaff.annotation.gtf.gz",
    gtf_by="gene_name"
)
rna = rna[:, rna.var.dropna().index.to_numpy()]
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()

atac = ad.read("data/10x-Multiome-Pbmc10k-ATAC.h5ad")
sc.pp.neighbors(atac, use_rep="X_spectral")
atac

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1])
atac.var["chromEnd"] = split.map(lambda x: x[2])
atac.var.head()

# graph = scglue.genomics.rna_anchored_prior_graph(rna, atac)
graph = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
graph

scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)
scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_spectral"
)

graph = graph.subgraph(itertools.chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
))

# This is patch for numpy >= 1.21 to fix a bug mentioned in
# https://github.com/lmcinnes/pynndescent/issues/163
# import pynndescent
# import numba
# import numpy as np
# @numba.njit(fastmath=True)
# def correct_alternative_cosine(ds):
#     result = np.empty_like(ds)
#     for i in range(ds.shape[0]):
#         result[i] = 1.0 - np.power(2.0, ds[i])
#     return result
# pynn_dist_fns_fda = pynndescent.distances.fast_distance_alternatives
# pynn_dist_fns_fda["cosine"]["correction"] = correct_alternative_cosine
# pynn_dist_fns_fda["dot"]["correction"] = correct_alternative_cosine

glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, graph,
    fit_kws={"directory": "glue"}
)

rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

combined = ad.concat([rna, atac])
sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)

sc.pl.umap(combined, color=["cell_type", "domain"], wspace=0.45)
