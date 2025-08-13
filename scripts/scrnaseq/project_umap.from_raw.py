#!/usr/bin/env python
# coding: utf-8

# In this script, I will run harmony and perform UMAP projection and leiden clustering

import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import leidenalg
import sys


# get input
adata_f = sys.argv[1]
batch_f = sys.argv[2]
out_f = sys.argv[3]


# read in adata
adata = sc.read_h5ad(adata_f)


# add batch information
batches = pd.read_csv("../../misc/lib_groups.txt", sep="\t")
adata.obs['BATCH'] = adata.obs[['library']].reset_index() \
    .merge(batches, left_on='library', right_on='LIBRARY') \
    .set_index('index').loc[adata.obs.index,'BATCH']


adata.layers['raw_counts'] = adata.X.copy()


# Normalize data
# Normalizing to 10000 total counts
sc.pp.normalize_total(adata, target_sum=1e4)
# Logarithmize the data
sc.pp.log1p(adata)


# highly variable genes
sc.pp.highly_variable_genes(adata, batch_key='BATCH')


# PCA
sc.tl.pca(adata, use_highly_variable=True, svd_solver='arpack')


# Run Harmony batch correction
sce.pp.harmony_integrate(adata, 'BATCH')


# Calculate neighbors
sc.pp.neighbors(adata, n_neighbors=50, n_pcs=30, use_rep='X_pca_harmony', key_added='neighbors_harmony')


# Leiden clustering
sc.tl.leiden(adata, neighbors_key="neighbors_harmony", flavor="igraph", n_iterations=2)


# UMAP projection
sc.tl.umap(adata, min_dist=0.25, neighbors_key='neighbors_harmony')


# Write data
adata.write(out_f)
