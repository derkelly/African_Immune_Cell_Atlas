#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scipy
from scipy import sparse
from scipy.sparse import csr_matrix
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import glob
import sys

counts_file = '/project/tishkofflab/projects/SingleCellRNA/test/test123/write/qc_020725_raw.h5ad'

adata = sc.read_h5ad(counts_file)

gene_counts = np.unique(adata.X.indices, return_counts=True)

idx_gte_10 = gene_counts[0][gene_counts[1] > 9]

adata_gte10 = adata[:,adata.var_names[idx_gte_10]]

# mtx_file = sparse.csr_matrix(adata.X)

# Write sparse matrix
PATH = "/project/tishkofflab/projects/SingleCellRNA/test/test123/write/scran/input/"
FILE = "counts_v1.mtx"

scipy.io.mmwrite(PATH + FILE, adata.X)

g_file = pd.Series(adata.var_names)

# Write gene names
FILE = "genes.csv"

g_file.to_csv(PATH + FILE, index=False, header=False)

bc_file = pd.Series(adata.obs_names)

FILE = "cells.csv"

bc_file.to_csv(PATH + FILE, index=False, header=False)
