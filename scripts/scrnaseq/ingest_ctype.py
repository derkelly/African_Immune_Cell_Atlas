#!/usr/bin/env python
# coding: utf-8

# In this script, I will ingest the cell type labelings of dataset A into dataset B.

import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import leidenalg
import sys


# get input
obs_col = sys.argv[1]
adata_Ai = sys.argv[2]
adata_Bi = sys.argv[3]
out_file = sys.argv[4]


# read in data
adata_A = sc.read_h5ad(adata_Ai)
adata_B = sc.read_h5ad(adata_Bi)


# if a raw_counts layer is in the object, replace X with this
if 'raw_counts' in adata_B.layers:
    adata_B.X = adata_B.layers['raw_counts'].copy()


# Normalize data
sc.pp.normalize_total(adata_B, target_sum=1e4)


# subset to shared genes
var_names = adata_A.var_names.intersection(adata_B.var_names)
adata_B = adata_B[:, var_names]


# Regress, and scale data
sc.pp.regress_out(adata_B, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata_B, max_value=10.0)


# perform ingestation
sc.tl.ingest(adata_B, adata_A, obs=obs_col)


# write output (cell index and cell type)
adata_B.obs[[obs_col]].to_csv(out_file)
