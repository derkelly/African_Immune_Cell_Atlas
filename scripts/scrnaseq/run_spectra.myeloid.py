#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import leidenalg
import glob
import sys
#import scvi
import gc
# import celltypist
# from celltypist import models
import gseapy as gp
import Spectra as spc
from Spectra import Spectra_util as spc_tl
from Spectra import K_est as kst
from Spectra import default_gene_sets
import cytopus as cp
import pickle


save_file = "write/lmtg_103024.h5ad"


adata = sc.read_h5ad(save_file)
adata_m = adata[adata.obs['Cluster_names']=="Myeloid cells"].copy()
sc.pp.highly_variable_genes(adata_m)

# generate annotations

# map cell types
ctype_map = {"Classical monocytes": "c-mono",
             "Non-classical monocytes" : "nc-mono",
             "cDCs" : "cDC",            
             "pDCs" : "p-DC"}

adata_m.obs['cytopus_names'] = adata_m.obs['SCluster_names'].apply(lambda x: ctype_map[x])

G = cp.KnowledgeBase()

celltype_of_interest = ['c-mono','nc-mono','cDC','p-DC']

global_celltypes = ['all-cells','leukocyte','M']

G.get_celltype_processes(celltype_of_interest,
                         global_celltypes = global_celltypes,
                         get_children=True,
                         get_parents =False)

annotations_m = G.celltype_process_dict


# run spectra

# fit the model
model = spc.est_spectra(adata=adata_m, 
    gene_set_dictionary=annotations_m, 
    use_highly_variable=True,
    cell_type_key="cytopus_names", 
    use_weights=True,
    lam=0.1, # varies depending on data and gene sets, try between 0.5 and 0.001
    delta=0.001, 
    kappa=None,
    rho=0.001, 
    use_cell_types=True,
    n_top_vals=50,
    label_factors=True, 
    overlap_threshold=0.2,
    clean_gs = True, 
    min_gs_num = 3,
    num_epochs=1000
)
# save adata
adata_m.write("write/lmtg_103024.wspectra.myeloid.h5ad")

# save model
with open('spectra_model.myeloid.pkl', 'wb') as file:
    pickle.dump(model, file)
