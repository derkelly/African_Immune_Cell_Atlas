#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import sys
import Spectra as spc
from Spectra import Spectra_util as spc_tl
from Spectra import K_est as kst
from Spectra import default_gene_sets
import cytopus as cp
import pickle


infile = sys.argv[1]
outfile = sys.argv[2]
modelout = sys.argv[3]

adata = sc.read_h5ad(infile)


# generate annotations

G = cp.KnowledgeBase()

celltype_of_interest = ['B','CD4-T','MAIT','M','NK','CD8-T','gdT']

global_celltypes = ['leukocyte']

G.get_celltype_processes(celltype_of_interest,
                         global_celltypes = global_celltypes,
                         get_children=True,
                         get_parents =False)

annotations = G.celltype_process_dict

# set empty sets for cell types not present in cytopus
annotations['DNT'] = {}
annotations['Progenitor'] = {}


# run spectra

# fit the model
model = spc.est_spectra(adata=adata, 
                        gene_set_dictionary=annotations, 
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
                        verbose=True,
                        num_epochs=4000
)

# save adata
adata.write(outfile)

# save model
with open(modelout, 'wb') as file:
    pickle.dump(model, file)
