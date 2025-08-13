#!/usr/bin/env python
# coding: utf-8

# In this script, I aim to read in libraries scRNA-seq data, add covariate information to cells, and perform the job of QC and data filtering, and save the data as hdf5 files for easier reading.

import numpy as np
import scipy
from scipy import sparse
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import sys
import memento
import gc
    
            
# set important variables
adata_f = sys.argv[1]
ctypes = sys.argv[2]
pops = sys.argv[3].split("_")
conds = sys.argv[4].split("_")
treatment_col = sys.argv[5]
out_f = sys.argv[6]


# read in adata and subset
adata = sc.read_h5ad(adata_f, backed='r')

adata_test = adata[(adata.obs['cytopus_names']==ctypes) & 
                   (adata.obs['Ethnicity'].isin(pops)) &
                   (adata.obs['COND'].isin(conds))].to_memory()


# read in metrics and convert "Sequence Saturation" column to ratio
metrics = pd.read_csv("/project/tishkofflab/projects/SingleCellRNA/sample_info/10x.metrics_summary.all.csv", thousands=",")
metrics['Sequencing Saturation'] = metrics['Sequencing Saturation'].str.rstrip('%').astype('float') / 100.0

# add sequencing saturation information to adata_test
adata_test.obs['seq_sat'] = adata_test.obs[['library']].reset_index() \
    .merge(metrics[['LIB','Sequencing Saturation']], left_on='library', right_on='LIB') \
    .set_index('index').loc[adata_test.obs.index,'Sequencing Saturation']
adata_test.obs['capture_rate'] = adata_test.obs['seq_sat'] * 0.25


# if the test is between conditions, use binary testing with replicates
if (treatment_col=="COND"):

    # test condition is the one that is not control
    test_cond = conds[conds != "CTL"]

    # make the condition numeric
    adata_test.obs[test_cond] = adata_test.obs['COND'].apply(lambda x: 0 if x == 'CTL' else 1)

    # set capture rate and get covariates
    memento.setup_memento(adata_test, q_column='capture_rate')
    memento.create_groups(adata_test, label_columns=[test_cond, 'TID'])
    memento.compute_1d_moments(adata_test, min_perc_group=.9)
    sample_meta = memento.get_groups(adata_test)
    sample_meta['TID'] = sample_meta['TID'].astype('category') # make sure to not confuse ourselves in case replicate labels are numbers
    
    treatment_df = sample_meta[[test_cond]]
    cov_df = pd.get_dummies(sample_meta['TID'].astype('category'))

    # run memento
    memento.ht_1d_moments(
        adata_test, 
        treatment=treatment_df,
        covariate=cov_df,
        num_boot=5000, 
        num_cpus=16)

    result_cond = memento.get_1d_ht_result(adata_test)

    result_cond.to_csv(out_f)


# if the test is between ethnic groups
elif (treatment_col=="Ethnicity"):

    test_pop = pops[1]

    # set up covariates
    sample_meta = adata_test.obs[['TID','SEX','AGE','BATCH','Ethnicity']].drop_duplicates().set_index('TID')
    sample_meta[test_pop] = sample_meta['Ethnicity'].apply(lambda x: 0 if x == pops[0] else 1)

    batch = pd.get_dummies(sample_meta['BATCH'].astype('category'))
    sex = pd.get_dummies(sample_meta['SEX'].astype('category'))
    cov_df = sample_meta[['AGE']].merge(sex['M'], left_index=True, right_index=True) \
                                 .merge(batch, left_index=True, right_index=True)

    # fills the place of the "snps" dataframe
    ethnicity_df = sample_meta[[test_pop]]

    # fills the place of "gene_snp_pairs," will test "is Fulani" for every gene
    gene_pop_pairs = pd.DataFrame({'gene': adata_test.var_names,
                                   'SNP': test_pop})


    result_pop = memento.run_eqtl(
        adata=adata_test,
        snps=ethnicity_df,
        cov=cov_df,
        gene_snp_pairs=gene_pop_pairs,
        donor_column='TID',
        num_cpu=16,
        num_blocks=1, # increase this if you run out of memory.
        num_boot=5000,
    )

    result_pop.to_csv(out_f)
