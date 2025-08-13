#!/usr/bin/env python
# coding: utf-8

# In this script, I aim to read in libraries scRNA-seq data pre-processed with SoupX, 
# performing the job of QC and data filtering, and saving the data as hdf5 files for 
# easier reading. Samples run with freemuxlet will be resolved at a later step.

import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd
import scanpy as sc
import anndata as ad
import glob
import pyarrow as pa
import sys
import scrublet


data_dir = '/project/tishkofflab/projects/SingleCellRNA/'


# # need to know whether to look for a freemuxlet file or a demuxlet file for each library
# muxlet_f = pd.read_csv(data_dir + "sample_info/muxlet_files.txt", sep="\t", names=("LIB","FILE"))
# muxlet_f = dict(zip(muxlet_f.LIB, muxlet_f.FILE))

# for now only use freemuxlet results 09/17/24


# subroutine to read in demuxlet or freemuxlet results
def get_demux(lib):





   
    demux_f = "/project/tishkofflab/projects/SingleCellRNA/RNA/RNA_demux/{}/{}_freemuxlet.clust1.samples.gz".format(lib, lib)
    demux_csv = pd.read_csv(demux_f,
                                delimiter="\t",
                                usecols=[1,4,5],
                                skiprows=[0],
                                names=['BARCODE','DROPLET_TYPE','BEST_GUESS'])

    demux_csv['LIBRARY'] = lib
    demux_csv[['S1', 'S2']] = demux_csv["BEST_GUESS"].apply(lambda x: pd.Series(str(x).split(","))[range(2)])
    demux_csv['BARCODE_LIB'] = demux_csv[['BARCODE','LIBRARY']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    demux_csv.index = demux_csv.BARCODE_LIB
    demux_csv.drop(['BEST_GUESS','BARCODE_LIB'], axis=1, inplace=True)

    return(demux_csv)



# subroutine for filtering Cell Ranger output
def calc_metrics(lib):
    
    demux_df = get_demux(lib)

    # # at CLUST if the muxlet result is a single value
    # demux_df.loc[demux_df["S1"].apply(lambda x: len(x)==1),"S1"] = demux_df[demux_df["S1"].apply(lambda x: len(x)==1)]["S1"].apply(lambda x: "CLUST" + x)
    # demux_df.loc[demux_df["S2"].apply(lambda x: len(x)==1),"S2"] = demux_df[demux_df["S2"].apply(lambda x: len(x)==1)]["S2"].apply(lambda x: "CLUST" + x)

    # read in adata
    adata = sc.read_10x_mtx("/project/tishkofflab/projects/SingleCellRNA/RNA/RNA_soupx/{}_soupx_filt".format(lib), var_names="gene_symbols", cache=False)
    adata.obs['index'] = [str(idx) + "_" + lib for idx in adata.obs.reset_index()['index']]
    adata.obs = adata.obs.set_index('index')
    
    # get intersecting barcodes between adata and demux results
    idx = list(set(adata.obs.index).intersection(set(demux_df.index)))

    # restrict to matching barcodes
    adata = adata[idx]
    demux_df = demux_df.loc[idx]

    # set relevant sample information
    adata.obs[['TYPE', 'LIBRARY', 'S1', 'S2']] = demux_df.loc[adata.obs.index, ['DROPLET_TYPE', 'LIBRARY', 'S1', 'S2']]  
    
    # annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # annotate the group of hemoglobin genes as 'hba'
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True)
    
    adata.obs['log10_total_counts'] = [np.log10(x) for x in adata.obs['total_counts']]
    
    return(adata)



################## MAIN ##################


lib = sys.argv[1]


# read in adata and calculate metrics
adata = calc_metrics(lib)


# create filtered dataset
adata_filt = adata[((adata.obs.TYPE == "SNG") &
                    (adata.obs.pct_counts_mt < 10) &
                    (adata.obs.total_counts > 500) & 
                    (adata.obs.n_genes_by_counts > 250))]


# run scrublet
sc.external.pp.scrublet(adata_filt, threshold=0.3)


# # set doublet calls from scrublet
# adata.obs['doublet_score'] = np.NaN
# adata.obs['predicted_doublet'] = False

# adata.obs.loc[adata_filt.obs.index,'doublet_score'] = adata_filt.obs.loc[adata_filt.obs.index,'doublet_score']
# adata.obs.loc[adata_filt.obs.index,'predicted_doublet'] = adata_filt.obs.loc[adata_filt.obs.index,'predicted_doublet']

adata_filt = adata_filt[~adata_filt.obs['predicted_doublet']]

# write results
adata_filt.write(data_dir + "/RNA/RNA_filt/{}_soupx_filt.h5ad".format(lib), compression="gzip")

