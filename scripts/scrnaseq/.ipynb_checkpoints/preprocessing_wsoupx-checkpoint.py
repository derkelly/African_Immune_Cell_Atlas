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


results_dir = 'hdf5_soupx'  # the directory that will store the analysis results

data_dir = '/project/tishkofflab/projects/SingleCellRNA/'


# dataframe of the directories containing demuxlet output and the first and last library in that dir
L5_L22 = [data_dir + "2022_05_05_NovaSeq_L5-L22/L5-L22_demuxlet/L" + 
           str(i) for i in range(5,23)]

L23_L43 = [data_dir + "Ning_NovASeq_L23-L43_0812/L23-L43_demuxlet/L" + 
           str(i) for i in range(23,44)]

L44_L64 = [data_dir + "Ning_NovASeq_L44-L64_0822/L44-L64_demuxlet/L" + 
           str(i) for i in range(44,65)]

L65_L85 = [data_dir + "Ning_NovASeq_L65-85_0829/L65-L85_demuxlet/L" + 
           str(i) for i in range(65,86)]

L86_L106 = [data_dir + "20230505_NovaSeq_L86-L106/L86-L106_demuxlet/L" +
           str(i) for i in range(86,107)]

L107_L127 = [data_dir + "tishkoff_L107-L127_NovaSeq_20230830/L107-L127_demuxlet/L" +
           str(i) for i in range(107,128)]

L128_L149 = [data_dir + "tishkoff_L128-L149_NovaSeq_20230912/L128-L149_demuxlet/L" +
           str(i) for i in range(128,150)]

L150_L157 = [data_dir + "tishkoff_L147_L150-157_MR13-20_Nova_20240112/L147_L150-L157_demuxlet/L" +
             str(i) for i in range(150,158)]

L158_L178 = [data_dir + "tishkoff_L158-L178_NovaSeq_20240119/L158-L178_demuxlet/L" +
             str(i) for i in range(158,179)]



paths = L5_L22 + L23_L43 + L44_L64 + L65_L85 + L86_L106 + L107_L127 + L128_L149 + L150_L157 + L158_L178
libs = ["L" + str(i) for i in range(5,179)]



demux_paths = dict(zip(libs, paths))
demux_paths['L147B'] = data_dir + "tishkoff_L147_L150-157_MR13-20_Nova_20240112/L147_L150-L157_demuxlet/L147"


# need to know whether to look for a freemuxlet file or a demuxlet file for each library
muxlet_f = pd.read_csv(data_dir + "sample_info/muxlet_files.txt", sep="\t", names=("LIB","FILE"))
muxlet_f = dict(zip(muxlet_f.LIB, muxlet_f.FILE))


# subroutine to read in demuxlet or freemuxlet results
def get_demux(lib):
    
    demux_f = demux_paths[lib] + "/" + muxlet_f[lib]
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
def filter_adata(lib):
    
    demux_df = get_demux(lib)

    # at CLUST if the muxlet result is a single value
    demux_df.loc[demux_df["S1"].apply(lambda x: len(x)==1),"S1"] = demux_df[demux_df["S1"].apply(lambda x: len(x)==1)]["S1"].apply(lambda x: "CLUST" + x)
    demux_df.loc[demux_df["S2"].apply(lambda x: len(x)==1),"S2"] = demux_df[demux_df["S2"].apply(lambda x: len(x)==1)]["S2"].apply(lambda x: "CLUST" + x)

    # read in adata
    adata = sc.read_10x_mtx("/project/tishkofflab/projects/SingleCellRNA/analysis23/soupx/soupx_{}_filt".format(lib), var_names="gene_symbols", cache=False)
    adata.obs['index'] = [str(idx) + "_" + lib for idx in adata.obs.reset_index()['index']]
    adata.obs = adata.obs.set_index('index')
    
    # set relevant sample information
    adata.obs[['TYPE', 'LIBRARY', 'S1', 'S2']] = demux_df[['DROPLET_TYPE', 'LIBRARY', 'S1', 'S2']]  
    
    # annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # annotate the group of hemoglobin genes as 'hba'
    adata.var['hba'] = adata.var_names.str.endswith(('HBB', 'HBA1', 'HBA2'))  
    sc.pp.calculate_qc_metrics(adata, qc_vars=['hba'], percent_top=None, log1p=False, inplace=True)
    
    adata.obs['log10_total_counts'] = [np.log10(x) for x in adata.obs['total_counts']]
    
    return(adata)


# calculate number of cells failing different QC metrics
# they are ranked as follows:
# 1. tagged as SNG using demuxlet or freemuxlet
# 2. is *not* predicted to be a doublet using scrublet
# 3. has more than 500 total UMIs
# 4. has more than 250 genes identified
# 5. has less than 20% mitochondrial reads
# 6. has less than 1% hemoglobin reads

# def calc_qc(adata):
    
#     d = {'demux': np.sum(adata.obs.TYPE=="DBL"),
#          'scrublet': np.sum(~(adata.obs.TYPE=="DBL") & 
#                             adata.obs.predicted_doublet),
#          'total_counts': np.sum(~((adata.obs.TYPE=="DBL") | 
#                                   adata.obs.predicted_doublet) & 
#                                 (adata.obs.total_counts<=500)),
#          'n_genes_by_counts': np.sum(~((adata.obs.TYPE=="DBL") | 
#                                        adata.obs.predicted_doublet | 
#                                        (adata.obs.total_counts<=500)) & 
#                                      (adata.obs.n_genes_by_counts<=250)),
#          'pct_counts_mt': np.sum(~((adata.obs.TYPE=="DBL") | 
#                                    adata.obs.predicted_doublet | 
#                                    (adata.obs.total_counts<=500) | 
#                                    (adata.obs.n_genes_by_counts<=250)) & 
#                                  (adata.obs.pct_counts_mt>=20)),
#          'pct_counts_hba': np.sum(~((adata.obs.TYPE=="DBL") | 
#                                     adata.obs.predicted_doublet | 
#                                     (adata.obs.total_counts<=500) | 
#                                     (adata.obs.n_genes_by_counts<=250) | 
#                                     (adata.obs.pct_counts_mt>=20)) & 
#                                   (adata.obs.pct_counts_hba>1))
#           }
    
#     qcs = pd.DataFrame(data=d, index=[0])
#     return(qcs)


################## MAIN ##################

lib = sys.argv[1]

adata = filter_adata(lib)


adata = adata[((adata.obs.TYPE == "SNG") &
               (adata.obs.pct_counts_mt < 20) &
               (adata.obs.pct_counts_hba < 1) & 
               (adata.obs.total_counts > 500) & 
               (adata.obs.n_genes_by_counts > 250))]


sc.external.pp.scrublet(adata, threshold=0.3)


# qcs = calc_qc(adata)
# qcs['lib'] = lib
# qcs.to_csv('qc/qc_filt.wsoupx.csv', sep='\t', mode='a', index=False, header=False)


adata[~adata.obs.predicted_doublet].write(results_dir + "/filtered/" + lib + "_wsoupx_filt.h5ad", compression="gzip")
