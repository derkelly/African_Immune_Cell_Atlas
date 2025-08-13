#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import pandas as pd


demux_f = sys.argv[1]
freemux_f = sys.argv[2]
out_f = sys.argv[3]


demux = pd.read_csv(demux_f, sep="\t")
freemux = pd.read_csv(freemux_f, sep="\t")


demux_llk_diff = demux[demux['DROPLET.TYPE']=="SNG"].groupby('BEST.GUESS')['DIFF.LLK.BEST.NEXT'].agg('median')
freemux_llk_diff = freemux[freemux['DROPLET.TYPE']=="SNG"].groupby('BEST.GUESS')['DIFF.LLK.BEST.NEXT'].agg('median')



if ((np.min(demux_llk_diff) > 0) | (np.min(demux_llk_diff) > np.min(freemux_llk_diff))):
    demux.to_csv(out_f + ".demux.txt", sep="\t", index=False)
else:
    freemux.to_csv(out_f + ".freemux.txt", sep="\t", index=False)
