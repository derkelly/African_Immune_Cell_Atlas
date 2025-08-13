{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "cad99223-3fc0-4dc2-bf69-c7144df755ed",
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import pandas as pd"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "8e648c41-1a4a-4131-8149-8f7f9f294a08",
   "metadata": {},
   "outputs": [],
   "source": [
    "demux_f = \"../../RNA/RNA_demux/L89/L89_demu.best\"\n",
    "freemux_f = \"../../RNA/RNA_demux/L89/L89_freemuxlet.clust1.samples.gz\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "9499ef3c-c859-4350-a6d4-d754caa43853",
   "metadata": {},
   "outputs": [],
   "source": [
    "demux = pd.read_csv(demux_f, sep=\"\\t\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "d8447886-92ed-44d8-8ae5-0e92382b0aa6",
   "metadata": {},
   "outputs": [],
   "source": [
    "demux[['S1','S2','FRAC']] = demux['BEST.GUESS'].str.split(',', expand=True)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "0a9a0def-1a33-4c75-9125-1d4444535ff1",
   "metadata": {},
   "outputs": [],
   "source": [
    "demux_cnts = demux[[\"S1\",\"DROPLET.TYPE\"]].groupby(\"S1\").value_counts().reset_index()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "01b51f90-06a9-4ec4-8afd-f3b5d128f525",
   "metadata": {},
   "outputs": [],
   "source": [
    "demux_cnts['frac'] = list(demux_cnts.groupby(\"S1\").apply(lambda x: x['count']/np.sum(x['count'])))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "c46c3999-5ad8-47ac-b4b8-24f36d579e5b",
   "metadata": {},
   "outputs": [],
   "source": [
    "if np.min(demux_cnts[demux_cnts['DROPLET.TYPE']==\"SNG\"]['frac']) > 0.5:\n",
    "    demux.to_csv(out_f + \".demux.txt\", sep=\"\\t\", index=False)\n",
    "else:\n",
    "    freemux.to_csv(out_f + \".freemux.txt\", sep=\"\\t\", index=False)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
