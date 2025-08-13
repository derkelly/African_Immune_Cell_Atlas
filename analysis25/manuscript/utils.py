#!/usr/bin/env python
# coding: utf-8
# utils.py


cond_cols = {'CTL': "#B3FFFF",
             'IFN': "#FF66FF",
             'LPS': "#FFA64D"}

pop_cols = {'CHG': "#AAFF00",
            'Fulani': "#E0508B",
            'Tikari': "#0073E6"}

ct_map = {'B cells': 'B',
          'CD4+ T cells': 'CD4T',
          'DN T cells': 'DNT',
          'MAIT cells': 'MAIT',
          'Myeloid cells': 'Myeloid',
          'NK cells': 'NK',
          'Progenitor cells': 'Progenitor',
          'TRAV1-2- CD8+ T cells': 'CD8T',
          'gd T cells': 'gdT'}

def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out
