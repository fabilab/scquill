"""Compress single cell data using biological cell metadata groups, including cell type."""

import gc
import numpy as np
import pandas as pd
import scanpy as sc


def _compress_biological(
    adata,
    celltype_column,
    additional_groupby_columns,
    measurement_type,
):
    """Compress at the cell type level"""
    features = adata.var_names
    nfeatures = len(features)

    # Metadata for groups
    groupby_columns = [celltype_column] + list(additional_groupby_columns)
    tmp = adata.obs[groupby_columns].copy()
    tmp["cell_count"] = 1
    # NOTE: "observbed" seems to be a little iffy incl. for the pandas folks who are defaulting it to True and deprecating it in the future
    # We recently changed this from False to True (17 July 2024)
    obs = tmp.groupby(groupby_columns, observed=True).sum().reset_index()
    ngroups = len(obs)

    # Names of rows
    obs_names = [
        "\t".join(str(y) for y in x) for _, x in obs[groupby_columns].iterrows()
    ]
    obs_names = pd.Index(
        obs_names,
        name="\t".join(groupby_columns),
    )

    # Molecular counts
    avg = np.zeros((nfeatures, ngroups), np.float32)
    if measurement_type == "gene_expression":
        frac = np.zeros_like(avg)
    idx = np.zeros(adata.n_obs, bool)
    for i, row in obs.iterrows():
        # Reconstruct indices of focal cells
        idx[:] = True
        for col in groupby_columns:
            idx &= adata.obs[col] == row[col]

        # Average across group
        Xidx = adata[idx].X
        avg[:, i] = np.asarray(Xidx.mean(axis=0))[0]
        if measurement_type == "gene_expression":
            frac[:, i] = np.asarray((Xidx > 0).mean(axis=0))[0]

    res = {
        "kind": "biological",
        "obs": obs,
        "obs_names": obs_names,
        "Xave": avg,
    }
    if measurement_type == "gene_expression":
        res["Xfrac"] = frac
    return res
