"""Read from files."""

import h5py
import hdf5plugin
import numpy as np
import pandas as pd
import anndata

from scquill.utils.types import _infer_dtype


def read_h5_to_anndata(
    h5_data,
    neighborhood,
    measurement_type,
):
    """Get an AnnData object in which each observation is an average."""
    if measurement_type not in h5_data["measurements"]:
        raise KeyError("Measurement type not found: {measurement_type}")

    me = h5_data["measurements"][measurement_type]
    compression = me.attrs.get("compression", None)

    if compression:
        try:
            import hdf5plugin
        except ImportError:
            raise ImportError(
                'You need the "hdf5plugin" package to decompress this approximation. You can install it e.g. via pip install hdf5plugin.'
            )

    var_names = me["var_names"].asstr()[:]

    # NOTE: there are a few older versions around
    if "groupby" in me:
        groupby_names = []
        groupby_dtypes = []
        n_levels = me["groupby"].attrs["n_levels"]
        for i in range(n_levels):
            groupby_names.append(me["groupby"]["names"].attrs[str(i)])
            groupby_dtypes.append(me["groupby"]["dtypes"].attrs[str(i)])
        groupby = "\t".join(groupby_names)
    else:
        gbykey = "grouped_by"
        keys = list(me[gbykey].keys())
        if len(keys) != 1:
            raise ValueError("Expected exactly one groupby key, found: {keys}")
        groupby = keys[0]
        groupby_names = me[gbykey][groupby]["names"].asstr()[:]
        groupby_dtypes = me[gbykey][groupby]["dtypes"].asstr()[:]
        n_levels = len(groupby_names)

    if neighborhood:
        neigroup = me["neighborhood"]
        Xave = neigroup["average"][:]
        if "quantisation" in me:
            quantisation = me["quantisation"][:]
            Xave = quantisation[Xave]

        groupby_order = me["obs_names"].asstr()[:]
        obs_names = neigroup["obs_names"].asstr()[:]
        ncells = neigroup["cell_count"][:]
        coords_centroid = neigroup["coords_centroid"][:]
        convex_hulls = []
        for ih in range(len(coords_centroid)):
            convex_hulls.append(neigroup["convex_hull"][str(ih)][:])

        if measurement_type == "gene_expression":
            Xfrac = neigroup["fraction"][:]
            adata = anndata.AnnData(
                X=Xave,
                layers={
                    "average": Xave,
                    "fraction": Xfrac,
                },
            )
        else:
            adata = anndata.AnnData(X=Xave)

        adata.obsm["X_ncells"] = ncells
        adata.obsm["X_umap"] = coords_centroid
        adata.uns["convex_hulls"] = convex_hulls

        adata.obs_names = pd.Index(obs_names, name="neighborhoods")
        adata.var_names = pd.Index(var_names, name="features")

    else:
        # NOTE: this is the historical way to store atlas approximations
        # when the organism has multiple tissues
        if "->" in groupby:
            read_fun = _read_data_from_stratified_group
        else:
            read_fun = _read_data_from_flat_group
        resdict = read_fun(me, groupby, groupby_names, groupby_dtypes, measurement_type)

        if measurement_type == "gene_expression":
            adata = anndata.AnnData(
                X=resdict["Xave"],
                obs=resdict["obs"],
                layers={
                    "average": resdict["Xave"],
                    "fraction": resdict["Xfrac"],
                },
            )
        else:
            adata = anndata.AnnData(
                X=resdict["Xave"],
                obs=resdict["obs"],
            )

        adata.var_names = pd.Index(var_names, name="features")

    adata.uns["approximation_groupby"] = {
        "names": groupby_names,
        "dtypes": groupby_dtypes,
    }
    if neighborhood:
        adata.uns["approximation_groupby"]["order"] = groupby_order
        adata.uns["approximation_groupby"]["cell_count"] = me["cell_count"][:]

    return adata


def _read_data_from_flat_group(
    me, groupby, groupby_names, groupby_dtypes, measurement_type
):
    # Data
    Xave = me["average"][:]
    if "quantisation" in me:
        quantisation = me["quantisation"][:]
        Xave = quantisation[Xave]
    if measurement_type == "gene_expression":
        Xfrac = me["fraction"][:]

    # Obs metadata
    obs_names = me["obs_names"].asstr()[:]
    obs = pd.DataFrame([], index=obs_names)
    for column, dtype in zip(groupby_names, groupby_dtypes):
        if _infer_dtype(dtype) == "S":
            obs[column] = me["obs"][column].asstr()[:]
        else:
            obs[column] = me["obs"][column][:]
    obs["cell_count"] = me["cell_count"][:]
    obs.index = pd.Index(obs.index, name=groupby)

    resdict = {
        "Xave": Xave,
        "obs": obs,
    }
    if measurement_type == "gene_expression":
        resdict["Xfrac"] = Xfrac
    return resdict


def _read_data_from_stratified_group(
    me, groupby, groupby_names, groupby_dtypes, measurement_type
):
    """Read data from a stratified group.

    NOTE: This is the historical way to store atlas approximations. It groups data by tissue (to have distinct neighborhoods)
    and then the obs_names inside each subgroup are the cell types. This function is therefore quite ad-hoc but it's ok for now.
    """
    # Iterate over the higher-level groups (tissues)
    Xave = []
    if measurement_type == "gene_expression":
        Xfrac = []
    obs = []

    # E.g. "tissues"
    groupby = "->".join(groupby_names)
    group = me["data"][groupby]
    for subgroupname, subgroup in group.items():
        # Data
        subXave = subgroup["average"][:]
        if "quantisation" in me:
            quantisation = me["quantisation"][:]
            subXave = quantisation[subXave]
        Xave.append(subXave)
        if measurement_type == "gene_expression":
            Xfrac.append(subgroup["fraction"][:])

        # Obs metadata (notice the tricky obs_names indexing)
        subobs_names = subgroup["obs_names"].asstr()[:]
        subobs = pd.DataFrame([], index=subobs_names)
        # Set the common column (typically tissue)
        subobs[groupby_names[0]] = subgroupname
        # Set the discriminative column (typicall cell type)
        subobs[groupby_names[1]] = subobs_names
        # NOTE: assuming there are not other columns, i.e. something like tissue->celltype\tage is not supported
        subobs["cell_count"] = subgroup["cell_count"][:]
        subobs.index = subobs[groupby_names].apply("->".join, axis=1)
        obs.append(subobs)

    # Concatenate the lists
    Xave = np.concatenate(Xave)
    obs = pd.concat(obs)
    obs.index = pd.Index(obs.index, name=groupby)
    resdict = {
        "Xave": Xave,
        "obs": obs,
    }

    if measurement_type == "gene_expression":
        Xfrac = np.concatenate(Xfrac)
        resdict["Xfrac"] = Xfrac

    return resdict
