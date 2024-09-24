"""Read from files."""

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

    if "quantisation" in me:
        quantisation = me["quantisation"][:]

    groupby_names = []
    groupby_dtypes = []
    n_levels = me["groupby"].attrs["n_levels"]
    for i in range(n_levels):
        groupby_names.append(me["groupby"]["names"].attrs[str(i)])
        groupby_dtypes.append(me["groupby"]["dtypes"].attrs[str(i)])
    groupby = "\t".join(groupby_names)

    if neighborhood:
        neigroup = me["neighborhood"]
        Xave = neigroup["average"][:]
        if "quantisation" in me:
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
        Xave = me["average"][:]
        if "quantisation" in me:
            Xave = quantisation[Xave]

        if measurement_type == "gene_expression":
            Xfrac = me["fraction"][:]
        obs_names = me["obs_names"].asstr()[:]
        # Add obs metadata
        obs = pd.DataFrame([], index=obs_names)
        for column, dtype in zip(groupby_names, groupby_dtypes):
            if _infer_dtype(dtype) == "S":
                obs[column] = me["obs"][column].asstr()[:]
            else:
                obs[column] = me["obs"][column][:]
        obs["cell_count"] = me["cell_count"][:]

        if measurement_type == "gene_expression":
            adata = anndata.AnnData(
                X=Xave,
                obs=obs,
                layers={
                    "average": Xave,
                    "fraction": Xfrac,
                },
            )
        else:
            adata = anndata.AnnData(
                X=Xave,
                obs=obs,
            )

        adata.var_names = pd.Index(var_names, name="features")
        adata.obs_names = pd.Index(obs_names, name=groupby)

    adata.uns["approximation_groupby"] = {
        "names": groupby_names,
        "dtypes": groupby_dtypes,
    }
    if neighborhood:
        adata.uns["approximation_groupby"]["order"] = groupby_order
        adata.uns["approximation_groupby"]["cell_count"] = me["cell_count"][:]

    return adata
