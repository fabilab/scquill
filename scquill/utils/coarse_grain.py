import numpy as np

# FIXME FIXME
def coarse_grain_anndata(adata, groupby):
    """Coarse grain an approximation further to broader grouping."""

    if 'approximation_groupby' not in adata.uns:
        raise KeyError(
            "missing .uns['approximation_groupby'] information, this does not look like an approximation",
        )

    if isinstance(groupby, str):
        groupby = [groupby]

    groupby = list(groupby)
    if len(groupby) == 0:
        raise ValueError(
            'groupby must be a sequence with at least one element',
        )

    groupby_original = list(adata.uns['approximation_groupby']['names'])
    for column in groupby:
        if column not in groupby_original:
            raise ValueError(
                f"Grouping not found in approximation: {column}",
            )

    # Trivial case
    if len(groupby_original) == 1:
        return adata.copy()

    indices_groupby = [groupby_original.index(gb) for gb in groupby]

    # Detect neighborhood
    neighborhood = 'X_ncells' in adata.obsm
    if neighborhood:
        # We only have to touch the X_ncells, the rest is grouping agnostic
        multiindex = (pd.Index(adata.uns['approximation_groupby']['order'])
                        .str
                        .split('\t', expand=True))
        multiindex.names = groupby_names
        ncells = pd.DataFrame(
            adata.obsm['X_ncells'],
            index=multiindex,
        )
        # FIXME: this probably needs more work
        ncells_cg = ncells.groupby(level=indices_groupby).sum()
        new_order = np.asarray(
            ['\t'.join(x) for x in ncells_cg.index],
        )

        adata_cg = adata.copy()
        adata_cg.obsm['X_ncells'] = ncells_cs.values
        adata_cg.uns['approximation_groupby'] = {
            'names': groupby,
            'dtypes': [adata.uns['approximation_groupby']['dtypes'][i] for i in indices_groupby],
            'order': new_order,
        }

    else:
        # If not neighborhood, we have to actually change the matrix X and, if present,
        # the layers (e.g. avg and fractions)
        ncells_with_meta = adata.obs[list(groupby_original) + ['cell_count']].copy()
        # FIXME: this probably needs more work
        ncells_cg = ncells_with_meta.groupby(groupby).sum('cell_count')
        new_order = np.asarray(
            ['\t'.join(str(x)) for x in ncells_cg.index],
        )

        X = np.zeros((len(new_order), adata.X.shape[1]), adata.X.dtype)
        if adata.layers:
            layers = {
                key: np.zeros((len(new_order), layer.shape[1]), layer.dtype) for key, layer in adata.layers.items()
            }
        for i in range(len(new_order)):
            groupid = ncells_cg.index[i]
            idx = ncells_with_meta[groupby[0]] == groupid[0]
            if len(groupid) > 1:
                for j in range(1, len(groupid)):
                    idx &= ncells_with_meta[groupby[j]] == groupid[j]
            idx = idx.values.nonzero()[0]
            X[i] = (X[idx] * ncells_with_meta.iloc[idx].values).sum(axis=0) / ncells_cg.iloc[i]
            if adata.layers:
                for key in layers:
                    Xl = adata.layers[key]
                    layers[key][i] = (Xl[idx] * ncells_with_meta.iloc[idx]).sum(axis=0) / ncells_cg.iloc[i]

        if adata.layers:
            adata_cg = anndata.AnnData(
                X=X,
                var=adata.var.copy(),
            )
        else:
            adata_cg = anndata.AnnData(
                X=X,
                layers=layers,
                var=adata.var.copy(),
            )

        # Obs names and metadata
        adata_cg.obs_names = new_order
        adata_cg.obs['cell_count'] = ncells_cg.values
        for gb in groupby:
            adata_cg.obs[gb] = ncells_cg.index.get_level_values(gb)

        adata_cg.uns['approximation_groupby'] = {
            'names': groupby,
            'dtypes': [adata.uns['approximation_groupby']['dtypes'][i] for i in indices_groupby],
        }

    return adata_cg

