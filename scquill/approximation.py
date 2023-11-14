import h5py
import pandas as pd
import anndata


class Approximation:
    """Access single cell approximations."""
    def __init__(
        self,
        filename,
    ):
        self.filename = filename
        self._adata_dict = {}

    def _infer_measurement_type(self, h5_data):
        measurement_types = list(h5_data['measurements'])
        if len(measurement_types) == 1:
            return measurement_types[0]
        elif len(measurement_types) == 0:
            raise KeyError("No measurements found in this approximation")
        else:
            raise KeyError(
                "Multiple measurement types found: {measurement_types}"
            )

    def _to_adata(
        self,
        neighborhood=False,
        measurement_type=None,
    ):
        """Get an AnnData object in which each observation is an average."""
        with h5py.File(self.filename) as h5_data:
            if measurement_type is None:
                measurement_type = self._infer_measurement_type(h5_data)

            if measurement_type not in h5_data['measurements']:
                raise KeyError(
                    "Measurement type not found: {measurement_type}"
                )

            me = h5_data['measurements'][measurement_type]
            compression = me.attrs['compression']

            if compression:
                try:
                    import hdf5plugin
                except ImportError:
                    raise ImportError(
                        "You need the \"hdf5plugin\" package to decompress this approximation. You can install it e.g. via pip install hdf5plugin."
                    )

            var_names = me['features'].asstr()[:]

            if 'quantisation' in me:
                quantisation = me['quantisation'][:]

            groupby_names = []
            groupby_dtypes = []
            n_levels = me['groupby'].attrs['n_levels']
            for i in range(n_levels):
                groupby_names.append(me['groupby']['names'].attrs[str(i)])
                groupby_dtypes.append(me['groupby']['dtypes'].attrs[str(i)])
            groupby = '>'.join(groupby_names)

            if neighborhood:
                neigroup = me['neighborhood']
                Xave = neigroup['average'][:]
                # TODO: quantisation

                if measurement_type == "gene_expression":
                    Xfrac = neigroup['fraction'][:]

                groupby_order = me['index'].asstr()[:]
                obs_names = neigroup['index'].asstr()[:]
                ncells = neigroup['cell_count'][:]
                coords_centroid = neigroup['coords_centroid'][:]
                convex_hulls = []
                for ih in range(len(coords_centroid)):
                    convex_hulls.append(
                        neigroup['convex_hull'][str(ih)][:]
                    )

                if measurement_type == "gene_expression":
                    adata = anndata.AnnData(
                        X=Xave,
                        layers={
                            'average': Xave,
                            'fraction': Xfrac,
                        }
                    )
                else:
                    adata = anndata.AnnData(X=Xave)

                adata.obs_names = pd.Index(obs_names, name='neighborhoods')
                adata.var_names = pd.Index(var_names, name='features')
                adata.obsm['X_ncells'] = ncells
                adata.obsm['X_umap'] = coords_centroid

                # TODO: Make a multiindex if needed
                adata.uns['convex_hulls'] = convex_hulls

            else:
                Xave = me['average'][:]
                # TODO: quantisation

                if measurement_type == "gene_expression":
                    Xfrac = me['fraction'][:]
                obs_names = me['index'].asstr()[:]
                ncells = me['cell_count'][:]

                if measurement_type == "gene_expression":
                    adata = anndata.AnnData(
                        X=Xave,
                        layers={
                            'average': Xave,
                            'fraction': Xfrac,
                        }
                    )
                else:
                    adata = anndata.AnnData(X=Xave)

                adata.var_names = pd.Index(var_names, name='features')
                adata.obs_names = pd.Index(obs_names, name=groupby)
                # Add obs metadata
                adata.obs['cell_count'] = ncells
                multiindex = adata.obs_names.str.split('\t', expand=True).to_frame()
                multiindex.columns = groupby_names
                for name in groupby_names:
                    adata.obs[name] = multiindex[name]

        adata.uns['approximation_groupby'] = {
            'names': groupby_names,
            'dtypes': groupby_dtypes,
        }
        if neighborhood:
            adata.uns['approximation_groupby']['order'] = groupby_order

        self._adata_dict[(measurement_type, groupby, neighborhood)] = adata

    def to_adata(
        self,
        groupby='celltype',
        neighborhood=False,
        measurement_type=None,
        ):

        if measurement_type is None:
            with h5py.File(self.filename) as h5_data:
                measurement_type = self._infer_measurement_type(h5_data)

        if (measurement_type, groupby, neighborhood) not in self._adata_dict:
            self._to_adata(
                neighborhood=neighborhood,
                measurement_type=measurement_type,
            )

        # FIXME: specify that it's a view somehow
        adata = self._adata_dict[(measurement_type, groupby, neighborhood)]

        # Coarse grain result if the approximation includes unnecessary metadata
        adata_cg = coarse_grain(adata, groupby)
        return adata_cg


def coarse_grain(adata, groupby):
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
        ncells_cg = ncells_with_meta.groupby(indices_groupby).sum()
        new_order = np.asarray(
            ['\t'.join(x) for x in ncells_cg.index],
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
            X[i] = (X[idx] * ncells_with_meta.iloc[idx]).sum(axis=0) / ncells_cg.iloc[i]
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
