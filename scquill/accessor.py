import h5py
import pandas as pd
import anndata


class Accessor:
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
        groupby='celltype',
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

            group = me[groupby]

            if neighborhood:
                neigroup = group['neighborhood']
                Xave = neigroup['average'][:]
                # TODO: quantisation

                if measurement_type == "gene_expression":
                    Xfrac = neigroup['fraction'][:]
                celltypes = neigroup['index'].asstr()[:]
                ncells = neigroup['cell_count'][:]
                coords_centroid = neigroup['coords_centroid'][:]
                obs_names = [f'neighborhood_{i+1}' for i in range(len(coords_centroid))]
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
                adata.uns['celltypes'] = pd.Index(celltypes, name='celltypes')
                adata.uns['convex_hulls'] = convex_hulls

            else:
                Xave = group['average'][:]
                # TODO: quantisation

                if measurement_type == "gene_expression":
                    Xfrac = group['fraction'][:]
                obs_names = group['index'].asstr()[:]
                ncells = group['cell_count'][:]

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

                adata.obs_names = pd.Index(obs_names, name='celltypes')
                adata.var_names = pd.Index(var_names, name='features')
                adata.obs['cell_count'] = ncells

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
                groupby=groupby,
                neighborhood=neighborhood,
                measurement_type=measurement_type,
            )

        # FIXME: specify that it's a view somehow
        return self._adata_dict[(measurement_type, groupby, neighborhood)]
