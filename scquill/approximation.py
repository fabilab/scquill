import h5py
import pandas as pd
import anndata

from .utils.coarse_grain import coarse_grain_anndata
from .utils.types import _infer_dtype


class Approximation:
    """Access single cell approximations."""

    def __init__(
        self,
    ):
        self._adata_dict = {}

    @classmethod
    def read_h5(
        cls,
        filename,
    ):
        """Lazy reader of approximation from file."""
        self = cls()
        self.approximation_dict = None
        self.filename = filename
        return self

    @classmethod
    def read_approximation_dict(
        cls,
        approximation_dict,
    ):
        """Lazy reader of approximation from dict."""
        self = cls()
        self.filename = None
        self.approximation_dict = approximation_dict
        return self

    def _infer_measurement_type(self, h5_data):
        measurement_types = list(h5_data["measurements"])
        if len(measurement_types) == 1:
            return measurement_types[0]
        elif len(measurement_types) == 0:
            raise KeyError("No measurements found in this approximation")
        else:
            raise KeyError("Multiple measurement types found: {measurement_types}")

    def _to_anndata(
        self,
        neighborhood=False,
        measurement_type=None,
    ):
        if self.approximation_dict:
            adata = self._appdict_to_anndata(
                neighborhood=neighborhood,
                measurement_type=measurement_type,
            )
        else:
            adata = self._h5file_to_anndata(
                neighborhood=neighborhood,
                measurement_type=measurement_type,
            )
        self._adata_dict[(measurement_type, neighborhood)] = adata

    def _h5file_to_anndata(
        self,
        neighborhood,
        measurement_type,
    ):
        """Get an AnnData object in which each observation is an average."""
        with h5py.File(self.filename) as h5_data:
            if measurement_type is None:
                measurement_type = self._infer_measurement_type(h5_data)

            if measurement_type not in h5_data["measurements"]:
                raise KeyError("Measurement type not found: {measurement_type}")

            me = h5_data["measurements"][measurement_type]
            compression = me.attrs["compression"]

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
                # TODO: quantisation

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

                # TODO: Make a multiindex if needed
                adata.uns["convex_hulls"] = convex_hulls

                adata.obs["cell_count"] = ncells
                adata.obs_names = pd.Index(obs_names, name="neighborhoods")
                adata.var_names = pd.Index(var_names, name="features")

            else:
                Xave = me["average"][:]
                # TODO: quantisation

                if measurement_type == "gene_expression":
                    Xfrac = me["fraction"][:]
                obs_names = me["obs_names"].asstr()[:]
                # Add obs metadata
                ncells = me["cell_count"][:]
                obs = pd.DataFrame([], index=obs_names)
                for column, dtype in zip(groupby_names, groupby_dtypes):
                    if _infer_dtype(dtype) == "S":
                        obs[column] = me["obs"][column].asstr()[:]
                    else:
                        obs[column] = me["obs"][column][:]
                obs["cell_count"] = ncells

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

        return adata

    def _appdict_to_anndata(
        neighborhood,
        measurement_type,
    ):
        compressed_atlas = self.approximation_dict

        if measurement_type is None:
            if len(compressed_atlas.keys()) == 1:
                measurement_type = list(compressed_atlas.keys())[0]
            else:
                raise ValueError(
                    "Multiple measurement types detected, which one would you like to look at?",
                )

        resd = compressed_atlas[measurement_type]
        var_names = resd["features"]

        if neighborhood:
            neid = resd["neighborhood"]
            Xave = neid["avg"].values
            obs_names = neid["avg"].index.values

            if measurement_type == "gene_expression":
                Xfrac = neid["frac"].values
                adata = anndata.AnnData(
                    X=Xave,
                    layers={
                        "average": Xave,
                        "fraction": Xfrac,
                    },
                )
            else:
                adata = anndata.AnnData(X=Xave)

            adata.obs["cell_count"] = neid["ncells"]
            adata.obsm["X_ncells"] = neid["ncells"]
            adata.obsm["X_umap"] = neid["coords_centroid"]
            adata.uns["convex_hulls"] = neid["convex_hull"]
            adata.obs_names = pd.Index(obs_names, name="neighborhoods")
            adata.var_names = pd.Index(var_names, name="features")

        else:
            Xave = resd["avg"].values
            if measurement_type == "gene_expression":
                Xfrac = resd["frac"].values
                adata = anndata.AnnData(
                    X=Xave,
                    layers={
                        "average": Xave,
                        "fraction": Xfrac,
                    },
                )
            else:
                adata = anndata.AnnData(
                    X=Xave,
                )

            groupby_names = resd["avg"].index.names
            groupby = "\t".join(groupby_names)
            obs_names = resd["avg"].index.map(lambda x: "\t".join(str(y) for y in x))
            obs = resd["obs"].copy()
            obs["cell_count"] = resd["ncells"]
            obs.index = obs_names
            adata.obs = obs
            adata.var_names = pd.Index(var_names, name="features")
            adata.obs_names = pd.Index(obs_names, name=groupby)

        adata.uns["approximation_groupby"] = {
            "names": groupby_names,
            "dtypes": groupby_dtypes,
        }
        if neighborhood:
            adata.uns["approximation_groupby"]["order"] = resd["avg"].index.values

        return adata

    def to_anndata(
        self,
        groupby="celltype",
        neighborhood=False,
        measurement_type=None,
    ):
        if measurement_type is None:
            with h5py.File(self.filename) as h5_data:
                measurement_type = self._infer_measurement_type(h5_data)

        if isinstance(groupby, str):
            groupby = [groupby]
        groupby = tuple(groupby)

        if (measurement_type, neighborhood) not in self._adata_dict:
            self._to_anndata(
                neighborhood=neighborhood,
                measurement_type=measurement_type,
            )

        # FIXME: specify that it's a view somehow
        adata = self._adata_dict[(measurement_type, neighborhood)]

        # Coarse grain result if the approximation includes unnecessary metadata
        adata_cg = coarse_grain_anndata(adata, groupby)
        return adata_cg
