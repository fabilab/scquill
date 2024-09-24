import os
import pathlib
import anndata

from .utils import (
    load_config,
    filter_cells,
    normalise_counts,
    correct_annotations,
    approximate_dataset,
    guess_measurement_type,
    guess_normalisation,
    guess_celltype_column,
    guess_celltype_order,
)
from .io import (
    write_to_h5,
)


class Compressor:
    def __init__(
        self,
        filename=None,
        adata=None,
        output_filename=None,
        celltype_column=None,
        celltype_order=None,
        additional_groupby_columns=tuple(),
        configuration=None,
        include_neighborhood=True,
        measurement_type="gene_expression",
    ):

        self.filename = filename
        self.adata = adata
        self.output_filename = output_filename
        self.celltype_column = celltype_column
        self.celltype_order = celltype_order
        self.additional_groupby_columns = additional_groupby_columns
        self.configuration = configuration
        self.include_neighborhood = include_neighborhood
        self.measurement_type = measurement_type

    def __call__(self):
        self.prepare()
        self.compress()
        self.store()

    def prepare(self):
        self._validate_constructor()
        self.load()
        self.preprocess()

    def _validate_constructor(self):
        self.configuration = self.configuration or {}
        self.include_neighborhood = bool(self.include_neighborhood)
        if self.output_filename:
            self.output_filename = pathlib.Path(self.output_filename)

    def _guess_annotation_info_if_needed(self):
        if self.celltype_column is None:
            self.celltype_column = guess_celltype_column(self.adata)

        if self.celltype_order is None:
            self.celltype_order = guess_celltype_order(
                self.adata,
                self.celltype_column,
            )

        if self.additional_groupby_columns is None:
            self.additional_groupby_columns = tuple()

    @property
    def additional_groupby_columns(self):
        return self._additonal_groupby_columns

    @additional_groupby_columns.setter
    def additional_groupby_columns(self, value):
        value = tuple(value)
        if self.celltype_column in value:
            raise ValueError(
                "The cell type column cannot also be an additional column",
            )
        self._additonal_groupby_columns = value

    def load(self):
        if (self.adata is None) and (self.filename is None):
            raise ValueError("Either filename or adata must be specified.")

        if self.adata is not None:
            return

        self.adata = anndata.read_h5ad(self.filename)

    def preprocess(self):
        config = self.configuration

        self._guess_annotation_info_if_needed()

        self.adata = filter_cells(self.adata, config)

        if "normalisation" not in config:
            config["normalisation"] = guess_normalisation(self.adata)
        if "measurement_type" not in config:
            config["measurement_type"] = guess_measurement_type(self.adata)

        self.adata = normalise_counts(
            self.adata,
            config["normalisation"],
            config["measurement_type"],
        )

        # self.adata = correct_annotations(
        #    self.adata,
        #    self.celltype_column,
        #    config,
        # )

    def compress(self):
        self.approximation = approximate_dataset(
            self.adata,
            self.celltype_column,
            self.additional_groupby_columns,
            measurement_type=self.measurement_type,
            include_neighborhood=self.include_neighborhood,
        )

    def store(self):
        if self.output_filename is not None:
            if pathlib.Path(self.output_filename).exists():
                os.remove(self.output_filename)
            write_to_h5(
                self.output_filename,
                self.approximation,
            )

    def to_anndata(self):
        """Export approximation to anndata object, excluding the neighborhood."""
        mt = self.measurement_type
        adata = anndata.AnnData(
            X=self.approximation[mt]["Xave"].T,
            obs=self.approximation[mt]["obs"],
            layers={
                "fraction_detected": self.approximation[mt]["Xfrac"].T,
            },
        )
        adata.obs_names = self.approximation[mt]["obs_names"]
        adata.var_names = self.approximation[mt]["var_names"]
        adata.uns["measurement_type"] = mt
        return adata

    @classmethod
    def from_anndata(
        cls,
        adata,
        output_filename=None,
        measurement_type=None,
    ):
        """Create a compressor from an anndata object."""
        self = cls(
            output_filename=output_filename,
            include_neighborhood=False,
        )

        if measurement_type is None:
            measurement_type = adata.uns.get("measurement_type", "gene_expression")

        self.measurement_type = measurement_type
        self.approximation = {
            measurement_type: {
                "Xave": adata.X.T,
                "Xfrac": adata.layers["fraction_detected"].T,
                "obs": adata.obs,
                "obs_names": adata.obs_names,
                "var_names": adata.var_names,
            }
        }

        return self
