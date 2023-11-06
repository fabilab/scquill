import pathlib
import anndata

from .utils import (
    load_config,
    filter_cells,
    normalise_counts,
    correct_annotations,
    approximate_dataset,
    store_approximation,
    guess_measurement_type,
    guess_normalisation,
    )



class Quill:
    def __init__(
        self,
        filename=None,
        adata=None,
        output_filename=None,
        celltype_column=None,
        configutation=None,
        include_neighborhood=True,
        ):

        self.filename = filename
        self.adata = adata
        self.output_filename = output_filename
        self.celltype_column = celltype_column
        self.configuration = configuration
        self.include_neighborhood = include_neighborhood

    def __call__(self):
        self._validate_constructor()
        self.load()
        self.preprocess()
        self.compress()
        self.store()

    def _validate_constructor():
        self.configuration = self.configuration or {}
        self.include_neighborhood = bool(self.include_neighborhood)

    def load(self):
        if (self.adata is None) and (self.filename is None):
            raise ValueError(
                "Either filename or adata must be specified."
            )

        if self.adata is not None:
            return

        self.adata = anndata.read(self.filename)

    def preprocess(self):
        config = self.configuration

        self.adata = filter_cells(self.adata, config)

        self.adata = normalise_counts(
            self.adata,
            config.get("normalisation", guess_normalisation(self.adata)),
            config.get("measurement_type", guess_measurement_type(self.adata)),
        )

        self.adata = correct_annotations(
            self.adata,
            self.celltype_column,
            config,
        )

    def compress(self):
        self.approximation = approximate_dataset(
            self.adata,
        )

    def store(self):
        config = self.configuration

        store_approximation(
            self.output_filename,
            self.approximation,
            config.get("measurement_type", guess_measurement_type(self.adata)),
        )
