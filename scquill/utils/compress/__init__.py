"""
Utility functions for the compression
"""

import gc
import numpy as np
import pandas as pd
import scanpy as sc


from .biological import _compress_biological
from .neighborhood import _compress_neighborhood


def approximate_dataset(
    adata,
    celltype_column,
    additional_groupby_columns=tuple(),
    measurement_type="gene_expression",
    include_neighborhood=True,
):
    """Compress atlas for one tissue after data is clean, normalised, and reannotated."""

    # Celltype averages
    biod = _compress_biological(
        adata,
        celltype_column,
        additional_groupby_columns,
        measurement_type,
    )
    features = adata.var_names
    result = {
        "var_names": features,
        "obs": biod["obs"],
        "obs_names": biod["obs_names"],
        "Xave": biod["Xave"],
    }
    if measurement_type == "gene_expression":
        result["Xfrac"] = biod["Xfrac"]

    if include_neighborhood:
        # Local neighborhoods
        neid = _compress_neighborhood(
            adata,
            biod["obs"],
            celltype_column,
            additional_groupby_columns,
            measurement_type,
        )
        result["neighborhood"] = neid

    compressed_atlas = {
        measurement_type: result,
    }

    return compressed_atlas
