'''
Utility functions for the compression
'''
import os
import gc
import pathlib
import yaml
import gzip
import numpy as np
import pandas as pd
import h5py

import scanpy as sc

from .config import load_config
from .preprocess import (
    filter_cells,
    normalise_counts,
    correct_annotations,
)
from .compress import (
    approximate_dataset,
    store_approximation,
)


def sanitise_gene_names(genes):
    genes_new = []
    for gene in genes:
        gene_new = gene.replace(',', ';')
        gene_new = gene_new.split(' ')[0]
        genes_new.append(gene_new)

    if len(set(genes_new)) != len(set(genes)):
        raise ValueError("Gene names are not unique after sanitisation.")

    return genes_new


# TODO
def guess_normalisation(adata):
    return None


# TODO
def guess_measurement_type(adata):
    return None
