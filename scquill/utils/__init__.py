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
)
from .store import (
    store_approximation,
)
from .heuristics import (
    guess_normalisation,
    guess_measurement_type,
    guess_celltype_column,
    guess_celltype_order,
)
