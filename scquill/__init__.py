import importlib.metadata

__version__ = importlib.metadata.version("scquill")

from .compressor import Compressor
from .approximation import Approximation

import scquill.pl as pl


__all__ = (
    "__version__",
    "Compressor",
    "Approximation",
    "pl",
)
