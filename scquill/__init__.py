import importlib.metadata

__version__ = importlib.metadata.version("scquill")

from .compressor import Compressor
from .approximation import ApproxAnnData

import scquill.pl as pl


__all__ = (
    "__version__",
    "Compressor",
    "ApproxAnnData",
    "pl",
)
