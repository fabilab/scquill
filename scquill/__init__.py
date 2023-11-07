import importlib.metadata

__version__ = importlib.metadata.version("scquill")

from .compressor import Compressor
from .accessor import Accessor

import scquill.plot as pl


__all__ = (
    "__version__",
    "Compressor",
    "Accessor",
    "pl",
)
