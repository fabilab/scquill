import importlib.metadata

__version__ = importlib.metadata.version("scquill")

from .compressor import Compressor
from .accessor import Accessor
