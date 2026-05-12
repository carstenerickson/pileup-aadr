"""pileup-aadr — one-shot user-BAM → AADR-site pseudohaploid genotypes."""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("pileup-aadr")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
