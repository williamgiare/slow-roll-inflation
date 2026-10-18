#========
# WG Edit
# Stable entry points shared by the command, validations and future notebooks.
#========
from importlib.metadata import version as _package_version

from .consistency import consistency
from .spectra import spectra
from .orders import ConsistencyOrder, UnsupportedOrderError

__version__ = _package_version("slow-roll-inflation")
__all__ = ['consistency', 'spectra', 'ConsistencyOrder', 'UnsupportedOrderError']
