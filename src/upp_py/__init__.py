"""Pure-Python post-processor for WRF-ARW, following UPP v4.1 algorithms.

The numeric core is NumPy. File IO uses xarray. Formulas, constants and
the order of operations follow EMC_post tag ``dtc_post_v4.1.0``
(``sorc/ncep_post.fd``). Output is not bit-for-bit with ``unipost.exe``:
GRIB packing, lookup-table quantization and a few documented
approximations remain.
"""

from upp_py.config import PostConfig, load_config
from upp_py.grib2 import write_grib2
from upp_py.post import dataset_from_state, post_wrfout, write_netcdf

__all__ = [
    "PostConfig",
    "dataset_from_state",
    "load_config",
    "post_wrfout",
    "write_grib2",
    "write_netcdf",
]
__version__ = "0.2.0"
