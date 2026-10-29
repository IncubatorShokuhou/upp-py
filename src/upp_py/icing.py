"""In-flight icing category from ``AVIATION.f`` ``CALICING``.

``MDL2P.f`` available-field 257 (GRIB pname ``TIPD``) calls ``CALICING``
on each requested isobaric level. The routine's comment says relative
humidity is a fraction, but the caller stores ``SAVRH`` from ``CALRH``
after multiplying by 100, and the threshold in the code is 70. Pass
percent.

Icing is 1 where omega is upward (Pa s-1, negative), temperature is
between 251 K and 273 K inclusive, and RH is at least 70%. Otherwise 0.
The comment also mentions equivalent potential vorticity and supercooled
large drops; the routine does not use them.

NCAR GFIP (``GFIP3.f`` ``icing_algo``, fields 450 and 480) is a separate
column algorithm. In v4.1 ``MDLFLD.f`` the isobaric write of
``icing_gfip`` is commented out. That algorithm is not ported.
"""

from __future__ import annotations

import numpy as np


def in_flight_icing(temperature: np.ndarray, rh_percent: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """Binary icing (1 or 0) from ``CALICING``.

    Parameters
    ----------
    temperature
        K.
    rh_percent
        Relative humidity in percent, as ``MDL2P.f`` passes ``SAVRH``.
    omega
        Vertical velocity in Pa s-1. Negative is ascent.
    """
    t = np.asarray(temperature, dtype=np.float64)
    rh = np.asarray(rh_percent, dtype=np.float64)
    w = np.asarray(omega, dtype=np.float64)
    hit = (w < 0.0) & (t <= 273.0) & (t >= 251.0) & (rh >= 70.0)
    return np.where(hit, 1.0, 0.0)
