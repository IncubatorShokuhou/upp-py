"""Absolute vorticity on an A-grid, ``CALVOR.f``.

For ``gridtype='A'`` (WRF-ARW mass points) and any model other than
RAP:

.. math::

    \\zeta_a = \\frac{\\partial v}{\\partial x} - \\frac{\\partial u}{\\partial y}
    + f + \\frac{u \\tan\\phi}{a}

with centered differences, ``\\Delta x = \\Delta x_{nom}/m`` (``DX = dxval/MSFT``
in ``INITPOST.F``) and ``f = 2\\Omega\\sin\\phi``. The one-point rim is
copied from the adjacent interior point; Fortran leaves it at the
zero it was initialized to.
"""

from __future__ import annotations

import numpy as np

from upp_py.constants import ERAD, OMEGA_EARTH


def absolute_vorticity(
    u: np.ndarray,
    v: np.ndarray,
    lat_deg: np.ndarray,
    dx_nominal: float,
    dy_nominal: float,
    mapfac: np.ndarray,
) -> np.ndarray:
    """Absolute vorticity (1/s) from grid-relative winds.

    ``u`` and ``v`` are 2-D (ny, nx) or 3-D (nz, ny, nx). ``lat_deg`` and
    ``mapfac`` are 2-D.
    """
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    squeeze = u.ndim == 2
    if squeeze:
        u = u[None]
        v = v[None]
    dx = dx_nominal / np.maximum(mapfac, 1.0e-6)
    dy = dy_nominal / np.maximum(mapfac, 1.0e-6)
    dvdx = (v[:, :, 2:] - v[:, :, :-2]) / (2.0 * dx[None, :, 1:-1])
    dudy = (u[:, 2:, :] - u[:, :-2, :]) / (2.0 * dy[None, 1:-1, :])
    # Align to the interior (1:-1, 1:-1).
    dvdx_i = dvdx[:, 1:-1, :]
    dudy_i = dudy[:, :, 1:-1]
    lat = np.deg2rad(lat_deg)
    f = 2.0 * OMEGA_EARTH * np.sin(lat)
    metric = u[:, 1:-1, 1:-1] * np.tan(lat[1:-1, 1:-1]) / ERAD
    zeta = dvdx_i - dudy_i + f[None, 1:-1, 1:-1] + metric
    out = np.empty_like(u)
    out[:] = np.nan
    out[:, 1:-1, 1:-1] = zeta
    out[:, 0, :] = out[:, 1, :]
    out[:, -1, :] = out[:, -2, :]
    out[:, :, 0] = out[:, :, 1]
    out[:, :, -1] = out[:, :, -2]
    return out[0] if squeeze else out
