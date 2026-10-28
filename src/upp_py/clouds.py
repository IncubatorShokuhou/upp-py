"""Low, mid and high cloud fraction, ``INITPOST.F`` (``COMMCODE`` / NCAR).

The DTC build defines ``COMMCODE``, so ARW does not post the model's
``CLDFRA``. It rebuilds a binary cloud mask (``qc + qi > 1e-7`` kg/kg)
and sets each point's fraction to the fraction of neighbors inside a
30 km box that are cloudy. Low / mid / high are the maximum of that
fraction in the layers

* low: ``p >= 642 hPa``
* mid: ``350 hPa <= p < 642 hPa``
* high: ``150 hPa <= p < 350 hPa``

Column total cloud cover is the max of those three, which is the max
over every level at or below 150 hPa.
"""

from __future__ import annotations

import numpy as np

from upp_py.constants import CLOUD_DEF_P, CLOUD_RADIUS_M, PTOP_HIGH, PTOP_LOW, PTOP_MID
from upp_py.wrf import WrfState


def _box_fraction(mask: np.ndarray, radius: int) -> np.ndarray:
    """Fraction of True neighbors in a square of half-width ``radius``.

    Edges use only the in-domain part of the window, matching the
    ``max(1, i-numr):min(im, i+numr)`` loops in ``INITPOST.F``.
    """
    cloudy = mask.astype(np.float64)
    # Integral image with a zero border so every window sum is four lookups.
    integ = np.pad(cloudy, 1).cumsum(0).cumsum(1)

    def _rect(y0, x0, y1, x1):
        return integ[y1, x1] - integ[y0, x1] - integ[y1, x0] + integ[y0, x0]

    ny, nx = mask.shape
    yy = np.arange(ny)[:, None]
    xx = np.arange(nx)[None, :]
    y0 = np.clip(yy - radius, 0, ny)
    x0 = np.clip(xx - radius, 0, nx)
    y1 = np.clip(yy + radius + 1, 0, ny)
    x1 = np.clip(xx + radius + 1, 0, nx)
    # integ is padded by 1, and the cumsum index equal to the exclusive
    # end in the original array is `end` (because of the pad) ... 
    # integ[k+1] includes original rows 0..k-1. A window [a, b) on the
    # original array sums integ[b+1] - integ[a] after the pad of 1? 
    # pad then cumsum: integ[i+1, j+1] = sum of cloudy[:i, :j].
    # sum over [y0:y1, x0:x1] = integ[y1, x1] - integ[y0, x1] - integ[y1, x0] + integ[y0, x0]
    # where integ index is the exclusive original end. With pad, the
    # cumsum at position (y1) in the padded array (index y1, since pad
    # shifts by 1 and cumsum[y1] includes padded rows 0..y1-1 which are
    # original 0..y1-2 plus the leading zero). Let me just use the
    # padded-cumsum convention carefully:
    #
    # padded[0] = 0, padded[1:] = data
    # integ[r] = sum(padded[:r]) along that axis after both cumsums,
    # so integ[y+1, x+1] = sum(data[:y, :x]).
    # sum(data[y0:y1, x0:x1]) = integ[y1+? 
    count = _rect(y0, x0, y1, x1)
    area = (y1 - y0) * (x1 - x0)
    return count / np.maximum(area, 1.0)


def cloud_fraction(state: WrfState) -> dict[str, np.ndarray]:
    """Low, mid, high and total cloud fraction in percent."""
    dx = state.dx if state.dx > 0 else 3000.0
    numr = int(np.rint(CLOUD_RADIUS_M / dx))
    numr = max(numr, 0)
    low = np.zeros(state.terrain.shape)
    mid = np.zeros_like(low)
    high = np.zeros_like(low)
    for k in range(state.nz):
        mask = (state.qc[k] + state.qi[k]) > CLOUD_DEF_P
        frac = _box_fraction(mask, numr) if numr else mask.astype(np.float64)
        p = state.p_mid[k]
        low = np.where(p >= PTOP_LOW, np.maximum(low, frac), low)
        mid = np.where((p < PTOP_LOW) & (p >= PTOP_MID), np.maximum(mid, frac), mid)
        high = np.where((p < PTOP_MID) & (p >= PTOP_HIGH), np.maximum(high, frac), high)
    total = np.maximum(np.maximum(low, mid), high)
    return {
        "low": np.clip(low, 0.0, 1.0) * 100.0,
        "mid": np.clip(mid, 0.0, 1.0) * 100.0,
        "high": np.clip(high, 0.0, 1.0) * 100.0,
        "total": np.clip(total, 0.0, 1.0) * 100.0,
    }
