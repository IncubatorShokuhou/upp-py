"""Aerosol fields read by ``INITPOST.F`` and scaled in ``MDLFLD.f``.

What v4.1 actually reads from a WRF-ARW file:

* ``QNCLOUD``, ``QNWFA`` and ``QNIFA`` only when ``imp_physics == 28``
  (aerosol-aware Thompson). ``MDLFLD.f`` writes them on model levels and
  sets values below ``1e-8`` to 0.
* ``smoke`` (WRF-Chem). The hybrid-level write is
  ``smoke * PMID / (RD * T)``, commented as ug kg-1 to ug m-3.
* ``TAOD5502D``, ``AERASY2D`` and ``AERSSA2D`` are copied through.
* ``DUST_1`` ... ``DUST_5`` are inside ``if(1==2)`` in ``INITPOST.F``, so
  ``unipost`` v4.1 never fills them. ``MDLFLD.f`` would multiply each bin
  by density the same way. This module applies that conversion when the
  arrays are in the file.

Sea salt, sulfate, organic carbon and soot appear in ``MDLFLD.f`` but
have no ``VarName`` read in the WRF ``INITPOST.F`` path. They are not
produced from wrfout.

GSD visibility can add ``aextc55`` from ``AOD3D_SMOKE / dz``. That
visibility routine is not ported; the 3-D smoke array is still available
as a mass concentration when the file contains ``smoke``.
"""

from __future__ import annotations

import numpy as np

from upp_py.constants import RD

_NUMBER_FLOOR = 1.0e-8


def clip_number_concentration(number: np.ndarray) -> np.ndarray:
    """``MDLFLD.f`` floor: values below ``1e-8`` are stored as 0."""
    arr = np.asarray(number, dtype=np.float64)
    return np.where(arr < _NUMBER_FLOOR, 0.0, arr)


def mass_concentration(mixing_ratio: np.ndarray, pressure: np.ndarray, temperature: np.ndarray) -> np.ndarray:
    """ug m-3 from ug kg-1, using ``PMID / (RD * T)`` as in the smoke write.

    Parameters
    ----------
    mixing_ratio
        ug kg-1, model levels, WRF bottom-up.
    pressure, temperature
        Pa and K on the same grid.
    """
    rho = np.asarray(pressure, dtype=np.float64) / (RD * np.asarray(temperature, dtype=np.float64))
    return np.asarray(mixing_ratio, dtype=np.float64) * rho
