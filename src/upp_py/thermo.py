"""Moisture diagnostics following UPP ``CALRH.f``, ``CALDWP.f``, ``DEWPOINT.f``.

Relative humidity uses the legacy saturation specific humidity

.. math::

    q_s = \\frac{P_{q0}}{p}\\exp\\left(A_2\\frac{T-A_3}{T-A_4}\\right)

with ``PQ0, A2, A3, A4`` from ``params.F``. Dewpoint inverts the same
vapor-pressure formula as ``DEWPOINT.f`` (Tuccillo 1990), including the
2000-entry table and linear interpolation in vapor pressure.
"""

from __future__ import annotations

import numpy as np

from upp_py.constants import (
    A2,
    A3,
    A4,
    CAPA,
    D001,
    EPS_MASS,
    H1M12,
    ONEPS,
    P1000,
    PQ0,
    RHMIN,
)

# DEWPOINT.f table (built once).
_NT = 2000
_RVP1 = 0.0001  # centibar
_RVP2 = 10.0
_RT3 = 273.16
_RVP3 = 0.611
_RA = 5.0065
_RB = 19.83923
_RAPB = _RA + _RB
_RTEST = 1.0e-6
_RDVP = (_RVP2 - _RVP1) / (_NT - 1.0)


def _build_dewpoint_table() -> np.ndarray:
    """Newton inversion of ``DEWPOINT.f`` for each table vapor pressure."""
    rlog3 = np.log(_RVP3)
    rgs = 1.0
    tdp = np.empty(_NT, dtype=np.float64)
    rvp = _RVP1 - _RDVP
    for nn in range(_NT):
        rvp = rvp + _RDVP
        rlvp = np.log(rvp) - rlog3 - _RAPB
        for _ in range(50):
            rn = _RA * np.log(rgs) - _RAPB * rgs - rlvp
            rd = (_RA / rgs) - _RAPB
            rch = rn / rd
            if abs(rch) < _RTEST:
                break
            rgs = rgs - rch
        tdp[nn] = _RT3 / rgs
    return tdp


_TDP = _build_dewpoint_table()
_A_TAB = 1.0 / _RDVP
_B_TAB = 1.0 - _A_TAB * _RVP1
_DNTM1 = float(_NT) - 0.01


def saturation_specific_humidity(temperature: np.ndarray, pressure: np.ndarray) -> np.ndarray:
    """Saturation specific humidity (kg/kg) from ``CALRH.f``.

    Parameters
    ----------
    temperature
        Air temperature (K).
    pressure
        Pressure (Pa).

    Returns
    -------
    numpy.ndarray
        ``q_s`` in kg/kg. This is the quantity UPP divides specific
        humidity by; it is not a WMO/Bolton saturation mixing ratio.
    """
    t = np.asarray(temperature, dtype=np.float64)
    p = np.asarray(pressure, dtype=np.float64)
    denom = t - A4
    # Guard the singularity at T=A4 (35.86 K), far below atmospheric values.
    denom = np.where(np.abs(denom) < 1.0e-6, np.sign(denom) * 1.0e-6 + 1.0e-6, denom)
    return (PQ0 / p) * np.exp(A2 * (t - A3) / denom)


def relative_humidity(
    pressure: np.ndarray,
    temperature: np.ndarray,
    specific_humidity: np.ndarray,
    *,
    percent: bool = False,
    adjust_q: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Relative humidity as in ``CALRH.f``.

    ``RH = q / q_s``, then clipped to ``[RHmin, 1]``. When the clip
    fires, Fortran also rewrites ``q``. Set ``adjust_q=True`` to return
    ``(rh, q_adjusted)``.

    Parameters
    ----------
    pressure, temperature, specific_humidity
        Pa, K, kg/kg.
    percent
        If True, return percent (0–100) as written to GRIB.
    adjust_q
        Also return specific humidity after the RH bounds check.
    """
    q = np.asarray(specific_humidity, dtype=np.float64)
    qs = saturation_specific_humidity(temperature, pressure)
    qs = np.maximum(qs, H1M12)
    rh = q / qs
    hi = rh > 1.0
    lo = rh < RHMIN
    rh = np.where(hi, 1.0, rh)
    rh = np.where(lo, RHMIN, rh)
    if adjust_q:
        q_out = np.where(hi | lo, rh * qs, q)
        rh_out = rh * 100.0 if percent else rh
        return rh_out, q_out
    return rh * 100.0 if percent else rh


def dewpoint(pressure: np.ndarray, specific_humidity: np.ndarray, temperature: np.ndarray) -> np.ndarray:
    """Dewpoint (K) from ``CALDWP.f`` + ``DEWPOINT.f``.

    Vapor pressure is ``e = p q / (eps + (1-eps) q)`` in Pa, converted
    to centibars, floored at ``1e-12`` Pa before the conversion, then
    looked up in the dewpoint table. The result is limited so it does
    not exceed the ambient temperature, matching ``CALDWP.f``.
    """
    p = np.asarray(pressure, dtype=np.float64)
    q = np.asarray(specific_humidity, dtype=np.float64)
    t = np.asarray(temperature, dtype=np.float64)
    evp = p * q / (EPS_MASS + ONEPS * q)
    evp_cb = np.maximum(H1M12, evp) * D001
    w1 = np.clip(_A_TAB * evp_cb + _B_TAB, 1.0, _DNTM1)
    jnt = np.floor(w1).astype(np.int64)
    jnt = np.clip(jnt, 1, _NT - 1)
    # Fortran is 1-based: TDP(JNT) is index jnt-1 in 0-based storage.
    td = _TDP[jnt - 1] + (w1 - jnt) * (_TDP[jnt] - _TDP[jnt - 1])
    return np.minimum(td, t)


def potential_temperature(temperature: np.ndarray, pressure: np.ndarray) -> np.ndarray:
    """Poisson potential temperature with UPP ``CAPA`` (not ``Rd/cp``)."""
    return np.asarray(temperature, dtype=np.float64) * (P1000 / np.asarray(pressure, dtype=np.float64)) ** CAPA


def temperature_from_theta(theta: np.ndarray, pressure: np.ndarray) -> np.ndarray:
    """Invert :func:`potential_temperature`."""
    return np.asarray(theta, dtype=np.float64) * (np.asarray(pressure, dtype=np.float64) / P1000) ** CAPA


def mixing_ratio_to_specific(mixing_ratio: np.ndarray) -> np.ndarray:
    """``q = r / (1+r)``, as ``INITPOST.F`` does for ``QVAPOR`` and ``Q2``."""
    r = np.asarray(mixing_ratio, dtype=np.float64)
    return r / (1.0 + r)


def fpvsnew(temperature: np.ndarray) -> np.ndarray:
    """Saturation vapor pressure (Pa), continuous form of ``fpvsnew``.

    ``CALRH_GFS.f`` tabulates the same Clausius–Clapeyron expression
    (water above 273.16 K, ice below 253.16 K, linear blend between).
    CAPE buoyancy (``CALCAPE.f``) calls ``fpvsnew`` rather than the
    ``PQ0`` formula. The table interpolation is about 6 decimal places;
    this evaluates the generating expression directly.
    """
    t = np.asarray(temperature, dtype=np.float64)
    ttp = 273.16
    psat = 610.78
    cvap = 1846.0
    cliq = 4185.5
    hvap = 2.5e6
    rv = 461.5
    csol = 2106.0
    hfus = 3.3358e5
    tliq = ttp
    tice = ttp - 20.0
    dldtl = cvap - cliq
    xponal = -dldtl / rv
    xponbl = -dldtl / rv + hvap / (rv * ttp)
    dldti = cvap - csol
    heati = hvap + hfus
    xponai = -dldti / rv
    xponbi = -dldti / rv + heati / (rv * ttp)
    # Clamp to the fpvs table range so extreme Newton steps stay defined.
    tt = np.clip(t, 180.0, 330.0)
    tr = ttp / tt
    pvl = psat * tr**xponal * np.exp(xponbl * (1.0 - tr))
    pvi = psat * tr**xponai * np.exp(xponbi * (1.0 - tr))
    w = (tt - tice) / (tliq - tice)
    return np.where(tt >= tliq, pvl, np.where(tt < tice, pvi, w * pvl + (1.0 - w) * pvi))


def specific_humidity_from_vapor_pressure(pressure: np.ndarray, vapor_pressure: np.ndarray) -> np.ndarray:
    """``q = eps e / (p - (1-eps) e)``, the form in ``CALCAPE.f``."""
    p = np.asarray(pressure, dtype=np.float64)
    e = np.minimum(np.asarray(vapor_pressure, dtype=np.float64), p * 0.99)
    e = np.maximum(e, 0.0)
    return EPS_MASS * e / (p - ONEPS * e)


def bolton_theta_e(pressure: np.ndarray, temperature: np.ndarray, specific_humidity: np.ndarray) -> np.ndarray:
    """Equivalent potential temperature from ``CALTHTE.f`` (Bolton 1980).

    Used only to choose the warmest boundary-layer parcel before
    ``CALCAPE`` (``MISCLN.f`` best-CAPE loop). The lifted parcel itself
    uses :func:`upp_py.cape.parcel_theta_e`.
    """
    p = np.asarray(pressure, dtype=np.float64)
    t = np.asarray(temperature, dtype=np.float64)
    q = np.maximum(np.asarray(specific_humidity, dtype=np.float64), 0.0)
    evp = p * q / (EPS_MASS + ONEPS * q)
    rmx = EPS_MASS * evp / np.maximum(p - evp, H1M12)
    ckapa = 0.2845 * (1.0 - 0.28 * rmx)
    rkapa = 1.0 / ckapa
    arg = np.maximum(H1M12, evp * 0.01)
    denom = 3.5 * np.log(t) - np.log(arg) - 4.805
    tlcl = 2840.0 / denom + 55.0
    fac = (P1000 / p) ** ckapa
    eterm = (3.376 / tlcl - 0.00254) * (rmx * 1.0e3 * (1.0 + 0.81 * rmx))
    return t * fac * np.exp(eterm)
