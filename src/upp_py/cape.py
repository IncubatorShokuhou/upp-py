"""CAPE, CIN, LCL and LFC following ``CALCAPE.f``, ``BNDLYR.f`` and ``CLDRAD.f``.

``CALCAPE`` lifts a parcel with a constant equivalent potential
temperature and sums

.. math::

    \\mathrm{CAPE,CIN} = \\sum g\\,\\ln\\frac{\\theta_{v,p}}{\\theta_{v,a}}\\,\\Delta z

from the LCL to the equilibrium level (the highest positively buoyant
level). Positive contributions go to CAPE, negative ones to CIN, which
is stored as a negative number (``CINS = min(CINS, 0)``).

The parcel theta-e is the one written in ``CALCAPE.f``:

.. math::

    \\theta_e = \\theta \\exp\\left(L/c_p\\, q\\, (10^5/p_{LCL})^{\\kappa}/\\theta\\right),

with ``\\kappa = CAPA`` and ``L/c_p = ELOCP`` from ``params.F``. ``p_LCL``
is the pressure where the dry adiabat's saturation specific humidity
(the ``PQ0`` formula, the generating expression of ``TABLE.f``) equals
the parcel ``q``. UPP looks this up in a 76×134 table; the root is
solved directly here, so CAPE is algorithmically the same formula
without the table's quantization.

Parcel choices, matching the GRIB products on the WRF control file:

* surface (``CLDRAD.f``): most-unstable in the lowest 10 hPa
  (``ITYPE=1``, ``DPBND=10e2``; the comment in that file says 70 mb).
* best / 180–0 hPa (``MISCLN.f``): the 30 hPa boundary-layer slab
  (``BNDLYR.f``, six slabs) with the largest Bolton theta-e
  (``CALTHTE.f``), lifted as ``ITYPE=2``.
* mixed layer / 90–0 hPa: pressure-weighted mean of the first three
  slabs, lifted as ``ITYPE=2``.
* most unstable / labelled 255–0 hPa: ``ITYPE=1``, ``DPBND=300 hPa``
  (the GRIB level is 255 hPa; the search depth in ``MISCLN.f`` is 300).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from upp_py.constants import (
    CAPA,
    D608,
    DPBND_LAYER,
    DPBND_MOST_UNSTABLE,
    DPBND_SURFACE,
    ELOCP,
    EPS_MASS,
    G,
    H1M12,
    NBND,
    ONEPS,
    P1000,
)
from upp_py.thermo import bolton_theta_e, fpvsnew, saturation_specific_humidity
from upp_py.wrf import WrfState


@dataclass
class CapeResult:
    """Column convective indices for one parcel definition."""

    cape: np.ndarray  # J/kg, >= 0
    cin: np.ndarray  # J/kg, <= 0
    lcl_pressure: np.ndarray  # Pa
    lcl_height_agl: np.ndarray  # m
    lfc_pressure: np.ndarray  # Pa, NaN if no LFC
    lfc_height_agl: np.ndarray  # m
    parcel_pressure: np.ndarray  # Pa
    parcel_temperature: np.ndarray  # K
    parcel_q: np.ndarray  # kg/kg


def _qsat_adiabat(theta: np.ndarray, pressure: np.ndarray) -> np.ndarray:
    """Saturation specific humidity along a dry adiabat of ``theta``."""
    t = theta * (pressure / P1000) ** CAPA
    return saturation_specific_humidity(t, pressure)


def lcl_pressure(theta: np.ndarray, specific_humidity: np.ndarray, pressure: np.ndarray) -> np.ndarray:
    """Pressure (Pa) where the dry adiabat saturates at the parcel ``q``."""
    q = np.maximum(np.asarray(specific_humidity, dtype=np.float64), H1M12)
    p_sfc = np.asarray(pressure, dtype=np.float64)
    theta = np.asarray(theta, dtype=np.float64)
    lo = np.full(p_sfc.shape, 500.0)
    hi = np.maximum(p_sfc, 600.0)
    # qs falls as pressure falls. If qs(hi) <= q the parcel is already
    # saturated (or the formula is below q) and the LCL is the parcel pressure.
    qs_hi = _qsat_adiabat(theta, hi)
    saturated = qs_hi <= q
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        qs = _qsat_adiabat(theta, mid)
        too_moist = qs > q  # still unsaturated: LCL is higher up
        hi = np.where(too_moist, mid, hi)
        lo = np.where(too_moist, lo, mid)
    p_lcl = 0.5 * (lo + hi)
    return np.where(saturated, p_sfc, p_lcl)


def parcel_theta_e(temperature: np.ndarray, pressure: np.ndarray, specific_humidity: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``CALCAPE`` theta-e and LCL pressure."""
    t = np.asarray(temperature, dtype=np.float64)
    p = np.asarray(pressure, dtype=np.float64)
    q = np.clip(np.asarray(specific_humidity, dtype=np.float64), 0.0, 0.1)
    theta = t * (P1000 / p) ** CAPA
    p_lcl = lcl_pressure(theta, q, p)
    apes = (P1000 / np.maximum(p_lcl, 1.0)) ** CAPA
    the = theta * np.exp(ELOCP * q * apes / np.maximum(theta, 1.0))
    return the, p_lcl


def _theta_e_saturated(temperature: np.ndarray, pressure: float | np.ndarray) -> np.ndarray:
    t = np.asarray(temperature, dtype=np.float64)
    p = pressure
    theta = t * (P1000 / p) ** CAPA
    qs = np.clip(saturation_specific_humidity(t, p), 0.0, 0.1)
    return theta * np.exp(ELOCP * qs * (P1000 / p) ** CAPA / np.maximum(theta, 1.0))


def temperature_along_moist_adiabat(theta_e: np.ndarray, pressure: np.ndarray, t_guess: np.ndarray) -> np.ndarray:
    """Invert theta-e(T, p, q_s(T, p)) = ``theta_e`` (the ``TTBLEX`` step).

    A short Newton loop with a one-sided difference. ``TTBLEX`` does the
    same inversion by spline on ``TABLE.f``; the root of the same
    expression is used here.
    """
    t = np.array(t_guess, dtype=np.float64, copy=True)
    target = np.asarray(theta_e, dtype=np.float64)
    for _ in range(10):
        te = _theta_e_saturated(t, pressure)
        te2 = _theta_e_saturated(t + 0.25, pressure)
        deriv = (te2 - te) / 0.25
        deriv = np.where(np.abs(deriv) < 1.0e-3, 1.0e-3, deriv)
        t = np.clip(t - (te - target) / deriv, 150.0, 400.0)
    return t


def _integrate(
    state: WrfState,
    t_parcel: np.ndarray,
    p_parcel: np.ndarray,
    q_parcel: np.ndarray,
) -> CapeResult:
    """Lift ``t_parcel, p_parcel, q_parcel`` with ``CALCAPE``'s sum."""
    ny, nx = state.p_sfc.shape
    nz = state.nz
    the, p_lcl = parcel_theta_e(t_parcel, p_parcel, q_parcel)

    # LCL model level: highest WRF k whose pressure is still >= p_LCL
    # (UPP: first level at or below the saturation pressure).
    p = state.p_mid
    ge = p >= p_lcl[None, :, :]
    k_idx = np.arange(nz)[:, None, None]
    k_lcl = np.max(np.where(ge, k_idx, 0), axis=0)

    # Parcel temperature on model levels at and above the LCL.
    tpar = np.zeros((nz, ny, nx), dtype=np.float64)
    # Guess: dry adiabat below, previous level above. One pass upward.
    t_guess = t_parcel * (state.p_mid / p_parcel[None, :, :]) ** CAPA
    for k in range(nz - 1, -1, -1):
        # Only levels with k >= k_lcl (at or above LCL) need a moist T,
        # but computing everywhere is simpler and the sum masks them.
        if k < nz - 1:
            guess = np.where(k >= k_lcl, tpar[k + 1], t_guess[k])
        else:
            guess = t_guess[k]
        tpar[k] = temperature_along_moist_adiabat(the, state.p_mid[k], guess)

    # Equilibrium level: highest level at or above the LCL where the
    # parcel is warmer than the environment. CALCAPE only tests levels
    # with L <= LCL and leaves IEQL at the ground if none qualify, which
    # makes the subsequent sum empty.
    above_lcl = k_idx >= k_lcl[None, :, :]
    warmer = (tpar > state.t) & above_lcl
    any_warm = np.any(warmer, axis=0)
    k_eql = np.max(np.where(warmer, k_idx, 0), axis=0)

    cape = np.zeros((ny, nx), dtype=np.float64)
    cin = np.zeros((ny, nx), dtype=np.float64)
    lfc_k = np.full((ny, nx), -1, dtype=np.int64)
    # Integrate from the LCL upward to the EL. UPP loops L from IEQL to
    # LCL, i.e. from the top of the buoyant column down; the sum itself
    # does not depend on order.
    for k in range(nz):
        in_layer = any_warm & (k >= k_lcl) & (k <= k_eql)
        if not np.any(in_layer):
            continue
        pres = state.p_mid[k]
        dz = state.z_int[k + 1] - state.z_int[k]
        gdz = G * dz
        esat = np.minimum(fpvsnew(tpar[k]), pres)
        qsat = EPS_MASS * esat / np.maximum(pres - ONEPS * esat, 1.0)
        tvp = tpar[k] * (1.0 + D608 * qsat)
        tva = state.t[k] * (1.0 + D608 * state.q[k])
        thetap = tvp * (P1000 / pres) ** CAPA
        thetaa = tva * (P1000 / pres) ** CAPA
        buoy = np.log(np.maximum(thetap, 1.0)) - np.log(np.maximum(thetaa, 1.0))
        pos = in_layer & (buoy > 0.0)
        neg = in_layer & (buoy < 0.0)
        cape = cape + np.where(pos, buoy * gdz, 0.0)
        cin = cin + np.where(neg, buoy * gdz, 0.0)
        first_pos = pos & (lfc_k < 0)
        lfc_k = np.where(first_pos, k, lfc_k)

    cape = np.maximum(cape, 0.0)
    cin = np.minimum(cin, 0.0)

    # LCL height: log-p interpolate interface height to p_LCL, then AGL.
    # Bracket p_LCL between the interface below and the one above.
    p_int = state.p_int
    z_int = state.z_int
    nzi = p_int.shape[0]
    gt = p_int >= p_lcl[None, :, :]
    ki = np.arange(nzi)[:, None, None]
    k_lo = np.max(np.where(gt, ki, 0), axis=0)
    k_hi = np.clip(k_lo + 1, 0, nzi - 1)
    jj, ii = np.indices((ny, nx))
    p_lo = p_int[k_lo, jj, ii]
    p_hi = p_int[k_hi, jj, ii]
    z_lo = z_int[k_lo, jj, ii]
    z_hi = z_int[k_hi, jj, ii]
    denom = np.log(np.maximum(p_lo, 1.0)) - np.log(np.maximum(p_hi, 1.0))
    denom = np.where(np.abs(denom) < 1.0e-8, 1.0e-8, denom)
    # If the bracket collapsed (LCL at the surface), stay on z_lo.
    same = k_hi == k_lo
    fact = (np.log(np.maximum(p_lcl, 1.0)) - np.log(np.maximum(p_lo, 1.0))) / denom
    z_lcl = np.where(same, z_lo, z_lo + (z_lo - z_hi) * fact)
    lcl_agl = z_lcl - state.z_int[0]

    lfc_p = np.full((ny, nx), np.nan)
    lfc_agl = np.full((ny, nx), np.nan)
    has = lfc_k >= 0
    if np.any(has):
        lfc_p = np.where(has, state.p_mid[lfc_k, jj, ii], np.nan)
        z_lfc = state.z_mid[lfc_k, jj, ii]
        lfc_agl = np.where(has, z_lfc - state.z_int[0], np.nan)

    return CapeResult(
        cape=cape,
        cin=cin,
        lcl_pressure=p_lcl,
        lcl_height_agl=lcl_agl,
        lfc_pressure=lfc_p,
        lfc_height_agl=lfc_agl,
        parcel_pressure=p_parcel,
        parcel_temperature=t_parcel,
        parcel_q=q_parcel,
    )


def _most_unstable(state: WrfState, depth_pa: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``ITYPE=1``: max theta-e in ``psfc-depth`` .. ``psfc``."""
    p = state.p_mid
    ps = state.p_mid[0]  # UPP uses PMID of the lowest mass level as PSFCK, not PINT
    # CALCAPE: PKL >= PSFCK-DPBND and PKL <= PSFCK. PSFCK = PMID(LM).
    in_layer = (p <= ps[None, :, :]) & (p >= (ps - depth_pa)[None, :, :])
    the = np.full(p.shape, -1.0e30)
    for k in range(state.nz):
        use = in_layer[k]
        if not np.any(use):
            continue
        the_k, _ = parcel_theta_e(state.t[k], state.p_mid[k], np.maximum(state.q[k], 0.0))
        the[k] = np.where(use, the_k, -1.0e30)
    k_best = np.argmax(the, axis=0)
    jj, ii = np.indices(ps.shape)
    return state.t[k_best, jj, ii], state.p_mid[k_best, jj, ii], np.maximum(state.q[k_best, jj, ii], 0.0)


def boundary_layer_means(state: WrfState, nlayers: int = NBND, dp: float = DPBND_LAYER):
    """Pressure-weighted means in ``nlayers`` slabs of thickness ``dp``.

    Mirrors ``BNDLYR.f``: interfaces start at the surface pressure and
    step up by ``DPBND`` (30 hPa). A model level contributes to a slab
    when its mid-level pressure lies between the slab interfaces.
    Weights are interface pressure thicknesses.
    """
    ny, nx = state.p_sfc.shape
    # PBINT(1) = PINT(LM+1) = surface; PBINT(n) = PBINT(n-1) - DPBND.
    top = np.empty((nlayers, ny, nx))
    bot = np.empty((nlayers, ny, nx))
    bot[0] = state.p_sfc
    for n in range(nlayers):
        if n:
            bot[n] = top[n - 1]
        top[n] = bot[n] - dp
    p_sum = np.zeros((nlayers, ny, nx))
    t_sum = np.zeros((nlayers, ny, nx))
    q_sum = np.zeros((nlayers, ny, nx))
    for k in range(state.nz):
        pm = state.p_mid[k]
        dp_k = state.p_int[k] - state.p_int[k + 1]
        dp_k = np.maximum(dp_k, 0.0)
        for n in range(nlayers):
            inside = (bot[n] >= pm) & (top[n] <= pm)
            w = np.where(inside, dp_k, 0.0)
            p_sum[n] += w
            t_sum[n] += state.t[k] * w
            q_sum[n] += state.q[k] * w
    tiny = p_sum <= 0.0
    t_mean = np.where(tiny, state.t[0], t_sum / np.maximum(p_sum, 1.0))
    q_mean = np.where(tiny, state.q[0], q_sum / np.maximum(p_sum, 1.0))
    p_mean = 0.5 * (bot + top)
    return p_mean, t_mean, q_mean


def cape_surface(state: WrfState) -> CapeResult:
    """Surface-based CAPE/CIN (``CLDRAD.f``, lowest 10 hPa)."""
    t, p, q = _most_unstable(state, DPBND_SURFACE)
    return _integrate(state, t, p, q)


def cape_most_unstable(state: WrfState) -> CapeResult:
    """Most-unstable CAPE in the lowest 300 hPa (``MISCLN.f`` ``ITYPE=1``)."""
    t, p, q = _most_unstable(state, DPBND_MOST_UNSTABLE)
    return _integrate(state, t, p, q)


def cape_best(state: WrfState) -> CapeResult:
    """Best CAPE: Bolton-theta-e max among six 30 hPa slabs (0–180 hPa)."""
    p_m, t_m, q_m = boundary_layer_means(state)
    the = np.stack([bolton_theta_e(p_m[n], t_m[n], q_m[n]) for n in range(p_m.shape[0])], axis=0)
    n_best = np.argmax(the, axis=0)
    jj, ii = np.indices(n_best.shape)
    return _integrate(state, t_m[n_best, jj, ii], p_m[n_best, jj, ii], np.maximum(q_m[n_best, jj, ii], 0.0))


def cape_mixed_layer(state: WrfState) -> CapeResult:
    """Mixed-layer CAPE: mean of the first three 30 hPa slabs (0–90 hPa)."""
    p_m, t_m, q_m = boundary_layer_means(state)
    t = (t_m[0] + t_m[1] + t_m[2]) / 3.0
    p = (p_m[0] + p_m[1] + p_m[2]) / 3.0
    q = (q_m[0] + q_m[1] + q_m[2]) / 3.0
    return _integrate(state, t, p, np.maximum(q, 0.0))
