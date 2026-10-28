"""Membrane (Mesinger) sea-level pressure, ``SLP_new.f`` / ``MEMSLP``.

UPP's WRF-ARW path calls ``MEMSLP`` (the log line is
``CALLING MEMSLP for A or B grid``). The routine:

1. Marks each standard pressure level above or below the ground.
2. Relaxes virtual temperature on underground points with an 8-point
   Jacobi smoother (``NRLX2 = 100``, coefficient ``OVERRC = 1.10``).
   Above-ground points that touch an underground neighbor are held at
   ``T(1+0.608 q)`` as the boundary condition. The stencil always uses
   the previous sweep, which is what makes the result independent of
   MPI decomposition.
3. Integrates the hydrostatic equation downward from the lowest
   above-ground pressure level. Where geopotential crosses zero,

   .. math::

      p_{sl} = p_1 \\exp\\left(\\frac{\\Phi_1}{R_d T}\\right).

4. Where the surface itself is already below the lowest output pressure
   (usually 1000 hPa) or the integration never reaches ``\\Phi = 0``,
   fall back to a fixed lapse ``-6.6 \\times 10^{-4}\\,\\mathrm{K\\,m^{-2}\\,s^{2}}``
   applied to the lowest layer, the "standard" reduction in the same
   subroutine. Flat points (``|\\Phi_s| < 1\\,\\mathrm{m^2\\,s^{-2}}``)
   keep the surface pressure.

``TPRES`` / ``QPRES`` / ``FIPRES`` are the isobaric fields produced by
:func:`upp_py.interp.interpolate_pressure_levels` on the full ``SPLDEF``
list in ``CTLBLK.f``.
"""

from __future__ import annotations

import numpy as np

from upp_py.constants import AD05, CFT0, D608, G, GAMMA, GI, NRLX_SLP, RD, SLOPE_SLP
from upp_py.interp import PressureLevel
from upp_py.wrf import WrfState


def _relax_level(t_init: np.ndarray, underground: np.ndarray, tv_bc: np.ndarray) -> np.ndarray:
    """100 Jacobi passes. ``underground`` points move; everyone else is fixed."""
    ttv = np.where(underground, t_init, tv_bc)
    # tv_bc is only meaningful on the above-ground rim; elsewhere keep t_init
    # for above-ground points that are not a boundary condition.
    fixed = ~underground
    ttv = np.where(fixed, tv_bc, t_init)
    for _ in range(NRLX_SLP):
        c = ttv
        a1 = c[1:-1, 0:-2]
        a2 = c[1:-1, 2:]
        a3 = c[0:-2, 1:-1]
        a4 = c[2:, 1:-1]
        a5 = c[0:-2, 0:-2]
        a6 = c[0:-2, 2:]
        a7 = c[2:, 0:-2]
        a8 = c[2:, 2:]
        center = c[1:-1, 1:-1]
        updated = AD05 * (4.0 * (a1 + a2 + a3 + a4) + a5 + a6 + a7 + a8) - CFT0 * center
        ttv = ttv.copy()
        inner_under = underground[1:-1, 1:-1]
        ttv[1:-1, 1:-1] = np.where(inner_under, updated, center)
    return ttv


def shuell_slp(state: WrfState) -> np.ndarray:
    """NCEP Shuell sea-level pressure (Pa), ``NGMSLP.f``.

    The grib2 WRF control file posts this as ``PRMSL`` (available-field
    index 105). It is the hydrostatic reduction

    .. math::

        p_{sl} = p_{sfc}\\exp(z_{sfc}/\\bar\\tau),\\quad
        \\tau = T_v R_d/g

    with a 6.5 K/km lapse from the lowest mass level and the Shuell cap
    at ``T = 290.66 K``. ``MEMSLP`` is a different field (``MSLET``).
    """
    z_sfc = state.z_int[0]
    z_mid = state.z_mid[0]
    psfc = state.p_sfc
    tvrt = state.t[0] * (1.0 + D608 * state.q[0])
    # NGMSLP.f after the 2006 tau-at-surface fix.
    tv_sfc = tvrt + (z_mid - z_sfc) * GAMMA
    tv_sl = tvrt + (z_mid - 0.0) * GAMMA
    tau_sfc = tv_sfc * RD * GI
    tau_sl = tv_sl * RD * GI
    taucr = RD * GI * 290.66
    const = 0.005 * G / RD
    both = (tau_sl > taucr) & (tau_sfc > taucr)
    only_sl = (tau_sl > taucr) & (tau_sfc <= taucr)
    tau_sl = np.where(only_sl, taucr, tau_sl)
    tau_sl = np.where(both, taucr - const * (tau_sfc - taucr) ** 2, tau_sl)
    tau_avg = 0.5 * (tau_sl + tau_sfc)
    fis = state.terrain * G
    pslp = psfc * np.exp(z_sfc / tau_avg)
    return np.where(np.abs(fis) > 1.0, pslp, psfc)


def membrane_slp(state: WrfState, levels: list[PressureLevel]) -> np.ndarray:
    """Mesinger membrane SLP (Pa) on the mass grid.

    Parameters
    ----------
    state
        WRF state. Surface geopotential is ``terrain * g``.
    levels
        Isobaric levels from top to bottom, including underground
        extrapolations. Must be the standard ``SPLDEF`` set (or any set
        that extends to 1000 hPa) so the downward integration has
        somewhere to land.

    Returns
    -------
    numpy.ndarray
        Sea-level pressure (Pa), shape ``(ny, nx)``.
    """
    if len(levels) < 2:
        raise ValueError("MEMSLP needs at least two pressure levels")
    ny, nx = state.p_sfc.shape
    nlev = len(levels)
    plev = np.array([lv.pressure for lv in levels], dtype=np.float64)
    # (nlev, ny, nx)
    tprs = np.stack([lv.temperature for lv in levels], axis=0)
    qprs = np.stack([lv.specific_humidity for lv in levels], axis=0)
    fprs = np.stack([lv.geopotential for lv in levels], axis=0)
    psfc = state.p_sfc
    fis = state.terrain * G
    # Above ground when the level pressure is less than the surface pressure.
    above = plev[:, None, None] < psfc[None, :, :]

    # Lowest above-ground level at each point (largest index still above ground).
    idx = np.arange(nlev)[:, None, None]
    lmho = np.max(np.where(above, idx, 0), axis=0)

    # Highest level index that still contains any underground point.
    lhmnt = nlev  # LSM+1 in 1-based, i.e. "no underground"
    for lev in range(nlev - 1, -1, -1):
        if np.any(~above[lev]):
            continue
        lhmnt = lev + 1
        break
    if lhmnt > nlev - 1:
        lhmnt = nlev  # nothing to relax

    t_relaxed = tprs.copy()
    if lhmnt < nlev:
        for lev in range(lhmnt, nlev):
            under = ~above[lev]
            # Above-ground points with any of 8 neighbors underground
            # become a virtual-temperature boundary condition.
            ht = above[lev].astype(np.float64)
            prod = (
                ht[1:-1, 0:-2]
                * ht[1:-1, 2:]
                * ht[0:-2, 1:-1]
                * ht[2:, 1:-1]
                * ht[0:-2, 0:-2]
                * ht[0:-2, 2:]
                * ht[2:, 0:-2]
                * ht[2:, 2:]
            )
            rim = np.zeros((ny, nx), dtype=bool)
            rim[1:-1, 1:-1] = above[lev, 1:-1, 1:-1] & (prod < 0.5)
            tv = tprs[lev] * (1.0 + D608 * qprs[lev])
            bc = np.where(rim, tv, tprs[lev])
            t_relaxed[lev] = _relax_level(tprs[lev], under, bc)

    pslp = psfc.copy()
    done = np.zeros((ny, nx), dtype=bool)

    flat = np.abs(fis) < 1.0
    pslp = np.where(flat, psfc, pslp)
    done = done | flat

    # Below-sea-level terrain: integrate model levels. Rare for ARW land
    # cases; kept so the branch in SLP_new.f is not dropped.
    below_sea = fis < -1.0
    if np.any(below_sea):
        # UPP walks model interfaces from the bottom until height > 0.
        # z_int[0] is the surface (negative), z_int[-1] is the lid.
        z_int = state.z_int
        p_int = state.p_int
        nzi = z_int.shape[0]
        for k in range(1, nzi):
            crossed = below_sea & ~done & (z_int[k] > 0.0) & (z_int[k - 1] <= 0.0)
            if not np.any(crossed):
                # Also the original test is ZINT(L) > 0 while scanning
                # downward; accept the first interface above sea level.
                crossed = below_sea & ~done & (z_int[k] > 0.0)
            if not np.any(crossed):
                continue
            # tem = 0.5*(T(L)+T(L-1))*(1+0.5*D608*(Q(L)+Q(L-1)))
            # L is the UPP level of this interface's lower mass level.
            # Use the mass level just below the interface when it exists.
            k_mass = min(k - 1, state.nz - 1)
            k_mass_up = min(k, state.nz - 1)
            tem = 0.5 * (state.t[k_mass] + state.t[k_mass_up]) * (
                1.0 + 0.5 * D608 * (state.q[k_mass] + state.q[k_mass_up])
            )
            p_est = p_int[k] / np.exp(-z_int[k] * G / (RD * tem))
            pslp = np.where(crossed, p_est, pslp)
            done = done | crossed

    # Mountain points underground at the lowest pressure level.
    # SLP_new.f gathers these only on the interior (I=2..IM-1, J=2..JM-1);
    # the rim uses the standard reduction below.
    under_bottom = ~above[-1]
    interior = np.zeros(under_bottom.shape, dtype=bool)
    interior[1:-1, 1:-1] = True
    kmm = under_bottom & ~done & interior
    if np.any(kmm):
        # Integrate downward from LMHO.
        gz1 = np.take_along_axis(fprs, lmho[None, :, :], axis=0)[0]
        p1 = plev[lmho]
        active = kmm.copy()
        for lev in range(nlev):
            # Only walk levels strictly below LMHO. Vectorised by masking.
            step = active & (lev > lmho) & (lev < nlev)
            if not np.any(step):
                continue
            t_here = t_relaxed[lev]
            t_prev = t_relaxed[lev - 1]
            tlyr = 0.5 * (t_here + t_prev)
            p2 = plev[lev]
            gz2 = gz1 + RD * tlyr * np.log(p1 / p2)
            hit = step & (gz2 <= 0.0)
            p_hit = p1 / np.exp(-gz1 / (RD * t_prev))
            pslp = np.where(hit, p_hit, pslp)
            done = done | hit
            active = active & ~hit
            # Advance the integrator only on points still walking.
            gz1 = np.where(step & ~hit, gz2, gz1)
            p1 = np.where(step & ~hit, p2, p1)
        # Still above sea level at the bottom of the pressure set.
        # SLP_new.f uses the geopotential just integrated with the relaxed
        # temperatures (FIPRES overwritten on the way down), not the
        # pre-relaxation underground height from MDL2P.
        remain = active & ~done
        if np.any(remain):
            tlyr = t_relaxed[-1] - 0.5 * gz1 * SLOPE_SLP
            p_est = p1 / np.exp(-gz1 / (RD * tlyr))
            pslp = np.where(remain, p_est, pslp)
            done = done | remain

    # Surface pressure already higher than the lowest isobaric level:
    # standard lapse reduction from the ground (SLP_new.f label 325).
    remain = ~done
    if np.any(remain):
        slope = SLOPE_SLP
        use_sfc = psfc > plev[-1]
        tvrt = state.t[0] * (1.0 + D608 * state.q[0])
        # DIS = ZINT(LM+1) - ZINT(LM) + 0.5*ZINT(LM+1)
        # WRF: z_int[0] is ZINT(LM+1), z_int[1] is ZINT(LM).
        dis = state.z_int[0] - state.z_int[1] + 0.5 * state.z_int[0]
        tlyr_sfc = tvrt - dis * G * slope
        p_from_sfc = psfc * np.exp(state.z_int[0] * G / (RD * tlyr_sfc))
        tlyr_bot = t_relaxed[-1] - 0.5 * fprs[-1] * slope
        p_from_bot = plev[-1] / np.exp(-fprs[-1] / (RD * tlyr_bot))
        pslp = np.where(remain & use_sfc, p_from_sfc, pslp)
        pslp = np.where(remain & ~use_sfc, p_from_bot, pslp)
    return pslp
