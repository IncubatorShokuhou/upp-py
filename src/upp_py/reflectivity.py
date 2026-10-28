"""Radar reflectivity for Thompson microphysics (``mp_physics=8``).

``MDLFLD.f`` (the branch dated 13 January 2009, after ``CALL PARAMR``)
follows Stoelinga: exponential size distributions for rain, snow and
graupel, with the sixth-moment reflectivity

.. math::

    Z_r = 720\\, N_{0r}\\, \\lambda_r^{-7},

and the ice species multiplied by ``alpha = 0.224`` and
``(\\rho_i/\\rho_w)^2``. Intercepts are the variable-``N_0`` fits in
``PARAMR.f`` (tanh transition for rain, Houze temperature dependence
for snow, the Roy Rasmussen graupel intercept). Bright-band boosting
is compiled off (``bb = 0``).

Composite reflectivity is the column-maximum linear ``Z`` plus a
convective ``Z = 300 R^{1.4}`` (``R`` in mm/h from the convective
step precipitation), then ``10 log10``. With no convective rain the
convective term is floored at 0.1 mm^6/m^3, so a hydrometeor-free
column is about −10 dBZ before reporting. The 1 km and 4 km AGL
products (REFD) write that clear-air / sub-zero result as 0 dBZ.
Those levels use the first mass level whose *interface* height above
the surface exceeds that depth, which is the
``DZ = ZINT(ll) - ZINT(lm+1)`` test in ``MDLFLD.f``.

Other microphysics schemes fall through to the Ferrier-style
``(q \\rho)^{1.75}`` branch in the same file (rain, snow, graupel
coefficients ``3.630803e-9``, ``2.185e-10``, ``1.033267e-9``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from upp_py.constants import (
    CONST1R,
    CONST2R,
    CONST_NG1,
    CONST_NG2,
    D608,
    DBZMIN,
    DELQR0,
    DGRAUPEL,
    DRAIN,
    DSNOW,
    GON,
    QR0,
    RD,
    REFL_ALPHA,
    REFL_R1,
    RON,
    SON,
    TOPG,
    TOPR,
    TOPS,
    ZR_A,
    ZR_B,
)
from upp_py.wrf import WrfState


@dataclass
class Reflectivity:
    """dBZ fields. Values below the floor are reported as ``DBZmin``."""

    dbz: np.ndarray  # (nz, ny, nx) 3-D equivalent reflectivity
    composite: np.ndarray  # value UPP writes as REFC
    agl_1km: np.ndarray
    agl_4km: np.ndarray
    diagnosed: np.ndarray  # Stoelinga column maximum, including convective Z


def _ze_thompson(state: WrfState) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Linear Z (mm^6/m^3) for rain, snow and graupel."""
    t = state.t
    rhod = state.p_mid / (RD * t * (state.q * D608 + 1.0))
    ze_r = np.full(t.shape, 1.0e-35)
    ze_s = np.full(t.shape, 1.0e-35)
    ze_g = np.full(t.shape, 1.0e-35)

    rain = np.maximum(REFL_R1, state.qr)
    use_r = state.qr >= 1.0e-6
    ronv = (CONST1R * np.tanh((QR0 - rain) / DELQR0) + CONST2R) / RON
    slor = (rhod * rain / (TOPR * ronv)) ** 0.25
    ze_r = np.where(use_r, (720.0 * ronv * RON * slor**7) * 1.0e18, ze_r)

    snow = np.maximum(REFL_R1, state.qs)
    use_s = state.qs >= 1.0e-6
    temp_c = np.minimum(-0.001, t - 273.15)
    sonv = np.minimum(2.0e8, 2.0e6 * np.exp(-0.12 * temp_c)) / SON
    slos = (rhod * snow / (TOPS * sonv)) ** 0.25
    ze_s = np.where(
        use_s,
        (720.0 * REFL_ALPHA * sonv * SON * slos**7 * (DSNOW / DRAIN) ** 2) * 1.0e18,
        ze_s,
    )

    graup = np.maximum(REFL_R1, state.qg)
    use_g = state.qg >= 1.0e-6
    rhoqg = rhod * graup
    gonv = CONST_NG1 * (rhoqg**CONST_NG2)
    gonv = np.clip(gonv, 1.0e4, GON) / GON
    slog = (rhoqg / (TOPG * gonv)) ** 0.25
    ze_g = np.where(
        use_g,
        (720.0 * REFL_ALPHA * gonv * GON * slog**7 * (DGRAUPEL / DRAIN) ** 2) * 1.0e18,
        ze_g,
    )
    return ze_r, ze_s, ze_g


def _ze_ferrier_like(state: WrfState) -> np.ndarray:
    """Non-Thompson ``(q rho)^1.75`` reflectivity, linear Z in mm^6/m^3."""
    rhod = state.p_mid / (RD * state.t * (state.q * D608 + 1.0))
    ze_r = (np.maximum(state.qr, 0.0) * rhod) ** 1.75 * 3.630803e-9 * 1.0e18
    ze_s = (np.maximum(state.qs, 0.0) * rhod) ** 1.75 * 2.18500e-10 * 1.0e18
    ze_g = (np.maximum(state.qg, 0.0) * rhod) ** 1.75 * 1.033267e-9 * 1.0e18
    # Schemes without ice (mp 1 or 3) put frozen condensate in snow only
    # when T < 0 is already partitioned by the reader.
    return ze_r + ze_s + ze_g


def _to_dbz(ze: np.ndarray) -> np.ndarray:
    out = np.full(ze.shape, DBZMIN, dtype=np.float64)
    good = ze > 0.0
    out = np.where(good, 10.0 * np.log10(ze), out)
    return np.maximum(out, DBZMIN)


def _reported_dbz(ze: np.ndarray) -> np.ndarray:
    """Stoelinga dBZ with clear air stored as 0.

    ``MDLFLD.f`` floors the convective contribution at 0.1 mm^6/m^3,
    and ``10 log10(0.1)`` is −10 dBZ. REFD treats that floor, and any
    other result below 0 dBZ, as clear sky.
    """
    dbz = _to_dbz(ze)
    return np.where(dbz < 0.0, 0.0, dbz)


def _level_at_agl(state: WrfState, ze: np.ndarray, depth: float) -> np.ndarray:
    """First mass level (from the ground) whose top interface exceeds ``depth``."""
    ny, nx = state.terrain.shape
    picked = np.full((ny, nx), np.nan)
    found = np.zeros((ny, nx), dtype=bool)
    z_sfc = state.z_int[0]
    for k in range(state.nz):
        agl = state.z_int[k + 1] - z_sfc
        hit = ~found & (agl > depth)
        if not np.any(hit):
            continue
        picked = np.where(hit, ze[k], picked)
        found = found | hit
        if np.all(found):
            break
    # If the column is shallower than `depth`, use the top level.
    picked = np.where(found, picked, ze[-1])
    return picked


def reflectivity(state: WrfState) -> Reflectivity:
    """Composite, 1 km and 4 km reflectivity, plus the 3-D field in dBZ.

    For Thompson on WRF-ARW (``MODELNAME=='NCAR'``), ``MDLFLD.f`` writes
    composite reflectivity as the column maximum of the model's
    ``REFL_10CM``, bounded to ``[DBZmin, 80]``. The Stoelinga diagnostic
    still runs (the unipost log prints ``dbze_max``) and is what 1 km
    and 4 km AGL fields use when the run is not RAP/HRRR. Both are
    returned: ``composite`` follows the GRIB write path, ``diagnosed``
    is the Stoelinga column maximum.
    """
    if state.mp_physics == 8:
        ze_r, ze_s, ze_g = _ze_thompson(state)
        ze = ze_r + ze_s + ze_g
    else:
        ze = _ze_ferrier_like(state)

    # Convective contribution. CPRATE is m per step; mm/h = m * 1000 / dt * 3600.
    dt = state.dt if state.dt > 0 else 1.0
    cuprate = (state.raincv / dt) * 3600.0  # mm/h ; RAINCV is already mm
    # INITPOST stores CPRATE = RAINCV_mm * 0.001 (meters), then
    # CUPRATE = 3.6e6/DT * CPRATE = RAINCV_mm / DT * 3600. Same mm/h.
    ze_conv = np.maximum(0.1, ZR_A * np.maximum(cuprate, 0.0) ** ZR_B)

    ze_max = np.max(ze, axis=0)
    # MDLFLD initializes ze_max to -10 and then takes max with ze_sum.
    # Any physical Z (positive) replaces -10, so the max of ze is correct
    # as long as ze has a positive floor. Our floor is 1e-35 * 1e18.
    diagnosed = _to_dbz(ze_max + ze_conv)
    model = getattr(state, "refl_10cm", None)
    if state.mp_physics == 8 and model is not None:
        # MDLFLD.f NCAR/NMM/FV3 branch: column max of REFL_10CM.
        composite = np.clip(np.max(model, axis=0), DBZMIN, 80.0)
    else:
        composite = diagnosed
    z1 = _level_at_agl(state, ze, 1000.0)
    z4 = _level_at_agl(state, ze, 4000.0)
    return Reflectivity(
        dbz=_to_dbz(ze),
        composite=composite,
        agl_1km=_reported_dbz(z1 + ze_conv),
        agl_4km=_reported_dbz(z4 + ze_conv),
        diagnosed=diagnosed,
    )
