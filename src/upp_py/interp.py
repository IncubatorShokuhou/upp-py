"""Vertical interpolation from model levels to pressure, following ``MDL2P.f``.

UPP indexes ``L = 1`` at the model top. This module keeps WRF order
(``k = 0`` at the ground) and applies the same formulas:

* Above the interface between the top two mass levels, copy the top
  mass-level value (no extrapolation of T/q/wind).
* Between and just below the lowest mass level but still above the
  ground, interpolate linearly in ``ln p``:

  .. math::

      f = f_L + (f_L - f_{L-1})\\,
      \\frac{\\ln p_t - \\ln p_L}{\\ln p_L - \\ln p_{L-1}}

  with the factor clipped to ``[-10, 10]`` (the NCAR/NMM branch).

* Below the ground (``p_t > p_{sfc}``), temperature follows a constant
  virtual-temperature lapse from the average of the 2nd and 3rd model
  layers above the ground (``RGAMOG`` in ``params.F``). Relative
  humidity is held at that layer's value. Wind and omega are taken from
  the lowest mass level. Height is integrated hydrostatically downward
  from the next pressure level above, so levels must be processed from
  the top of the atmosphere downward.

Height above ground is interpolated in ``ln p`` between interface
heights ``(PH+PHB)/g``, which is the non-GFS branch of ``MDL2P.f``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from upp_py.constants import D608, G, H1M12, RD, RGAMOG, RHMIN
from upp_py.thermo import saturation_specific_humidity
from upp_py.wrf import WrfState


@dataclass
class PressureLevel:
    """One isobaric surface."""

    pressure: float  # Pa
    height: np.ndarray  # m MSL
    temperature: np.ndarray  # K
    specific_humidity: np.ndarray  # kg/kg
    u: np.ndarray
    v: np.ndarray
    omega: np.ndarray  # Pa/s
    geopotential: np.ndarray  # m^2/s^2, the FPRS array MEMSLP consumes
    underground: np.ndarray  # bool


def _gather(field: np.ndarray, k: np.ndarray) -> np.ndarray:
    """Pick ``field[k[j,i], j, i]`` with ``field`` shaped (nz, ny, nx)."""
    nz, ny, nx = field.shape
    kk = np.clip(k, 0, nz - 1)
    jj, ii = np.indices((ny, nx))
    return field[kk, jj, ii]


def _logp_pair(p_mid: np.ndarray, field: np.ndarray, p_target: float) -> tuple[np.ndarray, np.ndarray]:
    """Log-p interpolation of ``field`` onto ``p_target``.

    Returns the interpolated field and a boolean mask that is True where
    ``p_target`` is above the interface between the top two mass levels
    (UPP copies the top level there).
    """
    nz = p_mid.shape[0]
    p_iface_top = 0.5 * (p_mid[-1] + p_mid[-2])
    above = p_target < p_iface_top
    gt = p_mid > p_target
    k_idx = np.arange(nz, dtype=np.int64)[:, None, None]
    k_below = np.max(np.where(gt, k_idx, 0), axis=0)
    k_above = np.clip(k_below + 1, 0, nz - 1)
    p_lo = _gather(p_mid, k_below)
    p_hi = _gather(p_mid, k_above)
    f_lo = _gather(field, k_below)
    f_hi = _gather(field, k_above)
    denom = np.log(p_lo) - np.log(np.maximum(p_hi, 1.0))
    denom = np.where(np.abs(denom) < 1.0e-6, 1.0e-6, denom)
    fact = (np.log(p_target) - np.log(p_lo)) / denom
    fact = np.clip(fact, -10.0, 10.0)
    val = f_lo + (f_lo - f_hi) * fact
    val = np.where(above, field[-1], val)
    return val, above


def _interface_height(state: WrfState, p_target: float) -> np.ndarray:
    """ln-p interpolation of interface height (m MSL)."""
    p_int = state.p_int
    z_int = state.z_int
    nzi = p_int.shape[0]
    # UPP searches interfaces from the first one below the model lid
    # (index nz-1 in WRF order) downward. The lid itself is the upper
    # bound of that bracket.
    n_search = nzi - 1
    gt = p_int[:n_search] > p_target
    k_idx = np.arange(n_search, dtype=np.int64)[:, None, None]
    k_lo = np.max(np.where(gt, k_idx, 0), axis=0)
    k_hi = k_lo + 1
    p_lo = _gather(p_int, k_lo)
    p_hi = _gather(p_int, k_hi)
    z_lo = _gather(z_int, k_lo)
    z_hi = _gather(z_int, k_hi)
    denom = np.log(p_lo) - np.log(np.maximum(p_hi, 1.0))
    denom = np.where(np.abs(denom) < 1.0e-6, 1.0e-6, denom)
    fact = (np.log(p_target) - np.log(np.maximum(p_lo, 1.0))) / denom
    return z_lo + (z_lo - z_hi) * fact


def _underground_tq(state: WrfState, p_target: float) -> tuple[np.ndarray, np.ndarray]:
    """Below-ground T and q from the 2nd/3rd layers, ``MDL2P.f``."""
    # WRF k=1 is the 2nd layer from the ground, k=2 the 3rd.
    t_l = 0.5 * (state.t[2] + state.t[1])
    q_l = 0.5 * (state.q[2] + state.q[1])
    p_l = state.p_int[2]
    qs = saturation_specific_humidity(t_l, p_l)
    rh = q_l / np.maximum(qs, H1M12)
    rh = np.clip(rh, RHMIN, 1.0)
    q_l = rh * qs
    tv = t_l * (1.0 + D608 * q_l)
    tv_blo = tv * (p_target / p_l) ** RGAMOG
    t_blo = tv_blo / (1.0 + D608 * q_l)
    qs_t = saturation_specific_humidity(t_blo, p_target)
    q_blo = np.maximum(1.0e-12, rh * qs_t)
    return t_blo, q_blo


def interpolate_pressure_levels(state: WrfState, pressures_pa: list[float] | tuple[float, ...]) -> list[PressureLevel]:
    """Interpolate model fields onto ``pressures_pa``.

    Parameters
    ----------
    state
        Mass-point WRF state from :func:`upp_py.wrf.open_wrfout`.
    pressures_pa
        Target pressures in Pa. They are processed from smallest to
        largest regardless of input order, because underground height
        depends on the level above. The returned list follows that
        top-to-bottom order.

    Returns
    -------
    list of PressureLevel
    """
    order = np.argsort(np.asarray(pressures_pa, dtype=np.float64))
    levels_pa = [float(pressures_pa[i]) for i in order]
    out: list[PressureLevel] = []
    prev: PressureLevel | None = None

    for p_t in levels_pa:
        t_int, above = _logp_pair(state.p_mid, state.t, p_t)
        q_int, _ = _logp_pair(state.p_mid, state.q, p_t)
        u_int, _ = _logp_pair(state.p_mid, state.u, p_t)
        v_int, _ = _logp_pair(state.p_mid, state.v, p_t)
        w_int, _ = _logp_pair(state.p_mid, state.omega, p_t)
        under = p_t > state.p_sfc

        # RH clamp on the log-p branch only (MDL2P.f, not the above-top copy).
        in_column = ~above & ~under
        qs = saturation_specific_humidity(t_int, p_t)
        rh = q_int / np.maximum(qs, H1M12)
        rh = np.clip(rh, RHMIN, 1.0)
        q_adj = np.where(in_column, rh * qs, q_int)

        if np.any(under):
            t_u, q_u = _underground_tq(state, p_t)
            t_int = np.where(under, t_u, t_int)
            q_adj = np.where(under, q_u, q_adj)
            u_int = np.where(under, state.u[0], u_int)
            v_int = np.where(under, state.v[0], v_int)
            w_int = np.where(under, state.omega[0], w_int)

        z = _interface_height(state, p_t)
        if prev is not None and np.any(under):
            t_prev = prev.temperature
            q_prev = prev.specific_humidity
            tv_layer = t_prev * (1.0 + D608 * q_prev) + t_int * (1.0 + D608 * q_adj)
            gz = prev.geopotential - RD * tv_layer * np.log(p_t / prev.pressure) / 2.0
            z = np.where(under, gz / G, z)
        gz_out = z * G
        out.append(
            PressureLevel(
                pressure=p_t,
                height=z,
                temperature=t_int,
                specific_humidity=q_adj,
                u=u_int,
                v=v_int,
                omega=w_int,
                geopotential=gz_out,
                underground=under,
            )
        )
        prev = out[-1]
    return out
