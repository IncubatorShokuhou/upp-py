"""Storm-relative helicity, ``CALHEL.f`` (Bunkers et al. 1998).

Storm motion is 7.5 m/s to the right of the 0–6 km mean wind, in the
direction perpendicular to the shear between the 0–0.5 km mean wind and
the 5.5–6 km mean wind. Means are arithmetic averages of mass-point
winds (not mass-weighted), which is the change Manikin made in 2003.

Helicity over a depth ``D`` (1 km or 3 km) is the Davies-Jones sum

.. math::

    H = \\sum\\left[
      (v-v_s)\\,\\widetilde{\\partial_z u}
      - (u-u_s)\\,\\widetilde{\\partial_z v}
    \\right]\\Delta z

with the vertical derivatives estimated from the levels above and below
the way ``CALHEL.f`` writes them. On an A-grid the four-point averages
in that routine collapse to the mass-point value, which is the WRF-ARW
case (``gridtype='A'`` in ``INITPOST.F``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from upp_py.wrf import WrfState


@dataclass
class HelicityResult:
    """0–1 km and 0–3 km helicity plus Bunkers storm motion."""

    hel_0_3km: np.ndarray  # m^2/s^2
    hel_0_1km: np.ndarray
    u_storm: np.ndarray  # m/s, grid-relative
    v_storm: np.ndarray
    u_shear_0_6: np.ndarray
    v_shear_0_6: np.ndarray


def storm_relative_helicity(state: WrfState) -> HelicityResult:
    """Compute Bunkers storm motion and 0–1 / 0–3 km SRH.

    Heights are AGL relative to the surface interface (``ZINT``), not
    the ``HGT`` field, matching ``HTSFC = FIS/g`` with ``FIS`` taken
    from the same geopotential UPP stored. For ARW, ``FIS = HGT * g``
    while ``ZINT`` at the ground is ``(PH+PHB)/g``; those two differ by
    centimeters, and ``CALHEL`` uses ``FIS/g``. This function uses
    ``terrain`` (``HGT``) so the AGL depths match the Fortran.
    """
    z_agl = state.z_mid - state.terrain[None, :, :]
    u = state.u
    v = state.v
    ny, nx = state.terrain.shape

    def _mean(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        count = mask.sum(axis=0)
        uu = np.where(mask, u, 0.0).sum(axis=0)
        vv = np.where(mask, v, 0.0).sum(axis=0)
        good = count > 0
        um = np.where(good, uu / np.maximum(count, 1), 0.0)
        vm = np.where(good, vv / np.maximum(count, 1), 0.0)
        return um, vm, good

    u6, v6, ok6 = _mean(z_agl <= 6000.0)
    u1, v1, ok1 = _mean(z_agl < 500.0)
    u5, v5, ok5 = _mean((z_agl >= 5500.0) & (z_agl < 6000.0))
    # If nothing landed in 5.5–6 km, use the first level in 6–7 km
    # encountered from the ground up (CALHEL loops L from LM down, i.e.
    # from the ground). The first one sets the mean and stops.
    need = ~ok5
    if np.any(need):
        for k in range(state.nz):
            hit = need & (z_agl[k] >= 6000.0) & (z_agl[k] < 7000.0)
            if not np.any(hit):
                continue
            u5 = np.where(hit, u[k], u5)
            v5 = np.where(hit, v[k], v5)
            ok5 = ok5 | hit
            need = ~ok5
            if not np.any(need):
                break

    ushr = u5 - u1
    vshr = v5 - v1
    denom = ushr**2 + vshr**2
    ok = ok6 & ok1 & ok5 & (denom != 0.0)
    scale = np.zeros_like(denom)
    np.divide(7.5, np.sqrt(np.maximum(denom, 1e-12)), out=scale, where=ok)
    u_storm = np.where(ok, u6 + scale * vshr, 0.0)
    v_storm = np.where(ok, v6 - scale * ushr, 0.0)

    def _helicity(depth: float) -> np.ndarray:
        heli = np.zeros((ny, nx), dtype=np.float64)
        # UPP L = 2..LM-1. WRF k = nz-L runs from 1 to nz-2.
        for k in range(1, state.nz - 1):
            dzabv = z_agl[k]
            use = dzabv < depth
            if not np.any(use):
                continue
            z1 = state.z_mid[k - 1]  # below
            z2 = state.z_mid[k]
            z3 = state.z_mid[k + 1]  # above
            dz = state.z_int[k + 1] - state.z_int[k]
            dz1 = z1 - z2
            dz2 = z2 - z3
            du1 = u[k - 1] - u[k]
            du2 = u[k] - u[k + 1]
            dv1 = v[k - 1] - v[k]
            dv2 = v[k] - v[k + 1]
            denom_z = dz1 + dz2
            safe = use & (np.abs(dz1) > 1.0e-3) & (np.abs(dz2) > 1.0e-3) & (np.abs(denom_z) > 1.0e-3)
            du = dz2 * (du1 / np.where(np.abs(dz1) < 1e-3, 1.0, dz1)) + dz1 * (
                du2 / np.where(np.abs(dz2) < 1e-3, 1.0, dz2)
            )
            dv = dz2 * (dv1 / np.where(np.abs(dz1) < 1e-3, 1.0, dz1)) + dz1 * (
                dv2 / np.where(np.abs(dz2) < 1e-3, 1.0, dz2)
            )
            inc = ((v[k] - v_storm) * du - (u[k] - u_storm) * dv) * dz / np.where(
                np.abs(denom_z) < 1e-3, 1.0, denom_z
            )
            heli = heli + np.where(safe, inc, 0.0)
        return heli

    return HelicityResult(
        hel_0_3km=_helicity(3000.0),
        hel_0_1km=_helicity(1000.0),
        u_storm=u_storm,
        v_storm=v_storm,
        u_shear_0_6=np.where(ok, ushr, 0.0),
        v_shear_0_6=np.where(ok, vshr, 0.0),
    )
