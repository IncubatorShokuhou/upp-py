"""Log-p interpolation and the flat-terrain limit of MEMSLP."""

import numpy as np

from upp_py.constants import G, SPL_PA
from upp_py.interp import interpolate_pressure_levels
from upp_py.mslp import membrane_slp
from upp_py.wrf import WrfState


def _tiny_state(**over):
    nz, ny, nx = 5, 4, 4
    z_int = np.array([0.0, 200.0, 1000.0, 3000.0, 8000.0, 16000.0])[:, None, None] * np.ones((1, ny, nx))
    p_mid = np.array([95000.0, 85000.0, 70000.0, 50000.0, 30000.0])[:, None, None] * np.ones((1, ny, nx))
    p_int = np.empty((nz + 1, ny, nx))
    p_int[0] = 100000.0
    p_int[-1] = 10000.0
    p_int[1:-1] = 0.5 * (p_mid[:-1] + p_mid[1:])
    t = np.linspace(290, 230, nz)[:, None, None] * np.ones((1, ny, nx))
    q = np.full((nz, ny, nx), 0.005)
    u = np.zeros((nz, ny, nx))
    u[:, :, :] = np.linspace(1, 10, nz)[:, None, None]
    base = dict(
        path="synthetic",
        valid_time=__import__("datetime").datetime(2018, 7, 23, 21),
        start_time=__import__("datetime").datetime(2018, 7, 23, 21),
        forecast_hour=0.0,
        mp_physics=8,
        dx=3000.0,
        dy=3000.0,
        map_proj=1,
        truelat1=30.0,
        truelat2=60.0,
        stand_lon=-120.0,
        cen_lat=40.0,
        cen_lon=-120.0,
        lat=np.full((ny, nx), 40.0),
        lon=np.full((ny, nx), -120.0),
        mapfac=np.ones((ny, nx)),
        terrain=np.zeros((ny, nx)),
        landmask=None,
        p_top=10000.0,
        p_mid=p_mid,
        p_sfc=np.full((ny, nx), 100000.0),
        p_int=p_int,
        z_int=z_int,
        z_mid=0.5 * (z_int[:-1] + z_int[1:]),
        t=t,
        q=q,
        u=u,
        v=np.zeros_like(u),
        w=np.zeros_like(u),
        omega=np.zeros_like(u),
        qc=np.zeros_like(u),
        qi=np.zeros_like(u),
        qr=np.zeros_like(u),
        qs=np.zeros_like(u),
        qg=np.zeros_like(u),
        th2=np.full((ny, nx), 300.0),
        q2_mixing=np.full((ny, nx), 0.008),
        t2=np.full((ny, nx), 290.0),
        q2=np.full((ny, nx), 0.008),
        p_shelter=np.full((ny, nx), 99900.0),
        u10=np.ones((ny, nx)),
        v10=np.ones((ny, nx)),
        rainc=np.zeros((ny, nx)),
        rainnc=np.zeros((ny, nx)),
        raincv=np.zeros((ny, nx)),
        rainncv=np.zeros((ny, nx)),
        dt=60.0,
    )
    base.update(over)
    return WrfState(**base)


def test_logp_recovers_a_level():
    state = _tiny_state()
    levels = interpolate_pressure_levels(state, [85000.0])
    # 850 hPa is a mass level, so T and u must match that level.
    np.testing.assert_allclose(levels[0].temperature, state.t[1], rtol=1e-5, atol=1e-3)
    np.testing.assert_allclose(levels[0].u, state.u[1], rtol=1e-5, atol=1e-4)


def test_logp_midpoint_between_levels():
    state = _tiny_state()
    # Geometric-mean pressure between 850 and 700 hPa.
    p_t = float(np.exp(0.5 * (np.log(85000.0) + np.log(70000.0))))
    levels = interpolate_pressure_levels(state, [p_t])
    expected = 0.5 * (state.t[1] + state.t[2])
    np.testing.assert_allclose(levels[0].temperature, expected, rtol=1e-4, atol=0.05)


def test_shuell_flat_is_surface_pressure():
    from upp_py.mslp import shuell_slp

    state = _tiny_state()
    np.testing.assert_allclose(shuell_slp(state), state.p_sfc)


def test_flat_terrain_slp_is_surface_pressure():
    state = _tiny_state()
    # |FIS| = 0 < 1, so MEMSLP returns PSFC without integrating.
    levels = interpolate_pressure_levels(state, list(SPL_PA))
    pslp = membrane_slp(state, levels)
    np.testing.assert_allclose(pslp, state.p_sfc)


def test_just_above_ground_stays_on_logp():
    state = _tiny_state()
    # Lowest mass level is 950 hPa, surface is 1000 hPa. 970 hPa is below
    # the mass level and still above the ground, so MDL2P keeps log-p.
    p_t = 97000.0
    levels = interpolate_pressure_levels(state, [p_t])
    p_lo, p_hi = state.p_mid[0, 0, 0], state.p_mid[1, 0, 0]
    fact = (np.log(p_t) - np.log(p_lo)) / (np.log(p_lo) - np.log(p_hi))
    expected = state.t[0] + (state.t[0] - state.t[1]) * fact
    np.testing.assert_allclose(levels[0].temperature, expected, rtol=1e-6, atol=1e-4)
    assert not np.any(levels[0].underground)


def test_below_ground_uses_second_and_third_layers():
    from upp_py.constants import D608, RGAMOG
    from upp_py.thermo import saturation_specific_humidity

    state = _tiny_state(p_sfc=np.full((4, 4), 80000.0))
    state.p_int[0] = state.p_sfc
    p_t = 85000.0
    levels = interpolate_pressure_levels(state, [70000.0, p_t])
    under = [lv for lv in levels if lv.pressure == p_t][0]
    assert np.all(under.underground)
    t_l = 0.5 * (state.t[2] + state.t[1])
    q_l = 0.5 * (state.q[2] + state.q[1])
    p_l = state.p_int[2]
    qs = saturation_specific_humidity(t_l, p_l)
    rh = np.clip(q_l / qs, 1.0e-6, 1.0)
    q_use = np.where((q_l / qs > 1.0) | (q_l / qs < 1.0e-6), rh * qs, q_l)
    tv = t_l * (1.0 + D608 * q_use)
    t_blo = tv * (p_t / p_l) ** RGAMOG / (1.0 + D608 * q_use)
    np.testing.assert_allclose(under.temperature, t_blo, rtol=1e-6, atol=1e-4)
    # Wind is copied from the lowest mass level.
    np.testing.assert_allclose(under.u, state.u[0])


def test_underground_temperature_is_finite():
    terrain = np.full((4, 4), 2500.0)
    # Raise the ground. Surface pressure well below 1000 hPa.
    state = _tiny_state(terrain=terrain, p_sfc=np.full((4, 4), 75000.0))
    state.z_int = state.z_int + 2500.0
    state.z_mid = 0.5 * (state.z_int[:-1] + state.z_int[1:])
    levels = interpolate_pressure_levels(state, [100000.0, 85000.0])
    # 1000 hPa is below 750 hPa surface.
    assert np.all(levels[1].underground) or levels[0].pressure == 85000.0
    under = [lv for lv in levels if lv.pressure == 100000.0][0]
    assert np.all(np.isfinite(under.temperature))
    assert np.all(under.temperature > 200.0)
    pslp = membrane_slp(state, interpolate_pressure_levels(state, list(SPL_PA)))
    assert np.all(np.isfinite(pslp))
    # Reduced pressure should exceed the surface pressure on a mountain.
    assert np.all(pslp > state.p_sfc - 1.0)
    _ = G
