"""CALVIS visibility and CALICING category."""

import numpy as np

from tests.test_interp_mslp import _tiny_state
from upp_py.aerosol import mass_concentration
from upp_py.config import PostConfig
from upp_py.constants import RD
from upp_py.icing import in_flight_icing
from upp_py.post import dataset_from_state
from upp_py.visibility import surface_visibility, visibility_m


def test_clear_air_visibility_is_the_24km_cap():
    state = _tiny_state()
    vis = surface_visibility(state)
    np.testing.assert_allclose(vis, 24135.0)


def test_cloud_water_lowers_visibility():
    state = _tiny_state()
    clear = visibility_m(state.q[0], state.qc[0], state.qr[0], state.qi[0], state.qs[0], state.t[0], state.p_mid[0])
    qc = np.zeros_like(state.qc)
    qc[0] = 1.0e-4
    cloudy = visibility_m(state.q[0], qc[0], state.qr[0], state.qi[0], state.qs[0], state.t[0], state.p_mid[0])
    assert np.all(cloudy < clear)
    assert np.all(cloudy < 24135.0)


def test_convective_rain_lowers_visibility():
    dry = _tiny_state()
    wet = _tiny_state(raincv=np.full((4, 4), 2.0))
    assert np.all(surface_visibility(wet) < surface_visibility(dry))


def test_calicing_thresholds():
    t = np.array([270.0, 270.0, 270.0, 250.0, 273.0, 274.0])
    rh = np.array([80.0, 80.0, 69.0, 90.0, 80.0, 90.0])
    w = np.array([-0.2, 0.1, -0.2, -0.2, -0.2, -0.2])
    got = in_flight_icing(t, rh, w)
    np.testing.assert_array_equal(got, [1.0, 0.0, 0.0, 0.0, 1.0, 0.0])


def test_smoke_mass_uses_dry_density():
    p = np.array([100000.0])
    t = np.array([280.0])
    q = np.array([2.0])  # ug/kg
    rho = p / (RD * t)
    np.testing.assert_allclose(mass_concentration(q, p, t), q * rho)


def test_posted_visibility_and_tipd_on_a_cold_ascent():
    state = _tiny_state()
    state.omega[:] = -0.5
    # 700 hPa is index 2 of the synthetic column, temperature is on the
    # linspace from 290 K to 230 K, so it sits inside 251-273 K.
    ds = dataset_from_state(state, PostConfig(fields=("VIS", "TIPD"), pressure_pa=(70000.0,)))
    np.testing.assert_allclose(ds["vis"].values, 24135.0)
    assert ds["tipd"].sizes["isobaric"] == 1
    # RH of q=0.005 at 700 hPa and ~260 K is high, so the category is 1.
    assert float(ds["tipd"].values.mean()) == 1.0


def test_qnwfa_skipped_unless_mp_physics_28():
    state = _tiny_state(qnwfa=np.ones((5, 4, 4)))
    ds = dataset_from_state(state, PostConfig(fields=("QNWFA",)))
    assert "qnwfa" not in ds
    assert "imp_physics==28" in ds.attrs["skipped"]
    state.mp_physics = 28
    posted = dataset_from_state(state, PostConfig(fields=("QNWFA",)))
    np.testing.assert_allclose(posted["qnwfa"].sel(model_level=0).values, 1.0)
    tiny = dataset_from_state(
        _tiny_state(mp_physics=28, qnwfa=np.full((5, 4, 4), 1.0e-12)),
        PostConfig(fields=("QNWFA",)),
    )
    np.testing.assert_allclose(tiny["qnwfa"].values, 0.0)
