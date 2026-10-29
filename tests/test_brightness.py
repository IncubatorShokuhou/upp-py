"""BRTMP calls crtm_py.brightness_temperature and skips without it."""

import numpy as np

from tests.test_interp_mslp import _tiny_state
from upp_py.brightness import call_kwargs, coeff_skip_reason, load_brightness_temperature, skip_reason
from upp_py.config import PostConfig
from upp_py.post import dataset_from_state
import upp_py.brightness as brightness


def test_missing_package_skip_names_crtm_py(monkeypatch):
    monkeypatch.setattr(brightness, "load_brightness_temperature", lambda: None)
    text = skip_reason()
    assert "crtm_py.brightness_temperature" in text
    assert "IncubatorShokuhou/crtm-py" in text
    assert "TSK" not in text
    ds = dataset_from_state(_tiny_state(), PostConfig(fields=("BRTMP",)))
    assert "brtmp" not in ds
    assert "crtm_py" in ds.attrs["skipped"]


def test_call_kwargs_follow_calrad_names():
    state = _tiny_state()
    state.xland = np.full((state.ny, state.nx), 2.0)
    state.tsk = np.full((state.ny, state.nx), 301.0)
    state.mp_physics = 8
    kwargs = call_kwargs(state, satellite="abi_g16", channels=(7, 14))
    assert kwargs["sensor"] == "abi_g16"
    assert kwargs["channels"] == (7, 14)
    assert kwargs["vertical_order"] == "bottom_up"
    assert kwargs["imp_physics"] == 8
    assert kwargs["include_scattering"] is True
    np.testing.assert_allclose(kwargs["sea_mask"], 1.0)
    np.testing.assert_allclose(kwargs["pressure_interface"], state.p_int)
    np.testing.assert_allclose(kwargs["specific_humidity"], state.q)
    np.testing.assert_allclose(kwargs["_qc"], state.qc)
    np.testing.assert_allclose(kwargs["skin_temperature"], 301.0)


def test_no_coefficients_are_skipped(monkeypatch):
    class _Module:
        def coeff_path(self):
            return None

        def brightness_temperature(self, **kwargs):
            raise AssertionError("synthetic forward should not run without coefficients")

    monkeypatch.setattr(brightness, "load_brightness_temperature", lambda: _Module())
    ds = dataset_from_state(_tiny_state(), PostConfig(fields=("BRTMP",)))
    assert "brtmp" not in ds
    assert "SpcCoeff" in ds.attrs["skipped"]
    assert "abi_g16" in coeff_skip_reason("abi_g16")


def test_grid_result_is_stored_as_channel_first(monkeypatch):
    state = _tiny_state()

    class _Result:
        optical_depth_model = "odps+cloudcoeff+ada"
        brightness_temperature = np.full((state.ny, state.nx, 2), 250.0)

    class _Module:
        def coeff_path(self):
            return "/tmp/coeffs"

        def Hydrometeors(self, **kwargs):
            return kwargs

        def surface_from_upp(self, **kwargs):
            return kwargs

        def brightness_temperature(self, **kwargs):
            assert kwargs["sensor"] == "abi_g17"
            assert kwargs["channels"] == (13, 14)
            assert kwargs["surface"]["ivegsrc"] == state.ivegsrc
            return _Result()

    monkeypatch.setattr(brightness, "load_brightness_temperature", lambda: _Module())
    ds = dataset_from_state(
        state,
        PostConfig(fields=("BRTMP",), brtmp_satellite="abi_g17", brtmp_channels=(13, 14)),
    )
    assert "skipped" not in ds.attrs
    assert ds["brtmp"].shape == (1, 2, state.ny, state.nx)
    np.testing.assert_allclose(ds["brtmp"].values, 250.0)
    assert list(ds["abi_channel"].values) == [13, 14]


def test_loader_is_none_when_the_package_is_absent():
    if load_brightness_temperature() is not None:
        return
    assert load_brightness_temperature() is None
