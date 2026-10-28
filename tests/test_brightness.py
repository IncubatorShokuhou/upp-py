"""BRTMP stays behind the optional upp_crtm import."""

import numpy as np

from tests.test_interp_mslp import _tiny_state
from upp_py.brightness import load_brightness_temperature, skip_reason
from upp_py.config import PostConfig
from upp_py.post import dataset_from_state


def test_package_is_absent_in_the_default_environment():
    assert load_brightness_temperature() is None
    text = skip_reason()
    assert "upp_crtm.brightness_temperature" in text
    assert "TSK" not in text


def test_request_without_package_is_a_skip():
    ds = dataset_from_state(_tiny_state(), PostConfig(fields=("BRTMP",)))
    assert "brtmp" not in ds
    assert "upp_crtm" in ds.attrs["skipped"]


def test_installed_callable_is_written(monkeypatch):
    import upp_py.brightness as brightness

    def fake(state, *, satellite, channels):
        assert satellite == "abi_g17"
        assert channels == (14,)
        return np.full((1, state.ny, state.nx), 221.5)

    monkeypatch.setattr(brightness, "load_brightness_temperature", lambda: fake)
    ds = dataset_from_state(
        _tiny_state(),
        PostConfig(fields=("BRTMP",), brtmp_satellite="abi_g17", brtmp_channels=(14,)),
    )
    assert "skipped" not in ds.attrs
    np.testing.assert_allclose(ds["brtmp"].values, 221.5)
    assert list(ds["abi_channel"].values) == [14]
    assert ds["brtmp"].attrs["satellite"] == "abi_g17"


def test_bad_shape_is_skipped(monkeypatch):
    import upp_py.brightness as brightness

    monkeypatch.setattr(
        brightness,
        "load_brightness_temperature",
        lambda: lambda state, *, satellite, channels: np.zeros((3, 3)),
    )
    ds = dataset_from_state(_tiny_state(), PostConfig(fields=("BRTMP",)))
    assert "brtmp" not in ds
    assert "shape" in ds.attrs["skipped"]
