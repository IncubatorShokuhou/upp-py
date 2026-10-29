"""GRIB2 smoke test. Skips when the eccodes extra is not installed."""

import numpy as np
import pytest
import xarray as xr

eccodes = pytest.importorskip("eccodes")

from upp_py.grib2 import write_grib2  # noqa: E402


def _sample_dataset():
    lat = np.array([[40.0, 40.0], [40.1, 40.1]])
    lon = np.array([[255.0, 255.1], [255.0, 255.1]])
    t = np.array([[[290.25, 291.5], [292.0, 293.75]]])
    apcp = np.array([[[0.0, 1.5], [2.25, 4.0]]])
    ds = xr.Dataset(
        data_vars={
            "t2m": (
                ("time", "y", "x"),
                t,
                {"grib2_shortName": "TMP", "long_name": "2 m temperature", "units": "K", "height": "2 m"},
            ),
            "apcp": (
                ("time", "y", "x"),
                apcp,
                {
                    "grib2_shortName": "APCP",
                    "long_name": "total precipitation",
                    "units": "kg m-2",
                    "grib_stat": "accumulation",
                },
            ),
            "crs": xr.DataArray(np.int32(0), attrs={"grid_mapping_name": "latitude_longitude"}),
        },
        coords={
            "time": np.array(["2018-07-23T22:00:00"], dtype="datetime64[ns]"),
            "y": [0, 1],
            "x": [0, 1],
            "latitude": (("y", "x"), lat),
            "longitude": (("y", "x"), lon),
            "forecast_hour": ("time", np.array([1.0])),
        },
        attrs={"forecast_reference_time": "2018-07-23T21:00:00Z", "DX": 1000.0, "DY": 1000.0},
    )
    return ds


def _messages(path):
    out = []
    with open(path, "rb") as handle:
        while True:
            gid = eccodes.codes_grib_new_from_file(handle)
            if gid is None:
                break
            out.append(
                {
                    "short": eccodes.codes_get(gid, "shortName"),
                    "values": np.array(eccodes.codes_get_values(gid), dtype=np.float64),
                    "stepType": eccodes.codes_get(gid, "stepType"),
                    "stepRange": eccodes.codes_get(gid, "stepRange"),
                    "level": eccodes.codes_get(gid, "level"),
                    "centre": eccodes.codes_get(gid, "centre"),
                }
            )
            eccodes.codes_release(gid)
    return out


def test_write_grib2_roundtrip(tmp_path):
    ds = _sample_dataset()
    path = tmp_path / "post.grib2"
    skipped = write_grib2(ds, str(path))
    assert skipped == []
    msgs = _messages(path)
    by_name = {m["short"]: m for m in msgs}
    assert set(by_name) == {"2t", "tp"}
    assert by_name["2t"]["level"] == 2
    assert by_name["2t"]["stepType"] == "instant"
    assert by_name["tp"]["stepType"] == "accum"
    assert by_name["tp"]["stepRange"] == "0-1"
    assert by_name["2t"]["centre"] in (7, "kwbc")
    np.testing.assert_allclose(by_name["2t"]["values"], ds["t2m"].values.ravel(), atol=0.02)
    np.testing.assert_allclose(by_name["tp"]["values"], ds["apcp"].values.ravel(), atol=0.02)


def test_lambert_temperature_roundtrip(tmp_path):
    ny, nx = 2, 3
    lat = np.array([[35.0, 35.1, 35.2], [35.3, 35.4, 35.5]])
    lon = np.array([[240.0, 240.2, 240.4], [240.0, 240.2, 240.4]])
    values = np.arange(ny * nx, dtype=np.float64).reshape(1, ny, nx) + 280.0
    ds = xr.Dataset(
        data_vars={
            "t": (
                ("time", "isobaric", "y", "x"),
                values[:, None, :, :],
                {"grib2_shortName": "TMP", "units": "K"},
            ),
            "crs": xr.DataArray(
                np.int32(1),
                attrs={
                    "grid_mapping_name": "lambert_conformal_conic",
                    "standard_parallel": (30.0, 60.0),
                    "longitude_of_central_meridian": -100.0,
                },
            ),
        },
        coords={
            "time": np.array(["2018-07-23T21:00:00"], dtype="datetime64[ns]"),
            "isobaric": ("isobaric", np.array([50000.0])),
            "y": np.arange(ny),
            "x": np.arange(nx),
            "latitude": (("y", "x"), lat),
            "longitude": (("y", "x"), lon),
            "forecast_hour": ("time", np.array([0.0])),
        },
        attrs={
            "forecast_reference_time": "2018-07-23T21:00:00Z",
            "DX": 3000.0,
            "DY": 3000.0,
        },
    )
    path = tmp_path / "lambert.grib2"
    skipped = write_grib2(ds, str(path))
    assert skipped == []
    msgs = _messages(path)
    assert len(msgs) == 1
    assert msgs[0]["short"] == "t"
    assert msgs[0]["level"] == 500
    np.testing.assert_allclose(msgs[0]["values"], values.ravel(), atol=0.05)


def test_unknown_short_name_is_reported(tmp_path):
    ds = _sample_dataset()
    ds["mystery"] = xr.DataArray(
        ds["t2m"].values,
        dims=("time", "y", "x"),
        attrs={"grib2_shortName": "NOT_A_PARAM"},
    )
    skipped = write_grib2(ds, str(tmp_path / "skip.grib2"))
    assert any("NOT_A_PARAM" in line for line in skipped)
    names = {m["short"] for m in _messages(tmp_path / "skip.grib2")}
    assert names == {"2t", "tp"}
