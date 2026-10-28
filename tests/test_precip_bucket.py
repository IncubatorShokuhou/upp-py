"""WRF rain-bucket totals and PREC_ACC interval fields."""

from datetime import datetime, timedelta

import numpy as np

from tests.test_interp_mslp import _tiny_state
from upp_py.config import PostConfig
from upp_py.post import dataset_from_state
from upp_py.precip import (
    bucket_total_mm,
    continuous_total_mm,
    window_accumulation_mm,
)


def _at(hour, rainc, *, i_rainc=None, bucket_mm=0.0, prec=None, rainnc=0.0):
    start = datetime(2018, 7, 23, 21)
    ny, nx = 4, 4
    kwargs = dict(
        valid_time=start + timedelta(hours=hour),
        start_time=start,
        forecast_hour=float(hour),
        rainc=np.full((ny, nx), rainc, dtype=np.float64),
        rainnc=np.full((ny, nx), rainnc, dtype=np.float64),
        bucket_mm=float(bucket_mm),
    )
    if i_rainc is not None:
        kwargs["i_rainc"] = np.full((ny, nx), i_rainc, dtype=np.float64)
        kwargs["i_rainnc"] = np.zeros((ny, nx))
    if prec is not None:
        kwargs["prec_acc_c"] = np.full((ny, nx), prec, dtype=np.float64)
        kwargs["prec_acc_nc"] = np.zeros((ny, nx))
        kwargs["prec_acc_dt"] = 3600.0
    return _tiny_state(**kwargs)


def test_continuous_total_without_counters_is_rainc_plus_rainnc():
    state = _at(2, rainc=3.0, rainnc=1.5)
    np.testing.assert_allclose(continuous_total_mm(state), 4.5)


def test_continuous_total_adds_bucket_mm_times_counter():
    # Remainder 2 mm after one dump of a 10 mm bucket -> 12 mm.
    state = _at(2, rainc=2.0, i_rainc=1, bucket_mm=10.0)
    np.testing.assert_allclose(continuous_total_mm(state), 12.0)


def test_ignored_counter_when_bucket_mm_is_not_positive():
    state = _at(2, rainc=2.0, i_rainc=3, bucket_mm=0.0)
    np.testing.assert_allclose(continuous_total_mm(state), 2.0)


def test_bucket_field_is_zero_at_forecast_hour_0():
    state = _at(0, rainc=0.0, prec=5.0)
    np.testing.assert_allclose(bucket_total_mm(state), 0.0)
    later = _at(1, rainc=4.0, prec=4.0)
    np.testing.assert_allclose(bucket_total_mm(later), 4.0)
    assert bucket_total_mm(_at(1, rainc=1.0, prec=None)) is None


def test_window_uses_counters_across_a_reset():
    series = [
        _at(1, rainc=6.0, i_rainc=0, bucket_mm=10.0),
        _at(2, rainc=2.0, i_rainc=1, bucket_mm=10.0),
    ]
    amounts, actual = window_accumulation_mm(series, hours=1.0)
    np.testing.assert_allclose(amounts[1], 6.0)
    assert actual[1] == 1.0
    # From the forecast start through hour 2 the continuous total is 12 mm.
    np.testing.assert_allclose(amounts[0], 6.0)
    full, _full_h = window_accumulation_mm(series, hours=5.0)
    np.testing.assert_allclose(full[1], 12.0)


def test_window_falls_back_to_prec_acc_when_the_raw_total_drops():
    # No I_RAINC. RAINC falls from 6 mm to 2 mm, which is a reset.
    # PREC_ACC at the later time holds the 6 mm that fell in that hour.
    series = [
        _at(1, rainc=6.0, prec=6.0),
        _at(2, rainc=2.0, prec=6.0),
    ]
    amounts, actual = window_accumulation_mm(series, hours=1.0)
    np.testing.assert_allclose(amounts[1], 6.0)
    assert actual[1] == 1.0


def test_dataset_writes_bucket_and_records_a_skip():
    state = _at(0, rainc=0.0, prec=3.0)
    ds = dataset_from_state(state, PostConfig(fields=("APCP", "APCP_BUCKET")))
    np.testing.assert_allclose(ds["apcp"].values, 0.0)
    np.testing.assert_allclose(ds["apcp_bucket"].values, 0.0)
    bare = _at(3, rainc=1.25)
    missing = dataset_from_state(bare, PostConfig(fields=("APCP", "APCP_BUCKET", "BRTMP")))
    assert "apcp_bucket" not in missing
    assert "APCP_BUCKET" in missing.attrs["skipped"]
    assert "BRTMP" in missing.attrs["skipped"]
    assert "upp_crtm" in missing.attrs["skipped"]
