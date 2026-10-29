"""Accumulated precipitation, including WRF buckets.

``INITPOST.F`` reads ``RAINC`` and ``RAINNC`` (mm) into ``CUPREC`` and
``ANCPRC`` after multiplying by 0.001. ``SURFCE.f`` writes accumulated
total precipitation (available-field index 87) as ``ACPREC * 1000``,
which is ``RAINC + RAINNC`` in mm, equal to kg m-2.

The same routine's comment (MEB, 2002-06-17) says accumulated fields are
assumed never to empty their buckets. WRF does empty them when
``bucket_mm > 0``: the remainder stays in ``RAINC``/``RAINNC`` and
``I_RAINC``/``I_RAINNC`` counts how many times ``bucket_mm`` was
subtracted. The continuous total below adds those counts back. Files
without the counters keep the UPP value ``RAINC + RAINNC``.

Interval buckets are a different pair of arrays. ``INITPOST.F`` reads
``PREC_ACC_C`` and ``PREC_ACC_NC`` (mm over ``PREC_ACC_DT`` seconds) and
sums them into ``PCP_BUCKET``. ``SURFCE.f`` available-field 434 writes
that sum, except at forecast hour 0, where the field is set to 0.
Convective and grid-scale buckets are fields 435 and 436 and use the
same forecast-hour-0 reset. ``SNOW_ACC_NC`` (field 437) is not reset
at hour 0.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np

from upp_py.wrf import WrfState

# Negative differences smaller than this are roundoff, not a bucket reset.
_RESET_TOL_MM = 1.0e-4


def precipitation_rate(state: WrfState) -> np.ndarray:
    """Instantaneous total precipitation rate (kg m-2 s-1).

    ``SURFCE.f`` available-field 167. ``PREC`` is ``(RAINCV+RAINNCV)*0.001``
    metres per model step and ``DTQ2`` is the model ``DT``, so the rate is
    ``PREC / DT * 1000``.
    """
    dt = float(state.dt)
    if dt <= 0.0:
        return np.full(state.rainc.shape, np.nan, dtype=np.float64)
    return (np.asarray(state.raincv, dtype=np.float64) + np.asarray(state.rainncv, dtype=np.float64)) * 0.001 / dt * 1000.0


def _plus_counter(amount: np.ndarray, counter: np.ndarray | None, bucket_mm: float) -> np.ndarray:
    total = np.asarray(amount, dtype=np.float64).copy()
    if bucket_mm > 0.0 and counter is not None:
        total = total + float(bucket_mm) * np.asarray(counter, dtype=np.float64)
    return total


def continuous_total_mm(state: WrfState) -> np.ndarray:
    """Run-total precipitation (mm) since the integration start.

    Equals ``RAINC + RAINNC`` when ``bucket_mm`` is not positive or the
    ``I_RAINC``/``I_RAINNC`` counters are absent, which is what
    ``SURFCE.f`` writes. Otherwise adds ``bucket_mm * (I_RAINC + I_RAINNC)``.
    """
    rainc = _plus_counter(state.rainc, state.i_rainc, state.bucket_mm)
    rainnc = _plus_counter(state.rainnc, state.i_rainnc, state.bucket_mm)
    return rainc + rainnc


def continuous_convective_mm(state: WrfState) -> np.ndarray:
    """Convective run total (mm), ``SURFCE.f`` field 33, plus ``I_RAINC``."""
    return _plus_counter(state.rainc, state.i_rainc, state.bucket_mm)


def continuous_gridscale_mm(state: WrfState) -> np.ndarray:
    """Grid-scale run total (mm), ``SURFCE.f`` field 34, plus ``I_RAINNC``."""
    return _plus_counter(state.rainnc, state.i_rainnc, state.bucket_mm)


def _zero_at_analysis(state: WrfState, amount: np.ndarray) -> np.ndarray:
    """``SURFCE.f`` sets bucket precip to 0 when ``IFHR == 0``."""
    if abs(float(state.forecast_hour)) < 1.0e-8:
        return np.zeros(state.rainc.shape, dtype=np.float64)
    return amount


def bucket_total_mm(state: WrfState) -> np.ndarray | None:
    """``PREC_ACC_C + PREC_ACC_NC`` (mm), or 0 at forecast hour 0.

    Returns ``None`` when neither accumulator is in the wrfout.
    ``SURFCE.f`` available-field 434.
    """
    if state.prec_acc_c is None and state.prec_acc_nc is None:
        return None
    c = 0.0 if state.prec_acc_c is None else np.asarray(state.prec_acc_c, dtype=np.float64)
    n = 0.0 if state.prec_acc_nc is None else np.asarray(state.prec_acc_nc, dtype=np.float64)
    return _zero_at_analysis(state, c + n)


def bucket_convective_mm(state: WrfState) -> np.ndarray | None:
    """``PREC_ACC_C`` (mm), 0 at forecast hour 0. Field 435."""
    if state.prec_acc_c is None:
        return None
    return _zero_at_analysis(state, np.asarray(state.prec_acc_c, dtype=np.float64))


def bucket_gridscale_mm(state: WrfState) -> np.ndarray | None:
    """``PREC_ACC_NC`` (mm), 0 at forecast hour 0. Field 436."""
    if state.prec_acc_nc is None:
        return None
    return _zero_at_analysis(state, np.asarray(state.prec_acc_nc, dtype=np.float64))


def bucket_snow_mm(state: WrfState) -> np.ndarray | None:
    """``SNOW_ACC_NC`` (mm). Field 437 does not zero the field at hour 0."""
    if state.snow_acc_nc is None:
        return None
    return np.asarray(state.snow_acc_nc, dtype=np.float64)


def bucket_graupel_mm(state: WrfState) -> np.ndarray | None:
    """``GRAUP_ACC_NC`` (mm). ``SURFCE.f`` available-field 775."""
    if state.graup_acc_nc is None:
        return None
    return np.asarray(state.graup_acc_nc, dtype=np.float64)


def _step_increment(
    previous: np.ndarray,
    current: np.ndarray,
    bucket_now: np.ndarray | None,
) -> np.ndarray:
    """Precipitation between two history times (mm).

    A drop in the continuous total means ``RAINC``/``RAINNC`` were reset
    and the counters were not in the file. On those points the PREC_ACC
    bucket at the later time is the accumulation since the last reset.
    """
    diff = np.asarray(current, dtype=np.float64) - np.asarray(previous, dtype=np.float64)
    if bucket_now is None:
        return np.where(diff >= -_RESET_TOL_MM, diff, np.nan)
    bucket = np.asarray(bucket_now, dtype=np.float64)
    return np.where(diff >= -_RESET_TOL_MM, diff, bucket)


def window_accumulation_mm(
    states: list[WrfState],
    hours: float,
) -> tuple[list[np.ndarray], list[float]]:
    """Trailing accumulation ending at each state.

    Parameters
    ----------
    states
        One wrfout time each, any order. Sorted by valid time here.
    hours
        Requested window length. The start is clamped to the forecast
        start, as ``SURFCE.f`` clamps ``ID(18)`` to 0.

    Returns
    -------
    amounts, actual_hours
        ``amounts[i]`` is mm (kg m-2) ending at the sorted state ``i``.
        ``actual_hours[i]`` is the length actually covered. It is shorter
        than ``hours`` when the window is clamped to the forecast start,
        and longer when the history times do not land on the requested
        start (the latest history time at or before that start is used).
        Points that cannot be closed are NaN.
    """
    if hours < 0.0:
        raise ValueError("accum hours must be >= 0")
    order = sorted(range(len(states)), key=lambda i: states[i].valid_time)
    ordered = [states[i] for i in order]
    totals = [continuous_total_mm(s) for s in ordered]
    buckets = [bucket_total_mm(s) for s in ordered]
    start: datetime = ordered[0].start_time
    steps: list[np.ndarray | None] = [None]
    for i in range(1, len(ordered)):
        steps.append(_step_increment(totals[i - 1], totals[i], buckets[i]))

    amounts: list[np.ndarray] = []
    actual: list[float] = []
    for i, state in enumerate(ordered):
        target = state.valid_time - timedelta(hours=float(hours))
        if target <= start:
            acc = np.array(totals[0], dtype=np.float64, copy=True)
            for k in range(1, i + 1):
                acc = acc + steps[k]
            amounts.append(acc)
            actual.append(max((state.valid_time - start).total_seconds() / 3600.0, 0.0))
            continue
        j = None
        for k in range(i, -1, -1):
            if ordered[k].valid_time <= target:
                j = k
                break
        if j is None:
            amounts.append(np.full(state.rainc.shape, np.nan, dtype=np.float64))
            actual.append(float("nan"))
            continue
        acc = np.zeros(state.rainc.shape, dtype=np.float64)
        for k in range(j + 1, i + 1):
            acc = acc + steps[k]
        amounts.append(acc)
        actual.append((state.valid_time - ordered[j].valid_time).total_seconds() / 3600.0)
    # Restore the caller's order.
    out_amounts: list[np.ndarray | None] = [None] * len(states)
    out_hours: list[float] = [0.0] * len(states)
    for pos, src in enumerate(order):
        out_amounts[src] = amounts[pos]
        out_hours[src] = actual[pos]
    return out_amounts, out_hours  # type: ignore[return-value]
