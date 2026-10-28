"""CALRH / CALDWP against the formulas in the Fortran."""

import numpy as np

from upp_py.constants import A2, A3, A4, PQ0, RHMIN
from upp_py.thermo import dewpoint, relative_humidity, saturation_specific_humidity


def test_saturation_matches_calrh_expression():
    t = np.array([280.0, 290.0, 300.0])
    p = np.array([85000.0, 100000.0, 50000.0])
    qs = saturation_specific_humidity(t, p)
    expected = PQ0 / p * np.exp(A2 * (t - A3) / (t - A4))
    np.testing.assert_allclose(qs, expected)


def test_rh_clips_and_rewrites_q():
    t = np.array([290.0])
    p = np.array([1.0e5])
    qs = saturation_specific_humidity(t, p)
    rh, q = relative_humidity(p, t, qs * 5, percent=False, adjust_q=True)
    assert rh[0] == 1.0
    np.testing.assert_allclose(q, qs)
    rh_lo, q_lo = relative_humidity(p, t, qs * 0, percent=True, adjust_q=True)
    np.testing.assert_allclose(rh_lo, RHMIN * 100)
    np.testing.assert_allclose(q_lo, RHMIN * qs)


def test_dewpoint_triple_point_and_not_above_t():
    # 6.11 hPa = 611 Pa = 0.611 centibar is the triple-point vapor pressure.
    # q from e = p q / (eps + (1-eps) q) inverted at p=1e5.
    from upp_py.constants import EPS_MASS, ONEPS

    e = 611.0
    p = 1.0e5
    q = EPS_MASS * e / (p - ONEPS * e)
    td = dewpoint(np.array([p]), np.array([q]), np.array([320.0]))
    assert abs(td[0] - 273.16) < 0.4
    td_cap = dewpoint(np.array([p]), np.array([q * 10]), np.array([280.0]))
    assert td_cap[0] <= 280.0 + 1e-6
