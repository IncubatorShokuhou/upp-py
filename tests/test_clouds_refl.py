"""Neighborhood cloud fraction and a zero-hydrometeor reflectivity floor."""

import numpy as np

from upp_py.clouds import _box_fraction
from upp_py.reflectivity import reflectivity
from tests.test_interp_mslp import _tiny_state


def test_box_fraction_full_and_edge():
    mask = np.zeros((5, 5), dtype=bool)
    mask[2, 2] = True
    frac = _box_fraction(mask, 0)
    np.testing.assert_allclose(frac[2, 2], 1.0)
    np.testing.assert_allclose(frac[0, 0], 0.0)
    frac1 = _box_fraction(mask, 1)
    # 3x3 window around the center contains 1 cloudy point.
    np.testing.assert_allclose(frac1[2, 2], 1.0 / 9.0)
    # Corner window is clipped: radius 1 at (0, 0) sees a 2x2, none cloudy.
    np.testing.assert_allclose(frac1[0, 0], 0.0)


def test_clear_air_composite_near_minus_10_dbz():
    state = _tiny_state()
    refl = reflectivity(state)
    # Raw Stoelinga column max still sees the 0.1 mm^6/m^3 floor.
    np.testing.assert_allclose(refl.diagnosed, -10.0, atol=0.05)
    # REFD writes clear air as 0 dBZ.
    np.testing.assert_allclose(refl.agl_1km, 0.0)
    np.testing.assert_allclose(refl.agl_4km, 0.0)
