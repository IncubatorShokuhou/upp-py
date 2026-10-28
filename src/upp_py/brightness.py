"""Optional GOES ABI brightness temperature.

UPP v4.1 computes this in ``CALRAD_WCLOUD_newcrtm.f`` (available-field
275, GRIB pname ``BRTMP``). The sensors used for GOES-16 and GOES-17
are ``abi_g16`` and ``abi_g17``. Each posts ABI infrared channels 7-16
(ten fields, ``IGET`` 927-936 and 937-946). Nadir ``abi_gr`` is
``IGET`` 958-967. The forward model is CRTM.

That model is not part of upp-py. Install the separate package
``upp-crtm`` (import name ``upp_crtm``) and call
``upp_crtm.brightness_temperature``. This module only locates that
function. It does not ship CRTM source, coefficient files, or a
``TSK**4`` substitute.

Expected signature::

    brightness_temperature(state, *, satellite: str, channels: tuple[int, ...]) -> numpy.ndarray

``state`` is a :class:`upp_py.wrf.WrfState`. The return is brightness
temperature in kelvin with shape ``(len(channels), ny, nx)``. Raise
``RuntimeError`` when coefficients or the CRTM library are missing;
upp-py copies that message into the dataset ``skipped`` attribute and
does not write ``brtmp``.
"""

from __future__ import annotations

import importlib
from collections.abc import Callable, Sequence

import numpy as np

from upp_py.wrf import WrfState

# CALRAD_WCLOUD_newcrtm.f sensorlist entries for ABI.
ABI_SATELLITES = ("abi_g16", "abi_g17", "abi_gr")
# Infrared channels 7-16. CALRAD loops ixchan = 1..10 over these.
ABI_IR_CHANNELS = (7, 8, 9, 10, 11, 12, 13, 14, 15, 16)

VIS_GSD_SKIP_REASON = (
    "CALVIS_GSD (field 410) is called from MDLFLD.f only when imp_physics==28, "
    "and it needs the cosine of the solar zenith from ZENSUN.f. "
    "Surface visibility here is CALVIS (field 180)."
)

GFIP_NOTE = (
    "NCAR GFIP icing_algo in GFIP3.f is not ported. "
    "MDLFLD.f computes it for fields 450 and 480, and the isobaric GRIB write is commented out. "
    "The isobaric icing product is CALICING / TIPD (field 257)."
)

_PACKAGE = "upp_crtm"
_ENTRY = "brightness_temperature"

BrightnessFn = Callable[..., np.ndarray]


class BrtmpUnavailable(RuntimeError):
    """The optional CRTM package cannot produce BRTMP."""


def skip_reason() -> str:
    """Why ``BRTMP`` is omitted when ``upp_crtm`` is not installed."""
    return (
        "CALRAD_WCLOUD_newcrtm.f BRTMP (GOES ABI IR channels 7-16) is computed by "
        f"the separate package {_PACKAGE}.{_ENTRY}. It is not installed. "
        "upp-py does not include CRTM source or coefficient files."
    )


def load_brightness_temperature() -> BrightnessFn | None:
    """Return ``upp_crtm.brightness_temperature``, or ``None`` if absent."""
    try:
        module = importlib.import_module(_PACKAGE)
    except ImportError:
        return None
    fn = getattr(module, _ENTRY, None)
    if not callable(fn):
        return None
    return fn


def compute_brtmp(
    state: WrfState,
    *,
    satellite: str = "abi_g16",
    channels: Sequence[int] = ABI_IR_CHANNELS,
) -> np.ndarray:
    """Call the optional package.

    Parameters
    ----------
    state
        One wrfout time.
    satellite
        CRTM sensor id from ``CALRAD_WCLOUD_newcrtm.f``:
        ``abi_g16``, ``abi_g17``, or ``abi_gr``.
    channels
        ABI channel numbers. The default is 7-16.

    Returns
    -------
    numpy.ndarray
        Kelvin, shape ``(len(channels), ny, nx)``.

    Raises
    ------
    BrtmpUnavailable
        The package is missing, the return shape is wrong, or the
        package raised ``RuntimeError``.
    """
    fn = load_brightness_temperature()
    if fn is None:
        raise BrtmpUnavailable(skip_reason())
    channel_list = tuple(int(c) for c in channels)
    try:
        tb = np.asarray(
            fn(state, satellite=satellite, channels=channel_list),
            dtype=np.float64,
        )
    except RuntimeError as exc:
        raise BrtmpUnavailable(str(exc)) from exc
    expected = (len(channel_list), state.ny, state.nx)
    if tb.shape != expected:
        raise BrtmpUnavailable(
            f"{_PACKAGE}.{_ENTRY} returned shape {tb.shape}, expected {expected}"
        )
    return tb
