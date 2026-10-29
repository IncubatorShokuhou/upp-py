"""Optional GOES ABI brightness temperature.

UPP v4.1 computes this in ``CALRAD_WCLOUD_newcrtm.f`` (available-field
275, GRIB pname ``BRTMP``). The sensors used for GOES-16 and GOES-17
are ``abi_g16`` and ``abi_g17``. Each posts ABI infrared channels 7-16
(ten fields, ``IGET`` 927-936 and 937-946). Nadir ``abi_gr`` is
``IGET`` 958-967. The forward model is CRTM.

That model is not part of upp-py. It lives in
https://github.com/IncubatorShokuhou/crtm-py (distribution ``crtm-py``,
import ``crtm_py``). This module passes a :class:`~upp_py.wrf.WrfState`
into :func:`crtm_py.brightness_temperature`, which builds the CALRAD
atmosphere, surface and geostationary zenith. upp-py does not ship CRTM
source or coefficient files, and it does not substitute ``TSK**4``.

Without the package, or without ``abi_*.SpcCoeff.bin`` and ODPS
``TauCoeff.bin``, the field is omitted and the reason is stored on the
dataset ``skipped`` attribute. A synthetic grey-gas result is not written
as ``BRTMP``.
"""

from __future__ import annotations

import importlib
from collections.abc import Sequence
from typing import Any

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

_PACKAGE = "crtm_py"
_ENTRY = "brightness_temperature"

class BrtmpUnavailable(RuntimeError):
    """The optional CRTM package cannot produce BRTMP."""


def skip_reason() -> str:
    """Why ``BRTMP`` is omitted when ``crtm_py`` is not installed."""
    return (
        "CALRAD_WCLOUD_newcrtm.f BRTMP (GOES ABI IR channels 7-16) is computed by "
        f"{_PACKAGE}.{_ENTRY} from https://github.com/IncubatorShokuhou/crtm-py. "
        "That package is not installed. upp-py does not include CRTM source or coefficient files."
    )


def coeff_skip_reason(satellite: str) -> str:
    """Why ``BRTMP`` is omitted when the ODPS coefficient files are absent."""
    return (
        f"crtm_py is installed but {satellite}.SpcCoeff.bin and an ODPS "
        f"{satellite}.TauCoeff.bin were not found. Set CRTM_PY_COEFF_DIR or run "
        "`crtm-py fetch`. A synthetic grey-gas brightness temperature is not written as BRTMP."
    )


def load_brightness_temperature() -> Any:
    """Return the ``crtm_py`` module, or ``None`` if it is not installed."""
    try:
        return importlib.import_module(_PACKAGE)
    except ImportError:
        return None


def _sea_mask(state: WrfState) -> np.ndarray:
    """UPP ``sm``: ``XLAND - 1`` (``INITPOST.F``). Else ``1 - LANDMASK``."""
    if state.xland is not None:
        return np.asarray(state.xland, dtype=np.float64) - 1.0
    if state.landmask is not None:
        return 1.0 - np.asarray(state.landmask, dtype=np.float64)
    return np.zeros((state.ny, state.nx), dtype=np.float64)


def call_kwargs(state: WrfState, *, satellite: str, channels: Sequence[int]) -> dict:
    """Arguments for :func:`crtm_py.brightness_temperature`.

    Vertical arrays stay in WRF bottom-up order. ``crtm_py`` flips them
    and applies the ``CALRAD`` absorber, ``EFFR`` and
    ``GEO_ZENITH_ANGLE`` steps. ``abi_g16`` uses the sub-satellite point
    (0 N, 75.2 W).
    """
    skin = state.tsk if state.tsk is not None else state.t[0]
    ivg = state.ivgtyp if state.ivgtyp is not None else np.ones((state.ny, state.nx), dtype=np.int64)
    return {
        "pressure_interface": state.p_int,
        "pressure": state.p_mid,
        "temperature": state.t,
        "specific_humidity": state.q,
        "sensor": satellite,
        "channels": tuple(int(c) for c in channels),
        "vertical_order": "bottom_up",
        "ozone": state.o3,
        "imp_physics": int(state.mp_physics),
        "latitude": state.lat,
        "longitude": state.lon,
        "skin_temperature": skin,
        "sea_mask": _sea_mask(state),
        "snow_fraction": state.snow_fraction if state.snow_fraction is not None else 0.0,
        "sea_ice": state.seaice if state.seaice is not None else 0.0,
        "vegetation_fraction": state.veg_fraction if state.veg_fraction is not None else 0.0,
        "u10": state.u10,
        "v10": state.v10,
        "ivgtyp": ivg,
        "ivegsrc": int(state.ivegsrc),
        "icu_physics": int(state.cu_physics),
        "hbot": state.hbot,
        "htop": state.htop,
        "cnvcfr": state.cnvcfr,
        "include_clouds": True,
        "include_scattering": True,
        "_qc": state.qc,
        "_qi": state.qi,
        "_qr": state.qr,
        "_qs": state.qs,
        "_qg": state.qg,
        "_qnr": state.qnr,
        "_qni": state.qni,
        "_f_rimef": state.f_rimef,
        "_snow_depth_mm": state.snow_depth_mm,
        "_smstot": state.smstot,
    }


def compute_brtmp(
    state: WrfState,
    *,
    satellite: str = "abi_g16",
    channels: Sequence[int] = ABI_IR_CHANNELS,
) -> np.ndarray:
    """GOES ABI brightness temperature in kelvin, shape ``(n_channels, ny, nx)``.

    Parameters
    ----------
    state
        One wrfout time, WRF bottom-up.
    satellite
        ``abi_g16``, ``abi_g17`` or ``abi_gr`` (``CALRAD`` sensor list).
    channels
        ABI band numbers. The default is 7-16.

    Raises
    ------
    BrtmpUnavailable
        ``crtm_py`` is missing, coefficient files are missing, the result
        is the synthetic grey gas, or the returned shape is wrong.
    """
    crtm_py = load_brightness_temperature()
    if crtm_py is None or not hasattr(crtm_py, _ENTRY):
        raise BrtmpUnavailable(skip_reason())
    root = crtm_py.coeff_path()
    if root is None:
        raise BrtmpUnavailable(coeff_skip_reason(satellite))
    kwargs = call_kwargs(state, satellite=satellite, channels=channels)
    hydro = crtm_py.Hydrometeors(
        qqw=kwargs.pop("_qc"),
        qqi=kwargs.pop("_qi"),
        qqr=kwargs.pop("_qr"),
        qqs=kwargs.pop("_qs"),
        qqg=kwargs.pop("_qg"),
        qqnr=kwargs.pop("_qnr"),
        qqni=kwargs.pop("_qni"),
        f_rimef=kwargs.pop("_f_rimef"),
    )
    skin = kwargs.pop("skin_temperature")
    snow_depth = kwargs.pop("_snow_depth_mm")
    soil = kwargs.pop("_smstot")
    surface = crtm_py.surface_from_upp(
        skin_temperature=skin,
        sea_mask=kwargs.pop("sea_mask"),
        ivgtyp=kwargs.pop("ivgtyp"),
        ivegsrc=kwargs.pop("ivegsrc"),
        snow_fraction=kwargs.pop("snow_fraction"),
        sea_ice=kwargs.pop("sea_ice"),
        vegetation_fraction=kwargs.pop("vegetation_fraction"),
        u10=kwargs.pop("u10"),
        v10=kwargs.pop("v10"),
        soil_moisture=np.nan if soil is None else soil,
        snow_depth=0.0 if snow_depth is None else snow_depth,
    )
    kwargs["hydrometeors"] = hydro
    kwargs["surface"] = surface
    kwargs["skin_temperature"] = skin
    try:
        result = crtm_py.brightness_temperature(**kwargs)
    except RuntimeError as exc:
        raise BrtmpUnavailable(str(exc)) from exc
    model = str(getattr(result, "optical_depth_model", ""))
    if model.startswith("synthetic"):
        raise BrtmpUnavailable(coeff_skip_reason(satellite))
    tb = np.asarray(result.brightness_temperature, dtype=np.float64)
    channel_list = kwargs["channels"]
    if tb.shape == (state.ny, state.nx, len(channel_list)):
        tb = np.moveaxis(tb, -1, 0)
    expected = (len(channel_list), state.ny, state.nx)
    if tb.shape != expected:
        raise BrtmpUnavailable(
            f"{_PACKAGE}.{_ENTRY} returned shape {tb.shape}, expected {expected}"
        )
    return tb
