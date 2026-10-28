"""Write a posted dataset to GRIB2 with the ECMWF eccodes Python bindings.

The optional extra ``grib`` installs ``eccodes`` (the wheel ships the C
library). Messages use NCEP as the originating centre (``kwbc``, centre
7) so short names match the tables wgrib2 uses for WRFPRS: ``prmsl``,
``tp``, ``refc``, ``vis``, ``tipd``. Packing is grid_simple at 24 bits.
It is not the decimal scale in ``params_grib2_tbl_new`` and it is not
bit-identical to ``unipost.exe``.

Horizontal grids are regular latitude-longitude or Lambert conformal,
taken from the dataset ``crs`` and ``latitude``/``longitude``. Lambert
longitudes are written on 0-360, which is what eccodes packs.

Accumulated precipitation uses ``stepType=accum`` and a ``stepRange``.
Run-total ``apcp`` spans forecast hour 0 to the valid hour. Bucket and
window fields use the length stored on the variable (``accum_length_hours``
or, for ``apcp_accum``, the per-time ``accum_hours_actual`` coordinate).
"""

from __future__ import annotations

import math
from datetime import datetime

import numpy as np
import xarray as xr

# wgrib2-style name on the DataArray -> eccodes shortName at centre kwbc.
_ECCODES_NAME = {
    "TMP": "t",
    "DPT": "dpt",
    "RH": "r",
    "HGT": "gh",
    "UGRD": "u",
    "VGRD": "v",
    "VVEL": "w",
    "SPFH": "q",
    "PRMSL": "prmsl",
    "MSLET": "mslet",
    "APCP": "tp",
    "ACPCP": "acpcp",
    "NCPCP": "ncpcp",
    "PRATE": "prate",
    "PWAT": "pwat",
    "CAPE": "cape",
    "CIN": "cin",
    "REFC": "refc",
    "REFD": "refd",
    "VIS": "vis",
    "TCDC": "tcc",
    "LCDC": "lcc",
    "MCDC": "mcc",
    "HCDC": "hcc",
    "HLCY": "hlcy",
    "TIPD": "tipd",
    "USTM": "ustm",
    "VSTM": "vstm",
    "ABSV": "absv",
    "PRES": "sp",
    "WEASD": "sdwe",
    "ASNOW": "asnow",
}

_HEIGHT_NAME = {
    ("TMP", "2 m"): "2t",
    ("DPT", "2 m"): "2d",
    ("RH", "2 m"): "2r",
    ("UGRD", "10 m"): "10u",
    ("VGRD", "10 m"): "10v",
}

_VAR_NAME = {
    "orog": "orog",
    "psfc": "sp",
}

_CLOUD_LEVEL = {
    "tcdc": "entireAtmosphere",
    "lcdc": "lowCloudLayer",
    "mcdc": "middleCloudLayer",
    "hcdc": "highCloudLayer",
}

# Parcel name -> (typeOfLevel, bottom, top). Pressure-from-ground values
# are Pa (GRIB surface type 108). 90 hPa, 180 hPa and 255 hPa match the
# UPP mixed-layer, best and most-unstable layer labels.
_PARCEL_LEVEL = {
    "surface": ("surface", None, None),
    "mixed_layer": ("pressureFromGroundLayer", 9000.0, 0.0),
    "best": ("pressureFromGroundLayer", 18000.0, 0.0),
    "most_unstable": ("pressureFromGroundLayer", 25500.0, 0.0),
}

_SKIP_VARS = {"crs"}


class Grib2Unavailable(RuntimeError):
    """eccodes is not installed."""


def _import_eccodes():
    try:
        from eccodes import (  # noqa: WPS433
            codes_grib_new_from_samples,
            codes_release,
            codes_set,
            codes_set_values,
            codes_write,
        )
    except ImportError as exc:  # pragma: no cover
        raise Grib2Unavailable(
            "GRIB2 output needs the eccodes package: pip install -e '.[grib]'"
        ) from exc
    return codes_grib_new_from_samples, codes_set, codes_set_values, codes_write, codes_release


def _as_pair(value) -> tuple[float, float]:
    if isinstance(value, (list, tuple, np.ndarray)):
        flat = np.asarray(value, dtype=np.float64).reshape(-1)
        if flat.size == 1:
            return float(flat[0]), float(flat[0])
        return float(flat[0]), float(flat[1])
    number = float(value)
    return number, number


def _lon360(lon: float) -> float:
    out = float(lon) % 360.0
    if out < 0.0:
        out += 360.0
    return out


def _ref_date(ds: xr.Dataset) -> tuple[int, int]:
    text = str(ds.attrs.get("forecast_reference_time", "")).replace("Z", "")
    if not text and "time" in ds.coords:
        raw = np.asarray(ds["time"].values).reshape(-1)[0]
        stamp = raw.astype("datetime64[s]").astype(datetime)
        return int(stamp.strftime("%Y%m%d")), int(stamp.strftime("%H%M"))
    if "T" in text:
        when = datetime.strptime(text[:19], "%Y-%m-%dT%H:%M:%S")
    else:
        when = datetime.strptime(text[:19], "%Y-%m-%d %H:%M:%S")
    return int(when.strftime("%Y%m%d")), int(when.strftime("%H%M"))


def _forecast_hours(ds: xr.Dataset) -> np.ndarray:
    ntime = int(ds.sizes["time"]) if "time" in ds.sizes else 1
    if "forecast_hour" in ds.coords:
        hours = np.asarray(ds["forecast_hour"].values, dtype=np.float64).reshape(-1)
        if hours.size == ntime:
            return hours
    if "forecast_hour" in ds.attrs:
        return np.full(ntime, float(ds.attrs["forecast_hour"]), dtype=np.float64)
    return np.zeros(ntime, dtype=np.float64)


def _grid_kwargs(ds: xr.Dataset) -> dict:
    lat = np.asarray(ds["latitude"].values, dtype=np.float64)
    lon = np.asarray(ds["longitude"].values, dtype=np.float64)
    if lat.ndim != 2:
        raise ValueError("GRIB2 write expects 2-D latitude and longitude")
    crs = ds["crs"].attrs if "crs" in ds else {}
    name = str(crs.get("grid_mapping_name", "latitude_longitude"))
    flip_j = lat[-1, 0] < lat[0, 0]
    flip_i = lon[0, -1] < lon[0, 0]
    if flip_j:
        lat = lat[::-1, :]
        lon = lon[::-1, :]
    if flip_i:
        lat = lat[:, ::-1]
        lon = lon[:, ::-1]
    ny, nx = lat.shape
    common = {
        "flip_j": flip_j,
        "flip_i": flip_i,
        "ny": ny,
        "nx": nx,
        "first_lat": float(lat[0, 0]),
        "first_lon": _lon360(float(lon[0, 0])),
    }
    if name == "lambert_conformal_conic":
        lat1, lat2 = _as_pair(crs.get("standard_parallel", (ds.attrs.get("TRUELAT1", 0.0), ds.attrs.get("TRUELAT2", 0.0))))
        lov = _lon360(float(crs.get("longitude_of_central_meridian", ds.attrs.get("STAND_LON", 0.0))))
        dx = float(ds.attrs.get("DX", crs.get("dx", 0.0)))
        dy = float(ds.attrs.get("DY", crs.get("dy", 0.0)))
        common.update(
            kind="lambert",
            latin1=lat1,
            latin2=lat2,
            lov=lov,
            lad=lat1,
            dx=dx,
            dy=dy,
        )
        return common
    dlon = float(lon[0, 1] - lon[0, 0]) if nx > 1 else 0.0
    dlat = float(lat[1, 0] - lat[0, 0]) if ny > 1 else 0.0
    common.update(
        kind="latlon",
        dlon=abs(dlon) if dlon != 0.0 else 0.0,
        dlat=abs(dlat) if dlat != 0.0 else 0.0,
        last_lat=float(lat[-1, -1]),
        last_lon=_lon360(float(lon[-1, -1])),
    )
    return common


def _orient(values: np.ndarray, grid: dict) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64)
    if grid["flip_j"]:
        out = out[::-1, :]
    if grid["flip_i"]:
        out = out[:, ::-1]
    return np.ascontiguousarray(out)


def _eccodes_short(var_name: str, short: str, height: str | None) -> str:
    if var_name in _VAR_NAME:
        return _VAR_NAME[var_name]
    if height and (short, height) in _HEIGHT_NAME:
        return _HEIGHT_NAME[(short, height)]
    if short not in _ECCODES_NAME:
        raise KeyError(short)
    return _ECCODES_NAME[short]


def _apply_grid(codes_set, gid, grid: dict) -> None:
    if grid["kind"] == "lambert":
        codes_set(gid, "gridType", "lambert")
        codes_set(gid, "Nx", int(grid["nx"]))
        codes_set(gid, "Ny", int(grid["ny"]))
        codes_set(gid, "latitudeOfFirstGridPointInDegrees", grid["first_lat"])
        codes_set(gid, "longitudeOfFirstGridPointInDegrees", grid["first_lon"])
        codes_set(gid, "LoVInDegrees", grid["lov"])
        codes_set(gid, "LaDInDegrees", grid["lad"])
        codes_set(gid, "Latin1InDegrees", grid["latin1"])
        codes_set(gid, "Latin2InDegrees", grid["latin2"])
        codes_set(gid, "DxInMetres", float(grid["dx"]))
        codes_set(gid, "DyInMetres", float(grid["dy"]))
        return
    codes_set(gid, "gridType", "regular_ll")
    codes_set(gid, "Ni", int(grid["nx"]))
    codes_set(gid, "Nj", int(grid["ny"]))
    codes_set(gid, "latitudeOfFirstGridPointInDegrees", grid["first_lat"])
    codes_set(gid, "longitudeOfFirstGridPointInDegrees", grid["first_lon"])
    codes_set(gid, "latitudeOfLastGridPointInDegrees", grid["last_lat"])
    codes_set(gid, "longitudeOfLastGridPointInDegrees", grid["last_lon"])
    codes_set(gid, "iDirectionIncrementInDegrees", float(grid["dlon"]) if grid["nx"] > 1 else 0.0)
    codes_set(gid, "jDirectionIncrementInDegrees", float(grid["dlat"]) if grid["ny"] > 1 else 0.0)


def _step_for(da: xr.DataArray, fhr: float, time_index: int, ds: xr.Dataset) -> tuple[str, str]:
    """Return ``(stepType, stepRange)`` in hours."""
    stat = str(da.attrs.get("grib_stat", "instant"))
    if stat != "accumulation":
        label = _hour_label(fhr)
        return "instant", label
    length = da.attrs.get("accum_length_hours")
    if da.name == "apcp_accum" and "accum_hours_actual" in ds.coords:
        length = float(np.asarray(ds["accum_hours_actual"].values).reshape(-1)[time_index])
        if math.isnan(length):
            length = 0.0
    if length is None:
        start = 0.0
    else:
        start = max(0.0, float(fhr) - float(length))
    if abs(float(fhr) - start) < 1.0e-8:
        return "instant", _hour_label(fhr)
    return "accum", f"{_hour_label(start)}-{_hour_label(fhr)}"


def _hour_label(hours: float) -> str:
    if abs(hours - round(hours)) < 1.0e-6:
        return str(int(round(hours)))
    return f"{hours:.4g}"


def _level_keys(var_name: str, da: xr.DataArray, sel: dict) -> dict:
    """eccodes level keys for one horizontal slice."""
    height = da.attrs.get("height")
    if height == "2 m":
        return {"typeOfLevel": "heightAboveGround", "level": 2}
    if height == "10 m":
        return {"typeOfLevel": "heightAboveGround", "level": 10}
    if var_name in _CLOUD_LEVEL:
        return {"typeOfLevel": _CLOUD_LEVEL[var_name]}
    if "isobaric" in sel:
        pa = float(da["isobaric"].values[sel["isobaric"]])
        return {"typeOfLevel": "isobaricInhPa", "level": int(round(pa / 100.0))}
    if "height_agl" in sel:
        metres = float(da["height_agl"].values[sel["height_agl"]])
        return {"typeOfLevel": "heightAboveGround", "level": int(round(metres))}
    if "hel_depth" in sel:
        depth = float(da["hel_depth"].values[sel["hel_depth"]])
        return {"typeOfLevel": "heightAboveGroundLayer", "topLevel": 0, "bottomLevel": depth}
    if "parcel" in sel:
        name = str(da["parcel"].values[sel["parcel"]])
        kind, bottom, top = _PARCEL_LEVEL.get(name, ("surface", None, None))
        keys = {"typeOfLevel": kind}
        if bottom is not None:
            keys["bottomLevel"] = bottom
            keys["topLevel"] = top
        return keys
    if "model_level" in sel:
        return {"typeOfLevel": "hybrid", "level": int(sel["model_level"]) + 1}
    if "dust_bin" in sel:
        return {"typeOfLevel": "hybrid", "level": 1}
    return {"typeOfLevel": "surface", "level": 0}


def _iter_slices(da: xr.DataArray, time_index: int):
    data = da.isel(time=time_index) if "time" in da.dims else da
    extra = [d for d in data.dims if d not in ("y", "x")]
    if not extra:
        yield {}, np.asarray(data.values)
        return
    ranges = [range(int(data.sizes[d])) for d in extra]
    # Local import keeps the module import light.
    import itertools

    for idx in itertools.product(*ranges):
        sel = dict(zip(extra, idx))
        yield sel, np.asarray(data.isel(sel).values)


def write_grib2(ds: xr.Dataset, path: str) -> list[str]:
    """Write ``ds`` to a GRIB2 file.

    Parameters
    ----------
    ds
        Dataset returned by :func:`upp_py.post_wrfout`. Every variable
        with a ``grib2_shortName`` that eccodes can encode is written.
        ``crs`` is skipped.
    path
        Output path. An existing file is replaced.

    Returns
    -------
    list of str
        Variables or slices that were not written, with the reason.
        Model-level aerosols and any short name missing from the kwbc
        tables land here. The NetCDF file is unaffected.

    Raises
    ------
    Grib2Unavailable
        The ``eccodes`` package is not installed.
    """
    new_from_samples, codes_set, codes_set_values, codes_write, codes_release = _import_eccodes()
    grid = _grid_kwargs(ds)
    data_date, data_time = _ref_date(ds)
    hours = _forecast_hours(ds)
    ntime = int(hours.size)
    skipped: list[str] = []
    with open(path, "wb") as handle:
        for name, da in ds.data_vars.items():
            if name in _SKIP_VARS or "grib2_shortName" not in da.attrs:
                continue
            if "y" not in da.dims or "x" not in da.dims:
                skipped.append(f"{name}: no horizontal grid")
                continue
            short = str(da.attrs["grib2_shortName"])
            try:
                ec_name = _eccodes_short(name, short, da.attrs.get("height"))
            except KeyError:
                skipped.append(f"{name}: no eccodes short name for {short}")
                continue
            for it in range(ntime):
                fhr = float(hours[it])
                try:
                    step_type, step_range = _step_for(da, fhr, it, ds)
                except Exception as exc:  # pragma: no cover
                    skipped.append(f"{name}: step {exc}")
                    continue
                for sel, values in _iter_slices(da, it):
                    gid = new_from_samples("GRIB2")
                    try:
                        codes_set(gid, "centre", "kwbc")
                        codes_set(gid, "shortName", ec_name)
                        _apply_grid(codes_set, gid, grid)
                        for key, val in _level_keys(name, da, sel).items():
                            codes_set(gid, key, val)
                        codes_set(gid, "dataDate", int(data_date))
                        codes_set(gid, "dataTime", int(data_time))
                        codes_set(gid, "significanceOfReferenceTime", 1)
                        codes_set(gid, "stepType", step_type)
                        if step_type == "accum":
                            codes_set(gid, "stepRange", step_range)
                        else:
                            codes_set(gid, "step", step_range)
                        field = _orient(values, grid)
                        flat = field.ravel(order="C")
                        if not np.isfinite(flat).all():
                            missing = 9.999e20
                            codes_set(gid, "bitmapPresent", 1)
                            flat = np.where(np.isfinite(flat), flat, missing)
                        codes_set(gid, "bitsPerValue", 24)
                        codes_set_values(gid, flat.astype(np.float64))
                        codes_write(gid, handle)
                    except Exception as exc:
                        where = ",".join(f"{k}={v}" for k, v in sel.items()) or "slice"
                        skipped.append(f"{name} t={it} {where}: {exc}")
                    finally:
                        codes_release(gid)
    return skipped
