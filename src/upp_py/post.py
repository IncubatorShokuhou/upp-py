"""Assemble an xarray dataset from a wrfout file.

This is the driver UPP splits across ``INITPOST``, ``MDLFLD``, ``MDL2P``,
``SURFCE``, ``CLDRAD`` and ``MISCLN``. Only the WRF-ARW core products
are filled; each variable's ``grib2_shortName`` attribute is the name
wgrib2 prints for the same quantity.
"""

from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import xarray as xr

from upp_py.cape import cape_best, cape_mixed_layer, cape_most_unstable, cape_surface
from upp_py.clouds import cloud_fraction
from upp_py.config import PostConfig
from upp_py.constants import GI, SPL_PA
from upp_py.helicity import storm_relative_helicity
from upp_py.interp import interpolate_pressure_levels
from upp_py.mslp import membrane_slp, shuell_slp
from upp_py.reflectivity import reflectivity
from upp_py.thermo import dewpoint, relative_humidity
from upp_py.vorticity import absolute_vorticity
from upp_py.wrf import WrfState, open_wrfout


def _time_coord(when: datetime) -> np.ndarray:
    if when.tzinfo is None:
        when = when.replace(tzinfo=timezone.utc)
    return np.array([np.datetime64(when.replace(tzinfo=None))], dtype="datetime64[ns]")


def _var(data, dims, short, long_name, units, **extra):
    attrs = {"grib2_shortName": short, "long_name": long_name, "units": units}
    attrs.update(extra)
    arr = np.asarray(data)
    if "time" in dims and arr.shape[0] != 1:
        arr = arr[None, ...] if arr.ndim == len(dims) - 1 else arr
    # data arrays are un-timed; add the time axis here.
    if dims[0] == "time" and (arr.ndim == len(dims) - 1):
        arr = arr[None, ...]
    return xr.DataArray(arr, dims=dims, attrs=attrs)


def precipitable_water(state: WrfState) -> np.ndarray:
    """Column water vapor (kg/m^2), ``CALPW.f`` with ``IDECID <= 1``.

    ``PW = sum q dp / g`` using interface pressure thicknesses. Cloud
    water is a separate column integral in that routine.
    """
    pw = np.zeros(state.terrain.shape, dtype=np.float64)
    for k in range(state.nz):
        dp = state.p_int[k] - state.p_int[k + 1]
        pw = pw + state.q[k] * np.maximum(dp, 0.0) * GI
    return np.maximum(pw, 0.0)


def _projection(state: WrfState) -> xr.DataArray:
    """CF grid mapping for Lambert, Mercator or polar stereographic."""
    proj = int(state.map_proj)
    attrs: dict = {"earth_radius": 6370000.0}
    if proj == 1:
        attrs.update(
            {
                "grid_mapping_name": "lambert_conformal_conic",
                "standard_parallel": (state.truelat1, state.truelat2),
                "longitude_of_central_meridian": state.stand_lon,
                "latitude_of_projection_origin": state.cen_lat,
            }
        )
    elif proj == 2:
        attrs.update(
            {
                "grid_mapping_name": "polar_stereographic",
                "straight_vertical_longitude_from_pole": state.stand_lon,
                "latitude_of_projection_origin": state.cen_lat,
                "standard_parallel": state.truelat1,
            }
        )
    elif proj == 3:
        attrs.update(
            {
                "grid_mapping_name": "mercator",
                "longitude_of_projection_origin": state.stand_lon,
                "standard_parallel": state.truelat1,
            }
        )
    else:
        attrs["grid_mapping_name"] = "latitude_longitude"
    attrs["wrf_map_proj"] = proj
    return xr.DataArray(np.int32(proj), attrs=attrs)


def post_wrfout(path: str, config: PostConfig | None = None, time_index: int = 0) -> xr.Dataset:
    """Post one wrfout time to an :class:`xarray.Dataset`.

    Parameters
    ----------
    path
        WRF-ARW history file.
    config
        Field list. ``None`` posts the default core set on the standard
        UPP pressure levels.
    time_index
        Which time record to read.
    """
    cfg = config or PostConfig()
    state = open_wrfout(path, time_index=time_index)
    return dataset_from_state(state, cfg)


def dataset_from_state(state: WrfState, cfg: PostConfig) -> xr.Dataset:
    """Build the output dataset from an already-read :class:`WrfState`."""
    ny, nx = state.ny, state.nx
    y = np.arange(ny)
    x = np.arange(nx)
    data_vars: dict[str, xr.DataArray] = {}
    coords = {
        "time": _time_coord(state.valid_time),
        "y": y,
        "x": x,
        "latitude": (("y", "x"), state.lat, {"units": "degrees_north", "standard_name": "latitude"}),
        "longitude": (("y", "x"), state.lon, {"units": "degrees_east", "standard_name": "longitude"}),
    }

    # Isobaric block. MSLP needs the full standard set even if the user
    # only asked for a few output levels, because MEMSLP walks down it.
    need_iso = any(
        cfg.wants(n) for n in ("HGT", "TMP", "DPT", "RH", "SPFH", "UGRD", "VGRD", "VVEL", "ABSV", "PRMSL", "MSLET")
    )
    levels_out = None
    if need_iso:
        slp_levels = list(dict.fromkeys(list(SPL_PA) + list(cfg.pressure_pa)))
        computed = {lv.pressure: lv for lv in interpolate_pressure_levels(state, slp_levels)}
        wanted = sorted(cfg.pressure_pa)
        levels_out = [computed[p] for p in wanted]
        coords["isobaric"] = (
            "isobaric",
            np.array(wanted, dtype=np.float64),
            {"units": "Pa", "positive": "down", "long_name": "pressure"},
        )

    if cfg.wants("PRMSL"):
        # post_avblfldidx 105 / NGMSLP.f. The grib2 short name is PRMSL.
        data_vars["prmsl"] = _var(
            shuell_slp(state), ("time", "y", "x"), "PRMSL", "Shuell mean sea level pressure", "Pa"
        )
    if cfg.wants("MSLET"):
        # post_avblfldidx 23 / MEMSLP in SLP_new.f.
        data_vars["mslet"] = _var(
            membrane_slp(state, [computed[p] for p in SPL_PA]),
            ("time", "y", "x"),
            "MSLET",
            "Mesinger membrane sea level pressure",
            "Pa",
        )

    def _stack(getter):
        return np.stack([getter(lv) for lv in levels_out], axis=0)

    if levels_out is not None:
        if cfg.wants("HGT"):
            data_vars["gh"] = _var(_stack(lambda lv: lv.height), ("time", "isobaric", "y", "x"), "HGT", "geopotential height", "m")
        if cfg.wants("TMP"):
            data_vars["t"] = _var(
                _stack(lambda lv: lv.temperature), ("time", "isobaric", "y", "x"), "TMP", "temperature", "K"
            )
        if cfg.wants("SPFH"):
            data_vars["q"] = _var(
                _stack(lambda lv: lv.specific_humidity),
                ("time", "isobaric", "y", "x"),
                "SPFH",
                "specific humidity",
                "kg kg-1",
            )
        if cfg.wants("UGRD"):
            data_vars["u"] = _var(_stack(lambda lv: lv.u), ("time", "isobaric", "y", "x"), "UGRD", "grid-relative u wind", "m s-1")
        if cfg.wants("VGRD"):
            data_vars["v"] = _var(_stack(lambda lv: lv.v), ("time", "isobaric", "y", "x"), "VGRD", "grid-relative v wind", "m s-1")
        if cfg.wants("VVEL"):
            data_vars["w"] = _var(
                _stack(lambda lv: lv.omega), ("time", "isobaric", "y", "x"), "VVEL", "vertical velocity (omega)", "Pa s-1"
            )
        if cfg.wants("RH") or cfg.wants("DPT"):
            rh = []
            td = []
            for lv in levels_out:
                rh_pct, _q = relative_humidity(lv.pressure, lv.temperature, lv.specific_humidity, percent=True, adjust_q=True)
                # q was already RH-clamped inside the interpolator; do not clip twice into the stored q.
                rh.append(relative_humidity(lv.pressure, lv.temperature, lv.specific_humidity, percent=True))
                td.append(dewpoint(lv.pressure, lv.specific_humidity, lv.temperature))
            if cfg.wants("RH"):
                data_vars["r"] = _var(np.stack(rh, 0), ("time", "isobaric", "y", "x"), "RH", "relative humidity", "%")
            if cfg.wants("DPT"):
                data_vars["dpt"] = _var(np.stack(td, 0), ("time", "isobaric", "y", "x"), "DPT", "dewpoint temperature", "K")
        if cfg.wants("ABSV"):
            absv = []
            for lv in levels_out:
                absv.append(
                    absolute_vorticity(lv.u, lv.v, state.lat, state.dx, state.dy, state.mapfac)
                )
            data_vars["absv"] = _var(
                np.stack(absv, 0), ("time", "isobaric", "y", "x"), "ABSV", "absolute vorticity", "s-1"
            )

    # Shelter and anemometer levels. Winds stay grid-relative, as UPP
    # writes U10/V10 without a Lambert rotation.
    if cfg.wants("TMP2M"):
        data_vars["t2m"] = _var(state.t2, ("time", "y", "x"), "TMP", "2 m temperature", "K", height="2 m")
    if cfg.wants("SPFH2M"):
        data_vars["q2m"] = _var(state.q2, ("time", "y", "x"), "SPFH", "2 m specific humidity", "kg kg-1", height="2 m")
    if cfg.wants("RH2M"):
        rh2 = relative_humidity(state.p_shelter, state.t2, state.q2, percent=True)
        data_vars["rh2m"] = _var(rh2, ("time", "y", "x"), "RH", "2 m relative humidity", "%", height="2 m")
    if cfg.wants("DPT2M"):
        td2 = dewpoint(state.p_shelter, state.q2, state.t2)
        data_vars["dpt2m"] = _var(td2, ("time", "y", "x"), "DPT", "2 m dewpoint temperature", "K", height="2 m")
    if cfg.wants("UGRD10M"):
        data_vars["u10"] = _var(state.u10, ("time", "y", "x"), "UGRD", "10 m grid-relative u wind", "m s-1", height="10 m")
    if cfg.wants("VGRD10M"):
        data_vars["v10"] = _var(state.v10, ("time", "y", "x"), "VGRD", "10 m grid-relative v wind", "m s-1", height="10 m")

    if cfg.wants("APCP"):
        # SURFCE.f: ACPREC = (RAINC+RAINNC)*0.001 m, GRIB multiplies by 1000 -> mm = kg/m2.
        apcp = state.rainc + state.rainnc
        data_vars["apcp"] = _var(apcp, ("time", "y", "x"), "APCP", "total accumulated precipitation", "kg m-2")
    if cfg.wants("PRATE"):
        dt = state.dt if state.dt > 0 else np.nan
        # PREC (m/step) / DT * 1000 = kg m-2 s-1. RAINCV+RAINNCV are mm/step.
        prate = (state.raincv + state.rainncv) / dt / 1000.0 * 1000.0
        # mm/step / dt = mm/s = kg/m2/s. The *1000/1000 cancels; keep it
        # readable as the SURFCE factor: (mm * 0.001) / dt * 1000.
        prate = (state.raincv + state.rainncv) * 0.001 / dt * 1000.0
        data_vars["prate"] = _var(prate, ("time", "y", "x"), "PRATE", "instantaneous precipitation rate", "kg m-2 s-1")

    if cfg.wants("PWAT"):
        data_vars["pwat"] = _var(precipitable_water(state), ("time", "y", "x"), "PWAT", "precipitable water", "kg m-2")

    need_refl = cfg.wants("REFC") or cfg.wants("REFD")
    if need_refl:
        refl = reflectivity(state)
        if cfg.wants("REFC"):
            data_vars["refc"] = _var(refl.composite, ("time", "y", "x"), "REFC", "composite reflectivity", "dBZ")
        if cfg.wants("REFD"):
            agl = []
            names = []
            for depth in cfg.reflectivity_agl_m:
                if abs(depth - 1000.0) < 1:
                    agl.append(refl.agl_1km)
                elif abs(depth - 4000.0) < 1:
                    agl.append(refl.agl_4km)
                else:
                    raise ValueError("REFD heights other than 1000 m and 4000 m are not implemented")
                names.append(depth)
            coords["height_agl"] = ("height_agl", np.array(names, dtype=np.float64), {"units": "m", "long_name": "height above ground"})
            data_vars["refd"] = _var(np.stack(agl, 0), ("time", "height_agl", "y", "x"), "REFD", "reflectivity above ground", "dBZ")

    need_cape = any(cfg.wants(n) for n in ("CAPE", "CIN", "LCL", "LFC"))
    if need_cape:
        parcels = {
            "surface": cape_surface(state),
            "mixed_layer": cape_mixed_layer(state),
            "best": cape_best(state),
            "most_unstable": cape_most_unstable(state),
        }
        # GRIB level tags: surface, 90-0 hPa (mixed), 180-0 hPa (best), 255-0 hPa (most unstable label).
        coords["parcel"] = (
            "parcel",
            np.array(list(parcels.keys()), dtype=object),
            {
                "grib_level": "surface | 90-0 hPa | 180-0 hPa | 255-0 hPa (search is 300 hPa)",
            },
        )
        if cfg.wants("CAPE"):
            data_vars["cape"] = _var(
                np.stack([parcels[k].cape for k in parcels], 0),
                ("time", "parcel", "y", "x"),
                "CAPE",
                "convective available potential energy",
                "J kg-1",
            )
        if cfg.wants("CIN"):
            data_vars["cin"] = _var(
                np.stack([parcels[k].cin for k in parcels], 0),
                ("time", "parcel", "y", "x"),
                "CIN",
                "convective inhibition",
                "J kg-1",
            )
        if cfg.wants("LCL"):
            data_vars["lcl_pressure"] = _var(
                np.stack([parcels[k].lcl_pressure for k in parcels], 0),
                ("time", "parcel", "y", "x"),
                "PLCL",
                "lifting condensation level pressure",
                "Pa",
            )
            data_vars["lcl_hgt"] = _var(
                np.stack([parcels[k].lcl_height_agl for k in parcels], 0),
                ("time", "parcel", "y", "x"),
                "HLCL",
                "LCL height above ground",
                "m",
            )
        if cfg.wants("LFC"):
            data_vars["lfc_pressure"] = _var(
                np.stack([parcels[k].lfc_pressure for k in parcels], 0),
                ("time", "parcel", "y", "x"),
                "PLFC",
                "level of free convection pressure",
                "Pa",
            )
            data_vars["lfc_hgt"] = _var(
                np.stack([parcels[k].lfc_height_agl for k in parcels], 0),
                ("time", "parcel", "y", "x"),
                "HLFC",
                "LFC height above ground",
                "m",
            )

    if any(cfg.wants(n) for n in ("HLCY", "USTM", "VSTM")):
        heli = storm_relative_helicity(state)
        if cfg.wants("HLCY"):
            coords["hel_depth"] = ("hel_depth", np.array([3000.0, 1000.0]), {"units": "m", "long_name": "helicity depth AGL"})
            data_vars["hlcy"] = _var(
                np.stack([heli.hel_0_3km, heli.hel_0_1km], 0),
                ("time", "hel_depth", "y", "x"),
                "HLCY",
                "storm relative helicity",
                "m2 s-2",
            )
        if cfg.wants("USTM"):
            data_vars["ustm"] = _var(heli.u_storm, ("time", "y", "x"), "USTM", "u component of storm motion", "m s-1")
        if cfg.wants("VSTM"):
            data_vars["vstm"] = _var(heli.v_storm, ("time", "y", "x"), "VSTM", "v component of storm motion", "m s-1")

    if any(cfg.wants(n) for n in ("TCDC", "LCDC", "MCDC", "HCDC")):
        cld = cloud_fraction(state)
        if cfg.wants("LCDC"):
            data_vars["lcdc"] = _var(cld["low"], ("time", "y", "x"), "LCDC", "low cloud fraction", "%")
        if cfg.wants("MCDC"):
            data_vars["mcdc"] = _var(cld["mid"], ("time", "y", "x"), "MCDC", "mid cloud fraction", "%")
        if cfg.wants("HCDC"):
            data_vars["hcdc"] = _var(cld["high"], ("time", "y", "x"), "HCDC", "high cloud fraction", "%")
        if cfg.wants("TCDC"):
            data_vars["tcdc"] = _var(cld["total"], ("time", "y", "x"), "TCDC", "total cloud fraction", "%")

    data_vars["orog"] = _var(state.terrain, ("time", "y", "x"), "HGT", "surface geopotential height", "m")
    data_vars["psfc"] = _var(state.p_sfc, ("time", "y", "x"), "PRES", "surface pressure", "Pa")
    data_vars["crs"] = _projection(state)

    ds = xr.Dataset(data_vars, coords=coords)
    for name, da in list(ds.data_vars.items()):
        if name == "crs":
            continue
        if "y" in da.dims and "x" in da.dims:
            da.attrs["grid_mapping"] = "crs"
            da.attrs["coordinates"] = "latitude longitude"
    ds.attrs["title"] = "upp-py WRF-ARW post"
    ds.attrs["source"] = f"upp-py post of {state.path}"
    ds.attrs["history"] = f"upp-py {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}"
    ds.attrs["forecast_reference_time"] = state.start_time.strftime("%Y-%m-%dT%H:%M:%SZ")
    ds.attrs["forecast_hour"] = float(state.forecast_hour)
    ds.attrs["mp_physics"] = int(state.mp_physics)
    ds.attrs["Conventions"] = "CF-1.8"
    ds.attrs["DX"] = float(state.dx)
    ds.attrs["DY"] = float(state.dy)
    return ds


def write_netcdf(ds: xr.Dataset, path: str) -> None:
    """Write the dataset. Object coordinates (parcel names) become strings."""
    encoding = {}
    out = ds.copy()
    if "parcel" in out.coords and out["parcel"].dtype == object:
        out = out.assign_coords(parcel=out["parcel"].astype(str))
    for name, da in out.data_vars.items():
        if np.issubdtype(da.dtype, np.floating):
            encoding[name] = {"zlib": True, "complevel": 4, "dtype": "float32"}
    out.to_netcdf(path, encoding=encoding)
