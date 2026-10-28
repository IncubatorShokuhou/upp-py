"""Assemble an xarray dataset from a wrfout file.

This is the driver UPP splits across ``INITPOST``, ``MDLFLD``, ``MDL2P``,
``SURFCE``, ``CLDRAD`` and ``MISCLN``. Only the WRF-ARW core products
are filled; each variable's ``grib2_shortName`` attribute is the name
wgrib2 prints for the same quantity.
"""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone

import numpy as np
import xarray as xr

from upp_py.aerosol import clip_number_concentration, mass_concentration
from upp_py.brightness import VIS_GSD_SKIP_REASON, compute_brtmp, skip_reason
from upp_py.cape import cape_best, cape_mixed_layer, cape_most_unstable, cape_surface
from upp_py.clouds import cloud_fraction
from upp_py.config import PostConfig
from upp_py.constants import GI, SPL_PA
from upp_py.helicity import storm_relative_helicity
from upp_py.icing import in_flight_icing
from upp_py.interp import interpolate_pressure_levels
from upp_py.mslp import membrane_slp, shuell_slp
from upp_py.precip import (
    bucket_convective_mm,
    bucket_graupel_mm,
    bucket_gridscale_mm,
    bucket_snow_mm,
    bucket_total_mm,
    continuous_convective_mm,
    continuous_gridscale_mm,
    continuous_total_mm,
    precipitation_rate,
    window_accumulation_mm,
)
from upp_py.reflectivity import reflectivity
from upp_py.thermo import dewpoint, relative_humidity
from upp_py.visibility import surface_visibility
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


def _ensure_model_level(coords: dict, nz: int) -> None:
    if "model_level" in coords:
        return
    coords["model_level"] = (
        "model_level",
        np.arange(nz, dtype=np.int32),
        {"long_name": "WRF mass level, 0 nearest the ground", "positive": "up"},
    )


def _add_aerosols(state: WrfState, cfg: PostConfig, data_vars: dict, coords: dict, notes: list[str]) -> None:
    """Model-level aerosols. Missing inputs are recorded in ``notes``."""
    number = (
        ("QNWFA", "qnwfa", state.qnwfa),
        ("QNIFA", "qnifa", state.qnifa),
        ("QNCLOUD", "qncloud", state.qncloud),
    )
    for short, key, arr in number:
        if not cfg.wants(short):
            continue
        if int(state.mp_physics) != 28:
            notes.append(
                f"{short}: INITPOST.F reads this only when imp_physics==28 "
                f"(this file is mp_physics={int(state.mp_physics)})."
            )
            continue
        if arr is None:
            notes.append(f"{short}: wrfout has no {short}.")
            continue
        _ensure_model_level(coords, state.nz)
        data_vars[key] = _var(
            clip_number_concentration(arr),
            ("time", "model_level", "y", "x"),
            short,
            f"{short} number concentration",
            "kg-1",
        )
    if cfg.wants("SMOKE"):
        smoke = state.aerosols.get("smoke")
        if smoke is None:
            notes.append("SMOKE: wrfout has no smoke (INITPOST.F VarName smoke).")
        else:
            _ensure_model_level(coords, state.nz)
            data_vars["smoke"] = _var(
                mass_concentration(smoke, state.p_mid, state.t),
                ("time", "model_level", "y", "x"),
                "SMOKE",
                "smoke mass concentration",
                "ug m-3",
            )
    if cfg.wants("DUST"):
        bins = []
        labels = []
        for ibin in range(1, 6):
            raw = state.aerosols.get(f"DUST_{ibin}")
            if raw is None:
                continue
            bins.append(mass_concentration(raw, state.p_mid, state.t))
            labels.append(ibin)
        if not bins:
            notes.append(
                "DUST: wrfout has no DUST_1..DUST_5. INITPOST.F v4.1 wraps those reads in if(1==2), "
                "so unipost does not fill them either. If the arrays are present, MDLFLD.f converts "
                "ug kg-1 to ug m-3 with density."
            )
        else:
            _ensure_model_level(coords, state.nz)
            coords["dust_bin"] = ("dust_bin", np.asarray(labels, dtype=np.int32), {"long_name": "GOCART dust bin"})
            data_vars["dust"] = _var(
                np.stack(bins, axis=0),
                ("time", "dust_bin", "model_level", "y", "x"),
                "DUST",
                "dust mass concentration",
                "ug m-3",
            )
    column = (
        ("AOTK", "aotk", "TAOD5502D", "aerosol optical thickness at 550 nm", "1"),
        ("AERASY", "aerasy", "AERASY2D", "aerosol asymmetry parameter", "1"),
        ("AERSSA", "aerssa", "AERSSA2D", "aerosol single-scattering albedo", "1"),
    )
    for short, key, wrf_name, long_name, units in column:
        if not cfg.wants(short):
            continue
        arr = state.aerosols.get(wrf_name)
        if arr is None:
            notes.append(f"{short}: wrfout has no {wrf_name}.")
            continue
        data_vars[key] = _var(arr, ("time", "y", "x"), short, long_name, units)


def _add_brtmp(state: WrfState, cfg: PostConfig, data_vars: dict, coords: dict, notes: list[str]) -> None:
    """Write ``brtmp`` when ``upp_crtm`` is installed. Otherwise record a skip."""
    try:
        tb = compute_brtmp(state, satellite=cfg.brtmp_satellite, channels=cfg.brtmp_channels)
    except Exception as exc:
        # BrtmpUnavailable, and any unexpected failure from the optional package.
        notes.append(f"BRTMP: {exc}")
        return
    coords["abi_channel"] = (
        "abi_channel",
        np.asarray(cfg.brtmp_channels, dtype=np.int32),
        {"long_name": "ABI channel number", "satellite": cfg.brtmp_satellite},
    )
    data_vars["brtmp"] = _var(
        tb,
        ("time", "abi_channel", "y", "x"),
        "BRTMP",
        f"brightness temperature {cfg.brtmp_satellite}",
        "K",
        satellite=cfg.brtmp_satellite,
    )


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


def post_wrfout(
    path: str,
    config: PostConfig | None = None,
    time_index: int = 0,
    time_indices: tuple[int, ...] | list[int] | None = None,
    accum_hours: float | None = None,
) -> xr.Dataset:
    """Post one or more wrfout times to an :class:`xarray.Dataset`.

    Parameters
    ----------
    path
        WRF-ARW history file.
    config
        Field list. ``None`` posts the default core set on the standard
        UPP pressure levels.
    time_index
        Which time record to read when ``time_indices`` is omitted.
    time_indices
        Several records from the same file, in any order. The result is
        sorted by valid time. ``apcp`` at each time is the continuous
        run total. Set ``accum_hours`` to also write ``apcp_accum``.
    accum_hours
        Overrides :attr:`PostConfig.accum_hours` when not ``None``.
    """
    cfg = config or PostConfig()
    if accum_hours is not None:
        cfg = replace(cfg, accum_hours=accum_hours)
    indices = (int(time_index),) if time_indices is None else tuple(int(i) for i in time_indices)
    if not indices:
        raise ValueError("time_indices is empty")
    states = [open_wrfout(path, time_index=i) for i in indices]
    order = sorted(range(len(states)), key=lambda i: states[i].valid_time)
    states = [states[i] for i in order]
    parts = [dataset_from_state(state, cfg) for state in states]
    if len(parts) == 1 and cfg.accum_hours is None:
        return parts[0]
    ds = _concat_times(parts)
    if cfg.accum_hours is not None and (cfg.wants("APCP") or cfg.wants("APCP_BUCKET")):
        ds = _add_window(ds, states, float(cfg.accum_hours))
    return ds


def _concat_times(parts: list[xr.Dataset]) -> xr.Dataset:
    crs = parts[0]["crs"]
    stripped = [part.drop_vars("crs") for part in parts]
    ds = xr.concat(stripped, dim="time", combine_attrs="override")
    ds["crs"] = crs
    return ds


def _add_window(ds: xr.Dataset, states: list[WrfState], hours: float) -> xr.Dataset:
    amounts, actual = window_accumulation_mm(states, hours)
    stacked = np.stack(amounts, axis=0)
    ds = ds.assign_coords(accum_hours_actual=("time", np.asarray(actual, dtype=np.float64)))
    ds["apcp_accum"] = _var(
        stacked,
        ("time", "y", "x"),
        "APCP",
        "precipitation accumulated over the requested window",
        "kg m-2",
        grib_stat="accumulation",
        accum_length_hours=float(hours),
    )
    # _var adds a time axis when the array is missing it. stacked already
    # has time as axis 0, and _var sees ndim == len(dims), so it keeps it.
    return ds


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
        "forecast_hour": ("time", np.array([state.forecast_hour], dtype=np.float64), {"units": "hours"}),
    }

    # Isobaric block. MSLP needs the full standard set even if the user
    # only asked for a few output levels, because MEMSLP walks down it.
    notes: list[str] = []
    need_iso = any(
        cfg.wants(n)
        for n in ("HGT", "TMP", "DPT", "RH", "SPFH", "UGRD", "VGRD", "VVEL", "ABSV", "PRMSL", "MSLET", "TIPD")
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
        if cfg.wants("TIPD"):
            tipd = []
            for lv in levels_out:
                rh_pct = relative_humidity(lv.pressure, lv.temperature, lv.specific_humidity, percent=True)
                tipd.append(in_flight_icing(lv.temperature, rh_pct, lv.omega))
            data_vars["tipd"] = _var(
                np.stack(tipd, 0),
                ("time", "isobaric", "y", "x"),
                "TIPD",
                "in-flight icing category",
                "1",
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
        # Continuous total. See upp_py.precip.continuous_total_mm.
        data_vars["apcp"] = _var(
            continuous_total_mm(state),
            ("time", "y", "x"),
            "APCP",
            "total precipitation since forecast start",
            "kg m-2",
            grib_stat="accumulation",
        )
    if cfg.wants("ACPCP"):
        data_vars["acpcp"] = _var(
            continuous_convective_mm(state),
            ("time", "y", "x"),
            "ACPCP",
            "convective precipitation since forecast start",
            "kg m-2",
            grib_stat="accumulation",
        )
    if cfg.wants("NCPCP"):
        data_vars["ncpcp"] = _var(
            continuous_gridscale_mm(state),
            ("time", "y", "x"),
            "NCPCP",
            "grid-scale precipitation since forecast start",
            "kg m-2",
            grib_stat="accumulation",
        )
    bucket_length_h = (state.prec_acc_dt / 3600.0) if state.prec_acc_dt else None
    if cfg.wants("APCP_BUCKET") or cfg.wants("ACPCP_BUCKET") or cfg.wants("NCPCP_BUCKET"):
        missing_prec = state.prec_acc_c is None and state.prec_acc_nc is None
        if missing_prec:
            notes.append(
                "APCP_BUCKET/ACPCP_BUCKET/NCPCP_BUCKET: wrfout has no PREC_ACC_C or PREC_ACC_NC "
                "(INITPOST.F VarName PREC_ACC_C / PREC_ACC_NC)."
            )
        else:
            if cfg.wants("APCP_BUCKET"):
                data_vars["apcp_bucket"] = _var(
                    bucket_total_mm(state),
                    ("time", "y", "x"),
                    "APCP",
                    "bucket total precipitation over PREC_ACC_DT",
                    "kg m-2",
                    grib_stat="accumulation",
                    accum_length_hours=bucket_length_h if bucket_length_h is not None else state.forecast_hour,
                    prec_acc_dt_s=float(state.prec_acc_dt),
                )
            if cfg.wants("ACPCP_BUCKET"):
                conv_b = bucket_convective_mm(state)
                if conv_b is None:
                    notes.append("ACPCP_BUCKET: wrfout has no PREC_ACC_C.")
                else:
                    data_vars["acpcp_bucket"] = _var(
                        conv_b,
                        ("time", "y", "x"),
                        "ACPCP",
                        "bucket convective precipitation over PREC_ACC_DT",
                        "kg m-2",
                        grib_stat="accumulation",
                        accum_length_hours=bucket_length_h if bucket_length_h is not None else state.forecast_hour,
                    )
            if cfg.wants("NCPCP_BUCKET"):
                grid_b = bucket_gridscale_mm(state)
                if grid_b is None:
                    notes.append("NCPCP_BUCKET: wrfout has no PREC_ACC_NC.")
                else:
                    data_vars["ncpcp_bucket"] = _var(
                        grid_b,
                        ("time", "y", "x"),
                        "NCPCP",
                        "bucket grid-scale precipitation over PREC_ACC_DT",
                        "kg m-2",
                        grib_stat="accumulation",
                        accum_length_hours=bucket_length_h if bucket_length_h is not None else state.forecast_hour,
                    )
    if cfg.wants("SNOW_BUCKET"):
        snow_b = bucket_snow_mm(state)
        if snow_b is None:
            notes.append("SNOW_BUCKET: wrfout has no SNOW_ACC_NC.")
        else:
            data_vars["snow_bucket"] = _var(
                snow_b,
                ("time", "y", "x"),
                "WEASD",
                "SNOW_ACC_NC bucket (UPP field 437 pname WEASD)",
                "kg m-2",
                grib_stat="accumulation",
                accum_length_hours=bucket_length_h if bucket_length_h is not None else state.forecast_hour,
            )
    if cfg.wants("GRAUP_BUCKET"):
        gr_b = bucket_graupel_mm(state)
        if gr_b is None:
            notes.append("GRAUP_BUCKET: wrfout has no GRAUP_ACC_NC.")
        else:
            data_vars["graup_bucket"] = _var(
                gr_b,
                ("time", "y", "x"),
                "GRAUP",
                "GRAUP_ACC_NC bucket (SURFCE.f field 775)",
                "kg m-2",
                grib_stat="accumulation",
            )
    if cfg.wants("PRATE"):
        data_vars["prate"] = _var(
            precipitation_rate(state),
            ("time", "y", "x"),
            "PRATE",
            "instantaneous precipitation rate",
            "kg m-2 s-1",
        )

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

    if cfg.wants("VIS"):
        data_vars["vis"] = _var(
            surface_visibility(state),
            ("time", "y", "x"),
            "VIS",
            "surface horizontal visibility",
            "m",
        )
    if cfg.wants("VIS_GSD"):
        notes.append("VIS_GSD: " + VIS_GSD_SKIP_REASON)
    if cfg.wants("BRTMP"):
        _add_brtmp(state, cfg, data_vars, coords, notes)

    _add_aerosols(state, cfg, data_vars, coords, notes)

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
    if notes:
        ds.attrs["skipped"] = " | ".join(notes)
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
