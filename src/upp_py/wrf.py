"""Read a WRF-ARW ``wrfout`` NetCDF the way ``INITPOST.F`` does.

Vertical index ``k = 0`` is the lowest mass level (WRF ``bottom_top``
order). UPP flips this so its ``L = 1`` is the model top; the arithmetic
below is written in WRF order and documented against the flipped index.

Destaggering (``INITPOST.F``):

* ``UH(i) = 0.5 * (U(i) + U(i+1))``
* ``VH(j) = 0.5 * (V(j) + V(j+1))``
* ``WH(k) = 0.5 * (W(k) + W(k+1))``

Hydrostatic pressure is ``P_HYD`` when the file has it (UPP v4.1 no
longer sums ``P+PB``). Temperature is ``(T_wrf + 300) * (p/1e5)**CAPA``.
Specific humidity is ``QVAPOR / (1 + QVAPOR)``. Geopotential height on
w levels is ``(PH + PHB) / g`` with ``g = 9.81``.

Surface interface pressure is dry ``MU + MUB + P_TOP`` plus the smoothed
column vapor pressure (Pyle/Chuang ARW fix in ``INITPOST.F``). Shelter
(2 m) pressure is ``PSHLTR = PSFC * exp(-0.068283 / T_lowest)`` from
``SURFCE.f``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
import xarray as xr

from upp_py.constants import (
    AD05,
    CAPA,
    CFT0,
    D608,
    G,
    NVAPOR_SMOOTH,
    RD,
    SHELTER_COEF,
    T0_THETA,
)
from upp_py.thermo import mixing_ratio_to_specific


@dataclass
class WrfState:
    """Mass-point state for one wrfout time, in WRF bottom-up order."""

    path: str
    valid_time: datetime
    start_time: datetime
    forecast_hour: float
    mp_physics: int
    dx: float
    dy: float
    map_proj: int
    truelat1: float
    truelat2: float
    stand_lon: float
    cen_lat: float
    cen_lon: float
    lat: np.ndarray
    lon: np.ndarray
    mapfac: np.ndarray
    terrain: np.ndarray
    landmask: np.ndarray | None
    p_top: float
    p_mid: np.ndarray  # (nz, ny, nx) Pa, P_HYD
    p_sfc: np.ndarray  # (ny, nx) moist surface pressure
    p_int: np.ndarray  # (nz+1, ny, nx) interface pressure, k=0 at ground
    z_int: np.ndarray  # (nz+1, ny, nx) m MSL
    z_mid: np.ndarray  # (nz, ny, nx) m MSL
    t: np.ndarray  # K
    q: np.ndarray  # kg/kg specific
    u: np.ndarray  # grid-relative m/s, mass point
    v: np.ndarray
    w: np.ndarray  # m/s mass point
    omega: np.ndarray  # Pa/s
    qc: np.ndarray
    qi: np.ndarray
    qr: np.ndarray
    qs: np.ndarray
    qg: np.ndarray
    th2: np.ndarray
    q2_mixing: np.ndarray
    t2: np.ndarray
    q2: np.ndarray
    p_shelter: np.ndarray
    u10: np.ndarray
    v10: np.ndarray
    rainc: np.ndarray  # mm
    rainnc: np.ndarray
    raincv: np.ndarray  # mm per model step (WRF RAINCV)
    rainncv: np.ndarray
    dt: float
    # Model diagnostic. NCAR + Thompson posts the column max of this
    # array as REFC (MDLFLD.f); None if the file does not have it.
    refl_10cm: np.ndarray | None = None
    attrs: dict = field(default_factory=dict)

    @property
    def nz(self) -> int:
        return int(self.t.shape[0])

    @property
    def ny(self) -> int:
        return int(self.t.shape[1])

    @property
    def nx(self) -> int:
        return int(self.t.shape[2])


def _scalar(ds: xr.Dataset, name: str, default=None):
    if name in ds.attrs:
        val = ds.attrs[name]
        if isinstance(val, bytes):
            val = val.decode()
        return val
    if name in ds.variables:
        val = ds[name].values
        if np.size(val) == 1:
            return val.reshape(-1)[0].item() if hasattr(val.reshape(-1)[0], "item") else val.reshape(-1)[0]
        return val
    if default is not None:
        return default
    raise KeyError(name)


def _parse_wrf_time(text: str) -> datetime:
    text = str(text).replace("T", "_")
    for fmt in ("%Y-%m-%d_%H:%M:%S", "%Y-%m-%d_%H:%M:%S.%f"):
        try:
            return datetime.strptime(text[:19] if fmt.endswith("%S") else text, fmt)
        except ValueError:
            continue
    return datetime.strptime(text[:19], "%Y-%m-%d_%H:%M:%S")


def _times_at(ds: xr.Dataset, index: int) -> datetime:
    raw = ds["Times"].values[index]
    if isinstance(raw, bytes):
        text = raw.decode()
    elif isinstance(raw, np.ndarray):
        text = "".join(chr(c) if isinstance(c, (int, np.integer)) else str(c) for c in raw)
        if not text or text.startswith("b"):
            text = bytes(raw).decode().strip("\x00")
    else:
        text = str(raw)
    text = text.strip().replace("\x00", "")
    return _parse_wrf_time(text)


def _get3d(ds: xr.Dataset, name: str, itime: int) -> np.ndarray | None:
    if name not in ds.variables:
        return None
    da = ds[name]
    if "Time" in da.dims:
        da = da.isel(Time=itime)
    data = np.asarray(da.values, dtype=np.float64)
    return data


def _optional3d(ds, name, itime, shape, fill=0.0) -> np.ndarray:
    data = _get3d(ds, name, itime)
    if data is None:
        return np.full(shape, fill, dtype=np.float64)
    return data


def _smooth_vapor(pvapor: np.ndarray) -> np.ndarray:
    """405-pass 8-point smoother from ``INITPOST.F`` (Jacobi, old field)."""
    orig = pvapor.copy()
    field = pvapor.copy()
    ny, nx = field.shape
    for _ in range(NVAPOR_SMOOTH):
        c = field
        new = AD05 * (
            4.0 * (c[1:-1, 0:-2] + c[1:-1, 2:] + c[0:-2, 1:-1] + c[2:, 1:-1])
            + c[0:-2, 0:-2]
            + c[0:-2, 2:]
            + c[2:, 0:-2]
            + c[2:, 2:]
        ) - CFT0 * c[1:-1, 1:-1]
        field = field.copy()
        field[1:-1, 1:-1] = new
    # Boundaries: original plus the interior increment of the adjacent row.
    out = field
    out[0, 1:-1] = orig[0, 1:-1] + (field[1, 1:-1] - orig[1, 1:-1])
    out[-1, 1:-1] = orig[-1, 1:-1] + (field[-2, 1:-1] - orig[-2, 1:-1])
    out[:, 0] = orig[:, 0] + (field[:, 1] - orig[:, 1])
    out[:, -1] = orig[:, -1] + (field[:, -2] - orig[:, -2])
    return out


def _enforce_monotonic_pressure(p_mid: np.ndarray) -> np.ndarray:
    """Match the non-monotonic ``PMID`` repairs in ``INITPOST.F``.

    UPP index increases downward. In WRF order, pressure must strictly
    decrease with ``k``. Interior bad levels are replaced by the mean of
    their neighbors; the lowest level uses the ``6/11`` extrapolation.
    """
    p = p_mid.copy()
    nz = p.shape[0]
    # Interior, scanning from the second-lowest upward (UPP ll = lm-1 ... 2).
    for k in range(1, nz - 1):
        # UPP ll corresponds to WRF k, ll-1 is above (k+1), ll+1 is below (k-1).
        bad = (p[k + 1] - p[k]) >= 0.0
        if np.any(bad):
            p[k] = np.where(bad, 0.5 * (p[k - 1] + p[k + 1]), p[k])
    if nz >= 3:
        fact = 6.0 / 11.0
        bad = (p[1] - p[0]) >= 0.0
        if np.any(bad):
            # PMID(LL)=PMID(ll-1) + fact*(PMID(ll-1)-PMID(ll-2))
            # ll = bottom = k 0, ll-1 = k 1, ll-2 = k 2
            repl = p[1] + fact * (p[1] - p[2])
            p[0] = np.where(bad, repl, p[0])
    return p


def open_wrfout(path: str, time_index: int = 0) -> WrfState:
    """Load one time from a WRF-ARW history file.

    Parameters
    ----------
    path
        Path to ``wrfout_d0*``.
    time_index
        Which ``Time`` record to post. UPP posts the time named in the
        ``itag`` file; the default is the first record.
    """
    ds = xr.open_dataset(path, decode_times=False, mask_and_scale=False)
    try:
        return _load(ds, path, time_index)
    finally:
        ds.close()


def _load(ds: xr.Dataset, path: str, itime: int) -> WrfState:
    t_pert = _get3d(ds, "T", itime)
    if t_pert is None:
        raise ValueError("wrfout is missing T (potential-temperature perturbation)")
    nz, ny, nx = t_pert.shape

    p_hyd = _get3d(ds, "P_HYD", itime)
    if p_hyd is None:
        pb = _get3d(ds, "PB", itime)
        pp = _get3d(ds, "P", itime)
        if pb is None or pp is None:
            raise ValueError("wrfout needs P_HYD or both P and PB")
        p_hyd = pb + pp
    theta = t_pert + T0_THETA
    # Convert theta -> T with the pressure *before* the monotonic repair,
    # which is the order in INITPOST.F.
    temperature = theta * (p_hyd * 1.0e-5) ** CAPA
    p_mid = _enforce_monotonic_pressure(p_hyd)

    qv = _optional3d(ds, "QVAPOR", itime, (nz, ny, nx))
    q = mixing_ratio_to_specific(qv)

    u_stag = _get3d(ds, "U", itime)
    v_stag = _get3d(ds, "V", itime)
    w_stag = _get3d(ds, "W", itime)
    if u_stag is None or v_stag is None or w_stag is None:
        raise ValueError("wrfout is missing U, V, or W")
    # After the Time axis is removed: U is (nz, ny, nx+1), V is (nz, ny+1, nx).
    u = 0.5 * (u_stag[..., :-1] + u_stag[..., 1:])
    v = 0.5 * (v_stag[:, :-1, :] + v_stag[:, 1:, :])
    w = 0.5 * (w_stag[:-1] + w_stag[1:])

    phb = _get3d(ds, "PHB", itime)
    ph = _get3d(ds, "PH", itime)
    if phb is None or ph is None:
        raise ValueError("wrfout is missing PH/PHB")
    z_int = (ph + phb) / G
    z_mid = 0.5 * (z_int[:-1] + z_int[1:])

    mu = np.asarray(ds["MU"].isel(Time=itime).values, dtype=np.float64)
    mub = np.asarray(ds["MUB"].isel(Time=itime).values, dtype=np.float64)
    p_top = float(np.asarray(ds["P_TOP"].values).reshape(-1)[0])
    p_dry = mu + mub + p_top

    # Column integrated vapor pressure (Pa). UPP loops from the top down;
    # QMEAN mixes a level with the one below it, except the lowest level.
    pvapor = np.zeros((ny, nx), dtype=np.float64)
    for k in range(nz - 1, -1, -1):
        dz = z_int[k + 1] - z_int[k]
        rho = p_mid[k] / (RD * temperature[k])
        if k > 0:
            qmean = 0.5 * (q[k] + q[k - 1])
        else:
            qmean = q[k]
        pvapor = pvapor + G * rho * dz * qmean
    pvapor = _smooth_vapor(pvapor)
    p_sfc = p_dry + pvapor
    # Keep the surface interface above the lowest mid-level pressure.
    p_sfc = np.maximum(p_sfc, p_mid[0] * 1.001)

    p_int = np.empty((nz + 1, ny, nx), dtype=np.float64)
    p_int[0] = p_sfc
    p_int[-1] = p_top
    if nz > 1:
        p_int[1:-1] = 0.5 * (p_mid[:-1] + p_mid[1:])

    omega = -w * p_mid * G / (RD * temperature * (1.0 + D608 * q))

    qc = _optional3d(ds, "QCLOUD", itime, (nz, ny, nx))
    qi = _optional3d(ds, "QICE", itime, (nz, ny, nx))
    qr = _optional3d(ds, "QRAIN", itime, (nz, ny, nx))
    qs = _optional3d(ds, "QSNOW", itime, (nz, ny, nx))
    qg = _optional3d(ds, "QGRAUP", itime, (nz, ny, nx))
    # WSM3 partitions cloud water/ice by temperature (INITPOST.F).
    mp = int(_scalar(ds, "MP_PHYSICS", 0))
    if mp == 3 and "QICE" not in ds.variables:
        cold = temperature < 273.15
        qi = np.where(cold, qc, 0.0)
        qc = np.where(cold, 0.0, qc)

    lat = np.asarray(ds["XLAT"].isel(Time=itime).values, dtype=np.float64)
    lon = np.asarray(ds["XLONG"].isel(Time=itime).values, dtype=np.float64)
    if "MAPFAC_M" in ds.variables:
        mapfac = np.asarray(ds["MAPFAC_M"].isel(Time=itime).values, dtype=np.float64)
    else:
        mapfac = np.ones((ny, nx), dtype=np.float64)
    terrain = np.asarray(ds["HGT"].isel(Time=itime).values, dtype=np.float64)
    land = None
    if "LANDMASK" in ds.variables:
        land = np.asarray(ds["LANDMASK"].isel(Time=itime).values, dtype=np.float64)

    th2 = np.asarray(ds["TH2"].isel(Time=itime).values, dtype=np.float64) if "TH2" in ds else theta[0]
    q2m = np.asarray(ds["Q2"].isel(Time=itime).values, dtype=np.float64) if "Q2" in ds else qv[0]
    u10 = np.asarray(ds["U10"].isel(Time=itime).values, dtype=np.float64) if "U10" in ds else u[0]
    v10 = np.asarray(ds["V10"].isel(Time=itime).values, dtype=np.float64) if "V10" in ds else v[0]
    rainc = np.asarray(ds["RAINC"].isel(Time=itime).values, dtype=np.float64) if "RAINC" in ds else np.zeros((ny, nx))
    rainnc = np.asarray(ds["RAINNC"].isel(Time=itime).values, dtype=np.float64) if "RAINNC" in ds else np.zeros((ny, nx))
    raincv = np.asarray(ds["RAINCV"].isel(Time=itime).values, dtype=np.float64) if "RAINCV" in ds else np.zeros((ny, nx))
    rainncv = np.asarray(ds["RAINNCV"].isel(Time=itime).values, dtype=np.float64) if "RAINNCV" in ds else np.zeros((ny, nx))

    # 2 m temperature: TH2 is potential temperature; SURFCE.f converts
    # with shelter pressure, not PSFC.
    t_low = temperature[0]
    p_shelter = p_sfc * np.exp(-SHELTER_COEF / t_low)
    t2 = th2 * (p_shelter * 1.0e-5) ** CAPA
    q2 = mixing_ratio_to_specific(q2m)

    valid = _times_at(ds, itime)
    start_attr = ds.attrs.get("START_DATE") or ds.attrs.get("SIMULATION_START_DATE")
    if isinstance(start_attr, bytes):
        start_attr = start_attr.decode()
    start = _parse_wrf_time(str(start_attr)) if start_attr else valid
    fhr = (valid - start).total_seconds() / 3600.0

    dx = float(_scalar(ds, "DX"))
    dy = float(_scalar(ds, "DY"))
    dt = abs(float(_scalar(ds, "DT", 0.0)))

    def _fattr(name, default=0.0):
        if name not in ds.attrs:
            return default
        return float(ds.attrs[name])

    return WrfState(
        path=path,
        valid_time=valid,
        start_time=start,
        forecast_hour=fhr,
        mp_physics=mp,
        dx=dx,
        dy=dy,
        map_proj=int(_scalar(ds, "MAP_PROJ", 1)),
        truelat1=_fattr("TRUELAT1"),
        truelat2=_fattr("TRUELAT2"),
        stand_lon=_fattr("STAND_LON"),
        cen_lat=_fattr("CEN_LAT") or _fattr("MOAD_CEN_LAT"),
        cen_lon=_fattr("CEN_LON"),
        lat=lat,
        lon=lon,
        mapfac=mapfac,
        terrain=terrain,
        landmask=land,
        p_top=p_top,
        p_mid=p_mid,
        p_sfc=p_sfc,
        p_int=p_int,
        z_int=z_int,
        z_mid=z_mid,
        t=temperature,
        q=q,
        u=u,
        v=v,
        w=w,
        omega=omega,
        qc=np.maximum(qc, 0.0),
        qi=np.maximum(qi, 0.0),
        qr=np.maximum(qr, 0.0),
        qs=np.maximum(qs, 0.0),
        qg=np.maximum(qg, 0.0),
        th2=th2,
        q2_mixing=q2m,
        t2=t2,
        q2=q2,
        p_shelter=p_shelter,
        u10=u10,
        v10=v10,
        rainc=rainc,
        rainnc=rainnc,
        raincv=raincv,
        rainncv=rainncv,
        dt=dt,
        refl_10cm=_get3d(ds, "REFL_10CM", itime),
        attrs={k: (v.decode() if isinstance(v, bytes) else v) for k, v in ds.attrs.items() if isinstance(v, (str, bytes, int, float, np.number))},
    )
