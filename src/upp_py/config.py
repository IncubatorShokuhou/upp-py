"""Field selection.

The YAML file is a short stand-in for ``postxconfig-NT-WRF.txt`` /
``wrf_cntrl.parm``: a list of GRIB2 short names and the isobaric levels
to post. It does not parse the Fortran control-file syntax.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from upp_py.constants import SPL_PA

# Short names this package knows how to compute. The set is the WRF-ARW
# core requested for upp-py, not the full UPP catalog.
KNOWN_FIELDS = (
    "PRMSL",
    "MSLET",
    "TMP",
    "DPT",
    "RH",
    "SPFH",
    "HGT",
    "UGRD",
    "VGRD",
    "VVEL",
    "ABSV",
    "TMP2M",
    "DPT2M",
    "RH2M",
    "SPFH2M",
    "UGRD10M",
    "VGRD10M",
    "REFC",
    "REFD",
    "CAPE",
    "CIN",
    "LCL",
    "LFC",
    "HLCY",
    "USTM",
    "VSTM",
    "PWAT",
    "APCP",
    "ACPCP",
    "NCPCP",
    "APCP_BUCKET",
    "ACPCP_BUCKET",
    "NCPCP_BUCKET",
    "SNOW_BUCKET",
    "GRAUP_BUCKET",
    "PRATE",
    "TCDC",
    "LCDC",
    "MCDC",
    "HCDC",
    "VIS",
    "VIS_GSD",
    "TIPD",
    "BRTMP",
    "QNWFA",
    "QNIFA",
    "QNCLOUD",
    "SMOKE",
    "DUST",
    "AOTK",
    "AERASY",
    "AERSSA",
)

DEFAULT_FIELDS = (
    "PRMSL",
    "MSLET",
    "HGT",
    "TMP",
    "DPT",
    "RH",
    "SPFH",
    "UGRD",
    "VGRD",
    "VVEL",
    "ABSV",
    "TMP2M",
    "DPT2M",
    "RH2M",
    "SPFH2M",
    "UGRD10M",
    "VGRD10M",
    "REFC",
    "REFD",
    "CAPE",
    "CIN",
    "LCL",
    "LFC",
    "HLCY",
    "USTM",
    "VSTM",
    "PWAT",
    "APCP",
    "ACPCP",
    "NCPCP",
    "APCP_BUCKET",
    "ACPCP_BUCKET",
    "NCPCP_BUCKET",
    "PRATE",
    "TCDC",
    "LCDC",
    "MCDC",
    "HCDC",
    "VIS",
    "TIPD",
)


@dataclass
class PostConfig:
    """What to write.

    Parameters
    ----------
    fields
        GRIB2 short names. Surface-based, mixed-layer and most-unstable
        CAPE/CIN are all written when ``CAPE`` or ``CIN`` is requested;
        they are distinguished by the ``parcel`` coordinate.
    pressure_pa
        Isobaric levels, Pa, any order.
    reflectivity_agl_m
        Heights for ``REFD``.
    accum_hours
        If set, also write ``apcp_accum``, the trailing accumulation of
        that many hours. ``None`` leaves ``apcp`` as the run total only.
        The window is applied by :func:`upp_py.post_wrfout` when more
        than one time is read, and also for a single time (clamped to
        the forecast start).
    """

    fields: tuple[str, ...] = DEFAULT_FIELDS
    pressure_pa: tuple[float, ...] = SPL_PA
    reflectivity_agl_m: tuple[float, ...] = (1000.0, 4000.0)
    accum_hours: float | None = None

    def wants(self, name: str) -> bool:
        return name in self.fields


def load_config(path: str | Path | None) -> PostConfig:
    """Load YAML. ``None`` returns :class:`PostConfig` defaults."""
    if path is None:
        return PostConfig()
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("reading a YAML config requires PyYAML") from exc
    with open(path, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}
    fields = tuple(raw.get("fields", DEFAULT_FIELDS))
    unknown = [f for f in fields if f not in KNOWN_FIELDS]
    if unknown:
        raise ValueError(f"unknown fields {unknown}; known: {', '.join(KNOWN_FIELDS)}")
    if "pressure_hpa" in raw:
        pressure = tuple(float(p) * 100.0 for p in raw["pressure_hpa"])
    elif "pressure_pa" in raw:
        pressure = tuple(float(p) for p in raw["pressure_pa"])
    else:
        pressure = SPL_PA
    agl = tuple(float(z) for z in raw.get("reflectivity_agl_m", (1000.0, 4000.0)))
    accum = raw.get("accum_hours")
    accum_hours = None if accum is None else float(accum)
    return PostConfig(fields=fields, pressure_pa=pressure, reflectivity_agl_m=agl, accum_hours=accum_hours)


def dump_example(path: str | Path) -> None:
    """Write a commented example config."""
    text = """# upp-py field list. Short names follow GRIB2.
# pressure_hpa is only used for isobaric fields (HGT TMP DPT RH SPFH UGRD VGRD VVEL ABSV).
# MSLP always integrates the full standard-level set internally.
fields:
  - PRMSL
  - HGT
  - TMP
  - DPT
  - RH
  - SPFH
  - UGRD
  - VGRD
  - VVEL
  - ABSV
  - TMP2M
  - DPT2M
  - RH2M
  - UGRD10M
  - VGRD10M
  - REFC
  - REFD
  - CAPE
  - CIN
  - LCL
  - LFC
  - HLCY
  - USTM
  - VSTM
  - PWAT
  - APCP
  - ACPCP
  - NCPCP
  - APCP_BUCKET
  - PRATE
  - VIS
  - TIPD
  - TCDC
  - LCDC
  - MCDC
  - HCDC
pressure_hpa: [250, 500, 700, 850, 925, 1000]
reflectivity_agl_m: [1000, 4000]
# accum_hours: 1
"""
    Path(path).write_text(text, encoding="utf-8")
