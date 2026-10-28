"""Satellite brightness temperature.

UPP v4.1 computes simulated brightness temperature in
``CALRAD_WCLOUD_newcrtm.f`` (available-field 275, pname ``BRTMP``).
The routine calls the CRTM forward model for GOES-16 and GOES-17 ABI
infrared channels 7-16 and reads CRTM spectral and transmittance
coefficients from a fix directory. Those coefficients are not in a
wrfout file.

upp-py does not wrap CRTM and does not substitute a clear-sky
Stefan-Boltzmann temperature for that product.
"""

BRTMP_SKIP_REASON = (
    "CALRAD_WCLOUD_newcrtm.f calls CRTM for GOES-16/17 ABI IR channels 7-16. "
    "The SpcCoeff and TauCoeff files are not in wrfout, and upp-py does not wrap CRTM."
)

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
