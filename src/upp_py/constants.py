"""Physical constants from UPP v4.1 ``params.F`` / ``params_mod``.

Values are copied from EMC_post tag ``dtc_post_v4.1.0``
(``sorc/ncep_post.fd/params.F``) so thermodynamic formulas stay on the
same constants the Fortran post uses. Do not replace them with CODATA
or ``scipy.constants`` values.
"""

from __future__ import annotations

# params.F
G = 9.81
RD = 287.04
ROG = RD / G
GI = 1.0 / G
CP = 1004.6
CAPA = 0.28589641  # used by CALCAPE / theta conversion, not exactly RD/CP
R_OVER_CP = RD / CP  # TABLE.f saturation lookup uses this exponent
P1000 = 1000.0e2
D608 = 0.608
GAMMA = 6.5e-3  # K/m, standard-atmosphere lapse used below ground
RGAMOG = RD * GAMMA / G
TFRZ = 273.15
EPS_MASS = 18.015 / 28.964  # FMW/FMD
ONEPS = 1.0 - EPS_MASS
ELWV = 2.50e6
ELIVW = 2.72e6  # params.F, ELOCP numerator
ELOCP = ELIVW / CP
ELIWV_TABLE = 2.683e6  # TABLE.f moist-adiabat constant (distinct from ELIVW)

# CALRH saturation specific humidity: QC = PQ0/p * exp(A2*(T-A3)/(T-A4))
PQ0 = 379.90516
A2 = 17.2693882
A3 = 273.16
A4 = 35.86
RHMIN = 1.0e-6

# DEWPOINT.f
H1M12 = 1.0e-12
D001 = 0.001  # Pa -> centibar

# MEMSLP successive over-relaxation (params.F)
OVERRC = 1.10
AD05 = OVERRC * 0.05
CFT0 = OVERRC - 1.0
NRLX_SLP = 100
NFILL_SLP = 0
SLOPE_SLP = -6.6e-4  # K/m^2 s^-2, i.e. dT/d(gz)

# shelter pressure reduction, SURFCE.f (2 m hydrostatic factor)
SHELTER_COEF = 0.068283

# Earth / Coriolis. ERAD is params.F; omega is the standard value UPP's
# coriolis parameter is built from (2*omega*sin(phi)).
ERAD = 6.371e6
OMEGA_EARTH = 7.2921e-5

# Standard isobaric levels (Pa), CTLBLK.f SPLDEF. Ordered top -> bottom.
SPL_PA = (
    200.0,
    500.0,
    700.0,
    1000.0,
    2000.0,
    3000.0,
    5000.0,
    7000.0,
    7500.0,
    10000.0,
    12500.0,
    15000.0,
    17500.0,
    20000.0,
    22500.0,
    25000.0,
    27500.0,
    30000.0,
    32500.0,
    35000.0,
    37500.0,
    40000.0,
    42500.0,
    45000.0,
    47500.0,
    50000.0,
    52500.0,
    55000.0,
    57500.0,
    60000.0,
    62500.0,
    65000.0,
    67500.0,
    70000.0,
    72500.0,
    75000.0,
    77500.0,
    80000.0,
    82500.0,
    85000.0,
    87500.0,
    90000.0,
    92500.0,
    95000.0,
    97500.0,
    100000.0,
)

# BNDLYR.f / CTLBLK.f
NBND = 6
DPBND_LAYER = 30.0e2  # Pa, each boundary-layer slab
# CLDRAD.f surface-based CAPE search depth (comment says 70 mb; the code uses 10 hPa)
DPBND_SURFACE = 10.0e2
# MISCLN.f most-unstable search (GRIB level is labelled 255 hPa)
DPBND_MOST_UNSTABLE = 300.0e2

# Cloud-fraction pressure tops (Pa), INITPOST.F
PTOP_LOW = 64200.0
PTOP_MID = 35000.0
PTOP_HIGH = 15000.0
CLOUD_DEF_P = 1.0e-7  # kg/kg qc+qi threshold
CLOUD_RADIUS_M = 30000.0

# Reflectivity floor, params.F
DBZMIN = -20.0

# Thompson / Stoelinga constants from PARAMR.f + MDLFLD.f (mp_physics=8)
PI = 3.141592653589793
RON = 8.0e6
RON2 = 1.0e9
RON_MIN = 2.0e7
SON = 2.0e7
GON = 4.0e6
DRAIN = 1000.0
DSNOW = 100.0
DGRAUPEL = 400.0
QR0 = 0.0002
DELQR0 = 0.5 * QR0
CONST1R = (RON2 - RON_MIN) * 0.5
CONST2R = (RON2 + RON_MIN) * 0.5
TOPR = PI * DRAIN * RON
TOPS = PI * DSNOW * SON
TOPG = PI * DGRAUPEL * GON
CONSTD = 1.0 / 0.52
CONSTGB = 12.0 / 13.0
CONST_NG1 = (1.57**CONSTD) * ((PI * DGRAUPEL) ** CONSTGB)
CONST_NG2 = -CONSTGB
REFL_ALPHA = 0.224
REFL_R1 = 1.0e-15
ZR_A = 300.0
ZR_B = 1.4

# Pyle/Chuang column-vapor smoother iteration count, INITPOST.F
NVAPOR_SMOOTH = 405

# WRF potential-temperature base used by INITPOST.F
T0_THETA = 300.0
