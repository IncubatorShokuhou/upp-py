"""Surface visibility from ``CALVIS.f`` (Stoelinga and Warner).

``MDLFLD.f`` calls ``CALVIS`` for available-field 180 on every
microphysics option other than the GSD branch. Inputs are the lowest
mass level. Hydrometeor mixing ratios are passed through. Water vapor
is the specific humidity ``Q`` produced by ``INITPOST.F`` (``QVAPOR``
divided by ``1 + QVAPOR``), even though the ``CALVIS`` comment says
mixing ratio. The call site passes that ``Q``, and so does this function.

Before the call, ``MDLFLD.f`` adds a convective contribution when
``CPRATE > 0`` and ``imp_physics /= 99``. ``CPRATE`` is ``RAINCV`` in
metres per step. The Zhao option 99 is the NMMB scheme, not WRF-ARW,
and is not applied here.

GSD visibility (``CALVIS_GSD``, field 410) runs only for
``imp_physics == 28`` and needs the cosine of the solar zenith from
``ZENSUN.f``. It is not computed.
"""

from __future__ import annotations

import numpy as np

from upp_py.constants import D608, RD
from upp_py.wrf import WrfState

# MDLFLD.f
RAINCON = 0.8333 * 1.1787e4
SNOCON = 0.94 * 1.4594e5

# CALVIS.f
_COEFLC = 144.7
_COEFLP = 2.24
_COEFFC = 327.8
_COEFFP = 10.36
_EXPONLC = 0.8800
_EXPONLP = 0.7500
_EXPONFC = 1.0000
_EXPONFP = 0.7776
_CONST1 = -np.log(0.02)
_RHOICE = 917.0
_RHOWAT = 1000.0
_VIS_MAX_M = 1.0e3 * 24.135


def convective_lowest_hydrometeors(state: WrfState) -> tuple[np.ndarray, np.ndarray]:
    """Lowest-level rain and snow mixing ratios after the convective term.

    Returns the unmodified ``QRAIN`` and ``QSNOW`` when ``RAINCV`` is 0.
    A missing ``SR`` is treated as 0, so the convective step is all rain.
    """
    qr = np.array(state.qr[0], dtype=np.float64, copy=True)
    qs = np.array(state.qs[0], dtype=np.float64, copy=True)
    dt = float(state.dt)
    if dt <= 0.0:
        return qr, qs
    cprate = np.asarray(state.raincv, dtype=np.float64) * 0.001
    sr = np.zeros(qr.shape, dtype=np.float64) if state.sr is None else np.asarray(state.sr, dtype=np.float64)
    active = cprate > 0.0
    if not np.any(active):
        return qr, qs
    t = state.t[0]
    p = state.p_mid[0]
    rdt = 1.0 / dt
    rainrate = np.maximum((1.0 - sr) * cprate * rdt, 0.0)
    # TERM1 and TERM2 in MDLFLD.f both use the lowest mass level when LMH = LM.
    tp = t / p
    term = tp**0.4167 * tp**0.5833 * rainrate**0.8333
    qr = np.where(active, qr + RAINCON * term, qr)
    snow = active & (sr > 0.0)
    if np.any(snow):
        snorate = np.maximum(sr * cprate * rdt, 0.0)
        term_s = tp**0.47 * tp**0.53 * snorate**0.94
        qs = np.where(snow, qs + SNOCON * term_s, qs)
    return qr, qs


def visibility_m(
    q_specific: np.ndarray,
    qc: np.ndarray,
    qr: np.ndarray,
    qi: np.ndarray,
    qs: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
) -> np.ndarray:
    """Horizontal visibility (m), capped at 24135 m.

    Extinction (km-1) for mass concentration ``C`` in g m-3 follows
    Kunkel (cloud water), Marshall-Palmer (rain), Rutledge and Hobbs
    (cloud ice) and Stallabrass (snow), as coded in ``CALVIS.f``.
    ``vis = -ln(0.02) / beta``, then converted to metres. The cap is
    24.135 km (Chuang, matching the observation ceiling).
    """
    q = np.asarray(q_specific, dtype=np.float64)
    temperature = np.asarray(temperature, dtype=np.float64)
    pressure = np.asarray(pressure, dtype=np.float64)
    tv = temperature * (1.0 + D608 * q)
    rho = pressure / (RD * tv)
    qclw = np.maximum(np.asarray(qc, dtype=np.float64), 0.0)
    qrain = np.maximum(np.asarray(qr, dtype=np.float64), 0.0)
    qice = np.maximum(np.asarray(qi, dtype=np.float64), 0.0)
    qsnow = np.maximum(np.asarray(qs, dtype=np.float64), 0.0)
    vovermd = (1.0 + q) / rho + (qclw + qrain) / _RHOWAT + (qice + qsnow) / _RHOICE
    # g m-3. The Fortran multiplies by 1000 after dividing by specific volume.
    scale = 1000.0 / vovermd
    conclc = np.maximum(qclw * scale, 0.0)
    conclp = np.maximum(qrain * scale, 0.0)
    concfc = np.maximum(qice * scale, 0.0)
    concfp = np.maximum(qsnow * scale, 0.0)
    beta = (
        _COEFFC * concfc**_EXPONFC
        + _COEFFP * concfp**_EXPONFP
        + _COEFLC * conclc**_EXPONLC
        + _COEFLP * conclp**_EXPONLP
        + 1.0e-10
    )
    return 1.0e3 * np.minimum(24.135, _CONST1 / beta)


def surface_visibility(state: WrfState) -> np.ndarray:
    """Field 180: ``CALVIS`` on the lowest mass level, after convective rain."""
    qr, qs = convective_lowest_hydrometeors(state)
    return visibility_m(state.q[0], state.qc[0], qr, state.qi[0], qs, state.t[0], state.p_mid[0])
