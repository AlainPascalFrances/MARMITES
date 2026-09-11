# -*- coding: utf-8 -*-
"""Phase-1 regression tests: legacy (v0.3, Decimal) vs new (v0.4, float64)
MARMITESsoil flux() on single-cell soil columns, plus mass-balance closure
tests on paths that were broken in the legacy code (EXF > 0).
"""
import importlib.util
import os
import sys
from decimal import Decimal, ROUND_HALF_EVEN

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.join(HERE, '..')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


legacy = _load('mmsoil_legacy', os.path.join(HERE, 'legacy', 'MARMITESsoil_v3_legacy.py'))
new = _load('mmsoil_new', os.path.join(TRUNK, 'MARMITESsoil', 'MARMITESsoil_v3.py'))

HNOFLO = -999.9
HDRY = 1.0e30


class _FakeProcess:
    def __init__(self, sy):
        self._sy = sy

    def float2array(self, _):
        return np.full((1, 1, 1), self._sy)


class _FakeMF:
    """Minimal stand-in for clsMF as used inside flux()."""

    def __init__(self, perlen=(1,), sy=0.05, wel_yn=1):
        self.perlen = list(perlen)
        self.hnoflo = HNOFLO
        self.hdry = HDRY
        self.wel_yn = wel_yn
        self.outcropL = np.ones((1, 1), dtype=int)
        self.sy_actual = sy
        self.cPROCESS = _FakeProcess(sy)


def _q(x):
    return Decimal(float(x)).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN)


def _column(nsl=2):
    """A loam column: parameters per soil layer (fractions), Tl in mm."""
    Sm = [0.41, 0.41][:nsl]
    Sfc = [0.30, 0.30][:nsl]
    Sr = [0.05, 0.05][:nsl]
    Ks = [50.0, 20.0][:nsl]
    Tl = [300.0, 700.0][:nsl]
    return Sm, Sfc, Sr, Ks, Tl


def _run_pair(Pe, PT, PE, Ssoil_ini_frac, EXF_ini=0.0, dgwt=5000.0,
              nsl=2,
              NVEG=2, VEGarea=(40.0, 30.0), LAI=(2.0, 1.5), Zr=(600.0, 1500.0),
              run_legacy=True):
    """Run legacy and new flux() with equivalent inputs; return dict pair.

    WP1d removed MMsoil's surface reservoir: ponding capacity and open-water
    evaporation moved to MODFLOW with the water (SFR for the channels, LAK for
    the charcas). The legacy side is therefore driven with Ssurf_max = 0 and
    Eosurf_max = 0, at which point it computes exactly what the new code does
    -- so this comparison is kept, and now proves the removal is equivalent to
    the old model with no pond capacity rather than a change of physics.
    """
    Eosurf_max = Ssurf_max = Ssurf_ini = 0.0
    Sm, Sfc, Sr, Ks, Tl = _column(nsl)
    TopSoil = 700_000.0  # mm elevation
    TopSoilLay = np.zeros(nsl, dtype=np.float32)
    BotSoilLay = np.zeros(nsl, dtype=np.float32)
    for l in range(nsl):
        TopSoilLay[l] = TopSoil if l == 0 else BotSoilLay[l - 1]
        BotSoilLay[l] = TopSoilLay[l] - Tl[l]
    Zr_elev = [TopSoil - z for z in Zr[:NVEG]]
    HEADSini = TopSoil - dgwt
    kTg_min = [0.0] * NVEG
    kTg_max = [0.9, 0.7][:NVEG]
    kT_f = [0.5] * NVEG
    kT_s = [1.0 / 10.0] * NVEG
    st = 'loam'
    Ssoil_ini_mm = [f * t for f, t in zip(Ssoil_ini_frac, Tl)]

    results = {}

    if run_legacy:
        cMF = _FakeMF()
        out = legacy.clsMMsoil(hnoflo=HNOFLO).flux(
            cMF, 1.0, Pe, np.array(PT, dtype=np.float64), PE, _q(Eosurf_max),
            [_q(z) for z in Zr_elev], np.array(VEGarea[:NVEG], dtype=np.float32),
            HEADSini, TopSoilLay, BotSoilLay,
            [_q(t) for t in Tl], nsl, [_q(x) for x in Sm], [_q(x) for x in Sfc],
            [_q(x) for x in Sr], [_q(x) for x in Ks], _q(Ssurf_max),
            [float(x) for x in Ssoil_ini_mm], Ssurf_ini, EXF_ini, dgwt, st,
            0, 0, 0, [_q(x) for x in kTg_min], [_q(x) for x in kTg_max],
            [_q(x) for x in kT_f], [_q(x) for x in kT_s], NVEG,
            np.array(LAI[:NVEG], dtype=np.float64))
        results['legacy'] = out

    cMF = _FakeMF()
    out = new.clsMMsoil(hnoflo=HNOFLO).flux(
        cMF, 1.0, Pe, np.array(PT, dtype=np.float64), PE,
        Zr_elev, np.array(VEGarea[:NVEG], dtype=np.float64),
        HEADSini, TopSoilLay.astype(np.float64), BotSoilLay.astype(np.float64),
        Tl, nsl, Sm, Sfc, Sr, Ks,
        Ssoil_ini_mm, EXF_ini, dgwt, st,
        0, 0, 0, kTg_min, kTg_max, kT_f, kT_s, NVEG,
        np.array(LAI[:NVEG], dtype=np.float64))
    results['new'] = out
    return results


KEYS = ['Eow', 'Ssurf', 'Ro', 'Rp', 'Esoil', 'Tsoil', 'Ssoil', 'Ssoil_pc',
        'Eg', 'Tg', 'HEADS_corr', 'dgwt_corr', 'SAT', 'Rexf', 'I']
ATOL = 5e-3  # mm -- legacy quantized at 1e-5 but rounding cascades


def _compare(res):
    lg, nw = res['legacy'], res['new']
    for k, a, b in zip(KEYS, lg, nw):
        if k == 'SAT':
            assert list(np.asarray(a, dtype=bool)) == list(np.asarray(b, dtype=bool)), k
            continue
        a = np.asarray(a, dtype=np.float64).ravel()
        b = np.asarray(b, dtype=np.float64).ravel()
        assert np.allclose(a, b, atol=ATOL, rtol=1e-4), \
            f'{k}: legacy={a} new={b}'


def _mass_balance(out, Pe, EXF_ini, Ssoil_ini_frac, Ssurf_ini=0.0, perlen=1.0):
    """Closure of the whole column: In - Out - dS ~ 0 (per day)."""
    (Eow, Ssurf, Ro, Rp, Esoil, Tsoil, Ssoil, _, _Eg, _Tg, _, _, _, Rexf, I) = out
    Sm, Sfc, Sr, Ks, Tl = _column(len(np.ravel(Ssoil)))
    Ssoil_ini_mm = np.array([f * t for f, t in zip(Ssoil_ini_frac, Tl)])
    dSsoil = (np.asarray(Ssoil, dtype=np.float64) - Ssoil_ini_mm) / perlen
    # soil column balance (I and EXF in, Rexf[0], Esoil, Tsoil, Rp[-1], dS out)
    mb = (float(I) + EXF_ini / perlen) - (float(Rexf[0]) + float(np.sum(Esoil))
                                          + float(np.sum(Tsoil)) + float(np.sum(dSsoil))
                                          + float(np.ravel(Rp)[-1]))
    # surface balance
    dSsurf = (float(Ssurf) - Ssurf_ini) / perlen
    mbsurf = Pe + float(Rexf[0]) - (float(Eow) + float(Ro) + float(I) + dSsurf)
    return mb, mbsurf


# --------------------------------------------------------------------- #
# regression: legacy vs new on paths that work in both
# --------------------------------------------------------------------- #

def test_dry_no_forcing():
    res = _run_pair(Pe=0.0, PT=[0.0, 0.0], PE=0.0, Ssoil_ini_frac=[0.10, 0.10])
    _compare(res)


def test_rain_infiltration_percolation():
    res = _run_pair(Pe=15.0, PT=[2.0, 3.0], PE=4.0, Ssoil_ini_frac=[0.35, 0.32])
    _compare(res)


def test_heavy_rain_runoff():
    res = _run_pair(Pe=120.0, PT=[1.0, 1.0], PE=2.0, Ssoil_ini_frac=[0.40, 0.40])
    _compare(res)


def test_shallow_wt_eg_tg():
    res = _run_pair(Pe=0.0, PT=[3.0, 4.0], PE=5.0,
                    Ssoil_ini_frac=[0.20, 0.25], dgwt=800.0)
    _compare(res)
    assert res['new'][KEYS.index('Eg')] > 0.0
    assert res['new'][KEYS.index('Tg')] > 0.0


def test_deep_wt_no_eg():
    res = _run_pair(Pe=0.0, PT=[3.0, 4.0], PE=5.0,
                    Ssoil_ini_frac=[0.20, 0.25], dgwt=90_000.0)
    _compare(res)
    assert res['new'][KEYS.index('Eg')] == 0.0


# --------------------------------------------------------------------- #
# the EXF > 0 path: legacy is broken (bug #1 of the review), new must
# run and close the mass balance
# --------------------------------------------------------------------- #

def test_legacy_exf_crashes():
    with pytest.raises(TypeError):
        _run_pair(Pe=2.0, PT=[1.0, 1.0], PE=2.0,
                  Ssoil_ini_frac=[0.38, 0.40], EXF_ini=30.0, dgwt=100.0)


def test_new_exf_mass_balance():
    res = _run_pair(Pe=2.0, PT=[1.0, 1.0], PE=2.0,
                    Ssoil_ini_frac=[0.38, 0.40], EXF_ini=30.0, dgwt=100.0,
                    run_legacy=False)
    mb, mbsurf = _mass_balance(res['new'], Pe=2.0, EXF_ini=30.0,
                               Ssoil_ini_frac=[0.38, 0.40])
    assert abs(mb) < 1e-6, f'soil column MB = {mb}'
    assert abs(mbsurf) < 1e-6, f'surface MB = {mbsurf}'


def test_new_exf_saturates_column():
    res = _run_pair(Pe=0.0, PT=[0.0, 0.0], PE=0.0,
                    Ssoil_ini_frac=[0.40, 0.40], EXF_ini=500.0, dgwt=50.0,
                    run_legacy=False)
    SAT = res['new'][KEYS.index('SAT')]
    assert bool(np.asarray(SAT)[-1]), 'bottom layer should saturate'


def test_mass_balance_all_regression_cases():
    for kw in (dict(Pe=0.0, PT=[0.0, 0.0], PE=0.0, Ssoil_ini_frac=[0.10, 0.10]),
               dict(Pe=15.0, PT=[2.0, 3.0], PE=4.0, Ssoil_ini_frac=[0.35, 0.32]),
               dict(Pe=120.0, PT=[1.0, 1.0], PE=2.0, Ssoil_ini_frac=[0.40, 0.40])):
        res = _run_pair(run_legacy=False, **kw)
        mb, mbsurf = _mass_balance(res['new'], Pe=kw['Pe'], EXF_ini=0.0,
                                   Ssoil_ini_frac=kw['Ssoil_ini_frac'])
        assert abs(mb) < 1e-6, f'{kw}: soil MB = {mb}'
        assert abs(mbsurf) < 1e-6, f'{kw}: surface MB = {mbsurf}'


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-v']))
