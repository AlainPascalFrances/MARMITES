# -*- coding: utf-8 -*-
"""MMsoil's stress period compiled with numba (marmites_soil_nb) against the
reference, clsMMsoil.step's own body, cell by cell (WP7 speed, 2026-10-11).

The toy is built to visit every branch the kernel has: soil zones of one,
two and three horizons; three vegetation types with dormant days; irrigated
fields (crop 0 included -- the reference indexes the crop tables with -1
there, Python's last element); cells part and wholly open water; dry cells,
exfiltration, rejected infiltration, a negative exfiltration; the CRR
cascade with sinks evaporating or routed; a float32 bottom grid, which numpy
compares in float32.
"""

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF6'):
    if os.path.join(TRUNK, p) not in sys.path:
        sys.path.insert(0, os.path.join(TRUNK, p))

import MARMITESsoil_v3 as ref                                   # noqa: E402
import marmites_soil_nb as nb                                   # noqa: E402
import marmites_crr as crrmod                                   # noqa: E402
from marmites_indices import INDEX_MM, INDEX_MM_SOIL            # noqa: E402

HNOFLO = 9999.999
needs_numba = pytest.mark.skipif(nb._jit is None,
                                 reason='numba is not installed')


class _MF:
    def __init__(self, nrow, ncol, nper, rng):
        self.nrow, self.ncol, self.nlay, self.nper = nrow, ncol, 1, nper
        self.perlen = [1.0] * nper
        self.perlen[3] = 2.0                       # a longer period inside
        self.nstp = [1] * nper
        self.hnoflo = HNOFLO
        self.hdry = None                           # MF6
        self.delr = [50.0] * ncol
        self.delc = [50.0] * nrow
        self.outcropL = np.ones((nrow, ncol), dtype=int)
        self.outcropL[rng.uniform(size=(nrow, ncol)) < 0.12] = 0


def _toy(seed, botm32=True, irr_column=True, sinks='evaporate',
         nrow=7, ncol=9, nper=14):
    rng = np.random.default_rng(seed)
    cMF = _MF(nrow, ncol, nper, rng)
    mm = ref.clsMMsoil(hnoflo=HNOFLO)
    cells = mm.build_cell_list(cMF)
    st = list(mm.paramEg)
    # three soil zones, 1 / 2 / 3 horizons
    _nsl = [1, 2, 3]
    _slprop = [np.array([1.0]), np.array([0.4, 0.6]),
               np.array([0.2, 0.3, 0.5])]
    _Sm = [[0.12], [0.35, 0.40], [0.22, 0.30, 0.42]]
    _Sfc = [[0.06], [0.20, 0.26], [0.15, 0.20, 0.28]]
    _Sr = [[0.02], [0.08, 0.10], [0.05, 0.07, 0.12]]
    _Ks = [[5.0], [3.0, 1.5], [40.0, 12.0, 2.0]]
    _Ssoil_ini = [[0.05], [0.25, 0.30], [0.10, 0.31, 0.40]]
    _st = [st[1], st[-1], st[4]]
    shp = (nrow, ncol)
    gridSOIL = rng.integers(1, 4, shp)
    gridMETEO = rng.integers(1, 3, shp)
    gridSOILthick = rng.uniform(0.3, 2.5, shp).astype(np.float32)
    TopSoil = (700.0 + rng.uniform(-8, 8, shp)) * 1000.0
    botm_l0 = TopSoil / 1000.0 - gridSOILthick - rng.uniform(4, 30, shp)
    if botm32:
        botm_l0 = botm_l0.astype(np.float32)
    nveg, nm, nf, nc = 3, 2, 2, 2
    gridVEGarea = np.zeros((nveg,) + shp, dtype=np.float32)
    gridVEGarea[0] = rng.uniform(0, 60, shp)
    gridVEGarea[1] = rng.uniform(0, 30, shp)
    gridVEGarea[2] = rng.uniform(0, 10, shp)
    P = rng.choice([0.0, 0.0, 0.3, 4.0, 18.0, 60.0], (nm, nper))
    Eo = rng.uniform(0.5, 6.0, (nm, nper))
    PT = rng.uniform(0.2, 5.0, (nm, nveg, nper))
    Pe = P[:, None, :] * rng.uniform(0.6, 1.0, (nm, nveg, nper))
    LAI = rng.choice([0.0, 0.5, 2.5], (nveg, nper))
    PE = rng.uniform(0.1, 4.0, (nm, 3, nper))
    Zr = [0.4, 15.0, 3.0]
    kTg_min, kTg_max = [0.05, 0.04, 0.0], [0.1, 0.86, 0.7]
    kT_f, kT_s = [0.5, 0.1, 0.3], [1.0, 0.05, 0.2]
    gridIRR = np.where(rng.uniform(size=shp) < 0.25,
                       rng.integers(1, nf + 1, shp), 0)
    crop = rng.integers(0, nc + 1, (nf, nper))
    P_irr = P[:, None, :] + rng.choice([0.0, 0.0, 12.0, 25.0], (nm, nf, nper))
    Pe_irr = P_irr * rng.uniform(0.7, 1.0, (nm, nf, nper))
    PT_irr = rng.uniform(1.0, 6.0, (nm, nf, nper))
    Zr_c, kTg_min_c, kTg_max_c = [0.6, 1.2], [0.05, 0.02], [0.1, 0.3]
    kT_f_c, kT_s_c = [0.5, 0.4], [1.0, 0.5]
    ctx = mm.build_context(
        cMF, cells, _nsl, 3, _st, _Sm, _Sfc, _Sr, _slprop, _Ssoil_ini,
        botm_l0, _Ks, gridSOIL, gridSOILthick, TopSoil, gridMETEO,
        INDEX_MM, INDEX_MM_SOIL, P, Eo, PT, Pe, PE, gridVEGarea, LAI, Zr,
        kTg_min, kTg_max, kT_f, kT_s, nveg, 1000.0, 1, P_irr, PT_irr,
        Pe_irr, crop, gridIRR, Zr_c, kTg_min_c, kTg_max_c, kT_f_c, kT_s_c)
    n = len(cells)
    ctx.gw_evt = True
    ctx.irr_column = irr_column
    fo = np.where(rng.uniform(size=n) < 0.15, rng.uniform(0.05, 0.95, n), 0.0)
    fo[rng.uniform(size=n) < 0.05] = 1.0
    ctx.f_open = fo
    # the cascade: 4-neighbours among the active cells
    pos = {(c[1], c[2]): c[0] for c in cells}
    nbrs = []
    for c in cells:
        i, j = c[1], c[2]
        nbrs.append(np.array([pos[q] for q in ((i - 1, j), (i + 1, j),
                                               (i, j - 1), (i, j + 1))
                              if q in pos], dtype=int))
    elev = np.array([TopSoil[c[1], c[2]] / 1000.0 for c in cells])
    xy = np.array([(c[2] * 50.0, -c[1] * 50.0) for c in cells])
    kind = np.where(fo >= 0.5, crrmod.KIND_SFR, crrmod.KIND_SOIL)
    crr = crrmod.CascadeNetwork(elev, nbrs, xy, kind, beta=0.9, sinks=sinks)
    forcing = []
    for t in range(nper):
        base = np.array([ctx.botm_l0[c[1], c[2]] for c in cells], float)
        soilbase = np.array([(TopSoil[c[1], c[2]] / 1000.0
                              - gridSOILthick[c[1], c[2]]) for c in cells])
        h = np.where(rng.uniform(size=n) < 0.1, base - 0.5,
                     soilbase - rng.uniform(-0.3, 6.0, n))
        exf = np.where(h > soilbase, rng.uniform(0.0, 3.0, n), 0.0)
        exf[rng.uniform(size=n) < 0.02] = -0.01
        rej = np.where(rng.uniform(size=n) < 0.1, rng.uniform(0, 2.0, n), 0.0)
        etu = rng.uniform(0.0, 1.5, n)
        # what UZF can take (MF6Coupler.perc_cap, 2026-10-18): a cap in
        # some cells, none in the others
        cap = np.where(rng.uniform(size=n) < 0.3, rng.uniform(0.0, 1.5, n),
                       np.inf)
        forcing.append((h, exf, rej, etu, cap))
    return mm, ctx, crr, forcing


def _run(mm, ctx, crr, forcing, compiled, monkeypatch):
    if compiled:
        monkeypatch.delenv('MM_SOIL_PYTHON', raising=False)
    else:
        monkeypatch.setenv('MM_SOIL_PYTHON', '1')
    ref._TALLIES.clear()
    state = mm.init_state(ctx)
    days = []
    for t, (h, exf, rej, etu, cap) in enumerate(forcing):
        out = mm.step(ctx, t, t, h, exf, state, rejinf_cell=rej,
                      etuzf_cell=etu, crr=crr, perc_cap=cap)
        days.append((out, state.Ssoil_ini.copy()))
    tallies = {k: dict(v) for k, v in ref._TALLIES.items()}
    ref._TALLIES.clear()
    return days, tallies


def _same(a, b, what):
    np.testing.assert_allclose(np.asarray(a, float), np.asarray(b, float),
                               rtol=1e-12, atol=1e-12, err_msg=what)


@needs_numba
@pytest.mark.parametrize('seed,botm32,irr_column,sinks', [
    (0, True, True, 'evaporate'), (1, False, False, 'route'),
    (2, True, False, 'route'), (3, False, True, 'evaporate')])
def test_the_compiled_period_is_the_reference_period(monkeypatch, capsys,
                                                     seed, botm32,
                                                     irr_column, sinks):
    mm, ctx, crr, forcing = _toy(seed, botm32, irr_column, sinks)
    a, ta = _run(mm, ctx, crr, forcing, False, monkeypatch)
    said_ref = capsys.readouterr().out
    b, tb = _run(mm, ctx, crr, forcing, True, monkeypatch)
    said_nb = capsys.readouterr().out
    assert getattr(ctx, '_nb_pack', None) is not None, \
        'the kernel did not run: the reference ran twice'
    for t, ((ra, sa), (rb, sb)) in enumerate(zip(a, b)):
        for k in ('MM', 'MM_S', 'perc', 'etg', 'petuzf', 'ro_deliver'):
            _same(rb[k], ra[k], 'day %d %s' % (t, k))
        for k in ('eg_pe', 'land', 'h0', 'tg_rate', 'tg_tip'):
            _same(rb['gw'][k], ra['gw'][k], 'day %d gw %s' % (t, k))
        shah = np.array([[p['y0'], p['b'], p['dll'], p['ext_d']]
                         for p in ra['gw']['shah']])
        _same(rb['gw']['shah'], shah, 'day %d gw shah' % t)
        for k in ('ro', 'runon', 'deliver', 'ecrr', 'ecrr_sink', 'routed'):
            _same(getattr(rb['crr'], k), getattr(ra['crr'], k),
                  'day %d crr %s' % (t, k))
        _same(sb, sa, 'day %d carried soil state' % t)
    assert ta.keys() == tb.keys()
    for k in ta:
        assert ta[k]['n'] == tb[k]['n'], k
        assert ta[k]['lo'] == pytest.approx(tb[k]['lo'], rel=1e-12), k
        assert ta[k]['hi'] == pytest.approx(tb[k]['hi'], rel=1e-12), k
    assert said_ref.count('EXFg < 0.0') == said_nb.count('EXFg < 0.0') > 0


@needs_numba
def test_the_compiled_period_without_a_cascade(monkeypatch):
    mm, ctx, _crr, forcing = _toy(5)
    a, _ = _run(mm, ctx, None, forcing, False, monkeypatch)
    b, _ = _run(mm, ctx, None, forcing, True, monkeypatch)
    for t, ((ra, sa), (rb, sb)) in enumerate(zip(a, b)):
        for k in ('MM', 'MM_S', 'perc', 'etg', 'petuzf'):
            _same(rb[k], ra[k], 'day %d %s' % (t, k))
        assert 'crr' not in rb and 'ro_deliver' not in rb
        _same(sb, sa, 'day %d state' % t)


def test_a_replaced_forcing_is_packed_again():
    """The pack is cached on the context; a forcing array put in its place
    (a test does exactly this with the LAI) must not be served stale."""
    mm, ctx, crr, _forcing = _toy(8)
    a = nb.pack(mm, ctx, crr)
    assert nb.pack(mm, ctx, crr) is a
    ctx.LAI_veg_zonesSP = np.asarray(ctx.LAI_veg_zonesSP) * 0.0
    b = nb.pack(mm, ctx, crr)
    assert b is not a and not b.LAI.any()


def test_the_uncoupled_path_stays_the_reference(monkeypatch):
    """No groundwater ET by EVT (the uncoupled run): the kernel declines."""
    mm, ctx, crr, forcing = _toy(6)
    ctx.gw_evt = False
    assert nb.pack(mm, ctx, crr) is None


@needs_numba
@pytest.mark.parametrize('compiled', [False, True])
def test_a_non_finite_storage_stops_the_period(monkeypatch, compiled):
    """A NaN in the carried soil water: both paths refuse, in words."""
    mm, ctx, crr, forcing = _toy(7)
    if compiled:
        monkeypatch.delenv('MM_SOIL_PYTHON', raising=False)
    else:
        monkeypatch.setenv('MM_SOIL_PYTHON', '1')
    state = mm.init_state(ctx)
    h, exf, rej, etu, _cap = forcing[0]
    mm.step(ctx, 0, 0, h, exf, state, rejinf_cell=rej, etuzf_cell=etu,
            crr=crr)
    state.Ssoil_ini[:, :] = np.nan
    h, exf, rej, etu, _cap = forcing[1]
    with pytest.raises(ref.MarmitesSoilError) as exc:
        mm.step(ctx, 1, 1, h, exf, state, rejinf_cell=rej, etuzf_cell=etu,
                crr=crr)
    assert 'Non-finite soil storage' in str(exc.value)
