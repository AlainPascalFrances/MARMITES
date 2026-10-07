# -*- coding: utf-8 -*-
"""Groundwater ET through EVT at the solved head -- the coupled run's only
groundwater-ET path since 2026-10-07 (the WEL route, et.gw_route = 'wel',
is gone).

MARMITES still computes Eg (Shah et al. 2007) and Tg (kTg, per vegetation
type while the water table is above its root tip); MF6 takes them through
two EVT packages at the head it solves for. Each day's curve
starts at the start-of-day head with MARMITES' rate and only falls with the
water table, so UZF's demand capped at what it leaves keeps total ET <= PET.
The API side was probed on libmf6 6.7 (E:/tmp_claude/helpers/evt_toy2.py):
SURFACE, RATE, DEPTH, PXDP and PETM written before the step are what MF6
takes, and SIMVALS is what it took.
"""
import importlib.util
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

import marmites_evt as me                                       # noqa: E402

SANDY_LOAM_FIELD = {'dll': 115.3, 'y0': 0.023, 'b': 0.013, 'ext_d': 1000.0}


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _g(depth, D, x, y):
    """MF6's segmented ET fraction ``depth`` below SURFACE."""
    return np.interp(np.asarray(depth) / D, np.r_[0.0, x, 1.0],
                     np.r_[1.0, y, 0.0])


# ------------------------------------------------------------- the curves
@pytest.mark.parametrize('d0', [0.0, 0.5, 1.0, 2.0, 5.0, 9.0])
def test_eg_follows_shah_below_the_start_of_day_head(d0):
    land, pe, nseg = 800.0, 0.003, 8
    s, R, D, x, y = me.eg_curve(pe, land, land - d0, SANDY_LOAM_FIELD, nseg)
    assert x.size == y.size == nseg - 1
    assert np.all(np.diff(x) > 0) and np.all(np.diff(y) <= 0)
    assert s == pytest.approx(land - d0)
    assert R == pytest.approx(pe * me.shah_f(d0, SANDY_LOAM_FIELD))
    dd = np.linspace(d0, 9.5, 200)
    evt = R * _g(dd - d0, D, x, y)
    assert np.max(np.abs(evt - pe * me.shah_f(dd, SANDY_LOAM_FIELD))) \
        < 0.02 * pe                                    # within 2 % of PE
    # zero at the extinction depth: no y0 step there
    assert R * _g(10.0 - d0, D, x, y) == pytest.approx(0.0)


def test_eg_above_the_land_is_measured_from_the_land():
    s, R, D, x, y = me.eg_curve(0.002, 800.0, 800.3, SANDY_LOAM_FIELD, 8)
    assert s == 800.0 and R == pytest.approx(0.002) and D == pytest.approx(10.0)


def test_eg_beyond_extinction_or_without_demand_takes_nothing():
    assert me.eg_curve(0.002, 800.0, 788.0, SANDY_LOAM_FIELD, 8)[1] == 0.0
    assert me.eg_curve(0.0, 800.0, 799.0, SANDY_LOAM_FIELD, 8)[1] == 0.0


def test_tg_ramps_out_type_by_type_at_the_root_tips():
    rates, tips = [0.001, 0.0005, 0.0002], [799.1, 784.5, 789.5]
    s, R, D, x, y = me.tg_curve(rates, tips, 799.5, 0.1, 8)
    assert s == 799.5 and R == pytest.approx(0.0017)
    assert D == pytest.approx(799.5 - 784.5)
    T = lambda h: R * _g(799.5 - h, D, x, y)               # noqa: E731
    assert T(799.5) == pytest.approx(0.0017)
    assert T(799.2) == pytest.approx(0.0017)                # above the ramp
    assert T(799.1) == pytest.approx(0.0007)                # grass is off
    assert T(789.5) == pytest.approx(0.0005)                # shrubs off
    assert T(784.5) == pytest.approx(0.0)
    # a type whose tip the start-of-day head is already below takes nothing
    assert me.tg_curve(rates, tips, 799.0, 0.1, 8)[1] == pytest.approx(0.0007)


def test_tg_needs_two_segments_per_type():
    with pytest.raises(ValueError, match='nseg'):
        me.tg_curve([1e-3] * 4, [790.0, 791.0, 792.0, 793.0], 799.0, 0.1, 6)


def test_the_curve_never_exceeds_the_start_of_day_rate():
    """Flat above SURFACE: a water table rising during the day takes no
    more than was reserved from UZF's demand."""
    s, R, D, x, y = me.eg_curve(0.003, 800.0, 798.0, SANDY_LOAM_FIELD, 8)
    assert np.all(y <= 1.0)
    s, R, D, x, y = me.tg_curve([1e-3], [790.0], 795.0, 0.1, 8)
    assert np.all(y <= 1.0)


# ------------------------------------------------------------- MMsoil
T_ = _load('t_runmmsoil_evt', os.path.join(HERE, 'test_runmmsoil.py'))
IX = T_.INDEX_MM


def _soil(coupled):
    """One MMsoil day: ``coupled`` hands Eg/Tg over to EVT (ctx.gw_evt, as
    MF6Coupler sets it); uncoupled, MMsoil takes them itself."""
    cMF = T_._FakeMF(nrow=1, ncol=2, nper=1, perlen=[1], sy=0.2)
    cMF.modelname = 'toy'
    cMF.outcropL[:] = 1
    inp = T_._build_inputs(cMF)
    mm = T_.new.clsMMsoil(hnoflo=T_.HNOFLO)
    cells = mm.build_cell_list(cMF)
    ctx = mm.build_context(
        cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'], inp['_Sm'],
        inp['_Sfc'], inp['_Sr'], inp['_slprop'], inp['_Ssoil_ini'],
        inp['botm_l0'], inp['_Ks'], inp['gridSOIL'], inp['gridSOILthick'],
        inp['TopSoil'], inp['gridMETEO'], T_.INDEX_MM, T_.INDEX_MM_S,
        inp['P_veg_zoneSP'], inp['Eo_zonesSP'], inp['PT_veg_zonesSP'],
        inp['Pe_veg_zonesSP'], inp['PE_zonesSP'], inp['gridVEGarea'],
        inp['LAI_veg_zonesSP'], inp['Zr'], inp['kTg_min'], inp['kTg_max'],
        inp['kT_f'], inp['kT_s'], inp['NVEG'], 1000.0, 0, None, None, None,
        None, None, None, None, None, None, None)
    ctx.gw_evt = bool(coupled)
    st = mm.init_state(ctx)
    out = mm.step(ctx, 0, 0, np.full(2, 698.0), np.zeros(2), st)
    return mm, ctx, out


def test_mmsoil_hands_over_the_potential_and_applies_nothing():
    _mm, ctx, own = _soil(False)
    _mm, ctx, evt = _soil(True)
    assert 'gw' not in own
    g = evt['gw']
    # nothing applied by MMsoil: no Eg/Tg in its own vector
    assert not evt['etg'].any()
    for key in ('iEg', 'iTg', 'iETg'):
        assert not evt['MM'][:, IX[key]].any()
    # the rest of the soil balance is untouched by who takes groundwater ET
    for key in ('iETsoil', 'iperc', 'iRo', 'iI', 'iPETuzf'):
        assert np.array_equal(evt['MM'][:, IX[key]], own['MM'][:, IX[key]])
    # what EVT starts the day with is what the uncoupled MMsoil takes itself
    # (no drawdown limit binds here: the table is 2 m down, 0.2 Sy)
    st_ = ctx._st[0]
    shah = T_.new.clsMMsoil(hnoflo=T_.HNOFLO).paramEg[st_]
    for k in range(ctx.ncell):
        r_eg = me.eg_curve(g['eg_pe'][k], g['land'][k], g['h0'][k], shah, 8)[1]
        r_tg = me.tg_curve(g['tg_rate'][k], g['tg_tip'][k], g['h0'][k], 0.1,
                           8)[1]
        assert r_eg * 1000.0 == pytest.approx(own['MM'][k, IX['iEg']],
                                              rel=1e-9)
        assert r_tg * 1000.0 == pytest.approx(own['MM'][k, IX['iTg']],
                                              rel=1e-9)
    assert own['MM'][:, IX['iETg']].sum() > 0.0          # there was some


# ------------------------------------------------------------- the coupler
coup = _load('marmites_coupler_evt', os.path.join(TRUNK, 'marmites_coupler.py'))
CM = _load('t_coupler_mock_evt', os.path.join(HERE, 'test_coupler_mock.py'))


def test_the_coupler_writes_caps_and_reads_back():
    """The mock API's EVT takes what is written -- as at a head above
    SURFACE, so SIMVALS = -RATE x area."""
    cpl0, _api, ctx = CM._setup(nper=3, heads0=698.0)
    mf6b = cpl0.mf6b
    mf6b.evt_nseg, mf6b.evt_ramp = 8, 0.1
    cpl = coup.MF6Coupler(cpl0.mm, ctx, cpl0.mm.init_state(ctx), mf6b,
                          conv_fact=1000.0)
    assert ctx.gw_evt, 'the coupler hands Eg/Tg over to EVT'
    api = CM.FakeApi('toy', ctx.cMF.nlay, ctx.cMF.nrow, ctx.cMF.ncol,
                     ctx.ncell, nuzf=ctx.ncell + 3, heads0=698.0)
    api.evt_area = cpl.area
    caps = []
    real = cpl._write_fluxes

    def spy(*a, **k):
        caps.append(k.get('etg_cap'))
        return real(*a, **k)
    cpl._write_fluxes = spy
    res = cpl.run(api)
    # no well anywhere: WEL is left for real pumping
    assert not np.any(api.Q) and not np.any(api.BOUND)
    # UZF's demand left room for exactly what EVT could take
    transient = [c for c in caps if c is not None]
    assert len(transient) == ctx.cMF.nper
    for c, r in zip(transient, api.evt_at_advance[1:]):
        assert np.allclose(c, r['EVT_EG'] + r['EVT_TG'])
    # ... and what MF6 took is in the books, Eg and Tg apart
    wb = res['wb_ts']
    assert wb[:, IX['iETg']].sum() > 0.0
    assert np.allclose(wb[:, IX['iETg']], wb[:, IX['iEg']] + wb[:, IX['iTg']])
    last = api.evt_at_advance[-1]
    eg = last['EVT_EG'] * 1000.0
    assert np.average(eg, weights=cpl.area) == pytest.approx(
        wb[-1, IX['iEg']], rel=1e-9)
    assert np.allclose(res['etg'][-1],
                       last['EVT_EG'] + last['EVT_TG'])


def test_evt_is_the_only_route():
    cpl, api, ctx = CM._setup(nper=2, heads0=698.0)
    res = cpl.run(api)
    assert set(cpl.p_evt) == {'eg', 'tg'}
    assert res['etg'].sum() > 0.0
    assert not hasattr(cpl, 'gw_route')
    assert 'p_evt_eg_sim' in coup.MF6Coupler.SUBSTEP_RATES


# ------------------------------------------------------- config and driver
def test_the_route_is_retired_and_the_curves_stay_panel_questions():
    import marmites_config as mcfg
    c = mcfg.RunConfig.from_dict({})
    assert not hasattr(c.et, 'gw_route')
    # an old file still loads: the key is dropped and the drop reported
    for old in ('wel', 'evt'):
        c = mcfg.RunConfig.from_dict({'et': {'gw_route': old}})
        assert any('et.gw_route is gone' in m for m in c.migrated), c.migrated
    with pytest.raises(mcfg.ConfigError, match='evt_nseg'):
        mcfg.RunConfig.from_dict({'et': {'evt_nseg': 2}})
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    for piece in ('b.evt_nseg = int(cfg.et.evt_nseg)',
                  'b.evt_ramp = float(cfg.et.evt_ramp)',
                  "DISCHARGE_PREFIXES = ('DRN', 'GHB', 'WEL', 'EVT')"):
        assert piece in src, piece
    assert 'gw_route' not in src


# ------------------------------------------------------------- the builder
def test_the_builder_makes_two_evt_packages_in_cell_order(tmp_path):
    lk = _load('t_lak_ponds_evt', os.path.join(HERE, 'test_lak_ponds.py'))
    make = lk._lamata_lak(tmp_path)
    b = make('evt')
    b.lak_shapefile = None
    b.evt_nseg = 8
    b.build()
    assert b.evt_packages == ['evt_eg', 'evt_tg']
    for pname in b.evt_packages:
        p = b.gwf.get_package(pname)
        assert p.nseg.get_data() == 8
        spd = p.stress_period_data.get_data(0)
        assert len(spd) == len(b.surf_cells)
        # one record per land column, in the cell order MMsoil uses
        assert [tuple(int(v) for v in r['cellid']) for r in spd[:3]] == \
            [tuple(int(v) for v in b._cellid(k, i, j))
             for (i, j, k) in b.surf_cells[:3]]
        assert not np.any(spd['rate'])
    # the ETg wells are gone: WEL is left for real pumping
    assert b.gwf.get_package('wel') is None


def test_a_root_tip_within_the_ramp_of_the_head_is_a_straight_line():
    """2026-10-05, the first EVT run: one type reaching the table, its tip
    0.05 m below the start-of-day head -- no inner breakpoint at all, and
    the curve builder failed on the empty array."""
    s, R, D, x, y = me.tg_curve([0.0, 0.002, 0.0], [700.0, 799.95, 700.0],
                                800.0, 0.1, 8)
    assert R == pytest.approx(0.002) and D == pytest.approx(0.05)
    assert x.size == y.size == 7
    assert np.all(np.diff(x) > 0) and np.all(np.diff(y) <= 0)
    assert y == pytest.approx(1.0 - x)                 # a straight ramp
    # and the same for Eg: a head right at the extinction depth's edge
    s, R, D, x, y = me.eg_curve(0.002, 800.0, 790.0001, SANDY_LOAM_FIELD, 8)
    assert x.size == 7 and np.all(np.diff(x) > 0)


@pytest.mark.parametrize('seed', range(5))
def test_any_cell_gets_a_valid_curve(seed):
    """Random heads, tips and rates around La Mata's: never an error, always
    nseg - 1 strictly increasing depths and non-increasing fractions."""
    rng = np.random.default_rng(seed)
    for _ in range(400):
        land = 800.0
        h0 = land - rng.uniform(-0.5, 12.0)
        tips = land - rng.uniform(0.0, 16.0, 3)
        rates = rng.choice([0.0, 1e-4, 1e-3], 3)
        for s, R, D, x, y in (me.tg_curve(rates, tips, h0, 0.1, 8),
                              me.eg_curve(rng.uniform(0, 3e-3), land, h0,
                                          SANDY_LOAM_FIELD, 8)):
            assert x.size == y.size == 7 and D > 0
            assert np.all(np.diff(x) > 0) and np.all((x > 0) & (x < 1))
            assert np.all(np.diff(y) <= 1e-15) and np.all((y >= 0) & (y <= 1))
