# -*- coding: utf-8 -*-
"""WP5 -- the CRR cascade: routing, reinfiltration, mass balance.

Daoud et al. (2022) Eq. 23 on the MARMITES cell list: runoff moves to the
lower face neighbours in proportion to the slope, the fractions sum to beta,
the rest evaporates; a cell with no lower neighbour is a sink. A soil column
takes the run-on by the infiltration law of the rain, a stream or pond cell
hands it to SFR / LAK. Cookbook 5.1-5.7.
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

import marmites_crr as crr                                     # noqa: E402
from marmites_crr import KIND_LAK, KIND_SFR, KIND_SOIL         # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


T = _load('t_runmmsoil_crr', os.path.join(HERE, 'test_runmmsoil.py'))
coup = _load('marmites_coupler_crr', os.path.join(TRUNK, 'marmites_coupler.py'))
IX = T.INDEX_MM


def _grid(z, kind=None, beta=1.0, sinks='evaporate', d=100.0):
    """A structured network on the elevations ``z`` (nrow, ncol)."""
    z = np.asarray(z, dtype=float)
    nrow, ncol = z.shape
    topo = crr.structured_topology([d] * ncol, [d] * nrow)
    k = (np.full(z.size, KIND_SOIL) if kind is None
         else np.asarray(kind, dtype=int).ravel())
    return crr.network_from_topology(topo, np.arange(z.size), z.ravel(), k,
                                     beta=beta, sinks=sinks)


def _own(own, cap=None):
    """A stand-in cell: its own runoff ``own`` [mm/d], and a soil that takes
    up to ``cap`` [mm/d] of the run-on. Records what each cell received."""
    seen = {}

    def fn(k, runon):
        seen[k] = runon
        took = 0.0 if cap is None else min(runon, cap[k])
        return own[k] + runon - took
    return fn, seen


PIT = [[5.0, 4.0, 5.0],
       [4.0, 1.0, 4.0],
       [5.0, 4.0, 5.0]]


# ------------------------------- the weights ---------------------------- #

def test_weights_sum_to_beta_and_follow_the_slope():
    z = [[9.0, 8.0, 7.0],
         [8.0, 6.0, 4.0],
         [7.0, 4.0, 1.0]]
    net = _grid(z, beta=0.85)
    for i in range(9):
        if net.sink[i]:
            continue
        assert net.alpha[i].sum() == pytest.approx(0.85, abs=1e-15)
    # the centre (6) drains to its lower neighbours 4 and 4 (right, down),
    # equally: same drop over the same distance
    c = 4
    assert sorted(net.receivers[c].tolist()) == [5, 7]
    assert net.alpha[c].tolist() == pytest.approx([0.425, 0.425])
    # cell 1 (8) -> 2 (7) drop 1, 4 (6) drop 2: weights 1:2
    r = dict(zip(net.receivers[1].tolist(), net.alpha[1].tolist()))
    assert r == pytest.approx({2: 0.85 / 3, 4: 0.85 * 2 / 3})


def test_only_face_neighbours_receive():
    """No flow through a corner -- MODFLOW connects faces only: a cell that
    touches the low one only at a corner, between two flat ones, is a sink."""
    z = [[9.0, 9.0], [9.0, 1.0]]
    net = _grid(z)
    assert net.sink.tolist() == [True, False, False, True]
    assert net.receivers[1].tolist() == [3]
    assert net.receivers[2].tolist() == [3]


def test_the_order_puts_every_receiver_after_its_donors():
    rng = np.random.default_rng(3)
    z = rng.uniform(0, 50, size=(5, 5))
    net = _grid(z)
    assert sorted(net.order.tolist()) == list(range(25))
    rank = np.empty(25, int)
    rank[net.order] = np.arange(25)
    for i in range(25):
        assert np.all(rank[net.receivers[i]] > rank[i])


def test_a_receiver_visited_before_its_donor_is_refused():
    net = _grid(PIT)
    net.receivers[4] = np.array([0])       # the pit "drains" up to a corner
    with pytest.raises(crr.CRRError, match='visited before'):
        net._check_order()


def test_bad_settings_are_refused():
    with pytest.raises(crr.CRRError, match='beta'):
        _grid(PIT, beta=0.0)
    with pytest.raises(crr.CRRError, match='sinks'):
        _grid(PIT, sinks='spill')
    with pytest.raises(crr.CRRError, match='elevation'):
        _grid([[1.0, np.nan]])


# -------------------------------- the sinks ----------------------------- #

def test_the_pit_is_the_one_sink_and_its_water_evaporates():
    net = _grid(PIT)
    assert net.sink.tolist() == [False] * 4 + [True] + [False] * 4
    assert any('1 topographic sink(s)' in s for s in net.summary())
    assert any('0.0 % in the streams' in s and '(100.0 % at sinks)' in s
               for s in net.summary(np.full(9, 1e4)))
    fn, _ = _own(np.full(9, 1.0))
    res = net.route(np.full(9, 1e4), fn)
    # every drop ends in the pit (beta = 1), and evaporates there
    assert res.ecrr_sink.sum() == pytest.approx(9 * 1.0 * 1e4 / 1000.0)
    assert res.ecrr.sum() == pytest.approx(res.ecrr_sink.sum())
    assert not res.deliver.any()


def test_route_sends_a_sink_to_the_nearest_open_water():
    kind = np.full(9, KIND_SOIL)
    kind[8] = KIND_SFR                     # a reach in the far corner
    net = _grid(PIT, kind=kind, sinks='route', beta=0.9)
    assert net.route_to[4] == 8
    fn, _ = _own(np.where(np.arange(9) == 4, 2.0, 0.0))
    res = net.route(np.full(9, 1e4), fn)
    assert res.routed[4] == pytest.approx(0.9 * 20.0)
    assert res.deliver[8] == pytest.approx(0.9 * 20.0)
    assert res.ecrr[4] == pytest.approx(0.1 * 20.0)


def test_an_edge_sink_is_reported_as_one():
    """A cell whose only way down leaves the catchment."""
    z = np.array([[3.0, 2.0, 1.0]])
    topo = crr.structured_topology([100.0] * 3, [100.0])
    # only the first two cells are active: the lowest one is outside
    net = crr.network_from_topology(topo, [0, 1], z.ravel()[:2],
                                    [KIND_SOIL, KIND_SOIL])
    assert net.sink.tolist() == [False, True]
    assert net.edge.tolist() == [True, True]
    assert any('1 of them on the catchment edge' in s for s in net.summary())


# ----------------------------- mass balance ----------------------------- #

@pytest.mark.parametrize('shape,beta', [((3, 3), 1.0), ((5, 5), 0.8),
                                        ((5, 5), 1.0)])
def test_closure_to_machine_precision(shape, beta):
    """Routed in = routed out + reinfiltrated + evaporated (cookbook 5.7),
    on a DEM with a known outlet (a stream cell) and a known sink."""
    rng = np.random.default_rng(7)
    nrow, ncol = shape
    yy, xx = np.mgrid[0:nrow, 0:ncol]
    z = 10.0 + 2.0 * xx + 1.0 * yy + rng.uniform(0, 0.5, size=shape)
    z[nrow // 2, ncol // 2] -= 20.0                     # the sink
    kind = np.full(z.shape, KIND_SOIL)
    kind[0, 0] = KIND_SFR                               # the outlet
    if nrow > 3:
        kind[0, 1] = KIND_LAK
    net = _grid(z, kind=kind, beta=beta)
    n = z.size
    area = rng.uniform(5e3, 2e4, size=n)
    own = rng.uniform(0.0, 5.0, size=n)
    cap = rng.uniform(0.0, 2.0, size=n)
    fn, seen = _own(own, cap)
    res = net.route(area, fn)
    out = res.ro.sum()
    dest = res.runon.sum() + res.deliver.sum() + res.ecrr.sum()
    assert abs(out - dest) <= 1e-12 * out
    # generated = reinfiltrated + delivered + evaporated
    reinf = sum(min(seen[k], cap[k]) * area[k] / 1000.0 for k in range(n))
    gen = float(np.sum(own * area / 1000.0))
    assert gen == pytest.approx(reinf + res.deliver.sum() + res.ecrr.sum(),
                                rel=1e-12)
    # the sink got water, and the outlet delivered
    assert res.ecrr_sink[(nrow // 2) * ncol + ncol // 2] > 0.0
    assert res.deliver[0] > 0.0


def test_beta_evaporates_on_every_hop():
    """A chain of four soil cells onto a reach, no infiltration: the reach
    gets r (beta + beta^2 + beta^3 + beta^4)."""
    z = [[5.0, 4.0, 3.0, 2.0, 1.0]]
    kind = [KIND_SOIL] * 4 + [KIND_SFR]
    b, r, a = 0.8, 2.0, 1e4
    net = _grid(z, kind=kind, beta=b)
    fn, _ = _own(np.array([r] * 4 + [0.0]))
    res = net.route(np.full(5, a), fn)
    q = r * a / 1000.0
    assert res.deliver[4] == pytest.approx(q * (b + b ** 2 + b ** 3 + b ** 4),
                                           rel=1e-14)
    assert res.ecrr.sum() == pytest.approx(4 * q - res.deliver[4], rel=1e-14)
    # the summary's connectivity line: (b + b^2 + b^3 + b^4) / 4 = 59.0 %
    assert any('would end 59.0 % in the streams' in s
               for s in net.summary(np.full(5, a)))


def test_open_water_takes_its_own_runoff_and_routes_nothing():
    z = [[1.0, 5.0, 3.0]]
    kind = [KIND_SOIL, KIND_LAK, KIND_SOIL]
    net = _grid(z, kind=kind)
    assert net.receivers[1].size == 0 and not net.sink[1]
    fn, _ = _own(np.array([0.0, 3.0, 0.0]))
    res = net.route(np.full(3, 1e4), fn)
    assert res.deliver.tolist() == pytest.approx([0.0, 30.0, 0.0])
    assert not res.runon.any()


# --------------------------- inside MMsoil ------------------------------ #

def _soil(nrow=3, ncol=3, P=60.0):
    cMF = T._FakeMF(nrow=nrow, ncol=ncol, nper=1, perlen=[1])
    cMF.modelname = 'toy'
    cMF.outcropL[:] = 1
    inp = T._build_inputs(cMF)
    inp['P_veg_zoneSP'][:] = P
    inp['Pe_veg_zonesSP'][:] = 0.9 * P
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    cells = mm.build_cell_list(cMF)
    ctx = mm.build_context(cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'],
                           inp['_Sm'], inp['_Sfc'], inp['_Sr'], inp['_slprop'],
                           inp['_Ssoil_ini'], inp['botm_l0'], inp['_Ks'],
                           inp['gridSOIL'], inp['gridSOILthick'], inp['TopSoil'],
                           inp['gridMETEO'], T.INDEX_MM, T.INDEX_MM_S,
                           inp['P_veg_zoneSP'], inp['Eo_zonesSP'],
                           inp['PT_veg_zonesSP'], inp['Pe_veg_zonesSP'],
                           inp['PE_zonesSP'], inp['gridVEGarea'],
                           inp['LAI_veg_zonesSP'], inp['Zr'], inp['kTg_min'],
                           inp['kTg_max'], inp['kT_f'], inp['kT_s'], inp['NVEG'],
                           1000.0, 0, None, None, None, None, None,
                           None, None, None, None, None)
    return mm, ctx, inp


def _state(mm, ctx, inp, wet):
    """Soil per cell [mm]: saturated where ``wet``, near wilting elsewhere."""
    st = mm.init_state(ctx)
    st.carried = True
    sm = np.asarray(inp['_Sm'][0])
    tl = 2000.0 * np.asarray(inp['_slprop'][0])
    for k in range(ctx.ncell):
        st.Ssoil_ini[k, :2] = (sm * tl) if wet[k] else (0.06 * tl)
    return st


def _step(mm, ctx, st, net=None, rej=None):
    h = np.full(ctx.ncell, 699.0)
    exf = np.zeros(ctx.ncell)
    return mm.step(ctx, 0, 0, h, exf, st, rejinf_cell=rej, crr=net)


def test_without_run_on_the_cascade_changes_nothing():
    """A flat field: every soil cell is a sink, nothing moves. Every flux,
    the bottom-up rejected-infiltration return included, is bit-identical to
    the run without the cascade; only the evaporated runoff is booked."""
    mm, ctx, inp = _soil()
    wet = np.ones(ctx.ncell, bool)
    rej = np.linspace(0.0, 4.0, ctx.ncell)
    a = _step(mm, ctx, _state(mm, ctx, inp, wet), rej=rej)
    net = _grid(np.zeros((3, 3)))
    assert net.sink.all()
    b = _step(mm, ctx, _state(mm, ctx, inp, wet), net, rej=rej)
    for k in ('perc', 'etg', 'petuzf', 'MM_S'):
        assert np.array_equal(a[k], b[k]), k
    keep = [v for key, v in IX.items() if key not in ('iEcrr', 'iETtot')]
    assert np.array_equal(a['MM'][:, keep], b['MM'][:, keep])
    assert a['MM'][:, IX['iRo']].min() > 0.0         # there was runoff
    # m3/d and back: equal to rounding
    np.testing.assert_allclose(b['MM'][:, IX['iEcrr']], a['MM'][:, IX['iRo']],
                               rtol=1e-14)
    np.testing.assert_allclose(b['MM'][:, IX['iETtot']],
                               a['MM'][:, IX['iETtot']]
                               + a['MM'][:, IX['iRo']], rtol=1e-14)
    # the off path never touches the new columns
    for key in ('iRunon', 'iReinf', 'iEcrr'):
        assert not a['MM'][:, IX[key]].any()


def test_run_on_reinfiltrates_top_down_and_the_rest_runs_on():
    """Saturated top row under heavy rain, dry cells below, a reach at the
    bottom: the dry columns take run-on by Eq. 1b, never beyond their
    capacity, and what passes them reaches the stream."""
    # a dry column's top layer holds ~210 mm more; the rain brings ~140 mm,
    # the run-on another ~140: about half of it reinfiltrates
    mm, ctx, inp = _soil(P=150.0)
    z = np.array([[9.0, 9.0, 9.0], [5.0, 5.0, 5.0], [1.0, 1.0, 1.0]])
    kind = np.full((3, 3), KIND_SOIL)
    kind[2, :] = KIND_SFR
    net = _grid(z, kind=kind)
    wet = np.array([True] * 3 + [False] * 6)
    base = _step(mm, ctx, _state(mm, ctx, inp, wet))
    st = _state(mm, ctx, inp, wet)
    out = _step(mm, ctx, st, net)
    MM, MMb = out['MM'], base['MM']
    mid = [3, 4, 5]
    assert np.all(MM[mid, IX['iRunon']] > 0.0)
    assert np.all(MM[mid, IX['iReinf']] > 0.0)
    assert not MM[[0, 1, 2, 6, 7, 8], IX['iRunon']].any()
    # infiltration grew by exactly the reinfiltration
    assert MM[mid, IX['iI']] - MMb[mid, IX['iI']] == pytest.approx(
        MM[mid, IX['iReinf']], abs=1e-9)
    # the surface balance of every column still closes
    assert np.abs(MM[:, IX['iMBsurf']]).max() < 1e-9
    # no column above its capacity (the top layer)
    sm = inp['_Sm'][0][0] * 2000.0 * inp['_slprop'][0][0]
    assert st.Ssoil_ini[:, 0].max() <= sm + 1e-9
    # the catchment: what left the soil surface for good went to the stream
    # (beta = 1, no sinks)
    area = ctx.geom.area
    net_ro = float(np.sum((MM[:, IX['iRo']] - MM[:, IX['iRunon']]) * area))
    deliv = float(np.sum(out['ro_deliver'] * area))
    assert net_ro == pytest.approx(deliv, rel=1e-12)
    assert not MM[:, IX['iEcrr']].any()
    # the stream cells take their own runoff and the cascade's
    assert np.all(out['ro_deliver'][6:] >= MM[6:, IX['iRo']] - 1e-12)
    assert not out['ro_deliver'][:6].any()


# ----------------------------- in the coupler --------------------------- #

def test_the_coupler_routes_and_reports(capsys):
    """The mock coupler of test_coupler_mock with the cascade on: the MM
    vector carries it, crr_ts closes period by period, the log reports it."""
    cm = _load('t_coupler_mock_crr', os.path.join(HERE, 'test_coupler_mock.py'))
    cpl, api, ctx = cm._setup(nper=3)
    ctx.P_veg_zoneSP[:] = 300.0
    ctx.Pe_veg_zonesSP[:] = 270.0
    nrow, ncol = ctx.cMF.nrow, ctx.cMF.ncol
    elev = (np.arange(nrow)[:, None] * 2.0 + np.arange(ncol)[None, :]
            + 700.0)[::-1, ::-1]
    cpl.crr = cpl._build_crr({'beta': 0.9, 'sinks': 'evaporate',
                              'elev': elev})
    res = cpl.run(api)
    ts = res['crr_ts']
    terms = [t.decode() for t in res['crr_terms']]
    c = dict((t, ts[:, i]) for i, t in enumerate(terms))
    assert c['runoff'].sum() > 0.0 and c['evap'].sum() > 0.0
    np.testing.assert_allclose(
        c['runoff'], c['runon'] + c['to_sfr'] + c['to_lak'] + c['evap'],
        rtol=1e-12)
    log = capsys.readouterr().out
    assert 'CRR cascade (mm/yr over the catchment)' in log
    # iEcrr is in the total ET
    wb = res['wb_ts']
    assert np.all(wb[:, IX['iEcrr']] >= 0.0) and wb[:, IX['iEcrr']].sum() > 0


def test_off_the_coupler_has_no_cascade_record():
    cm = _load('t_coupler_mock_crr2', os.path.join(HERE, 'test_coupler_mock.py'))
    cpl, api, ctx = cm._setup(nper=2)
    res = cpl.run(api)
    assert cpl.crr is None and 'crr_ts' not in res
    assert not res['wb_ts'][:, IX['iRunon']].any()


def test_delivery_uses_the_cascade_when_on():
    out = {'MM': np.array([[5.0] + [0.0] * 33]),
           'ro_deliver': np.array([7.0])}
    c = SimpleNamespace(crr=None, _iRo=0)
    assert coup.MF6Coupler._runoff_to_deliver(c, out).tolist() == [5.0]
    c.crr = object()
    assert coup.MF6Coupler._runoff_to_deliver(c, out).tolist() == [7.0]
