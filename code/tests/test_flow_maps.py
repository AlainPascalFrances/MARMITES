# -*- coding: utf-8 -*-
"""The maps added on 2026-10-10 (user): the streamflow along the network
and the exchange with the streams and the ponds -- since 2026-10-12 inside
the recharge maps: Rsw (the beds' leakage in), Exf (seepage + GW->SFR +
GW->LAK), Re = Rg + Rsw - Exf, Rn = Re - ETg."""
import os
import shutil
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import marmites_postprocess as PP                    # noqa: E402


def test_the_recharge_maps_put_the_streams_and_ponds_in():
    """User, 2026-10-12: MF6 hands the valley's groundwater to SFR and LAK
    directly. Rsw = the beds' leakage in; Exf = seepage + GW->SFR + GW->LAK
    (the OUTFLOW parts, kept apart from the inflow cell by cell); Re = Rg +
    Rsw - Exf; Rn = Re - ETg. Lateral drains and pumping are not in them."""
    nlay, n = 1, 3
    z = np.zeros((nlay, n))
    maps = {'UZF-GWRCH': z + [[2.0, 0.0, 1.0]],
            'DRN_SEEP': z + [[-0.5, 0.0, 0.0]],
            'EVT_EG': z + [[-0.2, -0.1, 0.0]],
            'EVT_TG': z + [[-0.1, 0.0, 0.0]],
            # cell 2: a stream gaining 3.0 overall, losing 0.5 at times
            'SFR': z + [[0.0, -2.5, 0.0]], 'SFR_IN': z + [[0.0, 0.5, 0.0]],
            'SFR_OUT': z + [[0.0, -3.0, 0.0]],
            'LAK': z + [[0.0, 0.0, 0.4]], 'LAK_IN': z + [[0.0, 0.0, 0.4]],
            'LAK_OUT': z,
            'DRN': z + [[-9.0, -9.0, -9.0]]}                   # lateral
    rm = PP.recharge_maps(maps, nlay)
    assert np.allclose(rm['Rsw'], [[0.0, 0.5, 0.4]])
    assert np.allclose(rm['EXFg'], [[0.5, 3.0, 0.0]])
    assert np.allclose(rm['Re'], [[1.5, -2.5, 1.4]])
    assert np.allclose(rm['Rn'], [[1.2, -2.6, 1.4]])
    # Re is the net exchange at the top, whatever the split
    top = sum(np.asarray(maps[k]) for k in ('UZF-GWRCH', 'DRN_SEEP', 'SFR',
                                           'LAK'))
    assert np.allclose(rm['Re'], top)
    # a digest made before the split: the cell's net, clipped
    old = {k: v for k, v in maps.items() if not k.endswith(('_IN', '_OUT'))}
    rm = PP.recharge_maps(old, nlay)
    assert np.allclose(rm['EXFg'], [[0.5, 2.5, 0.0]])
    assert np.allclose(rm['Re'], [[1.5, -2.5, 1.4]])
    assert PP.recharge_maps({'DRN': z}, nlay) == {}
    stems = [s for s, _l, _g in PP.RECHARGE_MAPS]
    assert stems == ['Rsw', 'EXFg', 'Re', 'Rn']
    assert not hasattr(PP, 'upper_face_flow'), 'the upper-face map is in Rn'
    assert 'SFR' not in [k for k, *_ in PP._AQ_MAPS], 'SFRg is in Rsw/Exf' 


def _mf6_exe():
    import mm_paths
    if os.path.isfile(mm_paths.MF6_EXE):
        return mm_paths.MF6_EXE
    return shutil.which('mf6')


@pytest.fixture(scope='module')
def stream_run(tmp_path_factory):
    """Three reaches in a row over a two-layer aquifer, fed 100 m3/d, with a
    streambed so tight the exchange is negligible: each reach hands on the
    100 m3/d, and the last one hands it out of the model."""
    flopy = pytest.importorskip('flopy')
    exe = _mf6_exe()
    if not exe:
        pytest.skip('no mf6')
    ws = str(tmp_path_factory.mktemp('stream'))
    name = 'strm'
    nlay, nrow, ncol = 2, 3, 4
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 2, 1.0)] * 2,
                          time_units='DAYS')
    flopy.mf6.ModflowIms(sim, complexity='SIMPLE',
                         linear_acceleration='BICGSTAB')
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True,
                               newtonoptions='NEWTON')
    flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=20.0,
                            delc=20.0, top=100.0, botm=[90.0, 50.0])
    flopy.mf6.ModflowGwfic(gwf, strt=95.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0, k33=0.1)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.1,
                            transient={0: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(1, 0, 0), 95.0]]})
    pkg = [(r, (0, 1, r + 1), 20.0, 1.0, 0.01, 99.0 - r, 0.5, 1e-9, 0.035,
            1 if r in (0, 2) else 2, 1.0, 0) for r in range(3)]
    con = [(0, -1), (1, 0, -2), (2, 1)]
    flopy.mf6.ModflowGwfsfr(gwf, nreaches=3, packagedata=pkg,
                            connectiondata=con,
                            perioddata={0: [(0, 'INFLOW', 100.0)]},
                            save_flows=True,
                            budget_filerecord='%s.sfr.cbc' % name)
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord='%s.cbc' % name,
                           saverecord=[('BUDGET', 'ALL')])
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    if not ok:
        pytest.skip('the stream model did not run')
    return ws, name, nlay, nrow, ncol


def test_each_reach_hands_on_its_flow(stream_run):
    ws, name, nlay, nrow, ncol = stream_run
    got = PP.sfr_reach_flows(ws, name)
    assert np.allclose(got['q'], 100.0, rtol=1e-3), got['q']
    # the GWF user node of each reach: layer 1, row 2, columns 2..4
    assert list(got['node']) == [1 * ncol + c + 1 for c in (1, 2, 3)]


def test_no_sfr_budget_no_streamflow(tmp_path):
    assert PP.sfr_reach_flows(str(tmp_path), 'none') is None


def test_the_map_pass_keeps_the_stream_inflow_and_outflow_apart(
        stream_run):
    """On a real MF6 budget: SFR = SFR_IN + SFR_OUT, the parts signed as
    MODFLOW's IN and OUT (+ into the aquifer)."""
    ws, name, nlay, nrow, ncol = stream_run
    cMF = SimpleNamespace(nlay=nlay, nrow=nrow, ncol=ncol)
    maps = PP._aquifer_map_pass(ws, name, cMF, 2, cache_dir=None,
                                verbose=False)
    assert {'SFR', 'SFR_IN', 'SFR_OUT'} <= set(maps)
    assert np.all(maps['SFR_IN'] >= 0.0) and np.all(maps['SFR_OUT'] <= 0.0)
    assert np.allclose(maps['SFR'], maps['SFR_IN'] + maps['SFR_OUT'])


def _pwb():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        '_pwb_flow', os.path.join(HERE, 'plot_water_budget.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_budget_table_shows_the_streams_and_compares_like_with_like():
    """User, 2026-10-10: the streams' groundwater never passes through the
    soil column, so the table lists it on its own rows, and adds it where
    the NWT run (no streams) counted it -- in Ro and in EXFg."""
    pwb = _pwb()
    labels = ['Ro', 'EXFg', 'P']
    sw = {'gw_to_sfr': 55.0, 'sfr_to_gw': 8.0, 'gw_to_lak': 0.0,
          'lak_to_gw': 0.6, 'q_outlet': 87.0}
    lines = pwb.surface_water_lines(sw, labels, [43.6, 8.2, 361.4],
                                    [89.2, 52.6, 361.7])
    row = {l.split()[0]: l.split() for l in lines if l.strip()}
    assert row['GW->SFR'][1:4] == ['55.0', '-', '-']
    assert row['Qout'][1] == '87.0'
    assert row['Ro+netSFR'][1:4] == ['90.6', '89.2', '1.4']
    assert row['EXFg+GWsw'][1:4] == ['63.2', '52.6', '10.6']
    # a new run alone: one value per row
    alone = pwb.surface_water_lines(sw, labels, [43.6, 8.2, 361.4], None)
    assert any(l.startswith('Ro+netSFR') and '90.6' in l for l in alone)
    assert pwb.surface_water_lines(None, labels, [1, 2, 3], None) == []


def test_the_stream_terms_come_from_the_listing(stream_run):
    ws, name, nlay, nrow, ncol = stream_run
    pwb = _pwb()
    assert pwb._gwf_name(ws) == name
    sw = pwb.surface_water_terms(ws, nrow * ncol * 400.0)
    assert sw is not None
    # a streambed of 1e-9 m/d: next to nothing crosses it
    assert abs(sw['gw_to_sfr']) < 1e-3 and abs(sw['sfr_to_gw']) < 1e-3
    assert sw['gw_to_lak'] is None and sw['q_outlet'] is None


def test_the_figures_take_the_streams_per_day(stream_run):
    """The NWT-comparison figures (01-04) carry the same rows as the table,
    per day: a series per period from the listing, the like-for-like totals
    built from them against the NWT term they stand for."""
    ws, name, nlay, nrow, ncol = stream_run
    pwb = _pwb()
    s = pwb.surface_water_series(ws, nrow * ncol * 400.0, nper=2)
    assert s is not None and len(s['gw_to_sfr']) == 2
    n = 3
    new = {'wb_ts': np.zeros((n, len(pwb.INDEX_MM))),
           'sw': {'gw_to_sfr': np.array([1.0, 2.0, 3.0]),
                  'sfr_to_gw': np.array([0.5, 0.5, 0.5]),
                  'gw_to_lak': None, 'lak_to_gw': None,
                  'q_outlet': np.array([4.0, 4.0, 4.0])}}
    new['wb_ts'][:, pwb.INDEX_MM['iRo']] = 1.0
    new['wb_ts'][:, pwb.INDEX_MM['iEXFg']] = 0.2
    ref = {'ts': np.ones((n, len(pwb.INDEX_MM)))}
    a, b = pwb._series(new, ref, 'Ro+netSFR')
    assert np.allclose(a, [1.5, 2.5, 3.5]) and np.allclose(b, 1.0)
    a, b = pwb._series(new, ref, 'EXFg+GWsw')
    assert np.allclose(a, [1.2, 2.2, 3.2])
    a, b = pwb._series(new, ref, 'GW->LAK')
    assert np.allclose(a, 0.0) and b is None
    assert np.allclose(pwb._series(new, ref, 'Qout')[0], 4.0)
