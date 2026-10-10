# -*- coding: utf-8 -*-
"""The maps added on 2026-10-10 (user): the streamflow along the network,
the exchange with the streams and the ponds, and the flow across the UPPER
face of each layer -- the companion of FLF, so the valley's upward flow can
be followed out of the aquifer."""
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


def test_the_upper_face_of_layer_1_is_the_net_exchange_at_the_top():
    """+ into the aquifer: recharge in, seepage / ET / a gaining stream out;
    layer 2's upper face is layer 1's lower face, plus what is attached to
    layer 2 itself (an outcrop cell, where layer 1 is absent)."""
    nlay, n = 2, 3
    z = np.zeros((nlay, n))
    maps = {'UZF-GWRCH': z + [[2.0, 0.0, 1.0], [0.0, 0.0, 0.5]],
            'DRN_SEEP': z + [[-0.5, 0.0, 0.0], [0.0, 0.0, 0.0]],
            'EVT_EG': z + [[-0.2, -0.1, 0.0], [0.0, 0.0, -0.1]],
            'SFR': z + [[0.0, -3.0, 0.0], [0.0, 0.0, 0.0]],
            'DRN': z + [[-9.0, -9.0, -9.0], [-9.0, -9.0, -9.0]],   # lateral
            'WEL': z + [[-7.0, 0.0, 0.0], [0.0, 0.0, 0.0]],       # pumping
            'FLF': z + [[1.0, -2.5, 0.0], [0.0, 0.0, 0.0]]}
    fuf = PP.upper_face_flow(maps, nlay)
    assert np.allclose(fuf[0], [1.3, -3.1, 1.0])
    assert np.allclose(fuf[1], [1.0, -2.5, 0.4])
    assert PP.upper_face_flow({'DRN': z}, nlay) is None


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


def test_layer_2s_upper_face_is_layer_1s_lower_face(stream_run):
    """On a real MF6 budget: nothing attached to layer 2 but a CHD, so what
    crosses its top is exactly what leaves layer 1 at the bottom."""
    ws, name, nlay, nrow, ncol = stream_run
    cMF = SimpleNamespace(nlay=nlay, nrow=nrow, ncol=ncol)
    maps = PP._aquifer_map_pass(ws, name, cMF, 2, cache_dir=None,
                                verbose=False)
    assert 'FLF' in maps and 'SFR' in maps
    fuf = PP.upper_face_flow(maps, nlay)
    flf = np.asarray(maps['FLF']).reshape(nlay, -1)
    assert np.allclose(fuf[1], flf[0])
    assert np.allclose(fuf[0], np.asarray(maps['SFR']).reshape(nlay, -1)[0])


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
