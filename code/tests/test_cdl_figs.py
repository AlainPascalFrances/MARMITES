# -*- coding: utf-8 -*-
"""The CdL figures (marmites_cdl_figs, marmites_cdl_out): what they read
from the written simulation -- the ponds' rim, bottom, stage and area, the
UZF land objects, the reaches and the stream's hand-offs to the pond -- the
stretch of stream each pond's section follows, and that every page is drawn.

A small two-layer DISV model in UTM coordinates with UZF, SFR and one pond
(LAK, its table and outlet) that the stream runs through by MVR, written by
flopy, never run.
"""
import datetime
import importlib.util
import json
import os
import sys
import types

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in ('', 'MARMITESutilities', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

flopy = pytest.importorskip('flopy')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CF = _load('marmites_cdl_figs',
           os.path.join(TRUNK, 'ppMF6', 'marmites_cdl_figs.py'))
CO = _load('marmites_cdl_out_m',
           os.path.join(TRUNK, 'ppMF6', 'marmites_cdl_out.py'))
GRID = _load('marmites_grid_m', os.path.join(TRUNK, 'marmites_grid.py'))

XLL, YLL, CS, NROW, NCOL = 739300.0, 4553050.0, 20.0, 5, 6
POND = (2, 2)                     # the pond's cell; the stream runs along row 3


def _ic(i, j):
    return i * NCOL + j


def _sim(ws):
    verts, cell2d, ncpl = GRID.disv_from_structured([CS] * NCOL, [CS] * NROW,
                                                    XLL, YLL)
    sim = flopy.mf6.MFSimulation(sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    m = flopy.mf6.ModflowGwf(sim, modelname='m')
    idm = np.ones((2, ncpl), dtype=int)
    idm[0, _ic(0, 5)] = 0                       # layer 1 pinched out here
    flopy.mf6.ModflowGwfdisv(m, nlay=2, ncpl=ncpl, nvert=len(verts),
                             vertices=verts, cell2d=cell2d, top=100.0,
                             botm=[80.0, 50.0], idomain=idm)
    flopy.mf6.ModflowGwfnpf(m, k=[[0.05] * ncpl, [0.5] * ncpl])
    flopy.mf6.ModflowGwfic(m, strt=97.0)
    # one land object per active top cell
    top_cells = [(0 if idm[0, c] else 1, c) for c in range(ncpl)]
    pk = [[n, cid, 1, -1, 0.25, 0.01 * (n + 1), 0.275, 0.285, 0.28, 3.5]
          for n, cid in enumerate(top_cells)]
    flopy.mf6.ModflowGwfuzf(m, nuzfcells=len(pk), packagedata=pk,
                            perioddata={0: [[n, 0.0002, 0.0, 1.5, 0.275, 0.0,
                                             0.0, 0.0] for n in range(len(pk))]})
    # the stream ends at reach 2, which hands its flow to the pond; the
    # pond spills into reach 3, which carries on (MVR, as the builder cuts
    # the reaches out of a pond's footprint)
    pdat = [[r, (0, _ic(3, 5 - r)), CS, 2.0, 0.01, 99.0 - 0.1 * r, 0.5, 0.1,
             0.035, 1 if r in (0, 2, 3, 5) else 2, 1.0, 0] for r in range(6)]
    conn = [[0, -1], [1, 0, -2], [2, 1], [3, -4], [4, 3, -5], [5, 4]]
    flopy.mf6.ModflowGwfsfr(m, pname='sfr', nreaches=6, packagedata=pdat,
                            connectiondata=conn, mover=True)
    flopy.mf6.ModflowGwflak(
        m, pname='lak', boundnames=True, nlakes=1, noutlets=1, ntables=1,
        mover=True,
        packagedata=[[0, 100.9, 1, 'pond7']],
        connectiondata=[[0, 0, (0, _ic(*POND)), 'embeddedv', 0.001, 0.0, 0.0,
                         0.0, 0.0]],
        outlets=[[0, 0, -1, 'manning', 101.0, 1.0, 0.03, 0.001]],
        tables=[[0, 'm.lak1.tab']])
    flopy.mf6.ModflowGwfmvr(m, maxmvr=2, maxpackages=2,
                            packages=[['sfr'], ['lak']],
                            perioddata={0: [['sfr', 2, 'lak', 0, 'FACTOR', 1.0],
                                            ['lak', 0, 'sfr', 3, 'FACTOR',
                                             1.0]]})
    sim.write_simulation(silent=True)
    with open(os.path.join(ws, 'm.lak1.tab'), 'w') as fh:
        fh.write('BEGIN dimensions\n  NROW  2\n  NCOL  3\nEND dimensions\n\n'
                 'BEGIN table\n  98.5  0.0  10.0\n  101.0  1000.0  400.0\n'
                 'END table\n')
    return ncpl


def _dataset(ds):
    os.makedirs(ds, exist_ok=True)
    ring = [[XLL, YLL], [XLL + NCOL * CS, YLL], [XLL + NCOL * CS,
            YLL + NROW * CS], [XLL, YLL + NROW * CS], [XLL, YLL]]
    fc = {'type': 'FeatureCollection', 'marmites': {'crs': 'EPSG:23029'},
          'features': [{'type': 'Feature', 'properties': {'id': 1},
                        'geometry': {'type': 'Polygon',
                                     'coordinates': [ring]}}]}
    with open(os.path.join(ds, 'inputBOUNDARY.geojson'), 'w') as fh:
        json.dump(fc, fh)
    # flopy numbers rows from the top: row 2 is the middle one
    cx = XLL + (POND[1] + 0.5) * CS
    cy = YLL + (NROW - POND[0] - 0.5) * CS
    pond = [[cx - 5, cy - 5], [cx + 5, cy - 5], [cx + 5, cy + 5],
            [cx - 5, cy + 5], [cx - 5, cy - 5]]
    fc['features'] = [{'type': 'Feature', 'properties': {'fid': 7},
                       'geometry': {'type': 'Polygon', 'coordinates': [pond]}}]
    with open(os.path.join(ds, 'inputPONDS.geojson'), 'w') as fh:
        json.dump(fc, fh)
    y = YLL + (NROW - 3 - 0.5) * CS
    with open(os.path.join(ds, 'inputSTREAM.csv'), 'w') as fh:
        fh.write('# test network\nseg_id,seq,x,y\n')
        for k, x in enumerate((XLL + NCOL * CS - 1.0, XLL + 1.0)):
            fh.write('0,%d,%r,%r\n' % (k, x, y))


def test_the_pond_reads_its_name_rim_bottom_and_stage(tmp_path):
    ws = str(tmp_path / 'ws')
    _sim(ws)
    gwf, mg = CF.load_model(ws, 'm')
    (lk,) = CF.lakes(gwf, mg, ws)
    assert lk['name'] == 'pond7'
    assert lk['cells'] == [_ic(*POND)]
    assert lk['strt'] == pytest.approx(100.9)
    assert lk['rim'] == pytest.approx(101.0)          # the outlet's invert
    assert lk['bottom'] == pytest.approx(98.5)        # the table's first stage
    assert lk['index'] == 0
    assert lk['area'] == pytest.approx(400.0)         # the table's at the rim


def test_the_uzf_land_objects_carry_their_layer_and_parameters(tmp_path):
    ws = str(tmp_path / 'ws')
    ncpl = _sim(ws)
    gwf, mg = CF.load_model(ws, 'm')
    u = CF.uzf_land(gwf, mg)
    assert u['layer'][_ic(0, 5)] == 2                 # layer 1 inactive there
    assert u['layer'][_ic(0, 0)] == 1
    assert u['vks'][0] == pytest.approx(0.01)
    assert np.allclose(u['extdp'], 1.5)
    assert np.isfinite(u['thts']).sum() == ncpl


def test_the_reaches_map_cell_to_streambed(tmp_path):
    ws = str(tmp_path / 'ws')
    _sim(ws)
    gwf, mg = CF.load_model(ws, 'm')
    s = CF.sfr_reaches(gwf, mg)
    assert sorted(s) == sorted(_ic(3, j) for j in range(NCOL))
    rno, rtp, rbth, _w = s[_ic(3, 5)]
    assert (rno, rtp, rbth) == (1, pytest.approx(99.0), pytest.approx(0.5))


def test_the_stream_hands_its_flow_to_the_pond_and_takes_its_spill(tmp_path):
    ws = str(tmp_path / 'ws')
    _sim(ws)
    gwf, _mg = CF.load_model(ws, 'm')
    down, up, inlets, outlets, into = CF.stream_links(gwf)
    assert inlets == {0: [2]} and outlets == {0: [3]} and into == {2: 0}
    assert down[0] == [1] and down[2] == [] and up[3] == []
    assert up[1] == [0] and down[3] == [4]


def test_a_pond_section_follows_the_stream_through_it(tmp_path):
    ws = str(tmp_path / 'ws')
    _sim(ws)
    gwf, _mg = CF.load_model(ws, 'm')
    links = CF.stream_links(gwf)
    rl = {r: CS for r in range(6)}
    # the whole stream: from the headwater to the pond, then to the outlet
    assert CF.stream_chain(0, links, length=1000.0, reach_len=rl) == (
        [0, 1, 2], [3, 4, 5])
    # one reach length each way: the reaches next to the pond
    assert CF.stream_chain(0, links, length=CS, reach_len=rl) == ([2], [3])


def test_the_main_stem_is_the_branch_with_most_reaches_above():
    # reach 4 takes reach 3 (a headwater) and reach 2 (below 0 and 1); the
    # pond's inlet is 4
    down = {0: [1], 1: [2], 2: [4], 3: [4], 4: []}
    up = {0: [], 1: [0], 2: [1], 3: [], 4: [2, 3]}
    links = (down, up, {0: [4]}, {}, {4: 0})
    ups, dns = CF.stream_chain(0, links, length=1e9)
    assert ups == [0, 1, 2, 4] and dns == []


def test_every_page_is_drawn(tmp_path):
    pytest.importorskip('matplotlib')
    ws, ds, out = (str(tmp_path / k) for k in ('ws', 'ds', 'out'))
    os.makedirs(out)
    ncpl = _sim(ws)
    _dataset(ds)
    import matplotlib.dates as mdates
    d0 = datetime.date(2000, 10, 1)
    days = [mdates.date2num(d0 + datetime.timedelta(days=k))
            for k in range(400)]
    cMF = types.SimpleNamespace(inputDate=days, nper=400,
                                elev=np.full(ncpl, 101.5))
    cover = np.zeros((2, ncpl))
    cover[0, ::2] = 80.0
    cover[1, 1::2] = 60.0
    ctx = types.SimpleNamespace(gridVEGarea=cover,
                                P_veg_zoneSP=np.full((1, 400), 1.5),
                                Eo_zonesSP=np.full((1, 400), 3.0))
    veg = [types.SimpleNamespace(name='grass', root_depth=0.5),
           types.SimpleNamespace(name='oak', root_depth=4.0)]
    cfg = types.SimpleNamespace(surface=types.SimpleNamespace(vegetation=veg))
    got = CF.cdl_figures(out, ws, 'm', ds, cMF=cMF, ctx=ctx, cfg=cfg,
                         verbose=False)
    names = sorted(os.path.basename(f) for f in got)
    assert names == ['IN_101_aquifer_params.png', 'IN_102_layer_surfaces.png',
                     'IN_103_ibound_uzf.png', 'IN_104_land_cover_extdp.png',
                     'IN_105_uzf_params.png', 'IN_106_forcing.png',
                     'IN_107_ponds_zoom_p1.png', 'IN_108_pond_sections_p1.png']
    assert all(os.path.getsize(f) > 0 for f in got)


def test_without_the_run_only_the_simulation_figures_are_drawn(tmp_path):
    pytest.importorskip('matplotlib')
    ws, out = str(tmp_path / 'ws'), str(tmp_path / 'out')
    os.makedirs(out)
    _sim(ws)
    got = CF.cdl_figures(out, ws, 'm', str(tmp_path / 'no_ds'), verbose=False)
    names = sorted(os.path.basename(f) for f in got)
    assert 'IN_104_land_cover_extdp.png' not in names    # needs ctx
    assert 'IN_108_pond_sections_p1.png' not in names    # needs cMF
    assert 'IN_101_aquifer_params.png' in names


# ------------------------------------------------------- output figures

def test_the_run_timing_is_drawn_from_the_run_record(tmp_path):
    pytest.importorskip('matplotlib')
    ws, out = str(tmp_path / 'ws'), str(tmp_path / 'out')
    _sim(ws)
    import matplotlib.dates as mdates
    d0 = datetime.date(2000, 10, 1)
    days = [mdates.date2num(d0 + datetime.timedelta(days=k))
            for k in range(30)]
    cMF = types.SimpleNamespace(inputDate=days)
    ctx = types.SimpleNamespace(P_veg_zoneSP=np.full((1, 30), 2.0))
    res = {'outer_iters': np.arange(30) % 7 + 5, 'iters_kind': b'linear',
           'sp_seconds': np.full(30, 0.5), 'ss_seconds': 3.0, 'ss_iters': 40}
    os.makedirs(out)
    written = []
    CO.fig_sp_timing(out, ws, 'm', str(tmp_path / 'ds'), cMF, ctx, res,
                     written)
    assert [os.path.basename(f) for f in written] == ['sp_run_timing.png']
    import pandas as pd
    t = pd.read_csv(os.path.join(out, 'sp_run_timing.csv'))
    assert len(t) == 30 and t['seconds'].sum() == pytest.approx(15.0)


def test_the_pond_recharge_page_and_its_numbers(tmp_path, monkeypatch):
    pytest.importorskip('matplotlib')
    ws, out = str(tmp_path / 'ws'), str(tmp_path / 'out')
    os.makedirs(out)
    _sim(ws)
    keys = ['from_stream', 'to_stream', 'outlet', 'runoff', 'rainfall',
            'evaporation', 'from_aquifer', 'to_aquifer', 'storage', 'other']
    rate = np.zeros((365, 2, len(keys)))
    rate[:, 0, keys.index('runoff')] = 10.0           # m3/d
    rate[:, 0, keys.index('from_aquifer')] = 2.0      # pond 7 gains
    rate[:, 1, keys.index('to_aquifer')] = -3.0       # pond 2 gives
    rate[:, :, keys.index('evaporation')] = -1.0
    perlen = np.ones(365)
    annual = np.einsum('p,pnc->nc', perlen, rate) / 365.0 * 365.25
    monkeypatch.setattr(CO, 'pond_volumes', lambda *a: (
        ['pond7', 'pond2'], keys, annual, rate, perlen))
    monkeypatch.setattr(CO, '_uzf_recharge', lambda *a: 100.0)
    written = []
    CO.fig_pond_recharge(out, ws, 'm', str(tmp_path / 'ds'), None, None,
                         written)
    assert [os.path.basename(f) for f in written] == [
        'pond_recharge_analysis.png']
    import pandas as pd
    t = pd.read_csv(os.path.join(out, 'pond_recharge_analysis.csv'),
                    index_col='pond')
    assert list(t.index) == ['pond2', 'pond7']        # in fid order
    # m3/yr, written to 0.1 m3
    assert t.loc['pond7', 'runoff'] == pytest.approx(3652.5, abs=0.05)
    assert t.loc['pond7', 'from_aquifer'] == pytest.approx(730.5, abs=0.05)
    assert t.loc['pond2', 'to_aquifer'] == pytest.approx(-1095.75, abs=0.05)


def test_no_steady_heads_page_without_a_steady_first_period(tmp_path):
    ws = str(tmp_path / 'ws')
    _sim(ws)
    with open(os.path.join(ws, 'm.sto'), 'w') as fh:
        fh.write('BEGIN PERIOD 1\n  TRANSIENT\nEND PERIOD\n')
    assert CO.steady_heads_at_points(ws, 'm', str(tmp_path)) is None
