# -*- coding: utf-8 -*-
"""The CdL input figures (marmites_cdl_figs): what they read from the
written simulation -- the ponds' rim, bottom and stage, the UZF land objects,
the reaches -- the section geometry, and that every page is drawn.

A small two-layer DISV model in UTM coordinates with UZF, SFR and one pond
(LAK, its table and outlet), written by flopy, never run.
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


CF = _load('marmites_cdl_figs_m',
           os.path.join(TRUNK, 'ppMF6', 'marmites_cdl_figs.py'))
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
    pdat = [[r, (0, _ic(3, 5 - r)), CS, 2.0, 0.01, 99.0 - 0.1 * r, 0.5, 0.1,
             0.035, 1 if r in (0, 5) else 2, 1.0, 0] for r in range(6)]
    conn = [[0, -1]] + [[r, r - 1, -(r + 1)] for r in range(1, 5)] + [[5, 4]]
    flopy.mf6.ModflowGwfsfr(m, pname='sfr', nreaches=6, packagedata=pdat,
                            connectiondata=conn)
    flopy.mf6.ModflowGwflak(
        m, pname='lak', boundnames=True, nlakes=1, noutlets=1, ntables=1,
        packagedata=[[0, 100.9, 1, 'pond7']],
        connectiondata=[[0, 0, (0, _ic(*POND)), 'embeddedv', 0.001, 0.0, 0.0,
                         0.0, 0.0]],
        outlets=[[0, 0, -1, 'manning', 101.0, 1.0, 0.03, 0.001]],
        tables=[[0, 'm.lak1.tab']])
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


def test_a_section_runs_across_the_nearest_stream():
    lines = [np.array([[0.0, 0.0], [100.0, 0.0]]),
             np.array([[0.0, 50.0], [0.0, 150.0]])]
    q, t = CF._nearest_on_lines(np.array([30.0, 10.0]), lines)
    assert np.allclose(q, [30.0, 0.0])
    assert abs(t[0]) == pytest.approx(1.0)            # along x: cut along y
    q, t = CF._nearest_on_lines(np.array([5.0, 120.0]), lines)
    assert np.allclose(q, [0.0, 120.0])
    assert abs(t[1]) == pytest.approx(1.0)


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
