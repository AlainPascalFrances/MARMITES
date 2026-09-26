# -*- coding: utf-8 -*-
"""The model map (IN_000_model_map.png): what it reads from the written
simulation, and that it is drawn.

A small DISV model in UTM coordinates with every feature the map shows --
an outlet drain, a GHB, an SFR network with a specified-inflow inlet and one
outlet, and a lake -- written by flopy, never run.
"""
import importlib.util
import json
import os
import sys

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


PP = _load('marmites_postprocess_m',
           os.path.join(TRUNK, 'ppMF6', 'marmites_postprocess.py'))
GRID = _load('marmites_grid_m', os.path.join(TRUNK, 'marmites_grid.py'))

XLL, YLL, CS, NROW, NCOL = 739300.0, 4553050.0, 50.0, 4, 5


def _ic(i, j):
    return i * NCOL + j


def _sim(ws, mvr=False):
    verts, cell2d, ncpl = GRID.disv_from_structured([CS] * NCOL, [CS] * NROW,
                                                    XLL, YLL)
    sim = flopy.mf6.MFSimulation(sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    m = flopy.mf6.ModflowGwf(sim, modelname='m')
    flopy.mf6.ModflowGwfdisv(m, nlay=1, ncpl=ncpl, nvert=len(verts),
                             vertices=verts, cell2d=cell2d, top=10.0,
                             botm=[0.0])
    flopy.mf6.ModflowGwfnpf(m)
    flopy.mf6.ModflowGwfic(m, strt=9.0)
    # the outlet drain on the west edge, a GHB on the east edge
    flopy.mf6.ModflowGwfdrn(m, pname='drn', stress_period_data={0: [
        [(0, _ic(3, 0)), 5.0, 1.0], [(0, _ic(2, 0)), 5.0, 1.0]]})
    # seepage drains everywhere else must NOT be drawn as outlet drains
    flopy.mf6.ModflowGwfdrn(m, pname='drn_seep', stress_period_data={0: [
        [(0, _ic(1, 2)), 9.0, 100.0]]})
    flopy.mf6.ModflowGwfghb(m, pname='ghb', stress_period_data={0: [
        [(0, _ic(1, 4)), 9.5, 1.0]]})
    # three reaches along row 3, east to west; reach 0 takes an inflow
    pdat = [[r, (0, _ic(3, 2 - r)), 50.0, 2.0, 0.01, 8.0 - 0.3 * r, 0.5, 0.1,
             0.035, 1 if r in (0, 2) else 2, 1.0, 0] for r in range(3)]
    conn = [[0, -1], [1, 0, -2], [2, 1]]
    flopy.mf6.ModflowGwfsfr(m, pname='sfr', nreaches=3, packagedata=pdat,
                            connectiondata=conn, mover=mvr,
                            perioddata={0: [[0, 'INFLOW', 25.0],
                                            [1, 'INFLOW', 0.0]]})
    flopy.mf6.ModflowGwflak(m, pname='lak', nlakes=1, noutlets=0, mover=mvr,
                            packagedata=[[0, 8.0, 1]],
                            connectiondata=[[0, 0, (0, _ic(0, 1)), 'vertical',
                                             0.001, 0.0, 0.0, 0.0, 0.0]])
    if mvr:
        # the last reach hands its flow to the pond: it ends there, but it
        # is not where the catchment drains
        flopy.mf6.ModflowGwfmvr(m, maxmvr=1, maxpackages=2,
                                packages=[['sfr'], ['lak']],
                                perioddata={0: [['sfr', 2, 'lak', 0,
                                                 'FACTOR', 1.0]]})
    sim.write_simulation(silent=True)
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
    cx, cy = XLL + 1.5 * CS, YLL + 3.5 * CS
    pond = [[cx - 5, cy - 5], [cx + 5, cy - 5], [cx + 5, cy + 5],
            [cx - 5, cy + 5], [cx - 5, cy - 5]]
    fc['features'] = [{'type': 'Feature', 'properties': {'id': 7},
                       'geometry': {'type': 'Polygon', 'coordinates': [pond]}}]
    with open(os.path.join(ds, 'inputPONDS.geojson'), 'w') as fh:
        json.dump(fc, fh)
    with open(os.path.join(ds, 'inputSTREAM.csv'), 'w') as fh:
        fh.write('# test network\nseg_id,seq,x,y\n')
        for k, x in enumerate((XLL + 2.9 * CS, XLL + 0.1 * CS)):
            fh.write('0,%d,%r,%r\n' % (k, x, YLL + 0.5 * CS))
    with open(os.path.join(ds, 'inputObs.txt'), 'w') as fh:
        fh.write('# Name X Y lay\nP1 %r %r 1\n#P2 %r %r 1\n'
                 % (XLL + 60, YLL + 110, XLL + 10, YLL + 10))


def test_the_features_come_from_the_written_packages(tmp_path):
    ws = str(tmp_path / 'ws')
    _sim(ws)
    f = PP.model_map_features(ws, 'm')
    assert f['mg'].grid_type == 'vertex'
    assert f['drn'] == sorted([_ic(3, 0), _ic(2, 0)])     # not drn_seep
    assert f['ghb'] == [_ic(1, 4)]
    assert f['sfr'] == sorted(_ic(3, j) for j in (0, 1, 2))
    assert f['nreaches'] == 3
    assert f['outlet'] == [_ic(3, 0)]                      # no downstream
    assert f['inlet'] == [_ic(3, 2)]                       # INFLOW > 0
    assert f['lak'] == [_ic(0, 1)]


def test_the_pond_footprint_comes_from_the_dataset_polygons(tmp_path):
    ws, ds = str(tmp_path / 'ws'), str(tmp_path / 'ds')
    _sim(ws)
    _dataset(ds)
    f = PP.model_map_features(ws, 'm', ds_ws=ds)
    assert f['lak_foot'] == [_ic(0, 1)]
    assert f['nlakes'] == 1


def test_a_reach_feeding_a_pond_is_not_an_outlet(tmp_path):
    ws = str(tmp_path / 'ws')
    _sim(ws, mvr=True)
    f = PP.model_map_features(ws, 'm')
    assert f['outlet'] == []


def test_no_simulation_no_map(tmp_path):
    assert PP.model_map_features(str(tmp_path), 'm') is None
    assert PP._fig_model_map(str(tmp_path), str(tmp_path), 'm',
                             str(tmp_path), verbose=False) == []


def test_the_map_is_drawn(tmp_path):
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('agg')
    ws, ds, out = (str(tmp_path / d) for d in ('ws', 'ds', 'out'))
    _sim(ws)
    _dataset(ds)
    os.makedirs(out)
    written = PP._fig_model_map(out, ws, 'm', ds, title='Test',
                                verbose=False)
    assert [os.path.basename(w) for w in written] == ['IN_000_model_map.png']
    assert os.path.getsize(written[0]) > 10000


def test_a_ring_centroid_survives_utm_coordinates():
    r = np.array([[739300.0, 4553050.0], [739302.0, 4553050.0],
                  [739302.0, 4553052.0], [739300.0, 4553052.0],
                  [739300.0, 4553050.0]])
    assert PP._ring_centroid(r) == pytest.approx((739301.0, 4553051.0),
                                                 abs=1e-9)
