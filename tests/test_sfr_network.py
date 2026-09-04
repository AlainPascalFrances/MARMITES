# -*- coding: utf-8 -*-
"""Routing the MARMITES channel map (PONDw) into an SFR reach network.

The routing is a priority flood grown upslope from the outlets, not greedy
steepest descent: the La Mata channel crosses flat stretches where every
neighbour has the same DEM value, and descent strands whole branches there.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..', 'trunk'))
DS = os.path.abspath(os.path.join(HERE, '..', 'DataSet_LaMata'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _load('marmites_sfr', os.path.join(TRUNK, 'ppMF6', 'marmites_sfr.py'))


def _asc(fn):
    a = np.loadtxt(fn, skiprows=6)
    return np.where(a <= -9999.0, 0.0, a)


def _lamata():
    fn = os.path.join(DS, 'inputPONDw.asc')
    if not os.path.exists(fn):
        pytest.skip('La Mata dataset not present')
    w = _asc(fn)
    dem = _asc(os.path.join(DS, 'MF_ws', 'elev_sinkfil.asc'))
    hm = _asc(os.path.join(DS, 'inputPONDhmax.asc'))
    dc = _asc(os.path.join(DS, 'MF_ws', 'drn_cond_l1.asc'))
    drn = [(int(i), int(j)) for i, j in zip(*np.where(dc > 0))]
    return w, dem, hm, drn


def _drains_to_outlet(net, c):
    seen = set()
    while c is not None:
        if c in net.outlets:
            return True
        if c in seen:
            return False          # cycle
        seen.add(c)
        c = net.recv[c]
    return False


# ------------------------------ routing ------------------------------- #

def test_simple_slope_routes_downhill():
    w = np.zeros((1, 5)); w[0, :] = 2.0
    dem = np.array([[104.0, 103.0, 102.0, 101.0, 100.0]])
    net = S.stream_network(w, dem, outlets=[(0, 4)])
    assert net.nreaches == 5
    for j in range(4):
        assert net.recv[(0, j)] == (0, j + 1)
    assert net.recv[(0, 4)] is None


def test_flat_stretch_is_routed_not_stranded():
    """Greedy steepest descent stops dead in a flat; the flood must not."""
    w = np.ones((1, 6))
    dem = np.array([[105.0, 104.0, 104.0, 104.0, 104.0, 100.0]])
    net = S.stream_network(w, dem, outlets=[(0, 5)])
    for j in range(5):
        assert _drains_to_outlet(net, (0, j)), 'cell (0,%d) is stranded' % j


def test_every_cell_reaches_an_outlet_and_no_cycles():
    w = np.zeros((5, 5))
    w[2, :] = 1.0
    w[:, 2] = 1.0                       # a cross: confluence in the middle
    dem = np.array([[10.0 + abs(i - 2) + abs(j - 2) for j in range(5)]
                    for i in range(5)])
    dem[2, 0] = 0.0
    net = S.stream_network(w, dem, outlets=[(2, 0)])
    assert net.nreaches == 9
    for c in net.cells:
        assert _drains_to_outlet(net, c)


def test_disconnected_component_is_reported_not_silently_dropped():
    w = np.zeros((1, 5))
    w[0, 0] = w[0, 1] = 1.0
    w[0, 3] = w[0, 4] = 1.0             # a separate island, no shared edge
    dem = np.array([[100.0, 101.0, 0.0, 102.0, 103.0]])
    with pytest.raises(ValueError, match='do not connect to any outlet'):
        S.stream_network(w, dem, outlets=[(0, 0)])


def test_outlet_must_be_a_stream_cell():
    w = np.zeros((1, 3)); w[0, 0] = w[0, 1] = 1.0
    dem = np.array([[100.0, 101.0, 102.0]])
    with pytest.raises(ValueError, match='not stream cells'):
        S.stream_network(w, dem, outlets=[(0, 2)])


def test_duplicate_outlets_are_collapsed():
    """A DRN cell list has one record per layer, so (i, j) repeats. Seeding
    the flood twice would duplicate reaches."""
    w = np.ones((1, 3))
    dem = np.array([[102.0, 101.0, 100.0]])
    net = S.stream_network(w, dem, drn_cells=[(0, 2), (0, 2), (0, 2)])
    assert net.nreaches == 3 and net.outlets == [(0, 2)]


def test_accumulation_counts_upstream_cells():
    w = np.ones((1, 4))
    dem = np.array([[103.0, 102.0, 101.0, 100.0]])
    net = S.stream_network(w, dem, outlets=[(0, 3)])
    assert net.acc[(0, 3)] == 4 and net.acc[(0, 0)] == 1


# --------------------------- reach attributes -------------------------- #

def _built(w, dem, **kw):
    net = S.stream_network(w, dem, outlets=kw.pop('outlets', None),
                           drn_cells=kw.pop('drn_cells', None))
    S.build_sfr(net, dem, delr=np.full(w.shape[1], 50.0),
                delc=np.full(w.shape[0], 50.0), verbose=False, **kw)
    return net


def test_reach_zero_is_a_headwater():
    """MF6 writes a downstream connection as -rno, and -0 loses its sign, so
    reach 0 must never be anyone's downstream target."""
    w = np.ones((1, 5))
    dem = np.array([[104.0, 103.0, 102.0, 101.0, 100.0]])
    net = _built(w, dem, outlets=[(0, 4)])
    targets = {-c for row in net.connectiondata for c in row[1:] if c < 0}
    assert 0 not in targets


def test_connection_count_matches_ncon_field():
    w = np.ones((1, 5))
    dem = np.array([[104.0, 103.0, 102.0, 101.0, 100.0]])
    net = _built(w, dem, outlets=[(0, 4)])
    for row in net.packagedata:
        assert row[9] == len(net.connectiondata[row[0]]) - 1


def test_downstream_connections_are_negative_upstream_positive():
    w = np.ones((1, 3))
    dem = np.array([[102.0, 101.0, 100.0]])
    net = _built(w, dem, outlets=[(0, 2)])
    mid = net.rno[(0, 1)]
    conns = net.connectiondata[mid][1:]
    assert net.rno[(0, 0)] in conns          # upstream, positive
    assert -net.rno[(0, 2)] in conns         # downstream, negative


def test_bed_is_downstream_monotonic():
    """DEM noise can leave a receiver higher than its own cell; SFR needs beds
    that never rise downstream (SFRmaker rule)."""
    w = np.ones((1, 4))
    dem = np.array([[103.0, 100.0, 102.0, 99.0]])   # (0,2) rises after (0,1)
    net = _built(w, dem, outlets=[(0, 3)])
    for c in net.cells:
        r = net.recv[c]
        if r is not None:
            assert net.reach_top[net.rno[r]] <= net.reach_top[net.rno[c]] + 1e-9
    assert net.nmono > 0, 'the rising bed should have been clamped'


def test_slope_never_below_the_floor():
    w = np.ones((1, 4))
    dem = np.array([[100.0, 100.0, 100.0, 100.0]])   # perfectly flat
    net = _built(w, dem, outlets=[(0, 3)], minslope=1e-4)
    assert min(net.reach_slope) >= 1e-4


def test_width_comes_from_the_channel_map():
    w = np.array([[1.5, 2.0, 3.0]])
    dem = np.array([[102.0, 101.0, 100.0]])
    net = _built(w, dem, outlets=[(0, 2)], pondw=w)
    assert net.reach_wid[net.rno[(0, 0)]] == pytest.approx(1.5)
    assert net.reach_wid[net.rno[(0, 2)]] == pytest.approx(3.0)


def test_bed_is_incised_below_land_surface_by_the_channel_depth():
    w = np.ones((1, 3))
    dem = np.array([[102.0, 101.0, 100.0]])
    hm = np.full((1, 3), 1.5)
    net = _built(w, dem, outlets=[(0, 2)], pondhmax=hm)
    assert net.reach_top[net.rno[(0, 0)]] == pytest.approx(100.5)


def test_diagonal_reach_is_longer_than_an_orthogonal_one():
    w = np.zeros((2, 2)); w[0, 0] = w[1, 1] = 1.0
    dem = np.array([[101.0, 0.0], [0.0, 100.0]])
    net = _built(w, dem, outlets=[(1, 1)])
    assert net.reach_len[net.rno[(0, 0)]] == pytest.approx(50.0 * np.sqrt(2.0))


def test_reach_descends_to_a_layer_deep_enough_for_the_streambed():
    w = np.ones((1, 2))
    dem = np.array([[101.0, 100.0]])
    # layer 0 is 0.2 m thick -- far too thin to contain a 0.5 m streambed
    botm = np.array([[[100.8, 99.8]], [[90.0, 89.0]]])
    net = S.stream_network(w, dem, outlets=[(0, 1)])
    S.build_sfr(net, dem, delr=np.full(2, 50.0), delc=np.full(1, 50.0),
                botm=botm, rbth=0.5, verbose=False)
    for row in net.packagedata:
        assert row[1][0] == 1, 'reach should have descended out of layer 0'


# ------------------------------ La Mata -------------------------------- #

def test_lamata_all_244_channel_cells_route_to_the_drn_outlets():
    w, dem, hm, drn = _lamata()
    net = S.stream_network(w, dem, drn_cells=drn)
    assert net.nreaches == 244
    assert len(net.outlets) == 6
    for c in net.cells:
        assert _drains_to_outlet(net, c)
    # one dominant branch: the catchment drains through a single trunk
    assert max(net.acc.values()) > 200


def test_lamata_reach_attributes_are_physical():
    w, dem, hm, drn = _lamata()
    net = S.stream_network(w, dem, drn_cells=drn)
    S.build_sfr(net, dem, pondhmax=hm, pondw=w, delr=np.full(60, 50.0),
                delc=np.full(65, 50.0), verbose=False)
    assert min(net.reach_wid) >= 1.5 and max(net.reach_wid) <= 3.0
    assert min(net.reach_len) > 0 and max(net.reach_len) <= 50.0 * np.sqrt(2.0) + 1e-6
    assert min(net.reach_slope) > 0
    assert 700.0 < min(net.reach_top) < max(net.reach_top) < 850.0
    assert len(net.packagedata) == len(net.connectiondata) == 244
