# -*- coding: utf-8 -*-
"""WP3 -- the SFR network as SFRmaker builds it (Leaf et al. 2021).

One outlet where the network leaves the catchment; each reach as long as the
channel MAPPED in its cell, its slope over the distance to the next reach and
its bed below the land surface; streambed parameters per segment; the pieces
a burn leaves touching only at a corner joined back; and the outlet written
as continuous observations. Synthetic cases where the answer is known by
construction, then the La Mata build.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


S = _load('marmites_sfr_o', os.path.join(TRUNK, 'ppMF6', 'marmites_sfr.py'))
TOPO = _load('marmites_topology_o',
             os.path.join(TRUNK, 'ppMF6', 'marmites_topology.py'))
GRID = _load('marmites_grid_o', os.path.join(TRUNK, 'marmites_grid.py'))


def _mesh(nrow=4, ncol=5, cs=50.0):
    """A structured grid AS a mesh: cells (ic, 0), rook adjacency only."""
    verts, cell2d, ncpl = GRID.disv_from_structured([cs] * ncol, [cs] * nrow,
                                                    0.0, 0.0)
    gp = {'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl, 'nlay': 1}
    return gp, TOPO.MeshTopology(gp)


def _ic(i, j, ncol=5):
    return (i * ncol + j, 0)


def _corner_case():
    """Row 3, columns 0-2, falls west to the outlet (3,0); a second piece
    climbs north from (2,3), which touches (3,2) only at a corner."""
    gp, topo = _mesh()
    w = np.zeros((20, 1))
    dem = np.full((20, 1), 800.0)
    for j, z in ((0, 700.0), (1, 701.0), (2, 702.0)):
        w[_ic(3, j)] = 2.0
        dem[_ic(3, j)] = z
    for i, z in ((2, 703.0), (1, 704.0), (0, 705.0)):
        w[_ic(i, 3)] = 2.0
        dem[_ic(i, 3)] = z
    return gp, topo, w, dem


# ------------------------------ the outlet ------------------------------ #

def test_the_outlet_is_the_lowest_stream_cell_on_the_catchment_edge():
    dem = np.array([[5.0, 4.0, 9.0],
                    [3.0, 1.0, 9.0],
                    [9.0, 2.0, 9.0]])
    cells = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 1)]

    def edge(c):
        return c[0] in (0, 2) or c[1] in (0, 2)
    # (1,1) is lowest overall but interior; (2,1) is the lowest on the edge
    assert S.catchment_outlet(cells, dem, edge) == (2, 1)


def test_no_stream_cell_on_the_edge_gives_no_outlet():
    dem = np.zeros((3, 3))
    assert S.catchment_outlet([(1, 1)], dem, lambda c: False) is None


def test_a_network_given_one_outlet_has_one_exit():
    gp, topo, w, dem = _corner_case()
    net = S.stream_network(w, dem, outlets=[_ic(3, 0)], topology=topo,
                           coords=lambda c: topo.xy[c[0]])
    assert net.outlets == [_ic(3, 0)]
    assert [c for c in net.cells if net.recv[c] is None] == [_ic(3, 0)]


# ------------------------------ bridging -------------------------------- #

def test_a_piece_touching_only_at_a_corner_is_an_error_without_coords():
    gp, topo, w, dem = _corner_case()
    with pytest.raises(Exception):
        S.stream_network(w, dem, outlets=[_ic(3, 0)], topology=topo)


def test_a_piece_touching_only_at_a_corner_is_joined_to_the_nearest_reach():
    gp, topo, w, dem = _corner_case()
    net = S.stream_network(w, dem, outlets=[_ic(3, 0)], topology=topo,
                           coords=lambda c: topo.xy[c[0]])
    assert net.nreaches == 6
    assert net.bridges == [(_ic(2, 3), _ic(3, 2))]
    assert net.recv[_ic(2, 3)] == _ic(3, 2)
    # the rest of the piece routes over shared faces into its lowest cell
    assert net.recv[_ic(0, 3)] == _ic(1, 3)
    assert net.recv[_ic(1, 3)] == _ic(2, 3)


def test_a_bridge_prefers_a_reach_that_is_not_higher():
    """The nearest routed cell is uphill of the piece; a lower one a little
    farther away is the one the water would reach."""
    gp, topo = _mesh()
    w = np.zeros((20, 1))
    dem = np.full((20, 1), 800.0)
    # routed: (3,0) outlet, (3,1), (3,2) HIGH, and (2,0) (1,0) going north
    for c, z in ((_ic(3, 0), 700.0), (_ic(3, 1), 701.0), (_ic(3, 2), 720.0),
                 (_ic(2, 0), 701.5), (_ic(1, 0), 702.0)):
        w[c], dem[c] = 2.0, z
    # the piece: (2,3) alone, at 710 -- nearest routed is (3,2) at 720
    w[_ic(2, 3)], dem[_ic(2, 3)] = 2.0, 710.0
    net = S.stream_network(w, dem, outlets=[_ic(3, 0)], topology=topo,
                           coords=lambda c: topo.xy[c[0]])
    (low, target), = net.bridges
    assert low == _ic(2, 3)
    assert target != _ic(3, 2)
    assert float(dem[target]) <= 710.0


# ------------------------- geometry from the map ------------------------ #

def _built(**kw):
    gp, topo, w, dem = _corner_case()
    net = S.stream_network(w, dem, outlets=[_ic(3, 0)], topology=topo,
                           coords=lambda c: topo.xy[c[0]])
    args = dict(pondw=w, delr=np.ones(1), delc=np.ones(20), verbose=False,
                spacing=lambda c, r: topo.distance(c[0], r[0]))
    args.update(kw)
    S.build_sfr(net, dem, **args)
    return net, topo, w, dem


def test_reach_length_is_the_channel_mapped_in_the_cell():
    chl = np.zeros((20, 1))
    for k, c in enumerate(zip(*np.where(_corner_case()[2] > 0))):
        chl[c] = 10.0 + k
    net, *_ = _built(reach_length=chl)
    for k, c in enumerate(net.cells):
        assert net.reach_len[k] == pytest.approx(chl[c])


def test_without_the_map_a_mesh_reach_is_not_a_cell_number_difference():
    """The legacy rule on the (ncpl, 1) proxy grid: |ic - ic'| metres."""
    net, topo, w, dem = _built()
    for k, c in enumerate(net.cells):
        r = net.recv[c]
        if r is not None:
            assert net.reach_len[k] == pytest.approx(topo.distance(c[0], r[0]))


def test_slope_is_taken_over_the_distance_to_the_next_reach():
    chl = np.full((20, 1), 3.0)            # short mapped reaches...
    net, topo, w, dem = _built(reach_length=chl, monotonic=False)
    k = net.rno[_ic(3, 1)]
    # ...but 1 m of fall over the 50 m between centroids, not over 3 m
    assert net.reach_slope[k] == pytest.approx(1.0 / 50.0)


def test_streambed_parameters_can_differ_per_cell():
    gp, topo, w, dem = _corner_case()
    rhk = np.where(w > 0, 0.3, 0.0)
    rhk[_ic(1, 3)] = 0.05
    man = np.where(w > 0, 0.04, 0.0)
    rbth = np.where(w > 0, 0.5, 0.0)
    rbth[_ic(3, 1)] = 1.2
    net, *_ = _built(rhk=rhk, man=man, rbth=rbth)
    for row, c in zip(net.packagedata, net.cells):
        assert row[6] == pytest.approx(rbth[c])
        assert row[7] == pytest.approx(rhk[c])
        assert row[8] == pytest.approx(man[c])


def test_monotonic_bed_can_be_switched_off():
    gp, topo, w, dem = _corner_case()
    dem = dem.copy()
    dem[_ic(3, 1)] = 699.0                 # a DEM pit below the outlet
    net_on = S.stream_network(w, dem, outlets=[_ic(3, 0)], topology=topo,
                              coords=lambda c: topo.xy[c[0]])
    S.build_sfr(net_on, dem, pondw=w, verbose=False, monotonic=True,
                spacing=lambda c, r: topo.distance(c[0], r[0]))
    net_off = S.stream_network(w, dem, outlets=[_ic(3, 0)], topology=topo,
                               coords=lambda c: topo.xy[c[0]])
    S.build_sfr(net_off, dem, pondw=w, verbose=False, monotonic=False,
                spacing=lambda c, r: topo.distance(c[0], r[0]))
    assert net_on.nmono >= 1 and net_off.nmono == 0
    k = net_off.rno[_ic(3, 0)]
    assert net_on.reach_top[k] < net_off.reach_top[k]
    # every reach the DEM's own bed, minus the same incision
    top_off = np.array(net_off.reach_top)
    z = np.array([float(dem[c]) for c in net_off.cells])
    assert np.ptp(z - top_off) == pytest.approx(0.0)


# ------------------------------ the burn -------------------------------- #

def test_a_projected_mesh_burns_onto_its_mesh_not_its_proxy_grid():
    """The model carries a projected mesh as ``mesh_gridprops`` with
    delr = delc = 1 on an (ncpl, 1) proxy: burning onto the proxy put 97
    mapped streams on 4 cells of La Mata's mesh."""
    mv = _load('marmites_vector_o', os.path.join(TRUNK, 'marmites_vector.py'))
    gp, _topo = _mesh()

    class _CMF:
        mesh_gridprops = gp
        delr = np.ones(1)
        delc = np.ones(20)
        xllcorner = yllcorner = 0.0
    tg = mv.TargetGrid.from_cMF(_CMF())
    assert tg.ncell == 20
    x0, y0, x1, y1 = tg.extent
    assert (x1 - x0, y1 - y0) == pytest.approx((250.0, 200.0))


# ------------------ La Mata: DRN at the outlet, and obs ------------------ #

def test_lamata_outlet_obs_and_only_the_outlet_drain_record_goes(tmp_path):
    flopy = pytest.importorskip('flopy')
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s1L.ini')):
        pytest.skip('La Mata dataset not present')
    import matplotlib
    matplotlib.use('agg')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    import marmites_channel as mch
    import marmites_config as cfgmod
    import marmites_vector as mv
    mf6mod = _load('marmites_mf6_o', os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'))
    XLL, YLL = 739300.0, 4553050.0
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s1L.ini',
                   xllcorner=XLL, yllcorner=YLL)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]
    b = mf6mod.clsMF6(c, top=np.asarray(c.elev, float),
                      botm=np.asarray(c.botm, float), sim_ws=str(tmp_path),
                      daily=True)
    b.verbose = False
    lines, seg_params = mch.read_stream_lines(
        os.path.join(DS, 'inputSTREAM.csv'),
        os.path.join(DS, 'inputSTREAM_param.csv'))
    vgrid = mv.TargetGrid.structured(c.delr, c.delc, c.xllcorner, c.yllcorner)
    present, seg_of_cell, ch_len = mch.burn_channel(lines, vgrid,
                                                    (c.nrow, c.ncol))
    act = np.asarray(c.outcropL) > 0
    b.sfr_pondw = np.where(act, present, 0.0)
    b.sfr_pondhmax = np.zeros_like(b.sfr_pondw)
    b.sfr_seg_of_cell = seg_of_cell
    b.sfr_seg_params = seg_params
    b.sfr_cell_length = ch_len
    b.cell_area = float(np.mean(c.delr)) * float(np.mean(c.delc))
    b.sfr_width_source = cfgmod.ParamSource(
        drainage={'w_min': 1.5, 'w_max': 3.0, 'power': 2.0})
    b.sfr_depth_source = cfgmod.ParamSource(value=1.0)
    b.build()
    b.write()

    net = b.sfr_net
    assert len(net.outlets) == 1
    # the reach lengths are the mapped channel, cell by cell
    for k, cc in enumerate(net.cells):
        assert net.reach_len[k] == pytest.approx(max(ch_len[cc], S.SFR_MINLEN))
    # DRN: all the legacy records but the outlet reach's own
    legacy = {tuple(b._cellid(l, i, j))
              for (l, i, j, _e, _c) in c.layer_row_column_elevation_cond[0]}
    drn = b.gwf.get_package('drn')
    kept = ({tuple(r['cellid']) for r in drn.stress_period_data.get_data(0)}
            if drn is not None else set())
    assert legacy - kept == legacy & b.sfr_outlet_cellids
    assert len(legacy - kept) <= 1
    # the outlet as continuous observations, reach numbers from 1
    rno = net.rno[net.outlets[0]] + 1
    obs_files = [f for f in os.listdir(tmp_path) if f.endswith('.sfr.obs')]
    assert len(obs_files) == 1
    txt = open(os.path.join(tmp_path, obs_files[0])).read().lower()
    assert b.sfr_obs_csv in txt
    for kind in ('ext-outflow', 'stage', 'sfr'):
        assert any(ln.split()[1:] == [kind, str(rno)]
                   for ln in txt.splitlines() if len(ln.split()) == 3)


# ------------------------- excision (WP4.1, ponds) ---------------------- #

def _chain():
    """Row 3 of the 4 x 5 mesh falling west to the outlet (3, 0), plus a
    tributary joining at (3, 2) from (2, 2). Built, so the rows are real."""
    gp, topo = _mesh()
    w = np.zeros((20, 1))
    dem = np.full((20, 1), 900.0)
    for j in range(5):
        w[_ic(3, j)], dem[_ic(3, j)] = 2.0, 700.0 + j
    w[_ic(2, 2)], dem[_ic(2, 2)] = 2.0, 710.0
    net = S.stream_network(w, dem, outlets=[_ic(3, 0)], topology=topo,
                           coords=lambda c: topo.xy[c[0]])
    S.build_sfr(net, dem, pondw=w, delr=np.ones(1), delc=np.ones(20),
                spacing=lambda c, r: topo.distance(c[0], r[0]), verbose=False)
    return net


def _consistent(net):
    n = net.nreaches
    assert [int(r[0]) for r in net.packagedata] == list(range(n))
    assert [int(r[0]) for r in net.connectiondata] == list(range(n))
    downs = set()
    for row, prow in zip(net.connectiondata, net.packagedata):
        conns = [int(v) for v in row[1:]]
        assert prow[9] == len(conns)
        assert sum(1 for v in conns if v < 0) <= 1
        assert all(0 <= abs(v) < n for v in conns)
        downs.update(-v for v in conns if v < 0)
    assert 0 not in downs                  # reach 0 is a headwater: no -0
    for c in net.cells:
        r = net.recv[c]
        k = net.rno[c]
        conns = [int(v) for v in net.connectiondata[k][1:]]
        if r is None:
            assert not any(v < 0 for v in conns)
        else:
            assert -net.rno[r] in conns


def test_excising_a_pond_cuts_the_stream_and_hands_over_both_ends():
    net = _chain()
    pond = _ic(3, 2)
    into, out_of = S.excise_reaches(net, [pond])
    # the main stem above and the tributary drain into the pond; the reach
    # below it takes the spill
    assert sorted(into) == sorted([(_ic(3, 3), pond), (_ic(2, 2), pond)])
    assert out_of == [(pond, _ic(3, 1))]
    assert net.nreaches == 5 and pond not in net.rno
    assert net.recv[_ic(3, 3)] is None and net.recv[_ic(2, 2)] is None
    assert net.outlets == [_ic(3, 0)] and net.recv[_ic(3, 0)] is None
    assert net.excised == [pond]
    _consistent(net)


def test_an_excised_network_keeps_its_reach_attributes():
    net = _chain()
    before = {c: (net.reach_len[net.rno[c]], net.reach_top[net.rno[c]])
              for c in net.cells}
    S.excise_reaches(net, [_ic(3, 2), _ic(3, 3)])
    for c in net.cells:
        k = net.rno[c]
        assert (net.reach_len[k], net.reach_top[k]) == before[c]
        assert net.packagedata[k][5] == pytest.approx(before[c][1])


def test_a_pond_on_the_outlet_is_refused():
    net = _chain()
    with pytest.raises(ValueError, match='outlet'):
        S.excise_reaches(net, [_ic(3, 0)])


def test_excising_nothing_changes_nothing():
    net = _chain()
    rows = [list(r) for r in net.connectiondata]
    assert S.excise_reaches(net, [_ic(0, 0)]) == ([], [])  # not a stream cell
    assert [list(r) for r in net.connectiondata] == rows


def test_lamata_the_stream_is_cut_through_the_soil_into_the_aquifer(tmp_path):
    """The stream's total depth below the land surface is the soil depth of
    its cell + the channel depth + the streambed thickness (user rule,
    2026-09-27): the streambed top one channel depth below the AQUIFER top,
    its bottom never ON it. Measured from the land surface instead, La
    Mata's 1.5 m of soil = 1.0 m channel + 0.5 m streambed put 76 % of the
    streambed bottoms exactly on the aquifer top, where the water table sits,
    and MF6 took up to 500 outer iterations a period."""
    pytest.importorskip('flopy')
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s1L.ini')):
        pytest.skip('La Mata dataset not present')
    import matplotlib
    matplotlib.use('agg')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    import marmites_channel as mch
    import marmites_config as cfgmod
    import marmites_vector as mv
    mf6mod = _load('marmites_mf6_depth', os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'))
    XLL, YLL = 739300.0, 4553050.0
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s1L.ini',
                   xllcorner=XLL, yllcorner=YLL)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]
    land = np.asarray(np.ma.getdata(c.elev), float)
    soil, depth, rbth = 1.5, 1.0, 0.5
    b = mf6mod.clsMF6(c, top=land - soil, botm=np.asarray(c.botm, float),
                      sim_ws=str(tmp_path), daily=True)
    b.verbose = False
    lines, seg_params = mch.read_stream_lines(
        os.path.join(DS, 'inputSTREAM.csv'),
        os.path.join(DS, 'inputSTREAM_param.csv'))
    vgrid = mv.TargetGrid.structured(c.delr, c.delc, c.xllcorner, c.yllcorner)
    present, seg_of_cell, ch_len = mch.burn_channel(lines, vgrid,
                                                    (c.nrow, c.ncol))
    act = np.asarray(c.outcropL) > 0
    b.sfr_pondw = np.where(act, present, 0.0)
    b.sfr_pondhmax = np.zeros_like(b.sfr_pondw)
    b.sfr_seg_of_cell, b.sfr_seg_params = seg_of_cell, None
    b.sfr_cell_length = ch_len
    b.cell_area = 2500.0
    b.sfr_rbth = rbth
    b.sfr_width_source = cfgmod.ParamSource(
        drainage={'w_min': 1.5, 'w_max': 3.0, 'power': 2.0})
    b.sfr_depth_source = cfgmod.ParamSource(value=depth)
    b.build()
    net = b.sfr_net
    tops = np.array([land[cc] - soil for cc in net.cells])
    rtp = np.array(net.reach_top)
    # the bed top one channel depth below the aquifer top -- lower only
    # where the downstream-monotonic rule smoothed it
    assert np.all(rtp <= tops - depth + 1e-6)
    assert np.median(tops - depth - rtp) == pytest.approx(0.0, abs=1e-6)
    # ... so no streambed bottom sits on the aquifer top any more
    assert np.all((rtp - rbth) <= tops - depth - rbth + 1e-6)
    assert np.min(tops - (rtp - rbth)) >= depth + rbth - 1e-6
