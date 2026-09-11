# -*- coding: utf-8 -*-
"""LAK for the La Mata charcas: sub-grid ponds as EMBEDDEDV lakes.

Every pond polygon in lm_ponds.shp is SMALLER than one 50 m cell (341-2036 m2
against 2500 m2) and two of them contain no cell centre at all, so a lake
cannot be built by excavating cells. Each pond therefore becomes a single
EMBEDDEDV lake inside one host cell, with its true area carried by a
stage-volume-area table instead of being inferred from the cell.
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

pytest.importorskip('shapefile')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


LK = _load('marmites_lak', os.path.join(TRUNK, 'ppMF6', 'marmites_lak.py'))

def _pond_shapefile():
    """Where the pond cartography actually lives.

    NOT in the repository: WP1 removed `example/LaMata/GIS/` because the repo
    holds only what MM and MF read directly. The canonical location is
    `$MM_DATA_ROOT/GIS` (mm_paths.GIS); the old in-dataset path is still tried
    so an unmigrated checkout keeps working.
    """
    try:
        import mm_paths
        cand = [os.path.join(str(mm_paths.GIS), 'lm_ponds.shp')]
    except Exception:
        cand = []
    cand.append(os.path.join(DS, 'GIS', 'lm_ponds.shp'))
    for c in cand:
        if os.path.exists(c):
            return c
    return cand[0]


SHP = _pond_shapefile()
XLL, YLL, CS, NROW, NCOL = 739300.0, 4553050.0, 50.0, 65, 60


def _grid_kwargs():
    return dict(xll=XLL, yll=YLL, delr=np.full(NCOL, CS), delc=np.full(NROW, CS),
                nrow=NROW, ncol=NCOL)


def _ponds():
    if not os.path.exists(SHP):
        pytest.skip('pond shapefile not present')
    return LK.read_pond_polygons(SHP)


# ------------------------------ polygons ------------------------------- #

def test_read_returns_twelve_ponds_with_ids():
    ponds = _ponds()
    assert len(ponds) == 12
    assert sorted(p.fid for p in ponds) == [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13]


def test_every_pond_is_smaller_than_one_cell():
    """The premise of the whole EMBEDDEDV design -- assert it, so that a future
    grid or shapefile change that breaks it is caught here."""
    for p in _ponds():
        assert 0 < p.area < CS * CS, 'pond %s is %.0f m2, no longer sub-grid' % (p.fid, p.area)


def test_polygon_area_matches_the_shoelace_area():
    p = _ponds()[0]
    pts = np.asarray(p.points)
    expect = 0.5 * abs(np.dot(pts[:, 0], np.roll(pts[:, 1], 1))
                       - np.dot(pts[:, 1], np.roll(pts[:, 0], 1)))
    # shapefile coordinates are ~7e5 m in float32, so the two summation orders
    # differ in the last bits; 1e-5 relative is far tighter than any real use
    assert p.area == pytest.approx(expect, rel=1e-5)


def test_centroid_of_a_square():
    sq = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]
    area, cen = LK._ring_area_centroid(sq)
    assert area == pytest.approx(100.0)
    assert cen == pytest.approx((5.0, 5.0))


# ---------------------------- cell assignment -------------------------- #

def test_every_pond_gets_a_distinct_host_cell():
    ponds = LK.assign_pond_cells(_ponds(), verbose=False, **_grid_kwargs())
    cells = [p.cell for p in ponds]
    assert all(c is not None for c in cells)
    assert len(set(cells)) == len(cells), 'EMBEDDEDV forbids two lakes in one cell'


def test_ponds_without_a_cell_centre_are_still_assigned():
    """Ponds 9 and 13 contain no cell centre; a centroid-based rule still
    places them, a 'centres inside the polygon' rule would drop them."""
    ponds = LK.assign_pond_cells(_ponds(), verbose=False, **_grid_kwargs())
    for fid in (9, 13):
        p = next(q for q in ponds if q.fid == fid)
        assert p.cell is not None


def test_host_cell_actually_contains_the_centroid():
    ponds = LK.assign_pond_cells(_ponds(), verbose=False, **_grid_kwargs())
    for p in ponds:
        i, j = p.cell
        x0, x1 = XLL + j * CS, XLL + (j + 1) * CS
        y1 = YLL + (NROW - i) * CS
        y0 = y1 - CS
        cx, cy = p.centroid
        assert x0 <= cx <= x1 and y0 <= cy <= y1, 'pond %s misplaced' % p.fid


def test_shared_host_cell_is_rejected():
    a = LK.PondLake(1, [(0, 0)], 100.0, (XLL + 10.0, YLL + 10.0))
    b = LK.PondLake(2, [(0, 0)], 100.0, (XLL + 20.0, YLL + 20.0))  # same cell
    with pytest.raises(ValueError, match='share a host cell'):
        LK.assign_pond_cells([a, b], verbose=False, **_grid_kwargs())


def test_inactive_host_falls_back_to_the_nearest_active_cell():
    idom = np.zeros((NROW, NCOL), dtype=int)
    idom[:, 30:] = 1
    p = LK.PondLake(1, [(0, 0)], 100.0, (XLL + 25.0, YLL + NROW * CS - 25.0))
    LK.assign_pond_cells([p], idomain=idom, verbose=False, **_grid_kwargs())
    assert idom[p.cell] > 0


# ------------------------------ lake table ----------------------------- #

def test_table_spans_bottom_to_above_the_rim():
    rows = LK.lake_table(100.0, 102.0, 1000.0)
    assert rows[0][0] == pytest.approx(100.0)
    assert rows[-1][0] > 102.0, 'the table must have headroom above the rim'
    stages = [r[0] for r in rows]
    assert stages == sorted(stages) and len(set(stages)) == len(stages)


def test_area_reaches_the_true_polygon_area_at_the_rim():
    rows = LK.lake_table(100.0, 102.0, 1234.0)
    at_rim = [r for r in rows if r[0] == pytest.approx(102.0)][0]
    assert at_rim[2] == pytest.approx(1234.0)


def test_area_tapers_toward_the_bed_so_the_lake_dries_smoothly():
    """A flat bottom steps the wetted area 0 -> full and a bone-dry lake
    degenerates numerically; the wedge bathymetry avoids that."""
    rows = LK.lake_table(100.0, 102.0, 1000.0)
    areas = [r[2] for r in rows]
    assert 0 < areas[0] < 0.05 * 1000.0
    assert areas == sorted(areas)


def test_volume_is_the_integral_of_area_over_stage():
    rows = LK.lake_table(100.0, 102.0, 1000.0)
    assert rows[0][1] == pytest.approx(0.0)
    for (s0, v0, a0, _), (s1, v1, a1, _) in zip(rows, rows[1:]):
        assert v1 - v0 == pytest.approx(0.5 * (a0 + a1) * (s1 - s0), abs=1e-3)


def test_bed_area_equals_surface_area():
    """EMBEDDEDV needs 4 columns; the wetted bed is the exchange area."""
    for r in LK.lake_table(100.0, 102.0, 900.0):
        assert len(r) == 4 and r[3] == pytest.approx(r[2])


# --------------------------- model integration ------------------------- #

def test_lamata_lak_mvr_model_builds_and_reloads(tmp_path):
    flopy = pytest.importorskip('flopy')
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s3L.ini')):
        pytest.skip('La Mata dataset not present')
    if not os.path.exists(SHP):
        pytest.skip('pond shapefile not present (%s)' % SHP)
    import matplotlib
    matplotlib.use('agg')
    from test_mf6_build import cmf as _cmf  # noqa: F401
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    mf6mod = _load('marmites_mf6', os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'))

    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s3L.ini',
                   xllcorner=XLL, yllcorner=YLL)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]

    def asc(fn):
        a = np.loadtxt(os.path.join(DS, fn), skiprows=6)
        return np.where(a <= -9999.0, 0.0, a)

    b = mf6mod.clsMF6(c, top=np.asarray(c.elev, float), botm=np.asarray(c.botm, float),
                      sim_ws=str(tmp_path), daily=True)
    b.verbose = False
    # WP1d: the network is the MAPPED hydrography burned onto the grid, not
    # inputSTREAMw.asc -- which was the alluvium footprint, and is retired.
    import marmites_channel as mch
    import marmites_config as cfgmod
    import marmites_vector as mv
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
    b.lak_shapefile = SHP
    b.lak_depth = np.full((c.nrow, c.ncol), 1.0)
    b.build()
    b.write()

    assert len(b.ponds) == 12
    assert sum(1 for p in b.ponds if p.on_channel) == 11
    # every lake has exactly one connection, as EMBEDDEDV requires
    lak = b.gwf.get_package('lak')
    conn = lak.connectiondata.get_data()
    assert len(conn) == 12
    for row in lak.packagedata.get_data():
        assert row['nlakeconn'] == 1
    # the stream is routed THROUGH the on-channel ponds
    mvr = b.gwf.get_package('mvr').perioddata.get_data(0)
    into = [r for r in mvr if str(r['pname1']).lower() == 'sfr']
    outof = [r for r in mvr if str(r['pname1']).lower() == 'lak']
    assert len(into) == 11 and len(outof) == 11
    # a pond never spills into the very reach that feeds it
    for p in b.ponds:
        if p.inlet_reach is not None:
            assert p.inlet_reach != p.outlet_reach

    sim = flopy.mf6.MFSimulation.load(sim_ws=str(tmp_path), verbosity_level=0)
    names = [pk.package_name for pk in sim.get_model().packagelist]
    assert 'lak' in names and 'mvr' in names and 'sfr' in names


def test_initial_stage_stays_between_bed_and_rim(tmp_path):
    """A lake started above its rim or below its bed blows up on step one."""
    pytest.importorskip('flopy')
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s3L.ini')):
        pytest.skip('La Mata dataset not present')
    if not os.path.exists(SHP):
        pytest.skip('pond shapefile not present (%s)' % SHP)
    import matplotlib
    matplotlib.use('agg')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    mf6mod = _load('marmites_mf6', os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'))
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s3L.ini',
                   xllcorner=XLL, yllcorner=YLL)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]
    b = mf6mod.clsMF6(c, top=np.asarray(c.elev, float), botm=np.asarray(c.botm, float),
                      sim_ws=str(tmp_path), daily=True)
    b.verbose = False
    b.lak_shapefile = SHP
    b.build()
    for p in b.ponds:
        assert p.bottom < p.strt <= p.rim, 'pond %s starts outside its own basin' % p.fid
