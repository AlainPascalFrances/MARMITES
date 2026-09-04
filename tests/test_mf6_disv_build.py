# -*- coding: utf-8 -*-
"""Phase-4 tests: clsMF6 builds a loadable DISV simulation from La Mata, and
the DISV model describes the same aquifer as the DIS model.

No mf6 binary required: build -> write -> flopy reload round-trip.
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

flopy = pytest.importorskip('flopy')
import matplotlib  # noqa: E402
matplotlib.use('agg')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mf6mod = _load('marmites_mf6_disv', os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'))
G = _load('marmites_grid_disv', os.path.join(TRUNK, 'marmites_grid.py'))


@pytest.fixture(scope='module')
def cmf():
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s3L.ini')):
        pytest.skip('La Mata dataset not present')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s3L.ini',
                   xllcorner=739300.0, yllcorner=4553050.0)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]
    return c


def _build(cmf, tmp, grid):
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp), grid=grid)
    b.build()
    return b


def test_disv_build_and_reload(cmf, tmp_path):
    b = _build(cmf, tmp_path, 'disv')
    assert b.grid == 'disv'
    assert b.ncpl == cmf.nrow * cmf.ncol == 3900
    assert len(b.vertices) == (cmf.nrow + 1) * (cmf.ncol + 1)
    b.write()
    files = os.listdir(tmp_path)
    name = cmf.modelname.lower()
    assert f'{name}.disv' in files
    assert f'{name}.dis' not in files
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(tmp_path), verbosity_level=0)
    gwf = sim.get_model(name)
    disv = gwf.get_package('disv')
    assert disv.nlay.get_data() == 6
    assert disv.ncpl.get_data() == 3900
    assert disv.nvert.get_data() == len(b.vertices)


def test_disv_cellids_use_icell2d(cmf, tmp_path):
    b = _build(cmf, tmp_path, 'disv')
    # UZF land cells: cellid must be (lay, icell2d) with icell2d = i*ncol+j
    for n, (i, j, k) in enumerate(b.surf_cells[:50]):
        cid = b.uzf_packagedata[n][1]
        assert len(cid) == 2, 'DISV cellid must be (lay, icell2d)'
        assert cid == (k, i * cmf.ncol + j)
    # WEL/DRN too
    wel = b.gwf.get_package('wel').stress_period_data.get_data(0)
    assert len(wel[0][0]) == 2
    drn = b.gwf.get_package('drn').stress_period_data.get_data(0)
    assert len(drn[0][0]) == 2


def test_disv_same_aquifer_as_dis(cmf, tmp_path):
    """The DISV model must describe the same aquifer: identical cell counts,
    identical per-cell top/botm/k when mapped through icell2d = i*ncol+j."""
    bd = _build(cmf, tmp_path / 'dis', 'dis')
    bv = _build(cmf, tmp_path / 'disv', 'disv')
    assert bd.ncell == bv.ncell
    assert bd.nuzfcells == bv.nuzfcells
    ncol = cmf.ncol
    # idomain
    id_d = bd.idomain
    id_v = bv._griddata(bd.idomain)
    for L in range(cmf.nlay):
        for (i, j, _k) in bd.surf_cells[:100]:
            assert id_d[L, i, j] == id_v[L, i * ncol + j]
    # NPF k / k33 and DISV top mapped consistently
    k_d = bd.gwf.get_package('npf').k.get_data()
    k_v = bv.gwf.get_package('npf').k.get_data()
    top_v = bv.gwf.get_package('disv').top.get_data()
    top_d = bd.gwf.get_package('dis').top.get_data()
    for (i, j, _k) in bd.surf_cells[:100]:
        ic = i * ncol + j
        assert np.isclose(k_d[0][i, j], np.ravel(k_v[0])[ic])
        assert np.isclose(top_d[i, j], np.ravel(top_v)[ic])


def test_disv_uzf_chaining_preserved(cmf, tmp_path):
    """Column chaining (land-first ordering + ivertcon) is grid-independent."""
    bd = _build(cmf, tmp_path / 'dis', 'dis')
    bv = _build(cmf, tmp_path / 'disv', 'disv')
    # same landflag / ivertcon structure, only cellids differ
    for n in range(0, bd.nuzfcells, 97):
        rd, rv = bd.uzf_packagedata[n], bv.uzf_packagedata[n]
        assert rd[0] == rv[0]            # iuzno
        assert rd[2] == rv[2]            # landflag
        assert rd[3] == rv[3]            # ivertcon
        assert np.isclose(rd[5], rv[5])  # vks
    assert all(r[2] == 1 for r in bv.uzf_packagedata[:bv.ncell])
    assert all(r[2] == 0 for r in bv.uzf_packagedata[bv.ncell:])


def test_disv_geometry_areas_match_dis(cmf, tmp_path):
    """Cell areas from the generated DISV vertices equal the DIS areas."""
    bv = _build(cmf, tmp_path, 'disv')
    cells = [(n, i, j, i * cmf.ncol + j) for n, (i, j, _k) in enumerate(bv.surf_cells)]
    g_dis = G.geometry_for(cmf, cells, grid='dis')
    g_disv = G.VertexGeometry.from_vertices(
        bv.vertices, bv.cell2d, np.array([c[3] for c in cells], dtype=int), nlay=cmf.nlay)
    assert np.allclose(g_dis.area, g_disv.area)
    assert np.allclose(g_disv.area, 50.0 * 50.0)   # La Mata 50 m cells
