# -*- coding: utf-8 -*-
"""grid.kind = 'imported' -- an existing MODFLOW 6 DISV grid taken as it is,
cells AND layers (2026-11-12, Casa de Lobos).

The fixture is a real binary grid file, written by MODFLOW 6 at test time
from a 6-cell, 3-layer toy model (_write_toy_grb; the repository holds no
.grb, by its own rule). Per cell, 1 = the layer exists, p = the file marks
it as a pass-through:

    cell   1    2    3    4    5    6
           111  1p1  p11  pp1  111  111
"""

import importlib.util
import os
import shutil
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6'), HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)



def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CFG = _load('marmites_config_imp', os.path.join(CODE, 'marmites_config.py'))
MESHES = _load('marmites_meshes_imp', os.path.join(CODE, 'marmites_meshes.py'))
MESH = _load('marmites_mesh_imp', os.path.join(CODE, 'marmites_mesh.py'))

pytest.importorskip('flopy')


def _write_toy_grb(exe, ws):
    """Run the toy model once; returns the path of its toy.disv.grb.

    The cells sit on the 4 x 5 grid of 50 m cells at (1000, 2000): cells
    1-5 are the 50 m columns over y 2000-2150, cell 6 the strip y 2150-2200
    above them. top 10 m; where a layer exists it is 2 m (layer 1), 3 m
    (layer 2) and down to -5 m (layer 3); an absent layer is written 0.1 m
    thick and as a pass-through, as CdL's grid does.
    """
    import flopy
    xs = [1000.0 + 50.0 * i for i in range(6)]
    verts, vid = [], {}
    for y in (2000.0, 2150.0, 2200.0):
        for x in xs:
            vid[(x, y)] = len(verts)
            verts.append([len(verts), x, y])
    cell2d = []
    for c in range(5):                       # clockwise from the upper left
        ring = [vid[(xs[c], 2150.0)], vid[(xs[c + 1], 2150.0)],
                vid[(xs[c + 1], 2000.0)], vid[(xs[c], 2000.0)]]
        cell2d.append([c, xs[c] + 25.0, 2075.0, 4] + ring)
    ring = [vid[(xs[0], 2200.0)], vid[(xs[5], 2200.0)]]         + [vid[(x, 2150.0)] for x in reversed(xs)]
    cell2d.append([5, 1125.0, 2175.0, len(ring)] + ring)
    idomain = np.where(EXISTS > 0, 1, -1)
    top = np.full(6, 10.0)
    botm = np.zeros((3, 6))
    upper = top.copy()
    for k, t in enumerate((2.0, 3.0, None)):
        botm[k] = (np.where(idomain[k] > 0, upper - t, upper - 0.1)
                   if t is not None else np.full(6, -5.0))
        upper = botm[k]
    sim = flopy.mf6.MFSimulation(sim_name='toy', sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname='toy')
    flopy.mf6.ModflowGwfdisv(gwf, nlay=3, ncpl=6, nvert=len(verts), top=top,
                             botm=botm, idomain=idomain, vertices=verts,
                             cell2d=cell2d)
    flopy.mf6.ModflowGwfic(gwf, strt=9.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(2, 0), 9.0]])
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    assert ok, 'MODFLOW 6 did not run the toy model'
    return os.path.join(ws, 'toy.disv.grb')


@pytest.fixture(scope='module')
def GRB(tmp_path_factory):
    import mm_paths
    exe = str(getattr(mm_paths, 'MF6_EXE', '') or '')
    if not os.path.exists(exe):
        pytest.skip('mf6.exe not found (%s): the grid-file fixture is '
                    'written by MODFLOW 6' % exe)
    return _write_toy_grb(exe, str(tmp_path_factory.mktemp('toy_grb')))

EXISTS = np.array([[1, 1, 0, 0, 1, 1],
                   [1, 0, 1, 0, 1, 1],
                   [1, 1, 1, 1, 1, 1]])
NROW, NCOL, NLAY = 4, 5, 3
XLL, YLL, CS = 1000.0, 2000.0, 50.0
HNOFLO = 9999.999


def _cfg(file='toy.disv.grb'):
    """A valid 3-layer configuration on an imported grid."""
    c = CFG.RunConfig.from_dict(
        {'grid': {'kind': 'imported', 'imported': {'file': 'toy.disv.grb'}},
         'layers': {'nlay': NLAY, 'k33': [2.0, 5.0, 10.0],
                    'sr': [0.2, 0.2, 0.1]}})
    c.grid.imported.file = file
    return c


def _about(c, key):
    return [m for m in c.problems() if key in m]


# ------------------------------------------------------------ configuration

def test_imported_is_a_grid_kind_and_a_mesh_kind():
    assert 'imported' in CFG.GRID_KINDS
    assert 'imported' in CFG.MESH_KINDS
    assert 'structured' not in CFG.MESH_KINDS and 'disv' not in CFG.MESH_KINDS


def test_an_imported_grid_needs_its_file():
    assert not _about(_cfg(), 'grid.imported.file')
    assert _about(_cfg(file=''), 'grid.imported.file')


def test_the_file_is_the_binary_grid_file_of_a_disv_model():
    got = _about(_cfg(file='cdl_gwf.disv'), 'grid.imported.file')
    assert got and '.disv.grb' in got[0]


def test_the_file_is_not_asked_of_another_kind():
    c = _cfg(file='')
    c.grid.kind = 'voronoi'
    assert not _about(c, 'grid.imported')


def test_ibound_is_not_required_with_an_imported_grid():
    """It is required of every other kind; here the grid file holds it."""
    c = _cfg()
    c.layers.ibound = CFG.VectorSource()
    assert not _about(c, 'layers.ibound')
    c.grid.kind = 'voronoi'
    assert _about(c, 'layers.ibound')


def test_the_block_survives_a_round_trip(tmp_path):
    c = _cfg(file='MF_ws/cdl_gwf.disv.grb')
    fn = c.write_toml(str(tmp_path / 'c.toml'))
    back = CFG.load_run_config(fn)
    assert back.grid_kind == 'imported'
    assert back.grid.imported.file == 'MF_ws/cdl_gwf.disv.grb'


# ------------------------------------------------------------ the grid file

def test_the_cells_come_out_as_the_file_has_them(GRB):
    gp = MESHES.read_grid_file(GRB)
    assert gp['ncpl'] == 6 and gp['nlay'] == 3
    assert len(gp['vertices']) == 18
    assert [r[0] for r in gp['cell2d']] == list(range(6))
    assert [r[3] for r in gp['cell2d']] == [4, 4, 4, 4, 4, 8]
    # open rings, as MODFLOW 6 wants them
    assert all(r[4] != r[-1] for r in gp['cell2d'])
    assert gp['cell2d'][0][1:3] == [1025.0, 2075.0]
    rep = MESHES.cell_size_report(gp)
    assert rep['area_min'] == pytest.approx(50.0 * 150.0)
    assert rep['area_max'] == pytest.approx(250.0 * 50.0)


def test_a_pass_through_in_the_file_is_an_absent_layer(GRB):
    gp = MESHES.read_grid_file(GRB)
    assert (np.asarray(gp['idomain']) == EXISTS).all()


def test_the_thickness_is_the_files_and_zero_where_absent(GRB):
    """The file writes an absent layer 0.1 m thick; it has NO thickness
    here, so the layers below are not pushed down by it."""
    th = np.asarray(MESHES.read_grid_file(GRB)['thick'])
    assert th[0].tolist() == pytest.approx([2, 2, 0, 0, 2, 2])
    assert th[1].tolist() == pytest.approx([3, 0, 3, 0, 3, 3])
    # layer 3 reaches -5 m from wherever the stack above it ends
    assert th[2].tolist() == pytest.approx([10, 12.9, 11.9, 14.8, 10, 10])


def test_a_missing_or_foreign_file_is_refused(tmp_path):
    with pytest.raises(MESHES.MeshBuildError):
        MESHES.read_grid_file(str(tmp_path / 'none.disv.grb'))
    bad = tmp_path / 'bad.disv.grb'
    bad.write_bytes(b'not a grid file at all' * 20)
    with pytest.raises(MESHES.MeshBuildError):
        MESHES.read_grid_file(str(bad))


# ------------------------------------------------------------ the producer

class _Stub:
    nlay = NLAY
    nrow, ncol = NROW, NCOL
    delr, delc = [CS] * NCOL, [CS] * NROW
    xllcorner, yllcorner = XLL, YLL


def test_build_mesh_reads_the_file_from_the_dataset(GRB, tmp_path):
    ds = tmp_path / 'ds'
    (ds / 'MF_ws').mkdir(parents=True)
    shutil.copy2(GRB, str(ds / 'MF_ws' / 'toy.disv.grb'))
    c = _cfg(file='MF_ws/toy.disv.grb')
    notes = []
    gp, info = MESHES.build_mesh(c, _Stub(), cache_dir=str(tmp_path / 'c'),
                                 dataset_dir=str(ds), warn=notes.append)
    assert info['kind'] == 'imported' and info['ncpl'] == 6
    assert not info['cached']
    assert any('2 / 4 / 6' not in n and '4 / 4 / 6' in n for n in notes)
    # the cache gives back the SAME thing, layers included
    gp2, info2 = MESHES.build_mesh(c, _Stub(), cache_dir=str(tmp_path / 'c'),
                                   dataset_dir=str(ds))
    assert info2['cached'] and info2['signature'] == info['signature']
    assert (np.asarray(gp2['idomain']) == np.asarray(gp['idomain'])).all()
    assert np.allclose(gp2['thick'], gp['thick'])
    assert gp2['cell2d'] == gp['cell2d']


def test_the_layer_count_must_be_the_files(GRB, tmp_path):
    shutil.copy2(GRB, str(tmp_path / 'toy.disv.grb'))
    c = _cfg(file='toy.disv.grb')
    stub = _Stub()
    stub.nlay = 2
    with pytest.raises(MESHES.MeshBuildError) as e:
        MESHES.build_mesh(c, stub, dataset_dir=str(tmp_path))
    assert 'layers.nlay' in str(e.value)


def test_the_signature_is_the_files_content(GRB, tmp_path):
    """Not its name, and not the rectangle the rasters are on."""
    fn = tmp_path / 'toy.disv.grb'
    shutil.copy2(GRB, str(fn))
    c = _cfg(file='toy.disv.grb')
    s1 = MESHES.mesh_signature(c, _Stub(), str(tmp_path))
    other = _Stub()
    other.ncol, other.delr = 7, [CS] * 7
    assert MESHES.mesh_signature(c, other, str(tmp_path)) == s1
    with open(str(fn), 'ab') as fh:
        fh.write(b'\0')
    assert MESHES.mesh_signature(c, _Stub(), str(tmp_path)) != s1


# ------------------------------------------------------------ the projection

class _Process:
    def __init__(self):
        self.nrow, self.ncol, self.nlay = NROW, NCOL, NLAY
        self.hnoflo = HNOFLO


class _FakeMF:
    """The model on the dataset grid, as props.placeholder_geometry leaves
    it: every layer active where the catchment is (here: not the north
    row), 1 m thick."""

    def __init__(self):
        self.nrow, self.ncol, self.nlay = NROW, NCOL, NLAY
        self.delr, self.delc = [CS] * NCOL, [CS] * NROW
        self.xllcorner, self.yllcorner = XLL, YLL
        self.hnoflo = HNOFLO
        self.ibound = np.ones((NLAY, NROW, NCOL), dtype=int)
        self.ibound[:, 0, :] = 0
        self.outcropL = (self.ibound[0] != 0).astype(int)
        self.iuzfbnd = self.outcropL.copy()
        elev = np.full((NROW, NCOL), 20.0)
        self.elev = np.ma.masked_array(elev, mask=(self.outcropL == 0))
        self.top = self.elev - 1.0                    # a 1 m soil column
        self.thick = np.ones((NLAY, NROW, NCOL))
        self.botm = np.stack([np.asarray(self.top) - (k + 1.0)
                              for k in range(NLAY)])
        self.strt = np.stack([elev - 2.0 for _ in range(NLAY)])
        self.hk_actual = [1.0, 2.0, 3.0]
        self.cPROCESS = _Process()


def _projected(GRB):
    notes = []
    gp = MESHES.normalise(MESHES.read_grid_file(GRB))
    m, _g, _p = MESH.project_model(_FakeMF(), gp, {}, warn=notes.append)
    return m, notes


def test_which_layers_exist_is_the_files_per_cell(GRB):
    m, _notes = _projected(GRB)
    ib = np.asarray(m.ibound)[:, :, 0]
    assert (ib[:, :5] == EXISTS[:, :5]).all()
    assert np.asarray(m.outcropL)[:5, 0].tolist() == [1, 1, 2, 3, 1]


def test_a_cell_the_model_does_not_reach_stays_out_and_is_said(GRB):
    """Cell 6 lies on the row the dataset grid has inactive: no land
    surface, no soil, no forcing there."""
    m, notes = _projected(GRB)
    assert (np.asarray(m.ibound)[:, 5, 0] == 0).all()
    assert any('1 cell(s) of the imported grid' in n for n in notes)


def test_the_bottoms_hang_from_the_models_own_top(GRB):
    """top = land surface - soil (19 m here), not the file's 10 m; below
    it the file's thicknesses, an absent layer taking none."""
    m, _notes = _projected(GRB)
    top = np.asarray(m.top)[:5, 0]
    botm = np.asarray(m.botm)[:, :5, 0]
    assert top.tolist() == pytest.approx([19.0] * 5)
    assert botm[0].tolist() == pytest.approx([17, 17, 19, 19, 17])
    assert botm[1].tolist() == pytest.approx([14, 17, 16, 19, 14])
    assert botm[2].tolist() == pytest.approx([4, 4.1, 4.1, 4.2, 4])
    th = np.asarray(m.thick)[:, :5, 0]
    assert th[1].tolist() == pytest.approx([3, 0, 3, 0, 3])


def test_a_mesh_without_layers_is_projected_as_before(GRB):
    """The voronoi and quadtree path: no idomain in the gridprops."""
    gp = MESHES.normalise(MESHES.read_grid_file(GRB))
    gp.pop('idomain'), gp.pop('thick')
    m, _g, _p = MESH.project_model(_FakeMF(), gp, {}, warn=lambda s: None)
    ib = np.asarray(m.ibound)[:, :, 0]
    assert (ib[:, :5] == 1).all() and (ib[:, 5] == 0).all()
    assert np.asarray(m.botm)[:, 0, 0].tolist() == pytest.approx([18, 17, 16])


# ------------------------------------------------------------ pass-through

def _mf6():
    return _load('marmites_mf6_imp',
                 os.path.join(CODE, 'ppMF6', 'marmites_mf6.py')).clsMF6


def test_an_absent_layer_between_two_that_exist_is_a_pass_through():
    """MODFLOW 6 cuts a column at idomain 0: the alluvium over an absent
    terrace would have no connection to the layer below."""
    idom = _mf6().pass_through(EXISTS.reshape(3, 6, 1))[:, :, 0]
    assert idom[:, 1].tolist() == [1, -1, 1]              # 1p1
    # nothing above (cells 3 and 4) or all present: unchanged
    assert idom[:, 2].tolist() == [0, 1, 1]
    assert idom[:, 3].tolist() == [0, 0, 1]
    assert idom[:, 0].tolist() == [1, 1, 1]
    assert int((idom == -1).sum()) == 1


def test_two_layers_never_need_a_pass_through():
    """La Mata: layer 1 absent over layer 2, nothing between."""
    two = np.array([[0, 1, 1], [1, 1, 0]]).reshape(2, 3, 1)
    assert (_mf6().pass_through(two) == two).all()


def test_a_gap_of_two_layers_is_passed_through_and_the_ends_are_not():
    col = np.array([0, 1, 0, 0, 1, 0]).reshape(6, 1, 1)
    assert _mf6().pass_through(col)[:, 0, 0].tolist() == [0, 1, -1, -1, 1, 0]


# ------------------------------------------------------------ the panel

def test_the_grid_panel_asks_for_the_file_and_nothing_of_another_kind():
    sys.path.insert(0, os.path.join(CODE, 'app'))
    from lib import schema
    assert schema.grid_fields('imported')[0] == 'grid.imported.file'
    assert schema.SUBPANELS['grid.imported'] == ('grid.kind', ('imported',))
    label, _units, help_ = schema.describe('grid.imported.file')
    assert '.disv.grb' in label
    assert 'thickness' in help_ and 'idomain' in help_
    assert 'imported' in schema.CHOICES['grid.kind']()
