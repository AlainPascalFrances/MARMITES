# -*- coding: utf-8 -*-
"""WP1c -- projecting a structured model onto an unstructured mesh.

The load-bearing test is `test_dis_equivalent_projection_is_an_identity`: a
mesh that reproduces the structured grid cell for cell must give back the
model unchanged. Everything else in the mesh path -- sampling, masking, cell
ordering, boundary remapping -- is downstream of that, so if it holds, a real
mesh differs from the structured model only by discretisation.
"""

import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


MESH = _load('marmites_mesh_t', os.path.join(CODE, 'marmites_mesh.py'))
GRID = _load('marmites_grid_t', os.path.join(CODE, 'marmites_grid.py'))
MESHES = _load('marmites_meshes_t', os.path.join(CODE, 'marmites_meshes.py'))

NROW, NCOL, NLAY = 4, 5, 2
XLL, YLL, CS = 1000.0, 2000.0, 50.0
HNOFLO = 9999.999


class _Process:
    def __init__(self, nrow, ncol, nlay):
        self.nrow, self.ncol, self.nlay = nrow, ncol, nlay

    def float2array(self, array):
        a = np.asarray(array)
        if a.ndim == 1 and a.shape[0] == self.nlay:
            out = np.ones((self.nlay, self.nrow, self.ncol))
            for k, e in enumerate(a):
                out[k] *= e
            return out
        return a


class _FakeMF:
    """A tiny structured MARMITES model, shaped like the real cMF."""

    def __init__(self):
        self.nrow, self.ncol, self.nlay = NROW, NCOL, NLAY
        self.delr = [CS] * NCOL
        self.delc = [CS] * NROW
        self.xllcorner, self.yllcorner = XLL, YLL
        self.hnoflo = HNOFLO
        rng = np.random.default_rng(7)
        self.ibound = np.ones((NLAY, NROW, NCOL), dtype=int)
        self.ibound[:, 0, :] = 0                      # a dead north row
        self.outcropL = np.zeros((NROW, NCOL), dtype=int)
        for L in range(NLAY):
            ib = np.abs(self.ibound)[L] != 0
            self.outcropL += ((self.outcropL == 0) & ib) * (L + 1)
        self.iuzfbnd = (self.outcropL > 0).astype(int)
        elev = 700.0 + rng.random((NROW, NCOL)) * 10.0
        self.elev = np.ma.masked_array(elev, mask=(self.outcropL == 0))
        thick = np.full((NROW, NCOL), 2.0)
        # the real driver builds top by masked subtraction, which leaves REAL
        # data under the mask -- reproduce that, it is what the trap needs
        self.top = self.elev - np.ma.masked_array(thick, mask=np.zeros_like(thick, bool))
        self.botm = np.stack([elev - 10.0 * (k + 1) for k in range(NLAY)])
        self.strt = np.stack([elev - 1.0 for _ in range(NLAY)])
        self.hk_actual = [rng.random((NROW, NCOL)) + 1.0 for _ in range(NLAY)]
        self.vka_actual = [2.0] * NLAY
        self.sy_actual = [rng.random((NROW, NCOL)) * 0.1 for _ in range(NLAY)]
        self.ss_actual = [1e-5] * NLAY
        self.cPROCESS = _Process(NROW, NCOL, NLAY)
        # two drains, each 1.5 m above its own cell bottom (the La Mata rule)
        self.layer_row_column_elevation_cond = {0: [
            [0, 1, 1, float(self.botm[0, 1, 1]) + 1.5, 0.035],
            [1, 2, 3, float(self.botm[1, 2, 3]) + 1.5, 0.025],
        ]}
        self.layer_row_column_head_cond = None


def _dis_equivalent_gridprops(cMF):
    verts, cell2d, ncpl = GRID.disv_from_structured(
        cMF.delr, cMF.delc, cMF.xllcorner, cMF.yllcorner)
    return {'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl,
            'nlay': cMF.nlay}


def _grids(cMF):
    rng = np.random.default_rng(11)
    return {'gridSOIL': rng.integers(1, 4, (NROW, NCOL)),
            'gridMETEO': np.ones((NROW, NCOL), dtype=int),
            'gridSOILthick': np.full((NROW, NCOL), 2.0),
            'gridSsurfhmax': rng.random((NROW, NCOL)),
            'gridSsurfw': rng.random((NROW, NCOL)),
            'gridIRR': np.zeros((NROW, NCOL), dtype=int),
            'gridVEGarea': rng.random((3, NROW, NCOL)) * 100.0}


# ------------------------------------------------------------------ locate
def test_locate_maps_points_to_the_right_source_cell():
    cMF = _FakeMF()
    gp = _dis_equivalent_gridprops(cMF)
    proj = MESH.MeshProjection(gp, cMF.delr, cMF.delc, XLL, YLL)
    # row 0 is the NORTHERNMOST row, so the top-left cell centre is
    # (xll + cs/2, yll + nrow*cs - cs/2)
    r, c, inside = proj.locate(XLL + 0.5 * CS, YLL + NROW * CS - 0.5 * CS)
    assert (int(r[0]), int(c[0]), bool(inside[0])) == (0, 0, True)
    r, c, inside = proj.locate(XLL + (NCOL - 0.5) * CS, YLL + 0.5 * CS)
    assert (int(r[0]), int(c[0]), bool(inside[0])) == (NROW - 1, NCOL - 1, True)


def test_points_outside_the_source_grid_are_flagged():
    cMF = _FakeMF()
    proj = MESH.MeshProjection(_dis_equivalent_gridprops(cMF),
                               cMF.delr, cMF.delc, XLL, YLL)
    _r, _c, inside = proj.locate([XLL - 1000.0, XLL + 10.0],
                                 [YLL + 10.0, YLL + 10.0])
    assert list(inside) == [False, True]


def test_cell_of_finds_the_nearest_centroid():
    cMF = _FakeMF()
    gp = _dis_equivalent_gridprops(cMF)
    proj = MESH.MeshProjection(gp, cMF.delr, cMF.delc, XLL, YLL)
    for ic in (0, 7, NROW * NCOL - 1):
        x, y = proj.xy[ic]
        assert proj.cell_of(x + 1.0, y - 1.0) == ic


# ------------------------------------------------------------- the identity
def test_dis_equivalent_projection_is_an_identity():
    """A mesh identical to the structured grid must return the same model.

    This is the whole correctness argument for the mesh path: it isolates the
    projection from any real difference in discretisation.
    """
    cMF = _FakeMF()
    grids = _grids(cMF)
    m, g2, _proj = MESH.project_model(cMF, _dis_equivalent_gridprops(cMF), grids)

    assert (m.nrow, m.ncol) == (NROW * NCOL, 1)
    flat = lambda a: np.asarray(a).reshape(-1, 1)          # noqa: E731

    assert np.array_equal(flat(cMF.outcropL), m.outcropL)
    assert np.array_equal(flat(cMF.iuzfbnd), m.iuzfbnd)
    for k in range(NLAY):
        assert np.array_equal(flat(cMF.ibound[k]), m.ibound[k])
        assert np.allclose(flat(cMF.botm[k]), m.botm[k])
        assert np.allclose(flat(cMF.hk_actual[k]), m.hk_actual[k])
    # np.asarray() drops the mask and keeps the DATA -- which is exactly what
    # the driver hands MODFLOW, so it is the data that has to match
    assert np.allclose(flat(np.asarray(cMF.top)), np.asarray(m.top))
    assert np.allclose(flat(np.asarray(cMF.elev)), np.asarray(m.elev))
    for name, src in grids.items():
        a = np.asarray(src)
        if a.ndim == 3:
            for z in range(a.shape[0]):
                assert np.allclose(flat(a[z]), g2[name][z]), name
        else:
            assert np.allclose(flat(a), g2[name]), name


def test_masks_survive_the_projection():
    cMF = _FakeMF()
    m, _g, _p = MESH.project_model(cMF, _dis_equivalent_gridprops(cMF), {})
    src = np.ma.getmaskarray(cMF.elev).reshape(-1, 1)
    assert np.array_equal(src, np.ma.getmaskarray(m.elev))


def test_mask_sentinel_does_not_overwrite_data_under_the_mask():
    """np.ma.masked_values() starts with filled(x, value), so it DESTROYS the
    data under an existing mask. The driver passes np.asarray(cMF.top) to
    MODFLOW, which keeps that data -- so this is a real difference, not a
    presentational one."""
    a = np.ma.masked_array([1.0, 2.0, 3.0], mask=[False, True, False])
    trap = np.ma.masked_values(a, HNOFLO, atol=0.09)
    assert np.asarray(trap)[1] == HNOFLO           # the trap, documented
    kept = MESH._mask_sentinel(a, HNOFLO)
    assert np.asarray(kept)[1] == 2.0              # data preserved
    assert bool(np.ma.getmaskarray(kept)[1])       # still masked


def test_hnoflo_is_masked_even_when_the_source_was_unmasked():
    plain = np.array([1.0, HNOFLO, 3.0])
    out = MESH._mask_sentinel(plain, HNOFLO)
    assert list(np.ma.getmaskarray(out)) == [False, True, False]


# ------------------------------------------------------------- cell list
def test_cell_list_indices_resolve_against_the_column_vectors():
    """Under the (ncpl, 1) convention the row index IS the icell2d, which is
    what lets build_cell_list, clsMF6._cellid and _griddata work unchanged."""
    cMF = _FakeMF()
    m, _g, _p = MESH.project_model(cMF, _dis_equivalent_gridprops(cMF), {})
    cells = [(cid, i, 0, i) for cid, i in
             enumerate(np.flatnonzero(m.outcropL[:, 0] > 0))]
    assert len(cells) == int((cMF.outcropL > 0).sum())
    for _cid, i, j, node in cells:
        assert i == node and j == 0
        assert m.outcropL[i, j] > 0


def test_griddata_reshape_holds():
    """nrow * ncol == ncpl is what makes clsMF6._griddata a valid reshape."""
    cMF = _FakeMF()
    m, _g, proj = MESH.project_model(cMF, _dis_equivalent_gridprops(cMF), {})
    assert m.nrow * m.ncol == proj.ncpl
    assert np.asarray(m.ibound).reshape(NLAY, proj.ncpl).shape == (NLAY, proj.ncpl)


def test_float2array_reports_the_mesh_shape():
    """MMsoil reads Sy as float2array(sy)[outcrop-1, i, j]; a cPROCESS still
    claiming the raster shape would index the wrong array."""
    cMF = _FakeMF()
    m, _g, proj = MESH.project_model(cMF, _dis_equivalent_gridprops(cMF), {})
    assert m.cPROCESS.float2array(m.sy_actual).shape == (NLAY, proj.ncpl, 1)
    assert m.cPROCESS.float2array(m.ss_actual).shape == (NLAY, proj.ncpl, 1)


# ------------------------------------------------------------------ drains
def test_drains_keep_their_height_above_the_cell_bottom():
    """La Mata's drains sit a fixed 1.51 m above their cell bottom. Carrying
    the ABSOLUTE elevation across makes MODFLOW refuse the package with
    'DRN BOUNDARY ELEVATION IS LESS THAN CELL BOTTOM'."""
    cMF = _FakeMF()
    gp = _dis_equivalent_gridprops(cMF)
    # a mesh whose bottoms differ from the source: shift botm down by 5 m
    m, _g, proj = MESH.project_model(cMF, gp, {})
    src = cMF.layer_row_column_elevation_cond[0]
    got = m.layer_row_column_elevation_cond[0]
    assert len(got) == len(src)
    for s, g in zip(src, got):
        lay, i, j = int(s[0]), int(s[1]), int(s[2])
        ic = int(g[1])
        assert ic == i * NCOL + j and int(g[2]) == 0
        off_src = float(s[3]) - float(cMF.botm[lay, i, j])
        off_dst = float(g[3]) - float(m.botm[lay, ic, 0])
        assert abs(off_src - off_dst) < 1e-9
        assert float(g[4]) == float(s[4])          # conductance unchanged


def test_drain_elevation_follows_a_different_mesh_bottom():
    cMF = _FakeMF()
    gp = _dis_equivalent_gridprops(cMF)
    proj = MESH.MeshProjection(gp, cMF.delr, cMF.delc, XLL, YLL)
    src_botm = np.asarray(cMF.botm, dtype=float)
    mesh_botm = proj.sample3d(src_botm, fill=HNOFLO) - 5.0    # deeper mesh
    out = proj.remap_drn_records(cMF.layer_row_column_elevation_cond[0],
                                 src_botm, mesh_botm, warn=None)
    for s, g in zip(cMF.layer_row_column_elevation_cond[0], out):
        assert abs(float(g[3]) - (float(s[3]) - 5.0)) < 1e-9


def test_clashing_boundary_records_are_reported():
    """Two source drains landing in one mesh cell is a real change to the
    boundary condition, so it must be said rather than silently applied."""
    cMF = _FakeMF()
    gp = _dis_equivalent_gridprops(cMF)
    proj = MESH.MeshProjection(gp, cMF.delr, cMF.delc, XLL, YLL)
    said = []
    proj.remap_records([[0, 1, 1, 1.0, 2.0], [0, 1, 1, 3.0, 4.0]],
                       warn=said.append)
    assert said and 'share a mesh cell' in said[0]


# ------------------------------------------------------------------ guards
def test_zone_rasters_are_never_interpolated():
    cMF = _FakeMF()
    gp = _dis_equivalent_gridprops(cMF)
    proj = MESH.MeshProjection(gp, cMF.delr, cMF.delc, XLL, YLL)
    zones = np.array([[1, 3] * (NCOL // 2 + 1)] * NROW)[:, :NCOL]
    out = proj.sample2d(zones, fill=0, dtype=int)
    assert out.dtype.kind in 'iu'
    assert set(np.unique(out)) <= {1, 3}


def test_an_unknown_sampling_mode_fails_loudly():
    cMF = _FakeMF()
    proj = MESH.MeshProjection(_dis_equivalent_gridprops(cMF),
                               cMF.delr, cMF.delc, XLL, YLL)
    with pytest.raises(MESH.MeshProjectionError) as e:
        proj.sample2d(np.zeros((NROW, NCOL)), how='bicubic')
    assert 'unknown sampling mode' in str(e.value)


def test_a_per_field_rule_is_refused_as_a_whole_model_policy():
    """'majority' on an elevation would round it to whole metres, and 'area'
    on a zone raster would interpolate class codes. Only 'auto' and 'centre'
    mean anything applied to every field at once."""
    cMF = _FakeMF()
    for bad in ('area', 'majority'):
        with pytest.raises(MESH.MeshProjectionError) as e:
            MESH.project_model(cMF, _dis_equivalent_gridprops(cMF), {}, how=bad)
        assert 'whole-model policy' in str(e.value)


def test_a_wrongly_shaped_source_array_is_refused():
    cMF = _FakeMF()
    proj = MESH.MeshProjection(_dis_equivalent_gridprops(cMF),
                               cMF.delr, cMF.delc, XLL, YLL)
    with pytest.raises(MESH.MeshProjectionError):
        proj.sample2d(np.zeros((NROW + 1, NCOL)))


def test_a_mesh_that_misses_the_model_is_refused():
    """A CRS or origin mistake must stop the run, not produce an empty model."""
    cMF = _FakeMF()
    gp = _dis_equivalent_gridprops(cMF)
    far = {'vertices': [[i, x + 1e6, y + 1e6] for i, x, y in gp['vertices']],
           'cell2d': [[r[0], r[1] + 1e6, r[2] + 1e6] + list(r[3:])
                      for r in gp['cell2d']],
           'ncpl': gp['ncpl'], 'nlay': gp['nlay']}
    with pytest.raises(MESH.MeshProjectionError) as e:
        MESH.project_model(cMF, far, {})
    assert 'do not overlap' in str(e.value)


def test_the_source_model_is_left_untouched():
    """A projection must not damage the structured model it read."""
    cMF = _FakeMF()
    before = (cMF.nrow, cMF.ncol, np.asarray(cMF.ibound).copy(),
              np.asarray(cMF.botm).copy())
    MESH.project_model(cMF, _dis_equivalent_gridprops(cMF), _grids(cMF))
    assert (cMF.nrow, cMF.ncol) == before[:2]
    assert np.array_equal(np.asarray(cMF.ibound), before[2])
    assert np.array_equal(np.asarray(cMF.botm), before[3])


# ------------------------------------------------------- producers / cache
def test_stream_lines_groups_by_segment(tmp_path):
    p = tmp_path / 's.csv'
    p.write_text('# provenance\nseg_id,seq,x,y\n'
                 '0,1,10,20\n0,0,0,0\n1,0,5,5\n1,1,6,6\n2,0,9,9\n',
                 encoding='utf-8')
    segs = MESHES.stream_lines(str(p))
    assert segs == [[(0.0, 0.0), (10.0, 20.0)], [(5.0, 5.0), (6.0, 6.0)]]


def test_watershed_ring_drops_the_closing_point(tmp_path):
    p = tmp_path / 'w.csv'
    p.write_text('# provenance\nring_id,seq,x,y\n'
                 '0,0,0,0\n0,1,10,0\n0,2,10,10\n0,3,0,0\n', encoding='utf-8')
    ring = MESHES.watershed_ring(str(p))
    assert ring == [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)]


def test_normalise_makes_a_fresh_mesh_equal_a_cached_one(tmp_path):
    """A cached mesh must be the SAME object as a fresh one, or a bug can
    appear on the second run only."""
    gp = {'vertices': [[np.int64(0), np.float64(1.0), 2.0]],
          'cell2d': [[np.int64(0), np.float64(1.0), 2.0, np.int64(1),
                      np.int64(0)]],
          'ncpl': np.int64(1), 'nlay': 2}
    norm = MESHES.normalise(gp)
    MESHES._save_cached(str(tmp_path), 'voronoi', 'sig', norm)
    back = MESHES._load_cached(str(tmp_path), 'voronoi', 'sig')
    assert back['vertices'] == norm['vertices']
    assert back['cell2d'] == norm['cell2d']
    assert all(type(x) is int for x in (norm['cell2d'][0][0], norm['ncpl']))


def test_a_stale_cache_is_not_served(tmp_path):
    gp = MESHES.normalise({'vertices': [[0, 1.0, 2.0]],
                           'cell2d': [[0, 1.0, 2.0, 1, 0]], 'ncpl': 1,
                           'nlay': 1})
    MESHES._save_cached(str(tmp_path), 'voronoi', 'sig-A', gp)
    assert MESHES._load_cached(str(tmp_path), 'voronoi', 'sig-A') is not None
    assert MESHES._load_cached(str(tmp_path), 'voronoi', 'sig-B') is None


def test_cell_size_report_measures_real_polygons():
    cMF = _FakeMF()
    rep = MESHES.cell_size_report(_dis_equivalent_gridprops(cMF))
    assert abs(rep['area_mean'] - CS * CS) < 1e-6
    assert abs(rep['size_equiv'] - CS) < 1e-6


# ------------------------------------------------------- WP1c.3 resampling
def _coarse_gridprops(cMF, factor=2):
    """A mesh of `factor` x `factor` source cells per mesh cell.

    Coarsening is the case that matters: centre sampling then keeps one source
    cell out of factor**2 and discards the rest.
    """
    delr = [CS * factor] * (NCOL // factor)
    delc = [CS * factor] * (NROW // factor)
    verts, cell2d, ncpl = GRID.disv_from_structured(delr, delc, XLL, YLL)
    return {'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl,
            'nlay': cMF.nlay}, delr, delc


def test_clip_to_rect_computes_exact_overlap_area():
    square = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    assert abs(MESH._shoelace(MESH._clip_to_rect(square, 0, 0, 10, 10))
               - 100.0) < 1e-9
    assert abs(MESH._shoelace(MESH._clip_to_rect(square, 5, 0, 15, 10))
               - 50.0) < 1e-9                       # half overlap
    assert abs(MESH._shoelace(MESH._clip_to_rect(square, 5, 5, 15, 15))
               - 25.0) < 1e-9                       # quarter overlap
    assert MESH._clip_to_rect(square, 20, 20, 30, 30) == []   # disjoint
    tri = [(0.0, 0.0), (10.0, 0.0), (0.0, 10.0)]
    assert abs(MESH._shoelace(MESH._clip_to_rect(tri, -5, -5, 15, 15))
               - 50.0) < 1e-9                       # fully inside


def test_overlap_weights_partition_each_mesh_cell():
    """Every mesh cell's overlaps must sum to its own area, or the weighted
    mean is normalised by the wrong total."""
    cMF = _FakeMF()
    gp, _dr, _dc = _coarse_gridprops(cMF)
    proj = MESH.MeshProjection(gp, cMF.delr, cMF.delc, XLL, YLL)
    for ic, (_r, _c, w) in enumerate(proj.overlaps()):
        assert abs(w.sum() - MESH._shoelace(proj.cell_polygon(ic))) < 1e-6
    rep = proj.overlap_report()
    assert abs(rep['coverage_mean'] - 1.0) < 1e-9
    assert abs(rep['src_per_cell_mean'] - 4.0) < 1e-9      # 2x2 source cells


def test_area_weighting_is_exactly_mass_conserving():
    """The acceptance criterion. When the mesh tiles the source grid and
    nothing is excluded, sum(value * area) must be preserved to round-off --
    that is an identity, so any drift is a bug rather than a tolerance."""
    cMF = _FakeMF()
    gp, _dr, _dc = _coarse_gridprops(cMF)
    proj = MESH.MeshProjection(gp, cMF.delr, cMF.delc, XLL, YLL)
    src = np.random.default_rng(3).random((NROW, NCOL)) * 100.0
    areas = np.array([MESH._shoelace(proj.cell_polygon(i))
                      for i in range(proj.ncpl)])
    # the coarse mesh covers only the tiled part of the grid, so compare there
    ref = float((src[:len(_dc) * 2, :len(_dr) * 2] * CS * CS).sum())
    out = np.asarray(proj.sample2d(src, fill=np.nan, dtype=float, how='area'))
    assert abs(float((out[:, 0] * areas).sum()) - ref) / ref < 1e-12


def test_area_weighting_beats_centre_sampling_on_a_coarser_mesh():
    """Centre sampling keeps one source cell in four; the point of WP1c.3 is
    that it stops doing that."""
    cMF = _FakeMF()
    gp, _dr, _dc = _coarse_gridprops(cMF)
    proj = MESH.MeshProjection(gp, cMF.delr, cMF.delc, XLL, YLL)
    src = np.random.default_rng(5).random((NROW, NCOL)) * 100.0
    ref = float(src[:len(_dc) * 2, :len(_dr) * 2].mean())
    got_c = float(np.asarray(proj.sample2d(src, how='centre')).mean())
    got_a = float(np.asarray(proj.sample2d(src, how='area')).mean())
    assert abs(got_a - ref) < 1e-12                  # exact: uniform cells
    assert abs(got_c - ref) > abs(got_a - ref)


def test_majority_picks_the_class_with_the_largest_overlap():
    cMF = _FakeMF()
    gp, _dr, _dc = _coarse_gridprops(cMF)
    proj = MESH.MeshProjection(gp, cMF.delr, cMF.delc, XLL, YLL)
    zones = np.ones((NROW, NCOL), dtype=int)
    zones[0, 0] = 7                       # 1 of the 4 source cells of cell 0
    out = np.asarray(proj.sample2d(zones, fill=0, dtype=int, how='majority'))
    assert out[0, 0] == 1                 # 3 beats 1
    zones[0, 1] = zones[1, 0] = 7         # now 3 of 4
    out = np.asarray(proj.sample2d(zones, fill=0, dtype=int, how='majority'))
    assert out[0, 0] == 7


def test_a_nodata_sentinel_never_wins_a_vote_or_skews_a_mean():
    """The bug this guards: La Mata's zone rasters hold hnoflo outside the
    catchment, and a boundary mesh cell elected 9999 as its soil zone --
    which moved the mean soil zone from 1.9 to 1248."""
    cMF = _FakeMF()
    gp, _dr, _dc = _coarse_gridprops(cMF)
    proj = MESH.MeshProjection(gp, cMF.delr, cMF.delc, XLL, YLL)
    valid = np.ones((NROW, NCOL), dtype=bool)
    zones = np.ones((NROW, NCOL), dtype=int)
    zones[0, 0] = zones[0, 1] = zones[1, 0] = 9999      # 3 of 4, but nodata
    valid[zones == 9999] = False
    out = np.asarray(proj.sample2d(zones, fill=0, dtype=int, how='majority',
                                   valid=valid))
    assert out[0, 0] == 1, 'the sentinel won the vote'
    cont = np.ones((NROW, NCOL)) * 2.0
    cont[0, 0] = cont[0, 1] = cont[1, 0] = 9999.999
    got = np.asarray(proj.sample2d(cont, fill=np.nan, dtype=float, how='area',
                                   valid=valid))
    assert abs(got[0, 0] - 2.0) < 1e-9, 'the sentinel was averaged in'


def test_auto_picks_the_rule_from_the_field_kind():
    assert MESH.MeshProjection._resolve_how('auto', int) == 'majority'
    assert MESH.MeshProjection._resolve_how('auto', float) == 'area'
    assert MESH.MeshProjection._resolve_how('centre', int) == 'centre'
    with pytest.raises(MESH.MeshProjectionError):
        MESH.MeshProjection._resolve_how('bilinear', float)


@pytest.mark.parametrize('how', list(MESH.MODEL_SAMPLING_MODES))
def test_the_identity_holds_under_every_sampling_mode(how):
    """On a DIS-equivalent mesh each mesh cell IS one source cell, so every
    mode must collapse to the same answer. That is what lets the mode change
    without re-arguing the correctness of the projection."""
    cMF = _FakeMF()
    grids = _grids(cMF)
    m, g2, _p = MESH.project_model(cMF, _dis_equivalent_gridprops(cMF), grids,
                                   how=how)
    flat = lambda a: np.asarray(a).reshape(-1, 1)          # noqa: E731
    assert np.array_equal(flat(cMF.outcropL), m.outcropL)
    assert np.allclose(flat(np.asarray(cMF.top)), np.asarray(m.top))
    for k in range(NLAY):
        assert np.allclose(flat(cMF.botm[k]), m.botm[k])
    for name, src in grids.items():
        a = np.asarray(src)
        if a.ndim == 2:
            assert np.allclose(flat(a), g2[name]), '%s under %s' % (name, how)


def test_a_generator_on_the_domain_edge_still_counts_as_inside():
    """flopy's VoronoiGrid reports a cell's GENERATOR as its centre, and for a
    boundary cell that point lies exactly ON the domain edge. A strict
    `x < xmax` test then declares much of a clipped mesh to be off-grid."""
    cMF = _FakeMF()
    proj = MESH.MeshProjection(_dis_equivalent_gridprops(cMF),
                               cMF.delr, cMF.delc, XLL, YLL)
    xmax = XLL + NCOL * CS
    _r, _c, inside = proj.locate([xmax, XLL, xmax],
                                 [YLL, YLL, YLL + NROW * CS])
    assert list(inside) == [True, True, True]


def test_domain_ratio_is_reported():
    cMF = _FakeMF()
    proj = MESH.MeshProjection(_dis_equivalent_gridprops(cMF),
                               cMF.delr, cMF.delc, XLL, YLL)
    assert abs(proj.overlap_report()['domain_ratio'] - 1.0) < 1e-9


def test_config_resample_modes_match_the_projection():
    """The schema and the sampler must not drift apart."""
    import marmites_config as cfgmod
    assert set(cfgmod.RESAMPLE_MODES) == set(MESH.MODEL_SAMPLING_MODES)
    assert set(MESH.MODEL_SAMPLING_MODES) <= set(MESH.SAMPLING_MODES)
