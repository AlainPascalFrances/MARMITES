# -*- coding: utf-8 -*-
"""Phase-4 tests: CellGeometry abstraction and the DISV-from-DIS grid.

Two things are proven here:
  1. StructuredGeometry reproduces the legacy delr/delc formulas exactly
     (the Phase-4 refactor must not change DIS results);
  2. a DISV grid built from the DIS geometry is geometrically identical --
     same cell areas, same widths for square cells -- which is the basis of
     the DISV-from-DIS equivalence acceptance test.
"""
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..', 'trunk'))
sys.path.insert(0, TRUNK)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


G = _load('marmites_grid', os.path.join(TRUNK, 'marmites_grid.py'))


class _CMF:
    """Minimal stand-in carrying only what the geometry factory reads."""

    def __init__(self, delr, delc, nlay=2, xll=1000.0, yll=2000.0):
        self.delr, self.delc = list(delr), list(delc)
        self.ncol, self.nrow = len(delr), len(delc)
        self.nlay = nlay
        self.xllcorner, self.yllcorner = xll, yll


def _cells(nrow, ncol, active=None):
    """MARMITES cell list (cid, i, j, node) over active (i,j)."""
    out, cid = [], 0
    for i in range(nrow):
        for j in range(ncol):
            if active is None or active(i, j):
                out.append((cid, i, j, i * ncol + j))
                cid += 1
    return out


# --------------------------------------------------------------------- #
# shoelace
# --------------------------------------------------------------------- #

def test_polygon_area_square_and_triangle():
    assert np.isclose(G.polygon_area([(0, 0), (0, 2), (2, 2), (2, 0)]), 4.0)
    assert np.isclose(G.polygon_area([(0, 0), (4, 0), (0, 3)]), 6.0)
    # orientation independent
    assert np.isclose(G.polygon_area([(0, 0), (2, 0), (2, 2), (0, 2)]), 4.0)
    assert np.isclose(G.polygon_area([(0, 0), (1, 1)]), 0.0)


# --------------------------------------------------------------------- #
# StructuredGeometry == legacy formulas
# --------------------------------------------------------------------- #

def test_structured_matches_legacy_formulas():
    delr = [10.0, 20.0, 30.0]        # column widths (x)
    delc = [5.0, 50.0]               # row heights (y)
    cmf = _CMF(delr, delc)
    cells = _cells(2, 3)
    g = G.StructuredGeometry.from_cells(cmf, cells)
    for cid, i, j, _node in cells:
        assert np.isclose(g.area[cid], delc[i] * delr[j]), (i, j)
        assert np.isclose(g.width[cid], delr[j]), (i, j)
    assert g.ncpl == 6 and g.nlay == 2 and g.ncell == 6
    assert g.kind == 'structured'


def test_structured_handles_inactive_cells():
    cmf = _CMF([10.0] * 4, [10.0] * 3)
    cells = _cells(3, 4, active=lambda i, j: (i + j) % 2 == 0)
    g = G.StructuredGeometry.from_cells(cmf, cells)
    assert g.ncell == len(cells) < 12
    assert np.allclose(g.area, 100.0)


# --------------------------------------------------------------------- #
# DISV grid built from DIS
# --------------------------------------------------------------------- #

def test_disv_from_structured_geometry():
    delr = [10.0, 20.0]
    delc = [5.0, 30.0]
    verts, cell2d, ncpl = G.disv_from_structured(delr, delc, xorigin=100.0, yorigin=200.0)
    assert ncpl == 4
    assert len(verts) == (2 + 1) * (2 + 1) == 9
    assert len(cell2d) == 4
    # icell2d == i*ncol + j (the MARMITES 'node')
    assert [r[0] for r in cell2d] == [0, 1, 2, 3]
    # row 0 is the northernmost: its centre y must exceed row 1's
    y0 = cell2d[0][2]
    y1 = cell2d[2][2]
    assert y0 > y1
    # polygon areas equal the structured areas
    for i in range(2):
        for j in range(2):
            ic = i * 2 + j
            ivs = cell2d[ic][4:]
            xy = [(verts[iv][1], verts[iv][2]) for iv in ivs]
            assert np.isclose(G.polygon_area(xy), delc[i] * delr[j]), (i, j)
    # origin honoured: south-west corner of the grid
    xs = [v[1] for v in verts]
    ys = [v[2] for v in verts]
    assert np.isclose(min(xs), 100.0)
    assert np.isclose(min(ys), 200.0)
    assert np.isclose(max(xs), 100.0 + sum(delr))
    assert np.isclose(max(ys), 200.0 + sum(delc))


def test_vertex_geometry_equals_structured_for_square_cells():
    """The Phase-4 correctness baseline: a DISV grid derived from a square
    DIS grid must give identical area AND width (sqrt(area) == delr)."""
    cmf = _CMF([50.0] * 5, [50.0] * 4, nlay=3)
    cells = _cells(4, 5, active=lambda i, j: not (i == 0 and j == 0))
    gs = G.geometry_for(cmf, cells, grid='dis')
    gv = G.geometry_for(cmf, cells, grid='disv')
    assert gv.kind == 'vertex'
    assert np.allclose(gs.area, gv.area)
    assert np.allclose(gs.width, gv.width)      # square cells
    assert gs.ncpl == gv.ncpl == 20


def test_vertex_width_is_sqrt_area_for_rectangular_cells():
    """Documented approximation: for non-square cells the vertex width is
    sqrt(area), not delr -- areas still agree exactly."""
    cmf = _CMF([10.0, 40.0], [10.0, 10.0])
    cells = _cells(2, 2)
    gs = G.geometry_for(cmf, cells, grid='dis')
    gv = G.geometry_for(cmf, cells, grid='disv')
    assert np.allclose(gs.area, gv.area)                 # areas identical
    assert np.allclose(gv.width, np.sqrt(gs.area))       # width differs
    assert not np.allclose(gs.width, gv.width)


def test_geometry_for_rejects_disu():
    cmf = _CMF([10.0], [10.0])
    cells = _cells(1, 1)
    try:
        G.geometry_for(cmf, cells, grid='disu')
    except ValueError as e:
        assert 'disu' in str(e).lower() or 'dis' in str(e).lower()
    else:
        raise AssertionError('DISU must be rejected (out of scope)')
