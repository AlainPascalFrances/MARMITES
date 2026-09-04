# -*- coding: utf-8 -*-
"""Phase-4 tests: the grid-agnostic plotting helper.

The point of the helper is that the SAME call renders a DIS and a DISV
model. These tests check the scatter logic exactly and that both grid types
produce a figure without touching the legacy imshow path.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..', 'trunk'))
sys.path.insert(0, TRUNK)

flopy = pytest.importorskip('flopy')
import matplotlib  # noqa: E402
matplotlib.use('agg')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P = _load('marmites_plot', os.path.join(TRUNK, 'marmites_plot.py'))
G = _load('marmites_grid_plot', os.path.join(TRUNK, 'marmites_grid.py'))


def _cells(nrow, ncol, skip=()):
    out, cid = [], 0
    for i in range(nrow):
        for j in range(ncol):
            if (i, j) in skip:
                continue
            out.append((cid, i, j, i * ncol + j))
            cid += 1
    return out


def test_scatter_places_values_at_nodes():
    cells = _cells(2, 3, skip={(0, 0)})
    vals = np.array([10., 20., 30., 40., 50.])
    arr = P.cells_to_grid_array(cells, vals, ncpl=6)
    assert np.isnan(arr[0])                    # inactive
    assert np.allclose(arr[1:], vals)


def test_scatter_rejects_wrong_length():
    cells = _cells(2, 2)
    with pytest.raises(ValueError):
        P.cells_to_grid_array(cells, np.ones(3), ncpl=4)


def test_plot_structured_grid(tmp_path):
    grid = flopy.discretization.StructuredGrid(
        delr=np.full(4, 50.0), delc=np.full(3, 50.0), nlay=1)
    cells = _cells(3, 4, skip={(0, 0)})
    vals = np.arange(len(cells), dtype=float)
    fn = str(tmp_path / 'dis.png')
    P.plot_cell_values(grid, cells, vals, label='test', fname=fn)
    assert os.path.exists(fn) and os.path.getsize(fn) > 0


def test_plot_vertex_grid(tmp_path):
    verts, cell2d, ncpl = G.disv_from_structured([50.0] * 4, [50.0] * 3)
    grid = flopy.discretization.VertexGrid(
        vertices=verts, cell2d=cell2d, nlay=1, ncpl=ncpl)
    cells = _cells(3, 4, skip={(0, 0)})
    vals = np.arange(len(cells), dtype=float)
    fn = str(tmp_path / 'disv.png')
    P.plot_cell_values(grid, cells, vals, label='test', fname=fn)
    assert os.path.exists(fn) and os.path.getsize(fn) > 0


def test_same_call_serves_both_grids(tmp_path):
    """The Phase-4 promise: one helper, both grid types, no branching by
    the caller."""
    cells = _cells(3, 4)
    vals = np.linspace(0.0, 1.0, len(cells))
    sgrid = flopy.discretization.StructuredGrid(
        delr=np.full(4, 50.0), delc=np.full(3, 50.0), nlay=1)
    verts, cell2d, ncpl = G.disv_from_structured([50.0] * 4, [50.0] * 3)
    vgrid = flopy.discretization.VertexGrid(
        vertices=verts, cell2d=cell2d, nlay=1, ncpl=ncpl)
    for name, g in (('s.png', sgrid), ('v.png', vgrid)):
        P.plot_cell_values(g, cells, vals, label='x', fname=str(tmp_path / name))
        assert os.path.getsize(str(tmp_path / name)) > 0
