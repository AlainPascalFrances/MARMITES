# -*- coding: utf-8 -*-
"""WP1c.7 -- the rasterising display adapter.

The figure suite takes ``(ndays, nlay, nrow, ncol)``. Rather than rewrite it,
mesh results are rasterised onto a display grid at the boundary. Two things
must hold: on a structured grid the adapter is the IDENTITY (so the regression
anchor cannot drift), and on a mesh every cell must actually reach a pixel.
"""

import importlib.util
import os
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


RAST = _load('marmites_rasterise_t',
             os.path.join(CODE, 'ppMF6', 'marmites_rasterise.py'))
GRID = _load('marmites_grid_t3', os.path.join(CODE, 'marmites_grid.py'))
MESH = _load('marmites_mesh_t3', os.path.join(CODE, 'marmites_mesh.py'))

NROW, NCOL, CS = 4, 5, 50.0
XLL, YLL = 1000.0, 2000.0


def _structured_gp(nrow=NROW, ncol=NCOL, cs=CS):
    verts, cell2d, ncpl = GRID.disv_from_structured([cs] * ncol, [cs] * nrow,
                                                    XLL, YLL)
    return {'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl, 'nlay': 2}


def _proj(gp=None):
    gp = gp or _structured_gp()
    return MESH.MeshProjection(gp, [CS] * NCOL, [CS] * NROW, XLL, YLL)


# --------------------------------------------------------------- identity
def test_a_dis_equivalent_mesh_rasterises_back_to_the_source_raster():
    """The load-bearing property: display a mesh that IS the source grid and
    the picture must be the original raster, pixel for pixel."""
    proj = _proj()
    dr = RAST.DisplayRaster.from_projection(proj, refine=1)
    assert (dr.nrow, dr.ncol) == (NROW, NCOL)
    src = np.arange(NROW * NCOL, dtype=float).reshape(NROW, NCOL)
    img = dr.field(src.reshape(-1))
    assert np.array_equal(img, src)
    assert dr.coverage()['cells_missing'] == 0


def test_every_pixel_maps_to_the_cell_it_sits_in():
    proj = _proj()
    dr = RAST.DisplayRaster.from_projection(proj, refine=1)
    for r in range(dr.nrow):
        for c in range(dr.ncol):
            assert dr.pix2cell[r, c] == r * NCOL + c


def test_from_cells_places_a_per_mm_cell_vector():
    proj = _proj()
    dr = RAST.DisplayRaster.from_projection(proj, refine=1)
    cells = [(0, 2, 0, 2), (1, 7, 0, 7)]        # icell2d 2 and 7
    img = dr.from_cells([10.0, 20.0], cells, nodata=-9999.0)
    assert img.reshape(-1)[2] == 10.0
    assert img.reshape(-1)[7] == 20.0
    assert (img == -9999.0).sum() == NROW * NCOL - 2


# ------------------------------------------------------------ refinement
def _two_scale_gp():
    """A 20 m cell beside two 10 m cells -- a mesh finer than the display."""
    v = {0: (0., 0.), 1: (20., 0.), 2: (20., 10.), 3: (20., 20.), 4: (0., 20.),
         5: (30., 0.), 6: (30., 10.), 7: (30., 20.)}
    return {'vertices': [[k, x, y] for k, (x, y) in v.items()],
            'cell2d': [[0, 10., 10., 5, 0, 1, 2, 3, 4],
                       [1, 25., 5., 4, 1, 5, 6, 2],
                       [2, 25., 15., 4, 2, 6, 7, 3]],
            'ncpl': 3, 'nlay': 1}


def test_a_cell_smaller_than_a_pixel_is_reported_as_missing():
    """The failure this catches is invisible in a picture: La Mata's quadtree
    refines to 12.5 m, and on the 50 m source grid 3828 of its 7728 cells fall
    between pixel centres and appear in NO figure at all."""
    gp = _two_scale_gp()
    proj = MESH.MeshProjection(gp, [30.0], [20.0], 0.0, 0.0)   # ONE 30x20 pixel
    dr = RAST.DisplayRaster(gp, proj._x_edge, proj._y_edge, refine=1)
    cov = dr.coverage()
    assert cov['cells_shown'] == 1
    assert cov['cells_missing'] == 2


def test_auto_refine_gives_every_cell_a_pixel():
    gp = _two_scale_gp()
    proj = MESH.MeshProjection(gp, [30.0], [20.0], 0.0, 0.0)
    k = RAST.DisplayRaster.auto_refine(proj)
    assert k > 1
    dr = RAST.DisplayRaster.from_projection(proj, refine=0)   # 0 = auto
    assert dr.coverage()['cells_missing'] == 0


def test_refinement_is_capped():
    assert RAST.DisplayRaster.MAX_REFINE >= 2
    gp = _two_scale_gp()
    proj = MESH.MeshProjection(gp, [30000.0], [20000.0], 0.0, 0.0)
    assert RAST.DisplayRaster.auto_refine(proj) <= RAST.DisplayRaster.MAX_REFINE


# --------------------------------------------------------------- adapter
class _StructMF:
    nlay, nrow, ncol = 2, NROW, NCOL
    delr = [CS] * NCOL
    delc = [CS] * NROW
    xllcorner, yllcorner = XLL, YLL
    hnoflo = 9999.999


def test_the_adapter_is_the_identity_on_a_structured_grid():
    """So the structured figures cannot drift when the mesh path changes."""
    DA = RAST.MapAdapter(_StructMF())
    assert not DA.on_mesh
    assert (DA.nrow, DA.ncol) == (NROW, NCOL)
    a = np.arange(2 * NROW * NCOL, dtype=float).reshape(2, NROW, NCOL)
    assert np.array_equal(DA.lay(a), a)
    cells = [(0, 1, 2, 1 * NCOL + 2)]
    g = DA.cells([5.0], cells, nodata=-1.0)
    assert g[1, 2] == 5.0
    assert np.allclose(DA.cell_area(), CS * CS)


def test_the_adapter_rasterises_on_a_mesh():
    proj = _proj()

    class _MeshMF:
        nlay, nrow, ncol = 2, proj.ncpl, 1
        delr, delc = [1.0], [1.0] * proj.ncpl
        xllcorner, yllcorner = XLL, YLL
        hnoflo = 9999.999
        mesh_proj = proj

    DA = RAST.MapAdapter(_MeshMF(), refine=1)
    assert DA.on_mesh
    assert (DA.nrow, DA.ncol) == (NROW, NCOL)
    a = np.arange(2 * proj.ncpl, dtype=float).reshape(2, proj.ncpl)
    out = DA.lay(a)
    assert out.shape == (2, NROW, NCOL)
    assert np.array_equal(out[0], a[0].reshape(NROW, NCOL))
    # cell_area is the MODEL cell's area, not the pixel's
    assert np.allclose(DA.cell_area(), CS * CS)


def test_cell_area_is_the_model_cell_not_the_pixel():
    """The m3/d -> mm/d conversion divides by the area the flux was computed
    in. Using the pixel area instead would scale every aquifer flux map by the
    refinement factor squared."""
    gp = _two_scale_gp()
    proj = MESH.MeshProjection(gp, [30.0], [20.0], 0.0, 0.0)

    class _MeshMF:
        nlay, nrow, ncol = 1, 3, 1
        delr, delc = [1.0], [1.0, 1.0, 1.0]
        xllcorner = yllcorner = 0.0
        hnoflo = 9999.999
        mesh_proj = proj

    DA = RAST.MapAdapter(_MeshMF(), refine=4)
    area = DA.cell_area()
    vals = np.unique(area[np.isfinite(area)])
    # the coarse cell is 400 m2, the two refined ones 100 m2 -- never the
    # pixel area, which is much smaller at refine=4
    assert set(np.round(vals).astype(int).tolist()) <= {100, 400}


def test_field_refuses_a_wrongly_sized_vector():
    proj = _proj()
    dr = RAST.DisplayRaster.from_projection(proj, refine=1)
    with pytest.raises(ValueError):
        dr.field(np.zeros(7))
