# -*- coding: utf-8 -*-
"""Quadtree (refined) DISV grids for MARMITES via GRIDGEN (Phase 4).

`marmites_grid.disv_from_structured` gives a vertex grid that is geometrically
identical to the DIS grid -- the correctness baseline. This module produces a
genuinely *refined* quadtree grid with USGS GRIDGEN (Lien & others, OFR
2014-1109) through flopy's wrapper, and rebuilds the MARMITES inputs on it.

Why a separate input layout
---------------------------
On a refined grid a cell has no (row, column). MARMITES' per-cell kernel
indexes its inputs as ``grid[i, j]``, so the refined path uses the degenerate
layout agreed in the design (code review section 5):

    cell list : (cid, cid, 0, icell2d)      -- i := cid, j := 0
    inputs    : column vectors of shape (ncell, 1)

Every existing ``grid[i, j]`` lookup then resolves to ``grid[cid, 0]`` and the
soil physics is reached unchanged -- ``flux()`` keeps its exact signature and
the frozen-oracle tests keep comparing like for like.

Requires the gridgen executable (not bundled):
    C:\\00MODFLOW\\gridgen.1.0.02\\bin\\gridgen_x64.exe   (Windows)
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import os

import numpy as np

from marmites_grid import VertexGeometry
from marmites_raster import sample_to_cells

__all__ = ['build_quadtree', 'cell_centres_from_gridprops', 'refined_cell_list',
           'resample_inputs', 'RefinedModel']


def build_quadtree(cMF, gridgen_exe, model_ws, refine_features=None,
                   level=1, layers=None, surface_interpolation='replicate'):
    """Run GRIDGEN and return DISV gridprops for a quadtree refinement.

    Parameters
    ----------
    cMF : parsed clsMF (supplies delr/delc/nlay/origin and top/botm shape).
    gridgen_exe : str      path to gridgen_x64(.exe).
    model_ws : str         scratch workspace for gridgen files.
    refine_features : list of feature specs, each
        ``(geometries, featuretype, level)`` where featuretype is
        'point' | 'line' | 'polygon'. Geometries are lists of coordinate
        sequences in model CRS. When None, no refinement is applied (the
        result is then equivalent to the base grid).
    level : int            default refinement level for features given
                           without their own level.
    layers : list[int]     layers to refine (default: all).

    Returns the flopy ``get_gridprops_disv()`` dictionary.
    """
    from flopy.discretization import StructuredGrid
    from flopy.utils.gridgen import Gridgen

    if not os.path.exists(gridgen_exe):
        raise FileNotFoundError('gridgen executable not found: %s' % gridgen_exe)
    nlay = int(cMF.nlay)
    sg = StructuredGrid(delr=np.asarray(cMF.delr, dtype=float),
                        delc=np.asarray(cMF.delc, dtype=float),
                        top=np.asarray(cMF.top, dtype=float),
                        botm=np.asarray(cMF.botm, dtype=float),
                        nlay=nlay,
                        xoff=float(getattr(cMF, 'xllcorner', 0.0)),
                        yoff=float(getattr(cMF, 'yllcorner', 0.0)))
    os.makedirs(model_ws, exist_ok=True)
    g = Gridgen(sg, model_ws=model_ws, exe_name=gridgen_exe,
                surface_interpolation=surface_interpolation)
    if layers is None:
        layers = list(range(nlay))
    for feat in (refine_features or []):
        if len(feat) == 3:
            geoms, ftype, lev = feat
        else:
            geoms, ftype = feat
            lev = level
        g.add_refinement_features(geoms, ftype, int(lev), layers)
    g.build(verbose=False)
    return g.get_gridprops_disv()


def cell_centres_from_gridprops(gridprops):
    """(ncpl, 2) cell-centre coordinates from DISV gridprops."""
    cell2d = gridprops['cell2d']
    xy = np.array([[float(r[1]), float(r[2])] for r in cell2d], dtype=float)
    return xy


def refined_cell_list(active_nodes):
    """MARMITES cell list for a refined grid: (cid, cid, 0, icell2d).

    ``i := cid`` and ``j := 0`` make the legacy ``[i, j]`` indexing resolve
    into (ncell, 1) column-vector inputs.
    """
    return [(cid, cid, 0, int(node)) for cid, node in enumerate(active_nodes)]


def _col(arr):
    """Shape a per-cell vector as an (ncell, 1) column vector."""
    a = np.asarray(arr)
    return a.reshape(-1, 1)


def resample_inputs(gridprops, active_nodes, raster_dir, rasters, hnoflo):
    """Sample MARMITES spatial inputs onto refined cells.

    rasters : dict name -> (filename, dtype). Returns dict name -> (ncell, 1).
    Zone rasters must be given dtype=int (nearest sampling, no interpolation).
    """
    xy_all = cell_centres_from_gridprops(gridprops)
    xy = xy_all[np.asarray(active_nodes, dtype=int)]
    out = {}
    for name, (fn, dt) in rasters.items():
        path = os.path.join(raster_dir, fn)
        vals = sample_to_cells(path, xy, dtype=dt, hnoflo=hnoflo)
        out[name] = _col(vals)
    return out


class RefinedModel:
    """Container tying a quadtree grid to MARMITES inputs.

    Attributes
    ----------
    gridprops   : DISV gridprops (vertices, cell2d, ncpl, top, botm, nlay)
    active_nodes: icell2d of every active soil column, in cell order
    cells       : MARMITES cell list (cid, cid, 0, icell2d)
    geom        : VertexGeometry over the refined polygons
    inputs      : dict of (ncell, 1) column vectors
    """

    def __init__(self, gridprops, active_nodes, inputs, nlay):
        self.gridprops = gridprops
        self.active_nodes = np.asarray(active_nodes, dtype=int)
        self.cells = refined_cell_list(self.active_nodes)
        self.inputs = inputs
        self.nlay = int(nlay)
        self.geom = VertexGeometry.from_vertices(
            gridprops['vertices'], gridprops['cell2d'],
            self.active_nodes, nlay=self.nlay)

    @property
    def ncell(self):
        return len(self.cells)

    @property
    def ncpl(self):
        return int(self.gridprops['ncpl'])

    def outcrop_column(self, layer_of_cell):
        """(ncell, 1) outcrop layer (1-based), the MARMITES convention."""
        return _col(np.asarray(layer_of_cell, dtype=int) + 1)

    def summary(self):
        return ('RefinedModel: ncpl=%d, active cells=%d, nlay=%d, '
                'cell area %.4g..%.4g m2'
                % (self.ncpl, self.ncell, self.nlay,
                   self.geom.area.min(), self.geom.area.max()))
