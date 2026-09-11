# -*- coding: utf-8 -*-
"""Grid-agnostic cell geometry for MARMITES (Phase 4).

MARMITES iterates a 1-D list of active soil columns (Phase 2). The only
spatial quantities the soil model and the coupler need per cell are:

  * ``area``  [L^2] -- cell plan area, used to convert MODFLOW volumetric
    fluxes (exfiltration m3/d, WEL m3/d) to/from MARMITES depths (mm/d);
  * ``width`` [L]   -- a characteristic horizontal cell dimension. It used to
    size the surface-storage geometry, which WP1d removed with the pond
    module; it is kept because a mesh cell still has no delr/delc.

This module supplies both for a structured (DIS) and a vertex (DISV) grid,
so the soil model never touches ``delr``/``delc`` again.

Width semantics
---------------
The legacy formulas use ``delr[j]`` (the cell's column width) as the
characteristic length. ``StructuredGeometry`` reproduces that exactly.
``VertexGeometry`` has no row/column notion, so it uses ``sqrt(area)`` --
identical for square cells (La Mata: 50 x 50 m) and the natural
generalization for irregular polygons. This is an approximation only for
non-square cells and only affects the ponding/open-water shape factor.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import numpy as np

__all__ = ['CellGeometry', 'StructuredGeometry', 'VertexGeometry',
           'disv_from_structured', 'polygon_area']


def polygon_area(xy):
    """Absolute polygon area by the shoelace formula (no shapely needed)."""
    a = np.asarray(xy, dtype=float)
    if a.shape[0] < 3:
        return 0.0
    x, y = a[:, 0], a[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


class CellGeometry:
    """Per-cell geometry indexed by MARMITES cell id (cid).

    Attributes
    ----------
    area, width : (ncell,) float arrays
    ncpl : int    number of map-view cells of the underlying grid
    nlay : int
    kind : 'structured' | 'vertex'
    """

    kind = 'abstract'

    def __init__(self, area, width, ncpl, nlay):
        self.area = np.asarray(area, dtype=float)
        self.width = np.asarray(width, dtype=float)
        self.ncpl = int(ncpl)
        self.nlay = int(nlay)
        if self.area.shape != self.width.shape:
            raise ValueError('area and width must have the same shape')

    @property
    def ncell(self):
        return int(self.area.size)

    def __repr__(self):
        return ('<%s ncell=%d ncpl=%d nlay=%d area[%.4g..%.4g]>'
                % (type(self).__name__, self.ncell, self.ncpl, self.nlay,
                   self.area.min() if self.ncell else 0,
                   self.area.max() if self.ncell else 0))


class StructuredGeometry(CellGeometry):
    """DIS geometry: area = delc[i]*delr[j], width = delr[j] (legacy)."""

    kind = 'structured'

    def __init__(self, delr, delc, i_arr, j_arr, nlay=1):
        delr = np.asarray(delr, dtype=float)
        delc = np.asarray(delc, dtype=float)
        i_arr = np.asarray(i_arr, dtype=int)
        j_arr = np.asarray(j_arr, dtype=int)
        area = delc[i_arr] * delr[j_arr]
        width = delr[j_arr]
        super().__init__(area, width, ncpl=delr.size * delc.size, nlay=nlay)
        self.delr, self.delc = delr, delc
        self.i_arr, self.j_arr = i_arr, j_arr

    @classmethod
    def from_cells(cls, cMF, cells):
        """Build from a clsMF configuration and the MARMITES cell list."""
        i_arr = np.array([c[1] for c in cells], dtype=int)
        j_arr = np.array([c[2] for c in cells], dtype=int)
        return cls(cMF.delr, cMF.delc, i_arr, j_arr, nlay=getattr(cMF, 'nlay', 1))


class VertexGeometry(CellGeometry):
    """DISV geometry: area from the cell polygon, width = sqrt(area)."""

    kind = 'vertex'

    def __init__(self, cell_areas, icell2d_arr, ncpl, nlay=1):
        cell_areas = np.asarray(cell_areas, dtype=float)      # (ncpl,)
        icell2d_arr = np.asarray(icell2d_arr, dtype=int)      # (ncell,)
        area = cell_areas[icell2d_arr]
        super().__init__(area, np.sqrt(area), ncpl=ncpl, nlay=nlay)
        self.icell2d_arr = icell2d_arr
        self.cell_areas = cell_areas

    @classmethod
    def from_vertices(cls, vertices, cell2d, icell2d_arr, nlay=1):
        """Build from MODFLOW DISV ``vertices``/``cell2d`` lists.

        vertices : [[iv, x, y], ...]
        cell2d   : [[icell2d, xc, yc, ncvert, iv1, iv2, ...], ...]
        """
        vxy = {int(v[0]): (float(v[1]), float(v[2])) for v in vertices}
        ncpl = len(cell2d)
        areas = np.zeros(ncpl, dtype=float)
        for rec in cell2d:
            ic = int(rec[0])
            ivs = [int(x) for x in rec[4:4 + int(rec[3])]]
            areas[ic] = polygon_area([vxy[iv] for iv in ivs])
        return cls(areas, icell2d_arr, ncpl=ncpl, nlay=nlay)

    @classmethod
    def from_modelgrid(cls, modelgrid, icell2d_arr, nlay=None):
        """Build from a flopy VertexGrid (uses get_cell_vertices)."""
        ncpl = int(np.atleast_1d(modelgrid.ncpl)[0])
        areas = np.array([polygon_area(modelgrid.get_cell_vertices(ic))
                          for ic in range(ncpl)], dtype=float)
        if nlay is None:
            nlay = int(getattr(modelgrid, 'nlay', 1))
        return cls(areas, icell2d_arr, ncpl=ncpl, nlay=nlay)


# --------------------------------------------------------------------- #
# DISV grid construction
# --------------------------------------------------------------------- #

def disv_from_structured(delr, delc, xorigin=0.0, yorigin=0.0):
    """Build DISV ``vertices``/``cell2d`` equivalent to a DIS grid.

    One cell2d per structured cell, with ``icell2d = i*ncol + j`` -- the same
    flat index the MARMITES cell list already stores as ``node``. This makes
    the DISV model geometrically identical to the DIS model, which is the
    Phase-4 correctness baseline (a DISV-from-DIS run must reproduce DIS).

    MODFLOW row 0 is the northernmost row, so y decreases with i; vertices
    are listed clockwise starting at the cell's north-west corner.

    Returns
    -------
    vertices : [[iv, x, y], ...]
    cell2d   : [[icell2d, xc, yc, 4, nw, ne, se, sw], ...]
    ncpl     : int
    """
    delr = np.asarray(delr, dtype=float)          # column widths (x)
    delc = np.asarray(delc, dtype=float)          # row heights (y)
    ncol, nrow = delr.size, delc.size

    x_edge = np.concatenate(([0.0], np.cumsum(delr))) + float(xorigin)
    y_top = float(yorigin) + float(np.sum(delc))
    y_edge = y_top - np.concatenate(([0.0], np.cumsum(delc)))   # y_edge[0] = north

    def iv(r, c):
        return r * (ncol + 1) + c

    vertices = [[iv(r, c), float(x_edge[c]), float(y_edge[r])]
                for r in range(nrow + 1) for c in range(ncol + 1)]

    cell2d = []
    for i in range(nrow):
        for j in range(ncol):
            ic = i * ncol + j
            xc = 0.5 * (x_edge[j] + x_edge[j + 1])
            yc = 0.5 * (y_edge[i] + y_edge[i + 1])
            cell2d.append([ic, float(xc), float(yc), 4,
                           iv(i, j), iv(i, j + 1), iv(i + 1, j + 1), iv(i + 1, j)])
    return vertices, cell2d, nrow * ncol


def geometry_for(cMF, cells, grid='dis', vertices=None, cell2d=None):
    """Convenience factory used by the drivers.

    grid='dis'  -> StructuredGeometry from cMF.delr/delc
    grid='disv' -> VertexGeometry; vertices/cell2d default to the
                   DIS-equivalent grid (disv_from_structured).
    """
    if grid == 'dis':
        return StructuredGeometry.from_cells(cMF, cells)
    if grid != 'disv':
        raise ValueError("grid must be 'dis' or 'disv' (DISU is out of scope)")
    if vertices is None or cell2d is None:
        vertices, cell2d, _ = disv_from_structured(
            cMF.delr, cMF.delc, getattr(cMF, 'xllcorner', 0.0),
            getattr(cMF, 'yllcorner', 0.0))
    icell2d_arr = np.array([c[3] for c in cells], dtype=int)   # node == i*ncol+j
    return VertexGeometry.from_vertices(vertices, cell2d, icell2d_arr,
                                        nlay=getattr(cMF, 'nlay', 1))
