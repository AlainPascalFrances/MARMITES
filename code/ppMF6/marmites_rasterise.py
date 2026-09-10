# -*- coding: utf-8 -*-
"""Rasterising display adapter for mesh results.  WP1c.7.

``MARMITESplot_v3.plotLAYER`` takes ``V`` shaped ``(ndays, nlay, nrow, ncol)``
and the post-processing carries ~100 more ``nrow``/``ncol`` references. None of
that is wrong -- it is a validated figure suite that draws real coordinate
axes, contour overlays, observation-point markers and a dozen map types.

So the mesh is converted TO A DISPLAY RASTER at the boundary, and the figure
suite is not touched at all. A per-cell vector becomes the ``(nrow, ncol)``
array it already expects; every figure then draws exactly as before.

Why the SOURCE grid is the display grid
---------------------------------------
The display raster defaults to the structured grid the model was projected
from -- La Mata's 50 m, 65 x 60. That is deliberate: it makes a mesh figure
directly comparable with the structured one, pixel for pixel, which is what
the WP1c.8 validation ladder needs. A finer raster can be asked for when the
mesh is finer than the source (the quadtree refines to 12.5 m), at the cost of
that comparability.

What this is NOT
----------------
Rasterising is for DISPLAY. It never feeds back into the model, and a figure
drawn this way shows a mesh cell as the block of pixels it covers, so a
Voronoi cell appears as a ragged blob rather than a polygon. True polygon maps
come free in the Streamlit app through flopy's ``PlotMapView`` and can be
added natively later; this adapter is what keeps the existing suite working in
the meantime.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import numpy as np

__all__ = ['DisplayRaster']


def _polygon_area(xy):
    a = np.asarray(xy, dtype=float)
    if a.shape[0] < 3:
        return 0.0
    x, y = a[:, 0], a[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _point_in_polygon(x, y, poly):
    inside = False
    x1, y1 = poly[-1]
    for x2, y2 in poly:
        if (y2 > y) != (y1 > y):
            if x < x2 + (y - y2) * (x1 - x2) / (y1 - y2):
                inside = not inside
        x1, y1 = x2, y2
    return inside


class DisplayRaster:
    """Maps mesh cells onto a regular raster, for plotting only.

    Attributes
    ----------
    nrow, ncol : int
    pix2cell : (nrow, ncol) int    icell2d of each pixel, -1 where none
    covered : (nrow, ncol) bool
    """

    def __init__(self, gridprops, x_edge, y_edge, refine=1):
        self.gridprops = gridprops
        self.ncpl = int(gridprops.get('ncpl', len(gridprops['cell2d'])))
        x_edge = np.asarray(x_edge, dtype=float)
        y_edge = np.asarray(y_edge, dtype=float)
        if int(refine) > 1:
            x_edge = _subdivide(x_edge, int(refine))
            y_edge = _subdivide(y_edge, int(refine))
        self.x_edge, self.y_edge = x_edge, y_edge
        self.ncol = x_edge.size - 1
        self.nrow = y_edge.size - 1
        self.xc = 0.5 * (x_edge[:-1] + x_edge[1:])
        self.yc = 0.5 * (y_edge[:-1] + y_edge[1:])
        self._vxy = {int(v[0]): (float(v[1]), float(v[2]))
                     for v in gridprops['vertices']}
        self.pix2cell = self._build()
        self.covered = self.pix2cell >= 0

    # ------------------------------------------------------------------ #
    @classmethod
    def from_projection(cls, proj, refine=0):
        """Display on the SOURCE structured grid of a MeshProjection.

        ``refine=0`` (the default) picks the factor automatically from the
        smallest mesh cell, so every cell gets at least one pixel centre. It
        matters: La Mata's quadtree refines to 12.5 m, and displayed on the
        50 m source grid **3828 of its 7728 cells fall between pixel centres
        and never appear in any figure at all**. ``refine=1`` forces the source
        resolution, which reproduces the structured picture exactly for a
        DIS-equivalent mesh.
        """
        if int(refine) <= 0:
            refine = cls.auto_refine(proj)
        return cls(proj.gridprops, proj._x_edge, proj._y_edge, refine=refine)

    # Beyond this the raster costs more than the figure is worth; a mesh that
    # needs it is better displayed as true polygons (the Streamlit app).
    MAX_REFINE = 8

    @staticmethod
    def auto_refine(proj):
        """Display refinement that gives the smallest mesh cell a pixel."""
        import numpy as _np
        areas = _np.array([_polygon_area([proj._vxy[int(iv)]
                                          for iv in rec[4:4 + int(rec[3])]])
                           for rec in proj.gridprops['cell2d']], dtype=float)
        areas = areas[areas > 0]
        if areas.size == 0:
            return 1
        src = min(float(_np.min(_np.diff(proj._x_edge))),
                  float(_np.min(-_np.diff(proj._y_edge))))
        need = src / float(_np.sqrt(areas.min()))
        return int(min(max(1, _np.ceil(need)), DisplayRaster.MAX_REFINE))

    def cell_polygon(self, ic):
        rec = self.gridprops['cell2d'][int(ic)]
        return [self._vxy[int(iv)] for iv in rec[4:4 + int(rec[3])]]

    def _build(self):
        """Which mesh cell each pixel centre falls in.

        Walks the MESH cells and fills the pixels inside each, rather than
        searching for every pixel: the work is then proportional to the number
        of pixels rather than to pixels x cells.
        """
        out = np.full((self.nrow, self.ncol), -1, dtype=int)
        # yc descends with the row index, so search the negated array
        neg_yc = -self.yc
        for ic in range(self.ncpl):
            poly = self.cell_polygon(ic)
            px = [p[0] for p in poly]
            py = [p[1] for p in poly]
            c0 = int(np.searchsorted(self.xc, min(px), 'left'))
            c1 = int(np.searchsorted(self.xc, max(px), 'right'))
            r0 = int(np.searchsorted(neg_yc, -max(py), 'left'))
            r1 = int(np.searchsorted(neg_yc, -min(py), 'right'))
            for r in range(max(r0, 0), min(r1, self.nrow)):
                y = self.yc[r]
                for c in range(max(c0, 0), min(c1, self.ncol)):
                    if out[r, c] < 0 and _point_in_polygon(self.xc[c], y, poly):
                        out[r, c] = ic
        return out

    # ------------------------------------------------------------------ #
    def field(self, values, nodata=np.nan):
        """A per-cell vector -> ``(nrow, ncol)``.

        ``values`` is indexed by icell2d and may be ``(ncpl,)`` or
        ``(ncpl, 1)``.
        """
        v = np.asarray(values, dtype=float).reshape(-1)
        if v.size != self.ncpl:
            raise ValueError('expected %d values (ncpl), got %d'
                             % (self.ncpl, v.size))
        idx = np.where(self.covered, self.pix2cell, 0)
        return np.where(self.covered, v[idx], nodata)

    def field_int(self, values, nodata=0):
        v = np.asarray(values).reshape(-1)
        idx = np.where(self.covered, self.pix2cell, 0)
        return np.where(self.covered, v[idx], nodata).astype(int)

    def stack(self, values, nodata=np.nan):
        """``(n, ncpl[, 1])`` -> ``(n, nrow, ncol)``; also handles
        ``(nt, nlay, ncpl)`` -> ``(nt, nlay, nrow, ncol)``."""
        a = np.asarray(values, dtype=float)
        if a.ndim >= 2 and a.shape[-1] == 1 and a.shape[-2] == self.ncpl:
            a = a[..., 0]
        if a.ndim == 1:
            return self.field(a, nodata)
        lead, out = a.shape[:-1], []
        for sub in a.reshape(-1, a.shape[-1]):
            out.append(self.field(sub, nodata))
        return np.asarray(out).reshape(lead + (self.nrow, self.ncol))

    def from_cells(self, values, cells, nodata=np.nan):
        """A per-MM-CELL vector (the active cells only) -> ``(nrow, ncol)``.

        ``cells`` is the MARMITES cell list, whose ``c[3]`` is the icell2d.
        """
        full = np.full(self.ncpl, np.nan, dtype=float)
        for v, c in zip(np.asarray(values, dtype=float).reshape(-1), cells):
            full[int(c[3])] = v
        out = self.field(full, nodata)
        return np.where(np.isnan(out), nodata, out)

    def coverage(self):
        """Fraction of pixels that landed on a mesh cell, and cells that got
        no pixel at all.

        A cell smaller than a display pixel can fall between pixel centres and
        vanish from every figure. On a coarse mesh that never happens; on a
        refined quadtree displayed at the source resolution it does, and the
        number is worth printing rather than discovering in a picture.
        """
        got = np.unique(self.pix2cell[self.covered])
        return {'pixels_covered': float(self.covered.mean()),
                'cells_shown': int(got.size),
                'cells_missing': int(self.ncpl - got.size)}


class MapAdapter:
    """Model-shaped arrays -> display-shaped arrays, for the native maps.

    On a structured grid every method is the identity, so the map code has one
    path and the regression anchor cannot drift. On a mesh the arrays are
    rasterised onto the display grid.

    The point of routing EVERY array through here is that a map mixes sources
    -- per-MM-cell result vectors, per-layer MF6 budgets, the ibound mask and
    the cell-area conversion -- and they must all end up on the same grid. One
    of them left in model shape is a broadcast error at best and a silently
    misaligned picture at worst.
    """

    def __init__(self, cMF, refine=0):
        self.cMF = cMF
        self.dr = None
        self.nrow, self.ncol = int(cMF.nrow), int(cMF.ncol)
        proj = getattr(cMF, 'mesh_proj', None)
        if proj is not None:
            self.dr = DisplayRaster.from_projection(proj, refine=refine)
            self.nrow, self.ncol = self.dr.nrow, self.dr.ncol

    @property
    def on_mesh(self):
        return self.dr is not None

    def report(self):
        if not self.on_mesh:
            return 'display: model grid %d x %d' % (self.nrow, self.ncol)
        c = self.dr.coverage()
        return ('display raster %d x %d (refine x%d): %d of %d mesh cells '
                'visible, %.1f%% of pixels covered'
                % (self.nrow, self.ncol,
                   (self.dr.x_edge.size - 1) // (self.cMF.mesh_proj.ncol_src),
                   c['cells_shown'], self.dr.ncpl + 0,
                   100.0 * c['pixels_covered']))

    def cells(self, values, cells, nodata=np.nan):
        """A per-MM-cell vector -> ``(nrow, ncol)``."""
        if not self.on_mesh:
            g = np.full((self.nrow, self.ncol), nodata, dtype=float)
            for v, c in zip(np.asarray(values, dtype=float).reshape(-1), cells):
                g[c[1], c[2]] = v
            return g
        return self.dr.from_cells(values, cells, nodata=nodata)

    def lay(self, arr, nodata=np.nan):
        """``(nlay, ...)`` model-shaped -> ``(nlay, nrow, ncol)``."""
        a = np.asarray(arr, dtype=float)
        nlay = a.shape[0]
        if not self.on_mesh:
            return a.reshape(nlay, self.nrow, self.ncol)
        return np.stack([self.dr.field(a[k].reshape(-1), nodata)
                         for k in range(nlay)])

    def lay_int(self, arr, nodata=0):
        a = np.asarray(arr)
        nlay = a.shape[0]
        if not self.on_mesh:
            return a.reshape(nlay, self.nrow, self.ncol)
        return np.stack([self.dr.field_int(a[k].reshape(-1), nodata)
                         for k in range(nlay)])

    def cell_area(self):
        """``(nrow, ncol)`` plan area of the MODEL cell each pixel belongs to.

        Not the pixel's own area: the m3/d -> mm/d conversion divides by the
        area of the cell the flux was computed in, and on a mesh that is the
        mesh cell, whatever the display resolution.
        """
        if not self.on_mesh:
            delr = np.asarray(self.cMF.delr, float)
            delc = np.asarray(self.cMF.delc, float)
            return delc[:, None] * delr[None, :]
        proj = self.cMF.mesh_proj
        areas = np.array([_polygon_area(self.dr.cell_polygon(i))
                          for i in range(self.dr.ncpl)], dtype=float)
        return self.dr.field(areas, nodata=np.nan)


def _subdivide(edge, k):
    """Split every interval of an edge array into ``k`` equal parts."""
    edge = np.asarray(edge, dtype=float)
    out = [edge[0]]
    for a, b in zip(edge[:-1], edge[1:]):
        for m in range(1, k + 1):
            out.append(a + (b - a) * m / k)
    return np.asarray(out, dtype=float)
