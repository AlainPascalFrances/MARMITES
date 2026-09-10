# -*- coding: utf-8 -*-
"""Project a structured MARMITES model onto an unstructured (DISV) mesh.

WP1c.1 -- the missing link between the mesh producers (quadtree today,
Voronoi in WP1c.2) and everything downstream: MMsoil, the coupler, the MF6
build and the post-processing.

The (ncpl, 1) convention
------------------------
A refined cell has no (row, column), and MARMITES' per-cell kernel indexes its
inputs as ``grid[i, j]``. Phase 4 solved that with a degenerate layout
(``marmites_gridgen.refined_cell_list``): cells ``(cid, cid, 0, icell2d)`` over
``(ncell, 1)`` column vectors, compacted to the ACTIVE cells.

This module uses the other degenerate layout, and it is the one that makes the
whole path collapse to a reshape:

    nrow := ncpl,  ncol := 1        arrays are (ncpl, 1) / (nlay, ncpl, 1)

Because ``i * ncol + j == i * 1 + 0 == i``, the row index **is** the icell2d.
Three things then work with no change at all:

  * ``MMsoil.build_cell_list`` walks nrow x ncol and emits
    ``(cid, icell2d, 0, icell2d)`` -- the correct node, for free;
  * ``clsMF6._cellid`` returns ``(k, i * ncol + j) == (k, icell2d)`` without
    needing its ``cell_nodes`` map;
  * ``clsMF6._griddata`` reshapes ``(nlay, nrow, ncol) -> (nlay, ncpl)``, which
    is exactly the MF6 DISV griddata layout, because ``nrow * ncol == ncpl``.

The compacted ``(ncell, 1)`` convention has none of those properties: it needs
an explicit node map and breaks the ``_griddata`` reshape. Both are correct;
this one is cheaper, so the driver uses it. ``RefinedModel`` keeps the older
convention for the standalone quadtree script and its tests.

Inactive mesh cells are carried through as full rows with ``ibound = 0``, so
they simply never enter the cell list -- the same way an inactive structured
cell behaves today.

Sampling  (WP1c.3)
------------------
``how='auto'`` is the default and the one to use: **area-weighted means** for
continuous fields, **largest-overlap class** for zone rasters. The weights are
EXACT overlap areas, from clipping each mesh cell polygon against each raster
cell it touches (Sutherland-Hodgman, no shapely).

Why it matters: a 100 m Voronoi cell covers about four 50 m raster cells, and
centre sampling keeps one of them and discards the rest. On La Mata's soil
thickness that is the difference between a resampled mean that tracks the
raster and one that drifts with wherever the centres happened to land.

Measured on the 100 m La Mata Voronoi mesh, error in the active-domain mean
against the raster: soil thickness 2.93% -> 0.06%, stream depth 14.47% ->
0.24%, stream width 14.55% -> 0.42%.

For ZONE rasters the comparison goes the other way and it is worth stating.
Majority vote is better PER CELL -- each cell gets the class that actually
dominates it -- but biased IN AGGREGATE, because a minority class that never
wins anywhere is absorbed by the dominant one. On La Mata's soil zones,
majority misassigns 7.64% of the active area against centre sampling's 3.66%,
and zone 1 drops from 13.1% of the catchment to 9.5%. Neither is simply
right; ``[grid] resample = "centre"`` is there for the case where preserving
zone proportions matters more than per-cell fidelity.

``how='centre'`` is also exact, and cheaper, whenever the mesh only REFINES
the raster (the quadtree case: every mesh cell lies inside one raster cell).

Both reduce to the same thing on a DIS-equivalent mesh, which is what keeps
the identity test valid under either mode.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import copy

import numpy as np

__all__ = ['MeshProjection', 'project_model', 'MeshProjectionError',
           'ZONE_FIELDS', 'SAMPLING_MODES', 'MODEL_SAMPLING_MODES']

# 'auto' is the one to use: it picks 'majority' for a zone raster and 'area'
# for a continuous field, so the rule follows the FIELD rather than the call
# site. 'centre' is kept because it is what WP1c.1 ran and therefore the
# comparison point for judging what the upgrade changed.
SAMPLING_MODES = ('auto', 'centre', 'area', 'majority')

# ...but only these two are coherent as a WHOLE-MODEL policy, and they are what
# `[grid] resample` offers. 'area' forced on a zone raster would interpolate
# class codes; 'majority' forced on an elevation would round it to whole
# metres. Those two are per-FIELD rules, chosen by 'auto'.
MODEL_SAMPLING_MODES = ('auto', 'centre')

# Field-by-field sampling rules. 'zone' means an integer class code, which may
# never be interpolated or averaged -- a soil zone halfway between 1 and 3 is
# not 2. Continuous fields may be averaged (WP1c.3).
ZONE_FIELDS = ('gridSOIL', 'gridMETEO', 'gridIRR', 'ibound', 'iuzfbnd',
               'outcropL')


class MeshProjectionError(Exception):
    """Raised when a model cannot be projected onto a mesh."""


def _shoelace(pts):
    """Absolute polygon area, no shapely needed."""
    a = np.asarray(pts, dtype=float)
    if a.shape[0] < 3:
        return 0.0
    x, y = a[:, 0], a[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _clip_to_rect(pts, xmin, ymin, xmax, ymax):
    """Sutherland-Hodgman clip of a polygon against an axis-aligned rectangle.

    Exact, and exact is what is wanted here: the clip window is a raster cell
    and the subject is a mesh cell, so the clipped area IS the overlap area
    used to weight the resampling. The algorithm requires a CONVEX clip
    window, which a rectangle is; every mesh this code produces (Voronoi
    cells, quadtree cells, DIS-equivalent rectangles) is convex too.
    """
    out = [(float(x), float(y)) for x, y in pts]
    for inside, cut in (
            (lambda p: p[0] >= xmin, lambda a, b: _cx(a, b, xmin, 0)),
            (lambda p: p[0] <= xmax, lambda a, b: _cx(a, b, xmax, 0)),
            (lambda p: p[1] >= ymin, lambda a, b: _cx(a, b, ymin, 1)),
            (lambda p: p[1] <= ymax, lambda a, b: _cx(a, b, ymax, 1))):
        if not out:
            return []
        buf, prev = [], out[-1]
        for cur in out:
            ci, pi = inside(cur), inside(prev)
            if ci:
                if not pi:
                    buf.append(cut(prev, cur))
                buf.append(cur)
            elif pi:
                buf.append(cut(prev, cur))
            prev = cur
        out = buf
    return out


def _point_in_polygon(x, y, poly):
    """Ray-casting point-in-polygon test (crossing number)."""
    inside = False
    n = len(poly)
    x1, y1 = poly[-1]
    for k in range(n):
        x2, y2 = poly[k]
        if (y2 > y) != (y1 > y):
            xin = x2 + (y - y2) * (x1 - x2) / (y1 - y2)
            if x < xin:
                inside = not inside
        x1, y1 = x2, y2
    return inside


def _cx(a, b, value, axis):
    """Intersection of segment a-b with the line coord[axis] == value."""
    d = b[axis] - a[axis]
    t = 0.0 if d == 0.0 else (value - a[axis]) / d
    other = 1 - axis
    p = [0.0, 0.0]
    p[axis] = value
    p[other] = a[other] + t * (b[other] - a[other])
    return (p[0], p[1])


def _mask_sentinel(arr, hnoflo, atol=0.09):
    """Mask ``hnoflo`` values WITHOUT disturbing the data under the mask.

    ``np.ma.masked_values(x, v)`` cannot be used here: it starts with
    ``filled(x, v)``, so every element already masked has its DATA replaced by
    ``v``. The driver hands MODFLOW ``np.asarray(cMF.top)``, which drops the
    mask and keeps the data, so that substitution is not cosmetic -- it
    silently changes the values MF6 receives at inactive cells. (``cMF.top``
    itself keeps real elevations under its mask, because numpy's masked
    subtraction copies the first operand's data there.)
    """
    data = np.ma.getdata(arr)
    m = np.ma.getmaskarray(arr) | np.isclose(data, hnoflo, rtol=0.0, atol=atol)
    return np.ma.array(data, mask=m)


def _edges(delr, delc, xllcorner, yllcorner):
    """Cell-edge coordinates of a structured grid.

    Returns ``(x_edge, y_edge)`` where ``x_edge`` increases west->east and
    ``y_edge`` DECREASES with the row index, because MODFLOW row 0 is the
    northernmost row.
    """
    delr = np.asarray(delr, dtype=float)
    delc = np.asarray(delc, dtype=float)
    x_edge = float(xllcorner) + np.concatenate(([0.0], np.cumsum(delr)))
    y_top = float(yllcorner) + float(np.sum(delc))
    y_edge = y_top - np.concatenate(([0.0], np.cumsum(delc)))
    return x_edge, y_edge


class MeshProjection:
    """Maps mesh cells onto the cells of a structured MARMITES model.

    Parameters
    ----------
    gridprops : dict
        DISV gridprops (``vertices``, ``cell2d``, ``ncpl``, ...).
    delr, delc : array-like
        Source structured grid spacings.
    xllcorner, yllcorner : float
        Source grid origin (lower-left), model CRS.

    Attributes
    ----------
    ncpl : int
    xy : (ncpl, 2)      mesh cell-centre coordinates
    row, col : (ncpl,)  source cell containing each mesh centre (clipped)
    inside : (ncpl,)    False where the centre falls outside the source grid
    """

    def __init__(self, gridprops, delr, delc, xllcorner, yllcorner):
        cell2d = gridprops['cell2d']
        self.gridprops = gridprops
        self.ncpl = int(gridprops.get('ncpl', len(cell2d)))
        if self.ncpl != len(cell2d):
            raise MeshProjectionError('ncpl=%d but %d cell2d records'
                                      % (self.ncpl, len(cell2d)))
        self.xy = np.array([[float(r[1]), float(r[2])] for r in cell2d],
                           dtype=float)
        self._vxy = {int(v[0]): (float(v[1]), float(v[2]))
                     for v in gridprops['vertices']}
        self._overlaps = None
        self.nrow_src = int(len(np.asarray(delc, dtype=float)))
        self.ncol_src = int(len(np.asarray(delr, dtype=float)))
        self._x_edge, self._y_edge = _edges(delr, delc, xllcorner, yllcorner)
        self.row, self.col, self.inside = self._locate(self.xy)

    # -- location ---------------------------------------------------------
    def _locate(self, xy):
        """Source (row, col) of each point; row/col clipped, `inside` is not.

        Uses the edge arrays rather than a division, so a non-uniform source
        grid works too -- La Mata is uniform 50 m but CdL need not be.
        """
        x = np.asarray(xy[:, 0], dtype=float)
        y = np.asarray(xy[:, 1], dtype=float)
        col = np.searchsorted(self._x_edge, x, side='right') - 1
        # y_edge decreases with the row index, so search the negated (thus
        # increasing) array.
        row = np.searchsorted(-self._y_edge, -y, side='right') - 1
        # The domain is CLOSED at its far edges. flopy's VoronoiGrid reports a
        # cell's "centre" as its GENERATOR, and for a boundary cell that point
        # lies exactly ON the domain boundary -- so a strict `x < xmax` test
        # declares a large share of a clipped Voronoi mesh to be outside the
        # model and fills those cells with nodata.
        xmin, xmax = self._x_edge[0], self._x_edge[-1]
        ymin, ymax = self._y_edge[-1], self._y_edge[0]
        eps = 1e-6 * max(xmax - xmin, ymax - ymin)
        inside = ((x >= xmin - eps) & (x <= xmax + eps)
                  & (y >= ymin - eps) & (y <= ymax + eps))
        return (np.clip(row, 0, self.nrow_src - 1),
                np.clip(col, 0, self.ncol_src - 1), inside)

    @staticmethod
    def _resolve_how(how, dtype):
        """'auto' picks the right rule from the FIELD KIND, not the caller.

        A zone raster carries class codes, so it can only ever be resolved by
        majority; a continuous field can only ever be area-averaged. Leaving
        that to each call site is how one field quietly gets interpolated.
        """
        if how == 'auto':
            return 'majority' if dtype is int else 'area'
        if how not in SAMPLING_MODES:
            raise MeshProjectionError(
                'unknown sampling mode %r (known: %s)'
                % (how, ', '.join(SAMPLING_MODES)))
        return how

    def locate(self, x, y):
        """Public point lookup: returns (row, col, inside) for scalars/arrays."""
        xy = np.column_stack([np.atleast_1d(x), np.atleast_1d(y)])
        return self._locate(xy)

    def _bboxes(self):
        if getattr(self, '_bbox', None) is None:
            bb = np.empty((self.ncpl, 4), dtype=float)
            for ic in range(self.ncpl):
                p = np.asarray(self.cell_polygon(ic), dtype=float)
                bb[ic] = (p[:, 0].min(), p[:, 1].min(),
                          p[:, 0].max(), p[:, 1].max())
            self._bbox = bb
        return self._bbox

    def cell_containing(self, x, y):
        """icell2d of the cell whose POLYGON contains (x, y).

        Falls back to the nearest centroid when no polygon contains the point
        -- which happens for a point exactly on a shared edge, and for one just
        outside the mesh boundary.

        Point-in-polygon rather than nearest-centroid throughout: the two agree
        on a true Voronoi mesh by construction, but NOT on a quadtree, where a
        small cell's centroid can be nearer to a point than the centroid of the
        large cell the point actually sits in.

        Ties are broken the way the legacy grid arithmetic breaks them. La
        Mata's observation point I1 sits exactly on a cell CORNER, where four
        cells contain it equally; ``cPROCESS.inputObs`` computes
        ``j = ceil(dx/cs) - 1`` and ``i = nrow - ceil(dy/cs)``, which puts such
        a point in the cell to the LEFT of a vertical edge and BELOW a
        horizontal one. Testing the point nudged by a micron in both negative
        directions reproduces that on any mesh -- so a DIS-equivalent mesh run
        places every observation in the same cell as the structured run, and
        the WP1c.8 comparison is not confounded by an arbitrary tie.
        """
        x, y = float(x), float(y)
        bb = self._bboxes()
        cand = np.flatnonzero((bb[:, 0] <= x) & (x <= bb[:, 2])
                              & (bb[:, 1] <= y) & (y <= bb[:, 3]))
        eps = 1e-6
        for probe in ((x - eps, y - eps), (x, y)):
            for ic in cand:
                if _point_in_polygon(probe[0], probe[1],
                                     self.cell_polygon(int(ic))):
                    return int(ic)
        return self.cell_of(x, y)

    def cell_of(self, x, y):
        """icell2d of the mesh cell whose CENTRE is nearest to (x, y).

        Used to move point features (drains, wells, observation points) onto
        the mesh. Nearest-centroid rather than point-in-polygon: it needs no
        geometry library, and for a Voronoi mesh the two are the same thing by
        construction -- a Voronoi cell IS the set of points nearest its seed.
        """
        d = ((self.xy[:, 0] - float(x)) ** 2
             + (self.xy[:, 1] - float(y)) ** 2)
        return int(np.argmin(d))

    # -- overlap weights (WP1c.3) -----------------------------------------
    def cell_polygon(self, ic):
        """Vertices of mesh cell ``ic`` as ``[(x, y), ...]``."""
        rec = self.gridprops['cell2d'][int(ic)]
        return [self._vxy[int(iv)] for iv in rec[4:4 + int(rec[3])]]

    def overlaps(self):
        """Per mesh cell, the source cells it covers and by how much.

        Returns a list of ``(rows, cols, weights)``; ``weights`` are overlap
        AREAS in m2. Computed once and cached, because it costs the same for
        every field and there are a dozen fields.
        """
        if self._overlaps is not None:
            return self._overlaps
        xe, ye = self._x_edge, self._y_edge
        out = []
        for ic in range(self.ncpl):
            poly = self.cell_polygon(ic)
            px = [p[0] for p in poly]
            py = [p[1] for p in poly]
            # candidate source cells: the raster window covering this polygon
            c0 = int(np.clip(np.searchsorted(xe, min(px), 'right') - 1,
                             0, self.ncol_src - 1))
            c1 = int(np.clip(np.searchsorted(xe, max(px), 'left'),
                             1, self.ncol_src))
            r0 = int(np.clip(np.searchsorted(-ye, -max(py), 'right') - 1,
                             0, self.nrow_src - 1))
            r1 = int(np.clip(np.searchsorted(-ye, -min(py), 'left'),
                             1, self.nrow_src))
            rr, cc, ww = [], [], []
            for r in range(r0, r1):
                ytop, ybot = ye[r], ye[r + 1]
                for c in range(c0, c1):
                    a = _shoelace(_clip_to_rect(poly, xe[c], ybot,
                                                xe[c + 1], ytop))
                    if a > 0.0:
                        rr.append(r)
                        cc.append(c)
                        ww.append(a)
            out.append((np.asarray(rr, dtype=int), np.asarray(cc, dtype=int),
                        np.asarray(ww, dtype=float)))
        self._overlaps = out
        return out

    def overlap_report(self):
        """How much of each mesh cell was actually covered by the source grid.

        A mean coverage below 1 means mesh cells hang over the edge of the
        raster; that is expected at the boundary and suspicious in bulk.
        """
        areas = np.array([_shoelace(self.cell_polygon(ic))
                          for ic in range(self.ncpl)], dtype=float)
        cov = np.array([w.sum() for _r, _c, w in self.overlaps()], dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            frac = np.where(areas > 0, cov / areas, 0.0)
        # How much of the SOURCE RECTANGLE the mesh actually tiles. Voronoi
        # cells do not close perfectly against a clipped boundary polygon, so
        # a mesh domain is typically a fraction of a percent smaller than the
        # structured one -- small, but it is a real difference in model area
        # and it is the reason a mass total over a field with a large nodata
        # sentinel will not balance to round-off.
        src_total = float((self._x_edge[-1] - self._x_edge[0])
                          * (self._y_edge[0] - self._y_edge[-1]))
        return {'coverage_mean': float(frac.mean()),
                'coverage_min': float(frac.min()),
                'domain_ratio': (float(areas.sum() / src_total)
                                 if src_total > 0 else float('nan')),
                'src_per_cell_mean': float(np.mean(
                    [len(w) for _r, _c, w in self.overlaps()]))}

    # -- sampling ---------------------------------------------------------
    def _sample_weighted(self, arr, fill, dtype, majority, valid=None):
        """Area-weighted mean, or largest-overlap class for a zone raster.

        ``valid`` marks the source cells that are PART OF THE MODEL. It must
        be passed for any field that carries a nodata sentinel outside the
        active domain, which is nearly all of them: La Mata's zone rasters
        hold hnoflo (9999.999) outside the catchment, and a boundary mesh cell
        overlapping four of those would otherwise elect 9999 as its soil zone
        by majority -- observed, and it moved the mean soil zone from 1.9 to
        1248. Averaging the same sentinel into an elevation gives a cell a
        mountain. Sniffing for the sentinel by value is not enough, because
        `fill` is not always the sentinel; the caller knows which cells are
        real, so it says so.

        Where a mesh cell overlaps nothing valid, the CENTRE value is kept and
        the result masked. That is what preserves the DIS-equivalent identity:
        a cell whose only source cell is masked still carries that cell's
        data, exactly as `how='centre'` would give.
        """
        data = np.ma.getdata(arr) if np.ma.isMaskedArray(arr) else np.asarray(arr)
        data = np.asarray(data, dtype=float)
        bad = np.zeros(data.shape, dtype=bool)
        if np.ma.isMaskedArray(arr):
            bad |= np.ma.getmaskarray(arr)
        bad |= ~np.isfinite(data)
        if valid is not None:
            bad |= ~np.asarray(valid, dtype=bool)
        vals = np.empty(self.ncpl, dtype=float)
        msk = np.zeros(self.ncpl, dtype=bool)
        centre = data[self.row, self.col]
        for ic, (rr, cc, ww) in enumerate(self.overlaps()):
            if rr.size:
                ok = ~bad[rr, cc]
                w, v = ww[ok], data[rr[ok], cc[ok]]
            else:
                w = np.empty(0)
                v = np.empty(0)
            if w.size == 0 or w.sum() <= 0.0:
                # Nothing valid underneath. If the cell still sits on the grid,
                # keep the centre datum (masked) so the DIS-equivalent identity
                # holds; if it is off the grid entirely, there is no datum.
                vals[ic] = centre[ic] if self.inside[ic] else fill
                msk[ic] = True
                continue
            if majority:
                classes, inv = np.unique(np.rint(v).astype(int),
                                         return_inverse=True)
                totals = np.bincount(inv, weights=w)
                vals[ic] = float(classes[int(np.argmax(totals))])
            else:
                vals[ic] = float(np.dot(w, v) / w.sum())
        # NOTE: `inside` is deliberately NOT applied to a cell that HAS valid
        # overlaps. A cell's value comes from what it covers, not from where
        # its nominal centre landed -- and a Voronoi centre is a generator on
        # the domain boundary often enough for the difference to matter.
        out = (np.rint(vals).astype(int) if dtype is int
               else vals.astype(float)).reshape(-1, 1)
        if np.ma.isMaskedArray(arr):
            return np.ma.array(out, mask=msk.reshape(-1, 1))
        return out

    def sample2d(self, arr, fill=np.nan, dtype=float, how='centre', valid=None):
        """Sample a source ``(nrow, ncol)`` array onto ``(ncpl, 1)``.

        ``fill`` is used where the mesh centre falls outside the source grid.

        A masked source gives a masked result, and the DATA UNDER THE MASK is
        carried through rather than replaced by ``fill``. That is not
        cosmetic: the driver hands MF6 ``np.asarray(cMF.top)``, which drops the
        mask, so overwriting the masked data would quietly change the values
        MODFLOW sees at inactive cells.
        """
        how = self._resolve_how(how, dtype)
        a = np.ma.getdata(arr) if np.ma.isMaskedArray(arr) else np.asarray(arr)
        if a.shape != (self.nrow_src, self.ncol_src):
            raise MeshProjectionError(
                'expected a (%d, %d) source array, got %s'
                % (self.nrow_src, self.ncol_src, (a.shape,)))
        if how in ('area', 'majority'):
            return self._sample_weighted(arr, fill, dtype,
                                         majority=(how == 'majority'),
                                         valid=valid)
        vals = np.asarray(a, dtype=float)[self.row, self.col]
        vals = np.where(self.inside, vals, fill)
        msk = None
        if np.ma.isMaskedArray(arr):
            # nearest-sample the mask itself, and mask anything off the grid:
            # a mesh cell with no source cell has no datum, masked or not.
            msk = np.ma.getmaskarray(arr)[self.row, self.col] | ~self.inside
        if dtype is int:
            bad = ~np.isfinite(vals)
            if bad.any():
                vals = np.where(bad, (fill if np.isfinite(fill) else 0), vals)
            out = np.rint(vals).astype(int).reshape(-1, 1)
        else:
            out = vals.astype(float).reshape(-1, 1)
        if msk is not None:
            return np.ma.array(out, mask=msk.reshape(-1, 1))
        return out

    def sample3d(self, arr, fill=np.nan, dtype=float, how='centre', valid=None):
        """Sample a source ``(nlay, nrow, ncol)`` array onto ``(nlay, ncpl, 1)``."""
        a = np.asarray(arr)
        if a.ndim != 3:
            raise MeshProjectionError('expected a 3-D array, got %s' % (a.shape,))
        masked = np.ma.isMaskedArray(arr)
        mask = np.ma.getmaskarray(arr) if masked else None
        layers = []
        for k in range(a.shape[0]):
            src = np.ma.array(a[k], mask=mask[k]) if masked else a[k]
            v = None if valid is None else (
                valid[k] if np.ndim(valid) == 3 else valid)
            layers.append(self.sample2d(src, fill=fill, dtype=dtype, how=how,
                                        valid=v))
        return np.stack(layers)

    def sample_layer_property(self, val, fill=np.nan, how='centre', valid=None):
        """Project a cMF layer property, which may be a list of scalars OR of
        arrays (``hk_actual`` is arrays, ``vka_actual`` is floats).

        Scalars are left alone: they are already grid-independent, and
        expanding them here would only make the MF6 build carry a full array
        for a constant.
        """
        if not isinstance(val, (list, tuple)):
            return val
        out = []
        for k, v in enumerate(val):
            a = np.asarray(v)
            if a.ndim == 2:
                vk = None if valid is None else (
                    valid[k] if np.ndim(valid) == 3 else valid)
                out.append(self.sample2d(v, fill=fill, how=how, valid=vk))
            else:
                out.append(v)
        return out

    def remap_drn_records(self, recs, botm_src, botm_mesh, warn=None):
        """Move DRN records onto the mesh, KEEPING THEIR HEIGHT ABOVE THE CELL
        BOTTOM rather than their absolute elevation.

        La Mata's twelve drains all sit exactly 1.510 m above their own cell's
        bottom: the elevation is defined relative to the layer, not in absolute
        metres. Carrying the absolute value across is therefore wrong twice
        over -- it loses the intent, and because a mesh cell's bottom differs
        from the source cell's, MODFLOW rejects the package outright with
        "DRN BOUNDARY ELEVATION IS LESS THAN CELL BOTTOM".

        Conductance is carried unchanged. It is a property of the drain, not of
        the cell; scaling it by an area ratio would be a modelling decision,
        and this function's job is transfer, not calibration.
        """
        bs = np.asarray(botm_src, dtype=float)
        bm = np.asarray(botm_mesh, dtype=float)
        out = self.remap_records(recs, warn=warn)
        offsets = []
        for src, dst in zip(recs, out):
            lay, i, j = int(src[0]), int(src[1]), int(src[2])
            ic = int(dst[1])
            off = float(src[3]) - float(bs[lay, i, j])
            dst[3] = float(bm[lay, ic, 0]) + off
            offsets.append(off)
        if warn and offsets:
            warn('%d drain(s) re-anchored to the mesh cell bottom, keeping '
                 'their height above it (%.3f..%.3f m).'
                 % (len(offsets), min(offsets), max(offsets)))
        return out

    def remap_records(self, recs, warn=None):
        """Move ``(layer, row, col, *rest)`` boundary records onto the mesh.

        DRN and GHB records carry a source (row, col). On the mesh they become
        ``(layer, icell2d, 0, *rest)`` so that the legacy consumers -- which
        read ``rec[1], rec[2]`` as (i, j) -- keep working under the (ncpl, 1)
        convention.

        A coarser mesh can put two source records in one cell. That is a real
        change to the boundary condition, not a rounding detail, so it is
        reported rather than silently summed.
        """
        out, seen, clashes = [], {}, []
        x_c = 0.5 * (self._x_edge[:-1] + self._x_edge[1:])
        y_c = 0.5 * (self._y_edge[:-1] + self._y_edge[1:])
        for rec in recs:
            lay, i, j = int(rec[0]), int(rec[1]), int(rec[2])
            ic = self.cell_of(x_c[j], y_c[i])
            key = (lay, ic)
            if key in seen:
                clashes.append((seen[key], (lay, i, j), ic))
            seen[key] = (lay, i, j)
            out.append([lay, ic, 0] + list(rec[3:]))
        if clashes and warn is not None:
            warn('%d boundary record(s) share a mesh cell with another: %s'
                 % (len(clashes),
                    '; '.join('(L%d r%d c%d) and (L%d r%d c%d) -> icell2d %d'
                              % (a[0], a[1], a[2], b[0], b[1], b[2], ic)
                              for a, b, ic in clashes[:5])))
        return out


# --------------------------------------------------------------------- #

class _MeshProcess:
    """cPROCESS proxy reporting the mesh shape.

    ``float2array`` expands per-layer scalars using ``self.nrow``/``self.ncol``,
    so it has to agree with the projected arrays or MMsoil's Sy lookup silently
    indexes the wrong shape.
    """

    def __init__(self, src, ncpl, nlay):
        # set through __dict__: __getattr__ below delegates everything else
        self.__dict__['_src'] = src
        self.__dict__['nrow'] = int(ncpl)
        self.__dict__['ncol'] = 1
        self.__dict__['nlay'] = int(nlay)

    def __getattr__(self, name):
        return getattr(self.__dict__['_src'], name)

    def float2array(self, array):
        if self.nlay < 2:
            if isinstance(array, list):
                return np.asarray(array).reshape((1, self.nrow, self.ncol))
            return np.asarray([array])
        a = np.asarray(array)
        if a.ndim == 1 and a.shape[0] == self.nlay:
            out = np.ones((self.nlay, self.nrow, self.ncol), dtype=float)
            for k, e in enumerate(a):
                out[k] *= e
            return out
        return a


def project_model(cMF, gridprops, grids, warn=None, how='auto'):
    """Re-express a structured MARMITES model on a DISV mesh.

    Parameters
    ----------
    cMF : parsed clsMF, already fully set up on the structured grid (top/botm
        adjusted for the soil column, outcropL computed).
    gridprops : dict     DISV gridprops from a mesh producer.
    grids : dict         the spatial inputs handed to ``build_context``
        (``gridSOIL``, ``gridSOILthick``, ``gridVEGarea``, ...).
    warn : callable      one-argument message sink (default: print).

    Returns
    -------
    (mesh_cMF, grids_out, proj)
        ``mesh_cMF`` is a shallow copy of ``cMF`` with nrow=ncpl, ncol=1 and
        every spatial array projected; ``grids_out`` mirrors ``grids``;
        ``proj`` is the :class:`MeshProjection` used, kept so that later steps
        (observation points, SFR, LAK) can move their own features.

    The source ``cMF`` is NOT modified: a failed projection must leave the
    structured model runnable.
    """
    if warn is None:
        def warn(msg):
            print('WARNING: %s' % msg)
    if how not in MODEL_SAMPLING_MODES:
        raise MeshProjectionError(
            'resampling policy %r is not a whole-model policy (use %s). '
            "'area' and 'majority' are per-field rules that 'auto' selects "
            'from the field kind.'
            % (how, ' or '.join(repr(x) for x in MODEL_SAMPLING_MODES)))
    hnoflo = float(getattr(cMF, 'hnoflo', 9999.999))
    proj = MeshProjection(gridprops, cMF.delr, cMF.delc,
                          getattr(cMF, 'xllcorner', 0.0),
                          getattr(cMF, 'yllcorner', 0.0))
    ncpl, nlay = proj.ncpl, int(cMF.nlay)

    m = copy.copy(cMF)
    m.nrow, m.ncol = ncpl, 1
    # One unit cell each way: the mesh's real geometry lives in `geom`
    # (VertexGeometry), and anything still reading delr/delc on a mesh is a
    # bug we want to see rather than a plausible wrong number.
    m.delr = [1.0]
    m.delc = [1.0] * ncpl
    m.cPROCESS = _MeshProcess(cMF.cPROCESS, ncpl, nlay)

    # ---- structural arrays. ibound decides what is active, so it is
    # projected first and everything else is made consistent with it.
    # ibound and iuzfbnd carry their own validity (a 0 IS the answer), so they
    # are resampled over every source cell. Everything else is resampled ONLY
    # from cells inside the model, because outside them the arrays hold hnoflo.
    m.ibound = proj.sample3d(cMF.ibound, fill=0, dtype=int, how=how)
    m.iuzfbnd = proj.sample2d(cMF.iuzfbnd, fill=0, dtype=int, how=how)
    src_active3 = np.abs(np.asarray(cMF.ibound)) != 0
    src_active = np.asarray(cMF.outcropL) > 0

    # outcropL is RECOMPUTED, not sampled: it is a derived quantity (the first
    # active layer), and sampling it independently of ibound is how a cell ends
    # up claiming to outcrop in a layer the projection made inactive.
    outcrop = np.zeros((ncpl, 1), dtype=int)
    for L in range(nlay):
        ib = np.abs(m.ibound[L]) != 0
        outcrop += ((outcrop == 0) & ib) * (L + 1)
    m.outcropL = outcrop

    # ---- elevations. sample2d carries the source mask AND the data under it;
    # _mask_sentinel additionally catches an unmasked hnoflo, for a cMF that
    # arrived as plain arrays, without overwriting anything.
    m.elev = _mask_sentinel(
        proj.sample2d(cMF.elev, fill=hnoflo, how=how, valid=src_active), hnoflo)
    m.top = _mask_sentinel(
        proj.sample2d(cMF.top, fill=hnoflo, how=how, valid=src_active), hnoflo)
    m.botm = proj.sample3d(np.asarray(cMF.botm), fill=hnoflo, how=how,
                           valid=src_active3)
    if getattr(cMF, 'strt', None) is not None:
        m.strt = proj.sample3d(np.asarray(cMF.strt), fill=hnoflo, how=how,
                               valid=src_active3)

    # ---- aquifer properties
    for name in ('hk_actual', 'vka_actual', 'ss_actual', 'sy_actual',
                 'vks_actual', 'thick'):
        val = getattr(cMF, name, None)
        if val is None:
            continue
        if isinstance(val, list):
            setattr(m, name, proj.sample_layer_property(
                val, fill=hnoflo, how=how, valid=src_active3))
            continue
        a = np.asarray(val)
        if a.ndim == 3:
            setattr(m, name, proj.sample3d(val, fill=hnoflo, how=how,
                                           valid=src_active3))
        elif a.ndim == 2:
            setattr(m, name, proj.sample2d(val, fill=hnoflo, how=how,
                                           valid=src_active))

    # ---- boundary records. DRN elevations are re-anchored to the receiving
    # cell's bottom (see remap_drn_records); GHB carries an absolute HEAD, a
    # boundary condition on the water table itself, so it is moved unchanged.
    _src_botm = np.asarray(cMF.botm, dtype=float)

    def _drn(v):
        return proj.remap_drn_records(v, _src_botm, m.botm, warn=warn)

    for name, fn in (('layer_row_column_elevation_cond', _drn),
                     ('layer_row_column_head_cond',
                      lambda v: proj.remap_records(v, warn=warn))):
        val = getattr(cMF, name, None)
        if not val:
            continue
        if isinstance(val, dict):
            setattr(m, name, {k: fn(v) for k, v in val.items()})
        else:
            setattr(m, name, [fn(v) for v in val])

    # ---- the MMsoil spatial inputs
    grids_out = {}
    for name, arr in (grids or {}).items():
        a = np.asarray(arr)
        is_zone = name in ZONE_FIELDS
        dt = int if is_zone else float
        fill = (0 if is_zone else hnoflo)
        if a.ndim == 3:
            # per-zone stacks such as gridVEGarea (nveg, nrow, ncol)
            grids_out[name] = np.stack(
                [proj.sample2d(a[k], fill=fill, dtype=dt, how=how,
                               valid=src_active) for k in range(a.shape[0])])
        elif a.ndim == 2:
            grids_out[name] = proj.sample2d(arr, fill=fill, dtype=dt, how=how,
                                            valid=src_active)
        else:
            grids_out[name] = arr

    nactive = int(np.count_nonzero(m.outcropL > 0))
    if nactive == 0:
        raise MeshProjectionError(
            'no active cell survived the projection: the mesh (ncpl=%d) and the '
            'source grid do not overlap. Check the mesh origin and CRS.' % ncpl)

    # Averaging top and botm independently can invert them, because a source
    # cell excluded from one weighted mean may be included in the other (top is
    # valid wherever the soil column outcrops; botm[k] wherever layer k
    # exists). MF6 would reject that later and less clearly.
    #
    # The condition checked is MF6's own, layer by layer and only where the
    # layer is ACTIVE -- an idomain=0 cell's geometry never reaches the solver,
    # so requiring it to be ordered would reject usable meshes.
    top_k = np.asarray(m.top)[:, 0]
    for k in range(nlay):
        bot_k = np.asarray(m.botm)[k][:, 0]
        live = (np.abs(m.ibound[k][:, 0]) != 0)
        bad = live & ~(top_k > bot_k)
        if bad.any():
            i0 = int(np.flatnonzero(bad)[0])
            raise MeshProjectionError(
                '%d active mesh cell(s) in layer %d ended up with the cell top '
                'at or below its bottom after resampling (first: icell2d %d, '
                'top %.3f, botm %.3f). This is a resampling artefact where the '
                'mesh straddles the model edge; try grid.resample = "centre", '
                'or a mesh that does not overhang the active domain.'
                % (int(bad.sum()), k, i0, float(top_k[i0]), float(bot_k[i0])))
        top_k = bot_k
    return m, grids_out, proj
