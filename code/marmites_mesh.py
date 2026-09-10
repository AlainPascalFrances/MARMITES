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

Sampling
--------
WP1c.1 samples at the mesh cell CENTRE, which is exact when the mesh is a
refinement of the source raster and lossy when a mesh cell is larger than a
raster cell. WP1c.3 replaces this with area-weighted means for continuous
fields and majority vote for zone rasters; ``how=`` is already threaded
through so that change lands in one place.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import copy

import numpy as np

__all__ = ['MeshProjection', 'project_model', 'MeshProjectionError',
           'ZONE_FIELDS']

# Field-by-field sampling rules. 'zone' means an integer class code, which may
# never be interpolated or averaged -- a soil zone halfway between 1 and 3 is
# not 2. Continuous fields may be averaged (WP1c.3).
ZONE_FIELDS = ('gridSOIL', 'gridMETEO', 'gridIRR', 'ibound', 'iuzfbnd',
               'outcropL')


class MeshProjectionError(Exception):
    """Raised when a model cannot be projected onto a mesh."""


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
        inside = ((col >= 0) & (col < self.ncol_src)
                  & (row >= 0) & (row < self.nrow_src))
        return (np.clip(row, 0, self.nrow_src - 1),
                np.clip(col, 0, self.ncol_src - 1), inside)

    def locate(self, x, y):
        """Public point lookup: returns (row, col, inside) for scalars/arrays."""
        xy = np.column_stack([np.atleast_1d(x), np.atleast_1d(y)])
        return self._locate(xy)

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

    # -- sampling ---------------------------------------------------------
    def sample2d(self, arr, fill=np.nan, dtype=float, how='centre'):
        """Sample a source ``(nrow, ncol)`` array onto ``(ncpl, 1)``.

        ``fill`` is used where the mesh centre falls outside the source grid.

        A masked source gives a masked result, and the DATA UNDER THE MASK is
        carried through rather than replaced by ``fill``. That is not
        cosmetic: the driver hands MF6 ``np.asarray(cMF.top)``, which drops the
        mask, so overwriting the masked data would quietly change the values
        MODFLOW sees at inactive cells.
        """
        if how != 'centre':
            raise MeshProjectionError(
                "sampling mode %r arrives with WP1c.3; only 'centre' is "
                "implemented" % (how,))
        a = np.ma.getdata(arr) if np.ma.isMaskedArray(arr) else np.asarray(arr)
        if a.shape != (self.nrow_src, self.ncol_src):
            raise MeshProjectionError(
                'expected a (%d, %d) source array, got %s'
                % (self.nrow_src, self.ncol_src, (a.shape,)))
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

    def sample3d(self, arr, fill=np.nan, dtype=float, how='centre'):
        """Sample a source ``(nlay, nrow, ncol)`` array onto ``(nlay, ncpl, 1)``."""
        a = np.asarray(arr)
        if a.ndim != 3:
            raise MeshProjectionError('expected a 3-D array, got %s' % (a.shape,))
        masked = np.ma.isMaskedArray(arr)
        mask = np.ma.getmaskarray(arr) if masked else None
        layers = []
        for k in range(a.shape[0]):
            src = np.ma.array(a[k], mask=mask[k]) if masked else a[k]
            layers.append(self.sample2d(src, fill=fill, dtype=dtype, how=how))
        return np.stack(layers)

    def sample_layer_property(self, val, fill=np.nan, how='centre'):
        """Project a cMF layer property, which may be a list of scalars OR of
        arrays (``hk_actual`` is arrays, ``vka_actual`` is floats).

        Scalars are left alone: they are already grid-independent, and
        expanding them here would only make the MF6 build carry a full array
        for a constant.
        """
        if not isinstance(val, (list, tuple)):
            return val
        out = []
        for v in val:
            a = np.asarray(v)
            if a.ndim == 2:
                out.append(self.sample2d(v, fill=fill, how=how))
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


def project_model(cMF, gridprops, grids, warn=None):
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
    m.ibound = proj.sample3d(cMF.ibound, fill=0, dtype=int)
    m.iuzfbnd = proj.sample2d(cMF.iuzfbnd, fill=0, dtype=int)

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
    m.elev = _mask_sentinel(proj.sample2d(cMF.elev, fill=hnoflo), hnoflo)
    m.top = _mask_sentinel(proj.sample2d(cMF.top, fill=hnoflo), hnoflo)
    m.botm = proj.sample3d(np.asarray(cMF.botm), fill=hnoflo)
    if getattr(cMF, 'strt', None) is not None:
        m.strt = proj.sample3d(np.asarray(cMF.strt), fill=hnoflo)

    # ---- aquifer properties
    for name in ('hk_actual', 'vka_actual', 'ss_actual', 'sy_actual',
                 'vks_actual', 'thick'):
        val = getattr(cMF, name, None)
        if val is None:
            continue
        if isinstance(val, list):
            setattr(m, name, proj.sample_layer_property(val, fill=hnoflo))
            continue
        a = np.asarray(val)
        if a.ndim == 3:
            setattr(m, name, proj.sample3d(val, fill=hnoflo))
        elif a.ndim == 2:
            setattr(m, name, proj.sample2d(val, fill=hnoflo))

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
            grids_out[name] = np.stack([proj.sample2d(a[k], fill=fill, dtype=dt)
                                        for k in range(a.shape[0])])
        elif a.ndim == 2:
            grids_out[name] = proj.sample2d(arr, fill=fill, dtype=dt)
        else:
            grids_out[name] = arr

    nactive = int(np.count_nonzero(m.outcropL > 0))
    if nactive == 0:
        raise MeshProjectionError(
            'no active cell survived the projection: the mesh (ncpl=%d) and the '
            'source grid do not overlap. Check the mesh origin and CRS.' % ncpl)
    return m, grids_out, proj
