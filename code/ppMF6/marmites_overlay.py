# -*- coding: utf-8 -*-
"""Vector layers onto the model grid: majority, area mean, area share.

The Soil panel asks for its inputs as POLYGON LAYERS -- the soil zones from
a shapefile's code column, the vegetation cover from its class column --
and until this module nothing could turn a polygon into a grid value: the
run read legacy rasters from fixed filenames instead, so the panel's
answers changed nothing. This is the missing step, three operations:

    majority(column)     the value covering the largest share of the cell
    area_mean(column)    the area-weighted mean of the value over the cell
    class_percent(...)   the share of the cell, in %, covered by each class

All three are EXACT AREA OVERLAYS -- a cell 37 % covered gets 37, not
whatever sat under its centre -- done with shapely 2's vectorised spatial
index, so a cell's intersections with the source polygons are found in
bulk rather than pair by pair.

GRID-AGNOSTIC: a "cell" is any polygon. ``structured_cells`` gives the
rectangles of a DIS grid; a mesh's polygons would do equally well.

The polygons come from the DATASET (the converter's GeoJSON, in model CRS),
never from the shapefiles directly: the run must not depend on the GIS
folder, which lives outside the repository.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import json
import os

import numpy as np

__all__ = ['OverlayError', 'Polygons', 'structured_cells', 'overlay',
           'majority', 'area_mean', 'class_percent']


class OverlayError(Exception):
    """A layer cannot be put onto the grid as asked."""


class Polygons(object):
    """The polygons of one layer, with the attribute columns asked for."""

    def __init__(self, geoms, attrs, source=''):
        self.geoms = geoms              # numpy array of shapely geometries
        self.attrs = attrs              # {column: list of values}
        self.source = source

    @classmethod
    def from_geojson(cls, path, columns):
        import shapely
        if not os.path.exists(path):
            raise OverlayError('%s is not there -- run the converter '
                               '(Soil panel, Cartography tab: Update '
                               'dataset), which writes it from the '
                               'cartography' % path)
        with open(path, encoding='utf-8') as fh:
            g = json.load(fh)
        feats = [f for f in g.get('features', []) if f.get('geometry')]
        if not feats:
            raise OverlayError('%s holds no polygons' % path)
        have = set(feats[0].get('properties') or {})
        missing = [c for c in columns if c not in have]
        if missing:
            raise OverlayError(
                '%s has no column %s -- it has %s. The converter exports '
                'the columns the configuration names; run it again on the '
                'Grid panel if the name changed.'
                % (os.path.basename(path), ', '.join(missing),
                   ', '.join(sorted(have))))
        geoms = shapely.from_geojson([json.dumps(f['geometry'])
                                      for f in feats])
        geoms = shapely.make_valid(geoms)
        attrs = {c: [f['properties'].get(c) for f in feats] for c in columns}
        return cls(np.asarray(geoms), attrs, source=path)


def structured_cells(xll, yll, delr, delc):
    """The cell rectangles of a DIS grid, row-major, row 0 at the NORTH.

    Returns a numpy array of shapely boxes of length nrow*ncol, so that
    ``result.reshape(nrow, ncol)`` puts cell (i, j) at [i, j].
    """
    import shapely
    delr = np.asarray(delr, dtype=float).ravel()
    delc = np.asarray(delc, dtype=float).ravel()
    xe = float(xll) + np.concatenate([[0.0], np.cumsum(delr)])
    ytop = float(yll) + float(delc.sum())
    ye = ytop - np.concatenate([[0.0], np.cumsum(delc)])      # decreasing
    x0 = np.tile(xe[:-1], delc.size)
    x1 = np.tile(xe[1:], delc.size)
    y1 = np.repeat(ye[:-1], delr.size)
    y0 = np.repeat(ye[1:], delr.size)
    return np.asarray(shapely.box(x0, y0, x1, y1))


# A polygon with more vertices than this is cut into TILE x TILE m tiles
# before the overlay. La Mata's vegetation layer holds one grass matrix of
# 284,640 vertices (the mean is 38): every mesh cell's 'intersects' test and
# intersection ran against all of it -- 655 s for the query and ~320 s for
# the intersections on the 15,915-cell Voronoi mesh. Intersection does not
# care how a polygon is split, so tiling is exact; each tile keeps the index
# of the polygon it came from.
TILE_VERTICES = 1000
TILE = 100.0


def _tiled(geoms):
    """``(pieces, source_index)``: small polygons as they are, big ones cut
    into tiles with GEOS's rectangle clip (one call per tile: shapely's
    clip_by_rect takes scalar bounds only)."""
    import shapely
    n = shapely.get_num_coordinates(geoms)
    big = np.nonzero(n > TILE_VERTICES)[0]
    if not big.size:
        return geoms, np.arange(len(geoms))
    small = np.nonzero(n <= TILE_VERTICES)[0]
    pieces, src = [geoms[small]], [small]
    for i in big:
        x0, y0, x1, y1 = geoms[i].bounds
        X0, Y0 = np.meshgrid(np.arange(x0, x1, TILE), np.arange(y0, y1, TILE))
        cut = np.array([shapely.clip_by_rect(geoms[i], x, y, x + TILE,
                                             y + TILE)
                        for x, y in zip(X0.ravel(), Y0.ravel())],
                       dtype=object)
        cut = cut[~shapely.is_empty(cut) & (shapely.area(cut) > 0.0)]
        # clip_by_rect may return an invalid ring on a degenerate edge;
        # an invalid operand would make the intersection below raise
        bad = ~shapely.is_valid(cut)
        if bad.any():
            cut[bad] = shapely.make_valid(cut[bad])
        pieces.append(cut)
        src.append(np.full(cut.size, i))
    return (np.concatenate([np.asarray(p, dtype=object) for p in pieces]),
            np.concatenate(src))


def overlay(cells, polys):
    """``(cell_index, polygon_index, area)`` of every non-empty intersection.

    One spatial-index query for the whole grid, then the intersections
    computed vectorised. Only pairs that actually overlap are returned; a
    polygon cut into tiles (TILE_VERTICES) may give several pairs for one
    cell, which every consumer sums.
    """
    import shapely
    geoms, src = _tiled(np.asarray(polys.geoms, dtype=object))
    tree = shapely.STRtree(geoms)
    ci, ti = tree.query(cells, predicate='intersects')
    if ci.size == 0:
        return ci, ti, np.zeros(0)
    area = shapely.area(shapely.intersection(cells[ci], geoms[ti]))
    keep = area > 0.0
    return ci[keep], src[ti][keep], area[keep]


def majority(cells, polys, column, fill, cast=int):
    """Per cell, the value of ``column`` covering the largest area.

    A cell no polygon touches gets ``fill``. Ties go to the value met first,
    which is deterministic for a given layer; they happen only where two
    polygons share a cell exactly half and half.
    """
    ci, pi, area = overlay(cells, polys)
    vals = polys.attrs[column]
    best = {}
    acc = {}
    for c, p, a in zip(ci.tolist(), pi.tolist(), area.tolist()):
        v = vals[p]
        if v is None:
            continue
        acc[(c, v)] = acc.get((c, v), 0.0) + a
    for (c, v), a in acc.items():
        if c not in best or a > best[c][1]:
            best[c] = (v, a)
    out = np.full(len(cells), fill, dtype=float)
    for c, (v, _a) in best.items():
        out[c] = float(cast(v))
    return out


def area_mean(cells, polys, column, fill):
    """Per cell, the area-weighted mean of ``column`` over the covered part.

    Weighted by the COVERED area, not the cell area: a cell half inside the
    layer gets the mean of the half that has a value, and the other half is
    not counted as a zero.
    """
    ci, pi, area = overlay(cells, polys)
    vals = np.array([np.nan if v is None else float(v)
                     for v in polys.attrs[column]])
    v = vals[pi]
    ok = np.isfinite(v)
    num = np.bincount(ci[ok], weights=area[ok] * v[ok], minlength=len(cells))
    den = np.bincount(ci[ok], weights=area[ok], minlength=len(cells))
    out = np.full(len(cells), fill, dtype=float)
    has = den > 0.0
    out[has] = num[has] / den[has]
    return out


def class_percent(cells, polys, column, classes, n):
    """Per class and cell, the percentage of the CELL covered by that class.

    ``classes`` maps a value of ``column`` to a class index, 1-based; values
    it does not list are ignored. Returns an array ``(n, len(cells))`` in %,
    against the cell's full area -- so a cell 37 % covered gets 37, and the
    uncovered rest is simply not vegetated.
    """
    import shapely
    ci, pi, area = overlay(cells, polys)
    cell_area = shapely.area(cells)
    out = np.zeros((n, len(cells)), dtype=float)
    vals = polys.attrs[column]
    idx = np.array([classes.get(str(vals[p]), 0) for p in pi.tolist()],
                   dtype=int) if pi.size else np.zeros(0, dtype=int)
    for k in range(1, n + 1):
        sel = idx == k
        if sel.any():
            out[k - 1] = np.bincount(ci[sel], weights=area[sel],
                                     minlength=len(cells))
    with np.errstate(invalid='ignore', divide='ignore'):
        out = np.where(cell_area > 0, 100.0 * out / cell_area, 0.0)
    # Overlapping source polygons can push a cell past 100 % -- the layer's
    # fault, not the overlay's, but the soil model refuses it, so say so.
    total = out.sum(axis=0)
    if np.any(total > 100.0 + 1e-6):
        worst = float(total.max())
        raise OverlayError(
            '%s: the classes cover up to %.1f %% of a cell -- polygons of the '
            'layer overlap. Repair the layer; the soil model refuses more '
            'than 100 %%.' % (os.path.basename(polys.source) or 'layer', worst))
    return out
