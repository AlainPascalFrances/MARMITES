# -*- coding: utf-8 -*-
"""WP1d -- wrapping VECTOR layers onto the model grid.

``marmites_mesh.py`` resamples RASTERS onto a mesh. This module is the other
half of the new input paradigm: it takes shapefiles -- points, lines and
polygons, with or without attribute columns -- and produces the per-cell
arrays the model wants, on whichever grid panel 1 defined (structured,
Voronoi or quadtree). Nothing here knows which of the three it is working on:
a grid is a list of convex cell polygons and that is all this module needs.

Why a new module rather than an extension of ``marmites_mesh``
--------------------------------------------------------------
``MeshProjection`` answers "which source RASTER cells does this mesh cell
overlap"; both sides are grids, and the source side is an axis-aligned
lattice, which is what makes its clipping fast. Here the source is an
arbitrary geometry collection with attributes, and the questions are
different in kind -- majority class, area fraction, area-weighted mean,
length of stream inside a cell, which cell holds this observation point. The
two share only the Sutherland-Hodgman idea, and even that runs the other way
round: there the clip window is a raster cell, here it is a MODEL cell.

Conventions
-----------
* Every array returned is shaped ``(nrow, ncol)``. On a mesh that is
  ``(ncpl, 1)`` -- the WP1c convention, so the row index IS the icell2d and
  nothing downstream has to know the difference.
* Sutherland-Hodgman requires a CONVEX clip window. Model cells are convex
  in all three grid kinds (rectangles, quadtree squares, Voronoi cells), so
  the subject geometry may be as irregular as it likes -- which is the point,
  since it comes from a digitised map.
* Shapefiles are read from ``DATA_ROOT/GIS`` and NEVER enter the repository.
  This module is called by the converter (WP1) and by the front-end preview;
  the model itself only ever sees the derived tables and arrays.

Areas are exact, not sampled. A cell that is 37 % alluvium gets 37, not the
class of whatever happened to sit under its centre.
"""

import os

import numpy as np

__all__ = ['VectorError', 'OVERLAY_MODES', 'Layer', 'TargetGrid',
           'overlay_polygons', 'burn_lines', 'locate_points',
           'coverage_report', 'write_geojson']

# 'majority'      the class covering the largest area of the cell
# 'area_fraction' percentage of the cell covered, 0..100  (what VEGarea wants)
# 'area_mean'     area-weighted mean of a numeric column
# 'presence'      1 where any feature touches the cell, else 0
OVERLAY_MODES = ('majority', 'area_fraction', 'area_mean', 'presence')

LINE_MODES = ('longest', 'length', 'presence')


class VectorError(Exception):
    """A vector layer cannot serve the use the configuration asks of it."""


# =====================================================================
#  geometry primitives
# =====================================================================

def _signed_area(pts):
    """Shoelace signed area: positive counter-clockwise, negative clockwise.

    The sign carries meaning here. ESRI writes outer rings clockwise and
    holes counter-clockwise, so summing ``-signed_area`` over the rings of a
    shape gives the net area with holes already subtracted -- no ring
    containment test needed.
    """
    a = np.asarray(pts, dtype=float)
    if a.shape[0] < 3:
        return 0.0
    x, y = a[:, 0], a[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _to_ccw(poly):
    """Return the polygon counter-clockwise, as the clipper below assumes."""
    p = [(float(x), float(y)) for x, y in poly]
    if len(p) > 2 and p[0] == p[-1]:
        p = p[:-1]
    return p if _signed_area(p) >= 0.0 else p[::-1]


def _edges_ccw(clip):
    """Edges of a CCW convex polygon as (ax, ay, bx, by) tuples."""
    return [(clip[k][0], clip[k][1],
             clip[(k + 1) % len(clip)][0], clip[(k + 1) % len(clip)][1])
            for k in range(len(clip))]


def _inside(px, py, e):
    """True when (px, py) is on the inner side of CCW edge ``e``."""
    return (e[2] - e[0]) * (py - e[1]) - (e[3] - e[1]) * (px - e[0]) >= 0.0


def _intersect(p, q, e):
    """Intersection of segment p-q with the infinite line through edge ``e``."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    ex, ey = e[2] - e[0], e[3] - e[1]
    den = ex * dy - ey * dx
    if den == 0.0:                      # parallel: keep the endpoint
        return (q[0], q[1])
    t = (ex * (p[1] - e[1]) - ey * (p[0] - e[0])) / den
    return (p[0] - dx * t, p[1] - dy * t)


def _clip_to_convex(subject, clip):
    """Sutherland-Hodgman clip of ``subject`` against the CONVEX ``clip``.

    ``subject`` may be non-convex and self-touching -- a digitised soil
    polygon usually is. ``clip`` must be convex and is taken CCW. Returns the
    clipped ring (possibly empty) as a list of (x, y).
    """
    out = [(float(x), float(y)) for x, y in subject]
    if len(out) > 2 and out[0] == out[-1]:
        out = out[:-1]
    for e in _edges_ccw(clip):
        if not out:
            return []
        buf, prev = [], out[-1]
        pi = _inside(prev[0], prev[1], e)
        for cur in out:
            ci = _inside(cur[0], cur[1], e)
            if ci:
                if not pi:
                    buf.append(_intersect(prev, cur, e))
                buf.append(cur)
            elif pi:
                buf.append(_intersect(prev, cur, e))
            prev, pi = cur, ci
        out = buf
    return out


def _clip_segment_to_convex(p, q, clip):
    """Cyrus-Beck: the portion of segment p-q inside a CCW convex polygon.

    Returns ``(length, (p_in, q_in))`` or ``(0.0, None)``. Used to measure how
    much stream lies in a cell, which is what decides which cell a reach
    belongs to when a segment crosses several.
    """
    dx, dy = q[0] - p[0], q[1] - p[1]
    if dx == 0.0 and dy == 0.0:
        return 0.0, None
    t0, t1 = 0.0, 1.0
    for e in _edges_ccw(clip):
        nx, ny = -(e[3] - e[1]), (e[2] - e[0])      # inward normal of a CCW edge
        den = nx * dx + ny * dy
        num = nx * (p[0] - e[0]) + ny * (p[1] - e[1])
        if den == 0.0:
            if num < 0.0:
                return 0.0, None                    # parallel and outside
            continue
        t = -num / den
        if den > 0.0:
            t0 = max(t0, t)
        else:
            t1 = min(t1, t)
        if t0 > t1:
            return 0.0, None
    a = (p[0] + dx * t0, p[1] + dy * t0)
    b = (p[0] + dx * t1, p[1] + dy * t1)
    return float(np.hypot(b[0] - a[0], b[1] - a[1])), (a, b)


def _point_in_convex(x, y, clip):
    return all(_inside(x, y, e) for e in _edges_ccw(clip))


# =====================================================================
#  the source: a shapefile
# =====================================================================

_KIND = {1: 'point', 11: 'point', 21: 'point', 8: 'point', 18: 'point',
         28: 'point',
         3: 'line', 13: 'line', 23: 'line',
         5: 'polygon', 15: 'polygon', 25: 'polygon'}


def _geojson_rings(gtype, coords):
    """Flatten a GeoJSON geometry into a list of coordinate rings/parts."""
    if gtype == 'LineString':
        return [coords]
    if gtype in ('MultiLineString', 'Polygon'):
        return list(coords)
    out = []                                   # MultiPolygon
    for poly in coords:
        out.extend(poly)
    return out


def write_geojson(path, kind, geometries, properties, crs='', provenance=None,
                  precision=2):
    """Write a derived layer into the dataset folder (WP1 tier A).

    ``geometries`` is one list of rings (or a single point) per feature, in
    MODEL CRS. Coordinates are rounded to ``precision`` decimals: these are
    metres, so two decimals is a centimetre -- far finer than any of the
    mapping behind them, and it roughly halves the file.

    The ``marmites`` member carries the provenance WP1 established: which
    shapefile this came from, its size and mtime, and the CRS as REPORTED by
    its .prj. A derived file with no provenance is not evidence of anything.
    """
    import json

    def rnd(xy):
        return [round(float(xy[0]), precision), round(float(xy[1]), precision)]

    feats = []
    for geom, props in zip(geometries, properties):
        if kind == 'point':
            g = {'type': 'Point', 'coordinates': rnd(geom)}
        elif kind == 'line':
            parts = [[rnd(p) for p in part] for part in geom]
            g = ({'type': 'LineString', 'coordinates': parts[0]} if len(parts) == 1
                 else {'type': 'MultiLineString', 'coordinates': parts})
        else:
            rings = [[rnd(p) for p in ring] for ring in geom]
            g = {'type': 'Polygon', 'coordinates': rings}
        feats.append({'type': 'Feature', 'geometry': g,
                      'properties': dict(props)})
    doc = {'type': 'FeatureCollection',
           'marmites': dict(provenance or {}, crs=crs),
           'features': feats}
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(doc, fh, separators=(',', ':'))
    return path


class _Shape:
    """The little of pyshp's shape interface this module uses."""

    __slots__ = ('points', 'parts', 'bbox')

    def __init__(self, points, parts, bbox=None):
        self.points = points
        self.parts = parts
        if bbox is None and points:
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            bbox = (min(xs), min(ys), max(xs), max(ys))
        self.bbox = bbox


class Layer:
    """One vector layer: a shapefile, or a derived GeoJSON in the dataset.

    Deliberately thin: no geometry engine, no CRS transformation. The CRS is
    REPORTED from the ``.prj`` rather than assumed, exactly as WP1 already
    does for the converter -- a silent reprojection is the one failure mode
    that produces a plausible-looking wrong answer.

    Two tiers, as WP1 established them. Shapefiles live in ``DATA_ROOT/GIS``
    and are read by the CONVERTER; the converter writes a derived GeoJSON
    into the dataset folder, carrying only the geometry and the columns the
    model was configured to use, and that is what a run reads. Both arrive
    here through the same class, so the wrapping code never knows which tier
    it was handed.
    """

    def __init__(self, path):
        if not os.path.exists(path):
            raise VectorError('vector layer not found: %s' % path)
        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0]
        if path.lower().endswith(('.geojson', '.json')):
            self._read_geojson(path)
        else:
            self._read_shapefile(path)

    # -- readers ----------------------------------------------------------
    def _read_shapefile(self, path):
        import shapefile                       # pyshp

        rdr = self._open(shapefile, path)
        self.shape_type = int(rdr.shapeType)
        self.kind = _KIND.get(self.shape_type)
        if self.kind is None:
            raise VectorError('%s: unsupported shape type %d'
                              % (self.name, self.shape_type))
        self.fields = [f[0] for f in rdr.fields[1:]]
        self.records = [list(r) for r in rdr.records()]
        self.shapes = list(rdr.shapes())
        self.bbox = tuple(float(v) for v in rdr.bbox)
        prj = os.path.splitext(path)[0] + '.prj'
        self.crs_wkt = (open(prj, encoding='utf-8', errors='replace').read()
                        if os.path.exists(prj) else '')
        rdr.close()

    def _read_geojson(self, path):
        import json

        with open(path, encoding='utf-8') as fh:
            doc = json.load(fh)
        feats = doc.get('features') or []
        self.provenance = doc.get('marmites') or {}
        self.crs_wkt = str(self.provenance.get('crs', ''))
        kinds, names = set(), []
        self.shapes, self.records = [], []
        for ft in feats:
            geom = ft.get('geometry') or {}
            gtype = geom.get('type', '')
            coords = geom.get('coordinates')
            if gtype == 'Point':
                kinds.add('point')
                pts, parts = [(float(coords[0]), float(coords[1]))], [0]
            elif gtype in ('LineString', 'MultiLineString', 'Polygon',
                           'MultiPolygon'):
                kinds.add('line' if 'Line' in gtype else 'polygon')
                rings = _geojson_rings(gtype, coords)
                pts, parts = [], []
                for r in rings:
                    parts.append(len(pts))
                    pts.extend((float(x), float(y)) for x, y in r)
            else:
                raise VectorError('%s: unsupported geometry %r'
                                  % (self.name, gtype))
            props = ft.get('properties') or {}
            for k in props:
                if k not in names:
                    names.append(k)
            self.shapes.append(_Shape(pts, parts))
            self.records.append(props)
        if len(kinds) > 1:
            raise VectorError('%s: mixed geometry types %s'
                              % (self.name, sorted(kinds)))
        self.kind = kinds.pop() if kinds else 'polygon'
        self.shape_type = {'point': 1, 'line': 3, 'polygon': 5}[self.kind]
        self.fields = names
        self.records = [[r.get(n) for n in names] for r in self.records]
        xs = [p[0] for s in self.shapes for p in s.points]
        ys = [p[1] for s in self.shapes for p in s.points]
        self.bbox = ((min(xs), min(ys), max(xs), max(ys))
                     if xs else (0.0, 0.0, 0.0, 0.0))

    @staticmethod
    def _open(shapefile, path):
        """Open the layer, surviving the encoding names GIS software writes.

        ArcGIS drops a ``.cpg`` saying ``ansi 1252``, which pyshp hands
        straight to ``bytes.decode`` and Python does not know -- so a file
        re-saved from ArcGIS stops opening, with a LookupError that says
        nothing about shapefiles. cp1252 is what that name means; latin-1 is
        the last resort because it cannot fail.
        """
        try:
            return shapefile.Reader(path)
        except (LookupError, UnicodeDecodeError):
            pass
        for enc in ('cp1252', 'latin-1'):
            try:
                return shapefile.Reader(path, encoding=enc,
                                        encodingErrors='replace')
            except (LookupError, UnicodeDecodeError):
                continue
        raise VectorError('%s: cannot decode the attribute table; the .cpg '
                          'file names an encoding Python does not know'
                          % os.path.basename(path))

    def __len__(self):
        return len(self.shapes)

    def require(self, *names):
        """Fail with a message naming what the layer DOES have."""
        miss = [n for n in names if n not in self.fields]
        if miss:
            raise VectorError(
                '%s: column(s) %s not in the layer, which has: %s'
                % (self.name, ', '.join(miss), ', '.join(self.fields) or '(none)'))

    def column(self, name, dtype=None):
        """One attribute column as a list, or a constant if ``name`` is None."""
        if name is None:
            return [1] * len(self.records)
        self.require(name)
        k = self.fields.index(name)
        vals = [r[k] for r in self.records]
        if dtype is not None:
            out = []
            for v in vals:
                try:
                    out.append(dtype(v))
                except (TypeError, ValueError):
                    out.append(None)
            return out
        return vals

    def rings(self, i):
        """Rings of feature ``i`` as lists of (x, y), holes included."""
        sh = self.shapes[i]
        pts = sh.points
        parts = list(sh.parts) + [len(pts)]
        return [[(float(x), float(y)) for x, y in pts[parts[k]:parts[k + 1]]]
                for k in range(len(parts) - 1) if parts[k + 1] - parts[k] >= 2]

    def feature_bbox(self, i):
        sh = self.shapes[i]
        bb = getattr(sh, 'bbox', None)
        if bb is not None and len(bb) == 4:
            return tuple(float(v) for v in bb)
        p = np.asarray(sh.points, dtype=float)
        return (p[:, 0].min(), p[:, 1].min(), p[:, 0].max(), p[:, 1].max())


# =====================================================================
#  the target: the model grid, whatever kind it is
# =====================================================================

class TargetGrid:
    """The model grid as a list of convex cell polygons, plus a bucket index.

    Built either from structured spacings or from DISV gridprops. The
    ``shape`` it reports is what every array returned by this module is
    reshaped to, so a caller never has to branch on the grid kind.
    """

    def __init__(self, polygons, shape):
        self.polygons = [_to_ccw(p) for p in polygons]
        self.ncell = len(self.polygons)
        self.shape = tuple(shape)
        if int(np.prod(self.shape)) != self.ncell:
            raise VectorError('shape %s does not hold %d cells'
                              % (self.shape, self.ncell))
        self._bbox = np.array(
            [(min(x for x, _ in p), min(y for _, y in p),
              max(x for x, _ in p), max(y for _, y in p))
             for p in self.polygons], dtype=float)
        self.area = np.array([abs(_signed_area(p)) for p in self.polygons],
                             dtype=float)
        self.extent = (float(self._bbox[:, 0].min()), float(self._bbox[:, 1].min()),
                       float(self._bbox[:, 2].max()), float(self._bbox[:, 3].max()))
        self._build_index()

    # -- constructors -----------------------------------------------------
    @classmethod
    def structured(cls, delr, delc, xllcorner, yllcorner):
        delr = np.asarray(delr, dtype=float)
        delc = np.asarray(delc, dtype=float)
        x = float(xllcorner) + np.concatenate(([0.0], np.cumsum(delr)))
        ytop = float(yllcorner) + float(delc.sum())
        y = ytop - np.concatenate(([0.0], np.cumsum(delc)))
        polys = []
        for i in range(len(delc)):                  # row 0 is the NORTHERNMOST
            for j in range(len(delr)):
                polys.append([(x[j], y[i + 1]), (x[j + 1], y[i + 1]),
                              (x[j + 1], y[i]), (x[j], y[i])])
        return cls(polys, (len(delc), len(delr)))

    @classmethod
    def from_gridprops(cls, gridprops):
        vxy = {int(v[0]): (float(v[1]), float(v[2]))
               for v in gridprops['vertices']}
        polys = []
        for rec in gridprops['cell2d']:
            ids = [int(v) for v in rec[4:]]
            if len(ids) > 2 and ids[0] == ids[-1]:
                ids = ids[:-1]
            polys.append([vxy[k] for k in ids])
        ncpl = int(gridprops.get('ncpl', len(polys)))
        return cls(polys, (ncpl, 1))                # the WP1c (ncpl, 1) convention

    @classmethod
    def from_cMF(cls, cMF):
        """Whatever grid the run is on -- DISV if it has gridprops, else DIS."""
        gp = getattr(cMF, 'gridprops', None)
        if gp:
            return cls.from_gridprops(gp)
        return cls.structured(cMF.delr, cMF.delc,
                              float(cMF.xllcorner), float(cMF.yllcorner))

    # -- candidate lookup -------------------------------------------------
    def _build_index(self):
        """Uniform bucket index over the cell bounding boxes.

        A quadtree would be tidier; a bucket grid is ~20 lines and fast enough
        for the 15 586 crown polygons of La Mata against a 1 000-cell mesh.
        """
        x0, y0, x1, y1 = self.extent
        n = max(1, int(np.sqrt(max(self.ncell, 1))))
        self._nb = n
        self._bx = (x1 - x0) / n if x1 > x0 else 1.0
        self._by = (y1 - y0) / n if y1 > y0 else 1.0
        self._x0, self._y0 = x0, y0
        buckets = {}
        for ic in range(self.ncell):
            b = self._bbox[ic]
            for bi in range(self._bi(b[0]), self._bi(b[2]) + 1):
                for bj in range(self._bj(b[1]), self._bj(b[3]) + 1):
                    buckets.setdefault((bi, bj), []).append(ic)
        self._buckets = buckets

    def _bi(self, x):
        return int(min(self._nb - 1, max(0, (x - self._x0) // self._bx)))

    def _bj(self, y):
        return int(min(self._nb - 1, max(0, (y - self._y0) // self._by)))

    def candidates(self, bbox):
        """Cells whose bounding box may intersect ``bbox`` (xmin,ymin,xmax,ymax)."""
        out = set()
        for bi in range(self._bi(bbox[0]), self._bi(bbox[2]) + 1):
            for bj in range(self._bj(bbox[1]), self._bj(bbox[3]) + 1):
                out.update(self._buckets.get((bi, bj), ()))
        bb = self._bbox
        return [ic for ic in out
                if not (bb[ic, 2] < bbox[0] or bb[ic, 0] > bbox[2]
                        or bb[ic, 3] < bbox[1] or bb[ic, 1] > bbox[3])]

    def cell_containing(self, x, y):
        """icell of the cell holding (x, y), or -1. Ties go to the lowest index."""
        for ic in sorted(self.candidates((x, y, x, y))):
            if _point_in_convex(x, y, self.polygons[ic]):
                return ic
        return -1

    def unravel(self, ic):
        """(row, col) of a flat cell index, in this grid's own shape."""
        return divmod(int(ic), self.shape[1])


# =====================================================================
#  the operations
# =====================================================================

def _feature_rings(layer, i):
    """Rings of feature ``i`` plus the sign that makes its OUTER ring positive.

    The ESRI spec writes outer rings clockwise and holes counter-clockwise,
    but plenty of real shapefiles do not, and trusting the spec turns a
    reversed file into zero coverage everywhere -- silently. So the
    orientation is taken from the feature's own largest ring, which is its
    outer one, and every other ring is read relative to that.

    Computed ONCE per feature: a big background polygon can be a candidate
    for every cell in the grid, and re-reading its vertices inside that loop
    is what makes a naive overlay unusable on 15 000 features.
    """
    rings = layer.rings(i)
    if not rings:
        return [], 1.0
    sgn = [_signed_area(r) for r in rings]
    k = int(np.argmax(np.abs(sgn)))
    return rings, (-1.0 if sgn[k] < 0.0 else 1.0)


def _area_in_cell(rings, flip, clip, bboxes=None):
    """Net area of one feature inside the convex cell ``clip``, holes out."""
    a = 0.0
    for n, ring in enumerate(rings):
        if bboxes is not None:
            b, c = bboxes[n], clip[1]
            if b[2] < c[0] or b[0] > c[2] or b[3] < c[1] or b[1] > c[3]:
                continue
        piece = _clip_to_convex(ring, clip[0])
        if len(piece) >= 3:
            a += flip * _signed_area(piece)
    return max(a, 0.0)


def overlay_polygons(layer, grid, field=None, how='majority', fill=0,
                     dtype=float, classes=None, select=None):
    """Wrap a polygon layer onto ``grid``.

    Parameters
    ----------
    layer : Layer
    grid : TargetGrid
    field : str or None
        Attribute column carrying the value. ``None`` means "the geometry
        itself is the value" -- valid for 'area_fraction' and 'presence'.
    how : {'majority', 'area_fraction', 'area_mean', 'presence'}
        'majority'      the class covering the most area of the cell
        'area_fraction' percent of the cell covered, 0..100
        'area_mean'     area-weighted mean of a numeric column
        'presence'      1 where any feature overlaps at all
    fill : scalar
        Value for cells no feature reaches.
    classes : sequence or None
        For 'area_fraction', restrict to features whose ``field`` value is in
        ``classes``; for 'majority', the allowed codes (others are ignored).
    select : callable or None
        ``select(record_value) -> bool``, applied to ``field``. Takes
        precedence over ``classes``; this is how "blank Species means grass"
        is expressed without teaching this module about vegetation.

    Returns
    -------
    (array, report) : ndarray shaped ``grid.shape``, and a dict of diagnostics.
    """
    if how not in OVERLAY_MODES:
        raise VectorError('unknown overlay mode %r (known: %s)'
                          % (how, ', '.join(OVERLAY_MODES)))
    if layer.kind != 'polygon':
        raise VectorError('%s is a %s layer; overlay_polygons needs polygons'
                          % (layer.name, layer.kind))
    if how in ('majority', 'area_mean') and field is None:
        raise VectorError('how=%r needs a field' % how)

    vals = layer.column(field)
    if how == 'area_mean':
        vals = layer.column(field, dtype=float)

    keep = np.ones(len(layer), dtype=bool)
    if select is not None:
        keep = np.array([bool(select(v)) for v in vals], dtype=bool)
    elif classes is not None:
        cset = set(classes)
        keep = np.array([v in cset for v in vals], dtype=bool)

    # area per (cell, value) for majority; running sums for the rest
    acc = [dict() for _ in range(grid.ncell)] if how == 'majority' else None
    tot = np.zeros(grid.ncell, dtype=float)
    wsum = np.zeros(grid.ncell, dtype=float)
    touched = np.zeros(grid.ncell, dtype=bool)
    n_used = 0

    for i in range(len(layer)):
        if not keep[i]:
            continue
        bb = layer.feature_bbox(i)
        cand = grid.candidates(bb)
        if not cand:
            continue
        n_used += 1
        v = vals[i]
        rings, flip = _feature_rings(layer, i)
        # Per-ring bounding boxes: a multipart feature -- the background
        # polygons of a crown map are exactly that -- otherwise clips every
        # one of its parts against every candidate cell.
        rbb = [(min(x for x, _ in r), min(y for _, y in r),
                max(x for x, _ in r), max(y for _, y in r)) for r in rings]
        for ic in cand:
            a = _area_in_cell(rings, flip,
                              (grid.polygons[ic], grid._bbox[ic]), rbb)
            if a <= 0.0:
                continue
            touched[ic] = True
            if how == 'majority':
                acc[ic][v] = acc[ic].get(v, 0.0) + a
            elif how == 'area_fraction':
                tot[ic] += a
            elif how == 'area_mean':
                if v is not None:
                    tot[ic] += a * v
                    wsum[ic] += a
            # 'presence' needs only `touched`

    out = np.full(grid.ncell, fill, dtype=dtype)
    if how == 'majority':
        allowed = None if classes is None else set(classes)
        for ic in range(grid.ncell):
            d = acc[ic]
            if allowed is not None:
                d = {k: a for k, a in d.items() if k in allowed}
            if d:
                out[ic] = max(d.items(), key=lambda kv: kv[1])[0]
    elif how == 'area_fraction':
        with np.errstate(invalid='ignore', divide='ignore'):
            frac = 100.0 * tot / np.where(grid.area > 0, grid.area, np.nan)
        out = np.where(np.isfinite(frac), frac, float(fill)).astype(dtype)
    elif how == 'area_mean':
        ok = wsum > 0.0
        out = out.astype(float)
        out[ok] = tot[ok] / wsum[ok]
        out = out.astype(dtype)
    else:                                   # presence
        out[touched] = 1

    report = {
        'layer': layer.name, 'how': how, 'field': field,
        'features': len(layer), 'features_used': n_used,
        'cells_touched': int(touched.sum()), 'cells': grid.ncell,
        'coverage_pct': 100.0 * float(touched.sum()) / max(grid.ncell, 1),
    }
    if how == 'area_fraction':
        report['mean_fraction_pct'] = float(np.nanmean(out))
        report['max_fraction_pct'] = float(np.nanmax(out))
        report['over_100_cells'] = int(np.sum(np.asarray(out, float) > 100.0 + 1e-6))
    return out.reshape(grid.shape), report


def burn_lines(layer, grid, field=None, how='longest', fill=0, dtype=float):
    """Wrap a line layer onto ``grid``.

    'longest'  the value of the feature contributing the most length to the
               cell -- the right rule for a per-segment attribute, because a
               cell crossed by two reaches takes the one that dominates it
    'length'   total length of the layer inside each cell [m]
    'presence' 1 where any line passes through

    Returns ``(array, report)``; the report carries ``cell_length`` (the per
    cell length in metres) whatever the mode, because the SFR builder wants it.
    """
    if how not in LINE_MODES:
        raise VectorError('unknown line mode %r (known: %s)'
                          % (how, ', '.join(LINE_MODES)))
    if layer.kind != 'line':
        raise VectorError('%s is a %s layer; burn_lines needs lines'
                          % (layer.name, layer.kind))
    vals = layer.column(field)
    best = np.zeros(grid.ncell, dtype=float)
    length = np.zeros(grid.ncell, dtype=float)
    pick = np.full(grid.ncell, -1, dtype=int)

    for i in range(len(layer)):
        for part in layer.rings(i):                # 'rings' = parts, for a line
            for k in range(len(part) - 1):
                p, q = part[k], part[k + 1]
                seg_bb = (min(p[0], q[0]), min(p[1], q[1]),
                          max(p[0], q[0]), max(p[1], q[1]))
                for ic in grid.candidates(seg_bb):
                    ln, _ = _clip_segment_to_convex(p, q, grid.polygons[ic])
                    if ln <= 0.0:
                        continue
                    length[ic] += ln
                    if ln > best[ic]:
                        best[ic], pick[ic] = ln, i

    out = np.full(grid.ncell, fill, dtype=dtype)
    if how == 'longest':
        for ic in np.nonzero(pick >= 0)[0]:
            out[ic] = vals[pick[ic]]
    elif how == 'length':
        out = length.astype(dtype)
    else:
        out[length > 0.0] = 1

    report = {
        'layer': layer.name, 'how': how, 'field': field,
        'features': len(layer),
        'cells_touched': int(np.sum(length > 0.0)), 'cells': grid.ncell,
        'total_length_m': float(length.sum()),
        'cell_length': length.reshape(grid.shape),
    }
    return out.reshape(grid.shape), report


def locate_points(layer, grid, name_field=None):
    """Cell of every point in ``layer``.

    Returns a list of dicts with ``name``, ``x``, ``y``, ``icell``, ``row``,
    ``col`` and ``inside``. Points that fall outside are REPORTED, not
    silently clipped -- on a Voronoi mesh clipped to the catchment, a
    piezometer just outside the boundary is a real thing to be told about.
    """
    if layer.kind != 'point':
        raise VectorError('%s is a %s layer; locate_points needs points'
                          % (layer.name, layer.kind))
    names = layer.column(name_field) if name_field else \
        ['pt%d' % (i + 1) for i in range(len(layer))]
    out = []
    for i in range(len(layer)):
        pt = layer.shapes[i].points[0]
        ic = grid.cell_containing(float(pt[0]), float(pt[1]))
        row, col = grid.unravel(ic) if ic >= 0 else (-1, -1)
        out.append({'name': str(names[i]), 'x': float(pt[0]), 'y': float(pt[1]),
                    'icell': int(ic), 'row': int(row), 'col': int(col),
                    'inside': ic >= 0})
    return out


def coverage_report(reports, indent='  '):
    """One line per overlay, for the run log and the front-end preview."""
    lines = []
    for r in reports:
        if 'total_length_m' in r:
            lines.append('%s%-22s %-9s %5d feature(s) -> %4d/%d cell(s), '
                         '%.0f m' % (indent, r['layer'], r['how'],
                                     r['features'], r['cells_touched'],
                                     r['cells'], r['total_length_m']))
        else:
            extra = ''
            if 'mean_fraction_pct' in r:
                extra = ', mean %.2f %% max %.2f %%' % (r['mean_fraction_pct'],
                                                        r['max_fraction_pct'])
                if r.get('over_100_cells'):
                    extra += ', %d cell(s) OVER 100 %%' % r['over_100_cells']
            lines.append('%s%-22s %-9s %5d feature(s) -> %4d/%d cell(s)%s'
                         % (indent, r['layer'], r['how'], r['features'],
                            r['cells_touched'], r['cells'], extra))
    return '\n'.join(lines)
