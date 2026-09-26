# -*- coding: utf-8 -*-
"""Build a MODFLOW 6 LAK package for the La Mata ponds (charcas).

The ponds are SUB-GRID
----------------------
``lm_ponds.shp`` holds 12 polygons of 341-2036 m2 against a 2500 m2 cell: every
pond is smaller than one 50 m cell, and two of them (ids 9 and 13) do not even
contain a cell centre. A lake cannot therefore be represented by excavating
cells -- there are no cells to excavate.

MODFLOW 6 has a connection type for exactly this case: ``EMBEDDEDV``. The lake
lives *inside* one host cell, has exactly one connection, and its real geometry
is supplied through a stage-volume-area table rather than being inferred from
the cell. This is also the design the CdL reference model converged on.

On a mesh the pond is RESOLVED
------------------------------
La Mata's Voronoi mesh is refined around the ponds: a pond covers 15-83
cells of ~25 m2. It stays one EMBEDDEDV lake -- connected through the cell
holding its centroid, the seeded cell -- but its FOOTPRINT (the cells whose
centre lies inside it) is what the stream is cut out of, so the stream runs
through the lake rather than past it (``pond_footprints``;
``marmites_mf6._build_ponds``). On the 50 m grid the footprint is the host.

Table shape
-----------
An excavated pond has a flat bottom, so the wetted area would jump from 0 to
the full footprint as a step and a bone-dry lake degenerates numerically. The
table therefore assumes a wedge bathymetry: the surface area grows linearly
from ~0 at the deepest point to the full polygon area at the rim, which lets
storage and fluxes vanish smoothly as the pond dries. ``barea`` (the bed
exchange area) is set equal to ``sarea``: the wetted bed IS the exchange area.

Evaporation
-----------
Written by the COUPLER each stress period, from the Eo forcing (WP1d). It
used to be left to MARMITES, and the package was given none of its own so the
pond would not evaporate twice; MARMITES no longer has a surface store to
evaporate from, so the evaporation followed the water here. The build still
specifies nothing -- the rate arrives through the API.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import numpy as np

__all__ = ['read_pond_polygons', 'pond_footprints', 'assign_pond_cells',
           'lake_table', 'PondLake']

LAK_BEDLEAK = 1e-3     # 1/d, lakebed leakance (clay-lined charca; El-Zehairy 2018)
LAK_SURFDEP = 0.05     # m, smooths the connection wetted area near the bed
POND_DEPTH = 1.5       # m, fallback pond depth when no depth map is supplied


class PondLake(object):
    """One pond: its polygon, its host cell and its lake geometry."""

    def __init__(self, fid, points, area, centroid):
        self.fid = fid
        self.points = points          # (n, 2) polygon vertices
        self.area = float(area)       # m2, true polygon area
        self.centroid = centroid      # (x, y)
        self.cell = None              # host (i, j): the EMBEDDEDV connection
        self.cells = []               # footprint: cells whose centre is inside
        self.klay = 0
        self.bottom = None            # lake bed elevation
        self.rim = None               # spill elevation (land surface)
        self.strt = None
        self.on_channel = False
        self.inlet_reaches = []       # SFR reaches handing their flow to the lake
        self.outlet_reaches = []      # SFR reaches the lake spills into

    @property
    def depth(self):
        return None if self.rim is None else self.rim - self.bottom

    def __repr__(self):
        return ('<PondLake fid=%s area=%.0f m2 cell=%s%s>'
                % (self.fid, self.area, self.cell,
                   ' on-channel' if self.on_channel else ''))


# --------------------------------------------------------------------- #
# geometry
# --------------------------------------------------------------------- #

def _ring_area_centroid(pts):
    """Signed-area centroid of a closed polygon ring.

    Taken relative to the first vertex: on absolute UTM coordinates the
    shoelace terms are ~1e12 and cancel (2026-09-24).
    """
    p = np.asarray(pts, dtype=float)
    if len(p) > 1 and not np.allclose(p[0], p[-1]):
        p = np.vstack([p, p[:1]])
    x0, y0 = p[0]
    x, y = p[:, 0] - x0, p[:, 1] - y0
    cross = x[:-1] * y[1:] - x[1:] * y[:-1]
    a2 = cross.sum()
    area = abs(a2) * 0.5
    if abs(a2) < 1e-12:                       # degenerate ring
        return area, (float(x0 + x.mean()), float(y0 + y.mean()))
    cx = float(((x[:-1] + x[1:]) * cross).sum() / (3.0 * a2))
    cy = float(((y[:-1] + y[1:]) * cross).sum() / (3.0 * a2))
    return area, (x0 + cx, y0 + cy)


def read_pond_polygons(path, id_field='id'):
    """Read pond polygons into :class:`PondLake` objects.

    Takes either tier: the derived ``inputPONDS.geojson`` the converter writes
    -- the one a run should use, since no shapefile belongs on the model path
    -- or a shapefile directly, which is what the GIS folder holds. Both
    arrive through ``marmites_vector.Layer``, so this does not branch.

    (Before WP1d this took a shapefile only, while ``[lak] source`` had been
    pointed at ``inputPONDS.csv``, so enabling LAK from the configuration
    failed on a missing .dbf.)
    """
    import os
    import sys

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if here not in sys.path:
        sys.path.insert(0, here)
    from marmites_vector import Layer, VectorError

    try:
        lay = Layer(path)
    except VectorError as exc:
        raise ValueError(str(exc)) from exc
    if lay.kind != 'polygon':
        raise ValueError('%s is a %s layer; the ponds must be polygons'
                         % (path, lay.kind))
    ids = lay.column(id_field) if id_field in lay.fields else None
    ponds = []
    for k in range(len(lay)):
        rings = [np.asarray(r, dtype=float) for r in lay.rings(k)
                 if len(r) >= 3]
        if not rings:
            continue
        # The outer ring is the largest by area; any interior ring is not
        # lake surface.
        pts = max(rings, key=lambda r: abs(_ring_area_centroid(r)[0]))
        area, cen = _ring_area_centroid(pts)
        ponds.append(PondLake(ids[k] if ids else k, pts, area, cen))
    if not ponds:
        raise ValueError('no polygons found in %s' % path)
    return ponds


def _inside(px, py, ring):
    """Even-odd point-in-polygon for many points against one ring."""
    r = np.asarray(ring, dtype=float)
    if len(r) > 1 and np.allclose(r[0], r[-1]):
        r = r[:-1]
    x0, y0 = r[0]
    x, y = r[:, 0] - x0, r[:, 1] - y0          # relative: UTM precision
    px = np.asarray(px, dtype=float) - x0
    py = np.asarray(py, dtype=float) - y0
    x2, y2 = np.roll(x, -1), np.roll(y, -1)
    inside = np.zeros(px.shape, dtype=bool)
    for a, b, c, d in zip(x, y, x2, y2):
        cross = (b > py) != (d > py)
        with np.errstate(divide='ignore', invalid='ignore'):
            xi = a + (py - b) * (c - a) / (d - b)
        inside ^= cross & (px < xi)
    return inside


def _cell_centres(grid, cells):
    """Area centroids of some cells of a ``marmites_vector.TargetGrid``."""
    out = []
    for ic in cells:
        _a, c = _ring_area_centroid(grid.polygons[ic])
        out.append(c)
    return np.asarray(out, dtype=float).reshape(-1, 2)


def pond_footprints(ponds, grid, active=None, verbose=True):
    """Give every pond its FOOTPRINT on the model grid and one HOST cell.

    The CdL design (cdl_gwf_model_fable_v2 §5b, 45 years converged): a pond
    owns the cells whose centre lies inside its polygon, and its one
    EMBEDDEDV connection goes to the cell holding its centroid -- on a mesh
    seeded at the pond centroids, the seeded cell. The footprint is what the
    stream is cut out of, so that it runs THROUGH the lake instead of past
    it. A pond on a coarse grid contains no cell centre at all (two La Mata
    ponds on the 50 m grid), so the host always belongs to the footprint.

    ``grid`` is a ``marmites_vector.TargetGrid`` -- structured or a mesh, so
    this never touches delr/delc, which on a projected mesh describe a
    (ncpl, 1) proxy grid of 1 m squares. ``active`` is a boolean mask in the
    grid's own shape.

    A pond whose centroid cell is inactive is hosted by the nearest active
    cell, provided the polygon reaches the active domain at all; a pond
    wholly outside it is not part of the model and is left out, reported.
    Sets ``p.cell`` (host, as (i, j)) and ``p.cells`` (footprint) and
    returns the ponds kept.
    """
    act = (np.ones(grid.ncell, dtype=bool) if active is None
           else np.asarray(active, dtype=bool).reshape(-1))
    if act.size != grid.ncell:
        raise ValueError('active mask holds %d cells, the grid %d'
                         % (act.size, grid.ncell))
    kept, dropped, used = [], [], {}
    for p in ponds:
        pts = np.asarray(p.points, dtype=float)
        bbox = (pts[:, 0].min(), pts[:, 1].min(),
                pts[:, 0].max(), pts[:, 1].max())
        cand = [ic for ic in grid.candidates(bbox) if act[ic]]
        foot = []
        if cand:
            cc = _cell_centres(grid, cand)
            foot = [ic for ic, ok in zip(cand, _inside(cc[:, 0], cc[:, 1], pts))
                    if ok]
        cx, cy = p.centroid
        host = grid.cell_containing(cx, cy)
        if host < 0 or not act[host]:
            touches = bool(foot) or any(
                (lambda c: c >= 0 and act[c])(grid.cell_containing(x, y))
                for x, y in pts)
            if not touches:
                dropped.append(p.fid)
                continue
            pool = foot or list(np.where(act)[0])
            cc = _cell_centres(grid, pool)
            host = int(pool[int(np.argmin((cc[:, 0] - cx) ** 2
                                          + (cc[:, 1] - cy) ** 2))])
        if host not in foot:
            foot.append(host)
        p.cell = tuple(int(v) for v in grid.unravel(host))
        p.cells = [tuple(int(v) for v in grid.unravel(ic)) for ic in sorted(foot)]
        used.setdefault(p.cell, []).append(p)
        kept.append(p)

    # An EMBEDDEDV lake must be the only lake connection in its cell, so two
    # ponds sharing a host would be an invalid model. Report it loudly rather
    # than letting MF6 fail with an opaque message.
    clashes = {c: [q.fid for q in v] for c, v in used.items() if len(v) > 1}
    if clashes:
        raise ValueError(
            'ponds share a host cell, which EMBEDDEDV forbids: %s. Refine the '
            'grid around them (grid.voronoi.refine_ponds) or merge the '
            'polygons.' % clashes)
    if verbose:
        if kept:
            nf = [len(p.cells) for p in kept]
            print('LAK: %d pond(s) on the grid, areas %.0f-%.0f m2; footprints '
                  '%d-%d cell(s), hosts %.0f-%.0f m2'
                  % (len(kept), min(p.area for p in kept),
                     max(p.area for p in kept), min(nf), max(nf),
                     min(grid.area[grid.shape[1] * p.cell[0] + p.cell[1]]
                         for p in kept),
                     max(grid.area[grid.shape[1] * p.cell[0] + p.cell[1]]
                         for p in kept)))
        if dropped:
            print('LAK: pond(s) %s lie wholly outside the active domain and '
                  'are not lakes of this model'
                  % ', '.join(str(f) for f in dropped))
    return kept


def assign_pond_cells(ponds, xll, yll, delr, delc, nrow, ncol, idomain=None,
                      verbose=True):
    """:func:`pond_footprints` on a structured grid given by its spacings.

    ``idomain`` may be (nrow, ncol) or (nlay, nrow, ncol); a cell is active
    when any layer is.
    """
    import os
    import sys

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if here not in sys.path:
        sys.path.insert(0, here)
    from marmites_vector import TargetGrid

    grid = TargetGrid.structured(delr, delc, xll, yll)
    active = None
    if idomain is not None:
        ib = np.abs(np.asarray(idomain))
        active = (ib.max(axis=0) > 0) if ib.ndim == 3 else (ib > 0)
    return pond_footprints(ponds, grid, active=active, verbose=verbose)


def lake_table(bottom, rim, area, nsteps=8, headroom=2.0, min_frac=0.01):
    """Stage / volume / sarea / barea rows for one embedded lake.

    A wedge bathymetry (area growing from ~0 at the bed to the full footprint
    at the rim) rather than a flat bottom, so the lake dries smoothly instead
    of stepping its wetted area to zero.
    """
    D = max(float(rim) - float(bottom), 0.1)
    amin = max(float(area) * float(min_frac), 1.0)
    stages = [bottom + D * k / nsteps for k in range(nsteps + 1)] + [rim + headroom]
    rows, vol, ps, pa = [], 0.0, None, None
    for s in stages:
        a = area if s >= rim else amin + (area - amin) * (s - bottom) / D
        if ps is not None:
            vol += 0.5 * (a + pa) * (s - ps)     # trapezoidal integral of sarea
        rows.append((round(float(s), 4), round(float(vol), 6),
                     round(float(a), 4), round(float(a), 4)))  # barea = sarea
        ps, pa = s, a
    return rows
