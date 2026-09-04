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

Table shape
-----------
An excavated pond has a flat bottom, so the wetted area would jump from 0 to
the full footprint as a step and a bone-dry lake degenerates numerically. The
table therefore assumes a wedge bathymetry: the surface area grows linearly
from ~0 at the deepest point to the full polygon area at the rim, which lets
storage and fluxes vanish smoothly as the pond dries. ``barea`` (the bed
exchange area) is set equal to ``sarea``: the wetted bed IS the exchange area.

Evaporation is left to MARMITES (E_ow in Eq. 1 of the paper), so the LAK
package is given no evaporation of its own -- otherwise the pond would
evaporate twice.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import numpy as np

__all__ = ['read_pond_polygons', 'assign_pond_cells', 'lake_table', 'PondLake']

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
        self.cell = None              # host (i, j)
        self.klay = 0
        self.bottom = None            # lake bed elevation
        self.rim = None               # spill elevation (land surface)
        self.strt = None
        self.on_channel = False
        self.inlet_reach = None       # SFR reach handing flow to the lake
        self.outlet_reach = None      # SFR reach the lake spills into

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
    """Signed-area centroid of a closed polygon ring."""
    p = np.asarray(pts, dtype=float)
    if len(p) > 1 and not np.allclose(p[0], p[-1]):
        p = np.vstack([p, p[:1]])
    x, y = p[:, 0], p[:, 1]
    cross = x[:-1] * y[1:] - x[1:] * y[:-1]
    a2 = cross.sum()
    area = abs(a2) * 0.5
    if abs(a2) < 1e-12:                       # degenerate ring
        return area, (float(x.mean()), float(y.mean()))
    cx = float(((x[:-1] + x[1:]) * cross).sum() / (3.0 * a2))
    cy = float(((y[:-1] + y[1:]) * cross).sum() / (3.0 * a2))
    return area, (cx, cy)


def read_pond_polygons(shp_fn, id_field='id'):
    """Read pond polygons from an ESRI shapefile into :class:`PondLake` objects."""
    try:
        import shapefile                       # pyshp
    except ImportError:                        # pragma: no cover
        raise ImportError(
            'reading %s needs pyshp; install it with "pip install pyshp" '
            '(flopy ships it as an optional dependency)' % shp_fn)
    r = shapefile.Reader(shp_fn)
    names = [f[0] for f in r.fields[1:]]
    ponds = []
    for k, sr in enumerate(r.iterShapeRecords()):
        pts = np.asarray(sr.shape.points, dtype=float)
        if len(pts) < 3:
            continue
        area, cen = _ring_area_centroid(pts)
        fid = sr.record[names.index(id_field)] if id_field in names else k
        ponds.append(PondLake(fid, pts, area, cen))
    if not ponds:
        raise ValueError('no polygons found in %s' % shp_fn)
    return ponds


def assign_pond_cells(ponds, xll, yll, delr, delc, nrow, ncol, idomain=None,
                      verbose=True):
    """Give every pond a single host cell (i, j).

    The host is the cell containing the pond centroid. That rule is used rather
    than "cells whose centre falls inside the polygon" because two La Mata
    ponds are small enough to contain no cell centre at all; a centroid always
    lands somewhere.
    """
    delr = np.asarray(delr, dtype=float)       # column widths (x)
    delc = np.asarray(delc, dtype=float)       # row heights (y)
    xedge = np.concatenate([[0.0], np.cumsum(delr)]) + float(xll)
    # rows run north -> south, so the top edge is yll + sum(delc)
    ytop = float(yll) + float(delc.sum())
    yedge = ytop - np.concatenate([[0.0], np.cumsum(delc)])

    used = {}
    for p in ponds:
        cx, cy = p.centroid
        j = int(np.clip(np.searchsorted(xedge, cx) - 1, 0, ncol - 1))
        i = int(np.clip(np.searchsorted(-yedge, -cy) - 1, 0, nrow - 1))
        if idomain is not None:
            ib = np.asarray(idomain)
            act = (np.abs(ib).max(axis=0) > 0) if ib.ndim == 3 else (np.abs(ib) > 0)
            if not act[i, j]:                  # nearest active cell instead
                ii, jj = np.where(act)
                xc = 0.5 * (xedge[:-1] + xedge[1:])
                yc = 0.5 * (yedge[:-1] + yedge[1:])
                d2 = (xc[jj] - cx) ** 2 + (yc[ii] - cy) ** 2
                n = int(np.argmin(d2))
                i, j = int(ii[n]), int(jj[n])
        p.cell = (i, j)
        used.setdefault((i, j), []).append(p)

    # An EMBEDDEDV lake must be the only lake connection in its cell, so two
    # ponds sharing a host would be an invalid model. Report it loudly rather
    # than letting MF6 fail with an opaque message.
    clashes = {c: [q.fid for q in v] for c, v in used.items() if len(v) > 1}
    if clashes:
        raise ValueError(
            'ponds share a host cell, which EMBEDDEDV forbids: %s. Refine the '
            'grid around them (--grid disv) or merge the polygons.' % clashes)
    if verbose:
        print('LAK: %d pond(s) assigned to host cells; areas %.0f-%.0f m2 '
              '(cell %.0f m2)'
              % (len(ponds), min(p.area for p in ponds),
                 max(p.area for p in ponds), float(delr[0] * delc[0])))
    return ponds


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
