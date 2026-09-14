# -*- coding: utf-8 -*-
"""Mesh producers: one function, selected by ``[grid] kind``.  WP1c.1.

    structured / disv   today's raster re-expressed as DISV
                        (``disv_from_structured``). KEPT PERMANENTLY: it is
                        the regression anchor and the MODFLOW-NWT comparison.
    quadtree            USGS GRIDGEN refinement of the same raster.
    voronoi             Triangle + flopy VoronoiGrid.       <-- WP1c.2

Every producer returns the SAME thing -- a MODFLOW 6 DISV ``gridprops`` dict --
so nothing downstream knows or cares which one ran. ``marmites_mesh`` then
projects the model onto it.

Caching
-------
A mesh takes seconds to minutes to build and is reused across runs, which is
exactly the situation that produced the CdL grid-design bug: a cached grid from
one design silently served a run configured for another. So the cache carries a
SIGNATURE of the inputs that determine the mesh, and a mismatch rebuilds rather
than warns. ``build_mesh(force=True)`` -- panel 1's Create grid -- forces it.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import hashlib
import json
import os

import numpy as np

__all__ = ['build_mesh', 'MeshBuildError', 'stream_lines', 'mesh_signature',
           'watershed_ring', 'normalise', 'cell_size_report', 'model_rectangle',
           'dataset_rectangle', 'rectangle_check', 'pond_rings', 'pond_seeds',
           'pond_cells', 'stream_cells', 'cells_touch', 'cell_polygons',
           'cell_centres', 'PRODUCER_VERSION']

# Identifies the MESH-PRODUCING BEHAVIOUR, not the module. Bump it on any
# change that would give a different mesh for the same configuration; it is
# part of the cache signature.
PRODUCER_VERSION = 5


class MeshBuildError(Exception):
    """Raised when a mesh cannot be produced."""


# --------------------------------------------------------------------- #
# refinement features
# --------------------------------------------------------------------- #

def stream_lines(csv_path):
    """Stream centre-lines as ``[[(x, y), ...], ...]``, one list per segment.

    Reads the WP1 ``inputSTREAM.csv`` (seg_id, seq, x, y), which is stored in
    model CRS and deliberately GRID-INDEPENDENT -- the same file feeds the
    quadtree refinement here, the Voronoi seeding of WP1c.2 and the SFR reach
    table of WP3.
    """
    if not os.path.exists(csv_path):
        raise MeshBuildError('stream table not found: %s' % csv_path)
    segs = {}
    with open(csv_path, encoding='utf-8') as fh:
        header = None
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if header is None:
                header = [c.strip() for c in line.split(',')]
                continue
            parts = line.split(',')
            rec = dict(zip(header, parts))
            segs.setdefault(int(rec['seg_id']), []).append(
                (int(rec['seq']), float(rec['x']), float(rec['y'])))
    out = []
    for sid in sorted(segs):
        pts = [(x, y) for _s, x, y in sorted(segs[sid])]
        if len(pts) >= 2:
            out.append(pts)
    if not out:
        raise MeshBuildError('no usable stream segment in %s' % csv_path)
    return out


# --------------------------------------------------------------------- #
# signature / cache
# --------------------------------------------------------------------- #

def mesh_signature(cfg, cMF):
    """Hash of everything that determines the mesh.

    Deliberately NARROW: the source grid geometry and the ``[grid]`` block.
    Changing the number of stress periods must not invalidate a mesh; changing
    the cell size must.
    """
    import dataclasses
    payload = {
        # Bump PRODUCER_VERSION whenever a producer would build a DIFFERENT
        # mesh from the same configuration -- retuning _TRI_AREA_PER_CELL did
        # exactly that, and without this the cache would have served the old
        # mesh forever. This is the CdL grid-design cache bug in miniature.
        'producer_version': PRODUCER_VERSION,
        'kind': cfg.grid_kind,
        'nrow': int(cMF.nrow), 'ncol': int(cMF.ncol), 'nlay': int(cMF.nlay),
        'delr': [float(x) for x in np.asarray(cMF.delr).ravel()],
        'delc': [float(x) for x in np.asarray(cMF.delc).ravel()],
        'xll': float(getattr(cMF, 'xllcorner', 0.0)),
        'yll': float(getattr(cMF, 'yllcorner', 0.0)),
        'voronoi': dataclasses.asdict(cfg.grid.voronoi),
        # The quadtree block was missing here until WP1d: changing
        # refine_level produces a genuinely different mesh, and without it in
        # the signature the cache served the old one. Same bug as the CdL
        # grid-design cache, one block further down.
        'quadtree': dataclasses.asdict(cfg.grid.quadtree),
    }
    blob = json.dumps(payload, sort_keys=True, default=str).encode('utf-8')
    return hashlib.sha256(blob).hexdigest()[:16]


def _cache_paths(cache_dir, kind):
    return (os.path.join(cache_dir, 'mesh_%s.json' % kind),
            os.path.join(cache_dir, 'mesh_%s.sig.json' % kind))


def _load_cached(cache_dir, kind, signature):
    grid_fn, sig_fn = _cache_paths(cache_dir, kind)
    if not (os.path.exists(grid_fn) and os.path.exists(sig_fn)):
        return None
    try:
        with open(sig_fn, encoding='utf-8') as fh:
            if json.load(fh).get('signature') != signature:
                return None
        with open(grid_fn, encoding='utf-8') as fh:
            gp = json.load(fh)
    except (ValueError, OSError):
        return None
    gp['vertices'] = [[int(v[0]), float(v[1]), float(v[2])] for v in gp['vertices']]
    gp['cell2d'] = [[int(r[0]), float(r[1]), float(r[2]), int(r[3])]
                    + [int(x) for x in r[4:]] for r in gp['cell2d']]
    for k in ('top', 'botm'):
        if k in gp and gp[k] is not None:
            gp[k] = np.asarray(gp[k], dtype=float)
    return gp


def normalise(gridprops):
    """Plain-Python vertices/cell2d, whatever the producer returned.

    flopy's VoronoiGrid hands back a mix of ``np.int64``/``np.float64`` and
    Python scalars. Normalising here rather than at the cache boundary matters
    for more than JSON: a fresh mesh and the same mesh reloaded from cache must
    be the SAME object, or a bug can appear on the second run only.
    """
    gp = dict(gridprops)
    gp['vertices'] = [[int(v[0]), float(v[1]), float(v[2])]
                      for v in gp['vertices']]
    gp['cell2d'] = [[int(r[0]), float(r[1]), float(r[2]), int(r[3])]
                    + [int(x) for x in r[4:]] for r in gp['cell2d']]
    for k in ('top', 'botm'):
        if gp.get(k) is not None:
            gp[k] = np.asarray(gp[k], dtype=float)
    for k in ('ncpl', 'nlay', 'nvert'):
        if gp.get(k) is not None:
            gp[k] = int(gp[k])
    return gp


def _save_cached(cache_dir, kind, signature, gridprops):
    os.makedirs(cache_dir, exist_ok=True)
    grid_fn, sig_fn = _cache_paths(cache_dir, kind)
    out = {}
    for k, v in gridprops.items():
        out[k] = v.tolist() if isinstance(v, np.ndarray) else v
    with open(grid_fn, 'w', encoding='utf-8') as fh:
        json.dump(out, fh)
    with open(sig_fn, 'w', encoding='utf-8') as fh:
        json.dump({'signature': signature, 'kind': kind,
                   'ncpl': int(gridprops.get('ncpl', len(gridprops['cell2d'])))},
                  fh, indent=2)


# --------------------------------------------------------------------- #
# producers
# --------------------------------------------------------------------- #

def _produce_structured(cfg, cMF, **_kw):
    from marmites_grid import disv_from_structured
    verts, cell2d, ncpl = disv_from_structured(
        cMF.delr, cMF.delc, getattr(cMF, 'xllcorner', 0.0),
        getattr(cMF, 'yllcorner', 0.0))
    return {'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl,
            'nlay': int(cMF.nlay)}


def _has_pyshp():
    """flopy's gridgen wrapper writes refinement features via pyshp."""
    import importlib.util
    return importlib.util.find_spec('shapefile') is not None


def _produce_quadtree(cfg, cMF, dataset_dir=None, model_ws=None, warn=None,
                      **_kw):
    """GRIDGEN quadtree, refined along the stream network and at the ponds.

    Refining on the streams rather than uniformly is what makes this a real
    end-to-end proof for WP1c.1: the mesh genuinely has no (row, column), so
    every projection and index path is exercised, not bypassed.

    The ponds go in as POLYGON features. GRIDGEN splits every cell a feature
    touches, so a footprint refines the cells it covers whether or not a
    mapped stream runs through it -- which is the point, since a charca off
    the network was being meshed at the background size. Unlike the Voronoi
    case there is nothing to reconcile: a quadtree refines by subdivision,
    features do not have to nest, and two features over the same cell simply
    take the deeper level.
    """
    import mm_paths
    from marmites_gridgen import build_quadtree
    exe = getattr(mm_paths, 'GRIDGEN_EXE', None)
    if not exe or not os.path.exists(exe):
        raise MeshBuildError(
            'gridgen executable not found (%s). Set $MM_GRIDGEN_EXE, or use '
            "grid.kind = 'structured'." % exe)
    feats = []
    csv = os.path.join(dataset_dir, 'inputSTREAM.csv') if dataset_dir else None
    if csv and not os.path.exists(csv):
        if warn:
            warn('no inputSTREAM.csv in %s: building an UNREFINED quadtree. '
                 'Run code/tools/gis_to_dataset.py first.' % dataset_dir)
    elif csv and not _has_pyshp():
        # flopy writes refinement features through a shapefile, and only pyshp
        # can write one for it. Degrade to an unrefined mesh with the reason
        # said out loud: an unrefined quadtree is geometrically the base grid,
        # so silently producing one would look like success.
        if warn:
            warn('pyshp is not installed, so the stream refinement is SKIPPED '
                 'and this quadtree is geometrically the base grid. Install it '
                 'with:  conda install -p C:\\miniconda3\\envs\\flopy pyshp')
    elif csv and not getattr(getattr(cfg.grid, 'quadtree', None),
                             'refine_streams', True):
        if warn:
            warn('grid.quadtree.refine_streams is off: building an UNREFINED '
                 'quadtree, which is geometrically the base grid.')
    elif csv:
        level = int(getattr(getattr(cfg.grid, 'quadtree', None),
                            'refine_level', 2))
        feats.append((stream_lines(csv), 'line', level))

    q = getattr(cfg.grid, 'quadtree', None)
    if dataset_dir and getattr(q, 'refine_ponds', False):
        rings = pond_rings(dataset_dir)
        if not rings:
            if warn:
                warn('grid.quadtree.refine_ponds is on but no pond footprint '
                     'was read from inputPONDS.geojson, so the ponds are NOT '
                     'refined. Run code/tools/gis_to_dataset.py.')
        elif not _has_pyshp():
            if warn and not csv:    # said once already when there is a csv
                warn('pyshp is not installed, so the pond refinement is '
                     'SKIPPED. Install it with:  conda install -p '
                     'C:\\miniconda3\\envs\\flopy pyshp')
        else:
            # 0 means "the same level as the streams", so the ponds follow
            # them unless there is a reason to split further.
            lev = int(getattr(q, 'pond_level', 0) or
                      getattr(q, 'refine_level', 2))
            # flopy reads a polygon as a list of RINGS, so each footprint is
            # wrapped in one, and pyshp appends a list to whatever it is
            # given -- so the ring has to arrive CLOSED and as a list, or it
            # fails on a tuple it cannot concatenate.
            feats.append(([[list(r) + [r[0]]] for r in rings],
                          'polygon', lev))
            if warn:
                warn('grid.quadtree: %d pond footprint(s) refined to level %d.'
                     % (len(rings), lev))
    ws = model_ws or os.path.join(os.getcwd(), '_gridgen')
    return build_quadtree(cMF, exe, ws, refine_features=feats,
                          layers=list(range(int(cMF.nlay))))


def rectangle_cell_size(cfg):
    """The cell size the model rectangle is snapped to, per grid kind.

    ``grid.cell_size`` for structured, disv and quadtree -- where it IS the
    cell, or the background GRIDGEN halves. ``grid.voronoi.cell_far`` for
    voronoi, where ``cell_size`` plays no part in the cells at all.
    """
    if cfg.grid_kind == 'voronoi':
        return float(cfg.grid.voronoi.cell_far)
    return float(cfg.grid.cell_size)


def model_rectangle(cfg, bbox=None):
    """The grid rectangle for this configuration (WP1d, panel 1).

    Returns ``(nrow, ncol, delr, delc, xllcorner, yllcorner)`` -- the six
    numbers every producer needs, and the only ones the mesh signature is
    keyed on besides the ``[grid]`` block itself.

    Derived from the CATCHMENT POLYGON's bounding box, which is the whole
    point of panel 1: the domain comes first and the grid is built inside it.
    The origin is snapped DOWN to a multiple of the cell size and the far edge
    UP, so the same polygon and cell size always give the same grid -- a
    rectangle that shifted with a re-exported shapefile would invalidate every
    cached mesh for no reason.

    ``grid.override.enable`` bypasses all of it and returns the origin and
    shape verbatim, which is how an existing grid is reproduced exactly.

    WP1d panel 1, D1 of §2A.7: on a VORONOI grid the rectangle is snapped to
    ``grid.voronoi.cell_far`` rather than ``grid.cell_size``. The Voronoi
    producer never uses ``cell_size`` for its cells -- ``cell_far`` is the
    size -- so leaving the snap on ``cell_size`` meant one number that did
    nothing visible and another that did the work.
    """
    import math

    ov = cfg.grid.override
    if ov.enable:
        if ov.nrow <= 0 or ov.ncol <= 0:
            raise MeshBuildError('grid.override needs nrow and ncol > 0')
        cs = float(cfg.grid.cell_size)
        return (int(ov.nrow), int(ov.ncol),
                np.full(int(ov.ncol), cs), np.full(int(ov.nrow), cs),
                float(ov.xllcorner), float(ov.yllcorner))
    if bbox is None:
        raise MeshBuildError(
            'no catchment bounding box: either set grid.override.enable, or '
            'give the polygon named by grid.boundary')
    cs = rectangle_cell_size(cfg)
    if cs <= 0:
        raise MeshBuildError('%s must be > 0'
                             % ('grid.voronoi.cell_far'
                                if cfg.grid_kind == 'voronoi'
                                else 'grid.cell_size'))
    b = float(cfg.grid.buffer)
    x0, y0, x1, y1 = (float(bbox[0]) - b, float(bbox[1]) - b,
                      float(bbox[2]) + b, float(bbox[3]) + b)
    if x1 <= x0 or y1 <= y0:
        raise MeshBuildError('the catchment bounding box is empty: %s' % (bbox,))
    x0 = math.floor(x0 / cs) * cs
    y0 = math.floor(y0 / cs) * cs
    ncol = int(math.ceil((x1 - x0) / cs))
    nrow = int(math.ceil((y1 - y0) / cs))
    return (nrow, ncol, np.full(ncol, cs), np.full(nrow, cs), x0, y0)


# --------------------------------------------------------------------- #
# the rectangle the dataset already carries                     WP1d, 1e
# --------------------------------------------------------------------- #
#
# STOPGAP, and it should stay one. A model is assembled from rasters written
# on ONE lattice -- soil zones, vegetation areas, and on the MODFLOW side
# elev, thickness, hk, ibound. `marmites_mesh.project_model` then resamples
# that assembled model onto the mesh, so the mesh and those rasters have to
# stand on the same ground. Panel 1 derives its rectangle from the catchment
# polygon; the committed rasters were frozen on whatever rectangle was
# current when somebody exported them. When the two disagree the mesh hangs
# over ground the rasters do not cover, and the projection fails late and
# obscurely -- top averaged over one subset of source cells and botm over
# another, until MF6's top > botm is violated.
#
# The cure is for the model's own rasters to be re-derived from the
# cartography onto whatever rectangle this panel produces. That belongs to
# the MODEL panel and is not built yet. Until it is, this says so HERE, where
# the rectangle is chosen, instead of leaving it to surface as a traceback.
#
# The DEM is excluded on purpose: it is kept at its own resolution and
# wrapped onto the cells at run time, so it is SUPPOSED to be on its own
# lattice (see marmites_dem).
RECT_EXCLUDE = ('inputdem.asc',)


def _rect_of(head):
    """The five numbers that place a raster: origin, shape, cell."""
    return (round(float(head['xllcorner']), 3),
            round(float(head['yllcorner']), 3),
            int(head['nrows']), int(head['ncols']),
            round(float(head['cellsize']), 6))


def _asc_files(dataset_dir, depth=2):
    """Every ESRI ASCII grid in the dataset tree.

    Two levels deep because the MODFLOW rasters live in a workspace
    subdirectory next to the MARMITES tables, and the model reads both.
    """
    out = []
    root = str(dataset_dir)
    for dirpath, dirnames, filenames in os.walk(root):
        rel = os.path.relpath(dirpath, root)
        level = 0 if rel == os.curdir else rel.count(os.sep) + 1
        if level >= depth:
            dirnames[:] = []
        for name in filenames:
            if name.lower().endswith('.asc') \
                    and name.lower() not in RECT_EXCLUDE:
                out.append(os.path.join(dirpath, name))
    return sorted(out)


def dataset_rectangle(dataset_dir):
    """The rectangle the model's rasters occupy.  ``(rect, names, others)``.

    ``rect`` is ``(xll, yll, nrow, ncol, cellsize)``, or None when the dataset
    holds no raster to compare against -- a brand new catchment, where there
    is nothing to disagree with yet. Rasters are GROUPED by the rectangle they
    declare and the largest group wins; ``others`` lists the groups that lost,
    because rasters disagreeing AMONG THEMSELVES is its own problem and hiding
    it behind a majority vote would be worse than saying it.
    """
    import marmites_dem as mdem

    groups = {}
    for path in _asc_files(dataset_dir):
        try:
            head = mdem.read_asc_header(path)
        except Exception:
            continue                      # not this function's to explain
        groups.setdefault(_rect_of(head), []).append(
            os.path.relpath(path, str(dataset_dir)))
    if not groups:
        return None, [], []
    ranked = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    (rect, names), others = ranked[0], ranked[1:]
    return rect, sorted(names), [(r, sorted(n)) for r, n in others]


def rectangle_check(cfg, bbox, dataset_dir):
    """Would the grid this panel builds stand on the dataset's rasters?

    Returns a dict whose ``status`` is one of

        'ok'        the derived rectangle is covered by the rasters
        'overhang'  it reaches past them -- the projection will fail
        'shifted'   covered, but the two lattices are not commensurate
        'none'      the dataset has no raster to compare against
        'error'     the rectangle itself cannot be derived

    with ``derived`` and ``source`` rectangles, the overhang per side in
    metres, and the rasters the source rectangle was read from.
    """
    out = {'status': 'none', 'derived': None, 'source': None, 'overhang': {},
           'names': [], 'others': [], 'detail': ''}
    try:
        nrow, ncol, delr, delc, xll, yll = model_rectangle(cfg, bbox)
    except MeshBuildError as exc:
        out['status'] = 'error'
        out['detail'] = str(exc)
        return out
    cs = float(delr[0])
    out['derived'] = (round(float(xll), 3), round(float(yll), 3),
                      int(nrow), int(ncol), round(cs, 6))

    rect, names, others = dataset_rectangle(dataset_dir)
    out['names'], out['others'] = names, others
    if rect is None:
        out['detail'] = ('the dataset carries no raster yet, so there is '
                         'nothing for this grid to disagree with')
        return out
    out['source'] = rect

    sx0, sy0, snr, snc, scs = rect
    sx1, sy1 = sx0 + snc * scs, sy0 + snr * scs
    dx0, dy0 = float(xll), float(yll)
    dx1, dy1 = dx0 + ncol * cs, dy0 + nrow * cs
    over = {'west': max(0.0, sx0 - dx0), 'east': max(0.0, dx1 - sx1),
            'south': max(0.0, sy0 - dy0), 'north': max(0.0, dy1 - sy1)}
    out['overhang'] = dict((k, round(v, 3)) for k, v in over.items()
                           if v > 1e-6)
    if out['overhang']:
        out['status'] = 'overhang'
        out['detail'] = ('the grid reaches %s past the rasters'
                         % ', '.join('%g m %s' % (v, k) for k, v
                                     in sorted(out['overhang'].items())))
        return out
    # Covered, but a derived cell whose edges fall inside a source cell is
    # resampled from fractions everywhere instead of cell on cell. Not fatal
    # -- worth saying, because it is the difference between a projection that
    # reproduces the legacy model and one that blurs it.
    off_x, off_y = (dx0 - sx0) % scs, (dy0 - sy0) % scs
    snapped = (min(off_x, scs - off_x) < 1e-6
               and min(off_y, scs - off_y) < 1e-6)
    if not snapped:
        out['status'] = 'shifted'
        out['detail'] = ('the grid sits inside the rasters but its origin is '
                         'offset %g m east and %g m north of their lattice'
                         % (round(off_x, 3), round(off_y, 3)))
        return out
    out['status'] = 'ok'
    out['detail'] = 'the grid stands entirely on the dataset rasters'
    return out


def grid_stub(cfg, bbox=None, nlay=1):
    """A cMF-shaped object carrying only what the mesh producers read.

    They use ``nrow``, ``ncol``, ``nlay``, ``delr``, ``delc``, the origin and
    -- GRIDGEN only -- ``top`` and ``botm``, so a mesh can be built, and
    previewed, without parsing the MODFLOW parameter file or the time
    discretisation.

    The elevations are PLACEHOLDERS, a unit-thick layer stack. GRIDGEN wants
    a StructuredGrid and a StructuredGrid wants elevations, but only the
    PLAN VIEW of the result is used here: the real top and botm come from the
    rasters when the model is built. Leaving them out is what made the
    quadtree producer fail with "SimpleNamespace has no attribute 'top'".
    """
    from types import SimpleNamespace

    nrow, ncol, delr, delc, xll, yll = model_rectangle(cfg, bbox)
    nlay = int(nlay)
    top = np.zeros((nrow, ncol), dtype=float)
    botm = np.stack([np.full((nrow, ncol), -(k + 1.0)) for k in range(nlay)])
    return SimpleNamespace(nrow=nrow, ncol=ncol, nlay=nlay,
                           delr=delr, delc=delc, top=top, botm=botm,
                           xllcorner=float(xll), yllcorner=float(yll))


def watershed_ring(csv_path):
    """Domain polygon as ``[(x, y), ...]`` from the WP1 ``inputWATERSHED.csv``.

    Only the first ring is used: the Voronoi domain is a single outer
    boundary. A closing point identical to the first is dropped, because
    Triangle wants an open ring and would otherwise see a zero-length segment.
    """
    if not os.path.exists(csv_path):
        raise MeshBuildError('watershed table not found: %s' % csv_path)
    rings = {}
    with open(csv_path, encoding='utf-8') as fh:
        header = None
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if header is None:
                header = [c.strip() for c in line.split(',')]
                continue
            rec = dict(zip(header, line.split(',')))
            rings.setdefault(int(rec['ring_id']), []).append(
                (int(rec['seq']), float(rec['x']), float(rec['y'])))
    if not rings:
        raise MeshBuildError('no ring in %s' % csv_path)
    pts = [(x, y) for _s, x, y in sorted(rings[min(rings)])]
    if len(pts) > 2 and abs(pts[0][0] - pts[-1][0]) < 1e-6 \
            and abs(pts[0][1] - pts[-1][1]) < 1e-6:
        pts = pts[:-1]
    if len(pts) < 3:
        raise MeshBuildError('watershed ring has %d point(s)' % len(pts))
    return pts


# Triangle's `-a` bounds the triangle area; a Voronoi cell is the dual of a
# VERTEX and so covers about two triangles, but Triangle's quality constraints
# make the real triangles smaller than the bound, so the two effects do not
# cancel. Measured on La Mata (cell_far=100 m, uniform): a bound of 0.433*s^2
# produced a mean cell area of 5328 m2, i.e. 1.23x the bound. Hence the
# constant below. It is a CALIBRATION, not an identity -- which is exactly why
# `cell_size_report` states the area actually achieved on every build.
_TRI_AREA_PER_CELL = 0.813


def _max_area_for(cell_size):
    """Triangle maximum triangle area for a target mean cell of ``cell_size``.

    ``cell_size`` is read as the side of the equivalent square, so
    ``cell_far = 100`` aims at a mean cell area of 10 000 m2.
    """
    return _TRI_AREA_PER_CELL * float(cell_size) ** 2


def _grid_rectangle(cMF):
    """The structured model's own extent as a ring, model CRS."""
    x0 = float(getattr(cMF, 'xllcorner', 0.0))
    y0 = float(getattr(cMF, 'yllcorner', 0.0))
    x1 = x0 + float(np.sum(np.asarray(cMF.delr, dtype=float)))
    y1 = y0 + float(np.sum(np.asarray(cMF.delc, dtype=float)))
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def _clip_to_grid(ring, cMF, warn=None):
    """Intersect the domain with the model grid extent.

    La Mata's ``Limite.shp`` is a 19.2 km2 rectangle while the MODFLOW grid
    covers 9.75 km2, so meshing the raw boundary spends two thirds of the
    cells outside the model, where the projection can only mark them inactive.
    Clipping is not cosmetic: it is the difference between 3584 cells with 893
    active and a mesh that actually resolves the catchment.
    """
    rect = _grid_rectangle(cMF)
    try:
        from shapely.geometry import Polygon
    except ImportError:
        if warn:
            warn('shapely is not available, so the domain cannot be clipped to '
                 'the model grid; using the grid extent itself.')
        return rect
    dom = Polygon(ring)
    grid = Polygon(rect)
    if not dom.is_valid:
        dom = dom.buffer(0)
    inter = dom.intersection(grid)
    if inter.is_empty:
        raise MeshBuildError(
            'the catchment boundary and the model grid do not overlap. Check '
            'that inputWATERSHED.csv and the .ini share a CRS.')
    if inter.geom_type == 'MultiPolygon':
        inter = max(inter.geoms, key=lambda g: g.area)
    if warn and inter.area < dom.area * 0.999:
        warn('catchment boundary clipped to the model grid extent: '
             '%.2f km2 of %.2f km2 kept.'
             % (inter.area / 1e6, dom.area / 1e6))
    pts = list(inter.exterior.coords)
    if len(pts) > 2 and pts[0] == pts[-1]:
        pts = pts[:-1]
    return [(float(x), float(y)) for x, y in pts]


def pond_rings(dataset_dir):
    """Pond footprints as ``[[(x, y), ...], ...]`` from ``inputPONDS.geojson``.

    The GeoJSON is what the converter writes for the outlines --
    ``inputPONDS.csv`` carries only the centroid, the area and the DEM
    statistics, which is enough to place a cell but not to draw a footprint or
    to size one. Both are in model CRS and grid-independent.

    Read with ``json`` alone: this is the model path, where the cookbook's D3
    rule keeps geopandas out. An unreadable or absent file is not an error --
    a catchment with no ponds is normal -- it is an empty list.
    """
    import json

    path = os.path.join(str(dataset_dir), 'inputPONDS.geojson')
    if not os.path.exists(path):
        return []
    try:
        with open(path, encoding='utf-8') as fh:
            doc = json.load(fh)
    except (ValueError, OSError):
        return []
    out = []
    for feat in doc.get('features', []):
        geom = feat.get('geometry') or {}
        kind = geom.get('type')
        if kind == 'Polygon':
            parts = [geom.get('coordinates', [[]])]
        elif kind == 'MultiPolygon':
            parts = geom.get('coordinates', [])
        else:
            continue
        for poly in parts:
            if not poly:
                continue
            ring = [(float(x), float(y)) for x, y in poly[0][:2000]]
            if len(ring) > 2 and ring[0] == ring[-1]:
                ring = ring[:-1]
            if len(ring) >= 3:
                out.append(ring)
    return out


def _ring_area_centre(ring):
    """``(area, x, y)`` of a closed ring by the shoelace formula.

    The CENTROID of the polygon, not the mean of its vertices: a rim mapped
    with many points on one side would drag the mean towards that side, and
    the seed has to sit where the cell should be centred.
    """
    a = cx = cy = 0.0
    n = len(ring)
    for i in range(n):
        x0, y0 = ring[i]
        x1, y1 = ring[(i + 1) % n]
        cross = x0 * y1 - x1 * y0
        a += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    if abs(a) < 1e-12:
        xs = [p[0] for p in ring]
        ys = [p[1] for p in ring]
        return 0.0, sum(xs) / n, sum(ys) / n
    return abs(a) / 2.0, cx / (3.0 * a), cy / (3.0 * a)


# A Voronoi cell can NEVER follow a pond outline: a face is the bisector of
# the segment between two generators, so a rim added as a constraint polygon
# is cut by the cells rather than honoured, and every cell along it straddles
# the rim half in and half out. CdL proved that on a 5922-cell rim mesh and
# settled on SEEDING instead (2026-07-04): one generator at the pond centre,
# plus a ring of helpers a little way out. The centre's cell is then bounded
# by the bisectors to its own ring, so it comes out pond-SIZED and
# pond-CENTRED and covers the footprint -- which is what LAK needs.
#
# Without the ring the centre generator simply takes the local background
# cell, because nothing near it competes. 2.5 x the pond radius with 6
# helpers is what CdL runs on.
POND_RING_FACTOR = 2.5
POND_RING_N = 6


def pond_seeds(dataset_dir, inside=None, factor=POND_RING_FACTOR,
               n_ring=POND_RING_N):
    """Generator nodes that give each pond a cell of its own.

    ``inside`` is an optional ``(x, y) -> bool`` test for the domain. A pond
    whose centre falls outside it is NOT meshed and is reported: Triangle
    discards a node outside the boundary polygon, so without this a pond
    beyond the catchment edge simply disappeared -- no cell, no message, and
    the LAK footprint later looking for one. Ring helpers outside are dropped
    on their own, since a generator beyond the boundary would pull a cell out
    of the domain.

    Returns ``(nodes, ponds, outside)``: the generator nodes, the ponds kept
    as ``(ring, area, cx, cy, radius)``, and the rings of those refused.
    """
    import math

    nodes, ponds, outside = [], [], []
    for ring in pond_rings(dataset_dir):
        area, cx, cy = _ring_area_centre(ring)
        if area <= 0.0:
            continue
        if inside is not None and not inside(cx, cy):
            outside.append(ring)
            continue
        r = math.sqrt(area / math.pi)
        ponds.append((ring, area, cx, cy, r))
        nodes.append((cx, cy))
        for k in range(int(n_ring)):
            th = 2.0 * math.pi * k / float(n_ring)
            px = cx + factor * r * math.cos(th)
            py = cy + factor * r * math.sin(th)
            if inside is None or inside(px, py):
                nodes.append((px, py))
    return nodes, ponds, outside


def cell_centres(gridprops):
    """``(ncpl, 2)`` array of cell centres from a DISV gridprops."""
    return np.asarray([(float(rec[1]), float(rec[2]))
                       for rec in gridprops['cell2d']], dtype=float)


def cell_polygons(gridprops):
    """``[[(x, y), ...], ...]``, one ring per cell, in cell2d order."""
    verts = {int(v[0]): (float(v[1]), float(v[2]))
             for v in gridprops['vertices']}
    return [[verts[int(i)] for i in rec[4:]] for rec in gridprops['cell2d']]


def pond_cells(gridprops, rings):
    """The cells each pond footprint owns.  ``[[icell2d, ...], ...]``.

    CENTRE-INSIDE, with the nearest cell as a fallback -- CdL's rule. A
    Voronoi cell belongs to the generator it surrounds, so "the centre is in
    the water" is the honest test of whether a cell is part of the lake; the
    fallback catches a pond smaller than the cell that covers it, which would
    otherwise own nothing at all and leave LAK with an empty footprint.

    This is the LAK footprint of WP4 read ahead of time. Panel 1 draws it so
    that a refinement can be SEEN to have worked, rather than taken on trust.
    """
    if not rings:
        return []
    cen = cell_centres(gridprops)
    out = []
    for ring in rings:
        inside = _inside_ring(ring)
        got = [i for i in range(cen.shape[0]) if inside(cen[i, 0], cen[i, 1])]
        if not got:
            _area, cx, cy = _ring_area_centre(ring)
            got = [int(np.argmin(np.hypot(cen[:, 0] - cx, cen[:, 1] - cy)))]
        out.append(got)
    return out


def stream_cells(gridprops, segments, step=None):
    """The cells a mapped stream runs through.  ``[icell2d, ...]``, sorted.

    By NEAREST CENTRE along a densely sampled centre-line. On a Voronoi mesh
    that is not an approximation -- the cell containing a point IS the cell
    whose generator is nearest -- and on a quadtree it is close enough for
    what this is for, which is showing where the SFR reaches will fall
    against the pond cells beside them.
    """
    if not segments:
        return []
    cen = cell_centres(gridprops)
    if step is None:
        # Half the smallest cell across, so no cell on the line is stepped
        # over. Bounded below, or a mesh with one sliver cell samples for
        # ever.
        polys = cell_polygons(gridprops)
        span = min((max(p[0] for p in poly) - min(p[0] for p in poly))
                   or 1e9 for poly in polys)
        step = max(span / 2.0, 1.0)
    pts = []
    for seg in segments:
        for k in range(len(seg) - 1):
            (x0, y0), (x1, y1) = seg[k], seg[k + 1]
            d = float(np.hypot(x1 - x0, y1 - y0))
            n = max(int(np.ceil(d / step)), 1)
            for j in range(n + 1):
                f = j / float(n)
                pts.append((x0 + f * (x1 - x0), y0 + f * (y1 - y0)))
    if not pts:
        return []
    pts = np.asarray(pts, dtype=float)
    hit = set()
    for a in range(0, pts.shape[0], 2048):
        chunk = pts[a:a + 2048]
        d = ((chunk[:, None, 0] - cen[None, :, 0]) ** 2
             + (chunk[:, None, 1] - cen[None, :, 1]) ** 2)
        hit.update(int(i) for i in np.argmin(d, axis=1))
    return sorted(hit)


def cells_touch(gridprops, a, b, tol=1e-6):
    """Do any cell of ``a`` and any cell of ``b`` share a vertex?

    The question WP4 will ask: a lake that does not touch the network it is
    supposed to drain into cannot be connected to it by a mover, whatever
    the package file says.
    """
    polys = cell_polygons(gridprops)

    def corners(cells):
        return set((round(x / tol), round(y / tol))
                   for i in cells for x, y in polys[i])
    return bool(corners(a) & corners(b))


def _produce_voronoi(cfg, cMF, dataset_dir=None, model_ws=None, warn=None,
                     **_kw):
    """Triangle + flopy VoronoiGrid over the catchment boundary.  WP1c.2.

    The domain comes from the WP1 ``inputWATERSHED.csv`` and the refinement
    from ``inputSTREAM.csv`` -- both grid-independent, both already in the
    repository, so the mesh is reproducible from the dataset alone with no GIS
    file in sight.
    """
    import mm_paths
    from flopy.utils.triangle import Triangle
    from flopy.utils.voronoi import VoronoiGrid

    v = cfg.grid.voronoi
    exe = mm_paths.TRIANGLE_EXE
    if not os.path.exists(exe):
        raise MeshBuildError(
            'triangle executable not found (%s). Set $MM_TRIANGLE_EXE.' % exe)
    if not dataset_dir:
        raise MeshBuildError('the voronoi producer needs the dataset folder')
    ring = watershed_ring(os.path.join(dataset_dir, 'inputWATERSHED.csv'))
    ring = _clip_to_grid(ring, cMF, warn)
    ws = model_ws or os.path.join(os.getcwd(), '_triangle')
    os.makedirs(ws, exist_ok=True)

    nodes, ponds, seeds = None, [], []
    if v.refine_ponds:
        # The seeds go to Triangle as NODES, not as regions or polygons: they
        # are generators, points the triangulation must contain, and it is
        # being a generator that gives the pond its cell.
        # With a pond SIZE asked for, the sizing zone is the ring of
        # generators and the helper ring would only add competitors inside
        # it -- the centre node is then the one thing near the pond.
        seeds, ponds, gone = pond_seeds(
            dataset_dir, inside=_inside_ring(ring),
            n_ring=0 if float(v.cell_pond) > 0.0 else POND_RING_N)
        if warn and gone:
            warn('%d pond(s) fall outside the catchment boundary and are NOT '
                 'meshed: %s.' % (len(gone), ', '.join(
                     '%.0f, %.0f' % _ring_area_centre(g)[1:] for g in gone)))
        if not seeds:
            if warn:
                warn('grid.voronoi.refine_ponds is on but no pond footprint was '
                     'read from inputPONDS.geojson: the mesh has no pond '
                     'cell. Run code/tools/gis_to_dataset.py.')
        else:
            nodes = np.asarray(seeds, dtype=float)

    # WHY maximum_area IS NOT SET when anything is refined.
    #
    # flopy passes it to Triangle as `-a<number>`, and Triangle reads `-a`
    # WITH a number as "this area, everywhere" -- per-region constraints in
    # the .poly are then not read at all. So every band's own maximum area
    # was written into the file and ignored: the grading that came out was
    # whatever the density of the band BOUNDARY segments happened to produce,
    # which is why a 5 m corridor was measuring 200 m2 cells and why
    # cell_pond changed nothing. A bare `-a` reads the regions, so the
    # background has to become a region like any other -- added below, once
    # the bands know where they are.
    refining = _will_refine(cfg, dataset_dir, ponds, warn)
    tri = Triangle(model_ws=ws, exe_name=exe, nodes=nodes,
                   maximum_area=None if refining
                   else _max_area_for(v.cell_far))
    tri.add_polygon(ring)

    if refining:
        # Only the ponds that are actually IN the domain: buffering one that
        # is not would put a band outside the catchment boundary.
        covered = _add_refinement_regions(tri, cfg, dataset_dir, warn,
                                          ponds=ponds)
        _add_background_region(tri, ring, covered, v.cell_far, warn)
    elif ponds and warn:
        warn('grid.voronoi.refine_ponds is on and nothing is refined, so '
             'each pond gets a cell CENTRED on it but at the background size: '
             'the graded sizes the refinement carries (cell_near_stream, '
             'grade_ratio) are what a pond would be refined to.')
    if ponds and warn:
        warn('pond seeding: %d pond(s), %d generator node(s).'
             % (len(ponds), len(nodes)))

    tri.build(verbose=False)
    gp = VoronoiGrid(tri).get_gridprops_vertexgrid()
    gp['nlay'] = int(cMF.nlay)
    return gp


def _inside_ring(ring):
    """``(x, y) -> bool`` for a simple polygon, by the crossing-number rule.

    Written out rather than reached for in shapely: it is a dozen lines, it
    is called a few dozen times, and the producer already degrades to an
    unrefined mesh when shapely is absent -- pond seeding should not be the
    one thing that makes it a hard dependency.
    """
    pts = list(ring)
    n = len(pts)

    def inside(x, y):
        hit = False
        j = n - 1
        for i in range(n):
            xi, yi = pts[i]
            xj, yj = pts[j]
            if (yi > y) != (yj > y):
                xc = xi + (y - yi) * (xj - xi) / (yj - yi)
                if x < xc:
                    hit = not hit
            j = i
        return hit
    return inside


# How far a pond's sizing zone stands off it, as a multiple of the pond
# radius. A Voronoi cell reaches half way to its neighbours, so a ring of
# generators at f x r gives the centre a cell of radius f x r / 2: at f = 2
# the pond's cell comes out at about the pond's own area, which is the
# one-cell-per-pond CdL builds.
POND_ZONE_FACTOR = 2.0
POND_ZONE_N = 12


def _will_refine(cfg, dataset_dir, ponds, warn):
    """Will anything actually be given a region of its own?

    Asked BEFORE Triangle is constructed, because the answer decides whether
    it is given a global maximum area -- and a global area silently disables
    every per-region one.
    """
    v = cfg.grid.voronoi
    if v.refine_ponds and float(v.cell_pond) > 0.0 and ponds:
        pass                       # the pond zones alone are enough
    elif not v.stream_refine:
        return False
    if not os.path.exists(os.path.join(dataset_dir, 'inputSTREAM.csv')) \
            and not (ponds and float(v.cell_pond) > 0.0):
        return False
    try:
        import shapely.geometry                        # noqa: F401
    except ImportError:
        return False
    return bool(v.graded_bands()) or bool(ponds and float(v.cell_pond) > 0.0)


def _add_background_region(tri, ring, covered, cell_far, warn):
    """The unrefined remainder, as a region of its own.

    With the bands carrying their own areas, the rest of the catchment has
    none unless it is seeded too -- and a part of a PSLG with no area
    constraint is meshed as coarsely as the quality rules allow, which on La
    Mata means a handful of enormous triangles.
    """
    from shapely.geometry import Polygon

    dom = Polygon(ring)
    if not dom.is_valid:
        dom = dom.buffer(0)
    rest = dom if covered is None else dom.difference(covered)
    if rest.is_empty:
        if warn:
            warn('grid.voronoi: the refinement covers the whole catchment, so '
                 'there is no background left to size.')
        return 0
    added = 0
    for geom in getattr(rest, 'geoms', [rest]):
        if geom.area <= 0.0:
            continue
        pt = geom.representative_point()
        tri.add_region((pt.x, pt.y), attribute=0,
                       maximum_area=_max_area_for(cell_far))
        added += 1
    return added


def _pond_footprints(rings, Polygon, tol=0.5):
    """The footprints themselves, as shapely polygons."""
    out = []
    for ring in rings:
        poly = Polygon(ring)
        if not poly.is_valid:
            poly = poly.buffer(0)
        poly = poly.simplify(tol)
        if not poly.is_empty and poly.area > 0.0:
            out.append(poly)
    return out


def _pond_zones(ponds, Polygon, factor=POND_ZONE_FACTOR, n=POND_ZONE_N):
    """A sizing zone standing off each pond, as a coarse ring of segments.

    NOT the footprint. Every vertex of a constraint polygon is a point the
    triangulation must contain, and in a Voronoi mesh a point is a GENERATOR:
    a rim added as a zone hands the pond four neighbours of its own, and the
    centre's cell -- bounded by the bisectors to them -- can only be a
    fraction of the water it is supposed to cover. Measured: a 1018 m2 pond
    came back with a 116 m2 cell and 29 cells inside its rim.

    A circle at ``factor`` x the pond radius, with few enough vertices to be
    the only generators near the pond, gives the centre a cell of about the
    pond's own size. What ``cell_pond`` then does is decide whether Triangle
    may put anything INSIDE that circle: at the pond's own width it may not,
    and the pond is one cell.
    """
    import math

    out = []
    for _ring, _area, cx, cy, r in ponds:
        rr = float(factor) * float(r)
        out.append(Polygon([(cx + rr * math.cos(2.0 * math.pi * k / n),
                             cy + rr * math.sin(2.0 * math.pi * k / n))
                            for k in range(int(n))]))
    return out


def _add_pond_zones(tri, zones, size, warn):
    """Each pond as a bounded region at ``size``.  Returns the count.

    Called only once the bands have swallowed the zones, so every rim here is
    strictly inside whatever band contains it and nothing crosses.
    """
    added = 0
    for poly in zones:
        for geom in getattr(poly, 'geoms', [poly]):
            ring = list(geom.exterior.coords)[:-1]
            if len(ring) < 3:
                continue
            tri.add_polygon(ring)
            pt = geom.representative_point()
            tri.add_region((pt.x, pt.y), attribute=100 + added,
                           maximum_area=_max_area_for(size))
            added += 1
    if warn and added:
        warn('grid.voronoi: %d pond zone(s) bounded at %g m.' % (added, size))
    return added


def _add_refinement_regions(tri, cfg, dataset_dir, warn, ponds=()):
    """Refine around the mapped features, graded outward.

    The features are the stream centre-lines and -- when
    ``grid.voronoi.refine_ponds`` is on -- the ponds, each given as
    ``(ring, area, cx, cy, radius)`` by :func:`pond_seeds`.

    With ``cell_pond`` left at 0 the footprints simply JOIN the geometry the
    bands are buffered from, so a pond is meshed at ``cell_near_stream`` and
    grades out with the corridor. One nested family of constraints, nothing
    crossing anything.

    With ``cell_pond`` set, each footprint also becomes a zone of its OWN,
    bounded by its rim and carrying its own maximum area -- the only way to
    make a pond cell coarser than the corridor around it. The rim would cross
    the band boundaries, since La Mata's charcas sit ON the mapped streams,
    so every band is first UNIONED with the footprints widened by a margin.
    That puts each rim strictly inside every band instead of across it: the
    constraints nest, which is the one thing Triangle insists on.

    Needs shapely to buffer the lines. The import is deliberately LOCAL: the
    repo-hygiene test forbids a top-level geometry import in model code, and a
    mesh producer that runs once per grid is exactly the place for an optional
    dependency rather than a hard one.

    """
    v = cfg.grid.voronoi
    csv = os.path.join(dataset_dir, 'inputSTREAM.csv')
    have_csv = os.path.exists(csv)
    if not have_csv and warn and v.stream_refine:
        warn('grid.voronoi.stream_refine is on but %s is missing, so there is '
             'no corridor.' % csv)
    try:
        from shapely.geometry import LineString, Polygon
        from shapely.ops import unary_union
    except ImportError:
        if warn:
            warn('shapely is not available, so the stream corridor is NOT '
                 'refined and the mesh is uniform at cell_far.')
        return None
    feats = ([LineString(p) for p in stream_lines(csv)]
             if have_csv and v.stream_refine else [])
    size_pond = float(v.cell_pond)
    zones = (_pond_zones(ponds, Polygon) if size_pond > 0.0
             else _pond_footprints([p[0] for p in ponds], Polygon))
    # The footprints are part of the geometry the bands are buffered from,
    # WHATEVER cell_pond says. With cell_pond at 0 that is the whole story --
    # a pond is meshed at cell_near_stream and grades out with the corridor.
    #
    # With cell_pond set, the rims are added as regions further down, and
    # they are added INSIDE the bands. Being part of the buffered geometry is
    # what makes that safe: band k then clears every rim by its own distance,
    # which GROWS with k, so no two bands ever follow the same arc. Padding
    # the bands by a constant instead -- the first thing tried -- made every
    # band trace the same pond bulge, and Triangle stops at "topological
    # inconsistency after splitting a segment" because the arcs coincide.
    feats.extend(zones)
    # Only to keep a BAND's seed out of a pond zone: the band's maximum area
    # would otherwise be the one that applies there.
    keep_out = unary_union(zones) if zones and size_pond > 0.0 else None
    if not feats:
        # Nothing mapped to refine along at all.
        return None
    if keep_out is not None and not (have_csv and v.stream_refine):
        # A pond size with no corridor: the zones ARE the refinement.
        _add_pond_zones(tri, zones, size_pond, warn)
        return keep_out
    lines = unary_union(feats)
    if warn and zones:
        warn('grid.voronoi: %d pond(s) %s.'
             % (len(zones),
                'given a sizing zone at %g x the pond radius, holding %g m '
                'cells' % (POND_ZONE_FACTOR, size_pond) if size_pond > 0.0
                else 'refined with the stream corridor, so a pond carries '
                     'cell_near_stream and grades out with it'))
    # Graded bands: the innermost carries cell_near_stream, each successive
    # band relaxes towards cell_far. Jumping straight from one size to the
    # other makes badly shaped cells along the seam.
    #
    # WP1d panel 1, D2: the bands are DERIVED from cell_near_stream, cell_far,
    # stream_buffer and grade_ratio -- see GridVoronoi.bands(). They used to be
    # a hand-written list unioned with stream_buffer, which let the outermost
    # band fall outside the corridor it was supposed to end at.
    bands = [(float(d), float(s)) for d, s in v.graded_bands()]
    if not bands:
        # No corridor is not the same as nothing to do: a pond SIZE stands on
        # its own, and the zones still have to be added or the size asked for
        # on the panel would quietly do nothing.
        if keep_out is not None:
            _add_pond_zones(tri, zones, size_pond, warn)
            return keep_out
        if warn:
            warn('grid.voronoi: no transition bands (cell_near_stream %g, '
                 'cell_far %g, stream_buffer %g) -- the mesh is uniform.'
                 % (v.cell_near_stream, v.cell_far, v.stream_buffer))
        return None
    added = 0
    prev = None
    for k, (dist, size) in enumerate(bands):
        # A buffer carries ~1 m of boundary detail at these radii, and every
        # vertex of it is a point Triangle must honour -- which is how a
        # corridor meant to hold 15 m cells came out at 3.7 m. Coarsen the
        # arc and simplify to a quarter of the cell this band carries: the
        # band moves by less than that, and the mesh is free again.
        poly = lines.buffer(dist, quad_segs=3).simplify(max(size / 4.0, 0.5))
        # The band's BOUNDARY has to go in as segments. A region is the
        # connected part of the triangulation containing its seed point, and
        # the parts are bounded by segments -- so without this every seed
        # point falls in the same region (the whole catchment) and one
        # maximum area silently applies everywhere. That is the bug that made
        # a "refined" mesh come out uniform: the sizes were computed, the
        # regions were seeded, and nothing bounded them.
        for geom in getattr(poly, 'geoms', [poly]):
            ring = list(geom.exterior.coords)[:-1]
            if len(ring) >= 3:
                tri.add_polygon(ring)
                added += 1
        # The seed must land in THIS band and not in the finer one inside it,
        # so it goes in the annulus rather than in the buffer.
        area = poly if prev is None else poly.difference(prev)
        if keep_out is not None:
            # ... and the band's own seed must not land in a pond zone, or
            # the band's maximum area would be the one applied there.
            area = area.difference(keep_out)
        if area.is_empty:
            prev = poly
            continue
        # EVERY disjoint part, not the largest: a stream network in two
        # catwalks makes an annulus in several pieces, and a piece with no
        # seed of its own is a piece with no area constraint at all -- it
        # would be meshed at whatever the quality rules allow.
        for part in getattr(area, 'geoms', [area]):
            if part.area <= 0.0:
                continue
            pt = part.representative_point()
            tri.add_region((pt.x, pt.y), attribute=k + 1,
                           maximum_area=_max_area_for(size))
        prev = poly
    if warn and not added:
        warn('grid.voronoi: the stream corridor produced no polygon, so the '
             'mesh is uniform at cell_far.')
    if keep_out is not None:
        _add_pond_zones(tri, zones, size_pond, warn)
    return prev if keep_out is None else unary_union([prev, keep_out])


_PRODUCERS = {
    'structured': _produce_structured,
    'disv': _produce_structured,
    'quadtree': _produce_quadtree,
    'voronoi': _produce_voronoi,
}


def build_mesh(cfg, cMF, cache_dir=None, dataset_dir=None, model_ws=None,
               warn=None, force=False):
    """Produce DISV gridprops for the configured ``[grid] kind``.

    ``force`` re-meshes even when the cache holds this signature. It is an
    ARGUMENT and not a configuration key (WP1d panel 1, D2 of §2A.7): a
    rebuild is something a modeller does once, from panel 1's *Create grid*,
    not a setting that stays true in the file long after the run that needed
    it -- which is how the CdL grid-design cache served the wrong mesh.

    Returns ``(gridprops, info)`` where ``info`` records what happened --
    the producer, the signature and whether the cache was used -- so the run
    log can state which mesh it actually ran on.
    """
    if warn is None:
        def warn(msg):
            print('WARNING: %s' % msg)
    kind = cfg.grid_kind
    producer = _PRODUCERS.get(kind)
    if producer is None:
        raise MeshBuildError('unknown grid.kind %r (known: %s)'
                             % (kind, ', '.join(sorted(_PRODUCERS))))
    sig = mesh_signature(cfg, cMF)
    cached = False
    gp = None
    if cache_dir and not force:
        gp = _load_cached(cache_dir, kind, sig)
        cached = gp is not None
    if gp is None:
        gp = normalise(producer(cfg, cMF, dataset_dir=dataset_dir,
                                model_ws=model_ws, warn=warn))
        if cache_dir:
            _save_cached(cache_dir, kind, sig, gp)
    gp.setdefault('ncpl', len(gp['cell2d']))
    gp.setdefault('nlay', int(cMF.nlay))
    info = {'kind': kind, 'signature': sig, 'cached': cached,
            'ncpl': int(gp['ncpl'])}
    info.update(cell_size_report(gp))
    return gp, info


def cell_size_report(gridprops):
    """Achieved cell areas, so a requested cell size can be CHECKED.

    ``_max_area_for`` converts a wanted cell size into a Triangle maximum
    triangle area through an equilateral-triangulation identity that is only
    approximate. Reporting what came out keeps that a starting point rather
    than an article of faith.
    """
    from marmites_grid import polygon_area
    vxy = {int(v[0]): (float(v[1]), float(v[2]))
           for v in gridprops['vertices']}
    areas = np.array([polygon_area([vxy[int(iv)] for iv in r[4:4 + int(r[3])]])
                      for r in gridprops['cell2d']], dtype=float)
    if areas.size == 0:
        return {}
    return {'area_mean': float(areas.mean()),
            'area_min': float(areas.min()),
            'area_max': float(areas.max()),
            'size_equiv': float(np.sqrt(areas.mean()))}
