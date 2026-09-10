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
than warns. ``[grid] rebuild = true`` forces it.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import hashlib
import json
import os

import numpy as np

__all__ = ['build_mesh', 'MeshBuildError', 'stream_lines', 'mesh_signature',
           'watershed_ring', 'normalise', 'cell_size_report',
           'PRODUCER_VERSION']

# Identifies the MESH-PRODUCING BEHAVIOUR, not the module. Bump it on any
# change that would give a different mesh for the same configuration; it is
# part of the cache signature.
PRODUCER_VERSION = 2


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
    """GRIDGEN quadtree, refined along the stream network.

    Refining on the streams rather than uniformly is what makes this a real
    end-to-end proof for WP1c.1: the mesh genuinely has no (row, column), so
    every projection and index path is exercised, not bypassed.
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
    elif csv:
        feats.append((stream_lines(csv), 'line', 2))
    ws = model_ws or os.path.join(os.getcwd(), '_gridgen')
    return build_quadtree(cMF, exe, ws, refine_features=feats,
                          layers=list(range(int(cMF.nlay))))


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

    tri = Triangle(model_ws=ws, exe_name=exe,
                   maximum_area=_max_area_for(v.cell_far))
    tri.add_polygon(ring)

    if v.stream_refine:
        _add_stream_regions(tri, cfg, dataset_dir, warn)
    if v.seed_ponds and warn:
        # A Voronoi cell can never follow a pond outline -- its faces bisect
        # the segment between two generators, so a constraint polygon is cut
        # rather than honoured. CdL's answer was to SEED a generator at the
        # pond centroid instead, which is a LAK concern; it lands with WP4.
        warn('grid.voronoi.seed_ponds is set, but pond seeding arrives with '
             'WP4 (LAK). The mesh is built without it.')

    tri.build(verbose=False)
    gp = VoronoiGrid(tri).get_gridprops_vertexgrid()
    gp['nlay'] = int(cMF.nlay)
    return gp


def _add_stream_regions(tri, cfg, dataset_dir, warn):
    """Refine a corridor along the stream network, graded outward.

    Needs shapely to buffer the lines. The import is deliberately LOCAL: the
    repo-hygiene test forbids a top-level geometry import in model code, and a
    mesh producer that runs once per grid is exactly the place for an optional
    dependency rather than a hard one.
    """
    v = cfg.grid.voronoi
    csv = os.path.join(dataset_dir, 'inputSTREAM.csv')
    if not os.path.exists(csv):
        if warn:
            warn('grid.voronoi.stream_refine is on but %s is missing: the mesh '
                 'will be uniform.' % csv)
        return
    try:
        from shapely.geometry import LineString
        from shapely.ops import unary_union
    except ImportError:
        if warn:
            warn('shapely is not available, so the stream corridor is NOT '
                 'refined and the mesh is uniform at cell_far.')
        return
    lines = unary_union([LineString(p) for p in stream_lines(csv)])
    # Graded bands: the innermost carries cell_near_stream, each successive
    # band relaxes towards cell_far. Jumping straight from one size to the
    # other makes badly shaped cells along the seam.
    bands = sorted(set([float(x) for x in v.trans_levels]
                       + [float(v.stream_buffer)]))
    n = len(bands)
    for k, dist in enumerate(bands):
        frac = k / float(max(n - 1, 1))
        size = v.cell_near_stream + frac * (v.cell_far - v.cell_near_stream)
        poly = lines.buffer(dist)
        pt = poly.representative_point()
        tri.add_region((pt.x, pt.y), attribute=k + 1,
                       maximum_area=_max_area_for(size))


_PRODUCERS = {
    'structured': _produce_structured,
    'disv': _produce_structured,
    'quadtree': _produce_quadtree,
    'voronoi': _produce_voronoi,
}


def build_mesh(cfg, cMF, cache_dir=None, dataset_dir=None, model_ws=None,
               warn=None):
    """Produce DISV gridprops for the configured ``[grid] kind``.

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
    if cache_dir and not cfg.grid.rebuild:
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
