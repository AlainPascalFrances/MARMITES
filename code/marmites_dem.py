# -*- coding: utf-8 -*-
"""The land surface, at its own resolution, wrapped onto the model grid.

WHY THIS EXISTS
---------------
The elevation used to reach the model as a raster already resampled to the
model's 50 m rectangle, which was then projected again onto whatever mesh the
run used. Two resamplings, and the first one threw away everything the 5 m
survey knew: on La Mata the round trip costs 0.65 m rms and 4.4 m at worst
against wrapping the original directly.

So the DEM stays at NATIVE RESOLUTION in the dataset -- grid-independent,
exactly like the vector layers -- and is wrapped onto the cells of whichever
grid panel 1 produced, structured, disv, voronoi or quadtree alike. The
wrapping is AREA-WEIGHTED: a cell's elevation is the mean of the DEM over the
cell's own polygon, not the pixel under its centre.

NO rasterio HERE. The converter reads the GIS raster, in whatever format GDAL
opens, and writes an ESRI ASCII grid into the dataset; this module reads that
with numpy alone, because the model path must keep running where rasterio is
not installed (cookbook decision D3).
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.1.0.dev0"

import hashlib
import json
import os

import numpy as np

__all__ = ['DEMError', 'read_asc', 'read_asc_header', 'write_asc',
           'wrap_to_grid', 'dem_path', 'HEADER_KEYS']

HEADER_KEYS = ('ncols', 'nrows', 'xllcorner', 'yllcorner', 'cellsize',
               'nodata_value')

# What the converter writes into the dataset, at the DEM's own resolution.
DATASET_DEM = 'inputDEM.asc'


class DEMError(Exception):
    """The elevation raster cannot serve the use asked of it."""


def dem_path(dataset_dir, name=DATASET_DEM):
    return os.path.join(str(dataset_dir), name)


def _header_from(fh, path):
    """The six header lines, lower-cased, off an already-open handle."""
    head = {}
    for _ in range(6):
        parts = fh.readline().split()
        if len(parts) != 2:
            raise DEMError('%s: expected six header lines' % path)
        head[parts[0].strip().lower()] = float(parts[1])
    missing = [k for k in HEADER_KEYS if k not in head]
    if missing:
        raise DEMError('%s: header has no %s' % (path, ', '.join(missing)))
    return head


def read_asc_header(path):
    """Where an ESRI ASCII grid SITS, without reading its body.

    The rectangle check on panel 1 asks this of every raster in a dataset;
    loading the arrays to compare six numbers would make the page crawl.
    """
    if not os.path.exists(path):
        raise DEMError('raster not found: %s' % path)
    with open(path, encoding='utf-8', errors='replace') as fh:
        return _header_from(fh, path)


def read_asc(path):
    """``(masked array, header)`` of an ESRI ASCII grid.

    The array is masked at NODATA, so a hole in the survey stays a hole
    instead of becoming an elevation of -9999 that averages into its
    neighbours.
    """
    if not os.path.exists(path):
        raise DEMError('elevation raster not found: %s' % path)
    with open(path, encoding='utf-8', errors='replace') as fh:
        head = _header_from(fh, path)
        try:
            arr = np.loadtxt(fh, dtype=float)
        except ValueError as exc:
            raise DEMError('%s: %s' % (path, exc))
    nrows, ncols = int(head['nrows']), int(head['ncols'])
    arr = np.atleast_2d(arr)
    if arr.shape != (nrows, ncols):
        raise DEMError('%s: header says %d x %d and the body is %d x %d'
                       % (path, nrows, ncols, arr.shape[0], arr.shape[1]))
    return np.ma.masked_values(arr, head['nodata_value'], atol=1e-6), head


def write_asc(path, arr, xll, yll, cellsize, nodata=-9999.0, fmt='%.2f'):
    """Write an ESRI ASCII grid.

    ``%.10g`` on the header numbers, not ``%g``: a UTM northing is seven
    digits and ``%g`` writes it as 4.55305e+06, which readers of this format
    do not parse -- and the file looks right until something fails on it.
    """
    arr = np.ma.filled(np.ma.asarray(arr, dtype=float), nodata)
    nrows, ncols = arr.shape
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('ncols         %d\n' % ncols)
        fh.write('nrows         %d\n' % nrows)
        fh.write('xllcorner     %.10g\n' % xll)
        fh.write('yllcorner     %.10g\n' % yll)
        fh.write('cellsize      %.10g\n' % cellsize)
        fh.write('NODATA_value  %.10g\n' % nodata)
        for r in range(nrows):
            fh.write(' '.join(fmt % v for v in arr[r]) + '\n')
    return path


def signature(head, gridprops):
    """Identifies a wrapped result: this DEM, on this mesh."""
    payload = {'dem': {k: float(head[k]) for k in HEADER_KEYS},
               'ncpl': int(gridprops.get('ncpl', len(gridprops['cell2d']))),
               'nvert': len(gridprops['vertices'])}
    blob = json.dumps(payload, sort_keys=True).encode('utf-8')
    return hashlib.sha256(blob).hexdigest()[:16]


def wrap_to_grid(path, gridprops, cache_dir=None, force=False, warn=None):
    """Area-weighted DEM elevation per model cell. ``(array (ncpl, 1), info)``.

    Works for every grid kind, because every kind arrives here as DISV
    gridprops -- a structured grid is its own polygons. The heavy lifting is
    ``marmites_mesh.MeshProjection``, which already answers "which source
    cells does this model cell overlap, and by how much"; the only new thing
    is that the source is the DEM's own fine lattice rather than the model's
    rectangle.

    Cached on (DEM header, mesh) because it is seconds, not milliseconds: 18 s
    for the 5 m La Mata DEM on a 5778-cell Voronoi mesh.
    """
    import marmites_mesh

    if warn is None:
        def warn(msg):
            print('WARNING: %s' % msg)

    dem, head = read_asc(path)
    sig = signature(head, gridprops)
    cache = (os.path.join(cache_dir, 'dem_%s.npy' % sig) if cache_dir else '')
    ncpl = int(gridprops.get('ncpl', len(gridprops['cell2d'])))
    if cache and not force and os.path.exists(cache):
        out = np.load(cache)
        if out.shape == (ncpl, 1):
            return np.ma.masked_values(out, -9999.0, atol=1e-6), {
                'signature': sig, 'cached': True, 'ncpl': ncpl,
                'cellsize': float(head['cellsize'])}

    cs = float(head['cellsize'])
    nrows, ncols = int(head['nrows']), int(head['ncols'])
    proj = marmites_mesh.MeshProjection(
        gridprops, np.full(ncols, cs), np.full(nrows, cs),
        float(head['xllcorner']), float(head['yllcorner']))
    out = proj.sample2d(dem, fill=-9999.0, how='area')
    arr = np.ma.masked_values(np.asarray(np.ma.filled(out, -9999.0)),
                              -9999.0, atol=1e-6).reshape(ncpl, 1)
    missing = int(np.ma.getmaskarray(arr).sum())
    if missing:
        warn('%d of %d cell(s) have no DEM under them; their elevation is '
             'NODATA. The raster does not cover the whole grid.'
             % (missing, ncpl))
    if cache:
        os.makedirs(cache_dir, exist_ok=True)
        np.save(cache, np.ma.filled(arr, -9999.0))
    return arr, {'signature': sig, 'cached': False, 'ncpl': ncpl,
                 'cellsize': cs, 'missing': missing}
