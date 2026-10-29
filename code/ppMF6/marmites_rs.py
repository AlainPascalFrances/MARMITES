# -*- coding: utf-8 -*-
"""Satellite NDVI for a catchment (user, 2026-10-29): Landsat 5, 7, 8 and 9
at 30 m and Sentinel-2 at 10 m, to compare with the model's wetness W
(marmites_wetimage).

The scenes come from Microsoft's Planetary Computer: its STAC catalogue
lists them, and their Cloud-Optimised GeoTIFFs are read over the
catchment's window only, with no account (a short-lived token the
catalogue hands out signs each file). Landsat is the USGS Collection 2
Level-2 surface reflectance; Sentinel-2 the ESA Level-2A.

What a fetch writes, in ``<RS>/<landsat7_8|sentinel2>/ndvi`` (the folders
the user named), on the model's CRS at the sensor's resolution:

    <platform>_<YYYY-MM-DD>_ndvi.tif   NDVI, float32, NaN where masked
    <platform>_<YYYY-MM-DD>_rgb.tif    the scene's true colour, uint8 --
                                       the visual support of the panel
    <platform>_<YYYY-MM-DD>.json       the scene and how it was read

Masked: clouds, their shadows, cirrus, snow and missing data (Landsat 7's
stripe gaps) -- the Landsat QA_PIXEL bits, the Sentinel-2 scene classes.
"""

import glob
import json
import os
from datetime import date

import numpy as np

__all__ = ['SENSORS', 'search', 'fetch', 'products', 'reflectance',
           'landsat_bad', 's2_bad']

STAC = 'https://planetarycomputer.microsoft.com/api/stac/v1/search'
SAS = 'https://planetarycomputer.microsoft.com/api/sas/v1/token/%s'

SENSORS = {
    'landsat': {
        'label': 'Landsat 5, 7, 8, 9 (30 m)', 'collection': 'landsat-c2-l2',
        'platforms': ('landsat-5', 'landsat-7', 'landsat-8', 'landsat-9'),
        'res': 30.0, 'folder': 'landsat7_8',
        'bands': {'red': 'red', 'nir': 'nir08', 'green': 'green',
                  'blue': 'blue', 'mask': 'qa_pixel'}},
    'sentinel2': {
        'label': 'Sentinel-2 (10 m)', 'collection': 'sentinel-2-l2a',
        'platforms': ('sentinel-2a', 'sentinel-2b', 'sentinel-2c'),
        'res': 10.0, 'folder': 'sentinel2',
        'bands': {'red': 'B04', 'nir': 'B08', 'green': 'B03', 'blue': 'B02',
                  'mask': 'SCL'}},
}


def folder(rs_root, sensor):
    """Where a sensor's products are kept: ``<RS>/<folder>/ndvi``."""
    return os.path.join(str(rs_root), SENSORS[sensor]['folder'], 'ndvi')


# ------------------------------------------------------------------ #
# the catalogue
# ------------------------------------------------------------------ #

def search(sensor, bbox_ll, start, end, max_cloud=None, session=None):
    """The scenes of ``sensor`` over ``bbox_ll`` (lon/lat: W, S, E, N)
    between ``start`` and ``end`` (dates), the oldest first: dicts of
    ``id, date, platform, cloud`` (the SCENE's cloud cover, %), its assets
    and its properties. ``max_cloud`` drops the cloudier scenes. Only the
    scenes whose footprint holds the whole box, one per platform and day
    (the least cloudy): Landsat's overlapping rows and Sentinel-2's tiles
    give several."""
    import requests
    from shapely.geometry import box, shape
    s = SENSORS[sensor]
    want = box(*bbox_ll)
    body = {'collections': [s['collection']], 'bbox': list(bbox_ll),
            'datetime': '%sT00:00:00Z/%sT23:59:59Z' % (start, end),
            'limit': 500}
    if max_cloud is not None:
        body['query'] = {'eo:cloud_cover': {'lte': float(max_cloud)}}
    post = (session or requests).post
    out, seen = [], set()
    while True:
        r = post(STAC, json=body, timeout=60)
        r.raise_for_status()
        js = r.json()
        for f in js.get('features', []):
            p = f.get('properties', {})
            plat = str(p.get('platform', '')).lower()
            if plat not in s['platforms'] or f['id'] in seen:
                continue
            seen.add(f['id'])
            geom = f.get('geometry')
            if geom and not shape(geom).contains(want):
                continue
            out.append({'id': f['id'], 'sensor': sensor,
                        'collection': s['collection'], 'platform': plat,
                        'date': str(p.get('datetime', ''))[:10],
                        'cloud': p.get('eo:cloud_cover'),
                        'assets': {k: v.get('href') for k, v in
                                   f.get('assets', {}).items()},
                        'properties': {k: p.get(k) for k in (
                            'datetime', 's2:processing_baseline',
                            'landsat:wrs_path', 'landsat:wrs_row',
                            's2:mgrs_tile', 'proj:epsg')}})
        nxt = [lk for lk in js.get('links', []) if lk.get('rel') == 'next']
        if not nxt:
            break
        body = nxt[0].get('body', body)
    best = {}
    for d in out:
        k = (d['platform'], d['date'])
        c = d['cloud'] if d['cloud'] is not None else 101.0
        if k not in best or c < (best[k]['cloud'] if best[k]['cloud']
                                 is not None else 101.0):
            best[k] = d
    return sorted(best.values(), key=lambda d: (d['date'], d['platform']))


_TOKENS = {}


def _signed(href, collection):
    """A catalogue file's URL with the read token appended (a local path,
    as in the tests, as it is)."""
    if not str(href).startswith('http'):
        return href
    if collection not in _TOKENS:
        import requests
        r = requests.get(SAS % collection, timeout=60)
        r.raise_for_status()
        _TOKENS[collection] = r.json()['token']
    return '%s%s%s' % (href, '&' if '?' in href else '?', _TOKENS[collection])


# ------------------------------------------------------------------ #
# reflectance and masks
# ------------------------------------------------------------------ #

# QA_PIXEL bits: 0 fill, 1 dilated cloud, 2 cirrus (8 and 9), 3 cloud,
# 4 cloud shadow, 5 snow
_LANDSAT_BAD = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5)
# Sentinel-2 scene classes: 0 no data, 1 saturated or defective, 3 cloud
# shadow, 8 and 9 cloud, 10 thin cirrus, 11 snow
_S2_BAD = (0, 1, 3, 8, 9, 10, 11)


def landsat_bad(qa):
    return (np.asarray(qa).astype(np.uint32) & _LANDSAT_BAD) > 0


def s2_bad(scl):
    return np.isin(np.asarray(scl), _S2_BAD)


def reflectance(sensor, dn, props=None):
    """Surface reflectance of the digital numbers; NaN where there is no
    data (0). Landsat C2 L2: DN x 2.75e-5 - 0.2. Sentinel-2 L2A: DN /
    10000, less 0.1 from processing baseline 04.00 (25 January 2022) on,
    when ESA added an offset of 1000 to every band."""
    dn = np.asarray(dn, dtype=float)
    if sensor == 'landsat':
        r = dn * 2.75e-5 - 0.2
    else:
        pb = str((props or {}).get('s2:processing_baseline') or '0')
        try:
            off = 0.1 if float(pb) >= 4.0 else 0.0
        except ValueError:
            off = 0.0
        r = dn / 10000.0 - off
    return np.where(dn > 0, r, np.nan)


def ndvi(red, nir):
    s = nir + red
    with np.errstate(invalid='ignore', divide='ignore'):
        v = (nir - red) / s
    return np.where(np.isfinite(v) & (s > 0), np.clip(v, -1, 1), np.nan)


# ------------------------------------------------------------------ #
# the fetch
# ------------------------------------------------------------------ #

def _read(href, crs, transform, shape):
    """One band of a (remote) GeoTIFF on the target grid, nearest
    neighbour (no mixing of a stripe gap's zeros into its neighbours)."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.vrt import WarpedVRT
    env = dict(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',
               CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif,.TIF,.tiff',
               GDAL_HTTP_MAX_RETRY='3', GDAL_HTTP_RETRY_DELAY='2')
    with rasterio.Env(**env):
        with rasterio.open(href) as src:
            with WarpedVRT(src, crs=crs, transform=transform,
                           width=shape[1], height=shape[0],
                           resampling=Resampling.nearest,
                           src_nodata=src.nodata if src.nodata is not None
                           else 0, nodata=0) as vrt:
                return vrt.read(1)


def _stretch(rgb, valid):
    """True colour to 0..255: ONE stretch for the three bands, from their
    2nd to their 98th percentile over the valid pixels -- a stretch per
    band upset the colour balance (Landsat 8 over La Mata came out
    purple)."""
    out = np.zeros(rgb.shape, np.uint8)
    m = valid[None] & np.isfinite(rgb)
    if not m.any():
        return out
    lo, hi = np.percentile(rgb[m], [2, 98])
    x = (rgb - lo) / max(hi - lo, 1e-6)
    return np.where(np.isfinite(x), np.clip(x * 255, 0, 255), 0).astype(np.uint8)


def fetch(scene, crs, bounds, dest):
    """Read ``scene`` (from ``search``) over ``bounds`` (xmin, ymin, xmax,
    ymax in ``crs``, the model's), write its NDVI, its true colour and its
    record into ``dest``; returns the record (a dict, with the file paths).
    The grid: the sensor's resolution, aligned on its multiples."""
    import rasterio
    from rasterio.transform import from_origin
    sensor = scene['sensor']
    s = SENSORS[sensor]
    res = s['res']
    x0 = np.floor(bounds[0] / res) * res
    y1 = np.ceil(bounds[3] / res) * res
    nx = int(np.ceil((bounds[2] - x0) / res))
    ny = int(np.ceil((y1 - bounds[1]) / res))
    T = from_origin(x0, y1, res, res)
    a = scene['assets']
    coll = scene['collection']
    band = {}
    for k, key in s['bands'].items():
        if key not in a:
            raise ValueError('scene %s has no %s band (%s)'
                             % (scene['id'], k, key))
        band[k] = _read(_signed(a[key], coll), crs, T, (ny, nx))
    bad = (landsat_bad(band['mask']) if sensor == 'landsat'
           else s2_bad(band['mask']))
    props = scene.get('properties', {})
    refl = {k: reflectance(sensor, band[k], props)
            for k in ('red', 'nir', 'green', 'blue')}
    v = ndvi(refl['red'], refl['nir'])
    v = np.where(bad, np.nan, v).astype(np.float32)
    valid = np.isfinite(v)
    rgb = _stretch(np.stack([refl['red'], refl['green'], refl['blue']]),
                   ~bad)
    os.makedirs(dest, exist_ok=True)
    stem = '%s_%s' % (scene['platform'], scene['date'])
    prof = dict(driver='GTiff', width=nx, height=ny, crs=crs, transform=T,
                compress='deflate')
    p_ndvi = os.path.join(dest, stem + '_ndvi.tif')
    with rasterio.open(p_ndvi, 'w', count=1, dtype='float32', nodata=np.nan,
                       **prof) as d:
        d.write(v, 1)
    p_rgb = os.path.join(dest, stem + '_rgb.tif')
    with rasterio.open(p_rgb, 'w', count=3, dtype='uint8', **prof) as d:
        d.write(rgb)
    rec = {'id': scene['id'], 'sensor': sensor, 'platform': scene['platform'],
           'date': scene['date'], 'collection': coll,
           'scene_cloud': scene.get('cloud'), 'res': res, 'crs': str(crs),
           'clear_share_window': float(valid.mean()),
           'properties': props, 'ndvi': os.path.basename(p_ndvi),
           'rgb': os.path.basename(p_rgb),
           'mask': ('QA_PIXEL fill, dilated cloud, cirrus, cloud, shadow, '
                    'snow' if sensor == 'landsat' else
                    'SCL 0, 1, 3, 8, 9, 10, 11')}
    with open(os.path.join(dest, stem + '.json'), 'w', encoding='utf-8') as fh:
        json.dump(rec, fh, indent=1)
    rec['dir'] = dest
    return rec


def products(*dirs):
    """The products fetched into ``dirs``: their records, the oldest
    first, each with its folder (``dir``)."""
    out = []
    for d in dirs:
        for p in glob.glob(os.path.join(str(d), '*.json')):
            try:
                with open(p, encoding='utf-8') as fh:
                    rec = json.load(fh)
            except (OSError, ValueError):
                continue
            if 'ndvi' in rec and os.path.exists(
                    os.path.join(os.path.dirname(p), rec['ndvi'])):
                rec['dir'] = os.path.dirname(p)
                out.append(rec)
    return sorted(out, key=lambda r: (r.get('date', ''), r.get('platform', '')))


def catchment_bounds(g, margin=100.0):
    """The active cells' extent (model CRS) with a margin, from a
    marmites_daymap grid."""
    act = ~g['off'].all(axis=0)
    xy = np.concatenate([g['polys'][c][:, :2] for c in np.flatnonzero(act)])
    return (float(xy[:, 0].min() - margin), float(xy[:, 1].min() - margin),
            float(xy[:, 0].max() + margin), float(xy[:, 1].max() + margin))


def lonlat_bbox(bounds, crs):
    """``bounds`` (model CRS) as a lon/lat box for the catalogue."""
    from pyproj import Transformer
    tr = Transformer.from_crs(crs, 'EPSG:4326', always_xy=True)
    xs = (bounds[0], bounds[2])
    ys = (bounds[1], bounds[3])
    pts = [tr.transform(x, y) for x in xs for y in ys]
    return (min(p[0] for p in pts), min(p[1] for p in pts),
            max(p[0] for p in pts), max(p[1] for p in pts))


def as_date(s):
    return date.fromisoformat(str(s)[:10])
