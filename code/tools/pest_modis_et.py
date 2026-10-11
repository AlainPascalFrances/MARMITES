# -*- coding: utf-8 -*-
"""MODIS MOD16A2GF actual ET over the catchment, per 8-day composite (WP7).

The ET observation group of the PEST++ calibration, as on CdL: the user
downloads the MOD16A2GF tiles (500 m, 8-day composites, ET in kg/m2 = mm
per composite) into a folder of the remote-sensing root (Home, ``rs``) --
``MODIS4LaMata`` -- and this turns them into ONE series: the catchment's
area-weighted mean ET of each composite.

    python code/tools/pest_modis_et.py --product MODIS4LaMata \\
        --config code/configs/lamata.toml --out modis_et_catchment.csv

- Each pixel weighs by the area it shares with the catchment polygon (the
  dataset's inputBOUNDARY.geojson, in the model CRS), computed in the MODIS
  sinusoidal projection the pixels are squares in.
- Fill values (32761..32767: water, urban, barren, no data) are left out,
  and the mean is over the valid area; ``valid_frac`` says how much of the
  catchment that was.
- A composite covers 8 days from its date, the last of a year up to
  31 December (5 or 6 days). A composite downloaded twice is read once.
- HDF4-EOS is read through rasterio's GDAL (its HDF4 driver), as
  ``HDF4_EOS:EOS_GRID:"file":MOD_Grid_MOD16A2:ET_500m``.

Geospatial libraries are used here, in a tool, never on the model path.
"""

import argparse
import datetime as dt
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
if CODE not in sys.path:
    sys.path.insert(0, CODE)

SUBDATASET = 'MOD_Grid_MOD16A2:ET_500m'
SCALE = 0.1                 # kg/m2/8day per count
VALID_MAX = 32700           # 32761..32767 are fills
_NAME = re.compile(r'MOD16A2GF?\.A(\d{4})(\d{3})\.(h\d+v\d+)\.')


def composite_period(year, doy):
    """(first day, number of days) of the composite dated year/doy: 8 days,
    the year's last one up to 31 December."""
    start = dt.date(year, 1, 1) + dt.timedelta(days=doy - 1)
    end_year = dt.date(year, 12, 31)
    n = min(8, (end_year - start).days + 1)
    return start, n


def list_composites(product_dir):
    """``[(year, doy, tile, path)]``, one per composite and tile, sorted; a
    composite found twice (two downloads) is kept once."""
    seen = {}
    for p in sorted(glob.glob(os.path.join(product_dir, '**', '*.hdf'),
                              recursive=True)):
        m = _NAME.search(os.path.basename(p))
        if not m:
            continue
        key = (int(m.group(1)), int(m.group(2)), m.group(3))
        seen.setdefault(key, p)
    return [k + (v,) for k, v in sorted(seen.items())]


def pixel_weights(transform, shape, polygon):
    """For a raster of ``shape`` with affine ``transform``, the pixels the
    ``polygon`` (same CRS) overlaps and the area each shares with it:
    ``(rows, cols, area)``, area in the CRS units squared."""
    import shapely
    from shapely.geometry import box
    xmin, ymin, xmax, ymax = polygon.bounds
    a, b, c, d, e, f = (transform.a, transform.b, transform.c, transform.d,
                        transform.e, transform.f)
    if b != 0 or d != 0:
        raise ValueError('a rotated raster is not supported')
    c0 = max(int(np.floor((xmin - c) / a)), 0)
    c1 = min(int(np.ceil((xmax - c) / a)), shape[1] - 1)
    r0 = max(int(np.floor((ymax - f) / e)), 0)
    r1 = min(int(np.ceil((ymin - f) / e)), shape[0] - 1)
    rows, cols, area = [], [], []
    shapely.prepare(polygon)
    for r in range(r0, r1 + 1):
        for k in range(c0, c1 + 1):
            px = box(c + k * a, f + (r + 1) * e, c + (k + 1) * a, f + r * e)
            if not polygon.intersects(px):
                continue
            s = polygon.intersection(px).area
            if s > 0.0:
                rows.append(r)
                cols.append(k)
                area.append(s)
    return np.array(rows, int), np.array(cols, int), np.array(area, float)


def weighted_mean(counts, area):
    """Area-weighted mean ET [mm per composite] over the VALID pixels, and
    the share of the area that was valid."""
    counts = np.asarray(counts, dtype=float)
    ok = (counts >= 0) & (counts <= VALID_MAX)
    tot = area.sum()
    if not ok.any() or tot <= 0:
        return np.nan, 0.0
    return (float(np.sum(counts[ok] * SCALE * area[ok]) / area[ok].sum()),
            float(area[ok].sum() / tot))


def catchment_polygon(boundary_geojson):
    """The catchment as one shapely polygon, from the dataset's boundary."""
    import json
    import shapely
    with open(boundary_geojson, encoding='utf-8') as fh:
        g = json.load(fh)
    geoms = shapely.from_geojson([json.dumps(f['geometry'])
                                  for f in g['features'] if f.get('geometry')])
    return shapely.make_valid(shapely.union_all(geoms))


def catchment_series(product_dir, polygon, epsg, start=None, end=None,
                     verbose=True):
    """The catchment-mean ET of every composite in ``product_dir`` that
    overlaps [start, end] (dates), as a DataFrame."""
    import pandas as pd
    import pyproj
    import rasterio
    from shapely.ops import transform as sh_transform
    comps = list_composites(product_dir)
    if not comps:
        raise SystemExit('no MOD16A2(GF) .hdf file under %s' % product_dir)
    tiles = sorted({c[2] for c in comps})
    if len(tiles) > 1:
        raise SystemExit('composites of %d tiles (%s): one tile is read; '
                         'the catchment must lie in one'
                         % (len(tiles), ', '.join(tiles)))
    rows = []
    cache = None
    for year, doy, tile, path in comps:
        first, n = composite_period(year, doy)
        last = first + dt.timedelta(days=n - 1)
        if start is not None and last < start:
            continue
        if end is not None and first > end:
            continue
        sds = 'HDF4_EOS:EOS_GRID:"%s":%s' % (path, SUBDATASET)
        with rasterio.open(sds) as src:
            if cache is None:
                tr = pyproj.Transformer.from_crs(
                    pyproj.CRS.from_epsg(int(epsg)),
                    pyproj.CRS.from_wkt(src.crs.to_wkt()), always_xy=True)
                poly = sh_transform(tr.transform, polygon)
                r, c, a = pixel_weights(src.transform, src.shape, poly)
                if not r.size:
                    raise SystemExit('the catchment does not overlap tile %s'
                                     % tile)
                cache = (r, c, a, (src.transform, src.shape))
                if verbose:
                    print('MODIS %s: the catchment overlaps %d pixel(s), '
                          '%.3f km2' % (tile, r.size, a.sum() / 1e6))
            elif (src.transform, src.shape) != cache[3]:
                raise SystemExit('%s is not on the grid of the first '
                                 'composite' % path)
            r, c, a, _ = cache
            win = rasterio.windows.Window(int(c.min()), int(r.min()),
                                          int(c.max() - c.min() + 1),
                                          int(r.max() - r.min() + 1))
            arr = src.read(1, window=win, masked=False)
            counts = arr[r - r.min(), c - c.min()]
        et, frac = weighted_mean(counts, a)
        rows.append({'start': first.isoformat(), 'end': last.isoformat(),
                     'ndays': n, 'et_mm': et,
                     'et_mmd': et / n if np.isfinite(et) else np.nan,
                     'valid_frac': frac, 'year': year, 'doy': doy,
                     'file': os.path.basename(path)})
    df = pd.DataFrame(rows)
    if verbose and len(df):
        print('MODIS ET: %d composite(s), %s .. %s; catchment mean %.2f '
              'mm/d, %.0f mm/yr; %d with no valid pixel'
              % (len(df), df['start'].iloc[0], df['end'].iloc[-1],
                 np.nanmean(df['et_mmd']),
                 np.nansum(df['et_mm']) / (df['ndays'].sum() / 365.25),
                 int(df['et_mm'].isna().sum())))
    return df


def main(argv=None):
    import marmites_config as mcfg
    import mm_paths
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--product', required=True,
                    help='the product folder under the remote-sensing root '
                         '(Home, rs), e.g. MODIS4LaMata, or a full path')
    ap.add_argument('--config', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    cfg = mcfg.load_run_config(a.config)
    prod = (a.product if os.path.isabs(a.product)
            else os.path.join(str(mm_paths.RS), a.product))
    ds = mm_paths.dataset_dir(cfg.paths.case)
    poly = catchment_polygon(os.path.join(str(ds), 'inputBOUNDARY.geojson'))
    dates = [l.split(',')[0].strip()[:10]
             for l in open(os.path.join(str(ds), 'inputDATE.txt'))
             if l.strip() and not l.startswith('#')]
    start = dt.date.fromisoformat(dates[0])
    end = dt.date.fromisoformat(dates[-1])
    df = catchment_series(prod, poly, cfg.grid.crs_epsg, start, end)
    df.to_csv(a.out, index=False, float_format='%.6g')
    print('written', a.out)


if __name__ == '__main__':
    main()
