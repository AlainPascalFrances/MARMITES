# -*- coding: utf-8 -*-
"""Satellite soil moisture against the model's saturation Se (user,
2026-10-29, branch swi-soil-moisture): a comparison now, a calibration
target later -- as MODIS ET is (group et_rs).

The products, one sub-folder each of the SOIL MOISTURE folder of panel 0
(mm_paths.SOIL_MOISTURE, by default RS/soil_moisture), downloaded by the
user:

    clms_swi_0.1deg  Copernicus Soil Water Index, global, daily, ~10 km,
                     ASCAT, 2007 on; root zone (T-values), % saturation
    clms_swi_1km     the same over Europe at 1 km (Sentinel-1 + ASCAT),
                     2015 on: the calibration target for Casa de Lobos
    clms_ssm_1km     Copernicus Surface Soil Moisture, Europe, 1 km,
                     Sentinel-1, every 1.5-4 days; top ~5 cm, % saturation

(SMOS-BEC L4 1 km, tried on 2026-10-30, is left out -- user: "skip
SMOS-BEC for La Mata and also for the suite": a 9 MB European map per pass
and day, no subset service reachable here, and La Mata often outside the
pass.)

A folder holds either the product's files (NetCDF or GeoTIFF) or the
Copernicus Browser's Statistical Info exports of it (user, 2026-10-30): one
CSV per layer and period, named by the Browser
``<dataset>-<layer>-<from>-<to>.csv`` -- SWI001..SWI100 for the T-values --
each the layer's daily mean over the area drawn there (the catchment), with
its spread. Exports of one layer over several periods join; on a day two of
them share they must agree (else they are of different areas).

The model side is the effective saturation of the soil column, Se = (S -
S_r) / (S_sat - S_r) (marmites_daymap.saturation), averaged over the
catchment -- each cell by its area -- and over each product pixel, by the
area the pixel shares with each cell. The scores and the figure cover the
run's whole hydrological years, as the ET scores do.

A Soil Water Index IS a degree of saturation of the root zone -- the
surface moisture filtered by an exponential of characteristic time T
(Wagner et al. 1999, the recursive form of Albergel et al. 2008) -- so it
is compared with 100 Se as it is. A surface product is filtered the same
way first (T, days), then scaled between its 2nd and 98th percentiles, as
is the model's series for it: its top few centimetres are not the soil
column, so only the timing and the anomalies are compared, not the level.

Writes ``<out>/_output/sm_products/smrs_<product>_T<T>.png`` (and .csv,
.json), listed in the Results calibration tab; Streamlit-free.
"""

import glob
import json
import os
import re
from datetime import date, timedelta

import numpy as np

import marmites_daymap as DM

__all__ = ['PRODUCTS', 'SUBDIR', 'available', 'read', 'exp_filter',
           'metrics', 'hydro_window', 'catchment_geojson', 'compare']

SUBDIR = 'sm_products'                       # under <out>/_output
FINE = 10.0                                  # m, the overlap of pixels and cells
_INK, _MUTED = '#222222', '#666666'
_MODEL, _SAT = '#2a78d6', '#e34948'

PRODUCTS = {
    'clms_swi_0.1deg': {
        'label': 'Copernicus Soil Water Index 0.1 deg (ASCAT), daily',
        'short': 'SWI 0.1 deg', 'kind': 'swi', 'unit': '% saturation',
        'depth': 'root zone (T-value)', 'crs': 'EPSG:4326',
        'tvalues': (1, 5, 10, 15, 20, 40, 60, 100)},
    'clms_swi_1km': {
        'label': 'Copernicus Soil Water Index 1 km Europe (Sentinel-1 + ASCAT)',
        'short': 'SWI 1 km', 'kind': 'swi', 'unit': '% saturation',
        'depth': 'root zone (T-value)', 'crs': 'EPSG:4326',
        'tvalues': (2, 5, 10, 15, 20, 40, 60, 100)},
    'clms_ssm_1km': {
        'label': 'Copernicus Surface Soil Moisture 1 km Europe (Sentinel-1)',
        'short': 'SSM 1 km', 'kind': 'surface', 'unit': '% saturation',
        'depth': 'top ~5 cm', 'crs': 'EPSG:4326', 'tvalues': ()},
}
T_DEFAULT = 10                               # days
VALID_MAX = 100.0              # % of saturation; a value beyond is a flag


def folder(sm_root, product):
    """A product's folder in the soil-moisture folder ``sm_root``."""
    return os.path.join(str(sm_root), product)


def _files(d):
    out = []
    for pat in ('*.nc', '*.nc4', '*.NC', '*.tif', '*.tiff', '*.TIF', '*.csv',
                '*.CSV'):
        out += glob.glob(os.path.join(str(d), '**', pat), recursive=True)
    return sorted(set(out))


def available(sm_root):
    """The products whose folder holds files: ``{product: n files}``."""
    out = {}
    for p in PRODUCTS:
        n = len(_files(folder(sm_root, p)))
        if n:
            out[p] = n
    return out


# ------------------------------------------------------------------ #
# reading
# ------------------------------------------------------------------ #

_DATE = re.compile(r'(?<!\d)(\d{4})(\d{2})(\d{2})(?!\d)')


def _file_date(path):
    for m in _DATE.finditer(os.path.basename(path)):
        try:
            return date(int(m.group(1)), int(m.group(2)), int(m.group(3)))
        except ValueError:
            continue
    return None


def _var_key(product, name):
    """The series a variable of a file is: a T-value (int) for an SWI, 's'
    for a surface moisture, None for anything else."""
    kind = PRODUCTS[product]['kind']
    n = str(name)
    if kind == 'swi':
        m = re.match(r'^SWI_?0*(\d{1,3})$', n, re.I)
        return int(m.group(1)) if m else None
    return 's' if re.match(r'^ssm$', n, re.I) else None


def _axes(ds):
    """The x/y (lon/lat) coordinate names of a dataset."""
    xs = [c for c in ('lon', 'longitude', 'x', 'X') if c in ds.coords or c in ds]
    ys = [c for c in ('lat', 'latitude', 'y', 'Y') if c in ds.coords or c in ds]
    if not xs or not ys:
        raise ValueError('no x/y coordinates in the file')
    return xs[0], ys[0]


def read(product, d):
    """All the files of ``product`` in ``d``: ``{'x', 'y' (pixel centres,
    the product's CRS), 'times' (dates), 'series' {key: (nt, ny, nx)},
    'crs'}``, NaN where there is no value. A file may hold one day or many;
    a day without a time coordinate takes the YYYYMMDD of its name.

    Copernicus Browser exports instead: ``{'aoi': True, 'x': None, 'y':
    None, 'times', 'series' {key: (nt, 1, 1)} the area's mean, 'spread'
    {key: (nt, 2)} its 10th and 90th percentiles, 'crs': None}``."""
    import xarray as xr
    files = _files(d)
    if not files:
        raise ValueError('no file in %s' % d)
    csvs = [f for f in files if f.lower().endswith('.csv')]
    if csvs:
        if len(csvs) < len(files):
            raise ValueError('%s holds both product files and Copernicus '
                             'Browser exports (CSV): keep one kind' % d)
        return _read_exports(product, csvs)
    frames = {}                               # date -> {key: (ny, nx)}
    grid = None
    crs = PRODUCTS[product]['crs']
    for f in files:
        if f.lower().endswith(('.tif', '.tiff')):
            got = _read_tif(product, f)
        else:
            ds = xr.open_dataset(f, mask_and_scale=True)
            try:
                got = _read_nc(product, ds, f)
            finally:
                ds.close()
        x, y, fcrs, per = got
        if fcrs:
            crs = fcrs
        if grid is None:
            grid = (x, y)
        elif x.shape != grid[0].shape or y.shape != grid[1].shape or \
                not (np.allclose(x, grid[0]) and np.allclose(y, grid[1])):
            raise ValueError('%s is on another grid than the files before it'
                             % os.path.basename(f))
        for t, vals in per.items():
            frames.setdefault(t, {}).update(vals)
    if not frames:
        raise ValueError('no %s variable found in the files of %s'
                         % (PRODUCTS[product]['short'], d))
    times = sorted(frames)
    keys = sorted({k for v in frames.values() for k in v},
                  key=lambda k: (isinstance(k, str), k))
    ny, nx = grid[1].size, grid[0].size
    series = {k: np.full((len(times), ny, nx), np.nan) for k in keys}
    for i, t in enumerate(times):
        for k, a in frames[t].items():
            series[k][i] = a
    return {'x': grid[0], 'y': grid[1], 'times': times, 'series': series,
            'crs': crs}


def _read_nc(product, ds, path):
    xn, yn = _axes(ds)
    x = np.asarray(ds[xn].values, float)
    y = np.asarray(ds[yn].values, float)
    fcrs = None
    for v in ds.variables.values():
        wkt = v.attrs.get('spatial_ref') or v.attrs.get('crs_wkt')
        if wkt:
            fcrs = wkt
            break
    per = {}
    tname = next((c for c in ('time', 't') if c in ds.dims), None)
    fdate = _file_date(path)
    for name, da in ds.data_vars.items():
        key = _var_key(product, name)
        if key is None:
            continue
        a = np.asarray(da.values, float)
        dims = list(da.dims)
        if tname and tname in dims:
            tt = [np.datetime64(t, 'D').astype(object) for t in
                  ds[tname].values]
            a = np.moveaxis(a, dims.index(tname), 0)
            dims = [tname] + [d_ for d_ in dims if d_ != tname]
        else:
            if fdate is None:
                raise ValueError('%s: no time coordinate, no YYYYMMDD in its '
                                 'name' % os.path.basename(path))
            tt = [fdate]
            a = a[None]
            dims = ['_t'] + dims
        # (t, y, x) whatever the order in the file
        sp = [d_ for d_ in dims[1:]]
        if sp == [xn, yn]:
            a = np.swapaxes(a, 1, 2)
        a = np.where((a >= 0) & (a <= VALID_MAX), a, np.nan)
        for i, t in enumerate(tt):
            per.setdefault(t, {})[key] = a[i]
    return x, y, fcrs, per


def _read_tif(product, path):
    import rasterio
    t = _file_date(path)
    if t is None:
        raise ValueError('%s: no YYYYMMDD in its name' % os.path.basename(path))
    with rasterio.open(path) as src:
        T = src.transform
        x = T.c + (np.arange(src.width) + 0.5) * T.a
        y = T.f + (np.arange(src.height) + 0.5) * T.e
        per = {t: {}}
        for b in range(1, src.count + 1):
            nm = (src.descriptions[b - 1] or
                  ('SWI_010' if PRODUCTS[product]['kind'] == 'swi' else 'ssm'))
            key = _var_key(product, nm)
            if key is None:
                continue
            a = src.read(b).astype(float)
            if src.nodata is not None:
                a = np.where(a == src.nodata, np.nan, a)
            sc = (src.scales or [1.0] * src.count)[b - 1]
            off = (src.offsets or [0.0] * src.count)[b - 1]
            a = a * sc + off
            per[t][key] = np.where((a >= 0) & (a <= VALID_MAX), a, np.nan)
        crs = src.crs.to_string() if src.crs else None
    return x, y, crs, per


def _read_export(product, path):
    """One Statistical Info export of the Copernicus Browser: ``(key,
    {date: (mean, p10, p90)})``, the key from the layer in the file's name.
    NaN on a day where more than half the area has no value, or where a
    flag (a value beyond the product's range) is in the statistics."""
    import pandas as pd
    name = os.path.basename(path)
    key = next((k for k in (_var_key(product, t) for t in
                            re.split(r'[-.]', name)) if k is not None), None)
    if key is None:
        raise ValueError('%s: no %s layer in its name (the Browser names an '
                         'export <dataset>-<layer>-<from>-<to>.csv)'
                         % (name, PRODUCTS[product]['short']))
    df = pd.read_csv(path)
    dc = next((c for c in df.columns if c.split('/')[-1] == 'date'), None)
    pre = dc.rsplit('/', 1)[0] + '/' if dc and '/' in dc else ''
    if dc is None or pre + 'mean' not in df.columns:
        raise ValueError('%s: not a Statistical Info export (no date and mean '
                         'columns)' % name)

    def col(s, fill):
        c = pre + s
        return (pd.to_numeric(df[c], errors='coerce').to_numpy(float)
                if c in df.columns else np.full(len(df), fill))
    mean, lo, hi = col('mean', np.nan), col('min', 0.0), col('max', 0.0)
    p10, p90 = col('p10', np.nan), col('p90', np.nan)
    n, nd = col('sampleCount', 1.0), col('noDataCount', 0.0)
    bad = (lo < 0) | (hi > VALID_MAX) | (2 * nd > n)
    days = [date.fromisoformat(str(s)[:10]) for s in df[dc]]
    per = {}
    for i, t in enumerate(days):
        v = ((np.nan,) * 3 if bad[i] else (mean[i], p10[i], p90[i]))
        if t in per and not _same(per[t][0], v[0]):
            raise ValueError('%s: two rows of %s differ' % (name, t))
        per[t] = v
    return key, per


def _same(a, b):
    return (np.isnan(a) and np.isnan(b)) or abs(a - b) <= 1e-6


def _read_exports(product, files):
    """The Browser exports of ``product``, joined layer by layer (see
    read)."""
    got = {}                                  # key -> {date: (m, p10, p90)}
    for f in files:
        key, per = _read_export(product, f)
        have = got.setdefault(key, {})
        for t, v in per.items():
            if t in have and not _same(have[t][0], v[0]):
                raise ValueError(
                    '%s and an export before it differ on %s (%.2f against '
                    '%.2f): exports of different areas?' % (
                        os.path.basename(f), t, v[0], have[t][0]))
            have[t] = v
    times = sorted({t for per in got.values() for t in per})
    series, spread = {}, {}
    for k, per in got.items():
        a = np.array([per.get(t, (np.nan,) * 3) for t in times], float)
        series[k] = a[:, 0].reshape(-1, 1, 1)
        spread[k] = a[:, 1:]
    return {'aoi': True, 'x': None, 'y': None, 'times': times,
            'series': series, 'spread': spread, 'crs': None}


# ------------------------------------------------------------------ #
# the model: Se over the catchment and over each pixel
# ------------------------------------------------------------------ #

def _pixel_polys(x, y, crs_from, crs_to):
    """Each pixel's outline (row-major, y then x) in ``crs_to``, from the
    regular grid of centres ``x``, ``y`` in ``crs_from``, each edge cut in
    8 so a reprojected pixel keeps its shape."""
    from pyproj import Transformer
    tr = Transformer.from_crs(crs_from, crs_to, always_xy=True)
    dx = np.median(np.diff(x)) if x.size > 1 else 0.1
    dy = np.median(np.diff(y)) if y.size > 1 else -0.1
    s = np.linspace(-0.5, 0.5, 9)
    out = []
    for yc in y:
        for xc in x:
            px = np.r_[xc + s * dx, np.full(9, xc + 0.5 * dx),
                       xc - s * dx, np.full(9, xc - 0.5 * dx)]
            py = np.r_[np.full(9, yc - 0.5 * dy), yc + s * dy,
                       np.full(9, yc + 0.5 * dy), yc - s * dy]
            out.append(np.column_stack(tr.transform(px, py)))
    return out


def overlap(g, pix_polys, fine=FINE):
    """``(npix, ncpl)`` m2 the pixels share with the active cells, from both
    cut on a ``fine`` grid over the catchment."""
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    from shapely.geometry import Polygon
    act = ~g['off'].all(axis=0)
    xy = np.concatenate([g['polys'][c][:, :2] for c in np.flatnonzero(act)])
    x0 = np.floor(xy[:, 0].min() / fine) * fine
    y1 = np.ceil(xy[:, 1].max() / fine) * fine
    nx = int(np.ceil((xy[:, 0].max() - x0) / fine))
    ny = int(np.ceil((y1 - xy[:, 1].min()) / fine))
    T = from_origin(x0, y1, fine, fine)
    cid = rasterize([(Polygon(p[:, :2]), c + 1) for c, p in
                     enumerate(g['polys']) if act[c]], out_shape=(ny, nx),
                    transform=T, fill=0, dtype='int32') - 1
    pid = rasterize([(Polygon(p), i + 1) for i, p in enumerate(pix_polys)],
                    out_shape=(ny, nx), transform=T, fill=0,
                    dtype='int32') - 1
    m = (cid >= 0) & (pid >= 0)
    n = len(pix_polys)
    w = np.zeros((n, g['ncpl']))
    np.add.at(w, (pid[m], cid[m]), fine * fine)
    return w


def catchment_geojson(info, g=None, tol=5.0):
    """The catchment -- the model's active cells -- as GeoJSON in lon/lat
    (WGS84), its outline simplified to ``tol`` m: the area to upload in the
    Copernicus Browser, whose Statistical Info then exports the mean over
    it."""
    from pyproj import Transformer
    from shapely.geometry import MultiPolygon, Polygon, mapping
    from shapely.geometry.polygon import orient
    from shapely.ops import transform, unary_union
    if not info.get('crs'):
        raise ValueError('the run needs its CRS')
    g = g or DM._grid(info['sim_ws'], info['name'])
    act = ~g['off'].all(axis=0)
    # a hair of buffer closes the slivers between neighbouring cells
    shp = unary_union([Polygon(g['polys'][c][:, :2]) for c in
                       np.flatnonzero(act)]).buffer(0.5).buffer(-0.5)
    shp = shp.simplify(tol)
    tr = Transformer.from_crs(info['crs'], 'EPSG:4326', always_xy=True)
    ll = transform(lambda x, y, z=None: tuple(np.round(tr.transform(x, y), 6)),
                   shp)
    parts = list(ll.geoms) if isinstance(ll, MultiPolygon) else [ll]
    ll = [orient(p, 1.0) for p in parts]          # RFC 7946: anticlockwise
    geom = mapping(ll[0] if len(ll) == 1 else MultiPolygon(ll))
    return json.dumps({'type': 'FeatureCollection', 'features': [
        {'type': 'Feature', 'properties': {'name': '%s catchment' % info['name']},
         'geometry': geom}]})


def model_se(info, g):
    """Se of every soil column every day: ``(nper, ncpl)``, NaN off the
    soil columns."""
    import h5py
    soil = info.get('soil')
    if not soil:
        raise ValueError("Se needs the run's soil horizons, from its "
                         "configuration copy")
    with h5py.File(info['h5'], 'r') as h:
        node = DM._nodes(h['cell_ij'][...], g)
        ss = np.asarray(h['ssoil'][...], float)          # (nper, ncell)
        tl = DM._thickness(h)
    ssat, sres = DM._column_storages(tl, soil)
    den = np.where(ssat > sres, ssat - sres, np.nan)
    se = np.full((ss.shape[0], g['ncpl']), np.nan)
    se[:, node] = np.clip((ss - sres[None]) / den[None], 0.0, 1.0)
    return se


def _wmean(v, w):
    """Weighted mean of ``v`` (..., n) with weights ``w`` (n,), NaN
    ignored; NaN where less than half the weight is valid."""
    ok = np.isfinite(v)
    sw = (ok * w).sum(axis=-1)
    s = np.where(ok, v, 0.0) @ w
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(sw >= 0.5 * w.sum(), s / sw, np.nan)


# ------------------------------------------------------------------ #
# filter and scores
# ------------------------------------------------------------------ #

def exp_filter(times, v, T):
    """The Soil Water Index of a surface series: the recursive exponential
    filter (Albergel et al. 2008), K_n = K_{n-1} / (K_{n-1} +
    exp(-(t_n - t_{n-1}) / T)), SWI_n = SWI_{n-1} + K_n (v_n - SWI_{n-1});
    on the days of ``times``, carried over the days without a value."""
    t = np.array([(d - times[0]).days for d in times], float)
    out = np.full(len(v), np.nan)
    k, s, tp = 1.0, np.nan, None
    for i, (ti, vi) in enumerate(zip(t, v)):
        if np.isfinite(vi):
            if tp is None:
                s, k = vi, 1.0
            else:
                k = k / (k + np.exp(-(ti - tp) / float(T)))
                s = s + k * (vi - s)
            tp = ti
        out[i] = s
    return out


def _anom(times, v, half=15):
    """``v`` less its day-of-year climatology (a +-``half``-day window over
    the years)."""
    doy = np.array([d.timetuple().tm_yday for d in times])
    clim = np.full(len(v), np.nan)
    for i, dd in enumerate(doy):
        dist = np.minimum(np.abs(doy - dd), 366 - np.abs(doy - dd))
        m = (dist <= half) & np.isfinite(v)
        if m.any():
            clim[i] = v[m].mean()
    return v - clim


def metrics(times, x, y):
    """Model ``x`` against product ``y`` (same unit): Pearson r, the r of
    their anomalies, the bias (y - x) and the unbiased RMSE, n days."""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return {'n': int(m.sum())}
    xa, ya = _anom(times, np.where(m, x, np.nan)), _anom(times, np.where(m, y, np.nan))
    ma = np.isfinite(xa) & np.isfinite(ya)
    d = (y[m] - y[m].mean()) - (x[m] - x[m].mean())
    return {'n': int(m.sum()), 'r': float(np.corrcoef(x[m], y[m])[0, 1]),
            'r_anom': float(np.corrcoef(xa[ma], ya[ma])[0, 1])
            if ma.sum() > 2 else None,
            'bias': float(np.mean(y[m] - x[m])), 'ubrmse': float(np.sqrt(np.mean(d ** 2)))}


def _scaled(v):
    """``v`` from its 2nd to its 98th percentile, as 0..100."""
    if not np.isfinite(v).any():
        return v
    lo, hi = np.nanpercentile(v, [2, 98])
    return np.clip((v - lo) / max(hi - lo, 1e-9), 0, 1) * 100.0


def hydro_window(dates, start_month=10):
    """The days of the run's whole hydrological years (each from day 1 of
    ``start_month``), as a mask: the window of the ET scores
    (marmites_postprocess._et_window); every day when the run holds no
    whole year."""
    m = int(start_month)
    nxt = [d + timedelta(days=1) for d in dates]
    first = next((i for i, d in enumerate(dates)
                  if d.month == m and d.day == 1), None)
    last = next((i for i in range(len(dates) - 1, -1, -1)
                 if nxt[i].month == m and nxt[i].day == 1), None)
    out = np.zeros(len(dates), bool)
    if first is None or last is None or last < first:
        out[:] = True
    else:
        out[first:last + 1] = True
    return out


# ------------------------------------------------------------------ #
# the comparison
# ------------------------------------------------------------------ #

def compare(out_dir, product, d, info=None, T=T_DEFAULT):
    """The run's catchment-mean Se against ``product`` read from folder
    ``d``, over the run's whole hydrological years; returns ``{'png',
    'csv', 'json', 'stats'}``. Raises ValueError with the reason when it
    cannot."""
    info = info or DM.outputs_of(out_dir)
    if info.get('why'):
        raise ValueError(info['why'])
    if not info.get('h5'):
        raise ValueError('the run needs its coupled results')
    P = PRODUCTS[product]
    prod = read(product, d)
    aoi = bool(prod.get('aoi'))
    g = DM._grid(info['sim_ws'], info['name'])
    if aoi:
        # the Browser's mean over the area drawn there -- the catchment --
        # is one pixel that covers it
        pp, W = [], None
        wp, keep = np.ones(1), np.ones(1, bool)
    else:
        if not info.get('crs'):
            raise ValueError('the run needs its CRS')
        pp = _pixel_polys(prod['x'], prod['y'], prod['crs'], info['crs'])
        W = overlap(g, pp)                                # (npix, ncpl)
        wp = W.sum(axis=1)                                # per pixel, m2
        keep = wp > 0
        if not keep.any():
            raise ValueError('no pixel of %s touches the catchment'
                             % P['short'])
    dates = list(info['dates'])
    hw = hydro_window(dates, info.get('hydro_year_start', 10))
    wd = [t for t, k in zip(dates, hw) if k]
    se = model_se(info, g)                                # (nper, ncpl)
    area = g['area']
    act = ~g['off'].all(axis=0) & np.isfinite(se).any(axis=0)
    se_c = _wmean(se[:, act], area[act])                  # catchment, (nper,)
    # the product on the run's days: each pixel by its share of the catchment
    pos = {t: i for i, t in enumerate(prod['times'])}
    ix = [pos.get(t) for t in dates]
    on = np.array([i is not None for i in ix])
    if not on[hw].any():
        raise ValueError('%s has no day in the whole hydrological years of '
                         'this run (%s to %s; the files: %s to %s)'
                         % (P['short'], wd[0], wd[-1], prod['times'][0],
                            prod['times'][-1]))
    flat = {k: v.reshape(len(prod['times']), -1)
            for k, v in prod['series'].items()}
    rows = [i for i in ix if i is not None]

    def catch(k):
        out = np.full(len(dates), np.nan)
        out[on] = _wmean(flat[k][rows][:, keep], wp[keep])
        return out

    def per_pixel(k):
        """Each pixel's own model series (its cells by the area they
        share) against its own value: r and days, for a product of
        several pixels (Casa de Lobos at 1 km)."""
        out = []
        for p in np.flatnonzero(keep):
            m = W[p] > 0
            xm = 100.0 * _wmean(se[:, m], W[p, m])
            ym = np.full(len(dates), np.nan)
            ym[on] = flat[k][rows, p]
            if P['kind'] != 'swi':
                ym = _scaled(exp_filter(dates, ym, T)[hw])
                xm = _scaled(xm[hw])
            else:
                xm, ym = xm[hw], ym[hw]
            out.append(dict(metrics(wd, xm, ym), pixel=int(p),
                            share=float(wp[p] / wp[keep].sum())))
        return out
    stats = {'product': product, 'label': P['label'], 'unit': P['unit'],
             'depth': P['depth'], 'run': os.path.basename(os.path.abspath(out_dir)),
             'window': [str(wd[0]), str(wd[-1])],
             'days_with_product': int(on[hw].sum())}
    if aoi:
        stats.update(source='Copernicus Browser export: the mean over the '
                     'area drawn there', pixels=None)
    else:
        stats.update(source='product files', pixels=int(keep.sum()),
                     pixel_shares=[float(s) for s in wp[keep] / wp[keep].sum()])
    x = 100.0 * se_c
    band = None
    if P['kind'] == 'swi':
        ts = sorted(k for k in flat if isinstance(k, int))
        if not ts:
            raise ValueError('no SWI T-value in the files')
        by_t = {t: metrics(wd, x[hw], catch(t)[hw]) for t in ts}
        T = T if T in ts else min(ts, key=lambda t: abs(t - T))
        y = catch(T)
        stats.update(T=T, by_T={str(t): m for t, m in by_t.items()},
                     comparison='100 Se against the SWI (both % of saturation)')
        xs, ys = x[hw], y[hw]
        if aoi:
            band = np.full((len(dates), 2), np.nan)
            band[on] = prod['spread'][T][rows]
            band = band[hw]
    else:
        if 's' not in flat:
            raise ValueError('no surface soil moisture variable in the files')
        y = exp_filter(dates, catch('s'), T)
        stats.update(T=T, comparison=(
            'the surface moisture filtered to T = %d days, then the model and '
            'the product each scaled between their 2nd and 98th percentiles'
            % T))
        xs, ys = _scaled(x[hw]), _scaled(y[hw])
    stats['fit'] = metrics(wd, xs, ys)
    if keep.sum() > 1:
        stats['per_pixel'] = per_pixel(T if P['kind'] == 'swi' else 's')
    dest = os.path.join(out_dir, '_output', SUBDIR)
    os.makedirs(dest, exist_ok=True)
    stem = 'smrs_%s_T%d' % (product, T)
    import pandas as pd
    df = pd.DataFrame({'date': wd, 'model_se_pct': x[hw], 'product': y[hw],
                       'model_compared': xs, 'product_compared': ys})
    if band is not None:
        df['product_p10'], df['product_p90'] = band[:, 0], band[:, 1]
    csv = os.path.join(dest, stem + '.csv')
    df.to_csv(csv, index=False, float_format='%.4f')
    png = os.path.join(dest, stem + '.png')
    _figure(png, g, pp, keep, wp, wd, xs, ys, stats, P, band)
    js = os.path.join(dest, stem + '.json')
    with open(js, 'w', encoding='utf-8') as fh:
        json.dump(stats, fh, indent=1)
    return {'png': png, 'csv': csv, 'json': js, 'stats': stats}


def _figure(png, g, pp, keep, wp, dates, xs, ys, st, P, band=None):
    """Four panels: the pixels over the catchment (or, for a Browser
    export, the catchment the area drawn there stands for); the two series,
    the product's spread over the area as a band; one against the other;
    the scores by T-value (an SWI) or by month."""
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    act = ~g['off'].all(axis=0)
    fig = plt.figure(figsize=(17, 10))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.1, 2.2, 1.1])
    a = fig.add_subplot(gs[:, 0])
    cells = [g['polys'][c][:, :2] for c in np.flatnonzero(act)]
    a.add_collection(PolyCollection(cells, facecolors='#dddddd',
                                    edgecolors='none'))
    if pp:
        sh = wp / max(wp.sum(), 1)
        for i, p in enumerate(pp):
            if keep[i]:
                a.plot(p[:, 0], p[:, 1], color=_SAT, lw=1.2)
                c = p.mean(axis=0)
                a.text(c[0], c[1], '%.0f %%' % (100 * sh[i]), color=_SAT,
                       ha='center', va='center', fontsize=9)
        xy = np.concatenate([pp[i] for i in np.flatnonzero(keep)])
        a.set_title('%d pixel(s) over the catchment\n(their share of it)'
                    % st['pixels'], fontsize=10, color=_INK)
    else:
        xy = np.concatenate(cells)
        a.set_title('The mean over the area drawn\nin the Copernicus Browser',
                    fontsize=10, color=_INK)
        a.text(0.5, 0.02, 'Statistical Info export: its outline\nis not in '
               'the file -- draw the catchment', transform=a.transAxes,
               ha='center', va='bottom', fontsize=8, color=_MUTED)
    a.set_xlim(xy[:, 0].min(), xy[:, 0].max())
    a.set_ylim(xy[:, 1].min(), xy[:, 1].max())
    a.set_aspect('equal', adjustable='box')
    a.tick_params(labelsize=7, colors=_MUTED)
    a.ticklabel_format(useOffset=False, style='plain')
    unit = ('% of saturation' if P['kind'] == 'swi'
            else 'scaled 2nd-98th percentile, %')
    a = fig.add_subplot(gs[0, 1])
    if band is not None:
        a.fill_between(dates, band[:, 0], band[:, 1], color=_SAT, alpha=0.18,
                       lw=0, label='%s, 10th-90th percentile over the area'
                       % P['short'])
    a.plot(dates, xs, color=_MODEL, lw=1.2, label='model: 100 Se, catchment mean')
    a.plot(dates, ys, color=_SAT, lw=1.0, label='%s%s' % (
        P['short'], ', T = %d d' % st['T'] if P['kind'] == 'swi'
        else ', filtered T = %d d' % st['T']))
    a.set_ylabel(unit, fontsize=9)
    a.legend(loc='upper right', fontsize=8)
    a.grid(alpha=0.3)
    a.tick_params(labelsize=8)
    f = st['fit']
    a.set_title('%s to %s: r %.2f, anomalies r %s, bias %+.1f, ubRMSE %.1f '
                '(%d days)' % (st['window'][0], st['window'][1],
                               f.get('r', np.nan), '%.2f' % f['r_anom']
                               if f.get('r_anom') is not None else '-',
                               f.get('bias', np.nan), f.get('ubrmse', np.nan),
                               f.get('n', 0)),
                fontsize=10, color=_INK)
    a = fig.add_subplot(gs[1, 1])
    m = np.isfinite(xs) & np.isfinite(ys)
    a.scatter(xs[m], ys[m], s=5, color=_MODEL, alpha=0.4, lw=0)
    lim = (0, 100)
    a.plot(lim, lim, color=_MUTED, lw=0.8, ls='--')
    a.set_xlim(*lim)
    a.set_ylim(*lim)
    a.set_xlabel('model [%s]' % unit, fontsize=9)
    a.set_ylabel('%s [%s]' % (P['short'], unit), fontsize=9)
    a.grid(alpha=0.3)
    a.tick_params(labelsize=8)
    a = fig.add_subplot(gs[:, 2])
    if P['kind'] == 'swi':
        ts = sorted(int(t) for t in st['by_T'])
        r = [st['by_T'][str(t)].get('r', np.nan) for t in ts]
        ra = [st['by_T'][str(t)].get('r_anom') or np.nan for t in ts]
        k = np.arange(len(ts))
        a.bar(k - 0.2, r, 0.4, color=_MODEL, label='r')
        a.bar(k + 0.2, ra, 0.4, color=_SAT, label='r of the anomalies')
        a.set_xticks(k)
        a.set_xticklabels(['T%d' % t for t in ts], fontsize=8)
        a.set_ylim(min(0, np.nanmin(r + ra) if len(ts) else 0), 1)
        a.set_title('By SWI T-value (days)', fontsize=10, color=_INK)
        a.legend(fontsize=8)
    else:
        mo = np.array([d_.month for d_ in dates])

        def by_month(v):
            out = []
            for i in range(1, 13):
                w = v[(mo == i) & np.isfinite(v)]
                out.append(w.mean() if w.size else np.nan)
            return out
        cx, cy = by_month(xs), by_month(ys)
        a.plot(range(1, 13), cx, 'o-', color=_MODEL, label='model')
        a.plot(range(1, 13), cy, 's-', color=_SAT, label=P['short'])
        a.set_xticks(range(1, 13))
        a.set_title('Mean by month', fontsize=10, color=_INK)
        a.legend(fontsize=8)
    a.grid(alpha=0.3)
    a.tick_params(labelsize=8)
    fig.suptitle('Soil moisture: the model\'s saturation Se against %s (%s, %s) '
                 '-- run %s' % (P['label'], P['depth'], P['unit'], st['run']),
                 fontsize=12, color=_INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(png, dpi=110, metadata={'MM-panels': '3'})
    plt.close(fig)
