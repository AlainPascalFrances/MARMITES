# -*- coding: utf-8 -*-
"""The model's wetness against a georeferenced image of the same day (user,
2026-10-29: "Make an option to produce the map wet_vs_rgb_20130615.png"),
a qualitative check of WHERE the model is wet.

Image: an RGB(A) GeoTIFF in any CRS -- a Google Earth export, a true-colour
Landsat scene -- reprojected onto the model's CRS at 2.5 m. Green is the
excess green ExG = 2g - r - b of the chromatic coordinates above its Otsu
threshold over the catchment, the ponds left out (open water is dark). In
a Mediterranean June the annual grass has dried off, so what stays green
is fed by groundwater or irrigated; the holm-oak crowns, dark and
grey-green, stay below the threshold.

Model: the wetness W of that day (marmites_daymap), on its cells.

Both on 50 m blocks. A block is image-wet when half of it is green, and the
model's wettest blocks over the same area are model-wet (equal areas, so
the comparison is of the pattern, not of a threshold): hit rate, false
alarms, Cohen's kappa, also within one block (50 m), and the Spearman rank
correlation of the green share against W, Se, theta and the exfiltration
-- with and without the irrigated fields, which are green when watered but
wet only to field capacity.

Writes ``<out>/_output/day_maps/wetimage_<image>_<YYYY-MM-DD>.png`` and its
``.json``; Streamlit-free, the Results panel only asks and shows.
"""

import glob
import json
import os
import re
from datetime import date, datetime

import numpy as np

import marmites_daymap as DM

__all__ = ['images', 'image_date', 'compare', 'PREFIX']

PREFIX = 'wetimage_'
_INK, _MUTED = '#222222', '#666666'
# the agreement classes: dry in both, wet in both, model only, image only
# (validated against colour-blind confusion, adjacent pairs)
_CLASSES = ('#e6e6e6', '#2a78d6', '#e34948', '#eda100')


def images(*folders):
    """The GeoTIFFs at the top of ``folders``, the newest first."""
    out = set()
    for d in folders:
        if d and os.path.isdir(str(d)):
            for pat in ('*.tif', '*.tiff', '*.TIF', '*.TIFF'):
                out.update(os.path.abspath(p) for p in
                           glob.glob(os.path.join(str(d), pat)))
    return sorted(out, key=os.path.getmtime, reverse=True)


def image_date(path):
    """The date a file name carries as YYYYMMDD (DR_RGB_20130615_...), or
    None."""
    for m in re.finditer(r'(?<!\d)(\d{4})(\d{2})(\d{2})(?!\d)',
                         os.path.basename(str(path))):
        try:
            return date(int(m.group(1)), int(m.group(2)), int(m.group(3)))
        except ValueError:
            continue
    return None


def _otsu(v, nb=256):
    """The threshold that splits ``v`` in two classes of least spread."""
    h, e = np.histogram(v, bins=nb)
    c = 0.5 * (e[1:] + e[:-1])
    w0 = np.cumsum(h)
    w1 = w0[-1] - w0
    s = np.cumsum(h * c)
    m0 = s / np.maximum(w0, 1)
    m1 = (s[-1] - s) / np.maximum(w1, 1)
    return float(c[np.argmax(w0 * w1 * (m0 - m1) ** 2)])


def _rgba(image, crs, transform, shape):
    """The image's bands 1-3 (and its alpha, else all opaque) on the pixel
    grid, uint8 (4, ny, nx)."""
    import rasterio
    from rasterio.warp import Resampling, reproject
    out = np.zeros((4,) + shape, np.uint8)
    with rasterio.open(str(image)) as src:
        if src.count < 3:
            raise ValueError('%s has %d band(s): an RGB image is needed'
                             % (os.path.basename(str(image)), src.count))
        for b in range(min(src.count, 4)):
            reproject(src.read(b + 1), out[b], src_transform=src.transform,
                      src_crs=src.crs, dst_transform=transform, dst_crs=crs,
                      resampling=Resampling.bilinear)
        if src.count < 4:
            out[3] = np.where(out[:3].sum(axis=0) > 0, 255, 0)
    return out


def _irrigated(layer, gis, crs, transform, shape):
    """The irrigated fields (the run's surface.irr_zones layer, a shapefile
    in DATA_ROOT/GIS) on the pixel grid, or None."""
    if not layer or not gis:
        return None
    p = os.path.join(str(gis), layer)
    if not os.path.exists(p):
        return None
    import geopandas as gpd
    from rasterio.features import rasterize
    df = gpd.read_file(p)
    if df.crs is not None and crs:
        df = df.to_crs(crs)
    geoms = [gm for gm in df.geometry if gm is not None and not gm.is_empty]
    if not geoms:
        return None
    return rasterize([(gm, 1) for gm in geoms], out_shape=shape,
                     transform=transform, fill=0, dtype='uint8').astype(bool)


def _spearman(a, b):
    """The rank correlation; NaN where one side is constant (a run without
    exfiltration), where it is not defined."""
    from scipy.stats import spearmanr
    if not np.ptp(a) or not np.ptp(b):
        return float('nan')
    return float(spearmanr(a, b).correlation)


def _scores(gfrac, wblk, m):
    """The equal-area comparison of the image-wet blocks (half green)
    with the model's wettest over the blocks ``m``."""
    from scipy.ndimage import binary_dilation
    img = gfrac >= 0.5
    p = float(img[m].mean())
    if p <= 0:
        raise ValueError('no block of the catchment is half green on the '
                         'image: nothing to compare')
    cut = float(np.quantile(wblk[m], 1 - p))
    mod = wblk >= cut
    a = int((img & mod & m).sum())
    b = int((~img & mod & m).sum())
    c = int((img & ~mod & m).sum())
    d = int((~img & ~mod & m).sum())
    n = a + b + c + d
    pe = ((a + b) * (a + c) + (c + d) * (b + d)) / n ** 2
    near = np.ones((3, 3), bool)
    return {'blocks': n, 'wet_share': p, 'model_cut': cut, 'both': a,
            'model_only': b, 'image_only': c, 'both_dry': d,
            'hit': a / max(a + c, 1), 'false_alarm': b / max(a + b, 1),
            'kappa': ((a + d) / n - pe) / (1 - pe) if pe < 1 else 0.0,
            'hit_50m': float((img & m & binary_dilation(mod & m, near)).sum()
                             / max((img & m).sum(), 1)),
            'false_alarm_50m': float((mod & m & ~binary_dilation(img & m, near))
                                     .sum() / max((mod & m).sum(), 1)),
            'spearman': _spearman(gfrac[m], wblk[m])}, img, mod


def compare(out_dir, image, day, info=None, gis=None, block=50.0, res=2.5):
    """Compare the run's wetness W on ``day`` with ``image``; returns
    ``{'png', 'json', 'scores'}``. Raises ValueError with the reason when
    it cannot."""
    import h5py
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    from shapely.geometry import Polygon
    info = info or DM.outputs_of(out_dir)
    if info.get('why'):
        raise ValueError(info['why'])
    if not info.get('crs'):
        raise ValueError("the run's configuration names no CRS (grid.crs_epsg)")
    if not info.get('h5'):
        raise ValueError('the coupled results are not in %s' % info['sim_ws'])
    if isinstance(day, datetime):
        day = day.date()
    dates = list(info['dates'])
    if day not in dates:
        raise ValueError('%s is not a day of this run (%s to %s)'
                         % (day, dates[0], dates[-1]))
    k, nper = dates.index(day), len(dates)
    sim_ws, name, crs = info['sim_ws'], info['name'], info['crs']
    g = DM._grid(sim_ws, name)
    active = ~g['off'].all(axis=0)
    # the pixel grid: whole blocks over the active cells
    xy = np.concatenate([g['polys'][c][:, :2] for c in np.flatnonzero(active)])
    x0 = np.floor(xy[:, 0].min() / block) * block
    y1 = np.ceil(xy[:, 1].max() / block) * block
    nbx = int(np.ceil((xy[:, 0].max() - x0) / block))
    nby = int(np.ceil((y1 - xy[:, 1].min()) / block))
    npx = int(round(block / res))
    shape = (nby * npx, nbx * npx)
    T = from_origin(x0, y1, res, res)
    cid = rasterize([(Polygon(p[:, :2]), c + 1) for c, p in enumerate(g['polys'])
                     if active[c]], out_shape=shape, transform=T, fill=0,
                    dtype='int32') - 1
    inside = cid >= 0
    at = np.maximum(cid, 0)
    rgba = _rgba(image, crs, T, shape)
    s = rgba[:3].astype(float).sum(axis=0)
    ok = (rgba[3] > 0) & (s > 0)
    r, gg, b = (np.where(s > 0, rgba[i] / np.where(s > 0, s, 1), 0.0)
                for i in range(3))
    exg = 2 * gg - r - b
    # the model on the day
    soil = info.get('soil')
    W = DM.day_values(sim_ws, name, k, nper, 'W', g=g, h5=info['h5'],
                      soil=soil)[0]
    SE = DM.saturation(info['h5'], k, g, soil)
    TH = DM.day_values(sim_ws, name, k, nper, 'theta', g=g, h5=info['h5'])[0]
    EXF = np.zeros(g['ncpl'])
    with h5py.File(info['h5'], 'r') as h:
        if 'exf' in h:
            EXF[DM._nodes(h['cell_ij'][...], g)] = h['exf'][k, :]
    ponds = np.zeros(g['ncpl'], bool)
    for cells in DM.pond_cells(sim_ws, name, g).values():
        ponds[cells] = True
    irr = _irrigated(info.get('irr_layer'), gis, crs, T, shape)
    use = inside & ok & ~ponds[at]
    if not use.any():
        raise ValueError('%s does not cover the catchment'
                         % os.path.basename(str(image)))
    v = exg[use]
    thr = _otsu(v[(v > -0.15) & (v < 0.35)])
    green = exg > thr

    def blocks(a, m):
        a = np.where(m, a, 0.0).reshape(nby, npx, nbx, npx)
        n = m.reshape(nby, npx, nbx, npx).sum(axis=(1, 3))
        return a.sum(axis=(1, 3)) / np.maximum(n, 1), n / npx ** 2

    def on_blocks(field):
        return blocks(np.nan_to_num(np.where(inside, field[at], np.nan)), use)[0]
    gfrac, cov = blocks(green.astype(float), use)
    wblk = on_blocks(W)
    valid = cov >= 0.5
    irr_blk = (blocks(irr.astype(float), inside)[0] if irr is not None
               else np.zeros_like(gfrac))
    core = valid & (irr_blk < 0.1)
    s_all, img_wet, mod_wet = _scores(gfrac, wblk, valid)
    s_core = _scores(gfrac, wblk, core)[0] if irr is not None else None
    alt = {nm: _spearman(gfrac[core], on_blocks(f)[core])
           for nm, f in (('W (Se + open water)', W), ('Se (soil only)', SE),
                         ('theta (% vol)', TH), ('exfiltration', EXF))}
    stats = {'image': os.path.abspath(str(image)), 'day': day.isoformat(),
             'period': k + 1, 'run': os.path.basename(os.path.abspath(out_dir)),
             'otsu': thr, 'green_share': float(green[use].mean()),
             'all': s_all, 'without_irrigated': s_core, 'spearman': alt,
             'irrigated_layer': info.get('irr_layer') if irr is not None else None}
    dest = os.path.join(out_dir, '_output', DM.DAY_DIR)
    os.makedirs(dest, exist_ok=True)
    stem = re.sub(r'[^\w.-]+', '_', os.path.splitext(os.path.basename(
        str(image)))[0])
    png = os.path.join(dest, '%s%s_%s.png' % (PREFIX, stem, day.isoformat()))
    _figure(png, rgba, inside, irr, cid, W, gfrac, wblk, valid, core,
            img_wet, mod_wet, stats, (x0, y1, res, npx, nbx, nby),
            'Model wetness W on %s (run %s, period %d) against the green of %s'
            % (day, stats['run'], k + 1, os.path.basename(str(image))))
    js = png[:-4] + '.json'
    with open(js, 'w', encoding='utf-8') as fh:
        json.dump(stats, fh, indent=1)
    return {'png': png, 'json': js, 'scores': stats}


def _figure(png, rgba, inside, irr, cid, W, gfrac, wblk, valid, core,
            img_wet, mod_wet, st, geo, suptitle):
    """Six panels: the image; its green per block; the model's W; the
    model's wettest blocks on the image; the agreement; green by W."""
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    x0, y1, res, npx, nbx, nby = geo
    ny, nx = inside.shape
    ext = (x0, x0 + nx * res, y1 - ny * res, y1)
    bext = (x0, x0 + nbx * npx * res, y1 - nby * npx * res, y1)
    xc = x0 + (np.arange(nx) + 0.5) * res
    yc = y1 - (np.arange(ny) + 0.5) * res
    xb = x0 + (np.arange(nbx) + 0.5) * npx * res
    yb = y1 - (np.arange(nby) + 0.5) * npx * res
    rows = np.flatnonzero(inside.any(axis=1))
    cols = np.flatnonzero(inside.any(axis=0))
    pad = npx
    xlim = (x0 + max(cols[0] - pad, 0) * res, x0 + min(cols[-1] + pad, nx) * res)
    ylim = (y1 - min(rows[-1] + pad, ny) * res, y1 - max(rows[0] - pad, 0) * res)
    rgb = np.moveaxis(rgba[:3], 0, -1).astype(float) / 255.0
    dim = np.where(inside[..., None], rgb, 0.55 * rgb + 0.45)

    def outline(ax, lw=0.9):
        ax.contour(xc, yc, inside.astype(float), levels=[0.5], colors=_INK,
                   linewidths=lw)
        if irr is not None:
            ax.contour(xc, yc, irr.astype(float), levels=[0.5], colors=_MUTED,
                       linewidths=0.8, linestyles='--')

    def frame(ax, title):
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect('equal', adjustable='box')
        ax.set_title(title, fontsize=10, color=_INK)
        ax.tick_params(labelsize=7, colors=_MUTED)
        ax.ticklabel_format(useOffset=False, style='plain')
        ax.set_xticks([t for t in ax.get_xticks()[::2]
                       if xlim[0] <= t <= xlim[1]])
        ax.set_xlim(*xlim)             # set_xticks widens the view to its ticks

    fig, ax = plt.subplots(2, 3, figsize=(18, 12.5))
    a = ax[0, 0]
    a.imshow(dim, extent=ext, interpolation='nearest')
    outline(a)
    frame(a, 'The image (catchment%s)' % ('; irrigated fields dashed'
                                          if irr is not None else ''))
    a = ax[0, 1]
    im = a.imshow(np.where(valid, gfrac, np.nan), extent=bext, cmap='Greens',
                  vmin=0, vmax=1, interpolation='nearest')
    outline(a)
    frame(a, 'Image: green share of each %g m block (ExG > %.3f, Otsu)'
          % (npx * res, st['otsu']))
    fig.colorbar(im, ax=a, shrink=0.75, label='green fraction [-]')
    a = ax[0, 2]
    im = a.imshow(np.where(inside, W[np.maximum(cid, 0)], np.nan), extent=ext,
                  cmap='Blues', vmin=0, vmax=1, interpolation='nearest')
    outline(a)
    frame(a, 'Model: wetness W = (1 - f_ow) Se + f_ow, on its cells')
    fig.colorbar(im, ax=a, shrink=0.75,
                 label='W [-]: 0 residual water, 1 saturated / open water')
    s = st['all']
    a = ax[1, 0]
    a.imshow(dim, extent=ext, interpolation='nearest')
    a.contour(xb, yb, (mod_wet & valid).astype(float), levels=[0.5],
              colors='#00e5ff', linewidths=1.4)
    outline(a, lw=0.6)
    frame(a, "The model's wettest %.0f %% of the blocks (cyan) on the image"
          % (100 * s['wet_share']))
    a = ax[1, 1]
    cls = np.full(valid.shape, np.nan)
    cls[valid & ~img_wet & ~mod_wet] = 0
    cls[valid & img_wet & mod_wet] = 1
    cls[valid & ~img_wet & mod_wet] = 2
    cls[valid & img_wet & ~mod_wet] = 3
    a.imshow(cls, extent=bext, cmap=ListedColormap(_CLASSES), vmin=-0.5,
             vmax=3.5, interpolation='nearest')
    outline(a)
    frame(a, 'Agreement (%g m blocks, equal wet areas): kappa %.2f, hit %.2f '
          '(%.2f within %g m)' % (npx * res, s['kappa'], s['hit'],
                                  s['hit_50m'], npx * res))
    a.legend(handles=[Patch(color=_CLASSES[1], label='wet in both (%d)' % s['both']),
                      Patch(color=_CLASSES[2], label='model only (%d)' % s['model_only']),
                      Patch(color=_CLASSES[3], label='image only (%d)' % s['image_only']),
                      Patch(color=_CLASSES[0], label='dry in both (%d)' % s['both_dry'])],
             loc='lower left', fontsize=8, frameon=True)
    a = ax[1, 2]
    w, gf = wblk[core], gfrac[core]
    q = np.quantile(w, np.linspace(0, 1, 11))
    k = np.clip(np.searchsorted(q, w, side='right') - 1, 0, 9)
    a.boxplot([gf[k == i] for i in range(10)], positions=range(1, 11),
              widths=0.6, showfliers=False,
              medianprops=dict(color=_CLASSES[1], lw=2),
              boxprops=dict(color=_MUTED), whiskerprops=dict(color=_MUTED),
              capprops=dict(color=_MUTED))
    a.plot(range(1, 11), [gf[k == i].mean() if (k == i).any() else np.nan
                          for i in range(10)], 'o', color=_INK, ms=4,
           label='mean')
    a.set_xticks(range(1, 11))
    a.set_xticklabels(['%.2f' % ((q[i] + q[i + 1]) / 2) for i in range(10)],
                      fontsize=7, rotation=40)
    a.set_xlabel('model wetness W, deciles of the blocks (class centre)',
                 fontsize=9)
    a.set_ylabel('image green fraction', fontsize=9)
    a.tick_params(labelsize=7)
    a.grid(alpha=0.3)
    which = ('irrigated fields out' if st['without_irrigated'] is not None
             else 'all blocks')
    a.text(0.03, 0.97, 'Spearman with the green fraction\n(%d blocks, %s):\n'
           % (core.sum(), which) + '\n'.join(
               '  %s: %.2f' % (k_, v) for k_, v in st['spearman'].items()),
           transform=a.transAxes, va='top', fontsize=8, color=_INK,
           bbox=dict(fc='white', ec='#cccccc'))
    a.set_title('Green by model wetness (%s)' % which, fontsize=10, color=_INK)
    a.legend(loc='upper right', fontsize=8)
    fig.suptitle(suptitle, fontsize=13, color=_INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    # a wide figure: the Results panel gives it a row of its own
    fig.savefig(png, dpi=110, metadata={'MM-panels': '3'})
    plt.close(fig)
