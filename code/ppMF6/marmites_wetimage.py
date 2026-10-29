# -*- coding: utf-8 -*-
"""The model's wetness W against the NDVI of a satellite scene of the same
day (user, 2026-10-29), by linear regression -- for Casa de Lobos, tried
on La Mata.

The scene is a product of marmites_rs: Landsat 5, 7, 8 or 9 at 30 m or
Sentinel-2 at 10 m, its NDVI and its true colour on the model's CRS. The
true colour is only the visual support; the comparison is NDVI = a + b W.

PER MODEL CELL (the user's choice): W is one value per cell, so the NDVI is
averaged over each cell -- the scene's pixels cut on a 2.5 m grid, so a
cell smaller than a pixel takes the pixel it lies in. A cell counts when
half of it is clear; the ponds (open water), the irrigated fields (green
because watered) and the cells of negative NDVI (water) are left out.

PER VEGETATION CLASS: the NDVI of a holm oak's crown, of a dried pasture
and of bare rock differ whatever the water, so one regression over all of
them would mostly fit the land cover. Each cell's class is its main
vegetation type in the run's own vegetation overlay (the cached one,
marmites_props.veg_cover) when it covers at least 10 % of the cell (the
panel can change it), else bare soil; one weighted
least squares per class of 30 cells or more, the weights the cells' areas
(half the La Mata cells are the small riparian ones), and one over all.
The residual, NDVI minus the class's fit, maps where the model is wetter
or drier than the vegetation says.

Writes ``<out>/_output/day_maps/wetndvi_<platform>_<YYYY-MM-DD>.png`` and
its ``.json``; Streamlit-free, the Results panel only asks and shows.
"""

import glob
import json
import os
from types import SimpleNamespace

import numpy as np

import marmites_daymap as DM

__all__ = ['compare', 'veg_classes', 'wls', 'PREFIX']

PREFIX = 'wetndvi_'
FINE = 2.5                                  # m, the cells cut on this grid
MIN_CELLS = 30                  # a class's own fit (17 pine cells gave noise)
MIN_COVER = 10.0                # %, a cell's main vegetation type, else bare
_INK, _MUTED = '#222222', '#666666'
# one colour per class, in a fixed order (validated against colour-blind
# confusion for adjacent pairs); bare soil last, in grey
_CLASS_COLOURS = ('#008300', '#2a78d6', '#eda100', '#6250d6', '#e34948',
                  '#eb6834')
_BARE = '#7f7f7f'


def wls(x, y, w):
    """Weighted least squares y = a + b x: ``{a, b, r2, rmse, n}``."""
    w = np.asarray(w, float)
    sw = w.sum()
    xm, ym = (w * x).sum() / sw, (w * y).sum() / sw
    sxx = (w * (x - xm) ** 2).sum()
    b = (w * (x - xm) * (y - ym)).sum() / sxx if sxx > 0 else 0.0
    a = ym - b * xm
    r = y - (a + b * x)
    syy = (w * (y - ym) ** 2).sum()
    return {'a': float(a), 'b': float(b),
            'r2': float(1 - (w * r ** 2).sum() / syy) if syy > 0 else 0.0,
            'rmse': float(np.sqrt((w * r ** 2).sum() / sw)), 'n': int(x.size)}


def veg_classes(info, g, ds_dir=None, min_cover=None):
    """``(names, cls)``: the vegetation types and bare soil, and each
    cell's class (``dominant``) -- from the run's own overlay, cached in
    its workspace when it ran; None when it cannot be found (the
    comparison is then one fit over all)."""
    if min_cover is None:
        min_cover = MIN_COVER
    if not info.get('config'):
        return None
    try:
        import marmites_config as mcfg
        import marmites_props as props
        cfg = mcfg.load_run_config(info['config'])
        if ds_dir is None:
            import mm_paths
            ds_dir = mm_paths.dataset_dir(cfg.paths.case)
        names = [str(v.name) for v in cfg.surface.vegetation]
        sig = None
        for p in glob.glob(os.path.join(info['sim_ws'], '_mesh',
                                        'mesh_*.sig.json')):
            with open(p, encoding='utf-8') as fh:
                js = json.load(fh)
            if int(js.get('ncpl', -1)) == g['ncpl']:
                sig = js.get('signature')
        if sig is None or not names:
            return None
        # the run's key (its mesh, files and class table): the cached
        # overlay, never recomputed here -- an empty cell list fails fast
        cov = props.veg_cover(cfg, SimpleNamespace(nrow=g['ncpl'], ncol=1),
                              ds_dir, len(names), verbose=False, cells=[],
                              cells_key=sig, cache_dir=os.path.join(
                                  info['sim_ws'], '_overlay'))
    except Exception:                                # noqa: BLE001
        return None
    return dominant(cov, names, min_cover)


def dominant(cov, names, min_cover=MIN_COVER):
    """``(names + ['bare soil'], class of each cell)``: its main vegetation
    type -- the largest cover (% of the cell) -- when that cover is
    ``min_cover`` or more, else bare soil. The uncovered rest is not a
    rival class: on La Mata it is the larger part of nine cells in ten
    (grass covers a median 26 % of a cell), so taking the largest of all
    made nearly every cell bare soil (user, 2026-10-29)."""
    cov = np.asarray(cov, float).reshape(len(names), -1)
    top = cov.max(axis=0)
    cls = np.where(top >= float(min_cover), cov.argmax(axis=0), len(names))
    return list(names) + ['bare soil'], cls


def _fine_cells(g, transform, shape, f):
    """Each 2.5 m pixel's cell (-1 outside the active cells), on the
    product's grid cut ``f`` times finer."""
    from rasterio.features import rasterize
    from rasterio.transform import Affine
    from shapely.geometry import Polygon
    t = transform * Affine.scale(1.0 / f)
    act = ~g['off'].all(axis=0)
    return rasterize([(Polygon(p[:, :2]), c + 1) for c, p in
                      enumerate(g['polys']) if act[c]],
                     out_shape=(shape[0] * f, shape[1] * f), transform=t,
                     fill=0, dtype='int32') - 1, t


def compare(out_dir, product, info=None, gis=None, ds_dir=None,
            min_cover=MIN_COVER):
    """NDVI of ``product`` (a record of marmites_rs.products) against the
    run's W on the scene's day; returns ``{'png', 'json', 'stats'}``.
    ``min_cover``: the % of a cell its main vegetation type must cover to
    class it, else bare soil. Raises ValueError with the reason when it
    cannot."""
    import rasterio
    from rasterio.features import rasterize
    info = info or DM.outputs_of(out_dir)
    if info.get('why'):
        raise ValueError(info['why'])
    if not info.get('h5'):
        raise ValueError('the coupled results are not in %s' % info['sim_ws'])
    import marmites_rs as RS
    day = RS.as_date(product['date'])
    dates = list(info['dates'])
    if day not in dates:
        raise ValueError('the scene of %s is not a day of this run (%s to %s)'
                         % (day, dates[0], dates[-1]))
    k, nper = dates.index(day), len(dates)
    sim_ws, name = info['sim_ws'], info['name']
    g = DM._grid(sim_ws, name)
    with rasterio.open(os.path.join(product['dir'], product['ndvi'])) as src:
        nd = src.read(1).astype(float)
        T, crs = src.transform, src.crs
    with rasterio.open(os.path.join(product['dir'], product['rgb'])) as src:
        rgb = src.read()
    if info.get('crs') and crs is not None and \
            crs.to_epsg() != int(str(info['crs']).split(':')[-1]):
        raise ValueError('the product is in %s, the run in %s: fetch it again '
                         'for this run' % (crs, info['crs']))
    res = abs(T.a)
    f = max(1, int(round(res / FINE)))
    cid, Tf = _fine_cells(g, T, nd.shape, f)
    inside = cid >= 0
    ndf = np.repeat(np.repeat(nd, f, axis=0), f, axis=1)
    ok = inside & np.isfinite(ndf)
    n = g['ncpl']
    tot = np.bincount(cid[inside], minlength=n).astype(float)
    nok = np.bincount(cid[ok], minlength=n).astype(float)
    sm = np.bincount(cid[ok], weights=ndf[ok], minlength=n)
    with np.errstate(invalid='ignore', divide='ignore'):
        NDc = np.where(nok > 0, sm / nok, np.nan)
        clear = np.where(tot > 0, nok / tot, 0.0)
    soil = info.get('soil')
    W = DM.day_values(sim_ws, name, k, nper, 'W', g=g, h5=info['h5'],
                      soil=soil)[0]
    ponds = np.zeros(n, bool)
    for cells in DM.pond_cells(sim_ws, name, g).values():
        ponds[cells] = True
    irr_cell = np.zeros(n, bool)
    irr_px = None
    layer = info.get('irr_layer')
    if layer and gis and os.path.exists(os.path.join(str(gis), layer)):
        import geopandas as gpd
        df = gpd.read_file(os.path.join(str(gis), layer))
        if df.crs is not None:
            df = df.to_crs(crs)
        geoms = [gm for gm in df.geometry if gm is not None and not gm.is_empty]
        if geoms:
            irr_px = rasterize([(gm, 1) for gm in geoms], out_shape=cid.shape,
                               transform=Tf, fill=0, dtype='uint8').astype(bool)
            ir = np.bincount(cid[inside & irr_px], minlength=n)
            irr_cell = np.where(tot > 0, ir / np.maximum(tot, 1), 0) >= 0.5
    vc = veg_classes(info, g, ds_dir=ds_dir, min_cover=min_cover)
    if vc is None:
        names, cls = ['all cells'], np.zeros(n, int)
    else:
        names, cls = vc
    area = g['area']
    use = ((tot > 0) & (clear >= 0.5) & np.isfinite(W) & np.isfinite(NDc)
           & ~ponds & ~irr_cell & (NDc >= 0))
    fits = {}
    for c, nm in enumerate(names):
        m = use & (cls == c)
        if m.sum() >= MIN_CELLS:
            fits[nm] = dict(wls(W[m], NDc[m], area[m]),
                            area_km2=float(area[m].sum() / 1e6))
    pooled = dict(wls(W[use], NDc[use], area[use]),
                  area_km2=float(area[use].sum() / 1e6)) if use.any() else None
    if pooled is None:
        raise ValueError('no cell of the catchment is clear on this scene')
    fit_a = np.full(n, pooled['a'])
    fit_b = np.full(n, pooled['b'])
    for c, nm in enumerate(names):
        if nm in fits:
            fit_a[cls == c], fit_b[cls == c] = fits[nm]['a'], fits[nm]['b']
    resid = np.where(use, NDc - (fit_a + fit_b * W), np.nan)
    act = ~g['off'].all(axis=0)
    stats = {'run': os.path.basename(os.path.abspath(out_dir)),
             'day': day.isoformat(), 'period': k + 1,
             'scene': {kk: product.get(kk) for kk in (
                 'id', 'sensor', 'platform', 'date', 'res', 'scene_cloud')},
             'clear_share_catchment': float((nok[act].sum() / max(tot[act].sum(), 1))),
             'cells': {'catchment': int(act.sum()), 'compared': int(use.sum()),
                       'cloudy_or_missing': int((act & (clear < 0.5)).sum()),
                       'ponds': int((act & ponds).sum()),
                       'irrigated': int((act & irr_cell).sum()),
                       'water_ndvi_below_0': int((act & (NDc < 0)).sum())},
             'classes': ('the main vegetation type of the run\'s overlay, '
                         'covering at least %g %% of the cell, else bare '
                         'soil' % min_cover if vc is not None
                         else 'none found: one fit over all'),
             'min_cover': float(min_cover),
             'class_cells': {nm: int((act & (cls == c)).sum())
                             for c, nm in enumerate(names)},
             'fits': fits, 'all_classes': pooled}
    dest = os.path.join(out_dir, '_output', DM.DAY_DIR)
    os.makedirs(dest, exist_ok=True)
    stem = '%s%s_%s' % (PREFIX, product['platform'], day.isoformat())
    png = os.path.join(dest, stem + '.png')
    _figure(png, g, rgb, nd, T, cid, f, irr_px, W, NDc, use, cls, names, fits,
            pooled, resid, stats)
    js = os.path.join(dest, stem + '.json')
    with open(js, 'w', encoding='utf-8') as fh:
        json.dump(stats, fh, indent=1)
    return {'png': png, 'json': js, 'stats': stats}


def _figure(png, g, rgb, nd, T, cid, f, irr_px, W, NDc, use, cls, names,
            fits, pooled, resid, st):
    """Six panels: the scene's true colour (support), its NDVI, the model's
    W, the NDVI per cell, NDVI against W per class, the residuals."""
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    from rasterio.transform import array_bounds
    ny, nx = nd.shape
    west, south, east, north = array_bounds(ny, nx, T)
    ext = (west, east, south, north)
    act = ~g['off'].all(axis=0)
    inside = cid >= 0
    xs = west + (np.arange(cid.shape[1]) + 0.5) * (east - west) / cid.shape[1]
    ys = north - (np.arange(cid.shape[0]) + 0.5) * (north - south) / cid.shape[0]
    xy = np.concatenate([g['polys'][c][:, :2] for c in np.flatnonzero(act)])
    pad = 60.0
    xlim = (xy[:, 0].min() - pad, xy[:, 0].max() + pad)
    ylim = (xy[:, 1].min() - pad, xy[:, 1].max() + pad)
    polys = [g['polys'][c][:, :2] for c in np.flatnonzero(act)]
    idx = np.flatnonzero(act)

    def outline(ax, lw=0.9):
        ax.contour(xs, ys, inside.astype(float), levels=[0.5], colors=_INK,
                   linewidths=lw)
        if irr_px is not None:
            ax.contour(xs, ys, irr_px.astype(float), levels=[0.5],
                       colors=_MUTED, linewidths=0.8, linestyles='--')

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

    def cells(ax, v, cmap, norm):
        pc = PolyCollection([polys[i] for i in range(len(idx))],
                            array=np.ma.masked_invalid(v[idx]), cmap=cmap,
                            norm=norm, edgecolors='none')
        ax.add_collection(pc)
        return pc

    sc = st['scene']
    lo, hi = np.nanpercentile(NDc[act], [2, 98])
    nnorm = Normalize(vmin=lo, vmax=hi)
    fig, ax = plt.subplots(2, 3, figsize=(18, 12.5))
    a = ax[0, 0]
    a.imshow(np.moveaxis(rgb, 0, -1), extent=ext, interpolation='nearest')
    outline(a)
    frame(a, '%s, %s: true colour (visual support)'
          % (sc['platform'].capitalize(), sc['date']))
    a = ax[0, 1]
    im = a.imshow(np.ma.masked_invalid(nd), extent=ext, cmap='YlGn', norm=nnorm,
                  interpolation='nearest')
    outline(a)
    frame(a, 'NDVI at %g m (clouds, gaps masked: %.0f %% of the catchment '
          'clear)' % (sc['res'], 100 * st['clear_share_catchment']))
    fig.colorbar(im, ax=a, shrink=0.75, label='NDVI [-]')
    a = ax[0, 2]
    pc = cells(a, W, 'Blues', Normalize(0, 1))
    outline(a)
    frame(a, 'Model: wetness W = (1 - f_ow) Se + f_ow, day %d' % st['period'])
    fig.colorbar(pc, ax=a, shrink=0.75,
                 label='W [-]: 0 residual water, 1 saturated / open water')
    a = ax[1, 0]
    pc = cells(a, np.where(use, NDc, np.nan), 'YlGn', nnorm)
    outline(a)
    frame(a, 'NDVI per model cell (%d cells compared; ponds, irrigated, '
          'cloudy out)' % st['cells']['compared'])
    fig.colorbar(pc, ax=a, shrink=0.75, label='NDVI [-]')
    a = ax[1, 1]
    x = np.linspace(0, 1, 11)
    for c, nm in enumerate(names):
        m = use & (cls == c)
        if not m.any():
            continue
        col = _BARE if nm == 'bare soil' else _CLASS_COLOURS[c % len(_CLASS_COLOURS)]
        a.scatter(W[m], NDc[m], s=4, color=col, alpha=0.35, lw=0)
        if nm in fits:
            fi = fits[nm]
            a.plot(x, fi['a'] + fi['b'] * x, color=col, lw=2,
                   label='%s: %.2f %+.2f W, R$^2$ %.2f, n %d'
                   % (nm, fi['a'], fi['b'], fi['r2'], fi['n']))
    a.plot(x, pooled['a'] + pooled['b'] * x, color=_INK, lw=1.2, ls='--',
           label='all: %.2f %+.2f W, R$^2$ %.2f, n %d'
           % (pooled['a'], pooled['b'], pooled['r2'], pooled['n']))
    a.set_xlim(0, 1)
    a.set_xlabel('model wetness W [-]', fontsize=9)
    a.set_ylabel('NDVI of the cell [-]', fontsize=9)
    a.tick_params(labelsize=7)
    a.grid(alpha=0.3)
    a.legend(loc='upper left', fontsize=8, frameon=True)
    a.set_title('NDVI = a + b W per vegetation class (main type >= %g %% of '
                'the cell), weighted by cell area' % st['min_cover'],
                fontsize=10, color=_INK)
    a = ax[1, 2]
    r = resid[act][np.isfinite(resid[act])]
    lim = float(np.percentile(np.abs(r), 98)) if r.size else 0.1
    pc = cells(a, resid, 'BrBG', Normalize(-lim, lim))
    outline(a)
    frame(a, 'Residual NDVI - fit: green, greener than its wetness predicts; '
          'brown, browner')
    fig.colorbar(pc, ax=a, shrink=0.75, label='NDVI residual [-]')
    fig.suptitle('Model wetness W (run %s, %s) against the NDVI of %s %s'
                 % (st['run'], st['day'], sc['platform'].capitalize(),
                    sc.get('id') or ''), fontsize=13, color=_INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    # a wide figure: the Results panel gives it a row of its own
    fig.savefig(png, dpi=110, metadata={'MM-panels': '3'})
    plt.close(fig)
