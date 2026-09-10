# -*- coding: utf-8 -*-
"""Post-processing for the coupled MARMITES / MODFLOW 6 run.

Reads the MODFLOW 6 output already on disk (heads, cell budget, UZF and SFR
budgets, the listing) and the coupled-run HDF5, and writes figures + tidy CSVs
into ``<sim_ws>/postproc/``. Mirrors the CdL post-processing
(code/SFR_LAK_CRR/postprocess_cdl.py) but for La Mata's structured DIS grid
and daily stress periods, so it needs none of the Voronoi/DISV machinery
(rasterio / shapely / geopandas / pyproj) that script depends on.

Figures / CSVs
--------------
  obs_heads          observed vs computed head time series at the piezometers
  mean_head_L*       time-mean head map per layer
  mean_depth_L*      time-mean depth-to-water map per layer (land surface - head)
  budget_compartment overall budget by compartment (m3/d) from the listing
  budget_yearly      yearly volumes (m3/yr) by compartment
  budget_uzf         UZF internal budget (recharge / seepage / ET / storage)
  budget_sfr         SFR internal budget (inflow / leakage / outflow / ...)
  storage_change_L*  per-layer storage change map

The steady-state stress period (kper 0) is excluded from every average and map.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import contextlib
import logging
import os
import sys
import time
import warnings

import numpy as np

__all__ = ['run_postproc', 'run_preproc', 'obs_points', 'obs_series',
           'budget_by_compartment', 'package_budget', 'layer_storage_change',
           'native_suite', 'COMPARTMENT']

# MARMITES input maps (soil-water-balance side) and the MODFLOW input maps
# (aquifer side). Both belong in the preprocessing channel: the MM script
# itself produces the soil/vegetation/meteo maps, and the aquifer maps come
# from the MF workspace. (label, filename-relative-to, colormap, integer?)
# compartment a listing budget term belongs to (same taxonomy as the CdL script)
COMPARTMENT = {
    'STO': 'aquifer (storage)', 'STORAGE': 'aquifer (storage)',
    'UZF': 'unsaturated (UZF->GW)', 'UZF-GWRCH': 'unsaturated (UZF->GW)',
    'SFR': 'surface (stream)', 'LAK': 'surface (lake)',
    'DRN': 'surface (drain)', 'DRN_SEEP': 'surface (seepage)',
    'WEL': 'ET (groundwater)', 'GHB': 'boundary (GHB)', 'RCH': 'recharge',
}
_COMP_COLORS = {
    'aquifer (storage)': 'tab:brown', 'unsaturated (UZF->GW)': 'tab:green',
    'surface (stream)': 'tab:blue', 'surface (lake)': 'tab:cyan',
    'surface (drain)': 'tab:purple', 'surface (seepage)': 'tab:red',
    'ET (groundwater)': 'tab:orange', 'boundary (GHB)': '0.4',
    'recharge': 'tab:olive', 'other': '0.7',
}


def _compartment_of(term):
    t = str(term).upper().strip()
    # longest key first, so DRN_SEEP wins over DRN and UZF-GWRCH over UZF
    for key in sorted(COMPARTMENT, key=len, reverse=True):
        if key in t:
            return COMPARTMENT[key]
    return 'other'


def _asc(fn):
    a = np.loadtxt(fn, skiprows=6)
    return np.where(a <= -9990.0, np.nan, a)


def _mkdir(sim_ws, sub='_output'):
    d = os.path.join(sim_ws, sub)
    os.makedirs(d, exist_ok=True)
    return d


# --------------------------------------------------------------------- #
# observation points
# --------------------------------------------------------------------- #

def obs_points(ds_ws, fn='inputObs.txt'):
    """Active observation points from inputObs.txt.

    Lines starting with '#' or '##' are disabled points and are skipped, as in
    the MARMITES input convention. Returns [{name, x, y, lay}, ...].
    """
    pts = []
    with open(os.path.join(ds_ws, fn)) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            p = s.split()
            if len(p) < 4:
                continue
            try:
                pts.append({'name': p[0], 'x': float(p[1]), 'y': float(p[2]),
                            'lay': int(p[3])})
            except ValueError:
                continue
    return pts


def obs_series(ds_ws, name, prefix='inputObsHEADS_'):
    """Observed head time series (date, head) for one point, or None."""
    import pandas as pd
    fn = os.path.join(ds_ws, '%s%s.txt' % (prefix, name))
    if not os.path.exists(fn):
        return None
    df = pd.read_csv(fn, sep=r'\s+', header=None, names=['date', 'head'],
                     engine='python')
    df['date'] = pd.to_datetime(df['date'], errors='coerce')
    return df.dropna().sort_values('date')


def _xy_to_ij(x, y, xll, yll, cs, nrow, ncol):
    # rows run north -> south: the top edge is yll + nrow*cs, row 0 is the north
    ytop = yll + nrow * cs
    j = int(np.clip((x - xll) // cs, 0, ncol - 1))
    i = int(np.clip((ytop - y) // cs, 0, nrow - 1))
    return i, j


# --------------------------------------------------------------------- #
# budgets
# --------------------------------------------------------------------- #

def budget_by_compartment(sim_ws, name):
    """Mean IN/OUT rate [m3/d] per compartment from the MF6 listing.

    Returns (df_terms, df_compartments) as pandas frames. Uses flopy's
    Mf6ListBudget so the numbers are exactly MF6's own accounting.
    """
    import pandas as pd
    import flopy
    lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, '%s.lst' % name))
    inc, cum = lst.get_dataframes()
    # drop the steady-state first period; average the incremental rates
    rate = inc.iloc[1:] if len(inc) > 1 else inc
    terms = []
    for col in rate.columns:
        if col in ('TOTAL_IN', 'TOTAL_OUT', 'IN-OUT', 'PERCENT_DISCREPANCY'):
            continue
        mean = float(rate[col].mean())
        if col.endswith('_OUT'):
            direction, base = 'OUT', col[:-4]
        elif col.endswith('_IN'):
            direction, base = 'IN', col[:-3]
        else:
            direction, base = '', col
        terms.append({'term': base, 'direction': direction, 'rate_m3d': mean,
                      'compartment': _compartment_of(base)})
    df = pd.DataFrame(terms)
    # net rate per compartment (IN positive, OUT negative)
    if not df.empty:
        df['signed'] = np.where(df['direction'] == 'OUT', -df['rate_m3d'], df['rate_m3d'])
        comp = df.groupby('compartment')['signed'].sum().reset_index()
        comp = comp.rename(columns={'signed': 'net_rate_m3d'}).sort_values('net_rate_m3d')
    else:
        comp = pd.DataFrame(columns=['compartment', 'net_rate_m3d'])
    return df, comp


def _subsample(seq, max_samples):
    """Evenly-spaced subsample of a list (all of it if short enough).

    Reading every one of ~1949 daily budget records is prohibitively slow and a
    time-mean is well approximated by an even subsample, so the .cbc means are
    taken over at most ``max_samples`` stress periods.
    """
    seq = list(seq)
    if max_samples is None or len(seq) <= max_samples:
        return seq
    idx = np.linspace(0, len(seq) - 1, max_samples).round().astype(int)
    return [seq[i] for i in sorted(set(idx))]


def package_budget(sim_ws, cbc_fn, kperkstp_skip=1, max_samples=120):
    """Mean rate [m3/d] of each budget term in a package .cbc file.

    Works for the UZF and SFR budget files (and the GWF cbc). Skips the first
    ``kperkstp_skip`` stress periods (the steady state) and averages over an
    even subsample of at most ``max_samples`` of the rest (set None for all).
    """
    import flopy
    cbc = flopy.utils.CellBudgetFile(os.path.join(sim_ws, cbc_fn))
    records = [r.strip().decode() if isinstance(r, bytes) else str(r).strip()
               for r in cbc.get_unique_record_names()]
    kk = cbc.get_kstpkper()
    keep = _subsample([k for k in kk if k[1] >= kperkstp_skip] or kk, max_samples)
    out = {}
    for rec in records:
        tot = 0.0
        n = 0
        for k in keep:
            try:
                data = cbc.get_data(kstpkper=k, text=rec)
            except Exception:
                continue
            if not data:
                continue
            arr = data[0]
            val = arr['q'].sum() if hasattr(arr, 'dtype') and arr.dtype.names and 'q' in arr.dtype.names \
                else np.asarray(arr, dtype=float).sum()
            tot += float(val)
            n += 1
        if n:
            out[rec] = tot / n
    return out


def layer_storage_change(sim_ws, name, nlay, nrow, ncol, kper_skip=1,
                         max_samples=120):
    """Time-mean storage-change map [m3/d] per layer from the GWF cbc."""
    import flopy
    cbc = flopy.utils.CellBudgetFile(os.path.join(sim_ws, '%s.cbc' % name))
    names = [r.strip().decode() if isinstance(r, bytes) else str(r).strip()
             for r in cbc.get_unique_record_names()]
    sto = [nm for nm in names if 'STO' in nm.upper()]
    kk = _subsample([k for k in cbc.get_kstpkper() if k[1] >= kper_skip]
                    or cbc.get_kstpkper(), max_samples)
    acc = np.zeros((nlay, nrow, ncol))
    n = 0
    for k in kk:
        step = np.zeros((nlay, nrow, ncol))
        got = False
        for nm in sto:
            try:
                d = cbc.get_data(kstpkper=k, text=nm, full3D=True)
            except Exception:
                d = None
            if d:
                step += np.asarray(d[0]).reshape(nlay, nrow, ncol)
                got = True
        if got:
            acc += step
            n += 1
    return acc / max(n, 1)


# --------------------------------------------------------------------- #
# main entry point
# --------------------------------------------------------------------- #

def run_postproc(sim_ws, ds_ws, name='lamatamm', dates=None,
                 xll=739300.0, yll=4553050.0, cs=50.0, verbose=True,
                 out_root=None):
    """Produce the full post-processing figure/CSV set into <out-dir>/_output/.

    Returns the list of files written. Each figure is guarded so a missing
    output file (e.g. no SFR in this run) skips that figure rather than
    aborting the whole report.
    """
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    import flopy

    out = _mkdir(out_root or sim_ws, '_output')
    written = []
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, verbosity_level=0)
    gwf = sim.get_model()
    mg = gwf.modelgrid
    # WP1c.7: a DISV model has no rows or columns. Everything here works on
    # the (ncpl, 1) shape MARMITES already uses for a mesh, so adopt it --
    # the CSV grids come out as one column per layer, and the head maps with
    # real axes are drawn by native_suite's rasterising path instead.
    locate = None
    if getattr(mg, 'nrow', None) is None:
        nlay = int(mg.nlay)
        nrow = int(np.atleast_1d(mg.ncpl)[0])
        ncol = 1
        # observation points cannot be found by grid arithmetic on a mesh;
        # flopy's own intersect() gives the icell2d, which IS the row here
        locate = lambda x, y: (int(mg.intersect(x, y)), 0)   # noqa: E731
    else:
        nlay, nrow, ncol = mg.nlay, mg.nrow, mg.ncol
    top = np.asarray(mg.top, dtype=float).reshape(nrow, ncol)

    hds = flopy.utils.HeadFile(os.path.join(sim_ws, '%s.hds' % name))
    kk = hds.get_kstpkper()
    kk_real = [k for k in kk if k[1] >= 1] or kk        # drop steady state
    # the mean map is taken over an even subsample; the obs series keeps the
    # full record (one cell, cheap) via _obs_series_from_hds
    kk_map = _subsample(kk_real, 200)
    H = np.array([hds.get_data(kstpkper=k) for k in kk_map])   # (nt, nlay, nrow, ncol)
    # a DISV read comes back (nlay, 1, ncpl); the model shape here is
    # (nlay, ncpl, 1), and the middle axis being 1 makes this a pure reshape
    H = H.reshape(H.shape[0], nlay, nrow, ncol)
    H = np.where(np.abs(H) > 1e29, np.nan, H)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)  # all-NaN layers
        mean_h = np.nanmean(H, axis=0)

    if dates is None:
        dates = _dates_from_dataset(ds_ws, len(kk_real))

    # --- obs vs computed heads ---------------------------------------- #
    try:
        f = _fig_obs_heads(hds, ds_ws, kk_real, dates, top, nlay, nrow, ncol,
                           xll, yll, cs, out, locate=locate)
        written += f
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   obs-vs-computed heads skipped: %r' % exc)

    # --- mean head + depth-to-water grids ----------------------------- #
    # CSV only: the figures these used to carry were plain imshow maps with
    # no coordinate frame, and native_suite draws the same two fields as
    # GWmap_head / MMmap_dgwt with the full plotLAYER axes.
    for L in range(nlay):
        for kind, grid in (('mean_head', mean_h[L]),
                           ('mean_depth', top - mean_h[L])):
            fn = os.path.join(out, '%s_L%d.csv' % (kind, L + 1))
            np.savetxt(fn, grid, delimiter=',')
            written += [fn]

    # --- budget by compartment + yearly ------------------------------- #
    try:
        written += _fig_budget(sim_ws, name, dates, out)
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   listing budget skipped: %r' % exc)

    # --- UZF / SFR internal budgets ----------------------------------- #
    for pkg, cbc in (('uzf', '%s.uzf.cbc' % name), ('sfr', '%s.sfr.cbc' % name)):
        path = os.path.join(sim_ws, cbc)
        if not os.path.exists(path):
            continue
        try:
            b = package_budget(sim_ws, cbc)
            if not b:
                continue
            s = pd.Series(b).sort_values()
            fig, ax = plt.subplots(figsize=(7, 0.4 * len(s) + 1.5))
            s.plot.barh(ax=ax, color='tab:green' if pkg == 'uzf' else 'tab:blue')
            ax.axvline(0, color='k', lw=0.6)
            ax.set_xlabel('mean rate [m3/d]')
            ax.set_title('%s internal budget' % pkg.upper())
            fn = os.path.join(out, 'budget_%s.png' % pkg)
            fig.savefig(fn, dpi=140, bbox_inches='tight')
            plt.close(fig)
            s.to_csv(os.path.join(out, 'budget_%s.csv' % pkg),
                     header=['mean_rate_m3d'])
            written += [fn]
        except Exception as exc:           # pragma: no cover
            if verbose:
                print('   %s budget skipped: %r' % (pkg, exc))

    # NOTE: the native plotLAYER head map is produced by native_suite's
    # _native_aquifer_maps(), per layer and from an exact mean over every
    # stress period, into this run's results folder. The older
    # native_layer_maps() that used to be called here wrote into the MODEL
    # workspace and predated the map corrections; it has been removed.

    # --- per-layer storage change ------------------------------------- #
    # likewise CSV only (m3/d, + = gain)
    try:
        sto = layer_storage_change(sim_ws, name, nlay, nrow, ncol)
        for L in range(nlay):
            fn = os.path.join(out, 'storage_change_L%d.csv' % (L + 1))
            np.savetxt(fn, sto[L], delimiter=',')
            written += [fn]
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   storage-change grid skipped: %r' % exc)

    if verbose:
        print('postproc: %d file(s) written to %s' % (len(written), out))
    return written


def run_preproc(sim_ws, ds_ws, name='lamatamm', mf_ws=None, verbose=True,
                out_root=None, cMF=None, ctx=None, res=None, trunk=None,
                gis_ws=None):
    """Input maps into <out-dir>/_input/.

    Every parameter field -- geometry, aquifer properties, UZF soil
    parameters, boundary packages, and the MARMITES soil / meteo / irrigation
    / vegetation zoning -- is drawn by the native ``plotLAYER`` as
    ``IN_<nnn>_<name>``, so all of them carry the MODFLOW index frame on the
    top and right and the projected coordinates on the bottom and left. The
    one figure that is not a parameter field is the site's general map,
    rebuilt from the GIS layers in ``gis_ws``.

    A second, plainer set of ``aq_*`` / ``mm_*`` imshow maps used to be drawn
    here as well. They duplicated the native set field for field, in a style
    borrowed from another project and without any coordinate frame, so they
    have been removed.
    """
    import matplotlib
    matplotlib.use('agg')

    out = _mkdir(out_root or sim_ws, '_input')
    mf_ws = mf_ws or os.path.join(ds_ws, 'MF_ws')
    written = []
    MMplot = _mmplot(trunk)

    # --- the site's general map, from the GIS layers ------------------ #
    try:
        written += _fig_general_map(out, cMF=cMF, gis_ws=gis_ws,
                                    verbose=verbose)
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   general map skipped: %r' % exc)

    # --- the native parameter-field maps ------------------------------ #
    if cMF is not None and ctx is not None and MMplot is not None:
        try:
            written += _native_input_maps(MMplot, out, cMF, ctx, res=res,
                                          verbose=verbose)
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   native input maps skipped: %r' % exc)
    if verbose:
        print('preproc: %d file(s) written to %s' % (len(written), out))
    return written


# Where the site's GIS layers live: in the WORKSPACE, never in the repo. The
# repo holds only what MM and MF read directly, and no shapefile is read by
# either -- this map is the one thing that touches them. So this is only a
# default: the caller can pass gis_ws, MARMITES_GIS_WS overrides it, and the
# figure is skipped when nothing is found. Forward slashes on purpose --
# Windows accepts them and they keep the literal free of escapes.
GIS_WS = os.environ.get('MARMITES_GIS_WS', 'E:/00code_ws/LAMATA_new/GIS')

# Soil_type.SoilType -> (face colour, hatch). chr(92) is a backslash: writing
# the hatch as an escaped literal here is unreadable.
_SOIL_STYLE = (
    ('Alluvium', '#7fc97f', '///'),
    ('Regolith', 'none', None),
    ('Outcrop', '#fdc086', chr(92) * 2),
)


def _fig_general_map(out, cMF=None, gis_ws=None, verbose=True):
    """The site's general map, rebuilt from the ArcMap GIS layers.

    A cartographic figure rather than a model one: soil types, irrigation
    plots, ponds, the hydrographic network, the catchment boundary and the
    monitoring / observation points over shaded relief and elevation
    contours. It follows the published La Mata figure
    (``GIS/LaMata_MM_MF_202109.png``, ArcMap 10.8) as closely as the
    available layers allow -- the orthophoto that figure uses as its
    background is not among them, so a hillshade off the model DEM stands in.

    A layer whose shapefile is absent is simply left out.

    Needs geopandas, which nothing else in this module does. Missing
    geopandas or a missing workspace skips the figure rather than failing.
    """
    try:
        import geopandas as gpd
    except Exception:                                # pragma: no cover
        if verbose:
            print('   general map skipped: geopandas not available')
        return []
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    # No fall-back onto ds_ws/GIS: the dataset folder in the repo is for what
    # MM and MF import, and putting shapefiles there to satisfy this figure is
    # exactly the mix-up to avoid.
    G = gis_ws or GIS_WS
    if not os.path.isdir(G):
        if verbose:
            print('   general map skipped: no GIS workspace at %r' % G)
        return []

    def rd(stem):
        fn = os.path.join(G, stem + '.shp')
        if not os.path.exists(fn):
            return None
        try:
            return gpd.read_file(fn)
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   general map: %s unreadable (%r)' % (stem, exc))
            return None

    soil, irr, pond = rd('Soil_type'), rd('Irr_Fields'), rd('lm_ponds')
    hydro, lim = rd('hydrography'), rd('lm_lim')
    mon, obs = rd('202109MonitPts'), rd('202109ObsPts')
    if all(g is None for g in (soil, irr, pond, hydro, lim, mon, obs)):
        if verbose:
            print('   general map skipped: no layers found in %r' % G)
        return []

    fig, ax = plt.subplots(figsize=(8.27, 8.27))
    handles = []
    K = 1000.0                                       # metres -> km
    ext = None

    def to_km(g):
        """The layer scaled to kilometres, so the axes read in km like every
        other map in the suite."""
        return g.set_geometry(g.geometry.scale(1 / K, 1 / K, origin=(0, 0)))

    # --- elevation: shaded relief plus labelled contours --------------- #
    if cMF is not None and getattr(cMF, 'xllcorner', None) is not None:
        from matplotlib.colors import LightSource
        from marmites_rasterise import MapAdapter
        # WP1c.7: a hillshade needs a raster. On a mesh the DEM is a per-cell
        # vector, and np.gradient on the (ncpl, 1) shape fails outright, so the
        # backdrop is rendered on the display raster.
        DA = MapAdapter(cMF)
        nrow, ncol = DA.nrow, DA.ncol
        dem = DA.lay(np.asarray(cMF.elev, dtype=float).reshape(1, -1))[0]
        if DA.on_mesh:
            dr = float(np.mean(np.diff(DA.dr.x_edge)))
            dc = float(np.mean(-np.diff(DA.dr.y_edge)))
            x0, y0 = float(DA.dr.x_edge[0]), float(DA.dr.y_edge[-1])
            # pixels the mesh does not cover would make the relief NaN
            if np.isnan(dem).any():
                dem = np.where(np.isnan(dem), np.nanmean(dem), dem)
        else:
            dr = float(np.mean(np.asarray(cMF.delr, dtype=float)))
            dc = float(np.mean(np.asarray(cMF.delc, dtype=float)))
            x0, y0 = float(cMF.xllcorner), float(cMF.yllcorner)
        ext = [x0 / K, (x0 + ncol * dr) / K, y0 / K, (y0 + nrow * dc) / K]
        # bilinear: the DEM is on the 50 m model grid, and nearest-neighbour
        # relief reads as a checkerboard at this scale
        ax.imshow(LightSource(azdeg=315, altdeg=45).hillshade(
                      dem, vert_exag=2.0, dx=dr, dy=dc),
                  cmap='Greys_r', extent=ext, origin='upper', alpha=0.4,
                  interpolation='bilinear', zorder=0)
        X = (x0 + (np.arange(ncol) + 0.5) * dr) / K
        Y = (y0 + (nrow - np.arange(nrow) - 0.5) * dc) / K
        lv = np.arange(np.floor(dem.min() / 10) * 10, dem.max() + 10, 10.0)
        CS = ax.contour(X, Y, dem, levels=lv, colors='#8c5a2b',
                        linewidths=0.6, zorder=2)
        ax.clabel(CS, CS.levels[::2], fmt='%d', fontsize=6, colors='#8c5a2b')
        handles.append(Line2D([], [], color='#8c5a2b', lw=0.8,
                              label='Elevation (m)'))

    # --- soil types ---------------------------------------------------- #
    if soil is not None and 'SoilType' in soil:
        for kind, face, hatch in _SOIL_STYLE:
            sel = soil[soil['SoilType'].astype(str) == kind]
            if not len(sel):
                continue
            to_km(sel).plot(ax=ax, facecolor=face, edgecolor='0.4', lw=0.3,
                            hatch=hatch, alpha=0.45, zorder=1)
            handles.append(Patch(facecolor=face, edgecolor='0.4',
                                 hatch=hatch, alpha=0.45, label=kind))

    # --- irrigation plots and ponds ------------------------------------ #
    if irr is not None and len(irr):
        # Irr_Fields holds the two pivots AND a frame polygon spanning the
        # whole sheet; hatched as-is the frame covers the entire map, so
        # anything larger than a fifth of the area is dropped
        if ext is not None:
            cap = 0.2 * (ext[1] - ext[0]) * (ext[3] - ext[2]) * K * K
            irr = irr[irr.geometry.area < cap]
        if len(irr):
            to_km(irr).plot(ax=ax, facecolor='#4292c6', edgecolor='#08519c',
                            lw=0.6, hatch='//', alpha=0.45, zorder=3)
            handles.append(Patch(facecolor='#4292c6', edgecolor='#08519c',
                                 hatch='//', alpha=0.45,
                                 label='Irrigation plot'))
    if pond is not None and len(pond):
        to_km(pond).plot(ax=ax, facecolor='#08306b', edgecolor='#08306b',
                         lw=0.6, zorder=4)
        handles.append(Patch(facecolor='#08306b', label='Ponds'))

    # --- hydrography and the catchment --------------------------------- #
    if hydro is not None and len(hydro):
        to_km(hydro).plot(ax=ax, color='#1f6fd0', lw=1.0, zorder=5)
        handles.append(Line2D([], [], color='#1f6fd0', lw=1.2,
                              label='Hydrography'))
    if lim is not None and len(lim):
        to_km(lim).boundary.plot(ax=ax, color='red', lw=2.0, ls='--', zorder=6)
        handles.append(Line2D([], [], color='red', lw=2.0, ls='--',
                              label='Catchment boundary'))

    # --- monitoring and observation points ----------------------------- #
    for g, colour, label in ((mon, '#33cc33', 'Monitoring points'),
                             (obs, '#33ddee', 'Observation points')):
        if g is None or not len(g):
            continue
        gx = to_km(g)
        ax.plot(gx.geometry.x, gx.geometry.y, 'o', ms=8, mfc=colour, mec='k',
                mew=0.8, ls='none', zorder=7)
        if 'Name' in gx:
            import matplotlib.patheffects as pe
            # P0 / SM / EC sit within ~100 m of each other, so the labels
            # are placed round the marker in turn instead of all up-right;
            # the halo keeps them readable over the relief
            box = ((7, 8, 'left'), (7, -14, 'left'),
                   (-7, 8, 'right'), (-7, -14, 'right'))
            for k, (nm, pt) in enumerate(zip(gx['Name'], gx.geometry)):
                dx, dy, ha = box[k % len(box)]
                ax.annotate(str(nm), xy=(pt.x, pt.y), xytext=(dx, dy),
                            textcoords='offset points', fontsize=8,
                            fontweight='bold', ha=ha, zorder=8,
                            path_effects=[pe.withStroke(linewidth=2,
                                                        foreground='white')])
        handles.append(Line2D([], [], color=colour, marker='o', ms=8, mec='k',
                              ls='none', label=label))

    # --- frame on the model grid, so it matches the other pages -------- #
    if ext is not None:
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
    ax.set_aspect('equal')
    ax.set_xlabel('X [km]', fontsize=10)
    ax.set_ylabel('Y [km]', fontsize=10)
    plt.setp(ax.get_yticklabels(), rotation=90, va='center')
    ax.set_title('La Mata catchment', fontsize=13)
    ax.legend(handles=handles, loc='upper right', fontsize=8, framealpha=0.9)
    crs = None
    for g in (lim, hydro, soil):
        if g is not None and g.crs is not None:
            crs = g.crs.to_string()
            break
    if crs:
        ax.annotate('Coordinate system: %s' % crs, xy=(0.01, 0.01),
                    xycoords='axes fraction', fontsize=7,
                    bbox=dict(fc='white', ec='0.5', lw=0.4))
    fn = os.path.join(out, 'IN_000_general_map.png')
    fig.savefig(fn, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return [fn]


# input fields drawn over the whole grid rather than over the active cells
_IN_NOMASK = ('ibound',)

# blue tones for the observed-vs-computed head figure, taken from the Blues
# ramp that every other water figure uses
_OBS_BLUE_LINE = '#2171b5'
_OBS_BLUE_MARK = '#08306b'

def _mmplot(trunk=None):
    """Import and return the native ``MARMITESplot_v3`` module (or None).

    Three entry points need it -- ``native_suite``, ``run_preproc`` and the
    network overlay -- and each used to carry its own copy of the sys.path
    dance.
    """
    if trunk is None:
        trunk = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    mmplot_dir = os.path.join(trunk, 'MARMITESutilities', 'MARMITESplot')
    if mmplot_dir not in sys.path:
        sys.path.insert(0, mmplot_dir)
    try:
        import MARMITESplot_v3 as MMplot
    except Exception:                                    # pragma: no cover
        return None
    return MMplot


def native_suite(out_dir, cMF, ctx, res, ds_ws=None, trunk=None, verbose=True,
                 sim_ws=None, sankey=True, sankey_full=True, sankey_min_flux=0.05,
                 map_days=6, sankey_obs_years=False):
    """Run the native MARMITESplot figures on the coupled run's in-memory data.

    Called from the runner where ``cMF``/``ctx``/``res`` exist. Produces the
    original MARMITES figures the NWT driver made:

      * plotTIMESERIES_CATCH  -- catchment water-balance time series
      * plotLAYER             -- per-layer maps of the time-mean fluxes
      * plotWBsankey          -- whole-catchment flux Sankey (core + full)

    The Sankey reads the aquifer per-layer fluxes from the MODFLOW cell budgets
    in ``sim_ws`` (defaults to the parent of ``out_dir``). ``sankey_min_flux``
    (mm/y) hides tiny flows on the core diagram; ``sankey_full`` also emits an
    all-flux version. Each figure is guarded so a missing input skips it rather
    than aborting.
    """
    import matplotlib
    matplotlib.use('agg')
    MMplot = _mmplot(trunk)
    nper = res['perc'].shape[0]
    if sim_ws is None:
        sim_ws = os.path.dirname(os.path.abspath(out_dir))
    name = str(getattr(cMF, 'modelname', 'lamatamm')).lower()
    written = []

    # NOTE: the catchment water-balance series (plotTIMESERIES_CATCH, written
    # as native_wb_catchment[_part2].png) is deliberately NOT produced. Its
    # panels came out empty -- the routine wants the legacy driver's combined
    # catchment-flux array, and reassembling that layout here never filled the
    # curves -- and the per-point series below carry the same fluxes.

    # --- water-balance Sankey (catchment core + full, then per-point) - #
    if sankey:
        # ONE pass over the cell budget covering the catchment AND every
        # observation cell. Reading it per target used to cost 9.3 min for 12
        # targets; the result is cached beside the figures as a small digest.
        agg = None
        try:
            targets = [[(c[1], c[2]) for c in ctx.cells]]      # 0 = catchment
            obs_ij = res.get('obs_ij')
            if obs_ij is not None:
                targets += [[(int(i), int(j))] for i, j in np.asarray(obs_ij)]
            agg = _aquifer_pass(sim_ws, name, cMF, targets, nper,
                                cache_dir=sim_ws, verbose=verbose)
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   aquifer pass failed (%r); falling back to per-target '
                      'reads' % exc)
        try:
            written += _native_sankey(MMplot, out_dir, cMF, ctx, res, sim_ws,
                                      name, min_flux=sankey_min_flux,
                                      full_diagram=sankey_full, verbose=verbose,
                                      agg=agg)
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   native Sankey skipped: %r' % exc)
        try:
            written += _native_sankey_obs(MMplot, out_dir, cMF, ctx, res, sim_ws,
                                          name, min_flux=sankey_min_flux,
                                          verbose=verbose, agg=agg,
                                          plot_years=sankey_obs_years)
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   per-point Sankey skipped: %r' % exc)

    # --- per-observation-point soil-column time series (Stage 2) ------- #
    try:
        written += _native_obs_timeseries(MMplot, out_dir, cMF, ctx, res, sim_ws,
                                          name, agg=agg if sankey else None,
                                          ds_ws=ds_ws, verbose=verbose)
    except Exception as exc:                         # pragma: no cover
        if verbose:
            print('   native obs time series skipped: %r' % exc)

    # --- result maps: the MM fluxes and the aquifer terms ------------- #
    try:
        written += _native_result_maps(MMplot, out_dir, cMF, ctx, res, sim_ws,
                                       name, ndays=map_days, verbose=verbose)
    except Exception as exc:                         # pragma: no cover
        if verbose:
            print('   native result maps skipped: %r' % exc)

    if verbose:
        print('native MARMITESplot: %d figure(s) -> %s' % (len(written), out_dir))
    return written


_MAP_FLUXES = (
    ('iP', 'P', 'rainfall', 'Blues'),
    ('iPe', 'Pe', 'effective rainfall', 'Blues'),
    ('iI', 'I', 'infiltration', 'Blues'),
    ('iperc', 'Rp', 'percolation', 'Blues'),
    ('iSsoil_pc', 'theta', 'soil moisture', 'Blues'),
    ('idgwt', 'dgwt', 'depth to water table', 'Blues'),
    ('iRo', 'Ro', 'runoff', 'Reds'),
    ('iEi', 'Ei', 'interception', 'Reds'),
    ('iEow', 'Eow', 'open-water evaporation', 'Reds'),
    ('iETsoil', 'ETsoil', 'soil ET', 'Reds'),
    ('iEg', 'Eg', 'groundwater evaporation', 'Reds'),
    ('iTg', 'Tg', 'groundwater transpiration', 'Reds'),
    ('iETg', 'ETg', 'groundwater ET', 'Reds'),
    ('iEXFg', 'EXFg', 'exfiltration', 'Reds'),
    ('iuzthick', 'uzthick', 'unsaturated thickness', 'YlOrBr'),
    ('idSsoil', 'dSsoil', 'change in soil storage', 'Blues'),
    ('idSsurf', 'dSsurf', 'change in surface storage', 'Blues'),
)


def _vrange(v):
    """(Vmin, Vmax) for plotLAYER, with a NEGLIGIBLE negative snapped to zero.

    plotLAYER switches to the diverging coolwarm_r ramp whenever Vmin < 0 <
    Vmax. A field that is physically one-signed but carries a -1e-9 rounding
    value would then waste half the colour range on a bound that is not real,
    so anything below 1e-6 of the positive range is treated as zero.
    """
    lo, hi = float(np.nanmin(v)), float(np.nanmax(v))
    tol = 1e-6 * max(abs(hi), abs(lo), 1.0)
    if -tol < lo < 0.0:
        lo = 0.0
    if 0.0 < hi < tol:
        hi = 0.0
    return lo, hi


def _obs4map(res):
    """Observation points for plotLAYER's ``points`` overlay: [lbl, i, j, lay],
    the four parallel lists the native routine expects."""
    if 'obs_ij' not in res:
        return None
    ij = np.asarray(res['obs_ij'])
    names = [n.decode() if isinstance(n, bytes) else str(n)
             for n in res.get('obs_names', [str(k) for k in range(len(ij))])]
    return [names, [int(v) for v in ij[:, 0]], [int(v) for v in ij[:, 1]],
            [0] * len(ij)]


def _hydro_year_index(DATE, ini_month):
    """Hydrological-year boundary indices, ported verbatim from the legacy
    driver (startMARMITES_v3.py ~449-483) so the Sankey aggregates exactly as
    the NWT post-processing did. Returns ``(HYindex, year_lst)``.

    ``DATE`` is the array of matplotlib date numbers (one per stress period).
    A run shorter than one hydrological year takes the driver's own short-run
    branch (``HYindex = [h0, h0, N-1, N-1]``), so the 'average' panel still
    renders (scaled to mm/y).
    """
    import matplotlib as mpl
    DATE = np.asarray(DATE, dtype=float)
    year_lst = []
    HYindex = []
    d0 = mpl.dates.num2date(DATE[0])
    if d0.month < ini_month or (d0.month == ini_month and d0.day == 1):
        year_lst.append(d0.year)
    else:
        year_lst.append(d0.year + 1)
    HYindex.append(int(np.argmax(DATE == mpl.dates.datestr2num(
        '%d-%d-01' % (year_lst[0], ini_month)))))
    if np.sum(DATE == mpl.dates.datestr2num(
            '%d-%d-01' % (year_lst[0] + 1, ini_month))) == 0:
        HYindex.append(HYindex[0])
        HYindex.append(len(DATE) - 1)
        HYindex.append(len(DATE) - 1)
        print('   Sankey: the run does not contain a full hydrological year; '
              'the whole-period panel is scaled to mm/y.')
    else:
        y = 0
        while DATE[-1] >= mpl.dates.datestr2num(
                '%d-%d-01' % (year_lst[y] + 1, ini_month)):
            if ini_month == 1:
                iniMonth_prev, add = 12, 0
            else:
                iniMonth_prev, add = ini_month - 1, 1
            indexend = int(np.argmax(DATE == mpl.dates.datestr2num(
                '%d-%d-30' % (year_lst[y] + add, iniMonth_prev))))
            if DATE[-1] >= mpl.dates.datestr2num(
                    '%d-%d-30' % (year_lst[y] + 2, iniMonth_prev)):
                year_lst.append(year_lst[y] + 1)
                HYindex.append(int(np.argmax(DATE == mpl.dates.datestr2num(
                    '%d-%d-01' % (year_lst[y] + 1, ini_month)))))
                indexend = int(np.argmax(DATE == mpl.dates.datestr2num(
                    '%d-%d-30' % (year_lst[y] + 2, iniMonth_prev))))
                y += 1
            else:
                break
        HYindex.append(indexend)
        HYindex.insert(0, 0)
        HYindex.append(len(DATE) - 1)
    return HYindex, year_lst


def _conv_fact(cMF):
    """MODFLOW volumetric-flux -> depth factor, keyed by the length unit
    (m3/d over m2 -> mm/d for meters)."""
    return {1: 304.8, 2: 1000.0, 3: 10.0}.get(int(getattr(cMF, 'lenuni', 2)), 1000.0)


# Budget records to harvest, as (key, cbc text, paknam2).
#
# The paknam2 filter matters: MODFLOW 6 writes BOTH drain packages under the
# single text 'DRN', distinguished only by the package name. On La Mata that is
# the 12-cell boundary drain AND the 1954-cell seepage face (~-2958 m3/d in
# layer 1). Reading text='DRN' alone silently returns whichever record comes
# first, so the seepage was being missed entirely and exfiltration fell back to
# the coupler's own `exf` array.
_AQ_RECORDS = (
    ('UZF-GWRCH', 'UZF-GWRCH', None),
    ('STO-SS', 'STO-SS', None),
    ('STO-SY', 'STO-SY', None),
    ('DRN', 'DRN', 'DRN'),                 # boundary drains
    ('DRN_SEEP', 'DRN', 'DRN_SEEP'),       # seepage face (seep='drn')
    ('WEL', 'WEL', None),                  # groundwater ET sink
)


def _grb_path(sim_ws, name):
    """The MF6 binary grid file, whichever discretisation wrote it (WP1c.7).

    A DIS model writes ``<name>.dis.grb`` and a DISV model ``<name>.disv.grb``.
    Looking only for the DIS name made every mesh run silently lose
    FLOW-JA-FACE: the file was simply absent, `have_flf` went False, and the
    inter-layer flow came out as zeros -- which is what made the Sankey fail
    with "Axis limits cannot be NaN or Inf" rather than draw a wrong number.
    """
    for ext in ('dis', 'disv'):
        cand = os.path.join(sim_ws, '%s.%s.grb' % (name, ext))
        if os.path.exists(cand):
            return cand
    return os.path.join(sim_ws, '%s.dis.grb' % name)


def _ja_down_index(grb_file, nlay, nrow, ncol):
    """Position in the FLOW-JA-FACE array of each cell's DOWNWARD connection.

    Computed once. This replaces ``flopy.mf6.utils.get_structured_faceflows``,
    which (a) re-parses the .grb on EVERY call -- 10.6 ms x nper, by far the
    dominant cost of the old per-stress-period reader -- and (b) whose
    documented ``ia``/``ja`` path is broken upstream (it does
    ``for n in range(grb.nodes)`` while ``grb`` is only bound when ``grb_file``
    is given; still present on flopy master 2026-09-07, PR #1968).

    Returns ``(src, pos, nodes)``: ``flf.ravel()[src] = flowja[pos]``.
    """
    from flopy.mf6.utils import MfGrdFile
    g = MfGrdFile(grb_file, verbose=False)
    ia = np.asarray(g.ia)
    ja = np.asarray(g.ja)
    ncpl = nrow * ncol
    nodes = nlay * ncpl
    src, pos = [], []
    for n in range(nodes - ncpl):
        lo, hi = ia[n], ia[n + 1]
        w = np.where(ja[lo:hi] == n + ncpl)[0]
        if w.size:
            src.append(n)
            pos.append(lo + w[0])
    return np.array(src, int), np.array(pos, int), nodes


def _aquifer_pass(sim_ws, name, cMF, targets, nper, cache_dir=None,
                  verbose=True):
    """Read the cell budget ONCE and reduce every record for EVERY target.

    ``targets`` is a list of ``(i, j)`` cell lists -- one entry per Sankey to be
    drawn (the catchment, then one per observation point). Reducing them all in
    a single sweep is what removes the old O(nper x ntarget) behaviour: the
    reader used to rescan the whole budget for each target, costing 9.3 min for
    12 targets where this takes seconds.

    Returns ``{record: (nper, ntarget, nlay) array in m3/d}`` plus the key
    ``'FLF'`` (flow across each layer's bottom face, + downward, flopy's sign
    convention). Cached as ``<cache_dir>/_aquifer_digest.npz`` -- the caller
    passes the MF6 WORKSPACE, not a results folder, because the digest is keyed
    on the cell budget's size and mtime: it belongs beside the file it derives
    from and is then reused by every later re-plot, whatever its run tag.
    """
    import flopy
    nlay, nrow, ncol = int(cMF.nlay), int(cMF.nrow), int(cMF.ncol)
    cbc_fn = os.path.join(sim_ws, '%s.cbc' % name)
    st = os.stat(cbc_fn)
    sig = '%d_%d_%d_%d' % (st.st_size, int(st.st_mtime), nper, len(targets))

    cache_fn = os.path.join(cache_dir, '_aquifer_digest.npz') if cache_dir else None
    if cache_fn and os.path.exists(cache_fn):
        try:
            z = np.load(cache_fn, allow_pickle=False)
            if str(z['sig']) == sig:
                if verbose:
                    print('   aquifer digest: reusing %s' % os.path.basename(cache_fn))
                return {k: z[k] for k in z.files if k != 'sig'}
        except Exception:
            pass                                   # stale/corrupt -> recompute

    # flat cell indices per target, so the reduction is a gather, not a mask
    flat = []
    for ij in targets:
        a = np.array([int(i) * ncol + int(j) for (i, j) in ij], int)
        flat.append(a)
    ntg = len(targets)

    cbc = flopy.utils.CellBudgetFile(cbc_fn)
    kk = cbc.get_kstpkper()
    off = max(len(kk) - nper, 0)                   # skip the steady SP0
    recs = set(r.strip() for r in cbc.get_unique_record_names(decode=True))
    out = {key: np.zeros((nper, ntg, nlay)) for key, _t, _p in _AQ_RECORDS}
    out['FLF'] = np.zeros((nper, ntg, nlay))

    have_flf = nlay > 1 and 'FLOW-JA-FACE' in recs
    grb = _grb_path(sim_ws, name)
    if have_flf and os.path.exists(grb):
        src, pos, nodes = _ja_down_index(grb, nlay, nrow, ncol)
    else:
        have_flf = False

    def _reduce(a, dest_k):
        """Reduce one stress period's (nlay, ncell) slab onto every target."""
        for t in range(ntg):
            dest_k[t] = a[:, flat[t]].sum(axis=1)

    def _fill(key, text, pak2, dest, flf_mode=False):
        """Fill (nper, ntg, nlay) for one record.

        Preferred path reads the WHOLE series in a single call -- one
        sequential sweep of the file instead of nper seeks. One record type is
        held at a time and freed before the next, so the peak is the largest
        single series (~350 MB for FLOW-JA-FACE on this model). Falls back to
        per-stress-period reads if that does not fit.
        """
        try:
            data = (cbc.get_data(text=text, paknam2=pak2) if flf_mode else
                    cbc.get_data(text=text, paknam2=pak2, full3D=True))
            if not data:
                return 'absent'
            if len(data) < nper + off:
                raise ValueError('%s: got %d records for %d stress periods'
                                 % (key, len(data), nper + off))
            for k in range(nper):
                rec = data[k + off]
                if flf_mode:
                    fj = np.asarray(rec, float).ravel()
                    a = np.zeros(nodes)
                    a[src] = -fj[pos]              # flopy negates; match it
                    a = a.reshape(nlay, -1)
                else:
                    a = np.ma.filled(np.asarray(rec, float), 0.0).reshape(nlay, -1)
                _reduce(a, dest[k])
            del data
            return 'bulk'
        except MemoryError:
            if verbose:
                print('   %s: series too large for one read, '
                      'falling back to per-stress-period' % key)
        for k in range(nper):
            kp = kk[k + off]
            if flf_mode:
                fj = np.asarray(cbc.get_data(text=text, paknam2=pak2,
                                             kstpkper=kp)[0], float).ravel()
                a = np.zeros(nodes)
                a[src] = -fj[pos]
                a = a.reshape(nlay, -1)
            else:
                d = cbc.get_data(text=text, paknam2=pak2, kstpkper=kp, full3D=True)
                if not d:
                    continue
                a = np.ma.filled(np.asarray(d[0], float), 0.0).reshape(nlay, -1)
            _reduce(a, dest[k])
        return 'per-SP'

    t0 = time.time()
    modes = {}
    for key, text, pak2 in _AQ_RECORDS:
        if text in recs:
            modes[key] = _fill(key, text, pak2, out[key])
    if have_flf:
        modes['FLF'] = _fill('FLF', 'FLOW-JA-FACE', None, out['FLF'],
                             flf_mode=True)
    if verbose and modes:
        slow = [r for r, m in modes.items() if m == 'per-SP']
        if slow:
            print('   per-SP fallback used for: %s' % ', '.join(slow))
    if verbose:
        print('   aquifer pass: %d SP x %d target(s) in %.1f s'
              % (nper, ntg, time.time() - t0))

    if cache_fn:
        try:
            os.makedirs(cache_dir, exist_ok=True)
            np.savez_compressed(cache_fn, sig=sig, **out)
        except Exception as exc:                   # pragma: no cover
            if verbose:
                print('   aquifer digest not cached: %r' % exc)
    return out


def _aquifer_layer_fluxes(sim_ws, name, cMF, ctx, res, sel_ij=None,
                          eg_series=None, tg_series=None,
                          agg=None, target=0):
    """Per-layer aquifer fluxes as **mm/d** time series (one value per transient
    stress period), keyed by the names plotWBsankey indexes: ``iRg_L, idSg_L,
    iEXFg_L, iEg_L, iTg_L, iWEL_L, iDRN_L, iFLF_L, iFRF_L, iFFF_L, iGHB_L,
    iCH_L`` for L = 1..nlay, plus ``idSu``.

    ``sel_ij`` selects the cells the balance is taken over: ``None`` = the whole
    catchment (values are the catchment-mean depth, the same basis as
    ``wb_ts``); a list of ``(i, j)`` = just those cells (one obs cell for a
    per-point Sankey). ``eg_series``/``tg_series`` give the Eg/Tg split ratio for
    the WEL sink; when omitted the catchment ``wb_ts`` ratio is used.

    Read straight from the MODFLOW 6 cell budgets (``<name>.cbc``). Sign
    convention matches the native Sankey's aquifer viewpoint:
      * ``Rg`` > 0 into the layer (UZF-GWRCH);
      * ``dSg`` > 0 = release from storage, a source (STO-SS + STO-SY);
      * ``EXFg`` and ``DRN`` < 0 = water leaving the aquifer;
      * ``Eg``/``Tg`` > 0 magnitudes (groundwater ET, drawn via the WEL sink and
        split by the Eg/Tg ratio -- in this MF6 model ETg *is* the WEL package,
        so WEL itself is not drawn: ``iWEL_L`` is left 0);
      * ``FLF[L]`` = flow across the bottom face of layer L, + downward.
    ``FRF``/``FFF`` are 0 (the legacy driver also zeroed the horizontal face
    flows); ``GHB``/``CH`` are 0 (no such packages here).
    """
    import flopy
    nlay = int(cMF.nlay)
    nper = int(res['perc'].shape[0])
    IX = dict(ctx.index)
    wb_ts = np.asarray(res['wb_ts'])
    from marmites_rasterise import model_cell_area

    # cell selection + depth conversion. A boolean (nrow, ncol) mask picks the
    # cells; the depth factor is conv_fact / (their total area) so the flux is a
    # mean depth over the selection (single-cell for an obs point).
    nrow, ncol = int(cMF.nrow), int(cMF.ncol)
    mask = np.zeros((nrow, ncol), bool)
    if sel_ij is None:
        for c in ctx.cells:
            mask[c[1], c[2]] = True
    else:
        for (i, j) in sel_ij:
            mask[int(i), int(j)] = True
    # NOT delc x delr: on a mesh those are the placeholder unit spacings, which
    # would make every cell 1 m2 (WP1c.7).
    area_sel = float(model_cell_area(cMF)[mask].sum())
    to_mm = _conv_fact(cMF) / area_sel if area_sel > 0 else 0.0

    # Volumetric (m3/d) per-layer totals for this target. ``agg`` comes from
    # _aquifer_pass(), which reads the budget ONCE for every target; without it
    # we fall back to a single-target pass so the function still works alone.
    if agg is None:
        agg = _aquifer_pass(sim_ws, name, cMF,
                            [sel_ij if sel_ij is not None
                             else [(c[1], c[2]) for c in ctx.cells]],
                            nper, cache_dir=None, verbose=False)
        target = 0

    def vol(rec):
        """(nper, nlay) m3/d for this target, zeros if the package is absent."""
        a = agg.get(rec)
        return np.zeros((nper, nlay)) if a is None else np.asarray(a)[:, target, :]

    eg = wb_ts[:, IX['iEg']] if eg_series is None else np.asarray(eg_series)
    tg = wb_ts[:, IX['iTg']] if tg_series is None else np.asarray(tg_series)
    egtot = eg + tg
    eg_frac = np.where(egtot > 0, eg / np.where(egtot == 0, 1.0, egtot), 0.5)

    Rg = vol('UZF-GWRCH') * to_mm
    dSg = (vol('STO-SS') + vol('STO-SY')) * to_mm
    DRN = vol('DRN') * to_mm                                  # <0 out
    FLF = vol('FLF') * to_mm
    wel = -vol('WEL') * to_mm                                 # >0 magnitude out
    Egl = wel * eg_frac[:, None]
    Tgl = wel * (1.0 - eg_frac[:, None])
    # exfiltration to the soil: prefer an explicit seepage-drain package,
    # else the coupler's captured exfiltration, assigned to the top layer
    EXF = vol('DRN_SEEP') * to_mm                             # <0 out, or zeros
    if not np.any(EXF) and 'exf' in res:
        exf = np.asarray(res['exf'])
        if sel_ij is None:
            EXF[:, 0] = -exf.mean(axis=1)
        else:
            sel = [_cell_pos(ctx, i, j) for (i, j) in sel_ij]
            sel = [p for p in sel if p is not None]
            if sel:
                EXF[:, 0] = -exf[:, sel].mean(axis=1)

    out = {}
    for L in range(nlay):
        out['iRg_%d' % (L + 1)] = Rg[:, L]
        out['idSg_%d' % (L + 1)] = dSg[:, L]
        out['iEXFg_%d' % (L + 1)] = EXF[:, L]
        out['iEg_%d' % (L + 1)] = Egl[:, L]
        out['iTg_%d' % (L + 1)] = Tgl[:, L]
        out['iWEL_%d' % (L + 1)] = np.zeros(nper)      # ETg drawn as Eg/Tg
        out['iDRN_%d' % (L + 1)] = DRN[:, L]
        out['iFLF_%d' % (L + 1)] = FLF[:, L]
        out['iFRF_%d' % (L + 1)] = np.zeros(nper)
        out['iFFF_%d' % (L + 1)] = np.zeros(nper)
        out['iGHB_%d' % (L + 1)] = np.zeros(nper)
        out['iCH_%d' % (L + 1)] = np.zeros(nper)
    # UZF unsaturated storage change, from the UZF mass balance so it closes:
    #   percolation in = recharge to GW out + dS_unsat
    perc = (wb_ts[:, IX['iperc']] if sel_ij is None
            else np.zeros(nper))              # per-point perc filled by caller
    out['idSu'] = perc - Rg.sum(axis=1)
    # per-layer active-cell and drain-cell counts (over the selection)
    ncell_MM = _active_cells_per_layer(cMF, nlay, mask)
    # a layer has drains for this target if its DRN volume is ever non-zero;
    # that is all the Sankey's 'drncells[L] > 0' guard needs
    drncells = [int(np.any(DRN[:, L] != 0.0)) for L in range(nlay)]
    return out, ncell_MM, drncells


def _cell_pos(ctx, i, j):
    """Position of cell (i, j) in the MM cell list, or None if inactive."""
    for p, c in enumerate(ctx.cells):
        if c[1] == i and c[2] == j:
            return p
    return None


def _active_cells_per_layer(cMF, nlay, mask=None):
    ib = np.abs(np.asarray(cMF.ibound))
    if mask is None:
        return [int((ib[L] != 0).sum()) for L in range(nlay)]
    return [int(((ib[L] != 0) & mask).sum()) for L in range(nlay)]


class _SankeyMF(object):
    """Thin wrapper over cMF supplying the plotting-metadata attributes the
    native Sankey reads, without mutating the real cMF."""
    def __init__(self, cMF, ncell_MM, drncells, dates):
        self._c = cMF
        nlay = int(cMF.nlay)
        self.wel_yn = 1               # groundwater ET drawn (as Eg/Tg)
        self.drn_yn = 1 if any(drncells) else 0
        self.ghb_yn = 0
        self.drncells = drncells
        self.ghbcells = [0] * nlay
        self.ncell_MM = ncell_MM
        self.inputDate = dates
        self.Mnlay = int(getattr(cMF, 'Mnlay', nlay))
        self.Mlay = list(getattr(cMF, 'Mlay', range(1, nlay + 1)))

    def __getattr__(self, name):
        return getattr(self._c, name)


def _assemble_flx(IX, IXS, mmv, mmsv, aq, nper):
    """Build the driver's ``(flx, flxIndex)`` from an MM flux table ``mmv``
    (nper, nidx), a soil table ``mmsv`` (nper, nsl, nidx_s) and the aquifer
    per-layer dict ``aq``. ``mmv`` is the catchment mean (``wb_ts``) or one obs
    cell's series (``mm_obs[:, p]``) -- the layout is identical either way."""
    flx, flxIndex = [], {}

    def put(nm, series):
        flxIndex[nm] = len(flx)
        flx.append(np.asarray(series, dtype=float))

    def mm(key):
        return mmv[:, IX[key]] if key in IX else np.zeros(nper)

    put('iP', mm('iP')); put('iEi', mm('iEi')); put('iPe', mm('iPe'))
    put('idSsurf', mm('idSsurf')); put('iRo', mm('iRo')); put('iEow', mm('iEow'))
    put('idSsoil', mm('idSsoil')); put('iEXFg', mm('iEXFg'))
    put('iI', mm('iI')); put('iSsurf', mm('iSsurf')); put('iperc', mm('iperc'))
    put('iETsoil', mm('iETsoil'))
    put('iEg', mm('iEg')); put('iTg', mm('iTg')); put('iETg', mm('iETg'))
    put('iEsoil', mmsv[:, :, IXS['iEsoil']].sum(axis=1))
    put('iTsoil', mmsv[:, :, IXS['iTsoil']].sum(axis=1))
    put('iExf_1', mmsv[:, 0, IXS['iExf']])                 # top soil layer
    for nm, series in aq.items():
        put(nm, series)
    return flx, flxIndex


def _sankey_dates(cMF, nper):
    """(DATE, HYindex, year_lst) for the Sankey, from cMF.inputDate (sliced to
    nper for a truncated run) or a synthetic daily fallback."""
    import matplotlib as mpl
    dates = getattr(cMF, 'inputDate', None)
    dates = None if dates is None else np.atleast_1d(dates)
    if dates is not None and len(dates) >= nper:
        dates = dates[:nper]
    else:
        import pandas as pd
        dates = mpl.dates.date2num(pd.date_range('2000-01-01', periods=nper, freq='D'))
    DATE = np.asarray(dates, float)
    ini_month = int(getattr(cMF, 'iniMonthHydroYear', 10))
    HYindex, year_lst = _hydro_year_index(DATE, ini_month)
    return DATE, HYindex, year_lst


def _render_sankey(MMplot, out_dir, DATE, flx, flxIndex, HYindex, year_lst,
                   smf, ncell_MM, ibound4Sankey, obspt, fntitle, treshold,
                   verbose, plot_years=True):
    """Call plotWBsankey once and collect the PNGs it wrote for this fntitle.

    "Ignoring fixed x/y limits to fulfill fixed data aspect with adjustable
    data limits" is silenced for the duration: nothing on our side causes it.
    matplotlib's own Sankey.finish() calls ax.axis([...]), which pins the limits
    and turns autoscaling off, and then asks for set_aspect('equal',
    adjustable='datalim') -- so it reports the conflict it just created, once
    per panel per axis (156 lines on a full La Mata run).

    Note it is emitted with _log.warning(), NOT the warnings module, so
    warnings.filterwarnings cannot touch it; only a logging filter can.
    """
    with _quiet_mpl_aspect():
        return _render_sankey_inner(MMplot, out_dir, DATE, flx, flxIndex,
                                    HYindex, year_lst, smf, ncell_MM,
                                    ibound4Sankey, obspt, fntitle, treshold,
                                    verbose, plot_years)


class _DropAspectNoise(logging.Filter):
    """Drop only matplotlib's fixed-aspect/fixed-limits complaint."""

    def filter(self, record):
        return not str(record.getMessage()).startswith('Ignoring fixed ')


@contextlib.contextmanager
def _quiet_mpl_aspect():
    log = logging.getLogger('matplotlib.axes._base')
    filt = _DropAspectNoise()
    log.addFilter(filt)
    try:
        yield
    finally:
        log.removeFilter(filt)


def _render_sankey_inner(MMplot, out_dir, DATE, flx, flxIndex, HYindex,
                         year_lst, smf, ncell_MM, ibound4Sankey, obspt,
                         fntitle, treshold, verbose, plot_years=True):
    written = []
    try:
        MMplot.plotWBsankey(out_dir, DATE, flx, flxIndex,
                            fn='%s_WBsankey' % fntitle, indexTime=HYindex,
                            year_lst=year_lst, cMF=smf, ncell_MM=ncell_MM,
                            obspt=obspt, fntitle=fntitle,
                            ibound4Sankey=ibound4Sankey, treshold=treshold,
                            plot_years=plot_years)
        written = [os.path.join(out_dir, f) for f in os.listdir(out_dir)
                   if f.startswith('_%s_WBsankey' % fntitle) and f.endswith('.png')]
        if verbose:
            print('   native Sankey (%s, treshold=%.3g): %d page(s)'
                  % (fntitle, treshold, len(written)))
    except Exception as exc:                             # pragma: no cover
        if verbose:
            print('   native Sankey (%s) skipped: %r' % (fntitle, exc))
    return written


def _native_sankey(MMplot, out_dir, cMF, ctx, res, sim_ws, name,
                   min_flux=0.05, full_diagram=True, verbose=True, agg=None):
    """Drive the native ``plotWBsankey`` for the whole catchment: MM terms from
    ``wb_ts``/``wb_ts_soil``, aquifer per-layer terms from the cell budgets.
    Renders a decluttered *core* diagram (fluxes below ``min_flux`` hidden) and,
    if ``full_diagram``, a full one (``treshold=0``). Returns the PNGs written."""
    IX = dict(ctx.index); IXS = dict(ctx.index_S)
    wb_ts = np.asarray(res['wb_ts'])
    wb_ts_soil = np.asarray(res['wb_ts_soil'])
    nper = wb_ts.shape[0]
    nlay = int(cMF.nlay)

    aq, ncell_MM, drncells = _aquifer_layer_fluxes(sim_ws, name, cMF, ctx, res,
                                                   agg=agg, target=0)
    flx, flxIndex = _assemble_flx(IX, IXS, wb_ts, wb_ts_soil, aq, nper)
    DATE, HYindex, year_lst = _sankey_dates(cMF, nper)
    smf = _SankeyMF(cMF, ncell_MM, drncells, DATE)
    ibound4Sankey = [1 if ncell_MM[L] > 0 else 0 for L in range(nlay)]

    written = []
    jobs = [('catchment', min_flux)]
    if full_diagram:
        jobs.append(('catchment_full', 0.0))
    for fntitle, tres in jobs:
        written += _render_sankey(MMplot, out_dir, DATE, flx, flxIndex, HYindex,
                                  year_lst, smf, ncell_MM, ibound4Sankey,
                                  'catchment', fntitle, tres, verbose)
    return written


def _native_sankey_obs(MMplot, out_dir, cMF, ctx, res, sim_ws, name,
                       min_flux=0.05, verbose=True, agg=None,
                       plot_years=False):
    """Per-observation-point Sankeys (the legacy ``flxObs_lst`` path).

    Needs the coupler's obs-cell capture: ``res['mm_obs']`` (nper, nobs, nidx),
    ``res['mms_obs']`` (nper, nobs, nsl, nidx_s), ``res['obs_ij']`` and
    ``res['obs_names']``. Builds a per-point ``flx`` (MM terms from that cell's
    captured series, aquifer terms from that single cell's cell-budget fluxes)
    and renders one core diagram per point. Returns the PNGs written."""
    if 'mm_obs' not in res:
        if verbose:
            print('   per-point Sankey skipped: run has no obs capture '
                  '(mm_obs); re-run after obs cells were resolved.')
        return []
    IX = dict(ctx.index); IXS = dict(ctx.index_S)
    mm_obs = np.asarray(res['mm_obs'])           # (nper, nobs, nidx)
    mms_obs = np.asarray(res['mms_obs'])         # (nper, nobs, nsl, nidx_s)
    obs_ij = np.asarray(res['obs_ij'])           # (nobs, 2)
    names = [n.decode() if isinstance(n, bytes) else str(n)
             for n in res.get('obs_names', [str(k) for k in range(len(obs_ij))])]
    nper = mm_obs.shape[0]
    nlay = int(cMF.nlay)
    DATE, HYindex, year_lst = _sankey_dates(cMF, nper)

    written = []
    for p in range(mm_obs.shape[1]):
        i, j = int(obs_ij[p, 0]), int(obs_ij[p, 1])
        mmv = mm_obs[:, p, :]
        mmsv = mms_obs[:, p, :, :]
        aq, ncell_MM, drncells = _aquifer_layer_fluxes(
            sim_ws, name, cMF, ctx, res, sel_ij=[(i, j)],
            eg_series=mmv[:, IX['iEg']], tg_series=mmv[:, IX['iTg']],
            agg=agg, target=p + 1)
        # the point's UZF storage change from its own percolation
        aq['idSu'] = mmv[:, IX['iperc']] - sum(aq['iRg_%d' % (L + 1)]
                                               for L in range(nlay))
        flx, flxIndex = _assemble_flx(IX, IXS, mmv, mmsv, aq, nper)
        smf = _SankeyMF(cMF, ncell_MM, drncells, DATE)
        ibound4Sankey = [1 if ncell_MM[L] > 0 else 0 for L in range(nlay)]
        written += _render_sankey(MMplot, out_dir, DATE, flx, flxIndex, HYindex,
                                  year_lst, smf, ncell_MM, ibound4Sankey,
                                  names[p], 'obs_%s' % names[p], min_flux,
                                  verbose, plot_years=plot_years)
    return written


# TeX labels, in INDEX_MM / INDEX_MM_SOIL order (from startMARMITES_v3.py).
_MM_TEX = [r'$P$', r'$PT$', r'$PE$', r'$Pe$', r'$S_{surf}$', r'$Ro$', r'$Exf_g$',
           r'$E_{ow}$', r'$MB_{soil}$', r'$E_I$', r'$E_o$', r'$E_g$', r'$T_g$',
           r'$\Delta S_{surf}$', r'$ET_g$', r'$ET_{soil}$', r'$\theta$',
           r'$\Delta S_{soil}$', r'$perc$', r'$h\/corr$', r'$d$', r'$thick_p$',
           r'$I$', r'$MB_{surf}$']
_MMS_TEX = ['E_{soil}', 'T_{soil}', '\\theta', 'R_{soil}', 'Exf',
            '\\Delta \\theta', 'S_{soil}', 'SAT', 'MB_{soil}']
_CLR_LST = ['darkgreen', 'firebrick', 'darkmagenta', 'goldenrod', 'green',
            'tomato', 'magenta', 'yellow']


def _is_vertex(hds):
    """True when a head file belongs to a DISV model (WP1c.7).

    flopy's ``get_ts`` wants ``(lay, row, col)`` for DIS and ``(lay, icell2d)``
    for DISV; handing it the wrong arity raises "Row index N out of range
    [0, 1)" -- which is what a mesh run used to do, because under the
    ``(ncpl, 1)`` convention the icell2d arrives as the ROW.
    """
    # MF6 writes a DISV head record with NROW=1 and NCOL=ncpl in its header
    # (verified: DISV -> nrow=1 ncol=989; DIS -> nrow=65 ncol=60).
    return (int(getattr(hds, 'nrow', 0) or 0) == 1
            and int(getattr(hds, 'ncol', 0) or 0) > 1)


def _cellid(vertex, lay, i, j):
    """Address for flopy's ``get_ts``.

    A DISV head record is STORED as (nlay, 1, ncpl), and flopy treats a file
    opened without a modelgrid as structured -- so the right address is the
    3-tuple ``(lay, 0, icell2d)`` in the file's own framing. That needs no
    modelgrid and no simulation load; under the ``(ncpl, 1)`` convention the
    icell2d arrives here as ``i``.
    """
    return (int(lay), 0, int(i)) if vertex else (int(lay), int(i), int(j))


def _obs_head_series(sim_ws, name, nlay, cells_ij, nper):
    """Head time series at each obs cell, per layer: (nobs, nlay, nper).

    Read straight from the MODFLOW 6 head file rather than any HDF5 copy.
    """
    import flopy
    hds = flopy.utils.HeadFile(os.path.join(sim_ws, '%s.hds' % name))
    kk = hds.get_kstpkper()
    off = max(len(kk) - nper, 0)                   # drop the steady SP0
    out = np.full((len(cells_ij), nlay, nper), np.nan)
    vertex = _is_vertex(hds)
    for p, (i, j) in enumerate(cells_ij):
        ts = hds.get_ts([_cellid(vertex, L, i, j) for L in range(nlay)])
        for L in range(nlay):
            out[p, L] = np.asarray(ts[off:off + nper, L + 1], float)
    return out


def _native_obs_timeseries(MMplot, out_dir, cMF, ctx, res, sim_ws, name,
                           agg=None, ds_ws=None, verbose=True):
    """Per-observation-point soil-column time series (native ``plotTIMESERIES``).

    Ports the observation loop of ``startMARMITES_v3.py`` (~1797-2130): the
    driver built, for ONE cell, a flux list holding every MM flux, every soil
    flux both summed and per soil layer, the per-layer heads and depths, the
    SATFLOW head, and the observed head / soil moisture / runoff series. The MM
    side now comes from the coupler's obs capture (``mm_obs``/``mms_obs``), the
    heads from the MF6 ``.hds`` and the recharge from the cell budget, so no
    legacy HDF5 is involved.
    """
    if 'mm_obs' not in res:
        if verbose:
            print('   obs time series skipped: run has no obs capture (mm_obs).')
        return []
    IX = dict(ctx.index)
    IXS = dict(ctx.index_S)
    mm_obs = np.asarray(res['mm_obs'])
    mms_obs = np.asarray(res['mms_obs'])
    obs_ij = np.asarray(res['obs_ij'])
    names = [n.decode() if isinstance(n, bytes) else str(n)
             for n in res.get('obs_names', [str(k) for k in range(len(obs_ij))])]
    nper, nobs = mm_obs.shape[0], mm_obs.shape[1]
    nlay = int(cMF.nlay)
    DATE, HYindex, _yr = _sankey_dates(cMF, nper)
    cMFd = _SankeyMF(cMF, [0] * nlay, [0] * nlay, DATE)   # supplies inputDate

    # observed series + SATFLOW parameters, from the native reader
    obs = {}
    try:
        obs, _ol, _oc, _ocl = cMF.cPROCESS.inputObs(
            inputObs_fn='inputObs.txt', inputObsHEADS_fn='inputObsHEADS',
            inputObsSM_fn='inputObsSM', inputObsRo_fn='inputObsRo',
            inputDate=DATE, _nslmax=int(ctx._nslmax), nlay=nlay)
    except Exception as exc:
        if verbose:
            print('   obs series unavailable (%r); plotting computed only' % exc)

    heads = _obs_head_series(sim_ws, name, nlay, obs_ij, nper)
    from marmites_rasterise import model_cell_area
    cell_area = model_cell_area(cMF)
    cf = _conv_fact(cMF)
    import MARMITESsoil_v3 as _MMsoil
    satflow = _MMsoil.SATFLOW()
    h_lbl = list(getattr(cMF, 'h_lbl', [str(L + 1) for L in range(nlay)]))
    written = []
    stats = []                 # (name, compCalibCritObs tuple) per obs point

    for p in range(nobs):
        i, j = int(obs_ij[p, 0]), int(obs_ij[p, 1])
        o = names[p]
        try:
            zone = int(ctx.gridSOIL[i, j]) - 1
            nsl = int(ctx._nsl[zone])
            l_high = int(cMF.outcropL[i, j]) - 1
            area = float(cell_area[i, j])
            flx, lbl, idx = [], [], {}

            def put(nm, series, tex):
                # keep masked arrays MASKED: the observed series carry the
                # hnoflo sentinel (9999.999) on days without a measurement, and
                # np.asarray() would silently turn those back into real values
                # and wreck the head axis.
                idx[nm] = len(flx)
                flx.append(series if np.ma.isMaskedArray(series)
                           else np.asarray(series, float))
                lbl.append(tex)

            # every MM flux, in index order
            for nm, k in sorted(IX.items(), key=lambda kv: kv[1]):
                put(nm, mm_obs[:, p, k],
                    _MM_TEX[k] if k < len(_MM_TEX) else nm)
            # soil fluxes summed over the soil layers, then per soil layer
            for nm, k in sorted(IXS.items(), key=lambda kv: kv[1]):
                put(nm, mms_obs[:, p, :, k].sum(axis=1), r'$%s$' % _MMS_TEX[k])
            for nm, k in sorted(IXS.items(), key=lambda kv: kv[1]):
                for l in range(nsl):
                    # the layer index goes INSIDE any existing subscript, so
                    # 'E_{soil}' becomes 'E_{soil,1}' and not the invalid
                    # 'E_{soil}_{1}' (same rule as the legacy driver)
                    tex = _MMS_TEX[k]
                    if '}' in tex:
                        tex = r'$%s,%d}$' % (tex[:tex.index('}')], l + 1)
                    else:
                        tex = r'$%s_{%d}$' % (tex, l + 1)
                    put('%s_%d' % (nm, l + 1), mms_obs[:, p, l, k], tex)
            # groundwater ET, attributed to the outcropping layer
            for nm in ('iEg', 'iTg'):
                for L in range(nlay):
                    put('%s_%d' % (nm, L + 1),
                        mm_obs[:, p, IX[nm]] if L == l_high else np.zeros(nper),
                        r'$%s_{g,%d}$' % ('E' if nm == 'iEg' else 'T', L + 1))
            # heads and depth to water, per layer, from the MF6 .hds
            for L in range(nlay):
                put('ih_%s' % h_lbl[L], heads[p, L], r'$h_{%s}$' % h_lbl[L])
                put('id_%s' % h_lbl[L], heads[p, L] - cMF.elev[i, j],
                    r'$d_{%s}$' % h_lbl[L])
            # aquifer terms at this cell, from the cell budget (mm/d)
            to_mm = cf / area
            if agg is not None:
                rg = np.asarray(agg['UZF-GWRCH'])[:, p + 1, :] * to_mm
                sg = ((np.asarray(agg['STO-SS'])[:, p + 1, :]
                       + np.asarray(agg['STO-SY'])[:, p + 1, :]) * to_mm)
            else:
                rg = sg = np.zeros((nper, nlay))
            put('iRg', rg[:, l_high], r'$Rg$')
            for L in range(nlay):
                put('iRg_%d' % (L + 1), rg[:, L], r'$Rg_{%d}$' % (L + 1))
                put('idSg_%d' % (L + 1), sg[:, L], r'$\Delta S_{g,%d}$' % (L + 1))
            put('idSu', mm_obs[:, p, IX['iperc']] - rg.sum(axis=1), r'$\Delta S_u$')
            # SATFLOW head from the same recharge, and the Picard-era corrections
            # (MF6 has no Picard loop, so hcorr is 0 and dcorr is just the depth)
            oo = obs.get(o, {})
            try:
                h_sf = satflow.runSATFLOW(rg[:, l_high], float(oo['hi']),
                                          float(oo['h0']), float(oo['RC']),
                                          float(oo['STO']))
            except Exception:
                h_sf = np.full(nper, np.nan)
            put('ih_SF', h_sf, r'$hSF$')
            put('id_SF', np.asarray(h_sf) - cMF.elev[i, j], r'$dSF$')
            put('idcorr', flx[idx['ihcorr']] - cMF.elev[i, j], r'$d\/corr$')
            # observed series
            oh = oo.get('obs_h')
            if oh is not None:
                v = np.ma.masked_values(np.asarray(oh)[0, :nper], cMF.hnoflo, atol=0.09)
                put('ihobs', v, r'$h\/obs$')
                put('idobs', v - cMF.elev[i, j], r'$d\/obs$')
            osm = oo.get('obs_SM')
            if osm is not None:
                for l in range(min(nsl, len(osm))):
                    put('iSobs_%d' % (l + 1),
                        np.ma.masked_values(np.asarray(osm[l])[:nper],
                                            cMF.hnoflo, atol=0.09),
                        r'$\theta_{%d}\/obs$' % (l + 1))
            oro = oo.get('obs_Ro')
            if oro is not None:
                put('iRoobs', np.ma.masked_values(np.asarray(oro)[0, :nper],
                                                  cMF.hnoflo, atol=0.09),
                    r'$Ro\/obs$')

            # The head panel must span EVERY head-like series it draws, not
            # just the MODFLOW ones: plotTIMESERIES also plots the SATFLOW head
            # and the observed record, and scaling to the MODFLOW heads alone
            # pushed hSF (739.6-749.5 m at P0) clean off the top of the axis.
            hser = [heads[p][np.isfinite(heads[p])], np.asarray(h_sf, float)]
            if 'ihobs' in idx:
                # masked entries are the hnoflo sentinel -> NaN, not values
                hser.append(np.ma.filled(np.ma.asarray(flx[idx['ihobs']],
                                                       dtype=float), np.nan))
            hs = np.concatenate([np.asarray(a, float).ravel() for a in hser])
            hs = hs[np.isfinite(hs)]
            hmax = float(np.nanmax(hs)) if hs.size else float(cMF.elev[i, j])
            hmin = float(np.nanmin(hs)) if hs.size else float(cMF.elev[i, j]) - 10.0
            fn = os.path.join(out_dir, '_0%s_ts.png' % o)
            MMplot.plotTIMESERIES(
                cMFd, i, j, flx, lbl, idx,
                ctx._Sm[zone], ctx._Sr[zone], fn,
                'Time series of fluxes at observation point %s' % o,
                'i = %d, j = %d, l_obs = %d, elev. = %.2f m' % (
                    i + 1, j + 1, l_high + 1, cMF.elev[i, j]),
                _CLR_LST, hmax, hmin, o, int(oo.get('lay', l_high)), nsl,
                float(cMF.elev[i, j]), int(getattr(cMF, 'iniMonthHydroYear', 10)),
                date_ini=DATE[HYindex[1]], date_end=DATE[HYindex[-2]])
            if os.path.exists(fn):
                written.append(fn)
            # groundwater-flux figure at the same point, from the same flux list
            fng = os.path.join(out_dir, '_0%s_tsGW.png' % o)
            try:
                MMplot.plotTIMESERIES_flxGW(
                    cMFd, flx, lbl, idx, fng,
                    'Groundwater fluxes at observation point %s' % o,
                    iniMonthHydroYear=int(getattr(cMF, 'iniMonthHydroYear', 10)),
                    date_ini=DATE[HYindex[1]], date_end=DATE[HYindex[-2]])
                # the native routine writes only the '_part3MF' variant
                stem = os.path.splitext(os.path.basename(fng))[0]
                written += [os.path.join(out_dir, f) for f in os.listdir(out_dir)
                            if f.startswith(stem) and f.endswith('.png')]
            except Exception as exc:                 # pragma: no cover
                if verbose:
                    print('   obs GW time series (%s) skipped: %r' % (o, exc))
            # calibration statistics for this point, over the hydro-year window
            a, b = HYindex[1], HYindex[-2]
            try:
                l_obs = int(oo.get('lay', l_high))
                stats.append((o, cMF.cPROCESS.compCalibCritObs(
                    mms_obs[a:b, p, 0:nsl, IXS['iSsoil_pc_s']],
                    np.ma.masked_values(heads[p, min(l_obs, nlay - 1), a:b],
                                        cMF.hnoflo, atol=0.09),
                    (np.asarray(oo['obs_SM'])[:, a:b]
                     if oo.get('obs_SM') is not None else []),
                    (np.asarray(oo['obs_h'])[0, a:b]
                     if oo.get('obs_h') is not None else None),
                    cMF.hnoflo, o, nsl, mm_obs[a:b, p, IX['ihcorr']])))
            except Exception as exc:                 # pragma: no cover
                if verbose:
                    print('   calib stats (%s) skipped: %r' % (o, exc))
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   obs time series (%s) skipped: %r' % (o, exc))
    if verbose:
        print('   native obs time series: %d figure(s) from %d point(s)'
              % (len(written), nobs))
    written += _native_calibcrit(MMplot, out_dir, cMF, stats, verbose=verbose)
    return written


def _aquifer_map_pass(sim_ws, name, cMF, nper, cache_dir=None, verbose=True):
    """Per-layer, per-CELL time-mean of each aquifer budget record.

    The Sankey pass reduces the budget to a few targets; the maps need the
    spatial field kept, so this is a separate single sweep that accumulates
    (nlay, nrow, ncol) sums. Returns ``{record: (nlay, nrow, ncol)}`` in m3/d,
    cached as ``<cache_dir>/_aquifer_maps.npz`` in the MF6 WORKSPACE rather
    than a results folder: the digest is keyed on the cell budget's size and
    mtime, so keeping it beside that file lets every later re-plot reuse it
    whatever its run tag.
    """
    import flopy
    nlay, nrow, ncol = int(cMF.nlay), int(cMF.nrow), int(cMF.ncol)
    cbc_fn = os.path.join(sim_ws, '%s.cbc' % name)
    st = os.stat(cbc_fn)
    sig = 'map_%d_%d_%d' % (st.st_size, int(st.st_mtime), nper)
    cache_fn = os.path.join(cache_dir, '_aquifer_maps.npz') if cache_dir else None
    if cache_fn and os.path.exists(cache_fn):
        try:
            z = np.load(cache_fn, allow_pickle=False)
            if str(z['sig']) == sig:
                if verbose:
                    print('   aquifer map digest: reusing %s'
                          % os.path.basename(cache_fn))
                return {k: z[k] for k in z.files if k != 'sig'}
        except Exception:
            pass
    cbc = flopy.utils.CellBudgetFile(cbc_fn)
    kk = cbc.get_kstpkper()
    off = max(len(kk) - nper, 0)
    recs = set(r.strip() for r in cbc.get_unique_record_names(decode=True))
    out = {}
    t0 = time.time()
    for key, text, pak2 in _AQ_RECORDS:
        if text not in recs:
            continue
        acc = np.zeros((nlay, nrow, ncol))
        n = 0
        try:
            data = cbc.get_data(text=text, paknam2=pak2, full3D=True)
            if not data or len(data) < nper + off:
                continue
            for k in range(nper):
                acc += np.ma.filled(np.asarray(data[k + off], float),
                                    0.0).reshape(nlay, nrow, ncol)
                n += 1
            del data
        except MemoryError:                          # pragma: no cover
            for k in range(nper):
                d = cbc.get_data(text=text, paknam2=pak2,
                                 kstpkper=kk[k + off], full3D=True)
                if d:
                    acc += np.ma.filled(np.asarray(d[0], float),
                                        0.0).reshape(nlay, nrow, ncol)
                    n += 1
        if n:
            out[key] = acc / n
    # vertical exchange, via the same precomputed JA index the Sankey uses
    grb = _grb_path(sim_ws, name)
    if nlay > 1 and 'FLOW-JA-FACE' in recs and os.path.exists(grb):
        try:
            src, pos, nodes = _ja_down_index(grb, nlay, nrow, ncol)
            data = cbc.get_data(text='FLOW-JA-FACE')
            acc = np.zeros(nodes)
            for k in range(nper):
                fj = np.asarray(data[k + off], float).ravel()
                a = np.zeros(nodes)
                a[src] = -fj[pos]                    # flopy's sign convention
                acc += a
            del data
            out['FLF'] = (acc / nper).reshape(nlay, nrow, ncol)
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   FLF map skipped: %r' % exc)
    if verbose:
        print('   aquifer map pass: %d record(s) over %d SP in %.1f s'
              % (len(out), nper, time.time() - t0))
    if cache_fn and out:
        try:
            os.makedirs(cache_dir, exist_ok=True)
            np.savez_compressed(cache_fn, sig=sig, **out)
        except Exception:                            # pragma: no cover
            pass
    return out


# Per-layer aquifer maps to draw: (record key, stem, colourbar label, sign)
# sign flips the packages MODFLOW reports as negative (water leaving the
# aquifer) so the map shows a positive magnitude.
_AQ_MAPS = (
    ('UZF-GWRCH', 'Rg', 'recharge to groundwater', +1.0),
    ('DRN_SEEP', 'EXFg', 'seepage to the surface', -1.0),
    ('DRN', 'DRN', 'boundary drainage', -1.0),
    ('WEL', 'ETg', 'groundwater ET', -1.0),
    ('FLF', 'FLF', 'flow across the lower face', +1.0),
)


def _native_result_maps(MMplot, out_dir, cMF, ctx, res, sim_ws, name,
                        ndays=0, verbose=True):
    """Every result map: the MM soil-water-balance fluxes (``MMmap_*``) and
    the aquifer terms read from the MF6 output (``GWmap_*``).

    Both families go through ONE ``draw``. They used to sit in two functions
    with two plotLAYER calls, and drifted: the MM maps passed ``nlay=1``,
    which puts plotLAYER in single-column mode and gave them a full-width
    panel on the sheet while every aquifer map got a half-width one. Here the
    MM fields are handed the same grid geometry as the aquifer ones and drawn
    as a single panel (they are surface fluxes, not per-layer), so the two
    sets come out the same size and in the same place on the page.

    Aquifer maps are time means -- the head from the ``.hds``, every budget
    term from the cell budget, one panel per layer, plus (when ``ndays`` is
    set) head maps on that many evenly spaced days. Volumetric budget terms
    are converted to mm/d per cell so they are comparable with the MM fluxes.
    """
    import flopy
    import matplotlib
    from marmites_rasterise import MapAdapter
    # WP1c.7: on a mesh every array below is rasterised onto a display grid,
    # so plotLAYER and the whole native suite are untouched. On the structured
    # grid the adapter is the identity.
    DA = MapAdapter(cMF)
    nlay = int(cMF.nlay)
    nrow, ncol = DA.nrow, DA.ncol
    if verbose and DA.on_mesh:
        print('   %s' % DA.report())
    nper = int(np.asarray(res['wb_ts']).shape[0])
    hnoflo = float(getattr(cMF, 'hnoflo', 9999.999))
    # Legacy convention: heads, recharge and storage in Blues; the terms that
    # take water OUT of the aquifer (exfiltration, drainage, ET) in Reds. The
    # MM fluxes carry their own colour in _MAP_FLUXES.
    cmap_in = matplotlib.colormaps['Blues']
    cmap_out = matplotlib.colormaps['Reds']
    pts = _obs4map(res)
    area = DA.cell_area()                           # (nrow, ncol) m2
    to_mm = _conv_fact(cMF) / area                  # m3/d -> mm/d, per cell
    ib = DA.lay_int(np.abs(np.asarray(cMF.ibound)))
    before = set(os.listdir(out_dir)) if os.path.isdir(out_dir) else set()

    def draw(V, stem, cblbl, unit, m3, nplot=None, cmap=None, days=None,
             dates=None, jd=None, prefix='GWmap'):
        """One plotLAYER page. ``m3`` is the (nlay, nrow, ncol) inactive-cell
        mask and ``nplot`` how many layer panels to draw (default: all)."""
        V = np.asarray(V, dtype=float)
        m = np.repeat(m3[None, :, :, :], V.shape[0], axis=0)
        VV = np.where(m, hnoflo, V)
        vals = V[~m]
        vals = vals[np.isfinite(vals)]
        if not vals.size:
            return
        lo, hi = _vrange(vals)
        try:
            MMplot.plotLAYER(
                days=days if days is not None else [0],
                str_per=days if days is not None else [0],
                Date=dates if dates is not None else 'NA',
                JD=jd if jd is not None else 'NA',
                ncol=ncol, nrow=nrow, nlay=nlay,
                nplot=nlay if nplot is None else nplot, V=VV,
                cmap=cmap if cmap is not None else cmap_in,
                CBlabel='%s [%s]' % (cblbl, unit), msg='',
                plt_title='%s_%s' % (prefix, stem), MM_ws=out_dir,
                interval_type='linspace', interval_num=5,
                Vmax=[hi], Vmin=[lo], fmt='%5.2f', points=pts, mask=m3,
                hnoflo=hnoflo, cMF=cMF)
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   %s map %s skipped: %r' % (prefix, stem, exc))

    # ---------------- MM soil-water-balance fluxes -------------------- #
    # Time-mean per-cell fluxes, following the legacy driver's output-map
    # conventions (startMARMITES_v3.py ~808-910): linspace/5, the model's own
    # hnoflo, the observation points overlaid.
    IX = dict(ctx.index)
    wb_map = np.asarray(res['wb_map'])              # (ncell, nidx)
    for key, stem, cblbl, cmname in _MAP_FLUXES:
        if key not in IX:
            continue
        g = DA.cells(wb_map[:, IX[key]], ctx.cells, nodata=hnoflo)
        unit = '-' if key in ('iSsoil_pc',) else (
            'm' if key in ('iuzthick', 'idgwt') else 'mm/d')
        # one surface layer, given the aquifer grid's shape so the page comes
        # out the same size, and drawn as the single panel it is
        mm = np.repeat(np.isclose(g, hnoflo, atol=0.09)[None, :, :], nlay,
                       axis=0)
        draw(np.repeat(g[None, None, :, :], nlay, axis=1), stem, cblbl, unit,
             mm, nplot=1, cmap=matplotlib.colormaps[cmname], prefix='MMmap')

    # ---------------- aquifer terms from the MF6 output --------------- #
    m3 = np.zeros((nlay, nrow, ncol), bool)
    for L in range(nlay):
        m3[L] = (ib[L] == 0)

    # heads: exact time mean over every stress period, then a time selection
    hds = flopy.utils.HeadFile(os.path.join(sim_ws, '%s.hds' % name))
    kk = hds.get_kstpkper()
    off = max(len(kk) - nper, 0)
    acc = None
    for k in range(nper):
        h = np.asarray(hds.get_data(kstpkper=kk[k + off]), dtype=float)
        h = np.where(np.abs(h) > 1e29, np.nan, h)
        acc = h if acc is None else acc + h
    draw(DA.lay(np.asarray(acc).reshape(nlay, -1) / nper)[None, :, :, :],
         'head', 'mean head', 'm', m3)
    if ndays and nper > 1:
        sel = np.unique(np.linspace(0, nper - 1, int(ndays)).astype(int))
        V = np.array([DA.lay(np.where(
            np.abs(np.asarray(hds.get_data(kstpkper=kk[s + off]), dtype=float))
            > 1e29, np.nan,
            np.asarray(hds.get_data(kstpkper=kk[s + off]), dtype=float)
        ).reshape(nlay, -1)) for s in sel])
        DATE, _hy, _yr = _sankey_dates(cMF, nper)
        import matplotlib as _mpl
        # JD must be a LIST here: plotLAYER indexes it per panel, so the
        # 'NA' placeholder used for single maps raises IndexError past i=1
        jd = [int(_mpl.dates.num2date(DATE[s]).timetuple().tm_yday) for s in sel]
        draw(V, 'head_series', 'head', 'm', m3,
             days=[int(s) for s in sel],
             dates=[float(DATE[s]) for s in sel], jd=jd)

    # budget terms, per layer, as mm/d
    maps = _aquifer_map_pass(sim_ws, name, cMF, nper, cache_dir=sim_ws,
                             verbose=verbose)
    for key, stem, cblbl, sgn in _AQ_MAPS:
        if key not in maps:
            continue
        V = (sgn * DA.lay(np.asarray(maps[key]).reshape(nlay, -1))
             * to_mm[None, :, :])[None, :, :, :]
        draw(V, stem, cblbl, 'mm/d', m3,
             cmap=cmap_out if sgn < 0 else cmap_in)
    # Effective and net recharge, exactly as the legacy driver derived them
    # (startMARMITES_v3.py ~2433-2462):
    #     Re = Rg + EXF          gross recharge net of exfiltration
    #     Rn = Re + WEL          ... and net of groundwater ET
    # The cbc already reports EXF and WEL negative (water leaving the aquifer),
    # so these are plain sums of the RAW records -- not the sign-flipped
    # magnitudes used for their individual maps. Both are signed fields, so
    # plotLAYER draws them with the remapped coolwarm_r ramp.
    if 'UZF-GWRCH' in maps:
        re = np.asarray(maps['UZF-GWRCH'], float).copy()
        if 'DRN_SEEP' in maps:
            re += np.asarray(maps['DRN_SEEP'], float)
        re = DA.lay(re.reshape(nlay, -1))
        draw((re * to_mm[None, :, :])[None, :, :, :],
             'Re', 'effective recharge (Rg + Exf)', 'mm/d', m3)
        if 'WEL' in maps:
            rn = re + DA.lay(np.asarray(maps['WEL'], float).reshape(nlay, -1))
            draw((rn * to_mm[None, :, :])[None, :, :, :],
                 'Rn', 'net recharge (Rg + Exf + ETg)', 'mm/d', m3)
    # storage change is the sum of the two storage records
    if 'STO-SS' in maps or 'STO-SY' in maps:
        s = (np.asarray(maps.get('STO-SS', 0.0))
             + np.asarray(maps.get('STO-SY', 0.0)))
        s = DA.lay(np.asarray(s, dtype=float).reshape(nlay, -1))
        draw((s * to_mm[None, :, :])[None, :, :, :],
             'dSg', 'release from groundwater storage', 'mm/d', m3)

    after = set(os.listdir(out_dir)) if os.path.isdir(out_dir) else set()
    got = [os.path.join(out_dir, f) for f in sorted(after - before)
           if f.endswith('.png')]
    if verbose:
        print('   native result maps: %d page(s)' % len(got))
    return got


def _native_input_maps(MMplot, out_dir, cMF, ctx, res=None, verbose=True):
    """The native INPUT maps, ported from startMARMITES_v3.py (~678-917).

    The legacy driver mapped every parameter field it had been given before
    running anything -- geometry, aquifer properties, the UZF soil parameters
    and the boundary packages -- as ``IN_<nnn>_<name>``. None of them were
    being produced. Conventions are the legacy ones: gist_rainbow_r (which is
    the INPUT-map colormap, unlike the results), 5 linspace intervals, the
    model's own hnoflo, and per-field number formats. elev / top / botm share
    one elevation scale so they can be read against each other.
    """
    import matplotlib
    from marmites_rasterise import MapAdapter
    # WP1c.7: normalise on the MODEL shape, draw on the DISPLAY shape.
    DA = MapAdapter(cMF)
    nlay = int(cMF.nlay)
    nrow_m, ncol_m = int(cMF.nrow), int(cMF.ncol)     # model
    nrow, ncol = DA.nrow, DA.ncol                     # display
    hnoflo = float(getattr(cMF, 'hnoflo', 9999.999))
    ib = DA.lay_int(np.abs(np.asarray(cMF.ibound)))
    mask = (ib == 0)                                  # (nlay, nrow, ncol)
    mask_all = mask.all(axis=0)                       # cells inactive everywhere
    cmap = matplotlib.colormaps['gist_rainbow_r']
    pts = _obs4map(res) if res is not None else None

    def arr(x):
        """Any parameter field as (nlay, nrow, ncol).

        The ini gives these in whatever form is shortest: a single scalar
        (uniform everywhere), one value per layer, one grid shared by all
        layers, or a full per-layer grid. The legacy driver leaned on
        cPROCESS.float2array plus a try/except; normalising here is explicit
        and covers all four, which is what the eps / thts / thti / thtr / vka
        fields need -- they are scalars on La Mata and were being dropped.
        """
        a = np.asarray(x, dtype=float)
        if a.ndim == 3:
            return DA.lay(a)
        if a.ndim == 2 and a.shape == (nrow_m, ncol_m):
            return DA.lay(np.repeat(a[None, :, :], nlay, axis=0))
        if a.size == 1:
            return np.full((nlay, nrow, ncol), float(a.ravel()[0]))
        if a.size == nlay:
            return np.repeat(a.reshape(nlay, 1, 1), nrow, axis=1).repeat(ncol, axis=2)
        raise ValueError('cannot map a field of shape %s onto (%d, %d, %d)'
                         % (a.shape, nlay, nrow_m, ncol_m))

    # transmissivity, as the legacy derived it (hk x thickness)
    thick = arr(cMF.thick)
    T = arr(cMF.hk_actual) * thick
    # aquifer top per layer: land surface for layer 1, the bottom above for the rest
    # arr() returns (nlay, nrow, ncol), so take the layer we want from each
    top_tmp = np.zeros((nlay, nrow, ncol))
    top_tmp[0] = arr(cMF.top)[0]
    _botm = arr(cMF.botm)
    for L in range(1, nlay):
        top_tmp[L] = _botm[L - 1]

    lst = [('elev', 'Elev.', arr(cMF.elev)),
           ('top', 'Aq. top - $top$', top_tmp),
           ('botm', 'Aq. bot. - $botm$', arr(cMF.botm)),
           ('thick', 'Aq. thick.', thick),
           ('strt', 'Init. heads - $strt$', arr(cMF.strt)),
           ('gridSOILthick', 'Soil thick.', arr(ctx.gridSOILthick)),
           ('gridSsurfhmax', 'Max. stream heigth', arr(ctx.gridSsurfhmax)),
           ('gridSsurfw', 'Stream width', arr(ctx.gridSsurfw)),
           ('hk', 'Horizontal hydraulic cond. - $hk$', arr(cMF.hk_actual)),
           ('Ss', 'Specific storage - $S_s$', arr(cMF.ss_actual)),
           ('Sy', 'Specific yield - $S_y$', arr(cMF.sy_actual)),
           ('vka', 'Vertical hydraulic cond. - $vka$', arr(cMF.vka_actual))]
    if T is not None:
        lst.insert(9, ('T', 'Transmissivity - $T$', T))
    if int(getattr(cMF, 'drn_yn', 0)) == 1:
        lst += [('drn_cond', 'Drain cond.', arr(cMF.drn_cond_array)),
                ('drn_elev', 'Drain elev.', arr(cMF.drn_elev_array))]
    if int(getattr(cMF, 'ghb_yn', 0)) == 1 and hasattr(cMF, 'ghb_cond_array'):
        lst += [('ghb_cond', 'GHB cond.', arr(cMF.ghb_cond_array)),
                ('ghb_head', 'GHB head', arr(cMF.ghb_head_array))]
    if int(getattr(cMF, 'uzf_yn', 0)) == 1:
        lst += [('eps', 'Epsilon - $eps$', arr(cMF.eps_actual)),
                ('thts', 'Sat. water content - $thts$', arr(cMF.thts_actual))]
        if int(getattr(cMF, 'iuzfopt', 0)) == 1:
            lst += [('vks', 'Sat. vert. hydraulic cond. - $vks$',
                     arr(cMF.vks_actual))]
        lst += [('thti', 'Initial water content - $thti$', arr(cMF.thti_actual)),
                ('thtr', 'Residual water content - $thtr$', arr(cMF.thtr_actual))]

    # the model footprint and the MARMITES zoning. These were only ever drawn
    # by the plain imshow maps that used to sit alongside this set; folding
    # them in here is what let those be removed.
    # `ib` came from DA.lay_int and is ALREADY display-shaped, so it must
    # not go through arr() again -- that would rasterise a picture.
    lst += [('ibound', 'Active cells - $ibound$', np.asarray(ib, dtype=float))]
    for attr, stem, cblbl in (('gridSOIL', 'SOILzones', 'Soil zone'),
                              ('gridMETEO', 'METEOzones', 'Meteo. zone'),
                              ('gridIRR', 'IRRzones', 'Irrigation zone')):
        g = getattr(ctx, attr, None)
        if g is not None:
            lst += [(stem, cblbl, arr(g))]
    veg = getattr(ctx, 'gridVEGarea', None)
    if veg is not None:
        veg = np.asarray(veg, dtype=float)
        for v in range(veg.shape[0]):
            lst += [('VEG%darea' % (v + 1), 'Veg. %d area [%%]' % (v + 1),
                     arr(veg[v]))]

    # elev / top / botm share one scale, so the three read against each other
    elev_all = np.concatenate([arr(cMF.elev).ravel(), top_tmp.ravel(),
                               arr(cMF.botm).ravel()])
    elev_all = elev_all[~np.isclose(elev_all, hnoflo, atol=0.09)]
    elev_lo, elev_hi = float(elev_all.min()), float(elev_all.max())

    before = set(os.listdir(out_dir)) if os.path.isdir(out_dir) else set()
    for n, (stem, cblbl, a) in enumerate(lst, start=1):
        try:
            V = a.reshape(1, nlay, nrow, ncol)
            # a field that does not vary by layer gets ONE panel, as the legacy
            # did, and is masked only where every layer is inactive
            if nlay > 1 and np.allclose(a[0], a[1:], equal_nan=True):
                m = np.repeat(mask_all[None, :, :], nlay, axis=0)
                nplot = 1
            else:
                m, nplot = mask, nlay
            if stem in _IN_NOMASK:
                # the footprint itself: masking it by the footprint would
                # leave a field that is 1 everywhere it is drawn
                m = np.zeros(m.shape, dtype=bool)
            good = V[0][~m]
            good = good[~np.isclose(good, hnoflo, atol=0.09)]
            good = good[np.isfinite(good)]
            if not good.size:
                continue
            if stem in ('Ss', 'hk', 'T', 'drn_cond', 'ghb_cond', 'vks'):
                fmt = '%5.e'
            elif stem in ('ibound', 'SOILzones', 'METEOzones', 'IRRzones'):
                fmt = '%5.0f'
            elif stem in ('Sy', 'thts', 'thti', 'thtr', 'gridSOILthick',
                          'gridSsurfhmax', 'gridSsurfw'):
                fmt = '%5.3f'
            else:
                fmt = '%5.1f'
            if stem in ('elev', 'top', 'botm'):
                lo, hi = elev_lo, elev_hi
            else:
                lo, hi = _vrange(good)
            if not hi > lo:
                # a uniform field (an unused zoning grid, say) would give a
                # zero-width colour scale and a degenerate BoundaryNorm
                hi = lo + 1.0
            # an integer field gets one colourbar tick per class, not five
            # linspace ones (1..3 in five steps prints "2" three times)
            nint = 5
            if fmt == '%5.0f':
                nint = max(2, min(5, int(round(hi - lo)) + 1))
            MMplot.plotLAYER(
                days=[0], str_per=[0], Date='NA', JD='NA', ncol=ncol, nrow=nrow,
                nlay=nlay, nplot=nplot, V=np.where(m, hnoflo, V), cmap=cmap,
                CBlabel=cblbl, msg='', plt_title='IN_%03d_%s' % (n, stem),
                MM_ws=out_dir, interval_type='linspace', interval_num=nint,
                Vmax=[hi], Vmin=[lo], fmt=fmt, points=pts, mask=m,
                hnoflo=hnoflo, cMF=cMF)
        except Exception as exc:                       # pragma: no cover
            if verbose:
                print('   input map %s skipped: %r' % (stem, exc))
    after = set(os.listdir(out_dir)) if os.path.isdir(out_dir) else set()
    got = [os.path.join(out_dir, f) for f in sorted(after - before)
           if f.endswith('.png')]
    if verbose:
        print('   native input maps: %d page(s) from %d field(s)'
              % (len(got), len(lst)))
    return got


def _native_calibcrit(MMplot, out_dir, cMF, stats, verbose=True):
    """Calibration-criteria figures (native ``plotCALIBCRIT``).

    ``stats`` is the list of ``(obs name, compCalibCritObs(...))`` tuples
    gathered per observation point. That native routine returns, in order:
    ``rmseHEADS, rmseHEADSc, rmseSM, rsrHEADS, rsrHEADSc, rsrSM,
    nseHEADS, nseHEADSc, nseSM, rHEADS, rHEADSc, rSM``. One figure is written
    per criterion (RMSE, RSR, NSE, r), each comparing soil moisture and heads
    across the points, exactly as the legacy driver did (~2211-2228).
    """
    if not stats:
        return []
    # unpack into the per-criterion lists plotCALIBCRIT expects, keeping each
    # observation name aligned with the points that actually produced a value
    sm = {k: [] for k in ('rmse', 'rsr', 'nse', 'r')}
    hd = {k: [] for k in ('rmse', 'rsr', 'nse', 'r')}
    hc = {k: [] for k in ('rmse', 'rsr', 'nse', 'r')}
    o_sm, o_hd, o_hc = [], [], []
    for o, t in stats:
        (rmseH, rmseHc, rmseS, rsrH, rsrHc, rsrS,
         nseH, nseHc, nseS, rH, rHc, rS) = t
        if rmseH is not None and rmseH != []:
            for k, v in zip(('rmse', 'rsr', 'nse', 'r'), (rmseH, rsrH, nseH, rH)):
                hd[k].append(v)
            o_hd.append(o)
        if rmseHc is not None and rmseHc != []:
            for k, v in zip(('rmse', 'rsr', 'nse', 'r'), (rmseHc, rsrHc, nseHc, rHc)):
                hc[k].append(v)
            o_hc.append(o)
        if rmseS is not None and rmseS != []:
            for k, v in zip(('rmse', 'rsr', 'nse', 'r'), (rmseS, rsrS, nseS, rS)):
                sm[k].append(v)
            o_sm.append(o)

    smmax = max([max(v) for v in sm['rmse']] or [None]) if sm['rmse'] else None
    hdmax = max([max(v) if hasattr(v, '__len__') else v
                 for v in hd['rmse']] or [None]) if hd['rmse'] else None
    jobs = [
        ('rmse', 'RMSE', 'Root mean square error', smmax, hdmax, 0,
         ['(m)', '(%wc)']),
        ('rsr', 'RSR', 'Root mean square error - observations standard '
                       'deviation ratio', None, None, 0, ['', '']),
        ('nse', 'NSE', 'Nash-Sutcliffe efficiency', 1.0, 1.0, None, ['', '']),
        ('r', 'r', "Pearson's correlation coefficient", 1.0, 1.0, -1.0,
         ['', '']),
    ]
    written = []
    for key, crit, title, smx, hmx, ymin, units in jobs:
        fn = os.path.join(out_dir, '__plt_calibcrit%s.png' % crit)
        try:
            MMplot.plotCALIBCRIT(
                calibcritSM=sm[key], calibcritSMobslst=o_sm,
                calibcritHEADS=hd[key], calibcritHEADSobslst=o_hd,
                calibcritHEADSc=hc[key], calibcritHEADScobslst=o_hc,
                plt_export_fn=fn,
                plt_title='Calibration criteria between simulated and observed '
                          'state variables\n%s' % title,
                calibcrit=crit, calibcritSMmax=smx, calibcritHEADSmax=hmx,
                ymin=ymin, units=units, hnoflo=cMF.hnoflo)
            # plotCALIBCRIT appends a page index, so collect by stem
            stem = os.path.splitext(os.path.basename(fn))[0]
            written += [os.path.join(out_dir, f) for f in os.listdir(out_dir)
                        if f.startswith(stem) and f.endswith('.png')]
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   calib criterion %s skipped: %r' % (crit, exc))
    if verbose:
        print('   native calibration criteria: %d figure(s) '
              '(heads at %d point(s), soil moisture at %d)'
              % (len(written), len(o_hd), len(o_sm)))
    return written


def _resolve_obs_cells_mesh(cMF, ctx, ds_ws, verbose=True):
    """WP1c.6: observation points -> mesh cells, by point-in-polygon.

    Reads the points with ``obs_points``, which applies the same '#'-disabled
    convention as the native reader, and keeps their file order so the obs
    capture is comparable across grids.
    """
    proj = cMF.mesh_proj
    pos = {int(c[1]): p for p, c in enumerate(ctx.cells)}
    obs_idx, obs_names, outside = [], [], []
    for pt in obs_points(ds_ws):
        ic = proj.cell_containing(pt['x'], pt['y'])
        p = pos.get(int(ic))
        if p is None:
            outside.append((pt['name'], ic))
            continue
        obs_idx.append(p)
        obs_names.append(pt['name'])
    if verbose:
        for nm, ic in outside:
            print('   obs point %s falls in mesh cell %d, which is not an '
                  'active MM cell; skipped.' % (nm, ic))
        print('resolve_obs_cells (mesh): %d observation cell(s) captured: %s'
              % (len(obs_idx), ', '.join(obs_names)))
    # Two points in one cell is not an error -- a coarse mesh can genuinely
    # merge nearby piezometers -- but it means their modelled series are
    # identical, which would otherwise look like a suspicious coincidence in
    # the calibration plots.
    if verbose and len(set(obs_idx)) < len(obs_idx):
        seen = {}
        for nm, p in zip(obs_names, obs_idx):
            seen.setdefault(p, []).append(nm)
        for p, nms in seen.items():
            if len(nms) > 1:
                print('   NOTE: %s share mesh cell %d, so their modelled '
                      'series are identical.' % (' and '.join(nms),
                                                 ctx.cells[p][1]))
    return obs_idx, obs_names


def resolve_obs_cells(cMF, ctx, ds_ws, verbose=True):
    """Map the enabled observation points (inputObs.txt) to MM cell-list
    positions. Returns ``(obs_idx, obs_names)`` for the coupler's obs capture.

    Reuses the native ``cPROCESS.inputObs`` reader (so the same points, grid
    mapping and disabled-line handling as the legacy driver), then finds each
    point's ``(i, j)`` in the MM cell list. Points whose cell is inactive are
    dropped with a warning.

    On a MESH (WP1c.6) the points are resolved GEOMETRICALLY from their
    coordinates instead. ``cPROCESS.inputObs`` derives ``(i, j)`` from
    ``nrow``/``ncol``/``cellsizeMF``, and under the ``(ncpl, 1)`` convention
    that grid is one cell wide -- so it would reject almost every point as
    lying outside the model before this function ever saw it."""
    if getattr(cMF, 'mesh_proj', None) is not None:
        return _resolve_obs_cells_mesh(cMF, ctx, ds_ws, verbose=verbose)
    try:
        obs, obs_list, _oc, _ocl = cMF.cPROCESS.inputObs(
            inputObs_fn='inputObs.txt', inputObsHEADS_fn='inputObsHEADS',
            inputObsSM_fn='inputObsSM', inputObsRo_fn='inputObsRo',
            inputDate=cMF.inputDate, _nslmax=int(ctx._nslmax), nlay=int(cMF.nlay))
    except Exception as exc:
        if verbose:
            print('resolve_obs_cells: inputObs read failed (%r); no obs capture.'
                  % exc)
        return [], []
    pos = {(c[1], c[2]): p for p, c in enumerate(ctx.cells)}
    obs_idx, obs_names = [], []
    for nm in obs_list:
        o = obs[nm]
        p = pos.get((int(o['i']), int(o['j'])))
        if p is None:
            if verbose:
                print('   obs point %s at (i=%d, j=%d) is not an active MM '
                      'cell; skipped.' % (nm, o['i'], o['j']))
            continue
        obs_idx.append(p)
        obs_names.append(nm)
    if verbose:
        print('resolve_obs_cells: %d observation cell(s) captured: %s'
              % (len(obs_idx), ', '.join(obs_names)))
    return obs_idx, obs_names


def _dates_from_dataset(ds_ws, n):
    import pandas as pd
    fn = os.path.join(ds_ws, 'inputDATE.txt')
    if not os.path.exists(fn):
        return pd.RangeIndex(n)
    ds = []
    with open(fn) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            ds.append(s.split(',')[0].strip())
    d = pd.to_datetime(ds, errors='coerce').dropna()
    if len(d) >= n:
        return d[:n]
    return pd.RangeIndex(n)


def _fig_obs_heads(hds, ds_ws, kk_real, dates, top, nlay, nrow, ncol,
                   xll, yll, cs, out, locate=None):
    import matplotlib.pyplot as plt
    import pandas as pd
    pts = obs_points(ds_ws)
    series = {p['name']: obs_series(ds_ws, p['name']) for p in pts}
    # EVERY monitoring point, not only the ones with an observed series: on
    # La Mata just 4 of the 11 active points in inputObs.txt have an
    # inputObsHEADS_* file, and the computed head is worth seeing at all of
    # them. A point without observations simply gets no markers.
    have = list(pts)
    if not have:
        return []
    dates = pd.to_datetime(dates)
    ncols = 2
    nrows_f = int(np.ceil(len(have) / ncols))
    fig, axes = plt.subplots(nrows_f, ncols, figsize=(12, 3 * nrows_f),
                             squeeze=False)
    tidy = []
    drawn = 0
    for ax, p in zip(axes.ravel(), have):
        try:
            i, j = (locate(p['x'], p['y']) if locate is not None
                    else _xy_to_ij(p['x'], p['y'], xll, yll, cs, nrow, ncol))
            L = min(max(p['lay'] - 1, 0), nlay - 1)
            # full-resolution series at this single cell (cheap via get_ts)
            ts = hds.get_ts(_cellid(_is_vertex(hds), L, i, j))
        except Exception:                  # pragma: no cover - point off grid
            ax.axis('off')
            continue
        comp = ts[:, 1]
        comp = np.where(np.abs(comp) > 1e29, np.nan, comp)
        # get_ts returns every saved step (incl. steady); align to transient
        if comp.shape[0] == len(kk_real) + 1:
            comp = comp[1:]
        n = min(len(dates), comp.shape[0])
        dd, comp = dates[:n], comp[:n]
        # blue tones throughout, as everywhere else water is plotted: the
        # computed head mid-blue, the observations the darkest blue of the
        # same ramp so they read as the same quantity
        ax.plot(dd, comp, '-', color=_OBS_BLUE_LINE, lw=1.0,
                label='computed L%d' % (L + 1))
        obs = series.get(p['name'])
        if obs is not None:
            ax.plot(obs['date'], obs['head'], 'o', ms=4, color=_OBS_BLUE_MARK,
                    mec=_OBS_BLUE_MARK, label='observed')
        ax.axhline(float(top[i, j]), color='0.6', ls=':', lw=0.8,
                   label='land surface')
        ax.set_title('%s  (%d,%d)' % (p['name'], i, j))
        ax.set_ylabel('head [m]')
        ax.legend(fontsize=7, loc='best')
        drawn += 1
        for d, v in zip(dd, comp):
            tidy.append((p['name'], d, float(v)))
    if not drawn:
        plt.close(fig)
        return []
    for ax in axes.ravel()[len(have):]:
        ax.axis('off')
    fig.tight_layout()
    fn = os.path.join(out, 'obs_heads.png')
    fig.savefig(fn, dpi=140, bbox_inches='tight')
    plt.close(fig)
    pd.DataFrame(tidy, columns=['point', 'date', 'computed_head_m']).to_csv(
        os.path.join(out, 'obs_heads_computed.csv'), index=False)
    return [fn]


def _fig_budget(sim_ws, name, dates, out):
    import matplotlib.pyplot as plt
    import pandas as pd
    df, comp = budget_by_compartment(sim_ws, name)
    written = []
    if not comp.empty:
        fig, ax = plt.subplots(figsize=(7, 0.5 * len(comp) + 1.5))
        colors = [_COMP_COLORS.get(c, '0.7') for c in comp['compartment']]
        ax.barh(comp['compartment'], comp['net_rate_m3d'], color=colors)
        ax.axvline(0, color='k', lw=0.6)
        ax.set_xlabel('net rate [m3/d]  (+ into aquifer)')
        ax.set_title('water budget by compartment (post spin-up mean)')
        fn = os.path.join(out, 'budget_compartment.png')
        fig.savefig(fn, dpi=140, bbox_inches='tight')
        plt.close(fig)
        comp.to_csv(os.path.join(out, 'budget_compartment.csv'), index=False)
        df.to_csv(os.path.join(out, 'budget_terms.csv'), index=False)
        written += [fn]
    return written
