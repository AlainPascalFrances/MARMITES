# -*- coding: utf-8 -*-
"""Water-budget figures for the coupled MARMITES-MODFLOW 6 run, with an
optional comparison against the previous MARMITES-MODFLOW-NWT (Picard) run.

Produces, in <out-dir>/figures_nwt_comparison/:

  01_wb_timeseries.png     daily catchment-mean fluxes (new vs reference)
  02_wb_cumulative.png     cumulative volumes over the simulation
  03_wb_totals.png         budget components as mm/yr, side-by-side bars
  04_wb_hydroyear.png      per-hydrological-year totals
  05_map_<flux>.png        time-mean maps: new | reference | difference
  06_heads.png             head time series and mean head map
  07_coupling.png          outer iterations, rejected infiltration

The comparison is a PLAUSIBILITY check, not a match: the MF6 run uses UZF6
(EPSILON clamped to 3.5 vs 2.0 in UZF1) and per-stress-period API coupling
instead of the whole-run Picard loop, so differences are expected and are
exactly what these figures are meant to show.

Usage:
    python tests/plot_water_budget.py
    python tests/plot_water_budget.py --no-reference        (new run only)
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))

# MODFLOW-NWT reference for the comparison panels. It lives in the legacy
# PhD/paper archive, NOT in the repository: it is 1.3 GB and cannot be
# regenerated (the NWT build path was removed from the code in Phase 1).
sys.path.insert(0, TRUNK)
import mm_paths as _mmp                              # noqa: E402
NWT_REF = str(_mmp.NWT_REF)   # WP1.4: ONE definition, in code/mm_paths.py
                              # ($MM_NWT_REF, or the legacy $MARMITES_NWT_REF)
# The dataset (its rasters give the 50 m grid the NWT reference ran on):
# panel 0's, mm_paths.dataset_dir -- the run passes its own to make_figures.
# It was this file's <repo>/example/LaMata, which the run no longer reads,
# and the figures were skipped once that folder was gone (2026-10-09).
DS = str(_mmp.dataset_dir('LaMata'))

import matplotlib  # noqa: E402
matplotlib.use('agg')
import matplotlib.pyplot as plt  # noqa: E402
import h5py  # noqa: E402
from marmites_indices import INDEX_MM  # noqa: E402
from marmites_coupler import MF6Coupler  # noqa: E402
sys.path.insert(0, os.path.join(TRUNK, 'MARMITESutilities', 'MARMITESplot'))
import MARMITESplot_v3 as MMplot  # noqa: E402

HNOFLO = 9999.999

# budget components to compare: label -> MM index key
FLUXES = [
    ('P',        'iP',       'Precipitation'),
    ('Pe',       'iPe',      'Effective precip.'),
    ('Ei',       'iEi',      'Interception'),
    ('Ro',       'iRo',      'Runoff'),
    ('I',        'iI',       'Infiltration'),
    ('Eow',      'iEow',     'Open-water evap.'),
    ('ETsoil',   'iETsoil',  'Soil ET'),
    ('Eg',       'iEg',      'Groundwater evap.'),
    ('Tg',       'iTg',      'Groundwater transp.'),
    ('ETg',      'iETg',     'Groundwater ET'),
    ('Rp',       'iperc',    'Percolation'),
    ('EXFg',     'iEXFg',    'GW exfiltration'),
    ('dSsoil',   'idSsoil',  'd(soil storage)'),
    ('dSsurf',   'idSsurf',  'd(surface storage)'),
]

# The streams and the ponds exchange with the aquifer DIRECTLY, so their
# water never passes through the soil column's Ro or EXFg (user,
# 2026-10-10): their own rows in the table and the figures, MF6 only, and
# the two totals that compare like with like with the NWT run, which routed
# no stream -- its valley groundwater exfiltrated into the soil column and
# partly ran off, so its Ro and EXFg hold what here enters the streams.
SW_TERMS = [
    ('GW->SFR', 'gw_to_sfr', 'groundwater into the streams'),
    ('SFR->GW', 'sfr_to_gw', 'stream losses to the aquifer'),
    ('GW->LAK', 'gw_to_lak', 'groundwater into the ponds'),
    ('LAK->GW', 'lak_to_gw', 'pond losses to the aquifer'),
    ('Qout', 'q_outlet', 'streamflow at the outlet'),
]
SW_TOTALS = [
    ('Ro+netSFR', 'Ro', 'Ro + GW->SFR - SFR->GW, against NWT Ro'),
    ('EXFg+GWsw', 'EXFg', 'EXFg + GW->SFR + GW->LAK, against NWT EXFg'),
]


def load_new(ws):
    fn = os.path.join(ws, MF6Coupler.RESULTS_H5)
    if not os.path.exists(fn):
        sys.exit('coupled results not found: %s\n(run tests/run_lamata_mf6.py first)' % fn)
    with h5py.File(fn, 'r') as h:
        d = {k: h[k][:] for k in ('heads', 'perc', 'etg', 'exf', 'rejinf',
                                  'outer_iters', 'cell_ij')}
        for k in ('wb_ts', 'wb_map', 'wb_ts_soil', 'wb_map_soil',
                  'cell_area', 'grid_shape'):
            d[k] = h[k][:] if k in h else None
        # what outer_iters counts: MF6's own step exposes only its LINEAR
        # solves (the coupler's res['iters_kind']); older runs: outer
        _kind = h['iters_kind'][()] if 'iters_kind' in h else b'outer'
        d['iters_kind'] = (_kind.decode() if isinstance(_kind, bytes)
                           else str(_kind))
    # written by the runner and never read: every map then fell back to
    # "grid size inferred from active cells" (7 warnings a run)
    if d['grid_shape'] is None:
        del d['grid_shape']
    # WP5: with the CRR cascade a cell's runoff includes the run-on of the
    # cells above it, counted again at every cell it crosses (2026-10-04:
    # "Ro 125.4" against 33 that left the soil). NET of the run-on, Ro is
    # what left the soil surface for good -- the NWT reference's meaning,
    # which routed nothing -- and the same as the Sankey's.
    ro, ron = INDEX_MM['iRo'], INDEX_MM['iRunon']
    for k in ('wb_ts', 'wb_map'):
        v = d.get(k)
        if v is not None and v.shape[-1] > ron:
            v[..., ro] = v[..., ro] - v[..., ron]
    d['_ws'] = ws
    return d


def load_reference(chunk=120, ndays=None):
    """Catchment-mean time series and time-mean maps from the NWT run.

    Read in day chunks: the reference MM array is ~370 MB. ``ndays`` limits
    it to the FIRST ndays -- the days the new run simulated: both start on
    the dataset's first day (2008-05-31), and a run with run.nsp = N covers
    days 0..N-1. Without it a 60-day summer run was compared with the NWT
    run's whole 1949-day record: rain 318 against 451 mm/yr, the same rain.
    """
    fn = NWT_REF
    if not os.path.exists(fn):
        print('plot_water_budget: NWT reference not found at\n  %s\n'
              '  -> comparison panels omitted. Set $MARMITES_NWT_REF to the '
              'legacy _h5_MM.h5 to restore them.' % fn)
        return None
    with h5py.File(fn, 'r') as h:
        MM = h['MM']
        nday, nrow, ncol, nidx = MM.shape
        if ndays is not None:
            nday = min(nday, int(ndays))
        ts = np.zeros((nday, nidx))
        acc = np.zeros((nrow, ncol, nidx))
        mask = None
        for a in range(0, nday, chunk):
            b = min(a + chunk, nday)
            blk = np.asarray(MM[a:b], dtype=np.float64)
            if mask is None:                       # active cells: not hnoflo
                mask = np.abs(blk[0, :, :, INDEX_MM['iP']] - HNOFLO) > 0.09
            m = mask[None, :, :, None]
            blk = np.where(m, blk, np.nan)
            ts[a:b] = np.nanmean(blk, axis=(1, 2))
            acc += np.nansum(blk, axis=0)
        acc /= nday
    return {'ts': ts, 'map': acc, 'mask': mask}


# Results root for this run, set by make_figures(). Figures belong with the
# run's results (<ws-root>/out_<stamp>_<tag>/figures), not in the model
# workspace and never in the repository.
_OUT_ROOT = None


def _fig(ws, name):
    d = os.path.join(_OUT_ROOT or ws, 'figures_nwt_comparison')
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)


def grid_shape(new, ref):
    """True (nrow, ncol) of the model grid.

    It cannot be inferred from the active cells alone: La Mata's outermost
    rows/columns are inactive, so max(i)+1, max(j)+1 gives 64x58 instead of
    65x60. Preference order: explicit dataset written by the runner, then the
    reference grid, then the (unreliable) max-index fallback.
    """
    gs = new.get('grid_shape')
    if gs is not None:
        return int(np.ravel(gs)[0]), int(np.ravel(gs)[1])
    if ref is not None:
        return ref['map'].shape[0], ref['map'].shape[1]
    ij = new['cell_ij']
    print('WARNING: grid size inferred from active cells; inactive border '
          'rows/columns will be missing from the maps.')
    return int(ij[:, 0].max()) + 1, int(ij[:, 1].max()) + 1


def _on_mesh(new):
    """A DISV run: its cells are (icell2d, 0) of an ncpl x 1 proxy grid."""
    gs = new.get('grid_shape')
    return gs is not None and int(np.ravel(gs)[1]) == 1         and int(np.ravel(gs)[0]) > 1


def _dataset_grid():
    """(xll, yll, nrow, ncol, cellsize) of the dataset's 50 m grid -- the
    grid the NWT reference was run on."""
    sys.path.insert(0, os.path.join(TRUNK, 'ppMF6'))
    import marmites_props as props
    rect, _names, _others = props.dataset_grid(DS)
    if rect is None:
        raise RuntimeError('no dataset raster declares the grid in %s' % DS)
    return rect


def _mesh_to_grid(new):
    """``(nrow, ncol, to_grid)`` putting a mesh run on the dataset grid.

    THE LINK TO THE NWT REFERENCE. On a mesh the maps were scattered by
    (icell2d, 0) onto an ncpl x 1 'grid' -- blank panels -- and the driver
    switched the reference off altogether, catchment series included. Each
    50 m cell now takes the AREA-weighted mean of the mesh cells overlapping
    it (exact polygon overlay; the polygons from MF6's own .grb), so the
    comparison is like-for-like with the 65 x 60 NWT run. Cached on ``new``.
    """
    if '_to_grid' in new:
        return new['_to_grid']
    import glob
    import flopy
    import shapely
    sys.path.insert(0, os.path.join(TRUNK, 'ppMF6'))
    import marmites_overlay as ov
    grb = sorted(glob.glob(os.path.join(new['_ws'], '*.disv.grb')))
    if not grb:
        raise RuntimeError('no .disv.grb in %s' % new['_ws'])
    mg = flopy.mf6.utils.MfGrdFile(grb[0], verbose=False).modelgrid
    ncol_m = int(np.ravel(new['grid_shape'])[1])
    icell = [int(i) * ncol_m + int(j) for i, j in new['cell_ij']]
    polys = ov.Polygons(np.array([shapely.Polygon(mg.get_cell_vertices(c))
                                  for c in icell], dtype=object), {})
    xll, yll, nrow, ncol, cs = _dataset_grid()
    cells = ov.structured_cells(xll, yll, [cs] * ncol, [cs] * nrow)
    ci, pi, area = ov.overlay(cells, polys)
    den = np.bincount(ci, weights=area, minlength=nrow * ncol)

    def to_grid(values):
        v = np.asarray(values, dtype=float)[pi]
        ok = np.isfinite(v)
        num = np.bincount(ci[ok], weights=v[ok] * area[ok],
                          minlength=nrow * ncol)
        wok = np.bincount(ci[ok], weights=area[ok], minlength=nrow * ncol)
        with np.errstate(invalid='ignore', divide='ignore'):
            g = np.where(wok > 0, num / wok, np.nan)
        return g.reshape(nrow, ncol)

    new['_to_grid'] = (nrow, ncol, to_grid, den.reshape(nrow, ncol))
    return new['_to_grid']


def _scatter(new, ref, values):
    """Per-cell values onto the full (nrow, ncol) grid, NaN elsewhere --
    on a mesh, the dataset grid through an area-weighted overlay."""
    if _on_mesh(new):
        return _mesh_to_grid(new)[2](values)
    nrow, ncol = grid_shape(new, ref)
    ij = new['cell_ij']
    g = np.full((nrow, ncol), np.nan)
    g[ij[:, 0], ij[:, 1]] = values
    return g


# --------------------------------------------------------------------- #

def _series(new, ref, k):
    """(MF6, NWT) per-day series [mm/d] of a FLUXES label, a stream / pond
    term (SW_TERMS, MF6 only: NWT is None) or a like-for-like total
    (SW_TOTALS, against the NWT term it compares with)."""
    fl = dict((a, b) for a, b, _ in FLUXES)
    if k in fl:
        idx = INDEX_MM[fl[k]]
        return (new['wb_ts'][:, idx],
                ref['ts'][:, idx] if ref is not None else None)
    sw = new.get('sw') or {}
    n = new['wb_ts'].shape[0]

    def z(key):
        v = sw.get(key)
        return np.zeros(n) if v is None else np.nan_to_num(np.asarray(v, float))

    if k == 'Ro+netSFR':
        a, b = _series(new, ref, 'Ro')
        return a + z('gw_to_sfr') - z('sfr_to_gw'), b
    if k == 'EXFg+GWsw':
        a, b = _series(new, ref, 'EXFg')
        return a + z('gw_to_sfr') + z('gw_to_lak'), b
    return z(dict((a, b) for a, b, _ in SW_TERMS)[k]), None


# what each like-for-like total stands against in the NWT run
_NWT_OF = dict((a, b) for a, b, _ in SW_TOTALS)


def plot_timeseries(new, ref, ws, dates=None):
    keys = ['P', 'Ro', 'ETsoil', 'ETg', 'Rp']
    sw = bool(new.get('sw'))
    if sw:
        keys += ['EXFg', 'streams']
    fig, axes = plt.subplots(len(keys), 1, figsize=(12, 11 + 4.4 * sw),
                             sharex=True)
    x = np.arange(new['wb_ts'].shape[0])
    for ax, k in zip(axes, keys):
        if k == 'streams':
            # the streams' own terms, MF6 only: NWT routed no stream
            for lab, col in (('GW->SFR', 'tab:blue'), ('SFR->GW', 'tab:orange'),
                             ('Qout', 'k')):
                ax.plot(x, _series(new, ref, lab)[0], lw=0.7, color=col,
                        label='%s (MF6)' % lab)
            ax.legend(fontsize=7, ncol=3, loc='upper right')
            ax.set_ylabel('streams\n(mm/d)', fontsize=9)
        else:
            a, b = _series(new, ref, k)
            ax.plot(x, a, lw=0.7, color='tab:blue', label='MF6 (API)')
            if b is not None:
                ax.plot(np.arange(len(b)), b, lw=0.7, color='tab:red',
                        alpha=0.75, label='NWT (Picard)')
            # like with like: what NWT counted in Ro and EXFg is here partly
            # in the streams
            tot = {'Ro': 'Ro+netSFR', 'EXFg': 'EXFg+GWsw'}.get(k) if sw else None
            if tot:
                ax.plot(x, _series(new, ref, tot)[0], lw=0.7,
                        color='tab:green', label='MF6 %s' % tot)
                ax.legend(fontsize=7, ncol=3, loc='upper right')
            ax.set_ylabel('%s\n(mm/d)' % k, fontsize=9)
        ax.grid(alpha=0.3)
        ax.tick_params(labelsize=8)
    axes[0].legend(fontsize=9, ncol=2)
    axes[-1].set_xlabel('day of simulation')
    fig.suptitle('La Mata - catchment-mean daily water-budget fluxes', fontsize=12)
    fig.tight_layout()
    fig.savefig(_fig(ws, '01_wb_timeseries.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)


def plot_cumulative(new, ref, ws):
    keys = ['P', 'Ro', 'ETsoil', 'ETg', 'Rp']
    if new.get('sw'):
        keys += ['Ro+netSFR', 'EXFg+GWsw', 'Qout']
    fig, ax = plt.subplots(figsize=(11, 6))
    cmap = plt.get_cmap('tab10')
    for c, k in enumerate(keys):
        a, b = _series(new, ref, k)
        ax.plot(np.cumsum(a), color=cmap(c), lw=1.6, label='%s (MF6)' % k)
        if b is not None:
            ax.plot(np.cumsum(b), color=cmap(c), lw=1.2, ls='--', alpha=0.8,
                    label='%s (NWT)' % _NWT_OF.get(k, k))
    ax.set_xlabel('day of simulation')
    ax.set_ylabel('cumulative depth (mm)')
    ax.set_title('Cumulative water-budget components (solid: MF6/API, dashed: NWT/Picard)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(_fig(ws, '02_wb_cumulative.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)


def plot_totals(new, ref, ws):
    """The totals as bars. Returns the soil-column rows only (labels, MF6,
    NWT) for the table: write_summary adds the streams' rows itself."""
    labels, newv, refv = [], [], []
    nday_new = new['wb_ts'].shape[0]
    for short, key, _long in FLUXES:
        idx = INDEX_MM[key]
        labels.append(short)
        newv.append(new['wb_ts'][:, idx].sum() / nday_new * 365.25)
        if ref is not None:
            refv.append(ref['ts'][:, idx].sum() / ref['ts'].shape[0] * 365.25)
    # the streams and the ponds, then the like-for-like totals, after them
    blab, bnew, bref = list(labels), list(newv), list(refv)
    if new.get('sw'):
        for k in [a for a, _b, _c in SW_TERMS] + [a for a, _b, _c in SW_TOTALS]:
            a, b = _series(new, ref, k)
            blab.append(k)
            bnew.append(float(np.sum(a)) / nday_new * 365.25)
            if ref is not None:
                bref.append(np.nan if b is None else
                            float(np.sum(b)) / len(b) * 365.25)
    x = np.arange(len(blab))
    w = 0.38 if ref is not None else 0.6
    fig, ax = plt.subplots(figsize=(12 + 0.5 * (len(blab) - len(labels)), 6))
    ax.bar(x - (w / 2 if ref is not None else 0), bnew, w, label='MF6 (API)',
           color='tab:blue')
    if ref is not None:
        ax.bar(x + w / 2, bref, w, label='NWT (Picard)', color='tab:red',
               alpha=0.85)
    if len(blab) > len(labels):
        sep = len(labels) - 0.5
        ax.axvline(sep, color='k', lw=0.8, ls='--')
        ax.text(sep + 0.15, 0.97, 'streams and ponds (MF6 only), then the '
                'totals that compare like with like with NWT Ro and EXFg',
                transform=ax.get_xaxis_transform(), fontsize=8, va='top')
    ax.set_xticks(x); ax.set_xticklabels(blab, rotation=45, ha='right')
    ax.set_ylabel('mm/year')
    ax.set_title('La Mata water budget - annual-equivalent rates')
    ax.grid(alpha=0.3, axis='y')
    ax.legend(loc='upper left')
    for i, v in enumerate(bnew):
        ax.annotate('%.0f' % v, (x[i] - (w / 2 if ref is not None else 0), v),
                    ha='center', va='bottom' if v >= 0 else 'top', fontsize=7)
    fig.tight_layout()
    fig.savefig(_fig(ws, '03_wb_totals.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)
    return labels, newv, (refv if ref is not None else None)


def plot_hydroyears(new, ref, ws, start_doy=273, days_per_year=365):
    """Totals per hydrological year (starting 1 October by default)."""
    keys = ['P', 'Ro', 'ETsoil', 'ETg', 'Rp']
    sw = bool(new.get('sw'))
    if sw:
        keys += ['Ro+netSFR', 'EXFg+GWsw', 'Qout']
    n = new['wb_ts'].shape[0]
    edges = list(range(0, n, days_per_year))
    fig, ax = plt.subplots(figsize=(11 + 2.0 * sw, 6))
    width = 0.8 / len(keys)
    xs = np.arange(len(edges))
    for c, k in enumerate(keys):
        a, b = _series(new, ref, k)
        tot = [a[s:min(s + days_per_year, n)].sum() for s in edges]
        ax.bar(xs + c * width, tot, width, label='%s (MF6)' % k)
        if b is not None:
            tr = [b[s:min(s + days_per_year, n)].sum() for s in edges]
            ax.plot(xs + c * width, tr, 'k_', ms=14, mew=1.5)
    ax.set_xticks(xs + 0.4); ax.set_xticklabels(['yr %d' % (i + 1) for i in xs])
    ax.set_ylabel('mm/year')
    ax.set_title('Per-year budget totals (bars: MF6/API; black dashes: '
                 'NWT/Picard%s)' % ('; Ro+netSFR and EXFg+GWsw against NWT Ro '
                                    'and EXFg' if sw else ''))
    ax.grid(alpha=0.3, axis='y'); ax.legend(fontsize=8, ncol=3)
    fig.tight_layout()
    fig.savefig(_fig(ws, '04_wb_hydroyear.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)


def plot_maps(new, ref, ws, which=('Rp', 'ETg', 'Eg', 'Tg', 'Ro', 'ETsoil')):
    for short in which:
        key = dict((a, b) for a, b, _ in FLUXES)[short]
        idx = INDEX_MM[key]
        gnew = _scatter(new, ref, new['wb_map'][:, idx])
        panels = [('MF6 (API)', gnew, None)]
        if ref is not None:
            gref = np.where(ref['mask'], ref['map'][:, :, idx], np.nan)
            panels.append(('NWT (Picard)', gref, None))
            panels.append(('difference (MF6 - NWT)', gnew - gref, 'diff'))
        fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 5))
        axes = np.atleast_1d(axes)
        vmax = np.nanmax(np.abs([np.nanmax(p[1]) for p in panels[:2]]))
        for ax, (title, arr, kind) in zip(axes, panels):
            if kind == 'diff':
                lim = np.nanmax(np.abs(arr)) or 1.0
                im = ax.imshow(arr, cmap='RdBu_r', vmin=-lim, vmax=lim)
            else:
                im = ax.imshow(arr, cmap='viridis', vmin=0, vmax=vmax)
            ax.set_title('%s\n%s' % (short, title), fontsize=10)
            ax.set_xticks([]); ax.set_yticks([])
            fig.colorbar(im, ax=ax, shrink=0.8, label='mm/d')
        fig.tight_layout()
        fig.savefig(_fig(ws, '05_map_%s.png' % short), dpi=200,
                    bbox_inches='tight')          # as every other map
        plt.close(fig)


def plot_heads(new, ws, ref=None):
    """Delegates to the recovered MARMITESplot module (single source of
    truth for MARMITES figures)."""
    if _on_mesh(new):
        # on the dataset grid, like the other maps; each 50 m cell weighted
        # by the mesh area it holds
        nrow, ncol, to_grid, den = _mesh_to_grid(new)
        H = np.array([to_grid(h) for h in np.asarray(new['heads'])])
        ii, jj = np.nonzero(den > 0)
        MMplot.plotHEADS(H[:, ii, jj], np.column_stack([ii, jj]),
                         (nrow, ncol), _fig(ws, '06_heads.png'),
                         area=den[ii, jj])
        return
    MMplot.plotHEADS(new['heads'], new['cell_ij'], grid_shape(new, ref),
                     _fig(ws, '06_heads.png'), area=new.get('cell_area'))


def plot_coupling(new, ws):
    """Delegates to the recovered MARMITESplot module."""
    MMplot.plotCOUPLING(new['outer_iters'], rejinf=new['rejinf'],
                        exf=new.get('exf'), plt_export_fn=_fig(ws, '07_coupling.png'),
                        area=new.get('cell_area'))


def _gwf_name(ws):
    """The GWF model's name, from mfsim.nam (a workspace can hold old
    listings of other models)."""
    try:
        for line in open(os.path.join(ws, 'mfsim.nam')):
            w = line.split()
            if len(w) >= 3 and w[0].upper().startswith('GWF6'):
                return w[2].lower()
    except OSError:
        pass
    return None


def surface_water_series(ws, area_m2, nper=None):
    """Per stress period, mm/d over the catchment: the streams' and the
    ponds' exchange with the aquifer (from the MF6 listing's cumulative
    volumes, one record per period -- as the runner's aquifer balance) and
    the outlet streamflow (the SFR observations). The steady first period is
    left out; ``nper`` keeps the last nper periods, to line up with the
    coupler's wb_ts. None for a run without streams or ponds, or without its
    files.
    """
    name = _gwf_name(ws)
    lst = os.path.join(ws, '%s.lst' % name) if name else None
    if not lst or not os.path.exists(lst) or not area_m2:
        return None
    try:
        import flopy
        sys.path.insert(0, os.path.join(TRUNK, 'ppMF6'))
        import marmites_postprocess as PP
        lb = flopy.utils.Mf6ListBudget(lst)
        cum = lb.get_dataframes(diff=False)[1]
        t = np.asarray(lb.get_times(), dtype=float)
        kp = np.array([int(k[1]) for k in lb.get_kstpkper()])
        steady = PP.steady_first(ws, name)
    except Exception:                                # pragma: no cover
        return None
    # the cumulative volume at the END of each period (its last record)
    last = {}
    for i, p in enumerate(kp):
        last[int(p)] = i
    rows = [last[p] for p in sorted(last)]
    tt = np.r_[0.0, t[rows]]
    dt = np.diff(tt)

    def rate(pfx, side):
        cols = [c for c in cum.columns if c.endswith('_' + side)
                and c.split('_')[0].split('-')[0].rstrip('0123456789') == pfx]
        if not cols:
            return None
        v = np.r_[0.0, sum(cum[c].to_numpy(float)[rows] for c in cols)]
        r = np.diff(v) / np.where(dt > 0, dt, np.nan) / area_m2 * 1000.0
        return r[1:] if steady and r.size > 1 else r

    out = {'gw_to_sfr': rate('SFR', 'OUT'), 'sfr_to_gw': rate('SFR', 'IN'),
           'gw_to_lak': rate('LAK', 'OUT'), 'lak_to_gw': rate('LAK', 'IN'),
           'q_outlet': None}
    if out['gw_to_sfr'] is None and out['gw_to_lak'] is None:
        return None
    try:
        per = PP.sfr_observations(ws, name)
        if per is not None and 'outflow' in per.columns:
            out['q_outlet'] = (-per['outflow'].to_numpy(float) / area_m2
                               * 1000.0)
    except Exception:                                # pragma: no cover
        pass
    if nper:
        for k, v in out.items():
            if v is not None:
                out[k] = v[-int(nper):] if len(v) >= nper else None
    return out


def surface_water_terms(ws, area_m2, nper=None):
    """The same over the run, in mm/yr -- annual-equivalent as the table's
    soil-column rows (the sum over the periods / their number * 365.25)."""
    s = surface_water_series(ws, area_m2, nper)
    if s is None:
        return None
    return {k: (None if v is None else
                float(np.nansum(v) / max(len(v), 1) * 365.25))
            for k, v in s.items()}


def surface_water_lines(sw, labels, newv, refv):
    """The table rows of surface_water_terms, and the two totals that
    compare like with like: NWT routed no stream, so its valley groundwater
    exfiltrated into the soil column and partly ran off -- its Ro and EXFg
    already hold what here enters the streams through their beds."""
    if not sw:
        return []
    val = dict(zip(labels, newv))
    ref = dict(zip(labels, refv)) if refv is not None else {}
    z = lambda k: sw.get(k) or 0.0                    # noqa: E731
    rows = [('GW->SFR', z('gw_to_sfr'), None, 'groundwater into the streams'),
            ('SFR->GW', z('sfr_to_gw'), None, 'stream losses to the aquifer'),
            ('GW->LAK', z('gw_to_lak'), None, 'groundwater into the ponds'),
            ('LAK->GW', z('lak_to_gw'), None, 'pond losses to the aquifer')]
    if sw.get('q_outlet') is not None:
        rows.append(('Qout', sw['q_outlet'], None, 'streamflow at the outlet'))
    if 'Ro' in val:
        rows.append(('Ro+netSFR', val['Ro'] + z('gw_to_sfr') - z('sfr_to_gw'),
                     ref.get('Ro'), 'Ro + GW->SFR - SFR->GW, against NWT Ro'))
    if 'EXFg' in val:
        rows.append(('EXFg+GWsw', val['EXFg'] + z('gw_to_sfr') + z('gw_to_lak'),
                     ref.get('EXFg'),
                     'EXFg + GW->SFR + GW->LAK, against NWT EXFg'))
    out = ['', 'streams and ponds -- they exchange with the aquifer directly, '
               'not through the soil column:']
    for lab, a, b, what in rows:
        if refv is None:
            out.append('%-10s %12.1f   %s' % (lab, a, what))
        elif b is None:
            out.append('%-10s %12.1f %12s %12s   %s' % (lab, a, '-', '-', what))
        else:
            out.append('%-10s %12.1f %12.1f %12.1f   %s'
                       % (lab, a, b, a - b, what))
    return out


def write_summary(labels, newv, refv, ws, new, sw=None):
    lines = ['La Mata water budget - annual-equivalent rates (mm/year), '
             'over the %d simulated day(s)%s' % (
                 new['wb_ts'].shape[0],
                 ', both runs' if refv is not None else ''), '']
    if refv is not None:
        lines.append('%-10s %12s %12s %12s' % ('flux', 'MF6 (API)', 'NWT (Picard)', 'diff'))
        lines.append('-' * 50)
        for lab, a, b in zip(labels, newv, refv):
            lines.append('%-10s %12.1f %12.1f %12.1f' % (lab, a, b, a - b))
    else:
        lines.append('%-10s %12s' % ('flux', 'MF6 (API)'))
        lines.append('-' * 26)
        for lab, a in zip(labels, newv):
            lines.append('%-10s %12.1f' % (lab, a))
    lines += surface_water_lines(sw, labels, newv, refv)
    lines += ['', 'MF6 %s iterations per stress period: mean %.1f, max %d'
              % (new.get('iters_kind', 'outer'), new['outer_iters'].mean(),
                 new['outer_iters'].max())]
    # area-weighted: a plain mean over mesh cells is the stream corridor's
    w = new.get('cell_area')
    rej = (float(np.average(new['rejinf'], axis=1, weights=w).mean())
           if w is not None else float(new['rejinf'].mean())) * 365.25
    lines.append('rejected infiltration returned to the soil column: %.1f mm/year' % rej)
    txt = '\n'.join(lines)
    with open(_fig(ws, '00_summary.txt'), 'w') as f:
        f.write(txt + '\n')
    print('\n' + txt)


def make_figures(ws, no_reference=False, verbose=True, out_dir=None,
                 dataset_dir=None):
    """Build the 01-07 water-budget figures for a coupled run.

    ``ws`` is the MODFLOW 6 workspace the results are READ from; ``out_dir`` is
    the run's results folder they are WRITTEN to (<ws-root>/out_<stamp>_<tag>/),
    keeping output out of the model workspace and out of the repository. When
    ``out_dir`` is None the figures land next to the model, as before.
    ``dataset_dir`` is the run's dataset (None: mm_paths.dataset_dir).

    Callable from the runner's --postproc as well as the CLI. Returns the
    figures directory, or None if the results file has no water-budget arrays.
    """
    global _OUT_ROOT, DS
    _OUT_ROOT = os.path.abspath(out_dir) if out_dir else None
    if dataset_dir:
        DS = str(dataset_dir)
    ws = os.path.abspath(ws)
    new = load_new(ws)
    if new['wb_ts'] is None:
        if verbose:
            print('plot_water_budget: results file has no wb_ts/wb_map; skipped.')
        return None
    ref = None
    if not no_reference:
        try:
            ref = load_reference(ndays=new['wb_ts'].shape[0])
        except Exception as exc:                # pragma: no cover
            if verbose:
                print('plot_water_budget: reference not loaded (%r); '
                      'new-run figures only.' % exc)
    if verbose:
        print('plot_water_budget: %d SPs x %d cells%s'
              % (new['wb_ts'].shape[0], new['perc'].shape[1],
                 '' if ref is None else '   |   NWT reference loaded'))
    _area = new.get('cell_area')
    _area = float(np.sum(_area)) if _area is not None else None
    new['sw'] = surface_water_series(ws, _area, nper=new['wb_ts'].shape[0])
    plot_timeseries(new, ref, ws)
    plot_cumulative(new, ref, ws)
    labels, newv, refv = plot_totals(new, ref, ws)
    plot_hydroyears(new, ref, ws)
    plot_maps(new, ref, ws)
    plot_heads(new, ws, ref)
    plot_coupling(new, ws)
    sw = (None if not new.get('sw') else
          {k: (None if v is None else
               float(np.nansum(v) / max(len(v), 1) * 365.25))
           for k, v in new['sw'].items()})
    write_summary(labels, newv, refv, ws, new, sw=sw)
    figdir = os.path.join(_OUT_ROOT or ws, 'figures_nwt_comparison')
    if verbose:
        print('Figures written to %s' % figdir)
    return figdir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ws', default=os.path.join(str(_mmp.WS_ROOT), 'MF6_ws'),
                    help='the MODFLOW 6 workspace (default: mm_paths.WS_ROOT/MF6_ws)')
    ap.add_argument('--out-dir', default=None,
                    help='results folder to write figures into')
    ap.add_argument('--no-reference', action='store_true')
    a = ap.parse_args()
    if make_figures(a.ws, a.no_reference, out_dir=a.out_dir) is None:
        sys.exit('This results file predates the water-budget output.\n'
                 'Re-run tests/run_lamata_mf6.py to record wb_ts/wb_map.')


if __name__ == '__main__':
    main()
