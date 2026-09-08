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
    python tests/plot_water_budget.py --mode iterative
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..', 'trunk'))
DS = os.path.abspath(os.path.join(HERE, '..', 'DataSet_LaMata'))

# MODFLOW-NWT reference for the comparison panels. It lives in the legacy
# PhD/paper archive, NOT in the repository: it is 1.3 GB and cannot be
# regenerated (the NWT build path was removed from the code in Phase 1).
NWT_REF = os.environ.get(
    'MARMITES_NWT_REF',
    os.path.join('E:' + os.sep, '00code_ws', 'LaMata_new_PhD_artigo_2s3L',
                 '_h5_MM.h5'))
sys.path.insert(0, TRUNK)

import matplotlib  # noqa: E402
matplotlib.use('agg')
import matplotlib.pyplot as plt  # noqa: E402
import h5py  # noqa: E402
from marmites_indices import INDEX_MM  # noqa: E402
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


def load_new(ws, mode):
    fn = os.path.join(ws, '_coupled_%s.h5' % mode)
    if not os.path.exists(fn):
        sys.exit('coupled results not found: %s\n(run tests/run_lamata_mf6.py first)' % fn)
    with h5py.File(fn, 'r') as h:
        d = {k: h[k][:] for k in ('heads', 'perc', 'etg', 'exf', 'rejinf',
                                  'outer_iters', 'cell_ij')}
        for k in ('wb_ts', 'wb_map', 'wb_ts_soil', 'wb_map_soil'):
            d[k] = h[k][:] if k in h else None
    return d


def load_reference(chunk=120):
    """Catchment-mean time series and time-mean maps from the NWT run.

    Read in day chunks: the reference MM array is ~370 MB.
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


def _scatter(new, ref, values):
    """Per-cell values onto the full (nrow, ncol) grid, NaN elsewhere."""
    nrow, ncol = grid_shape(new, ref)
    ij = new['cell_ij']
    g = np.full((nrow, ncol), np.nan)
    g[ij[:, 0], ij[:, 1]] = values
    return g


# --------------------------------------------------------------------- #

def plot_timeseries(new, ref, ws, dates=None):
    keys = ['P', 'Ro', 'ETsoil', 'ETg', 'Rp']
    fig, axes = plt.subplots(len(keys), 1, figsize=(12, 11), sharex=True)
    x = np.arange(new['wb_ts'].shape[0])
    for ax, k in zip(axes, keys):
        idx = INDEX_MM[dict((a, b) for a, b, _ in FLUXES)[k]]
        ax.plot(x, new['wb_ts'][:, idx], lw=0.7, color='tab:blue', label='MF6 (API)')
        if ref is not None:
            ax.plot(np.arange(ref['ts'].shape[0]), ref['ts'][:, idx], lw=0.7,
                    color='tab:red', alpha=0.75, label='NWT (Picard)')
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
    fig, ax = plt.subplots(figsize=(11, 6))
    cmap = plt.get_cmap('tab10')
    for c, k in enumerate(keys):
        idx = INDEX_MM[dict((a, b) for a, b, _ in FLUXES)[k]]
        ax.plot(np.cumsum(new['wb_ts'][:, idx]), color=cmap(c), lw=1.6, label='%s (MF6)' % k)
        if ref is not None:
            ax.plot(np.cumsum(ref['ts'][:, idx]), color=cmap(c), lw=1.2, ls='--',
                    alpha=0.8, label='%s (NWT)' % k)
    ax.set_xlabel('day of simulation')
    ax.set_ylabel('cumulative depth (mm)')
    ax.set_title('Cumulative water-budget components (solid: MF6/API, dashed: NWT/Picard)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(_fig(ws, '02_wb_cumulative.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)


def plot_totals(new, ref, ws):
    labels, newv, refv = [], [], []
    nday_new = new['wb_ts'].shape[0]
    for short, key, _long in FLUXES:
        idx = INDEX_MM[key]
        labels.append(short)
        newv.append(new['wb_ts'][:, idx].sum() / nday_new * 365.25)
        if ref is not None:
            refv.append(ref['ts'][:, idx].sum() / ref['ts'].shape[0] * 365.25)
    x = np.arange(len(labels))
    w = 0.38 if ref is not None else 0.6
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(x - (w / 2 if ref is not None else 0), newv, w, label='MF6 (API)', color='tab:blue')
    if ref is not None:
        ax.bar(x + w / 2, refv, w, label='NWT (Picard)', color='tab:red', alpha=0.85)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('mm/year')
    ax.set_title('La Mata water budget - annual-equivalent rates')
    ax.grid(alpha=0.3, axis='y')
    ax.legend()
    for i, v in enumerate(newv):
        ax.annotate('%.0f' % v, (x[i] - (w / 2 if ref is not None else 0), v),
                    ha='center', va='bottom' if v >= 0 else 'top', fontsize=7)
    fig.tight_layout()
    fig.savefig(_fig(ws, '03_wb_totals.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)
    return labels, newv, (refv if ref is not None else None)


def plot_hydroyears(new, ref, ws, start_doy=273, days_per_year=365):
    """Totals per hydrological year (starting 1 October by default)."""
    keys = ['P', 'Ro', 'ETsoil', 'ETg', 'Rp']
    n = new['wb_ts'].shape[0]
    edges = list(range(0, n, days_per_year))
    fig, ax = plt.subplots(figsize=(11, 6))
    width = 0.8 / len(keys)
    xs = np.arange(len(edges))
    for c, k in enumerate(keys):
        idx = INDEX_MM[dict((a, b) for a, b, _ in FLUXES)[k]]
        tot = [new['wb_ts'][a:min(a + days_per_year, n), idx].sum() for a in edges]
        ax.bar(xs + c * width, tot, width, label='%s (MF6)' % k)
        if ref is not None:
            tr = [ref['ts'][a:min(a + days_per_year, n), idx].sum() for a in edges]
            ax.plot(xs + c * width, tr, 'k_', ms=14, mew=1.5)
    ax.set_xticks(xs + 0.4); ax.set_xticklabels(['yr %d' % (i + 1) for i in xs])
    ax.set_ylabel('mm/year')
    ax.set_title('Per-year budget totals (bars: MF6/API; black dashes: NWT/Picard)')
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
        fig.savefig(_fig(ws, '05_map_%s.png' % short), dpi=140, bbox_inches='tight')
        plt.close(fig)


def plot_heads(new, ws, ref=None):
    """Delegates to the recovered MARMITESplot module (single source of
    truth for MARMITES figures)."""
    MMplot.plotHEADS(new['heads'], new['cell_ij'], grid_shape(new, ref),
                     _fig(ws, '06_heads.png'))


def plot_coupling(new, ws):
    """Delegates to the recovered MARMITESplot module."""
    MMplot.plotCOUPLING(new['outer_iters'], rejinf=new['rejinf'],
                        exf=new.get('exf'), plt_export_fn=_fig(ws, '07_coupling.png'))


def write_summary(labels, newv, refv, ws, new):
    lines = ['La Mata water budget - annual-equivalent rates (mm/year)', '']
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
    lines += ['', 'MF6 outer iterations: mean %.1f, max %d'
              % (new['outer_iters'].mean(), new['outer_iters'].max())]
    rej = new['rejinf'].mean() * 365.25
    lines.append('rejected infiltration returned to the soil column: %.1f mm/year' % rej)
    txt = '\n'.join(lines)
    with open(_fig(ws, '00_summary.txt'), 'w') as f:
        f.write(txt + '\n')
    print('\n' + txt)


def make_figures(ws, mode='lagged', no_reference=False, verbose=True,
                 out_dir=None):
    """Build the 01-07 water-budget figures for a coupled run.

    ``ws`` is the MODFLOW 6 workspace the results are READ from; ``out_dir`` is
    the run's results folder they are WRITTEN to (<ws-root>/out_<stamp>_<tag>/),
    keeping output out of the model workspace and out of the repository. When
    ``out_dir`` is None the figures land next to the model, as before.

    Callable from the runner's --postproc as well as the CLI. Returns the
    figures directory, or None if the results file has no water-budget arrays.
    """
    global _OUT_ROOT
    _OUT_ROOT = os.path.abspath(out_dir) if out_dir else None
    ws = os.path.abspath(ws)
    new = load_new(ws, mode)
    if new['wb_ts'] is None:
        if verbose:
            print('plot_water_budget: results file has no wb_ts/wb_map; skipped.')
        return None
    ref = None
    if not no_reference:
        try:
            ref = load_reference()
        except Exception as exc:                # pragma: no cover
            if verbose:
                print('plot_water_budget: reference not loaded (%r); '
                      'new-run figures only.' % exc)
    if verbose:
        print('plot_water_budget: %d SPs x %d cells%s'
              % (new['wb_ts'].shape[0], new['perc'].shape[1],
                 '' if ref is None else '   |   NWT reference loaded'))
    plot_timeseries(new, ref, ws)
    plot_cumulative(new, ref, ws)
    labels, newv, refv = plot_totals(new, ref, ws)
    plot_hydroyears(new, ref, ws)
    plot_maps(new, ref, ws)
    plot_heads(new, ws, ref)
    plot_coupling(new, ws)
    write_summary(labels, newv, refv, ws, new)
    figdir = os.path.join(_OUT_ROOT or ws, 'figures_nwt_comparison')
    if verbose:
        print('Figures written to %s' % figdir)
    return figdir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ws', default=os.environ.get(
        'MARMITES_WS_ROOT', os.path.join('E:' + os.sep, '00code_ws', 'LaMata_MM-MF6')) + os.sep + 'MF6_ws')
    ap.add_argument('--out-dir', default=None,
                    help='results folder to write figures into')
    ap.add_argument('--mode', default='lagged')
    ap.add_argument('--no-reference', action='store_true')
    a = ap.parse_args()
    if make_figures(a.ws, a.mode, a.no_reference, out_dir=a.out_dir) is None:
        sys.exit('This results file predates the water-budget output.\n'
                 'Re-run tests/run_lamata_mf6.py to record wb_ts/wb_map.')


if __name__ == '__main__':
    main()
