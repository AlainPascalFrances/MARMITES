# -*- coding: utf-8 -*-
"""Post-processing for the coupled MARMITES / MODFLOW 6 run.

Reads the MODFLOW 6 output already on disk (heads, cell budget, UZF and SFR
budgets, the listing) and the coupled-run HDF5, and writes figures + tidy CSVs
into ``<sim_ws>/postproc/``. Mirrors the CdL post-processing
(trunk/SFR_LAK_CRR/postprocess_cdl.py) but for La Mata's structured DIS grid
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

import os
import sys
import time
import warnings

import numpy as np

__all__ = ['run_postproc', 'run_preproc', 'obs_points', 'obs_series',
           'budget_by_compartment', 'package_budget', 'layer_storage_change',
           'native_layer_maps', 'native_suite', 'COMPARTMENT']

# MARMITES input maps (soil-water-balance side) and the MODFLOW input maps
# (aquifer side). Both belong in the preprocessing channel: the MM script
# itself produces the soil/vegetation/meteo maps, and the aquifer maps come
# from the MF workspace. (label, filename-relative-to, colormap, integer?)
_MM_INPUT_MAPS = [
    ('soil zones', 'inputSOILzones.asc', 'tab20', True),
    ('meteo zones', 'inputMETEOzones.asc', 'tab10', True),
    ('irrigation zones', 'inputIRRzones.asc', 'tab10', True),
    ('soil thickness [m]', 'inputSOILthick.asc', 'YlOrBr', False),
    ('pond width [m]', 'inputPONDw.asc', 'Blues', False),
    ('pond depth max [m]', 'inputPONDhmax.asc', 'Blues', False),
    ('vegetation 1 area', 'inputVEG1area.asc', 'Greens', False),
    ('vegetation 2 area', 'inputVEG2area.asc', 'Greens', False),
    ('vegetation 3 area', 'inputVEG3area.asc', 'Greens', False),
]

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


def _mkdir(sim_ws, sub='postproc'):
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
    """Produce the full post-processing figure/CSV set into <sim_ws>/postproc/.

    Returns the list of files written. Each figure is guarded so a missing
    output file (e.g. no SFR in this run) skips that figure rather than
    aborting the whole report.
    """
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    import flopy

    out = _mkdir(out_root or sim_ws, 'postproc')
    written = []
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, verbosity_level=0)
    gwf = sim.get_model()
    mg = gwf.modelgrid
    nlay, nrow, ncol = mg.nlay, mg.nrow, mg.ncol
    top = np.asarray(mg.top, dtype=float).reshape(nrow, ncol)

    hds = flopy.utils.HeadFile(os.path.join(sim_ws, '%s.hds' % name))
    kk = hds.get_kstpkper()
    kk_real = [k for k in kk if k[1] >= 1] or kk        # drop steady state
    # the mean map is taken over an even subsample; the obs series keeps the
    # full record (one cell, cheap) via _obs_series_from_hds
    kk_map = _subsample(kk_real, 200)
    H = np.array([hds.get_data(kstpkper=k) for k in kk_map])   # (nt, nlay, nrow, ncol)
    H = np.where(np.abs(H) > 1e29, np.nan, H)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)  # all-NaN layers
        mean_h = np.nanmean(H, axis=0)

    if dates is None:
        dates = _dates_from_dataset(ds_ws, len(kk_real))

    # --- obs vs computed heads ---------------------------------------- #
    try:
        f = _fig_obs_heads(hds, ds_ws, kk_real, dates, top, nlay, nrow, ncol,
                           xll, yll, cs, out)
        written += f
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   obs-vs-computed heads skipped: %r' % exc)

    # --- mean head + depth-to-water maps ------------------------------ #
    for L in range(nlay):
        for kind, arr, cmap, lbl in (
                ('mean_head', mean_h[L], 'viridis', 'mean head [m]'),
                ('mean_depth', top - mean_h[L], 'viridis_r',
                 'mean depth to water [m]')):
            fig, ax = plt.subplots(figsize=(6, 6))
            im = ax.imshow(arr, cmap=cmap)
            ax.set_title('%s  layer %d' % (lbl, L + 1))
            fig.colorbar(im, ax=ax, shrink=0.8, label=lbl)
            fn = os.path.join(out, '%s_L%d.png' % (kind, L + 1))
            fig.savefig(fn, dpi=140, bbox_inches='tight')
            plt.close(fig)
            np.savetxt(os.path.join(out, '%s_L%d.csv' % (kind, L + 1)),
                       arr, delimiter=',')
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

    # --- native MARMITES plotLAYER head maps -------------------------- #
    try:
        written += native_layer_maps(sim_ws, name=name, verbose=verbose)
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   native plotLAYER skipped: %r' % exc)

    # --- per-layer storage change map --------------------------------- #
    try:
        sto = layer_storage_change(sim_ws, name, nlay, nrow, ncol)
        for L in range(nlay):
            fig, ax = plt.subplots(figsize=(6, 6))
            vmax = np.nanmax(np.abs(sto[L])) or 1.0
            im = ax.imshow(sto[L], cmap='RdBu', vmin=-vmax, vmax=vmax)
            ax.set_title('mean storage change  layer %d' % (L + 1))
            fig.colorbar(im, ax=ax, shrink=0.8, label='m3/d (+ gain)')
            fn = os.path.join(out, 'storage_change_L%d.png' % (L + 1))
            fig.savefig(fn, dpi=140, bbox_inches='tight')
            plt.close(fig)
            written += [fn]
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   storage-change map skipped: %r' % exc)

    if verbose:
        print('postproc: %d file(s) written to %s' % (len(written), out))
    return written


def run_preproc(sim_ws, ds_ws, name='lamatamm', mf_ws=None, verbose=True,
                out_root=None):
    """Input maps into <sim_ws>/preproc/: MARMITES soil/veg/meteo maps AND the
    MODFLOW aquifer maps (top, per-layer K / Ss / Sy / thickness, ibound+UZF
    footprint, ponds+stream overlay).
    """
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import flopy

    out = _mkdir(out_root or sim_ws, 'preproc')
    mf_ws = mf_ws or os.path.join(ds_ws, 'MF_ws')
    written = []

    def _one(arr, title, cmap, fn, integer=False):
        fig, ax = plt.subplots(figsize=(6, 6))
        a = np.asarray(arr, dtype=float)
        im = ax.imshow(a, cmap=cmap)
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.8)
        path = os.path.join(out, fn)
        fig.savefig(path, dpi=140, bbox_inches='tight')
        plt.close(fig)
        written.append(path)

    # --- MARMITES input maps (soil-water-balance side) ---------------- #
    for label, fn, cmap, integer in _MM_INPUT_MAPS:
        p = os.path.join(ds_ws, fn)
        if os.path.exists(p):
            _one(_asc(p), 'MM input: %s' % label, cmap,
                 'mm_%s.png' % fn.replace('input', '').replace('.asc', '').lower())

    # --- MODFLOW aquifer maps ----------------------------------------- #
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, verbosity_level=0)
    gwf = sim.get_model()
    mg = gwf.modelgrid
    nlay, nrow, ncol = mg.nlay, mg.nrow, mg.ncol
    top = np.asarray(mg.top, dtype=float).reshape(nrow, ncol)
    botm = np.asarray(mg.botm, dtype=float).reshape(nlay, nrow, ncol)
    _one(top, 'aquifer top / land surface [m]', 'terrain', 'aq_top.png')

    npf = gwf.get_package('npf')
    if npf is not None:
        k = np.asarray(npf.k.array).reshape(nlay, nrow, ncol)
        k33 = np.asarray(npf.k33.array).reshape(nlay, nrow, ncol)
        for L in range(nlay):
            _one(k[L], 'K horizontal L%d [m/d]' % (L + 1), 'viridis',
                 'aq_k_L%d.png' % (L + 1))
            _one(k33[L], 'K vertical L%d [m/d]' % (L + 1), 'viridis',
                 'aq_k33_L%d.png' % (L + 1))
            _one(botm[L], 'bottom L%d [m]' % (L + 1), 'terrain',
                 'aq_botm_L%d.png' % (L + 1))

    sto = gwf.get_package('sto')
    if sto is not None:
        ss = np.asarray(sto.ss.array).reshape(nlay, nrow, ncol)
        sy = np.asarray(sto.sy.array).reshape(nlay, nrow, ncol)
        for L in range(nlay):
            _one(ss[L], 'specific storage L%d [1/m]' % (L + 1), 'magma',
                 'aq_ss_L%d.png' % (L + 1))
            _one(sy[L], 'specific yield L%d' % (L + 1), 'cividis',
                 'aq_sy_L%d.png' % (L + 1))

    # ibound / idomain + UZF footprint
    idom = np.asarray(mg.idomain).reshape(nlay, nrow, ncol) \
        if mg.idomain is not None else np.ones((nlay, nrow, ncol))
    _one(idom[0], 'idomain / active cells L1', 'Greys', 'aq_idomain_L1.png')

    # ponds + stream overlay
    try:
        written += _fig_network_overlay(sim_ws, ds_ws, gwf, top, out)
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   network overlay skipped: %r' % exc)

    if verbose:
        print('preproc: %d file(s) written to %s' % (len(written), out))
    return written


def _fig_network_overlay(sim_ws, ds_ws, gwf, top, out):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 7))
    im = ax.imshow(top, cmap='terrain', alpha=0.8)
    fig.colorbar(im, ax=ax, shrink=0.8, label='elevation [m]')
    pondw = os.path.join(ds_ws, 'inputPONDw.asc')
    if os.path.exists(pondw):
        w = _asc(pondw)
        ii, jj = np.where(w > 0)
        ax.plot(jj, ii, 's', ms=2, color='tab:blue', label='stream cells')
    sfr = gwf.get_package('sfr')
    if sfr is not None:
        cd = sfr.packagedata.get_data()
        ri = [rec['cellid'][-2] for rec in cd]
        cj = [rec['cellid'][-1] for rec in cd]
        ax.plot(cj, ri, '.', ms=3, color='navy', label='SFR reaches')
    lak = gwf.get_package('lak')
    if lak is not None:
        cd = lak.connectiondata.get_data()
        li = [rec['cellid'][-2] for rec in cd]
        lj = [rec['cellid'][-1] for rec in cd]
        ax.plot(lj, li, 'o', ms=6, mfc='none', color='tab:red', label='ponds (LAK)')
    ax.set_title('stream network and ponds')
    ax.legend(fontsize=8, loc='best')
    ax.invert_yaxis()
    fn = os.path.join(out, 'network_overlay.png')
    fig.savefig(fn, dpi=140, bbox_inches='tight')
    plt.close(fig)
    return [fn]


# labels for the MM flux indices, used by the native catchment plot
_MM_LABEL = {
    'iP': 'P', 'iPe': 'Pe', 'iPT': 'PT', 'iPE': 'PE', 'iRo': 'Ro',
    'iI': 'I', 'iEXFg': 'EXFg', 'iEow': 'Eow', 'iEi': 'Ei', 'iEg': 'Eg',
    'iTg': 'Tg', 'iETg': 'ETg', 'iETsoil': 'ETsoil', 'iperc': 'Rp',
    'idSsurf': 'dSsurf', 'iSsurf': 'Ssurf', 'idSsoil': 'dSsoil',
    'iSsoil_pc': 'Ssoil%', 'iMB': 'MB', 'iEo': 'Eo', 'ihcorr': 'hcorr',
    'idgwt': 'dgwt', 'iuzthick': 'uzthick', 'iMBsurf': 'MBsurf',
}


def native_suite(out_dir, cMF, ctx, res, ds_ws=None, trunk=None, verbose=True,
                 sim_ws=None, sankey=True, sankey_full=True, sankey_min_flux=0.05):
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
    if trunk is None:
        trunk = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    mmplot_dir = os.path.join(trunk, 'MARMITESutilities', 'MARMITESplot')
    if mmplot_dir not in sys.path:
        sys.path.insert(0, mmplot_dir)
    import MARMITESplot_v3 as MMplot
    IX = dict(ctx.index)
    nper = res['perc'].shape[0]
    if sim_ws is None:
        sim_ws = os.path.dirname(os.path.abspath(out_dir))
    name = str(getattr(cMF, 'modelname', 'lamatamm')).lower()
    written = []

    # --- catchment water-balance time series -------------------------- #
    # plotTIMESERIES_CATCH expects the driver's COMBINED catchment-flux array:
    # the per-cell MM fluxes (wb_ts) plus the soil-layer fluxes summed over
    # layers (wb_ts_soil) under extra index keys (iEsoil, iTsoil, iRsoil), plus
    # a few derived ones (iRg == recharge). Reassemble that layout here.
    try:
        fn = _native_catchment_series(MMplot, out_dir, cMF, ctx, res, nper)
        if fn:
            written.append(fn)
    except Exception as exc:                         # pragma: no cover
        if verbose:
            print('   native catchment WB series skipped: %r' % exc)

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
                                cache_dir=out_dir, verbose=verbose)
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
                                          verbose=verbose, agg=agg)
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

    # --- per-layer maps of the time-mean fluxes ----------------------- #
    try:
        written += _native_flux_maps(MMplot, out_dir, cMF, ctx, res, verbose=verbose)
    except Exception as exc:                         # pragma: no cover
        if verbose:
            print('   native flux maps skipped: %r' % exc)

    if verbose:
        print('native MARMITESplot: %d figure(s) -> %s' % (len(written), out_dir))
    return written


def _native_catchment_series(MMplot, out_dir, cMF, ctx, res, nper):
    """Assemble the driver's combined catchment-flux array and call
    plotTIMESERIES_CATCH. Returns the file path or None."""
    IX = dict(ctx.index)
    IXS = dict(ctx.index_S)
    wb_ts = np.asarray(res['wb_ts'])                # (nper, nMM)
    wb_ts_soil = np.asarray(res['wb_ts_soil'])      # (nper, nsl, nsoil)

    comb = dict(IX)                                 # start from the MM indices
    rows = [wb_ts[:, i] for _, i in sorted(IX.items(), key=lambda kv: kv[1])]
    # append the soil-layer fluxes, summed over layers, under the keys the
    # native routine references
    for key, skey in (('iEsoil', 'iEsoil'), ('iTsoil', 'iTsoil'),
                      ('iRsoil', 'iRsoil')):
        if skey in IXS and key not in comb:
            comb[key] = len(rows)
            rows.append(wb_ts_soil[:, :, IXS[skey]].sum(axis=1))
    # derived/aliased indices the routine also uses
    if 'iRg' not in comb:                           # groundwater recharge ~ perc
        comb['iRg'] = len(rows)
        rows.append(wb_ts[:, IX['iperc']] if 'iperc' in IX else np.zeros(nper))
    # observation overlays are optional (obs_catch=None); give them safe slots
    for key in ('iRoobs', 'ihobs', 'idobs', 'ih_SF', 'idcorr', 'idSg_1'):
        if key not in comb:
            comb[key] = len(rows)
            rows.append(np.zeros(nper))
    flx = np.vstack(rows)                           # (nflux, nper)
    flxLbl = [''] * len(comb)
    for k, i in comb.items():
        flxLbl[i] = _MM_LABEL.get(k, k[1:] if k.startswith('i') else k)

    dates = getattr(cMF, 'inputDate', None)
    if dates is None or len(np.atleast_1d(dates)) != nper:
        import pandas as pd
        dates = pd.date_range('2000-01-01', periods=nper, freq='D')
        cMF = _ShimDate(cMF, dates)
    d0, d1 = _to_ordinal(dates[0]), _to_ordinal(dates[-1])
    fn = os.path.join(out_dir, 'native_wb_catchment.png')
    MMplot.plotTIMESERIES_CATCH(
        cMF, flx, flxLbl, fn, 'catchment_wb',
        float(np.nanmax(flx)), float(np.nanmin(flx)),
        int(getattr(cMF, 'iniMonthHydroYear', 10)), d0, d1, comb)
    return fn if os.path.exists(fn) else None


class _ShimDate(object):
    """Thin wrapper adding inputDate to a cMF that lacks it (fallback only)."""
    def __init__(self, cMF, dates):
        self._c = cMF
        self.inputDate = dates

    def __getattr__(self, name):
        return getattr(self._c, name)


def _to_ordinal(d):
    try:
        import matplotlib.dates as mdates
        return mdates.date2num(d)
    except Exception:
        return float(getattr(d, 'toordinal', lambda: 0)())


# Time-mean MM flux maps to draw, as (index key, file/label stem, colourbar).
_MAP_FLUXES = (
    ('iP', 'P', 'rainfall'),
    ('iPe', 'Pe', 'effective rainfall'),
    ('iRo', 'Ro', 'runoff'),
    ('iI', 'I', 'infiltration'),
    ('iEi', 'Ei', 'interception'),
    ('iEow', 'Eow', 'open-water evaporation'),
    ('iETsoil', 'ETsoil', 'soil ET'),
    ('iEg', 'Eg', 'groundwater evaporation'),
    ('iTg', 'Tg', 'groundwater transpiration'),
    ('iETg', 'ETg', 'groundwater ET'),
    ('iperc', 'Rp', 'percolation'),
    ('iEXFg', 'EXFg', 'exfiltration'),
    ('idSsoil', 'dSsoil', 'change in soil storage'),
    ('idSsurf', 'dSsurf', 'change in surface storage'),
    ('iSsoil_pc', 'theta', 'soil moisture'),
    ('iuzthick', 'uzthick', 'unsaturated thickness'),
    ('idgwt', 'dgwt', 'depth to water table'),
)


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


def _native_flux_maps(MMplot, out_dir, cMF, ctx, res, verbose=True):
    """plotLAYER maps of the time-mean MM fluxes.

    Follows the conventions the legacy driver used for its input/output maps
    (startMARMITES_v3.py ~808-910): ``Date='NA'``/``JD='NA'`` because these are
    time-MEAN maps rather than a day in a series, ``interval_type='linspace'``
    with 5 intervals, the model's own ``hnoflo``, and the observation points
    overlaid. It previously passed ``msg='arange'`` and a ``pref_plt_title``
    that duplicated the title, giving files named ``native_native_map_*``.
    """
    import matplotlib
    IX = dict(ctx.index)
    wb_map = np.asarray(res['wb_map'])              # (ncell, nidx)
    cells = ctx.cells
    nrow, ncol = int(cMF.nrow), int(cMF.ncol)
    hnoflo = float(getattr(cMF, 'hnoflo', 9999.999))
    cmap = matplotlib.colormaps['gist_rainbow_r']
    pts = _obs4map(res)
    before = set(os.listdir(out_dir)) if os.path.isdir(out_dir) else set()
    for key, stem, cblbl in _MAP_FLUXES:
        if key not in IX:
            continue
        grid = np.full((1, 1, nrow, ncol), hnoflo)
        for n, c in enumerate(cells):
            grid[0, 0, c[1], c[2]] = wb_map[n, IX[key]]
        mask = np.isclose(grid[0], hnoflo, atol=0.09)
        vals = grid[0][~mask]
        if not vals.size:
            continue
        unit = '-' if key in ('iSsoil_pc',) else (
            'm' if key in ('iuzthick', 'idgwt') else 'mm/d')
        try:
            MMplot.plotLAYER(
                days=[0], str_per=[0], Date='NA', JD='NA', ncol=ncol, nrow=nrow,
                nlay=1, nplot=1, V=grid, cmap=cmap,
                CBlabel='%s [%s]' % (cblbl, unit), msg='',
                plt_title='MMmap_%s' % stem, MM_ws=out_dir,
                interval_type='linspace', interval_num=5,
                Vmax=[float(vals.max())], Vmin=[float(vals.min())],
                fmt='%5.2f', points=pts, mask=mask, hnoflo=hnoflo)
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   flux map %s skipped: %r' % (stem, exc))
    after = set(os.listdir(out_dir)) if os.path.isdir(out_dir) else set()
    return [os.path.join(out_dir, f) for f in sorted(after - before)]


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
    convention). Cached to ``<cache_dir>/_aquifer_digest.npz``, keyed on the
    budget file's size and mtime, so re-plotting a run is instant.
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
    grb = os.path.join(sim_ws, '%s.dis.grb' % name)
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
    delr = np.asarray(cMF.delr, float)
    delc = np.asarray(cMF.delc, float)

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
    area_sel = float((delc[:, None] * delr[None, :])[mask].sum())
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
                   verbose):
    """Call plotWBsankey once and collect the PNGs it wrote for this fntitle."""
    written = []
    try:
        MMplot.plotWBsankey(out_dir, DATE, flx, flxIndex,
                            fn='%s_WBsankey' % fntitle, indexTime=HYindex,
                            year_lst=year_lst, cMF=smf, ncell_MM=ncell_MM,
                            obspt=obspt, fntitle=fntitle,
                            ibound4Sankey=ibound4Sankey, treshold=treshold)
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
                       min_flux=0.05, verbose=True, agg=None):
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
                                  names[p], 'obs_%s' % names[p], min_flux, verbose)
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


def _obs_head_series(sim_ws, name, nlay, cells_ij, nper):
    """Head time series at each obs cell, per layer: (nobs, nlay, nper).

    Read straight from the MODFLOW 6 head file rather than any HDF5 copy.
    """
    import flopy
    hds = flopy.utils.HeadFile(os.path.join(sim_ws, '%s.hds' % name))
    kk = hds.get_kstpkper()
    off = max(len(kk) - nper, 0)                   # drop the steady SP0
    out = np.full((len(cells_ij), nlay, nper), np.nan)
    for p, (i, j) in enumerate(cells_ij):
        ts = hds.get_ts([(L, int(i), int(j)) for L in range(nlay)])
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
    delr = np.asarray(cMF.delr, float)
    delc = np.asarray(cMF.delc, float)
    cf = _conv_fact(cMF)
    import MARMITESsoil_v3 as _MMsoil
    satflow = _MMsoil.SATFLOW()
    h_lbl = list(getattr(cMF, 'h_lbl', [str(L + 1) for L in range(nlay)]))
    written = []

    for p in range(nobs):
        i, j = int(obs_ij[p, 0]), int(obs_ij[p, 1])
        o = names[p]
        try:
            zone = int(ctx.gridSOIL[i, j]) - 1
            nsl = int(ctx._nsl[zone])
            l_high = int(cMF.outcropL[i, j]) - 1
            area = delr[j] * delc[i]
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

            hs = heads[p][np.isfinite(heads[p])]
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
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   obs time series (%s) skipped: %r' % (o, exc))
    if verbose:
        print('   native obs time series: %d/%d point(s)' % (len(written), nobs))
    return written


def resolve_obs_cells(cMF, ctx, ds_ws, verbose=True):
    """Map the enabled observation points (inputObs.txt) to MM cell-list
    positions. Returns ``(obs_idx, obs_names)`` for the coupler's obs capture.

    Reuses the native ``cPROCESS.inputObs`` reader (so the same points, grid
    mapping and disabled-line handling as the legacy driver), then finds each
    point's ``(i, j)`` in the MM cell list. Points whose cell is inactive are
    dropped with a warning."""
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


def native_layer_maps(sim_ws, name='lamatamm', trunk=None, verbose=True):
    """Render the mean-head map per layer with the native MARMITES plotLAYER.

    This reuses the original MARMITES figure style (the same routine the NWT
    driver used for its head maps) on the coupled MF6 output, so the native
    plotting suite is exercised alongside the cdl-style maps. Returns the files
    written, or [] if plotLAYER is unavailable.
    """
    import flopy
    if trunk is None:
        trunk = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    mmplot_dir = os.path.join(trunk, 'MARMITESutilities', 'MARMITESplot')
    if mmplot_dir not in sys.path:
        sys.path.insert(0, mmplot_dir)
    try:
        import MARMITESplot_v3 as MMplot
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   native plotLAYER unavailable: %r' % exc)
        return []
    out = _mkdir(sim_ws, 'postproc')
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, verbosity_level=0)
    mg = sim.get_model().modelgrid
    nlay, nrow, ncol = mg.nlay, mg.nrow, mg.ncol
    hds = flopy.utils.HeadFile(os.path.join(sim_ws, '%s.hds' % name))
    kk = [k for k in hds.get_kstpkper() if k[1] >= 1] or hds.get_kstpkper()
    H = np.array([hds.get_data(kstpkper=k) for k in _subsample(kk, 200)])
    H = np.where(np.abs(H) > 1e29, np.nan, H)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)  # all-NaN layers
        mean_h = np.nanmean(H, axis=0).reshape(1, nlay, nrow, ncol)
    mask = ~np.isfinite(mean_h[0])
    import matplotlib
    cmap = matplotlib.colormaps['viridis']       # plotLAYER wants a cmap object
    MMplot.plotLAYER(
        days=[0], str_per=[0], Date=['mean'], JD=[0], ncol=ncol, nrow=nrow,
        nlay=nlay, nplot=nlay, V=np.nan_to_num(mean_h, nan=-999.9),
        cmap=cmap, CBlabel='mean head [m]', msg='arange',
        plt_title='mean_head_native', MM_ws=out, mask=mask, hnoflo=-999.9,
        pref_plt_title='native')
    files = [os.path.join(out, f) for f in os.listdir(out)
             if f.startswith('native_mean_head_native')]
    if verbose:
        print('native plotLAYER: %d head map page(s) -> %s' % (len(files), out))
    return files


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
                   xll, yll, cs, out):
    import matplotlib.pyplot as plt
    import pandas as pd
    pts = obs_points(ds_ws)
    series = {p['name']: obs_series(ds_ws, p['name']) for p in pts}
    have = [p for p in pts if series.get(p['name']) is not None]
    if not have:
        return []
    dates = pd.to_datetime(dates)
    ncols = 2
    nrows_f = int(np.ceil(len(have) / ncols))
    fig, axes = plt.subplots(nrows_f, ncols, figsize=(12, 3 * nrows_f),
                             squeeze=False)
    tidy = []
    for ax, p in zip(axes.ravel(), have):
        i, j = _xy_to_ij(p['x'], p['y'], xll, yll, cs, nrow, ncol)
        L = min(max(p['lay'] - 1, 0), nlay - 1)
        # full-resolution series at this single cell (cheap via get_ts)
        ts = hds.get_ts((L, i, j))
        comp = ts[:, 1]
        comp = np.where(np.abs(comp) > 1e29, np.nan, comp)
        # get_ts returns every saved step (incl. steady); align to transient
        if comp.shape[0] == len(kk_real) + 1:
            comp = comp[1:]
        n = min(len(dates), comp.shape[0])
        dd, comp = dates[:n], comp[:n]
        ax.plot(dd, comp, '-', color='#5c3211', lw=1.0, label='computed L%d' % (L + 1))
        obs = series[p['name']]
        ax.plot(obs['date'], obs['head'], 'o', ms=4, color='tab:red', label='observed')
        ax.axhline(float(top[i, j]), color='0.6', ls=':', lw=0.8, label='land surface')
        ax.set_title('%s  (%d,%d)' % (p['name'], i, j))
        ax.set_ylabel('head [m]')
        ax.legend(fontsize=7, loc='best')
        for d, v in zip(dd, comp):
            tidy.append((p['name'], d, float(v)))
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
