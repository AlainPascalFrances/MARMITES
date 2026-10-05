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
  outlet_streamflow  streamflow at the catchment outlet (m3/d and mm/yr)
                     against the gauge, <ro_prefix>_catchment.txt (WP3.6)
  budget_sfr_ts      the stream network's budget over time, in mm
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
           'native_suite', 'COMPARTMENT', 'sfr_observations',
           'obs_streamflow', 'model_map_features']

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

# ---- where the measured series are: the State variables panel ------------
# These were LITERALS in five places of this module -- the point table, the
# head, soil-moisture and runoff prefixes, and the name column of the point
# layer -- so the panel that names them, [obs], changed nothing; the audit
# behind the cookbook's Appendix B found it. They are one table now, filled
# from the configuration by use_observations() before a run reads anything.
# The defaults are the La Mata names the literals were.
#
# The ACTUAL-ET series (obs.aet_prefix) is not here: no reader takes one yet.
OBS = {'table': 'inputObs.txt', 'heads': 'inputObsHEADS',
       'sm': 'inputObsSM', 'ro': 'inputObsRo', 'name_column': 'Name'}


def _tick_years(cMF):
    """The time-series tick density, from [postproc] via cMF.

    A record shorter than the first bound gets quarterly minor ticks, one
    shorter than the second half-yearly. The plotting routines have always
    taken these; nothing passed them, so the panel's values never arrived.
    """
    return {'maxYearsTickTrimester': int(getattr(cMF, 'maxYearsTickTrimester',
                                                 5)),
            'maxYearsTickSemester': int(getattr(cMF, 'maxYearsTickSemester',
                                                10))}


def use_observations(cfg):
    """Take the observation files from the State variables panel.

    Called once, by the driver, before anything reads them. Returns what is
    now in use, for the log.
    """
    o = cfg.obs
    OBS.update(table=o.table, heads=o.heads_prefix, sm=o.sm_prefix,
               ro=o.ro_prefix, name_column=o.name_column)
    return dict(OBS)


def obs_points(ds_ws, fn=None):
    """Active observation points from the point table ([obs] table).

    Lines starting with '#' or '##' are disabled points and are skipped, as in
    the MARMITES input convention. Returns [{name, x, y, lay}, ...].
    """
    fn = fn or OBS['table']
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


def obs_series(ds_ws, name, prefix=None):
    """Observed head time series (date, head) for one point, or None.

    ``prefix`` defaults to the panel's head prefix, joined to the point name
    with '_' as the files are named (inputObsHEADS_P1.txt).
    """
    import pandas as pd
    prefix = prefix if prefix is not None else OBS['heads'] + '_'
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

def steady_first(sim_ws, name=None):
    """Is stress period 1 of the model in ``sim_ws`` steady-state?

    Read from the model's own STO file -- its first PERIOD block says
    STEADY-STATE or TRANSIENT -- rather than assumed: a periodic spin-up
    cycle starts from the previous cycle's heads and has no steady period,
    and dropping kper 0 there would drop the run's first day. No STO file
    (a model without storage) is steady throughout.
    """
    import glob
    import re
    cands = ([os.path.join(sim_ws, '%s.sto' % name)] if name else []) +         sorted(glob.glob(os.path.join(sim_ws, '*.sto')))
    fn = next((c for c in cands if os.path.exists(c)), None)
    if fn is None:
        return True
    with open(fn, encoding='utf-8', errors='replace') as fh:
        txt = fh.read().upper()
    m = re.search(r'BEGIN\s+PERIOD\s+(\d+)(.*?)END\s+PERIOD', txt, re.S)
    if m is None or int(m.group(1)) != 1:
        return True          # MF6's default for a period never set: steady
    return 'STEADY-STATE' in m.group(2)


def budget_by_compartment(sim_ws, name):
    """Mean IN/OUT rate [m3/d] per compartment from the MF6 listing.

    Returns (df_terms, df_compartments) as pandas frames. Uses flopy's
    Mf6ListBudget so the numbers are exactly MF6's own accounting.
    """
    import pandas as pd
    import flopy
    lst = flopy.utils.Mf6ListBudget(os.path.join(sim_ws, '%s.lst' % name))
    inc, cum = lst.get_dataframes()
    # drop the steady-state first period, if there is one; average the
    # incremental rates
    rate = (inc.iloc[1:] if len(inc) > 1 and steady_first(sim_ws, name)
            else inc)
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


def period_steps(times, kk, steady=False):
    """Group a binary output's records into its STRESS PERIODS.

    ``times`` and ``kk`` are the file's record times [d] and (kstp, kper),
    one per record (flopy ``get_times`` / ``get_kstpkper``). MF6 saves a
    record at every TIME STEP, and adaptive time stepping cuts a period into
    as many as it needs -- 1354 records for 365 periods on 2026-10-05 -- so
    reading the last nper records as nper periods shifted every series in
    time (the CdL trap, cookbook WP6.4). Returns one list per stress period,
    in order, of (record index, weight), the weights the steps' shares of
    the period (summing to 1): a RATE is the period's time-weighted mean, a
    STATE its LAST record. The steady first period is left out when
    ``steady``.
    """
    t = np.asarray(times, dtype=float).ravel()
    kp = np.array([int(k[1]) for k in kk])
    dt = np.diff(np.r_[0.0, t])
    groups = {}
    for i, per in enumerate(kp):
        groups.setdefault(int(per), []).append(i)
    pers = sorted(groups)
    if steady and len(pers) > 1:
        pers = pers[1:]
    out = []
    for per in pers:
        idx = groups[per]
        w = dt[idx]
        tot = float(w.sum())
        w = w / tot if tot > 0 else np.full(len(idx), 1.0 / len(idx))
        out.append([(int(i), float(x)) for i, x in zip(idx, w)])
    return out


def period_end_kk(hds, steady=False):
    """The (kstp, kper) of the LAST step of every stress period."""
    kk = hds.get_kstpkper()
    return [kk[g[-1][0]] for g in period_steps(hds.get_times(), kk, steady)]


def package_budget(sim_ws, cbc_fn, kperkstp_skip=None, max_samples=120):
    """Mean rate [m3/d] of each budget term in a package .cbc file.

    Works for the UZF and SFR budget files (and the GWF cbc). Skips the first
    ``kperkstp_skip`` stress periods (the steady state; None = 1 when the
    model HAS a steady first period, 0 when it does not) and averages over an
    even subsample of at most ``max_samples`` of the rest (set None for all).
    """
    import flopy
    if kperkstp_skip is None:
        kperkstp_skip = 1 if steady_first(sim_ws) else 0
    cbc = flopy.utils.CellBudgetFile(os.path.join(sim_ws, cbc_fn))
    records = [r.strip().decode() if isinstance(r, bytes) else str(r).strip()
               for r in cbc.get_unique_record_names()]
    kk = cbc.get_kstpkper()
    # each record weighted by its STEP LENGTH: under ATS a hard day leaves
    # dozens of short steps, and counting each as much as a whole day pulled
    # the "mean rate" toward those days
    dt = dict(zip(kk, np.diff(np.r_[0.0, np.asarray(cbc.get_times(),
                                                   dtype=float)])))
    keep = _subsample([k for k in kk if k[1] >= kperkstp_skip] or kk, max_samples)
    out = {}
    for rec in records:
        tot = 0.0
        n = 0.0
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
            w = float(dt.get(k, 1.0)) or 1.0
            tot += float(val) * w
            n += w
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
    # THE END OF EVERY STRESS PERIOD, not every saved step: ATS saves one
    # per time step (period_steps)
    kk_real = period_end_kk(hds, steady=steady_first(sim_ws, name))
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
    for pkg, cbc in (('uzf', '%s.uzf.cbc' % name), ('sfr', '%s.sfr.cbc' % name),
                     ('lak', '%s.lak.cbc' % name)):
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

    # --- the stream at the outlet, and its budget over time (WP3.6) -- #
    try:
        written += _fig_stream(sim_ws, name, ds_ws, out, mg, dates=dates,
                               verbose=verbose)
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   streamflow figures skipped: %r' % exc)

    # --- the ponds: stages and the LAK budget over time (WP4.8) ------- #
    try:
        written += _fig_lakes(sim_ws, name, ds_ws, out, dates=dates,
                              verbose=verbose)
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   lake figures skipped: %r' % exc)

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


# --------------------------------------------------------------------- #
# the stream (WP3.6)
# --------------------------------------------------------------------- #

def _tdis_perlen(sim_ws):
    """Stress-period lengths [d] from the simulation's TDIS file, or None."""
    import glob
    import re
    fn = None
    nam = os.path.join(sim_ws, 'mfsim.nam')
    if os.path.exists(nam):
        with open(nam, encoding='utf-8', errors='replace') as fh:
            m = re.search(r'^\s*TDIS6\s+(\S+)', fh.read(), re.I | re.M)
        if m:
            fn = os.path.join(sim_ws, m.group(1))
    if fn is None or not os.path.exists(fn):
        fn = next(iter(sorted(glob.glob(os.path.join(sim_ws, '*.tdis')))),
                  None)
    if fn is None:
        return None
    with open(fn, encoding='utf-8', errors='replace') as fh:
        m = re.search(r'BEGIN\s+PERIODDATA(.*?)END\s+PERIODDATA', fh.read(),
                      re.I | re.S)
    if m is None:
        return None
    per = []
    for line in m.group(1).splitlines():
        p = line.split('#')[0].split()
        if not p:
            continue
        try:
            per.append(float(p[0]))
        except ValueError:                 # OPEN/CLOSE: not read here
            return None
    return np.asarray(per) if per else None


def sfr_observations(sim_ws, name, perlen=None):
    """The SFR continuous observations as ONE ROW PER STRESS PERIOD.

    MF6 writes <name>.obs.sfr.csv at every time step, and adaptive time
    stepping cuts a period into a parameter-dependent number of them. Each
    row here is the time-weighted mean of a period's steps -- the rate that
    moves the period's volume -- so there are as many rows as periods
    whatever the stepping (the CdL trap, cookbook WP6.4). The steady first
    period, when there is one, is dropped. Columns: the observation names
    in lower case, and ``perlen`` [d]. None when the run wrote no file.
    """
    import pandas as pd
    fn = os.path.join(sim_ws, '%s.obs.sfr.csv' % name)
    if not os.path.exists(fn):
        return None
    df = pd.read_csv(fn)
    df.columns = [str(c).strip().lower() for c in df.columns]
    t = df['time'].to_numpy(dtype=float)
    dt = np.diff(np.r_[0.0, t])
    if perlen is None:
        perlen = _tdis_perlen(sim_ws)
    if perlen is None:                     # no TDIS to read: a row per step
        perlen = dt
    perlen = np.asarray(perlen, dtype=float)
    ends = np.cumsum(perlen)
    # a step belongs to the period it ENDS in; the tolerance absorbs the
    # rounding of the accumulated time, far below any ATS minimum step
    kper = np.searchsorted(ends + 1e-6, t, side='left')
    keep = kper < len(ends)
    k, w = kper[keep], dt[keep]
    wsum = np.bincount(k, weights=w, minlength=len(ends))
    out = {}
    for c in df.columns:
        if c == 'time':
            continue
        v = df[c].to_numpy(dtype=float)[keep]
        with np.errstate(invalid='ignore', divide='ignore'):
            out[c] = (np.bincount(k, weights=v * w, minlength=len(ends))
                      / np.where(wsum > 0, wsum, np.nan))
    res = pd.DataFrame(out)
    res['perlen'] = perlen
    if len(res) > 1 and steady_first(sim_ws, name):
        res = res.iloc[1:].reset_index(drop=True)
    return res


def obs_streamflow(ds_ws, fn=None):
    """Observed catchment streamflow [mm/d] as a date-indexed Series, or None.

    ``<ro_prefix>_catchment.txt`` (date, value), read as the legacy driver
    reads it: the gauge discharge over the catchment area, which it
    compared with the catchment-mean MM runoff in mm/d.
    """
    import pandas as pd
    fn = fn or os.path.join(ds_ws, '%s_catchment.txt' % OBS['ro'])
    if not os.path.exists(fn):
        return None
    df = pd.read_csv(fn, sep=r'\s+', header=None, usecols=[0, 1],
                     names=['date', 'q'], engine='python', comment='#')
    df['date'] = pd.to_datetime(df['date'], errors='coerce')
    df['q'] = pd.to_numeric(df['q'], errors='coerce')
    df = df.dropna()
    df = df[df['q'] > -9000.0]             # hnoflo marks a gap
    if df.empty:
        return None
    return df.set_index('date')['q'].sort_index()


def _active_area(mg):
    """The model's active footprint [m2]: the catchment the outlet drains."""
    idm = getattr(mg, 'idomain', None)
    if getattr(mg, 'nrow', None) is None:
        verts = np.asarray(mg.verts, dtype=float)
        area = []
        for iv in mg.iverts:
            xy = verts[list(iv)]
            # relative to the first vertex: on absolute UTM coordinates the
            # shoelace cancels catastrophically (2026-09-24)
            x, y = (xy - xy[0]).T
            area.append(0.5 * abs(np.dot(x, np.roll(y, -1))
                                  - np.dot(y, np.roll(x, -1))))
        area = np.asarray(area)
    else:
        area = np.outer(np.asarray(mg.delc, dtype=float),
                        np.asarray(mg.delr, dtype=float)).ravel()
    if idm is None:
        return float(area.sum())
    act = (np.asarray(idm).reshape(int(mg.nlay), -1) > 0).any(axis=0)
    return float(area[act].sum())


def _period_mean(series, starts, perlen):
    """Mean of a date-indexed series over each period [start, start+perlen)."""
    ot = series.index.values
    ov = series.to_numpy(dtype=float)
    t0 = np.asarray(starts, dtype='datetime64[ns]')
    t1 = t0 + (np.asarray(perlen, dtype=float) * 86400e9).astype(
        'timedelta64[ns]')
    a = np.searchsorted(ot, t0, side='left')
    b = np.searchsorted(ot, t1, side='left')
    csum = np.r_[0.0, np.cumsum(ov)]
    cnt = b - a
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(cnt > 0, (csum[b] - csum[a]) / np.maximum(cnt, 1),
                        np.nan)


def _fit(sim, obs):
    """NSE, r, volume bias [%] and n over the periods both series exist."""
    sim, obs = np.asarray(sim, dtype=float), np.asarray(obs, dtype=float)
    m = np.isfinite(sim) & np.isfinite(obs)
    s, o = sim[m], obs[m]
    if len(o) < 2 or np.ptp(o) == 0:
        return None
    return {'nse': 1.0 - np.sum((s - o) ** 2) / np.sum((o - o.mean()) ** 2),
            'r': float(np.corrcoef(s, o)[0, 1]),
            'bias': (100.0 * (s.sum() - o.sum()) / o.sum() if o.sum()
                     else np.nan),
            'n': int(len(o)), 'mask': m}


# the network budget from the stream's point of view, + into it. MF6 reports
# ext-inflow, from-mvr and evaporation positive, to-mvr and ext-outflow
# negative -- and the aquifer exchange POSITIVE WHEN THE STREAM LOSES: the
# reach solver negates the leakage, qgwf = cond * (h_sfr - h_gw), and routes
# qd = qsrc - qgwf downstream (gwf-sfr-steady.f90). Checked on a five-reach
# model: 500 in - 0.1 evaporated - 1.633 'sfr' = 498.267 out.
_SFR_TERMS = (('net_inflow', 1.0, 'MM runoff in', 'tab:blue'),
              ('net_from_mvr', 1.0, 'from the ponds', 'tab:cyan'),
              ('net_leakage', -1.0, 'from the aquifer (net)', 'tab:brown'),
              ('net_to_mvr', 1.0, 'to the ponds', 'tab:olive'),
              ('net_evaporation', -1.0, 'open-water evaporation', 'tab:red'),
              ('net_outflow', 1.0, 'outflow at the outlet', 'tab:purple'))


def _fig_stream(sim_ws, name, ds_ws, out, mg, dates=None, verbose=True):
    """Outlet hydrograph against the gauge, and the SFR budget over time.

    Both from the SFR continuous observations -- one small csv, exact at
    every step -- rather than the SFR budget file. Nothing is drawn for a
    run without SFR.
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    per = sfr_observations(sim_ws, name)
    if per is None or 'outflow' not in per.columns:
        return []
    written = []
    n = len(per)
    area = _active_area(mg)
    if dates is None or len(dates) < n:
        dates = _dates_from_dataset(ds_ws, n)
    is_dt = isinstance(dates, pd.DatetimeIndex)
    x = pd.DatetimeIndex(dates[:n]) if is_dt else np.arange(n)
    mmyr = 1000.0 * 365.25 / area          # m3/d -> mm/yr over the catchment
    q = -per['outflow'].to_numpy(dtype=float)       # leaving: > 0
    tab = pd.DataFrame({'q_sim_m3d': q, 'q_sim_mmd': q * 1000.0 / area},
                       index=x)
    if 'stage' in per.columns:
        tab['stage_m'] = per['stage'].to_numpy()
    if 'leakage' in per.columns:
        # MF6's sign: > 0 when the outlet reach loses water to the aquifer
        tab['to_aquifer_m3d'] = per['leakage'].to_numpy()
    obs = obs_streamflow(ds_ws) if is_dt else None
    fit = None
    if obs is not None:
        tab['q_obs_mmd'] = _period_mean(obs, x, per['perlen'])
        tab['q_obs_m3d'] = tab['q_obs_mmd'] * area / 1000.0
        fit = _fit(tab['q_sim_m3d'], tab['q_obs_m3d'])
    fn = os.path.join(out, 'outlet_streamflow.csv')
    tab.to_csv(fn, index_label='date' if is_dt else 'period')
    written.append(fn)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True)
    for ax, log in zip(axes, (False, True)):
        ax.plot(x, q, color='tab:blue', lw=0.9, label='simulated, SFR outlet')
        if obs is not None:
            ax.plot(x, tab['q_obs_m3d'], ls='none', marker='.', ms=2.5,
                    color='k', label='observed, %s_catchment' % OBS['ro'])
        ax.set_ylabel('streamflow [m$^3$/d]')
        sec = ax.secondary_yaxis('right', functions=(lambda v: v * mmyr,
                                                     lambda v: v / mmyr))
        sec.set_ylabel('[mm/yr]')
        ax.grid(alpha=0.3)
        if log:
            pos = q[np.isfinite(q) & (q > 0)]
            if 'q_obs_m3d' in tab:
                o = tab['q_obs_m3d'].to_numpy()
                pos = np.r_[pos, o[np.isfinite(o) & (o > 0)]]
            if len(pos):
                ax.set_yscale('log')
                ax.set_ylim(max(pos.min(), pos.max() * 1e-5) * 0.8,
                            pos.max() * 1.5)
            ax.set_title('the same, log scale: baseflow and recessions',
                         fontsize=9)
    axes[0].legend(loc='upper right', fontsize=8)
    head = ('Streamflow at the catchment outlet (catchment %.2f km$^2$); '
            'simulated mean %.0f mm/yr' % (area / 1e6, np.nanmean(q) * mmyr))
    if fit is not None:
        o = tab['q_obs_m3d'].to_numpy()[fit['mask']]
        s = q[fit['mask']]
        head += ('\nover the %d period(s) observed: NSE %.2f, r %.2f, '
                 'volume bias %+.0f %%, simulated %.0f vs observed %.0f mm/yr'
                 % (fit['n'], fit['nse'], fit['r'], fit['bias'],
                    s.mean() * mmyr, o.mean() * mmyr))
    elif obs is not None:
        head += '\nno observed day falls inside the run'
    axes[0].set_title(head, fontsize=10)
    fn = os.path.join(out, 'outlet_streamflow.png')
    fig.savefig(fn, dpi=140, bbox_inches='tight')
    plt.close(fig)
    written.append(fn)

    # --- the network budget over time, in mm over the catchment -------- #
    terms = [t for t in _SFR_TERMS if t[0] in per.columns]
    if terms:
        plen = per['perlen'].to_numpy(dtype=float)
        vol = pd.DataFrame({lab: sgn * per[c].to_numpy() * plen * 1000.0 / area
                            for c, sgn, lab, _col in terms}, index=x)
        by, unit = vol, 'per stress period'
        if is_dt and n > 62:
            by = vol.groupby(x.to_period('M')).sum()
            by.index = by.index.to_timestamp()
            unit = 'per month'
        years = float(plen.sum()) / 365.25
        fig, ax = plt.subplots(figsize=(11, 5))
        xb = np.arange(len(by))
        pos = np.zeros(len(by))
        neg = np.zeros(len(by))
        for (_c, _s, lab, col) in terms:
            v = by[lab].to_numpy()
            tot = vol[lab].sum() / years if years > 0 else np.nan
            up, dn = np.where(v > 0, v, 0.0), np.where(v < 0, v, 0.0)
            ax.bar(xb, up, bottom=pos, color=col, width=0.85,
                   label='%s  (%+.1f mm/yr)' % (lab, tot))
            ax.bar(xb, dn, bottom=neg, color=col, width=0.85)
            pos += up
            neg += dn
        ax.plot(xb, by.sum(axis=1).to_numpy(), color='k', lw=0.8, marker='.',
                ms=3, label='closure (sum of the terms)')
        ax.axhline(0, color='k', lw=0.5)
        ax.set_ylabel('mm over the catchment, %s' % unit)
        if is_dt:
            ix = np.unique(np.linspace(0, len(by) - 1,
                                       min(len(by), 12)).astype(int))
            ax.set_xticks(ix)
            ax.set_xticklabels([by.index[i].strftime('%Y-%m') for i in ix],
                               rotation=45, ha='right', fontsize=8)
        ax.set_title('Stream network budget (SFR): + into the stream, '
                     '- out of it', fontsize=10)
        ax.legend(fontsize=8, loc='best')
        ax.grid(alpha=0.3, axis='y')
        fn = os.path.join(out, 'budget_sfr_ts.png')
        fig.savefig(fn, dpi=140, bbox_inches='tight')
        plt.close(fig)
        written.append(fn)
        rates = pd.DataFrame({c: per[c].to_numpy() for c, *_r in terms},
                             index=x)
        fn = os.path.join(out, 'budget_sfr_period.csv')
        rates.to_csv(fn, index_label='date' if is_dt else 'period')
        written.append(fn)
    if verbose:
        msg = '   outlet streamflow: mean %.0f m3/d (%.0f mm/yr)' % (
            np.nanmean(q), np.nanmean(q) * mmyr)
        if fit is not None:
            msg += '; against the gauge NSE %.2f, bias %+.0f %% (%d periods)' \
                % (fit['nse'], fit['bias'], fit['n'])
        print(msg)
    return written


def _lake_table_file(fn):
    """(bed, rim) of one LAK stage table (``marmites_lak.lake_table``: the
    first row is the bed, the last the rim plus a headroom row)."""
    rows, inside = [], False
    with open(fn) as fh:
        for line in fh:
            t = line.strip().upper()
            if t.startswith('BEGIN TABLE'):
                inside = True
            elif t.startswith('END TABLE'):
                break
            elif inside and t and not t.startswith('#'):
                rows.append(float(line.split()[0]))
    if not rows:
        raise ValueError('no table rows in %s' % fn)
    return rows[0], rows[-2] if len(rows) > 1 else rows[0]


def _lake_names(sim_ws, name, n):
    """The boundnames of the LAK packagedata, or lake 1..n."""
    names = []
    fn = os.path.join(sim_ws, '%s.lak' % name)
    if os.path.exists(fn):
        inside = False
        with open(fn) as fh:
            for line in fh:
                t = line.strip().upper()
                if t.startswith('BEGIN PACKAGEDATA'):
                    inside = True
                elif t.startswith('END PACKAGEDATA'):
                    break
                elif inside and t and not t.startswith('#'):
                    parts = line.split()
                    names.append(parts[3].strip("'\"") if len(parts) > 3
                                 else 'lake %s' % parts[0])
    return names if len(names) == n else ['lake %d' % (k + 1) for k in range(n)]


def lake_series(sim_ws, name):
    """Each lake's stage at the end of every stress period, its bed and rim,
    and the LAK budget per period [m3/d, + into the lakes, sub-step means].
    None for a run without LAK. The steady period, if any, is dropped."""
    import glob
    import re
    import flopy
    import pandas as pd
    sfn = os.path.join(sim_ws, '%s.lak.stage' % name)
    if not os.path.exists(sfn):
        return None
    st = flopy.utils.HeadFile(sfn, text='STAGE', precision='double')
    last = {}
    for ks, kp in st.get_kstpkper():
        last[kp] = ks
    pers = sorted(last)
    if steady_first(sim_ws, name) and len(pers) > 1:
        pers = pers[1:]
    stage = np.array([np.asarray(st.get_data(kstpkper=(last[p], p)),
                                 float).ravel() for p in pers])
    # MF6 writes a DRY lake's stage as its 1e30 no-data value
    stage = np.where(np.abs(stage) > 1e29, np.nan, stage)
    n = stage.shape[1]
    tabs = sorted(glob.glob(os.path.join(sim_ws, '%s.lak*.tab' % name)),
                  key=lambda f: int(re.findall(r'lak(\d+)\.tab$', f)[0]))
    geo = [_lake_table_file(f) for f in tabs[:n]]
    bed = np.array([g[0] for g in geo]) if len(geo) == n else None
    rim = np.array([g[1] for g in geo]) if len(geo) == n else None
    budget = None
    cfn = os.path.join(sim_ws, '%s.lak.cbc' % name)
    if os.path.exists(cfn):
        cb = flopy.utils.CellBudgetFile(cfn, precision='double')
        kk = cb.get_kstpkper()
        t = np.asarray(cb.get_times(), float)
        dt = np.diff(np.concatenate([[0.0], t]))
        terms = [r.decode().strip() if isinstance(r, bytes) else r.strip()
                 for r in cb.get_unique_record_names()]
        terms = [r for r in terms if r != 'FLOW-JA-FACE']
        rate = {r: np.zeros(len(pers)) for r in terms}
        span = np.zeros(len(pers))
        where = {p: k for k, p in enumerate(pers)}
        for i, (ks, kp) in enumerate(kk):
            if kp not in where:
                continue
            k = where[kp]
            span[k] += dt[i]
            for r in terms:
                d = cb.get_data(kstpkper=(ks, kp), text=r)
                if d:
                    rate[r][k] += float(np.sum(d[0]['q'])) * dt[i]
        span = np.where(span > 0, span, 1.0)
        budget = pd.DataFrame({r: v / span for r, v in rate.items()})
    return {'periods': pers, 'stage': stage, 'bed': bed, 'rim': rim,
            'names': _lake_names(sim_ws, name, n), 'budget': budget}


def _fig_lakes(sim_ws, name, ds_ws, out, dates=None, verbose=True):
    """The ponds (WP4.8): each lake's stage against its bed and rim, and
    the LAK budget over time. Nothing for a run without LAK."""
    import matplotlib.pyplot as plt
    import pandas as pd
    ls = lake_series(sim_ws, name)
    if ls is None:
        return []
    written = []
    stage, n = ls['stage'], ls['stage'].shape[1]
    nper = stage.shape[0]
    if dates is None or len(dates) < nper:
        dates = _dates_from_dataset(ds_ws, nper)
    is_dt = isinstance(dates, pd.DatetimeIndex)
    x = pd.DatetimeIndex(dates[:nper]) if is_dt else np.arange(nper)
    tab = pd.DataFrame(stage, index=x, columns=ls['names'])
    fn = os.path.join(out, 'lake_stage.csv')
    tab.to_csv(fn, index_label='date' if is_dt else 'period')
    written.append(fn)

    ncol = min(4, n)
    nrow = int(np.ceil(n / float(ncol)))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.4 * nrow),
                             sharex=True, squeeze=False)
    for L in range(nrow * ncol):
        ax = axes.flat[L]
        if L >= n:
            ax.set_visible(False)
            continue
        ax.plot(x, stage[:, L], color='tab:blue', lw=1.0)
        title = ls['names'][L]
        if ls['bed'] is not None:
            ax.axhline(ls['bed'][L], color='saddlebrown', lw=0.8, ls='--')
            ax.axhline(ls['rim'][L], color='k', lw=0.8, ls=':')
            d = np.nan_to_num(stage[:, L] - ls['bed'][L], nan=0.0)
            title += ': %.2f-%.2f m deep, rim %.2f m' % (
                d.min(), d.max(), ls['rim'][L] - ls['bed'][L])
            if np.isnan(stage[:, L]).any():
                title += ', dry %d period(s)' % int(np.isnan(stage[:, L]).sum())
        ax.set_title(title, fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
        if is_dt:
            for lab in ax.get_xticklabels():
                lab.set_rotation(45)
                lab.set_ha('right')
    fig.suptitle('Lake stages (LAK) [m]: dashed the bed, dotted the rim (the '
                 'outlet sill)', fontsize=10)
    fig.tight_layout()
    fn = os.path.join(out, 'lake_stage.png')
    fig.savefig(fn, dpi=140, bbox_inches='tight')
    plt.close(fig)
    written.append(fn)

    bud = ls['budget']
    if bud is not None and len(bud):
        bud.index = x
        fn = os.path.join(out, 'budget_lak_period.csv')
        bud.to_csv(fn, index_label='date' if is_dt else 'period')
        written.append(fn)
        by, unit = bud, 'm$^3$/d, per stress period'
        if is_dt and nper > 62:
            by = bud.groupby(x.to_period('M')).mean()
            by.index = by.index.to_timestamp()
            unit = 'm$^3$/d, monthly mean'
        cmap = plt.get_cmap('tab10')
        fig, ax = plt.subplots(figsize=(11, 5))
        xb = np.arange(len(by))
        pos = np.zeros(len(by))
        neg = np.zeros(len(by))
        for c, term in enumerate(by.columns):
            v = by[term].to_numpy()
            up, dn = np.where(v > 0, v, 0.0), np.where(v < 0, v, 0.0)
            ax.bar(xb, up, bottom=pos, width=0.85, color=cmap(c % 10),
                   label='%s  (%+.1f m$^3$/d)' % (term, bud[term].mean()))
            ax.bar(xb, dn, bottom=neg, width=0.85, color=cmap(c % 10))
            pos += up
            neg += dn
        ax.plot(xb, by.sum(axis=1).to_numpy(), color='k', lw=0.8, marker='.',
                ms=3, label='closure (sum of the terms)')
        ax.axhline(0, color='k', lw=0.5)
        ax.set_ylabel(unit)
        if is_dt:
            ix = np.unique(np.linspace(0, len(by) - 1,
                                       min(len(by), 12)).astype(int))
            ax.set_xticks(ix)
            ax.set_xticklabels([by.index[i].strftime('%Y-%m') for i in ix],
                               rotation=45, ha='right', fontsize=8)
        ax.set_title('Lake budget (LAK, all %d lakes): + into the lakes, - out '
                     'of them (STORAGE - = the lakes filling)' % n, fontsize=10)
        ax.legend(fontsize=8, loc='best')
        ax.grid(alpha=0.3, axis='y')
        fn = os.path.join(out, 'budget_lak_ts.png')
        fig.savefig(fn, dpi=140, bbox_inches='tight')
        plt.close(fig)
        written.append(fn)
    if verbose:
        msg = '   lakes: %d' % n
        if ls['bed'] is not None:
            d = np.nan_to_num(stage - ls['bed'][None, :], nan=0.0)
            msg += ('; depth above the bed %.2f..%.2f m, %d of them dry at '
                    'some point' % (d.min(), d.max(),
                                    int(np.sum(d.min(axis=0) < 0.01))))
        if bud is not None and len(bud):
            msg += '; budget (m3/d, + in): ' + ', '.join(
                '%s %+.1f' % (t, bud[t].mean()) for t in bud.columns)
        print(msg)
    return written


def run_preproc(sim_ws, ds_ws, name='lamatamm', mf_ws=None, verbose=True,
                out_root=None, cMF=None, ctx=None, res=None, trunk=None,
                gis_ws=None, input_maps=True):
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

    # --- the model map: grid, streams, ponds, boundary packages ------- #
    try:
        written += _fig_model_map(out, sim_ws, name, ds_ws,
                                  title=os.path.basename(os.path.normpath(ds_ws)),
                                  verbose=verbose)
    except Exception as exc:               # pragma: no cover
        if verbose:
            print('   model map skipped: %r' % exc)

    # --- the native parameter-field maps ------------------------------ #
    # [postproc] input_maps: drawn whatever the Plots panel said, until now.
    if not input_maps:
        if verbose:
            print('   input parameter maps: off on the Plots panel')
    elif cMF is not None and ctx is not None and MMplot is not None:
        try:
            written += _native_input_maps(MMplot, out, cMF, ctx, res=res,
                                          verbose=verbose)
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   native input maps skipped: %r' % exc)
    if verbose:
        print('preproc: %d file(s) written to %s' % (len(written), out))
    return written


# --------------------------------------------------------------------- #
# the model map: grid, catchment, streams, ponds, boundary packages
# --------------------------------------------------------------------- #

def _geojson_polygons(fn):
    """``([(properties, [exterior ring, ...]), ...], crs)`` from a GeoJSON file
    the converter wrote; rings as (n, 2) arrays. ``([], None)`` when absent."""
    import json
    if not fn or not os.path.exists(fn):
        return [], None
    with open(fn, encoding='utf-8') as fh:
        d = json.load(fh)
    out = []
    for f in d.get('features', []):
        g = f.get('geometry') or {}
        polys = ([g.get('coordinates', [])] if g.get('type') == 'Polygon'
                 else g.get('coordinates', []) if g.get('type') == 'MultiPolygon'
                 else [])
        rings = [np.asarray(p[0], dtype=float)[:, :2] for p in polys if p]
        out.append((f.get('properties') or {}, rings))
    return out, (d.get('marmites') or {}).get('crs')


def _ring_centroid(r):
    """Area centroid of a ring, taken relative to its first vertex (UTM)."""
    x0, y0 = r[0]
    x, y = r[:, 0] - x0, r[:, 1] - y0
    cr = x * np.roll(y, -1) - np.roll(x, -1) * y
    a = cr.sum() / 2.0
    if abs(a) < 1e-12:
        return float(r[:, 0].mean()), float(r[:, 1].mean())
    return (float(x0 + ((x + np.roll(x, -1)) * cr).sum() / (6.0 * a)),
            float(y0 + ((y + np.roll(y, -1)) * cr).sum() / (6.0 * a)))


def model_map_features(sim_ws, name=None, ds_ws=None):
    """What the model map draws, read from the WRITTEN simulation.

    Returns a dict: ``mg`` (the flopy grid), ``node`` (cellid -> cell
    number), and the cell numbers of the SFR reaches, the SFR outlets and
    inlets (reaches given a specified inflow), the LAK host cells, and the
    DRN (the outlet drain, not the seepage drains) and GHB cells. A reach
    that hands its flow to a pond through MVR ends there but is not an
    outlet. With ``ds_ws`` and a LAK package, ``lak_foot`` holds the pond
    footprints -- the cells the builder cut the stream out of, recomputed
    from the dataset's pond polygons by the same rule. None when ``sim_ws``
    holds no simulation.
    """
    import flopy
    if not os.path.exists(os.path.join(sim_ws, 'mfsim.nam')):
        return None
    sim = flopy.mf6.MFSimulation.load(
        sim_ws=sim_ws, verbosity_level=0,
        load_only=['dis', 'disv', 'drn', 'ghb', 'sfr', 'lak', 'mvr'])
    gwf = (sim.get_model(name) if name and name in sim.model_names
           else sim.get_model())
    mg = gwf.modelgrid
    vertex = mg.grid_type == 'vertex'
    ncol = 1 if vertex else int(mg.ncol)

    def node(cid):
        c = [int(v) for v in cid]
        return c[1] if vertex else c[1] * ncol + c[2]

    pk = {p.package_name.lower(): p for p in gwf.packagelist}
    f = {'mg': mg, 'node': node, 'sfr': [], 'outlet': [], 'inlet': [],
         'lak': [], 'lak_foot': [], 'drn': [], 'ghb': [], 'nreaches': 0,
         'nlakes': 0}
    to_mvr = set()
    if 'mvr' in pk:
        recs = pk['mvr'].perioddata.get_data(0)
        for r in (recs if recs is not None else []):
            if str(r['pname1']).lower() == 'sfr':
                to_mvr.add(int(r['id1']))
    for key in ('drn', 'ghb'):
        if key in pk:
            spd = pk[key].stress_period_data.get_data(0)
            if spd is not None and len(spd):
                f[key] = sorted({node(c) for c in spd['cellid']})
    if 'sfr' in pk:
        s = pk['sfr']
        pdat = s.packagedata.get_data()
        cell_of = {int(r['ifno']): node(r['cellid']) for r in pdat}
        f['sfr'] = sorted(set(cell_of.values()))
        f['nreaches'] = len(cell_of)
        cd = s.connectiondata.get_data()
        ics = [n for n in cd.dtype.names if n != 'ifno']
        for r in cd:
            v = np.array([r[n] for n in ics], dtype=float)
            v = v[np.isfinite(v)]
            # MF6 writes a downstream connection as a NEGATIVE reach number,
            # and reach 0 is a headwater, so -0 never occurs
            if not (v < 0).any() and int(r['ifno']) not in to_mvr:
                f['outlet'].append(cell_of[int(r['ifno'])])
        spd = s.perioddata.get_data(0) if s.perioddata.has_data() else None
        if spd is not None:
            for r in spd:
                if str(r['sfrsetting']).lower() == 'inflow':
                    try:
                        if float(r['sfrsetting_data']) > 0:
                            f['inlet'].append(cell_of[int(r['ifno'])])
                    except (TypeError, ValueError):
                        pass       # a time series: named, so not an inlet here
    if 'lak' in pk:
        cd = pk['lak'].connectiondata.get_data()
        if cd is not None and len(cd):
            f['lak'] = sorted({node(c) for c in cd['cellid']})
        f['nlakes'] = int(pk['lak'].nlakes.get_data() or 0)
        fn = os.path.join(ds_ws, 'inputPONDS.geojson') if ds_ws else None
        if fn and os.path.exists(fn):
            f['lak_foot'] = _pond_footprint_cells(mg, fn)
    return f


def _pond_footprint_cells(mg, fn):
    """Cell numbers of the pond footprints on a flopy grid, by the
    builder's own rule (marmites_lak.pond_footprints)."""
    import marmites_lak as LK
    import marmites_vector as mv
    ncell = int(np.asarray(mg.xcellcenters).size)
    polys = []
    for n in range(ncell):
        v = [tuple(map(float, xy[:2])) for xy in mg.get_cell_vertices(n)]
        if len(v) > 3 and v[0] == v[-1]:
            v = v[:-1]
        polys.append(v)
    shape = ((ncell, 1) if mg.grid_type == 'vertex'
             else (int(mg.nrow), int(mg.ncol)))
    grid = mv.TargetGrid(polys, shape)
    idm = getattr(mg, 'idomain', None)
    act = (None if idm is None else
           (np.asarray(idm).reshape(int(mg.nlay), -1) > 0).any(axis=0))
    ponds = LK.pond_footprints(LK.read_pond_polygons(fn), grid, active=act,
                               verbose=False)
    return sorted({i * shape[1] + j for p in ponds for i, j in p.cells})


def _fig_model_map(out, sim_ws, name, ds_ws, title=None, verbose=True):
    """IN_000_model_map.png -- the model on one page, in the CdL symbology.

    The grid, the catchment boundary (black), the observation points
    (magenta diamonds), the mapped stream network (blue), the ponds
    (outlined dark blue) and the cells the packages occupy: SFR light blue, LAK orange.
    Boundary packages follow cdl_gwf_model_fable_v2 §14b -- the SHAPE is the
    package (DRN square, SFR triangle, GHB circle) and the COLOUR the
    direction (red out, blue in). The SFR outlet's cell number and the pond
    IDs are written in orange. Everything is read from the written
    simulation and the dataset's own files, never from a shapefile.
    """
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PatchCollection
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.patches import Polygon as MplPoly
    f = model_map_features(sim_ws, name, ds_ws=ds_ws)
    if f is None:
        return []
    mg = f['mg']
    xc = np.asarray(mg.xcellcenters, dtype=float).ravel()
    yc = np.asarray(mg.ycellcenters, dtype=float).ravel()
    ncell = len(xc)

    def patches(nodes):
        return [MplPoly(np.asarray(mg.get_cell_vertices(int(n)))[:, :2],
                        closed=True) for n in nodes]

    x0, x1, y0, y1 = mg.extent
    h = float(np.clip(10.0 * (y1 - y0) / max(x1 - x0, 1.0), 6.0, 14.0))
    fig, ax = plt.subplots(figsize=(10, h))
    ax.add_collection(PatchCollection(patches(range(ncell)), facecolor='none',
                                      edgecolor='0.8', lw=0.2, zorder=1))
    handles, labels = [], []

    def key(h_, lab):
        handles.append(h_)
        labels.append(lab)

    # the catchment boundary
    bnd, crs = _geojson_polygons(os.path.join(ds_ws, 'inputBOUNDARY.geojson'))
    for _p, rings in bnd:
        for r in rings:
            ax.plot(r[:, 0], r[:, 1], color='k', lw=1.2, zorder=5)
    if bnd:
        key(Line2D([0], [0], color='k', lw=1.2), 'catchment boundary')

    # the mapped stream network
    fn = os.path.join(ds_ws, 'inputSTREAM.csv')
    if os.path.exists(fn):
        try:
            import marmites_channel as mch
            lines, _par = mch.read_stream_lines(fn)
            ax.add_collection(LineCollection(
                [np.asarray(v) for v in lines.values()], colors='tab:blue',
                lw=1.3, zorder=4))
            key(Line2D([0], [0], color='tab:blue', lw=1.3), 'stream network')
        except Exception as exc:           # pragma: no cover
            if verbose:
                print('   model map: stream network skipped: %r' % exc)

    halo = [pe.withStroke(linewidth=2.2, foreground='white')]
    if f['sfr']:
        # a stroked edge as well as the fill: a mesh refined along the
        # channel has reach cells a few metres wide, which a bare fill would
        # leave hidden under the stream line
        ax.add_collection(PatchCollection(patches(f['sfr']),
                                          facecolor='deepskyblue', alpha=0.45,
                                          edgecolor='deepskyblue', lw=2.5,
                                          zorder=2))
        key(Patch(facecolor='deepskyblue', alpha=0.45), 'SFR reach cells')
    if f['lak_foot']:
        ax.add_collection(PatchCollection(patches(f['lak_foot']),
                                          facecolor='orange', alpha=0.75,
                                          edgecolor='none', zorder=3))
        key(Patch(facecolor='orange', alpha=0.75), 'LAK cells (pond footprint)')
    if f['lak']:
        # the one EMBEDDEDV connection of each lake
        ax.add_collection(PatchCollection(patches(f['lak']),
                                          facecolor='orange', alpha=0.9,
                                          edgecolor='crimson', lw=0.9,
                                          zorder=3.5))
        key(Patch(facecolor='orange', edgecolor='crimson'),
            'LAK connection (host cell)')

    # the ponds, OUTLINED in dark blue as in CdL (a fill would hide the
    # LAK cells beneath), and their IDs in orange
    ponds, _c = _geojson_polygons(os.path.join(ds_ws, 'inputPONDS.geojson'))
    for props, rings in ponds:
        for r in rings:
            ax.plot(r[:, 0], r[:, 1], color='navy', lw=1.0, zorder=4)
        if rings:
            pid = props.get('id', props.get('fid', ''))
            cx, cy = _ring_centroid(rings[0])
            ax.annotate('%s' % pid, (cx, cy), color='darkorange', fontsize=7.5,
                        fontweight='bold', ha='center', va='bottom',
                        xytext=(0, 6), textcoords='offset points', zorder=9,
                        path_effects=halo)
    if ponds:
        key(Line2D([0], [0], color='navy', lw=1.0), 'ponds (ID in orange)')

    # boundary packages: shape = package, colour = direction
    M_DRN, M_SFR, M_GHB = 's', '^', 'o'
    C_IN, C_OUT = 'royalblue', 'red'
    S_BAND, S_PT = 27, 82
    if f['ghb']:
        ax.scatter(xc[f['ghb']], yc[f['ghb']], marker=M_GHB, s=S_BAND, c=C_IN,
                   edgecolors='k', linewidths=0.4, zorder=6)
    if f['drn']:
        ax.scatter(xc[f['drn']], yc[f['drn']], marker=M_DRN, s=S_BAND,
                   c=C_OUT, edgecolors='k', linewidths=0.4, zorder=6)
    if f['inlet']:
        ax.scatter(xc[f['inlet']], yc[f['inlet']], marker=M_SFR, s=S_PT + 20,
                   c=C_IN, edgecolors='k', linewidths=0.9, zorder=8)
    for n in f['outlet']:
        ax.scatter(xc[n], yc[n], marker=M_SFR, s=S_PT + 20, c=C_OUT,
                   edgecolors='k', linewidths=0.9, zorder=8)
        ax.annotate('outlet: cell %d' % n, (xc[n], yc[n]), color='darkorange',
                    fontsize=8, fontweight='bold', xytext=(9, -14),
                    textcoords='offset points', zorder=9, path_effects=halo)

    def shape(m):
        return Line2D([0], [0], marker=m, color='w', markerfacecolor='none',
                      markeredgecolor='k', markeredgewidth=1.1, markersize=7,
                      linestyle='none')
    if f['sfr']:
        key(shape(M_SFR), 'SFR  (inflow / outlet)')
    if f['ghb']:
        key(shape(M_GHB), 'GHB  (inflow)')
    if f['drn']:
        key(shape(M_DRN), 'DRN  (outflow)')
    if f['sfr'] or f['ghb'] or f['drn']:
        key((Patch(facecolor=C_OUT, edgecolor='k'),
             Patch(facecolor=C_IN, edgecolor='k')),
            'red: outflow;  blue: inflow')

    # observation points
    try:
        pts = obs_points(ds_ws)
    except (OSError, ValueError):
        pts = []
    for p in pts:
        ax.scatter(p['x'], p['y'], marker='D', s=48, c='magenta',
                   edgecolors='k', linewidths=0.8, zorder=7)
        ax.annotate(p['name'], (p['x'], p['y']), textcoords='offset points',
                    xytext=(5, 4), fontsize=8, fontweight='bold', zorder=7,
                    path_effects=halo)
    if pts:
        key(Line2D([0], [0], marker='D', color='w', markerfacecolor='magenta',
                   markeredgecolor='k', markersize=6, linestyle='none'),
            'obs points')

    ax.legend(handles, labels, loc='upper right', fontsize=6, framealpha=0.95,
              labelspacing=0.35, handlelength=1.4, handleheight=1.0,
              handletextpad=0.5, borderpad=0.4,
              handler_map={tuple: HandlerTuple(ndivide=None)})
    kind = 'DISV' if mg.grid_type == 'vertex' else 'DIS'
    head = '%s model -- %d cells (%s)' % (title or name, ncell, kind)
    if f['nreaches']:
        head += ', %d SFR reaches' % f['nreaches']
    if ponds:
        head += ', %d ponds' % len(ponds)
        if f['nlakes']:
            head += ' (%d lakes)' % f['nlakes']
    ax.set_title(head)
    if f['outlet']:
        ax.text(0.01, 0.01, 'cell numbers from 0, as in the run log; the MF6 '
                'files count from 1', transform=ax.transAxes, fontsize=6.5,
                color='0.4', zorder=9)
    ax.ticklabel_format(useOffset=False, style='plain')
    ax.set_xlabel('X (m%s)' % (', ' + crs if crs else ''))
    ax.set_ylabel('Y (m)')
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect('equal')
    fn = os.path.join(out, 'IN_000_model_map.png')
    fig.savefig(fn, dpi=150, bbox_inches='tight')
    plt.close(fig)
    if verbose:
        print('   model map: %s' % fn)
    return [fn]


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
        # backdrop is rendered on the display raster -- the one map that
        # keeps a raster: a hillshade and DEM contours need a regular grid.
        DA = MapAdapter(cMF, raster=True)
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
        if OBS['name_column'] in gx:
            import matplotlib.patheffects as pe
            # P0 / SM / EC sit within ~100 m of each other, so the labels
            # are placed round the marker in turn instead of all up-right;
            # the halo keeps them readable over the relief
            box = ((7, 8, 'left'), (7, -14, 'left'),
                   (-7, 8, 'right'), (-7, -14, 'right'))
            for k, (nm, pt) in enumerate(zip(gx[OBS['name_column']],
                                                 gx.geometry)):
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
                 map_days=6, sankey_obs_years=False, obs_series=True,
                 result_maps=True):
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
    # [postproc] obs_series and result_maps: these were drawn whatever the
    # Plots panel said.
    if not obs_series:
        if verbose:
            print('   per-point time series: off on the Plots panel')
    else:
        try:
            written += _native_obs_timeseries(
                MMplot, out_dir, cMF, ctx, res, sim_ws, name,
                agg=agg if sankey else None, ds_ws=ds_ws, verbose=verbose)
        except Exception as exc:                     # pragma: no cover
            if verbose:
                print('   native obs time series skipped: %r' % exc)

    # --- result maps: the MM fluxes and the aquifer terms ------------- #
    if not result_maps:
        if verbose:
            print('   result maps: off on the Plots panel')
    else:
        try:
            written += _native_result_maps(MMplot, out_dir, cMF, ctx, res,
                                           sim_ws, name, ndays=map_days,
                                           verbose=verbose)
        except Exception as exc:                     # pragma: no cover
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
    # WP5, the CRR cascade -- drawn only when it ran (all zero otherwise)
    ('iRunon', 'Runon', 'run-on from upslope (CRR)', 'Blues'),
    ('iReinf', 'Reinf', 'reinfiltrated run-on (CRR)', 'Blues'),
    ('iEcrr', 'Ecrr', 'runoff evaporated by the cascade (CRR)', 'Reds'),
)
_CRR_MAPS = ('iRunon', 'iReinf', 'iEcrr')


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


def _obs4map(res, DA):
    """Observation points for plotLAYER's ``points`` overlay: [lbl, y, x,
    lay], already in plotLAYER's frame -- the map adapter places them (one
    cell off on the structured grid before, and on a mesh at x = 0, y =
    icell2d, which stretched every map into a strip)."""
    if res is None or 'obs_ij' not in res:
        return None
    ij = np.asarray(res['obs_ij'])
    names = [n.decode() if isinstance(n, bytes) else str(n)
             for n in res.get('obs_names', [str(k) for k in range(len(ij))])]
    return DA.points(names, [(int(i), int(j)) for i, j in ij])


def _hydro_year_index(DATE, ini_month):
    """Hydrological-year boundary indices, ported verbatim from the legacy
    driver (startMARMITES_v3.py ~449-483) so the Sankey aggregates exactly as
    the NWT post-processing did. Returns ``(HYindex, year_lst)``.

    ``DATE`` is the array of matplotlib date numbers (one per stress period).

    A RUN WITHOUT A FULL HYDROLOGICAL YEAR is drawn as a whole: ``HYindex =
    [0, 0, N-1, N-1]`` and no year of its own (``year_lst`` empty). The
    legacy branch started the panel at the run's first 1 October instead
    (``[h0, h0, N-1, N-1]``) and called it "the average of 1 hydrological
    year": La Mata's run of 31 May 2008 .. 30 May 2009 showed October-May
    scaled to a year -- P 401 mm/y where the run had 353, groundwater ET 7
    where it had 38 (2026-09-25).
    """
    import matplotlib as mpl
    import matplotlib.dates  # noqa: F401  -- a submodule: mpl.dates needs it
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
        print('   Sankey: the run does not contain a full hydrological year; '
              'the panel is the whole run, scaled to mm/y.')
        return [0, 0, len(DATE) - 1, len(DATE) - 1], []
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
    # ... or, on et.gw_route = 'evt', the two EVT packages (Eg, Tg apart)
    ('EVT_EG', 'EVT', 'EVT_EG'),
    ('EVT_TG', 'EVT', 'EVT_TG'),
    ('GHB', 'GHB', None),                  # head-dependent boundary
    # the streams and the ponds exchange with the aquifer directly, both
    # ways (+ into the aquifer): missing here, a losing stream's seepage and
    # a gaining stream's baseflow read as aquifer storage change
    ('SFR', 'SFR', None),
    ('LAK', 'LAK', None),
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
    # the RECORDS read are part of the key: a digest made before a record
    # was added would otherwise be reused without it
    # ... and the METHOD: digests made before the periods were read as
    # periods (period_steps, 2026-10-05) hold the last nper TIME steps
    sig = 'p2_%d_%d_%d_%d_%s' % (st.st_size, int(st.st_mtime), nper,
                                 len(targets),
                                 '+'.join(k for k, _t, _p in _AQ_RECORDS))

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
    # each stress period's records and their shares of it (period_steps):
    # the period's rate is the time-weighted mean of its steps' rates
    pgroups = period_steps(cbc.get_times(), kk,
                           steady=steady_first(sim_ws, name))[-nper:]
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
        def slab(rec):
            if flf_mode:
                fj = np.asarray(rec, float).ravel()
                a = np.zeros(nodes)
                a[src] = -fj[pos]                  # flopy negates; match it
                return a.reshape(nlay, -1)
            return np.ma.filled(np.asarray(rec, float), 0.0).reshape(nlay, -1)

        try:
            data = (cbc.get_data(text=text, paknam2=pak2) if flf_mode else
                    cbc.get_data(text=text, paknam2=pak2, full3D=True))
            if not data:
                return 'absent'
            if len(data) < len(kk):
                raise ValueError('%s: got %d records for %d time steps'
                                 % (key, len(data), len(kk)))
            # each period's rate: its steps' rates weighted by their length
            for k, grp in enumerate(pgroups):
                a = sum(w * slab(data[r]) for r, w in grp)
                _reduce(a, dest[k])
            del data
            return 'bulk'
        except MemoryError:
            if verbose:
                print('   %s: series too large for one read, '
                      'falling back to per-stress-period' % key)
        for k, grp in enumerate(pgroups):
            a = None
            for r, w in grp:
                d = (cbc.get_data(text=text, paknam2=pak2, kstpkper=kk[r])
                     if flf_mode else
                     cbc.get_data(text=text, paknam2=pak2, kstpkper=kk[r],
                                  full3D=True))
                if not d:
                    continue
                s_ = w * slab(d[0])
                a = s_ if a is None else a + s_
            if a is not None:
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
    flows); ``CH`` is 0 (no such package here). ``GHB``, and ``SFR`` and
    ``LAK`` (``iSFR_L``, ``iLAK_L``: the streams and the ponds), are signed,
    + into the aquifer -- seepage from them -- and - out of it, groundwater
    discharging into them.
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
    GHB = vol('GHB') * to_mm                                  # + in, - out
    SFR = vol('SFR') * to_mm                                  # + in, - out
    LAK = vol('LAK') * to_mm                                  # + in, - out
    wel = -vol('WEL') * to_mm                                 # >0 magnitude out
    # on et.gw_route = 'evt' MF6 itself keeps Eg and Tg apart; the WEL is
    # then zero but for a steady first period's mean
    Egl = wel * eg_frac[:, None] - vol('EVT_EG') * to_mm
    Tgl = wel * (1.0 - eg_frac[:, None]) - vol('EVT_TG') * to_mm
    # exfiltration to the soil: prefer an explicit seepage-drain package,
    # else the coupler's captured exfiltration, assigned to the top layer
    EXF = vol('DRN_SEEP') * to_mm                             # <0 out, or zeros
    if not np.any(EXF) and 'exf' in res:
        exf = np.asarray(res['exf'])
        if sel_ij is None:
            # the catchment: area-weighted, never a plain mean over cells
            _g = getattr(ctx, 'geom', None)
            _w = (np.asarray(_g.area, dtype=float) if _g is not None
                  else None)
            EXF[:, 0] = -(np.average(exf, axis=1, weights=_w)
                          if _w is not None and _w.size == exf.shape[1]
                          else exf.mean(axis=1))
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
        out['iGHB_%d' % (L + 1)] = GHB[:, L]
        out['iSFR_%d' % (L + 1)] = SFR[:, L]
        out['iLAK_%d' % (L + 1)] = LAK[:, L]
        out['iCH_%d' % (L + 1)] = np.zeros(nper)
    # UZF unsaturated storage change, from the UZF mass balance so it closes:
    #   percolation in = recharge to GW out + dS_unsat
    perc = (wb_ts[:, IX['iperc']] if sel_ij is None
            else np.zeros(nper))              # per-point perc filled by caller
    out['idSu'] = perc - Rg.sum(axis=1) - _uzf_out(wb_ts, IX, sel_ij is None)
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


def _ghb_cells(aq, nlay):
    """Per layer, does this target have a GHB flow at all? (the Sankey's
    ``ghbcells[L] > 0`` guard, as ``drncells`` is for the drains)"""
    return [int(np.any(np.asarray(aq.get('iGHB_%d' % (L + 1), 0.0)) != 0.0))
            for L in range(nlay)]


class _SankeyMF(object):
    """Thin wrapper over cMF supplying the plotting-metadata attributes the
    native Sankey reads, without mutating the real cMF."""
    def __init__(self, cMF, ncell_MM, drncells, dates, ghbcells=None):
        self._c = cMF
        nlay = int(cMF.nlay)
        self.wel_yn = 1               # groundwater ET drawn (as Eg/Tg)
        self.drn_yn = 1 if any(drncells) else 0
        self.drncells = drncells
        self.ghbcells = list(ghbcells) if ghbcells is not None else [0] * nlay
        self.ghb_yn = 1 if any(self.ghbcells) else 0
        self.ncell_MM = ncell_MM
        self.inputDate = dates

    def __getattr__(self, name):
        return getattr(self._c, name)


def _uzf_out(mmv, IX, use=True):
    """What leaves UZF other than recharge: its actual ET and the
    infiltration it rejected (WP2). The storage change was perc - Rg, so it
    silently absorbed both -- 58.7 mm/yr on 2026-09-24, 36 of it rejected
    infiltration."""
    if not use:
        return 0.0
    v = np.asarray(mmv)
    out = np.zeros(v.shape[0])
    # ... and what MF6's ET routine removed beyond the demand, within its
    # wave tolerance: not ET, but water gone from UZF all the same
    for k in ('iETuzf', 'iRejInf', 'iETuzf_num'):
        if k in IX:
            out = out + v[:, IX[k]]
    return out


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
        # a column the run's table does not have -- a run older than the
        # index -- is zero, as it was in that run
        if key in IX and IX[key] < np.shape(mmv)[1]:
            return mmv[:, IX[key]]
        return np.zeros(nper)

    put('iP', mm('iP')); put('iEi', mm('iEi')); put('iPe', mm('iPe'))
    # WP5: with the CRR cascade a cell's runoff is partly the run-on of the
    # cells above it, which another cell already counted as runoff. NET of
    # it, Ro is what leaves the soil surface for good -- to SFR, LAK, or the
    # cascade's evaporation -- and the surface balance Pe + Exf_1 = I + Ro
    # still closes (I includes the reinfiltration). Zero run-on: unchanged.
    put('idSsurf', mm('idSsurf')); put('iRo', mm('iRo') - mm('iRunon'))
    put('iEow', mm('iEow'))
    put('idSsoil', mm('idSsoil')); put('iEXFg', mm('iEXFg'))
    put('iI', mm('iI')); put('iSsurf', mm('iSsurf')); put('iperc', mm('iperc'))
    put('iETsoil', mm('iETsoil'))
    put('iEg', mm('iEg')); put('iTg', mm('iTg')); put('iETg', mm('iETg'))
    put('iEsoil', mmsv[:, :, IXS['iEsoil']].sum(axis=1))
    put('iTsoil', mmsv[:, :, IXS['iTsoil']].sum(axis=1))
    put('iExf_1', mmsv[:, 0, IXS['iExf']])                 # top soil layer
    # WP2: the deep unsaturated zone's actual ET, and the rejected
    # infiltration UZF returns to the soil column
    put('iETuzf', mm('iETuzf'))
    put('iRejInf', mm('iRejInf'))
    # WP2 row 2: the open water by package (iEow is their sum)
    put('iEow_sfr', mm('iEow_sfr'))
    put('iEow_lak', mm('iEow_lak'))
    # WP5: the cascade's terms, for the Sankey's CRR arms (WP6)
    put('iRunon', mm('iRunon'))
    put('iReinf', mm('iReinf'))
    put('iEcrr', mm('iEcrr'))
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
                            plot_years=plot_years,
                            # [postproc] wb_unit, set on cMF by
                            # props.apply_plot_settings
                            per_day=(str(getattr(smf, 'plt_WB_unit', 'year'))
                                     == 'day'))
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
    smf = _SankeyMF(cMF, ncell_MM, drncells, DATE,
                    ghbcells=_ghb_cells(aq, nlay))
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
        aq['idSu'] = (mmv[:, IX['iperc']] - sum(aq['iRg_%d' % (L + 1)]
                                                for L in range(nlay))
                      - _uzf_out(mmv, IX))
        flx, flxIndex = _assemble_flx(IX, IXS, mmv, mmsv, aq, nper)
        smf = _SankeyMF(cMF, ncell_MM, drncells, DATE,
                        ghbcells=_ghb_cells(aq, nlay))
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
    # the END of every stress period (period_steps): the last nper saved
    # steps were the last nper TIME steps, not periods, under ATS
    rows = [g[-1][0] for g in period_steps(
        hds.get_times(), hds.get_kstpkper(),
        steady=steady_first(sim_ws, name))][-nper:]
    out = np.full((len(cells_ij), nlay, nper), np.nan)
    vertex = _is_vertex(hds)
    for p, (i, j) in enumerate(cells_ij):
        ts = hds.get_ts([_cellid(vertex, L, i, j) for L in range(nlay)])
        for L in range(nlay):
            out[p, L, :len(rows)] = np.asarray(ts[rows, L + 1], float)
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
            inputObs_fn=OBS['table'], inputObsHEADS_fn=OBS['heads'],
            inputObsSM_fn=OBS['sm'], inputObsRo_fn=OBS['ro'],
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
            put('idSu', mm_obs[:, p, IX['iperc']] - rg.sum(axis=1)
                - _uzf_out(mm_obs[:, p, :], IX), r'$\Delta S_u$')
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
                date_ini=DATE[HYindex[1]], date_end=DATE[HYindex[-2]],
                **_tick_years(cMF))
            if os.path.exists(fn):
                written.append(fn)
            # groundwater-flux figure at the same point, from the same flux list
            fng = os.path.join(out_dir, '_0%s_tsGW.png' % o)
            try:
                MMplot.plotTIMESERIES_flxGW(
                    cMFd, flx, lbl, idx, fng,
                    'Groundwater fluxes at observation point %s' % o,
                    iniMonthHydroYear=int(getattr(cMF, 'iniMonthHydroYear', 10)),
                    date_ini=DATE[HYindex[1]], date_end=DATE[HYindex[-2]],
                    **_tick_years(cMF))
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
    sig = 'map_p2_%d_%d_%d' % (st.st_size, int(st.st_mtime), nper)
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
    # the time mean over the periods, each period its steps' time-weighted
    # mean (period_steps) -- not the mean of the last nper saved steps
    pgroups = period_steps(cbc.get_times(), kk,
                           steady=steady_first(sim_ws, name))[-nper:]
    np_ = max(len(pgroups), 1)
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
            if not data or len(data) < len(kk):
                continue
            for grp in pgroups:
                for r, w in grp:
                    acc += w * np.ma.filled(np.asarray(data[r], float),
                                            0.0).reshape(nlay, nrow, ncol)
                n += 1
            del data
        except MemoryError:                          # pragma: no cover
            for grp in pgroups:
                for r, w in grp:
                    d = cbc.get_data(text=text, paknam2=pak2,
                                     kstpkper=kk[r], full3D=True)
                    if d:
                        acc += w * np.ma.filled(np.asarray(d[0], float),
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
            for grp in pgroups:
                for r, w in grp:
                    fj = np.asarray(data[r], float).ravel()
                    a = np.zeros(nodes)
                    a[src] = -fj[pos]                # flopy's sign convention
                    acc += w * a
            del data
            out['FLF'] = (acc / np_).reshape(nlay, nrow, ncol)
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
    ('EVT_EG', 'Eg', 'groundwater evaporation (EVT)', -1.0),
    ('EVT_TG', 'Tg', 'groundwater transpiration (EVT)', -1.0),
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
    pts = _obs4map(res, DA)
    geo = DA.plot_geometry()      # the DISPLAY grid's axes
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
                hnoflo=hnoflo, cMF=geo, polys=DA.polys)
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
        if key not in IX or IX[key] >= wb_map.shape[1]:
            continue
        if key in _CRR_MAPS and not np.any(wb_map[:, IX[key]]):
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
    # the END of every stress period (period_steps), not the last nper saved
    # steps -- ATS saves one per time step
    kk = period_end_kk(hds, steady=steady_first(sim_ws, name))[-nper:]
    nper_h = len(kk)
    acc = None
    for k in range(nper_h):
        h = np.asarray(hds.get_data(kstpkper=kk[k]), dtype=float)
        h = np.where(np.abs(h) > 1e29, np.nan, h)
        acc = h if acc is None else acc + h
    draw(DA.lay(np.asarray(acc).reshape(nlay, -1) / nper_h)[None, :, :, :],
         'head', 'mean head', 'm', m3)
    if ndays and nper_h > 1:
        sel = np.unique(np.linspace(0, nper_h - 1, int(ndays)).astype(int))
        V = np.array([DA.lay(np.where(
            np.abs(np.asarray(hds.get_data(kstpkper=kk[s]), dtype=float))
            > 1e29, np.nan,
            np.asarray(hds.get_data(kstpkper=kk[s]), dtype=float)
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
        if 'WEL' in maps or 'EVT_EG' in maps:
            et = sum(np.asarray(maps[k], float) for k in
                     ('WEL', 'EVT_EG', 'EVT_TG') if k in maps)
            rn = re + DA.lay(np.asarray(et, float).reshape(nlay, -1))
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
    pts = _obs4map(res, DA)
    geo = DA.plot_geometry()      # the DISPLAY grid's axes

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
                          ):
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
                hnoflo=hnoflo, cMF=geo, polys=DA.polys)
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
            inputObs_fn=OBS['table'], inputObsHEADS_fn=OBS['heads'],
            inputObsSM_fn=OBS['sm'], inputObsRo_fn=OBS['ro'],
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
        # get_ts returns every SAVED STEP (and the steady period): keep the
        # end of each stress period, as kk_real lists them
        _pos = {k: r for r, k in enumerate(hds.get_kstpkper())}
        comp = ts[[_pos[k] for k in kk_real], 1]
        comp = np.where(np.abs(comp) > 1e29, np.nan, comp)
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
