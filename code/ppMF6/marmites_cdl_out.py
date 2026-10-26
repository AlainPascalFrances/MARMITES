# -*- coding: utf-8 -*-
"""The CdL model's output figures, for a MARMITES run (user, 2026-10-26),
drawn after the run into <run>/_output/:

  pond_recharge_analysis.png    each pond's water over the run, from MF6's
                                own LAK budget: the runoff it captures, what
                                it gives the aquifer and takes from it, its
                                evaporation; the seasons of those against the
                                rainfall; the ponds in the catchment's water
                                balance. (The CdL figure was an offline
                                estimate, its LAK having been retired; here
                                the ponds are simulated.)
  ss_obs_vs_computed_heads.png  a run with a steady first period: the heads
                                that period ends at against the observed
                                ones, point by point
  sp_run_timing.png             wall-clock time and MF6 iterations per stress
                                period, with the rainfall and the streamflow
                                at the outlet

Each figure is skipped with a note when its data are not there.
"""

import os

import numpy as np

import marmites_postprocess as PP

DAYS_PER_YEAR = 365.25
# the pond-budget components as drawn: (key, label, colour); the keys and
# colours are the LAK budget figures'
_POND_TERMS = (('runoff', 'runoff captured (RUNOFF)', '#2ca02c'),
               ('to_aquifer', 'given to the aquifer (LAK -> GW)', '#c49c94'),
               ('from_aquifer', 'taken from the aquifer (GW -> LAK)',
                '#8c564b'),
               ('evaporation', 'evaporation (EVAPORATION)', '#d62728'))


def _text(v):
    return v.decode() if isinstance(v, bytes) else str(v)


def _dates(cMF, ds_ws, n):
    """``n`` dates of the transient periods: the run's own, else the
    dataset's."""
    import matplotlib.dates as mdates
    import pandas as pd
    d = getattr(cMF, 'inputDate', None) if cMF is not None else None
    if d is not None and len(d) >= n:
        return pd.DatetimeIndex(pd.to_datetime(
            mdates.num2date(np.asarray(d, dtype=float)[:n])).tz_localize(None))
    return PP._dates_from_dataset(ds_ws, n)


def _rain(ctx, n):
    """Catchment rainfall [mm/d] per day, the mean over the meteo zones; None
    without ctx."""
    P = getattr(ctx, 'P_veg_zoneSP', None) if ctx is not None else None
    if P is None:
        return None
    P = np.asarray(P, dtype=float)
    return P.reshape(P.shape[0], -1).mean(axis=0)[:n]


def _bars(ax, x, y, width, **kw):
    """Bars per period -- or, over a long run, the same as a filled step
    profile: a day is less than a pixel wide there and its bar vanishes."""
    if len(y) > 300:
        ax.fill_between(x, 0, y, step='mid', lw=0, **kw)
    else:
        ax.bar(x, y, width=width, lw=0, **kw)


def _save(fig, out, fn, written):
    import matplotlib.pyplot as plt
    path = os.path.join(out, fn)
    fig.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    written.append(path)


# --------------------------------------------------------------------- #
# the ponds and the recharge
# --------------------------------------------------------------------- #

def pond_volumes(sim_ws, name):
    """Each pond's mean annual volumes [m3/yr] per LAK component over the
    run, with the per-period rates they come from: ``(names, keys, annual
    [pond, component], rate [period, pond, component] m3/d, perlen [d])``;
    None for a run without LAK."""
    bl = PP.lake_budget_by_lake(sim_ws, name)
    if bl is None:
        return None
    perlen = PP._tdis_perlen(sim_ws)
    if perlen is None:
        perlen = np.ones(max(bl['periods']) + 1)
    perlen = np.asarray(perlen, dtype=float)[bl['periods']]
    rate = bl['rate']
    days = float(perlen.sum())
    annual = np.einsum('p,pnc->nc', perlen, rate) / days * DAYS_PER_YEAR
    return bl['names'], bl['keys'], annual, rate, perlen


def _uzf_recharge(sim_ws, name):
    """The recharge to the water table [m3/d], MF6's UZF -> GW over the run;
    None when the listing does not say."""
    try:
        _df, comp = PP.budget_by_compartment(sim_ws, name)
    except Exception:
        return None
    row = comp[comp['compartment'].astype(str).str.startswith('unsaturated')]
    return float(row['net_rate_m3d'].sum()) if len(row) else None


def _active_area(sim_ws, name):
    """The model's area [m2]: every cell active in some layer."""
    from marmites_cdl_figs import _layered, _polys, load_model
    gwf, mg = load_model(sim_ws, name)
    if mg is None:
        return None, None, None
    polys = _polys(mg)
    act = (_layered(mg.idomain, mg) > 0).any(axis=0)
    a = np.array([0.5 * abs(np.dot(p[:, 0], np.roll(p[:, 1], -1))
                            - np.dot(p[:, 1], np.roll(p[:, 0], -1)))
                  for p in polys])
    return float(a[act].sum()), gwf, mg


def fig_pond_recharge(out, sim_ws, name, ds_ws, cMF, ctx, written):
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.patches import Patch
    pv = pond_volumes(sim_ws, name)
    if pv is None:
        return
    names, keys, annual, rate, perlen = pv
    col = {k: keys.index(k) for k in keys}
    # the ponds in their fid order
    def fid(nm):
        s = ''.join(ch for ch in str(nm) if ch.isdigit())
        return int(s) if s else 0
    order = sorted(range(len(names)), key=lambda i: fid(names[i]))
    labels = [str(names[i]).replace('pond', '') for i in order]
    vals = np.array([[annual[i, col[k]] for k, _l, _c in _POND_TERMS]
                     for i in order])
    net_gw = vals[:, 1] + vals[:, 2]          # + = taken from the aquifer
    gains = net_gw > 0

    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.35,
                          wspace=0.35)
    ax = fig.add_subplot(gs[0, :])
    x = np.arange(len(order))
    w = 0.2
    for j, (_k, lab, c) in enumerate(_POND_TERMS):
        ax.bar(x + (j - 1.5) * w, np.abs(vals[:, j]), width=w, color=c,
               label=lab, edgecolor='k', lw=0.3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    for t, g in zip(ax.get_xticklabels(), gains):
        t.set_color('#d95f02' if g else '#1f78b4')
        t.set_fontweight('bold')
    ax.set_xlabel('pond (fid): blue -- gives the aquifer more than it takes '
                  '(net recharge); orange -- takes more than it gives (gains '
                  'groundwater)', fontsize=9)
    ax.set_ylabel('m3/yr')
    ax.grid(alpha=0.3, axis='y')
    ax.set_title('Per pond: mean annual volumes over the run (MF6 LAK budget, '
                 '%d d)' % int(round(perlen.sum())))

    # the seasons: every pond together, per calendar month
    dates = _dates(cMF, ds_ws, rate.shape[0])
    axm = fig.add_subplot(gs[1, :2])
    if isinstance(dates, pd.DatetimeIndex):
        vol = rate.sum(axis=1) * perlen[:, None]          # m3 per period
        df = pd.DataFrame(vol, index=dates, columns=keys)
        per_month = df.resample('MS').sum()
        n_year = per_month.groupby(per_month.index.month).size()
        clim = per_month.groupby(per_month.index.month).sum().div(n_year,
                                                                  axis=0)
        m = np.arange(1, 13)
        clim = clim.reindex(m, fill_value=0.0)
        for j, (k, lab, c) in enumerate(_POND_TERMS):
            axm.bar(m + (j - 1.5) * w, np.abs(clim[k].to_numpy()), width=w,
                    color=c, edgecolor='k', lw=0.3)
        axm.set_xticks(m)
        axm.set_xlabel('month')
        axm.set_ylabel('m3/month, all ponds')
        axm.grid(alpha=0.3, axis='y')
        rain = _rain(ctx, len(dates))
        if rain is not None:
            pr = pd.Series(rain, index=dates[:len(rain)]).resample('MS').sum()
            pclim = pr.groupby(pr.index.month).mean().reindex(m)
            ax2 = axm.twinx()
            ax2.plot(m, pclim.to_numpy(), 'k.-', label='rainfall')
            ax2.set_ylabel('rainfall (mm/month)')
            ax2.legend(fontsize=8, loc='upper center')
        axm.set_title('Seasons (monthly means over the run)')
    else:
        axm.axis('off')

    # the ponds in the catchment's water balance
    axt = fig.add_subplot(gs[1, 2])
    axt.axis('off')
    area, gwf, mg = _active_area(sim_ws, name)
    rch = _uzf_recharge(sim_ws, name)
    tot = vals.sum(axis=0)
    lines = ['WATER BALANCE CONTEXT', '-' * 30]
    if area:
        lines.append('model area:          %8.2f km2' % (area / 1e6))
    rain = _rain(ctx, len(perlen))
    if rain is not None:
        lines.append('mean annual rainfall:%8.0f mm/yr'
                     % (np.mean(rain) * DAYS_PER_YEAR))
    rch_yr = None
    if rch is not None and area:
        rch_yr = rch * DAYS_PER_YEAR
        lines += ['recharge to the water', '  table (UZF -> GW): %8.0f mm/yr'
                  % (rch_yr / area * 1000.0),
                  '                     %8.3f hm3/yr' % (rch_yr / 1e6)]
    pond_area = 0.0
    if gwf is not None:
        from marmites_cdl_figs import lakes
        pond_area = sum(lk['area'] for lk in lakes(gwf, mg, sim_ws))
    lines += ['', 'THE %d PONDS (%.0f m2)' % (len(order), pond_area),
              '-' * 30,
              'runoff captured:     %8.0f m3/yr' % tot[0],
              'given to the aquifer:%8.0f m3/yr' % abs(tot[1]),
              'taken from it:       %8.0f m3/yr' % tot[2],
              'evaporated:          %8.0f m3/yr' % abs(tot[3]),
              'net to the aquifer:  %+8.0f m3/yr' % -(tot[1] + tot[2])]
    if rch_yr:
        lines.append('  = %+.2f %% of the recharge'
                     % (-100.0 * (tot[1] + tot[2]) / rch_yr))
    lines += ['', '%d pond(s) gain groundwater,' % int(gains.sum()),
              '%d give the aquifer more' % int((~gains).sum())]
    axt.text(0.0, 1.0, '\n'.join(lines), va='top', ha='left',
             family='monospace', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='0.96', edgecolor='0.6'))
    fig.legend(handles=[Patch(color=c, label=lab)
                        for _k, lab, c in _POND_TERMS],
               loc='lower center', ncol=4, fontsize=8)
    fig.suptitle('Ponds and recharge -- what each pond captures, gives the '
                 'aquifer and takes from it', fontsize=12, fontweight='bold')
    _save(fig, out, 'pond_recharge_analysis.png', written)
    pd.DataFrame(vals, index=['pond%s' % s for s in labels],
                 columns=[k for k, _l, _c in _POND_TERMS]).to_csv(
        os.path.join(out, 'pond_recharge_analysis.csv'),
        index_label='pond', float_format='%.1f')


# --------------------------------------------------------------------- #
# the steady state against the observed heads
# --------------------------------------------------------------------- #

def steady_heads_at_points(sim_ws, name, ds_ws):
    """At each observation point with an observed head series: the head the
    steady first period ends at and the observed mean, min and max.
    ``[(name, layer, computed, obs mean, obs min, obs max)]``; None when the
    run has no steady first period."""
    import flopy
    if not PP.steady_first(sim_ws, name):
        return None
    from marmites_cdl_figs import load_model
    _gwf, mg = load_model(sim_ws, name)
    hds = flopy.utils.HeadFile(os.path.join(sim_ws, '%s.hds' % name))
    kk = [k for k in hds.get_kstpkper() if int(k[1]) == 0]
    h = np.asarray(hds.get_data(kstpkper=kk[-1]), dtype=float).reshape(
        int(mg.nlay), -1)
    vertex = mg.grid_type == 'vertex'
    out = []
    for p in PP.obs_points(ds_ws):
        s = PP.obs_series(ds_ws, p['name'])
        if s is None or not len(s):
            continue
        v = s['head'].to_numpy(dtype=float)
        v = v[np.isfinite(v) & (np.abs(v) < 9000.0)]
        if not v.size:
            continue
        try:
            c = mg.intersect(p['x'], p['y'])
        except Exception:
            continue
        node = int(c) if vertex else int(c[0]) * int(mg.ncol) + int(c[1])
        L = min(max(int(p['lay']) - 1, 0), int(mg.nlay) - 1)
        hc = h[L, node]
        if not np.isfinite(hc) or abs(hc) > 1e29:
            continue
        out.append((p['name'], L + 1, float(hc), float(v.mean()),
                    float(v.min()), float(v.max())))
    return out


def fig_ss_heads(out, sim_ws, name, ds_ws, res, written, verbose=True):
    import matplotlib.pyplot as plt
    import pandas as pd
    pts = steady_heads_at_points(sim_ws, name, ds_ws)
    if pts is None:
        if verbose:
            print('   SS heads figure: this run has no steady first period '
                  '(it started from saved heads)')
        return
    if not pts:
        return
    df = pd.DataFrame(pts, columns=['point', 'layer', 'computed_m',
                                    'obs_mean_m', 'obs_min_m', 'obs_max_m'])
    err = df['computed_m'] - df['obs_mean_m']
    rmse, bias = float(np.sqrt(np.mean(err ** 2))), float(err.mean())
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    lo = min(df['obs_min_m'].min(), df['computed_m'].min())
    hi = max(df['obs_max_m'].max(), df['computed_m'].max())
    pad = 0.05 * (hi - lo + 1.0)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], 'k--', lw=1,
            label='1:1')
    om, cm = df['obs_mean_m'].to_numpy(), df['computed_m'].to_numpy()
    ax.errorbar(om, cm,
                xerr=np.vstack([om - df['obs_min_m'].to_numpy(),
                                df['obs_max_m'].to_numpy() - om]),
                fmt='o', ms=9, color='#8c564b', ecolor='0.55', capsize=3,
                label='observation points (observed range)')
    for p, lay, x, y in zip(df['point'], df['layer'], om, cm):
        ax.annotate('%s (L%d)' % (p, lay), (float(x), float(y)),
                    xytext=(6, 4), textcoords='offset points', fontsize=9,
                    fontweight='bold')
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_aspect('equal')
    ax.grid(alpha=0.3)
    ax.set_xlabel('OBSERVED head (m a.s.l.), mean of the record')
    ax.set_ylabel('COMPUTED steady-state head (m a.s.l.)')
    ax.legend(loc='upper left')
    rch = res.get('steady_recharge') if res is not None else None
    ax.set_title('Steady state: observed vs computed heads%s\nRMSE = %.2f m'
                 '   bias = %+.2f m   (n = %d)'
                 % (' (recharge %.0f mm/yr)' % (float(rch) * 1000.0
                                               * DAYS_PER_YEAR)
                    if rch is not None else '', rmse, bias, len(df)))
    _save(fig, out, 'ss_obs_vs_computed_heads.png', written)
    df.to_csv(os.path.join(out, 'ss_obs_vs_computed_heads.csv'), index=False,
              float_format='%.3f')


# --------------------------------------------------------------------- #
# the run's time per stress period
# --------------------------------------------------------------------- #

def fig_sp_timing(out, sim_ws, name, ds_ws, cMF, ctx, res, written):
    import matplotlib.pyplot as plt
    import pandas as pd
    if res is None or 'outer_iters' not in res:
        return
    it = np.asarray(res['outer_iters'], dtype=float)
    n = it.size
    sec = res.get('sp_seconds')
    sec = (np.asarray(sec, dtype=float)[:n] if sec is not None
           and np.size(sec) and np.nansum(sec) > 0 else None)
    kind = _text(res.get('iters_kind', 'outer'))
    dates = _dates(cMF, ds_ws, n)
    is_dt = isinstance(dates, pd.DatetimeIndex)
    x = dates if is_dt else np.arange(1, n + 1)
    fig, (ax, axq) = plt.subplots(2, 1, figsize=(14, 8), sharex=True,
                                  gridspec_kw={'height_ratios': [2.2, 1.0]})
    ss = res.get('ss_seconds')
    width = 1.0 if is_dt else 0.8
    if sec is not None:
        _bars(ax, x, sec, width, color='tab:red')
        ax.set_ylabel('wall-clock time per period (s)')
        y = sec
        total = float(sec.sum()) + (float(ss) if ss is not None else 0.0)
        head = 'total %dh%02dm over %d periods' % (total // 3600,
                                                  (total % 3600) // 60, n)
        ax2 = ax.twinx()
        ax2.plot(x, it, color='0.25', lw=0.5)
        ax2.set_ylabel('MF6 %s iterations' % kind, fontsize=8)
    else:
        _bars(ax, x, it, width, color='tab:red')
        ax.set_ylabel('MF6 %s iterations per period' % kind)
        y = it
        head = ('%d periods -- this run did not record its time per period '
                '(runs before 2026-10-26): iterations only' % n)
    if ss is not None and sec is not None:
        head += '; the steady first period took %.0f s (%s iterations)' % (
            float(ss), res.get('ss_iters', '?'))
    for k in np.argsort(y)[::-1][:3]:
        lab = 'SP%d' % (k + 1)
        if is_dt:
            lab += ' %s' % dates[k].strftime('%Y-%m-%d')
        ax.annotate(lab, (x[k], y[k]), xytext=(0, 3),
                    textcoords='offset points', ha='center', fontsize=7)
    ax.set_title('Run time per stress period (%s)' % head, fontsize=10)
    ax.grid(alpha=0.3, axis='y')
    # what the slow periods coincide with: the rain and the flow at the outlet
    q = None
    try:
        sf = PP.sfr_observations(sim_ws, name)
        if sf is not None:
            c = next((c for c in sf.columns if c.startswith('outflow')), None)
            if c is not None:
                q = np.abs(sf[c].to_numpy(dtype=float))[:n]
    except Exception:
        q = None
    if q is not None:
        axq.plot(x[:q.size], np.where(q > 0, q, np.nan), color='tab:blue',
                 lw=0.8, label='streamflow at the outlet')
        axq.set_yscale('log')
        axq.set_ylabel('outlet flow (m3/d)')
    rain = _rain(ctx, n)
    if rain is not None:
        axr = axq.twinx()
        _bars(axr, x[:rain.size], rain, width, color='0.55',
              label='rainfall')
        axr.set_ylim(max(rain.max(), 1.0) * 2.5, 0)
        axr.set_ylabel('rainfall (mm/d)', fontsize=8)
    axq.grid(alpha=0.3)
    axq.set_xlabel('date' if is_dt else 'stress period')
    _save(fig, out, 'sp_run_timing.png', written)
    if sec is not None:
        pd.DataFrame({'period': np.arange(1, n + 1), 'seconds': sec,
                      'iterations': it.astype(int)},
                     index=dates if is_dt else None).to_csv(
            os.path.join(out, 'sp_run_timing.csv'),
            index_label='date' if is_dt else None, float_format='%.3f')


# --------------------------------------------------------------------- #
# all of them
# --------------------------------------------------------------------- #

def cdl_output_figures(out, sim_ws, name, ds_ws, cMF=None, ctx=None,
                       res=None, verbose=True):
    """The CdL output figures into ``out``; the files written. A figure that
    cannot be drawn is skipped with a note."""
    import matplotlib
    matplotlib.use('agg')
    os.makedirs(out, exist_ok=True)
    written = []
    jobs = [('pond recharge', lambda: fig_pond_recharge(
                out, sim_ws, name, ds_ws, cMF, ctx, written)),
            ('SS heads', lambda: fig_ss_heads(out, sim_ws, name, ds_ws, res,
                                              written, verbose=verbose)),
            ('run timing', lambda: fig_sp_timing(out, sim_ws, name, ds_ws,
                                                 cMF, ctx, res, written))]
    for what, job in jobs:
        try:
            job()
        except Exception as exc:                 # pragma: no cover
            if verbose:
                print('   %s figure skipped: %r' % (what, exc))
    if verbose:
        print('CdL output figures: %d file(s)' % len(written))
    return written
