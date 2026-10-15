# -*- coding: utf-8 -*-
"""The observations of the La Mata PEST++ calibration (WP7): what is
compared, and the model's value of each, from a run's observation exports.

TARGETS are built once, at setup, from the observed records and the dates
of the run (``build_targets``): one row per observation with its name,
group, the model quantity it is compared with, the observed value, a weight
and a measurement standard deviation. COLLECTION (``collect``) reads a
run's exports (obs_heads.csv, obs_sm.csv, obs_sfr.csv, obs_et.csv, one row
per stress period -- the WP6.4 contract) and gives the model's value of
every target, in the targets' order: what the forward run writes for the
instruction file.

Groups (docs/MM-MF6/WP7_PESTpp_LaMata_draft.md, decisions of 2026-10-11):

    hds_p0   P0's head, near-daily, measured           (heads 35 %)
    hds_c    C1..C5, W1 heads, quarterly: SOFT data, a synthetic annual
             cycle and a constant -- each weighs 25 % of a P0 head
             (user, 2026-10-12), whatever share of phi that gives
    sm_p0    P0's soil moisture, 2 horizons             (sm 35 %, half)
    sm_sm    SM's soil moisture, 2 horizons             (sm 35 %, half)
    q_log    outlet flow, log10(Q + 0.01 mm/d), daily   (flow 20 %: 10)
    q_dry    the gauge's DRY days (Q <= 0.005 mm/d): the model's flow
             against zero, in mm/d -- La Mata has no summer flow, and the
             pilot kept ~0.05 mm/d of baseflow (2026-10-13)  (flow 20 %: 5)
    q_year   outlet volume per hydrological year, over
             the days the gauge has                     (flow 20 %: 5)
    et_rs    MODIS MOD16A2GF catchment ET per composite (ET 10 %)
    et_pt    measured ET at the observation points -- the dataset's
             <obs.aet_prefix>_<point>.txt: the EC tower, lysimeters --
             daily, from the exports' obs_et_mmd        (0: validation)

A target before the calibration window (the warm-up year) is kept, with
weight 0, so it is reported and never fitted.
"""

import datetime as dt
import os

import numpy as np

__all__ = ['GROUP_SHARE', 'NOISE', 'build_targets', 'collect',
           'phi_by_group', 'balance_weights']

# share of the initial phi per group (the user's 35/35/20/10, split);
# the soft data (SOFT) take no share of their own
GROUP_SHARE = {'hds_p0': 0.35, 'hds_c': 0.0, 'sm_p0': 0.175,
               'sm_sm': 0.175, 'q_log': 0.10, 'q_dry': 0.05, 'q_year': 0.05,
               'et_rs': 0.10, 'et_pt': 0.0}
# soft data: weighted at a fraction of a measured group's weight per
# observation (user, 2026-10-12: at most 25 % of the measured heads')
SOFT = {'hds_c': ('hds_p0', 0.25)}
# measurement noise (standard deviation) per group: absolute, or 'rel' x
# the observed value with a floor
NOISE = {'hds_p0': 0.1, 'hds_c': 0.4, 'sm_p0': 0.01, 'sm_sm': 0.01,
         'q_log': 0.1, 'q_dry': 0.01, 'q_year': ('rel', 0.10, 1.0),
         'et_rs': ('rel', 0.20, 1.0), 'et_pt': ('rel', 0.20, 0.1)}
Q_EPS = 0.01                 # mm/d, so a dry day has a finite log
Q_DRY = 0.005                # mm/d: the gauge's dry day


def _sd(group, obs):
    n = NOISE[group]
    if isinstance(n, tuple):
        return max(n[1] * abs(float(obs)), n[2])
    return float(n)


def _hydro_year(d, start_month=10):
    return d.year + 1 if d.month >= start_month else d.year


def build_targets(exports_dir, calib_start, modis=None, hy_month=10):
    """The targets, from a run's exports (which carry the observed values
    beside the simulated ones) plus the MODIS table.

    ``calib_start``: the first day fitted (a ``date``); targets before it
    get weight 0. ``modis``: DataFrame of pest_modis_et (start, end, et_mm).
    Returns a DataFrame with obsnme, obgnme, kind, key, obsval, weight (1 or
    0, balanced later), standard_deviation."""
    import pandas as pd
    rows = []
    cs = pd.Timestamp(calib_start)

    def add(name, group, kind, key, obs, date):
        if obs is None or not np.isfinite(obs):
            return
        fitted = GROUP_SHARE[group] > 0 or group in SOFT
        rows.append({'obsnme': name, 'obgnme': group, 'kind': kind,
                     'key': key, 'obsval': float(obs),
                     'weight': 1.0 if (pd.Timestamp(date) >= cs
                                       and fitted) else 0.0,
                     'standard_deviation': _sd(group, obs)})

    h = pd.read_csv(os.path.join(exports_dir, 'obs_heads.csv'))
    for r in h.itertuples():
        g = 'hds_p0' if r.point.upper() == 'P0' else 'hds_c'
        add('h_%s_%04d' % (r.point.lower(), r.period), g, 'head',
            '%s|%d' % (r.point, r.period), r.obs_head_m, r.date)
    s = pd.read_csv(os.path.join(exports_dir, 'obs_sm.csv'))
    for r in s.itertuples():
        p = r.point.upper()
        if p not in ('P0', 'SM'):
            continue
        add('sm_%s_%d_%04d' % (p.lower(), r.soil_layer, r.period),
            'sm_%s' % p.lower(), 'sm', '%s|%d|%d' % (p, r.soil_layer,
                                                     r.period),
            r.obs_theta, r.date)
    q = pd.read_csv(os.path.join(exports_dir, 'obs_sfr.csv'))
    for r in q.itertuples():
        if np.isfinite(r.obs_q_mmd):
            add('q_%04d' % r.period, 'q_log', 'qlog', str(r.period),
                np.log10(max(r.obs_q_mmd, 0.0) + Q_EPS), r.date)
        if np.isfinite(r.obs_q_mmd) and r.obs_q_mmd <= Q_DRY:
            add('qd_%04d' % r.period, 'q_dry', 'qlin', str(r.period), 0.0,
                r.date)
    qd = q.dropna(subset=['obs_q_mmd']).copy()
    qd['hy'] = [_hydro_year(pd.Timestamp(d), hy_month) for d in qd['date']]
    for hy, grp in qd.groupby('hy'):
        # fitted only when the whole hydrological year is in the window
        add('qv_%d' % hy, 'q_year', 'qyear',
            ';'.join(str(p) for p in grp['period']),
            float(grp['obs_q_mmd'].sum()),
            pd.Timestamp(int(hy) - 1, hy_month, 1))
    e = pd.read_csv(os.path.join(exports_dir, 'obs_et.csv'))
    cat = e[e['site'] == 'catchment']
    if modis is not None and len(cat):
        days = pd.to_datetime(cat['date']).dt.date.tolist()
        per = dict(zip(days, cat['period']))
        for r in modis.itertuples():
            a = dt.date.fromisoformat(str(r.start))
            b = dt.date.fromisoformat(str(r.end))
            ps = [per.get(a + dt.timedelta(days=k))
                  for k in range((b - a).days + 1)]
            if any(p is None for p in ps) or not np.isfinite(r.et_mm):
                continue                 # a composite the run does not cover
            add('etm_%d%03d' % (r.year, r.doy), 'et_rs', 'etsum',
                ';'.join(str(p) for p in ps), r.et_mm, a)
    # the measured ET at the points, read by the run itself (2026-10-15:
    # was an Excel sheet read here); the catchment's is MODIS, above
    for r in e[e['site'] != 'catchment'].itertuples():
        add(r.obsnme, 'et_pt', 'etpt', '%s|%d' % (r.site, r.period),
            r.obs_et_mmd, r.date)
    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError('no observation found in %s' % exports_dir)
    if df['obsnme'].duplicated().any():
        raise ValueError('duplicate observation names: %s' % df.loc[
            df['obsnme'].duplicated(), 'obsnme'].tolist()[:5])
    return df


def collect(exports_dir, targets):
    """The model's value of every target, in the targets' order: a
    DataFrame obsnme, value. A value the run did not produce is a large
    finite sentinel, so PEST++ sees a bad fit, never a parse error."""
    import pandas as pd
    h = pd.read_csv(os.path.join(exports_dir, 'obs_heads.csv'))
    # MMsoil's CORRECTED head: MF6's plus the saturated soil above the soil
    # base (user, 2026-10-13) -- MF6's top is the soil base, so its head
    # cannot follow a water table that rises into the soil, as P0's does in
    # winter (observed up to the land surface). Equal to MF6's head where
    # the soil is not saturated; MF6's where a run has no corrected head.
    sim_h = h['sim_hcorr_m'].where(np.isfinite(h['sim_hcorr_m']),
                                   h['sim_head_m'])
    hv = dict(zip(h['point'].astype(str) + '|' + h['period'].astype(str),
                  sim_h))
    s = pd.read_csv(os.path.join(exports_dir, 'obs_sm.csv'))
    sv = dict(zip(s['point'].astype(str).str.upper() + '|'
                  + s['soil_layer'].astype(str) + '|'
                  + s['period'].astype(str), s['sim_theta']))
    q = pd.read_csv(os.path.join(exports_dir, 'obs_sfr.csv'))
    qv = dict(zip(q['period'], q['sim_q_mmd']))
    e = pd.read_csv(os.path.join(exports_dir, 'obs_et.csv'))
    cat = e[e['site'] == 'catchment']
    etc = dict(zip(cat['period'], cat['sim_et_total_mmd']))
    etp = dict(zip(e['site'].astype(str) + '|' + e['period'].astype(str),
                   e['sim_et_total_mmd']))
    out = []
    for r in targets.itertuples():
        k, key = r.kind, str(r.key)
        if k == 'head':
            v = hv.get(key)
        elif k == 'sm':
            v = sv.get(key)
        elif k == 'qlog':
            x = qv.get(int(key))
            v = None if x is None else np.log10(max(float(x), 0.0) + Q_EPS)
        elif k == 'qlin':
            x = qv.get(int(key))
            v = None if x is None else max(float(x), 0.0)
        elif k == 'qyear':
            xs = [qv.get(int(p)) for p in key.split(';')]
            v = None if any(x is None for x in xs) else float(np.sum(xs))
        elif k == 'etsum':
            xs = [etc.get(int(p)) for p in key.split(';')]
            v = None if any(x is None for x in xs) else float(np.sum(xs))
        elif k == 'etpt':
            v = etp.get(key)
        else:
            raise ValueError('unknown target kind %r' % k)
        if v is None or not np.isfinite(v):
            v = 1.0e10
        out.append((r.obsnme, float(v)))
    return pd.DataFrame(out, columns=['obsnme', 'value'])


def phi_by_group(targets, sim, weights=None):
    """Sum of squared weighted residuals per group."""
    w = targets['weight'] if weights is None else weights
    r = (targets['obsval'].to_numpy() - sim['value'].to_numpy()) * w
    return (pd_series(r ** 2, targets['obgnme'])).groupby(level=0).sum()


def pd_series(values, index):
    import pandas as pd
    return pd.Series(np.asarray(values, float), index=np.asarray(index))


def balance_weights(targets, sim, total=1.0e4):
    """Weights scaled per group so each carries GROUP_SHARE of the initial
    phi (``total``); weight-0 targets stay 0. Returns the new weights."""
    w = targets['weight'].to_numpy(dtype=float).copy()
    r = targets['obsval'].to_numpy() - sim['value'].to_numpy()
    g = targets['obgnme'].to_numpy()
    for grp, share in GROUP_SHARE.items():
        if grp in SOFT:
            continue
        m = (g == grp) & (w > 0)
        if not m.any() or share <= 0:
            w[g == grp] = 0.0
            continue
        phi = float(np.sum((r[m] * w[m]) ** 2))
        if phi <= 0:
            continue
        w[m] *= np.sqrt(share * total / phi)
    # soft data: a fraction of the measured group's weight per observation
    for grp, (ref, frac) in SOFT.items():
        m = (g == grp) & (w > 0)
        mr = (g == ref) & (w > 0)
        if m.any():
            w[m] = frac * (float(np.median(w[mr])) if mr.any() else 0.0)
    return w
