# -*- coding: utf-8 -*-
"""Figures and a summary of a La Mata pestpp-ies run (WP7): the prior Monte
Carlo or an IES iteration, from the master folder.

    python code/tools/pest_lamata_report.py [--iter 0]

Writes into ``<pest_ws>/master/report_<iter>/``:
  phi_groups.png    phi of each group per realisation (box), the base
  ens_<group>.png   the ensemble's envelope (5-95 %, median) against the
                    observations, through time, per group and site
  summary.txt       realisations run / failed, phi per group, the
                    parameters on their bounds
"""

import argparse
import datetime as dt
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

START = dt.date(2008, 5, 31)          # period 1 of the La Mata run


def _period_of(kind, key):
    k = str(key)
    if kind == 'head':
        return int(k.split('|')[1])
    if kind == 'sm':
        return int(k.split('|')[2])
    if kind in ('qlog', 'etec'):
        return int(k)
    return int(k.split(';')[0])        # qyear, etsum: the first day


def _site_of(kind, key, name):
    k = str(key)
    if kind == 'head':
        return k.split('|')[0]
    if kind == 'sm':
        return '%s horizon %s' % tuple(k.split('|')[:2])
    return {'qlog': 'outlet log10(Q+0.01)', 'qyear': 'outlet volume per year',
            'etsum': 'MODIS catchment ET per composite',
            'etec': 'EC tower ET'}.get(kind, name)


def report(master, targets_csv, it=0, out=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    import pyemu
    pst = pyemu.Pst(os.path.join(master, 'lamata.pst'))
    oe_fn = os.path.join(master, 'lamata.%d.obs.jcb' % it)
    pe_fn = os.path.join(master, 'lamata.%d.par.jcb' % it)
    # the ensembles' DataFrames (pyemu wraps them)
    oe = pyemu.ObservationEnsemble.from_binary(pst, oe_fn)._df
    pe = pyemu.ParameterEnsemble.from_binary(pst, pe_fn)._df
    oe.index = oe.index.astype(str)
    pe.index = pe.index.astype(str)
    tg = pd.read_csv(targets_csv).set_index('obsnme')
    out = out or os.path.join(master, 'report_%d' % it)
    os.makedirs(out, exist_ok=True)
    lines = []
    failed = oe.index[(oe.values >= 1.0e9).any(axis=1)].tolist()
    lines.append('iteration %d: %d realisation(s) returned, %d failed '
                 '(sentinel values) %s' % (it, len(oe), len(failed),
                                           failed[:10]))
    # phi per group
    w = pst.observation_data.loc[oe.columns, 'weight'].to_numpy()
    obs = pst.observation_data.loc[oe.columns, 'obsval'].to_numpy()
    grp = pst.observation_data.loc[oe.columns, 'obgnme'].to_numpy()
    ok = oe.drop(index=failed, errors='ignore')
    res2 = ((ok.values - obs) * w) ** 2
    phi = pd.DataFrame({g: res2[:, grp == g].sum(axis=1)
                        for g in sorted(set(grp))}, index=ok.index)
    phi['total'] = phi.sum(axis=1)
    lines.append('phi per group (realisations that ran):')
    lines.append(phi.describe().T[['mean', 'min', '50%', 'max']]
                 .round(1).to_string())
    fig, ax = plt.subplots(figsize=(10, 4.5))
    cols = [c for c in phi.columns if c != 'total' and phi[c].max() > 0]
    ax.boxplot([phi[c] for c in cols], tick_labels=cols)
    if 'base' in phi.index:
        ax.plot(range(1, len(cols) + 1), [phi.loc['base', c] for c in cols],
                'r*', ms=10, label='base')
        ax.legend()
    ax.set_yscale('log')
    ax.set_ylabel('phi (weighted)')
    ax.set_title('La Mata, iteration %d: phi per group over %d realisations'
                 % (it, len(phi)))
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'phi_groups.png'), dpi=120)
    plt.close(fig)
    # parameters on their bounds
    pdat = pst.parameter_data.loc[pe.columns]
    lo = (pe.values <= pdat['parlbnd'].to_numpy() * (1 + 1e-6)).mean(axis=0)
    hi = (pe.values >= pdat['parubnd'].to_numpy() * (1 - 1e-6)).mean(axis=0)
    gb = pd.DataFrame({'pargp': pdat['pargp'].to_numpy(), 'lo': lo,
                       'hi': hi}).groupby('pargp').mean()
    lines.append('share of the ensemble on a bound, per group:')
    lines.append((gb * 100).round(0).astype(int).astype(str).add(' %')
                 .to_string())
    # envelopes against the observations
    tg = tg.loc[[c for c in oe.columns if c in tg.index]]
    tg['period'] = [_period_of(k, key) for k, key in zip(tg['kind'],
                                                         tg['key'])]
    tg['date'] = [START + dt.timedelta(days=int(p) - 1) for p in tg['period']]
    tg['site'] = [_site_of(k, key, n) for k, key, n in
                  zip(tg['kind'], tg['key'], tg.index)]
    q = ok[tg.index].quantile([0.05, 0.5, 0.95])
    for g in sorted(tg['obgnme'].unique()):
        sub = tg[tg['obgnme'] == g]
        sites = list(dict.fromkeys(sub['site']))
        fig, axs = plt.subplots(len(sites), 1, figsize=(11, 2.6 * len(sites)),
                                squeeze=False, sharex=True)
        for ax, s in zip(axs[:, 0], sites):
            t = sub[sub['site'] == s].sort_values('date')
            n = t.index
            ax.fill_between(t['date'], q.loc[0.05, n], q.loc[0.95, n],
                            color='tab:blue', alpha=0.3, label='5-95 %')
            ax.plot(t['date'], q.loc[0.5, n], color='tab:blue', lw=1,
                    label='median')
            if 'base' in ok.index:
                ax.plot(t['date'], ok.loc['base', n], color='k', lw=0.8,
                        label='base')
            fit = t['weight'] > 0
            ax.plot(t['date'][fit], t['obsval'][fit], 'r.', ms=3,
                    label='observed (fitted)')
            ax.plot(t['date'][~fit], t['obsval'][~fit], '.', color='0.6',
                    ms=3, label='observed (weight 0)')
            ax.set_title(s, fontsize=9)
            ax.grid(alpha=0.3)
        h, lb = axs[0, 0].get_legend_handles_labels()
        fig.legend(h, lb, fontsize=7, ncol=5, loc='lower center')
        fig.suptitle('La Mata, iteration %d: %s' % (it, g))
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        fig.savefig(os.path.join(out, 'ens_%s.png' % g), dpi=110)
        plt.close(fig)
    with open(os.path.join(out, 'summary.txt'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print('report written to %s' % out)
    return out


def main(argv=None):
    import mm_paths
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--iter', type=int, default=0)
    ap.add_argument('--master', default=None)
    a = ap.parse_args(argv)
    pest_ws = str(mm_paths.PEST_WS or '')
    master = a.master or os.path.join(pest_ws, 'master')
    report(master, os.path.join(master, 'pest_targets.csv'), it=a.iter)


if __name__ == '__main__':
    main()
