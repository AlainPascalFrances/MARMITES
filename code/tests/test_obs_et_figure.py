# -*- coding: utf-8 -*-
"""The ET figure (user, 2026-10-28: "where are the ET graphs for
calibration, i.e. computed ET vs MODIS one and EC one?").

From the observation export obs_et.csv: every point with a measured ET
series day by day, the catchment against MODIS summed per 8-day composite
-- what the calibration fits -- over whole hydrological years only, a
composite the window cuts left out.
"""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

pd = pytest.importorskip('pandas')
import matplotlib                                       # noqa: E402
matplotlib.use('agg')
import matplotlib.figure                                # noqa: E402

import marmites_postprocess as PP                       # noqa: E402

DATES = pd.date_range('2008-05-31', '2010-10-15', freq='D')
LO, HI = pd.Timestamp('2008-10-01'), pd.Timestamp('2010-09-30')


def _write(out, modis=True):
    n = len(DATES)
    sim = 1.0 + 0.5 * np.sin(np.arange(n) / 58.0)
    ec = np.full(n, np.nan)
    ec[200:260] = sim[200:260]                  # 60 days, a perfect fit
    rows = []
    for site, obs in (('EC', ec),
                      ('catchment', 1.25 * sim if modis
                       else np.full(n, np.nan))):
        rows.append(pd.DataFrame({
            'obsnme': ['et_%s_%04d' % (site.lower(), k + 1)
                       for k in range(n)],
            'site': site, 'period': np.arange(1, n + 1),
            'date': DATES.strftime('%Y-%m-%d'), 'sim_et_total_mmd': sim,
            'obs_et_mmd': obs}))
    # a point with no measured series gets no row
    rows.append(rows[0].assign(site='P0', obs_et_mmd=np.nan))
    pd.concat(rows).to_csv(os.path.join(out, 'obs_et.csv'), index=False)


@pytest.fixture
def spy(monkeypatch):
    seen = {}
    real = matplotlib.figure.Figure.savefig

    def save(fig, fn, *a, **k):
        seen[os.path.basename(str(fn))] = fig
        return real(fig, fn, *a, **k)
    monkeypatch.setattr(matplotlib.figure.Figure, 'savefig', save)
    return seen


def _whole_composites(lo, hi):
    n = 0
    for y in range(lo.year, hi.year + 1):
        for k in range(46):
            a = pd.Timestamp(y, 1, 1) + pd.Timedelta(days=8 * k)
            b = min(a + pd.Timedelta(days=7), pd.Timestamp(y, 12, 31))
            if a.year == y and a >= lo and b <= hi:
                n += 1
    return n


def test_the_modis_composites_start_on_1_january():
    d = pd.to_datetime(['2008-01-01', '2008-01-08', '2008-01-09',
                        '2008-12-31', '2009-01-01'])
    assert PP.modis_composites(d) == [(2008, 0), (2008, 0), (2008, 1),
                                      (2008, 45), (2009, 0)]


def test_ec_by_day_and_the_catchment_by_composite(tmp_path, spy):
    _write(str(tmp_path))
    files = PP._fig_obs_et(str(tmp_path), hydro_year_start=10, verbose=False)
    assert [os.path.basename(f) for f in files] == ['obs_et.png']
    ax = spy['obs_et.png'].axes
    titles = [a.get_title() for a in ax]
    # two rows (EC, the catchment) of a series and a scatter; P0 has none
    assert len(ax) == 4
    assert titles[0].startswith('EC:') and '60 days' in titles[0]
    assert 'NSE 1.00' in titles[1] and 'bias +0 %' in titles[1]
    assert titles[2].startswith('catchment:') and 'MODIS' in titles[2]
    assert '%d composites' % _whole_composites(LO, HI) in titles[2]
    assert 'bias -20 %' in titles[3]          # simulated = MODIS / 1.25
    # whole hydrological years only
    x = pd.to_datetime(ax[0].get_lines()[0].get_xdata())
    assert x.min() == LO and x.max() == HI
    assert 'whole hydrological years' in spy['obs_et.png']._suptitle.get_text()


def test_without_modis_the_catchment_row_says_so(tmp_path, spy):
    _write(str(tmp_path), modis=False)
    PP._fig_obs_et(str(tmp_path), verbose=False)
    ax = spy['obs_et.png'].axes
    t = ax[2].get_title()
    assert 'no MODIS' in t and '_catchment.txt' in t
    assert not ax[3].get_visible() or not ax[3].axison


def test_no_export_no_figure(tmp_path):
    assert PP._fig_obs_et(str(tmp_path), verbose=False) == []


def test_the_run_draws_it():
    src = open(os.path.join(CODE, 'ppMF6', 'marmites_postprocess.py'),
               encoding='utf-8').read()
    body = src[src.index('def run_postproc('):]
    body = body[:body.index('\ndef ')]
    assert '_fig_obs_et(out, hydro_year_start=hydro_year_start' in body
