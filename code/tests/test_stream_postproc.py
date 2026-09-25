# -*- coding: utf-8 -*-
"""WP3.6 -- the stream in the post-processing.

The SFR continuous observations collapse to one row per stress period (the
time-weighted mean of the ATS steps, the steady period dropped), the gauge
is read as the legacy driver read it, and the outlet hydrograph and the
network budget are drawn from them. Synthetic output where the answer is
known; no MODFLOW run.
"""
import importlib.util
import os
import sys
import types

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in ('', 'MARMITESutilities', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

pd = pytest.importorskip('pandas')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


PP = _load('marmites_postprocess_s',
           os.path.join(TRUNK, 'ppMF6', 'marmites_postprocess.py'))


def _sim(ws, perlen, steady=True, rows=(), cols=('outflow',)):
    """A model workspace with just what the reader looks at."""
    os.makedirs(ws, exist_ok=True)
    with open(os.path.join(ws, 'mfsim.nam'), 'w') as fh:
        fh.write('BEGIN timing\n  TDIS6  m.tdis\nEND timing\n')
    with open(os.path.join(ws, 'm.tdis'), 'w') as fh:
        fh.write('BEGIN options\n  TIME_UNITS days\nEND options\n'
                 'BEGIN dimensions\n  NPER %d\nEND dimensions\n'
                 'BEGIN perioddata\n' % len(perlen))
        for p in perlen:
            fh.write('  %r 1 1.0\n' % float(p))
        fh.write('END perioddata\n')
    with open(os.path.join(ws, 'm.sto'), 'w') as fh:
        fh.write('BEGIN period 1\n  %s\nEND period 1\n'
                 % ('STEADY-STATE' if steady else 'TRANSIENT'))
        if steady:
            fh.write('BEGIN period 2\n  TRANSIENT\nEND period 2\n')
    with open(os.path.join(ws, 'm.obs.sfr.csv'), 'w') as fh:
        fh.write('time,%s\n' % ','.join(c.upper() for c in cols))
        for r in rows:
            fh.write(','.join('%r' % float(v) for v in r) + '\n')
    return ws


# --------------------------- one row per period ------------------------- #

def test_ats_steps_collapse_to_their_time_weighted_period_mean(tmp_path):
    ws = _sim(str(tmp_path), [1.0, 1.0, 1.0], steady=True,
              rows=[(1.0, 100.0),                   # steady period
                    (1.25, 10.0), (2.0, 20.0),      # day 1, two ATS steps
                    (3.0, 30.0)])                   # day 2, one step
    per = PP.sfr_observations(ws, 'm')
    assert len(per) == 2                            # the steady day is gone
    assert per['outflow'].tolist() == pytest.approx([17.5, 30.0])
    assert per['perlen'].tolist() == [1.0, 1.0]


def test_the_row_count_does_not_depend_on_the_stepping(tmp_path):
    a = _sim(str(tmp_path / 'a'), [1.0] * 4, steady=False,
             rows=[(1, 5), (2, 5), (3, 5), (4, 5)])
    b = _sim(str(tmp_path / 'b'), [1.0] * 4, steady=False,
             rows=[(0.1, 5), (0.5, 5), (1, 5), (2, 5), (2.3, 5), (2.31, 5),
                   (3, 5), (4, 5)])
    pa, pb = PP.sfr_observations(a, 'm'), PP.sfr_observations(b, 'm')
    assert len(pa) == len(pb) == 4
    assert pb['outflow'].tolist() == pytest.approx([5.0] * 4)


def test_without_a_steady_period_the_first_day_is_kept(tmp_path):
    ws = _sim(str(tmp_path), [1.0, 1.0], steady=False,
              rows=[(1.0, 7.0), (2.0, 9.0)])
    assert PP.sfr_observations(ws, 'm')['outflow'].tolist() == [7.0, 9.0]


def test_no_sfr_no_observations(tmp_path):
    assert PP.sfr_observations(str(tmp_path), 'm') is None


# ------------------------------ the gauge ------------------------------- #

def test_the_gauge_is_read_and_its_gaps_skipped(tmp_path):
    fn = tmp_path / 'inputObsRo_catchment.txt'
    fn.write_text('2008-06-11\t0.5\n2008-06-12\t-9999\n2008-06-13\t0.25\n')
    s = PP.obs_streamflow(str(tmp_path))
    assert list(s.index.strftime('%Y-%m-%d')) == ['2008-06-11', '2008-06-13']
    assert s.tolist() == [0.5, 0.25]


def test_a_period_takes_the_mean_of_the_observed_days_inside_it():
    s = pd.Series(np.arange(10.0),
                  index=pd.date_range('2008-01-01', periods=10, freq='D'))
    m = PP._period_mean(s, pd.DatetimeIndex(['2008-01-01', '2008-01-05',
                                             '2008-03-01']), [4.0, 6.0, 1.0])
    assert m[0] == pytest.approx(1.5)          # days 0..3
    assert m[1] == pytest.approx(6.5)          # days 4..9
    assert np.isnan(m[2])                      # nothing observed


# ------------------------------ the figures ----------------------------- #

def _grid(nrow=2, ncol=5, cs=100.0):
    return types.SimpleNamespace(nlay=1, nrow=nrow, ncol=ncol,
                                 delr=np.full(ncol, cs), delc=np.full(nrow, cs),
                                 idomain=np.ones((1, nrow, ncol), dtype=int))


def test_outlet_hydrograph_against_a_perfect_gauge(tmp_path):
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('agg')
    area = 2 * 5 * 100.0 * 100.0               # 0.1 km2
    days = pd.date_range('2008-06-10', periods=6, freq='D')
    q = np.array([100.0, 200.0, 400.0, 300.0, 150.0, 120.0])   # m3/d out
    leak = np.full(6, 5.0)                     # the stream loses 5 m3/d
    evap = np.full(6, 1.0)
    inflow = q + leak + evap                   # closes exactly
    rows = [(0.0 + 1.0, 0, 0, 0, 0)]           # the steady day
    rows += [(2.0 + k, -q[k], inflow[k], evap[k], leak[k])
             for k in range(6)]
    ws = _sim(str(tmp_path / 'ws'), [1.0] * 7, steady=True, rows=rows,
              cols=('outflow', 'net_inflow', 'net_evaporation',
                    'net_leakage'))
    ds = tmp_path / 'ds'
    ds.mkdir()
    (ds / 'inputDATE.txt').write_text(
        '\n'.join(d.strftime('%Y-%m-%d') for d in days) + '\n')
    # the gauge in mm/d over the catchment, from day 2 on
    (ds / 'inputObsRo_catchment.txt').write_text(''.join(
        '%s\t%.12g\n' % (d.strftime('%Y-%m-%d'), qq * 1000.0 / area)
        for d, qq in list(zip(days, q))[1:]))
    out = tmp_path / 'out'
    out.mkdir()
    written = PP._fig_stream(ws, 'm', str(ds), str(out), _grid(),
                             verbose=False)
    names = {os.path.basename(f) for f in written}
    assert {'outlet_streamflow.png', 'outlet_streamflow.csv',
            'budget_sfr_ts.png', 'budget_sfr_period.csv'} <= names
    tab = pd.read_csv(out / 'outlet_streamflow.csv', index_col=0,
                      parse_dates=True)
    assert tab['q_sim_m3d'].tolist() == pytest.approx(q.tolist())
    assert np.isnan(tab['q_obs_m3d'].iloc[0])
    assert tab['q_obs_m3d'].iloc[1:].tolist() == pytest.approx(q[1:].tolist())
    fit = PP._fit(tab['q_sim_m3d'], tab['q_obs_m3d'])
    assert fit['n'] == 5 and fit['nse'] == pytest.approx(1.0)
    assert fit['bias'] == pytest.approx(0.0)
    # the network budget closes: in - evaporation - to the aquifer - out = 0
    terms = [t for t in PP._SFR_TERMS
             if t[0] in ('net_inflow', 'net_evaporation', 'net_leakage')]
    rates = pd.read_csv(out / 'budget_sfr_period.csv', index_col=0)
    closure = sum(sgn * rates[c] for c, sgn, _l, _k in terms) - q
    assert np.allclose(closure, 0.0)


def test_the_catchment_area_of_a_mesh_is_its_active_cells(tmp_path):
    # two unit squares far out in UTM, the second inactive
    verts = np.array([[739300.0, 4553050.0], [739301.0, 4553050.0],
                      [739301.0, 4553051.0], [739300.0, 4553051.0],
                      [739302.0, 4553050.0], [739302.0, 4553051.0]])
    mg = types.SimpleNamespace(nlay=1, verts=verts,
                               iverts=[[0, 3, 2, 1], [1, 2, 5, 4]],
                               idomain=np.array([[1, 0]]))
    assert PP._active_area(mg) == pytest.approx(1.0, abs=1e-9)
