# -*- coding: utf-8 -*-
"""Each pond's water budget per hydrological year (user, 2026-10-09).

One figure with every pond -- a group of bars per year, a bar per pond --
and one with all of them summed, each stacked by flux: from and to the
stream, the outlet flow the mover did not pass on, MARMITES runoff,
rainfall, evaporation, from and to the aquifer, storage and the rest. The
numbers come from the LAK budget file, pond by pond, as sub-step means.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in ('', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

flopy = pytest.importorskip('flopy')
pd = pytest.importorskip('pandas')
import matplotlib  # noqa: E402
matplotlib.use('agg')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PP = _load('_pp_lake_years', os.path.join(TRUNK, 'ppMF6',
                                          'marmites_postprocess.py'))
TERMS = ('GWF', 'RAINFALL', 'EVAPORATION', 'RUNOFF', 'EXT-INFLOW',
         'WITHDRAWAL', 'EXT-OUTFLOW', 'STORAGE', 'CONSTANT', 'FROM-MVR',
         'TO-MVR')


def _rec(nodes, qs):
    a = np.zeros(len(qs), dtype=[('node', 'i4'), ('node2', 'i4'),
                                 ('q', 'f8')])
    a['node'], a['node2'], a['q'] = nodes, nodes, qs
    return a


class _FakeLakCbc:
    """Two ponds, two daily periods; period 1 solved in two half steps.
    Pond 1 has two aquifer connections, one gaining and one losing."""
    STEPS = [(0, 0), (1, 0), (0, 1)]
    TIMES = [0.5, 1.0, 2.0]

    def __init__(self, *a, **k):
        pass

    def get_kstpkper(self):
        return list(self.STEPS)

    def get_times(self):
        return list(self.TIMES)

    def get_unique_record_names(self):
        return [t.rjust(16).encode() for t in TERMS]

    def get_data(self, text=None, kstpkper=None):
        out = []
        for i in range(len(self.STEPS)):
            if text == 'GWF':                       # pond 1: +2 and -6 m3/d
                out.append(_rec([1, 1, 2], [2.0, -6.0, -1.0 * (i + 1)]))
            elif text == 'EVAPORATION':
                out.append(_rec([1, 2], [-1.0, -2.0]))
            elif text == 'FROM-MVR':
                out.append(_rec([1, 2], [100.0 * (i + 1), 0.0]))
            elif text == 'TO-MVR':
                out.append(_rec([1, 2], [-95.0 * (i + 1), 0.0]))
            else:
                out.append(_rec([1, 2], [0.0, 0.0]))
        return out


@pytest.fixture
def lak_ws(tmp_path, monkeypatch):
    (tmp_path / 'toy.lak.cbc').write_bytes(b'')
    (tmp_path / 'toy.sto').write_text('BEGIN PERIOD 1\n  TRANSIENT\n'
                                      'END PERIOD 1\n')
    monkeypatch.setattr(flopy.utils, 'CellBudgetFile', _FakeLakCbc)
    return str(tmp_path)


def test_each_pond_is_kept_apart_and_the_aquifer_split_by_sign(lak_ws):
    by = PP.lake_budget_by_lake(lak_ws, 'toy')
    k = {c: i for i, c in enumerate(by['keys'])}
    r = by['rate']
    assert r.shape == (2, 2, len(PP.LAK_COMPONENTS))
    assert by['names'] == ['lake 1', 'lake 2']
    # pond 1 gains 2 and loses 6 through its two connections: kept apart
    assert r[0, 0, k['from_aquifer']] == pytest.approx(2.0)
    assert r[0, 0, k['to_aquifer']] == pytest.approx(-6.0)
    assert r[0, 1, k['from_aquifer']] == 0.0
    # period 1 = two half steps at 100 and 200 m3/d: their TIME mean
    assert r[0, 0, k['from_stream']] == pytest.approx(150.0)
    assert r[1, 0, k['from_stream']] == pytest.approx(300.0)
    assert r[0, 1, k['to_aquifer']] == pytest.approx(-1.5)
    assert r[0, 1, k['evaporation']] == pytest.approx(-2.0)


def test_hydrological_years_start_on_the_panel_month():
    dates = pd.date_range('2008-05-31', '2009-11-02', freq='D')
    yrs = PP.hydro_years(dates, 10)
    assert [y for y, _ix in yrs] == ['2007/08', '2008/09', '2009/10']
    assert [len(ix) for _y, ix in yrs] == [123, 365, 33]
    assert dates[yrs[1][1][0]] == pd.Timestamp('2008-10-01')
    assert [y for y, _ix in PP.hydro_years(dates, 1)] == ['2008', '2009']


def test_the_volumes_are_rates_times_the_period_lengths():
    rate = np.zeros((4, 1, len(PP.LAK_COMPONENTS)))
    rate[:, 0, 0] = [1.0, 2.0, 3.0, 4.0]
    by = {'rate': rate, 'names': ['p'], 'keys': [c[0] for c in
                                                 PP.LAK_COMPONENTS]}
    dates = pd.to_datetime(['2008-09-29', '2008-09-30', '2008-10-01',
                            '2008-10-02'])
    years, days, vol = PP.lake_budget_years(by, dates, [1, 1, 2, 2], 10)
    assert years == ['2007/08', '2008/09']
    assert list(days) == [2.0, 4.0]
    assert vol[:, 0, 0] == pytest.approx([3.0, 14.0])


def test_the_figures_and_the_table(lak_ws, tmp_path):
    by = PP.lake_budget_by_lake(lak_ws, 'toy')
    dates = pd.to_datetime(['2008-09-30', '2008-10-01'])
    out = tmp_path / 'out'
    out.mkdir()
    files = PP._fig_lake_budget_years(by, dates, [1.0, 1.0], 10, str(out),
                                      verbose=False)
    names = sorted(os.path.basename(f) for f in files)
    assert names == ['lake_budget_years.csv', 'lake_budget_years_by_pond.png',
                     'lake_budget_years_total.png']
    tab = pd.read_csv(os.path.join(str(out), 'lake_budget_years.csv'))
    assert set(tab['pond']) == {'lake 1', 'lake 2', 'all ponds'}
    assert set(tab['component']) == {c[0] for c in PP.LAK_COMPONENTS}
    one = tab[(tab.hydro_year == '2008/09') & (tab.component == 'from_stream')]
    assert dict(zip(one.pond, one.m3)) == {'lake 1': 300.0, 'lake 2': 0.0,
                                           'all ponds': 300.0}


def test_the_run_draws_them_with_the_panel_year():
    src = open(os.path.join(TRUNK, 'ppMF6', 'marmites_postprocess.py'),
               encoding='utf-8').read()
    body = src[src.index('def _fig_lakes('):]
    body = body[:body.index('\ndef ')]
    assert '_fig_lake_budget_years(' in body and 'hydro_year_start' in body
    drv = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert "hydro_year_start=int(getattr(cMF, 'iniMonthHydroYear'" in drv
