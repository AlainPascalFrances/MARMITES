# -*- coding: utf-8 -*-
"""The Sankey's time windows (2026-09-25).

The one-year La Mata run (31 May 2008 .. 30 May 2009) contains no full
hydrological year, and its "whole period" Sankey showed October..May
scaled to a year -- P 401 mm/y where the run had 353, groundwater ET 7
where it had 38 -- under the caption "Average of the 1 hydrological
year(s)". Every window also dropped its last day and summed in float16.
"""

import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

mpl = pytest.importorskip('matplotlib')
pytest.importorskip('matplotlib.dates')
pd = pytest.importorskip('pandas')


def _pp():
    spec = importlib.util.spec_from_file_location(
        'marmites_postprocess_hy', os.path.join(CODE, 'ppMF6',
                                                'marmites_postprocess.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _dates(start, n):
    return mpl.dates.date2num(pd.date_range(start, periods=n, freq='D'))


def test_a_run_without_a_full_year_is_drawn_whole():
    hy, years = _pp()._hydro_year_index(_dates('2008-05-31', 365), 10)
    assert hy == [0, 0, 364, 364]
    assert years == [], 'no hydrological year of its own'


def test_full_years_still_run_october_to_september():
    D = _dates('2008-05-31', 365 * 3)
    hy, years = _pp()._hydro_year_index(D, 10)
    day = lambda s: int(np.argmax(D == mpl.dates.datestr2num(s)))
    assert years == [2008, 2009]
    assert hy[1] == day('2008-10-01') and hy[2] == day('2009-10-01')
    assert hy[-2] == day('2010-09-30')


def test_every_window_is_inclusive_and_in_double_precision():
    src = open(os.path.join(CODE, 'MARMITESutilities', 'MARMITESplot',
                            'MARMITESplot_v3.py'), encoding='utf-8').read()
    body = src[src.index('def plotWBsankey'):]
    body = body[:body.index('\ndef ', 10)]
    assert 'np.float16(' not in body
    assert '[i:indexend]' not in body
    assert '[i:indexend + 1]' in body
    # the last year ends ON its 30 September
    assert 'else indexTime[-2])' in body
    # a whole-run panel says what it is
    assert 'Whole run, %s to %s (%d days)' in body


def test_summing_the_windows_covers_every_day_once():
    """The rule as plotWBsankey applies it: the per-year windows tile the
    whole-period window exactly."""
    D = _dates('2008-05-31', 365 * 3)
    it, _years = _pp()._hydro_year_index(D, 10)
    days = []
    for k in range(1, len(it) - 2):
        i = it[k]
        end = it[k + 1] - 1 if k + 1 < len(it) - 2 else it[-2]
        days += list(range(i, end + 1))
    assert days == list(range(it[1], it[-2] + 1))
