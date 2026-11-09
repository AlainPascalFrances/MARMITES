# -*- coding: utf-8 -*-
"""Dry and wet seasons started by the rain (user, 2026-11-09;
MARMITESsurf/rain_seasons.py): La Mata's grass greens up with the first
autumn rain of 10 mm after 1 September and is wilted after 30 days with
less than 10 mm after 1 May, instead of on fixed julian days."""
import os
import sys
from datetime import date, datetime, timedelta

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'MARMITESsurf')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rain_seasons as RS                            # noqa: E402
import marmites_config as mcfg                       # noqa: E402

mdates = pytest.importorskip('matplotlib.dates')


def _hourly(first, last, rain_by_day):
    """An hourly record (matplotlib day numbers) from ``first`` to ``last``
    (dates), each day's rain spread over its 24 hours."""
    n = (last - first).days + 1
    t0 = mdates.date2num(datetime(first.year, first.month, first.day))
    t = t0 + np.arange(n * 24) / 24.0
    r = np.zeros(n * 24)
    for d, mm in rain_by_day.items():
        k = (d - first).days
        r[k * 24:(k + 1) * 24] = mm / 24.0
    return t, r


def _yday(d):
    return d.timetuple().tm_yday


def test_the_wet_season_starts_with_the_first_rain_of_the_threshold():
    """A 15 mm day on 31 August and a 9.9 mm day in September do not start
    it; the 12 mm of 20 October does, its green-up lasting the transition;
    a year without such a day keeps the fixed day; every step of a year
    carries its year's day."""
    rain = {date(2009, 8, 31): 15.0, date(2009, 9, 10): 9.9,
            date(2009, 10, 20): 12.0, date(2010, 9, 5): 30.0,
            date(2011, 9, 3): 5.0}
    t, r = _hourly(date(2009, 6, 1), date(2011, 12, 31), rain)
    jw, rows = RS.wet_season_start(t, r, 10.0, 9, 1, 300, 20)
    assert rows[0] == (2009, date(2009, 10, 20), pytest.approx(12.0),
                       _yday(date(2009, 10, 20)) + 20)
    assert rows[1][1] == date(2010, 9, 5)
    assert rows[2] == (2011, None, rows[2][2], 300.0)
    assert np.isnan(rows[2][2])
    years = np.array([mdates.num2date(x).year for x in t])
    assert set(jw[years == 2009]) == {_yday(date(2009, 10, 20)) + 20}
    assert set(jw[years == 2011]) == {300.0}


def test_the_rain_of_a_day_is_its_total_and_zones_are_averaged():
    """10 mm spread over a day's hours is a 10 mm day; with two zones,
    20 mm in one and none in the other is a 10 mm day too."""
    t, r = _hourly(date(2009, 8, 1), date(2009, 12, 31),
                   {date(2009, 9, 15): 10.0})
    _jw, rows = RS.wet_season_start(t, r, 10.0, 9, 1, 300, 20)
    assert rows[0][1] == date(2009, 9, 15)
    _jw, rows = RS.wet_season_start(t, np.vstack([2 * r, 0 * r]), 10.0, 9, 1,
                                    300, 20)
    assert rows[0][1] == date(2009, 9, 15)


def test_a_record_starting_after_the_date_keeps_the_fixed_day():
    """Rain before the record is unknown: the first 'first rain' cannot be
    told, so that year keeps j_wet (and j_dry)."""
    t, r = _hourly(date(2009, 10, 1), date(2009, 12, 31),
                   {date(2009, 10, 5): 25.0})
    _jw, rows = RS.wet_season_start(t, r, 10.0, 9, 1, 300, 20)
    assert rows == [(2009, None, rows[0][2], 300.0)]
    _jd, rows = RS.dry_season_start(t, r, 10.0, 30, 5, 1, 150, 30)
    assert rows[0][1] is None and rows[0][3] == 150.0


def test_the_grass_is_wilted_after_thirty_days_under_ten_mm():
    """Rain every 4 days until 20 May (3 mm each), then dry: the first 30
    days from 1 May totalling under 10 mm run from 9 May (the rains of 12,
    16 and 20 May, 9 mm) to 7 June, when the grass is wilted; the wet-to-dry
    transition starts 30 days before. The rains of March-April do not
    count."""
    first = date(2010, 3, 1)
    rain = {first + timedelta(days=k): 3.0 for k in range(0, 81, 4)}
    t, r = _hourly(first, date(2010, 12, 31), rain)
    last_rain = max(rain)
    assert last_rain == date(2010, 5, 20)
    jd, rows = RS.dry_season_start(t, r, 10.0, 30, 5, 1, 150, 30)
    (y, end, mm, j_d), = rows
    # the window of 30 days whose total drops under 10 mm: it may hold
    # the last three rains (9 mm) at most
    assert y == 2010 and mm < 10.0
    win = [d for d in rain if end - timedelta(days=29) <= d <= end]
    assert sum(rain[d] for d in win) == pytest.approx(mm)
    assert len(win) == 3 and end == date(2010, 6, 7)
    assert j_d == _yday(end) - 30
    assert set(jd) == {_yday(end) - 30}


def test_a_wet_spring_keeps_the_fixed_dry_day():
    first = date(2010, 3, 1)
    rain = {first + timedelta(days=k): 5.0 for k in range(0, 306, 7)}
    t, r = _hourly(first, date(2010, 12, 31), rain)
    _jd, rows = RS.dry_season_start(t, r, 10.0, 30, 5, 1, 150, 30)
    assert rows[0][1] is None and rows[0][3] == 150.0


def test_the_configuration_checks_the_dates_and_amounts():
    v = mcfg.Vegetation()
    assert (v.wet_rain_mm, v.dry_rain_mm) == (0.0, 0.0), 'off by default'
    assert (v.wet_rain_from, v.dry_rain_from, v.dry_rain_days) == (
        '09-01', '05-01', 30)
    assert mcfg.month_day('09-01') == (9, 1)
    assert mcfg.month_day('02-29') is None
    assert mcfg.month_day('Sept 1') is None
