# -*- coding: utf-8 -*-
"""Dry and wet seasons started by the rain (user, 2026-11-09).

A vegetation type's seasonal values (LAI, height, albedo, shelter) changed
on fixed julian days, j_dry and j_wet, the same every year. La Mata's annual
grass follows the rain instead: it greens up with the first autumn rains and
wilts when spring rain stops. The fixed day 300 left the first 25-125 mm of
autumn rain on bare soil, and the soil re-wetted weeks before the probes did
(scratch diagnosis of the 2026-10-28 run).

With a threshold set, each calendar year:

* the WET season's transition starts on the first day from a date (1
  September) whose rain reaches ``mm``; its wet values are reached the
  dry-to-wet transition later;
* the DRY season is reached -- the grass wilted -- at the end of the first
  ``ndays`` days in a row from a date (1 May) whose rain totals less than
  ``mm``; the wet-to-dry transition leads up to that day.

A year whose record holds no such day -- or starts after the date -- keeps
the fixed julian day. The rain is the daily total of the meteo record, the
mean over its zones.
"""

from datetime import date

import numpy as np

__all__ = ['wet_season_start', 'dry_season_start']


def _daily(datenum, rain):
    """The record's days (dates), their rain (mm) and each step's day."""
    import matplotlib.dates as mdates
    t = np.asarray(datenum, dtype=float)
    r = np.asarray(rain, dtype=float)
    if r.ndim == 2:
        r = r.mean(axis=0)
    days, inv = np.unique(np.floor(t), return_inverse=True)
    daily = np.bincount(inv, weights=r, minlength=days.size)
    return [mdates.num2date(d).date() for d in days], daily, inv


def _per_step(when, inv, by_year):
    years = np.array([d.year for d in when])[inv]
    return np.array([by_year[y][2] for y in years.tolist()], dtype=float)


def wet_season_start(datenum, rain, mm, month, day, j_wet, trans):
    """The julian day each time step's WET season is reached, and each
    year's start: ``(J_w (nt,), [(year, rain day or None, its rain mm,
    J_w)])``.

    ``datenum``: matplotlib day numbers of the record's steps (any step
    length; a day is ``floor(datenum)``). ``rain``: mm per step, (nt,) or
    (nzone, nt) -- zones are averaged. ``mm``: the rain of the day that
    starts the season. ``month``/``day``: the date the search starts each
    year. ``j_wet``: the fixed julian day of a year without such a day.
    ``trans``: the dry-to-wet transition, days.
    """
    when, daily, inv = _daily(datenum, rain)
    by_year = {}
    for y in sorted({d.year for d in when}):
        start = date(y, int(month), int(day))
        hit = None
        if when[0] <= start:
            hit = next((k for k, d in enumerate(when) if d.year == y and
                        d >= start and daily[k] >= mm), None)
        by_year[y] = ((None, float('nan'), float(j_wet)) if hit is None else
                      (when[hit], float(daily[hit]),
                       float(when[hit].timetuple().tm_yday + trans)))
    return (_per_step(when, inv, by_year),
            [(y,) + by_year[y] for y in sorted(by_year)])


def dry_season_start(datenum, rain, mm, ndays, month, day, j_dry, trans):
    """The julian day each time step's DRY season STARTS -- the wet-to-dry
    transition, which ends, the vegetation wilted, on the last day of the
    first ``ndays`` days in a row from ``month``/``day`` whose rain totals
    less than ``mm`` -- and each year's: ``(J_d (nt,), [(year, that last
    day or None, the spell's rain mm, J_d)])``. ``j_dry``: the fixed julian
    day of a year without such a spell; ``trans``: the wet-to-dry
    transition, days. Other arguments as wet_season_start.
    """
    when, daily, inv = _daily(datenum, rain)
    n = int(ndays)
    csum = np.concatenate([[0.0], np.cumsum(daily)])
    by_year = {}
    for y in sorted({d.year for d in when}):
        start = date(y, int(month), int(day))
        hit = None
        if when[0] <= start:
            for k, d in enumerate(when):
                end = k + n - 1
                if d.year != y or d < start or end >= len(when):
                    continue
                if when[end].year != y or (when[end] - d).days != n - 1:
                    continue                  # past the year, or a gap
                if csum[end + 1] - csum[k] < mm:
                    hit = end
                    break
        if hit is None:
            by_year[y] = (None, float('nan'), float(j_dry))
        else:
            by_year[y] = (when[hit], float(csum[hit + 1] - csum[hit + 1 - n]),
                          float(when[hit].timetuple().tm_yday - trans))
    return (_per_step(when, inv, by_year),
            [(y,) + by_year[y] for y in sorted(by_year)])
