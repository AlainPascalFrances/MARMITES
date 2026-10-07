# -*- coding: utf-8 -*-
"""Groundwater ET as MF6 EVT curves, evaluated at the SOLVED head -- the
coupled run's only groundwater-ET path since 2026-10-07.

MARMITES computes groundwater evaporation Eg (Shah et al. 2007) and
transpiration Tg (the kTg function of soil moisture, per vegetation type,
while the water table is above the root tip). Until 2026-10-07 the coupled
run could also hand the result over as a fixed WEL rate at the head of the
PREVIOUS day (et.gw_route = 'wel'): whatever the head did during the day the
well kept pumping, and MMsoil approximated the drawdown itself, with Sy.
That route is gone; the uncoupled MMsoil keeps its own Eg/Tg and drawdown.

Here MMsoil decides how much each source COULD take that day, from the
demand, soil moisture, roots and Shah curve, and MF6 takes it at the head it
solves for, through two EVT packages (one for Eg, one for Tg, so the budget
keeps them apart). Each day's curve:

  * starts at the START-OF-DAY head (SURFACE) with MARMITES' rate: above it
    the rate is flat, so within the day EVT can only fall, never exceed what
    was reserved -- and UZF's demand capped at what that leaves keeps
    ETsoil + ETuzf + ETg <= PE + PT exactly;
  * falls as the water table does: Eg along Shah's curve to its extinction
    depth, Tg type by type as the head passes each root tip, over a ramp of
    width ``ramp`` (a step would be a kink MF6's Newton solve stalls on).

All lengths in m, rates in m/d per unit of CELL area. The curves are
MF6's segmented form: NSEG segments, PXDP the depth fractions of DEPTH and
PETM the ET fractions of RATE at the NSEG - 1 inner points, both strictly
monotone.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import numpy as np

__all__ = ['shah_f', 'eg_curve', 'tg_curve', 'flat_curve', 'NSEG_MIN']

NSEG_MIN = 3
_EPS = 1e-6          # the smallest step between two points [fraction]


def shah_f(d, p):
    """Shah et al. (2007) Eq. 17 as MARMITES uses it: the fraction of the
    potential evaporation at a water-table depth ``d`` [m] below the land
    surface; ``p`` holds Shah's dll, y0, b, ext_d in cm (cm^-1 for b). At
    most 1, zero beyond the extinction depth."""
    d = np.asarray(d, dtype=float) * 100.0               # m -> cm
    y0, b, dll, ext = p['y0'], p['b'], p['dll'], p['ext_d']
    f = np.minimum(y0 + np.exp(-b * (d - dll)), 1.0)
    f = np.where(d <= dll, 1.0, f)
    return np.where(d >= ext, 0.0, f)


def _flat(surface, nseg):
    """A curve that takes nothing (rate 0): valid, monotone, inert."""
    x = np.linspace(0.0, 1.0, nseg + 1)[1:-1]
    return (float(surface), 0.0, 1.0, x, np.zeros(nseg - 1))


def flat_curve(rate, surface, bottom, nseg):
    """A STEADY period's groundwater ET: ``rate`` [m/d] wherever the water
    table stands above ``bottom`` [m], tapering to zero over the last
    segment above it; full above ``surface``.

    A steady first period (a cold start) is driven by the MEAN forcing --
    per-cell mean recharge and mean ETg from an earlier run -- and has no
    start-of-day head to anchor a day's curve to. The WEL route drew that
    mean as a fixed rate with AUTO_FLOW_REDUCE; this is the same rate,
    reduced as the cell dries, through the EVT package the transient days
    use (2026-10-07). Returns (SURFACE, RATE, DEPTH, PXDP, PETM)."""
    depth = float(surface) - float(bottom)
    if not rate > 0.0 or not depth > 0.0:
        return _flat(surface, nseg)
    x = np.linspace(0.0, 1.0, nseg + 1)[1:-1]
    return (float(surface), float(rate), depth, x, np.ones(nseg - 1))


def _pad(x, y, nseg):
    """Exactly ``nseg - 1`` inner points on the curve through (0, 1),
    (x, y), (1, 0): too few are filled in along the longest segments (the
    curve does not change), never too many (the caller sizes nseg)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = (x > _EPS) & (x < 1.0 - _EPS)
    x, y = x[keep], y[keep]
    o = np.argsort(x)
    x, y = x[o], y[o]
    # no inner point at all is a valid curve -- one type whose root tip is
    # within its ramp of the start-of-day head: a straight line (2026-10-05)
    u = np.ones(x.size, dtype=bool)
    u[1:] = np.diff(x) > _EPS
    x, y = x[u], y[u]
    if x.size > nseg - 1:
        raise ValueError('%d breakpoints need nseg >= %d (et.evt_nseg)'
                         % (x.size, x.size + 1))
    while x.size < nseg - 1:
        X = np.r_[0.0, x, 1.0]
        Y = np.r_[1.0, y, 0.0]
        k = int(np.argmax(np.diff(X)))
        xm = 0.5 * (X[k] + X[k + 1])
        ym = 0.5 * (Y[k] + Y[k + 1])
        x = np.insert(x, k, xm)
        y = np.insert(y, k, ym)
    # strictly decreasing PETM is not required, non-increasing is: clip
    # the rounding of an interpolated midpoint
    y = np.minimum.accumulate(np.clip(y, 0.0, 1.0))
    return x, y


def eg_curve(pe, land, h0, p, nseg, last=0.05):
    """Eg: ``pe`` the potential [m/d] (what the soil and the deep
    unsaturated zone left of PE), ``land`` the land surface and ``h0`` the
    start-of-day head [m], ``p`` Shah's parameters for the soil.

    RATE = pe * f(d0), d0 = land - h0, flat above h0; below, f(d) / f(d0)
    to the extinction depth. Shah's curve stops at y0 > 0 there; it is
    brought to zero over the last ``last`` of the depth instead of with a
    step. Returns (SURFACE, RATE, DEPTH, PXDP, PETM)."""
    d0 = float(land) - float(h0)
    ext = float(p['ext_d']) / 100.0
    dll = float(p['dll']) / 100.0
    f0 = float(shah_f(max(d0, 0.0), p))
    if pe <= 0.0 or f0 <= 0.0 or d0 >= ext:
        return _flat(h0, nseg)
    # a head above the land surface evaporates at the full rate already:
    # the curve is measured from the land then
    surface = float(h0) if d0 >= 0.0 else float(land)
    d0 = max(d0, 0.0)
    D = ext - d0
    d_end = ext - last * D            # where the curve starts its last ramp
    pts = []
    if d0 < dll < d_end:
        pts.append(dll)
    start = max(d0, dll)
    n_in = max(nseg - 2 - len(pts), 0)
    if n_in and d_end > start:
        # the knots that keep a piecewise-linear exponential closest
        # everywhere equidistribute sqrt|f''| -- equal steps of
        # exp(-b (d - dll) / 2): close together where the curve bends
        b = float(p['b']) * 100.0                          # 1/cm -> 1/m
        u = np.linspace(np.exp(-b * (start - dll) / 2.0),
                        np.exp(-b * (d_end - dll) / 2.0), n_in + 2)[1:-1]
        pts.extend((dll - 2.0 * np.log(u) / b).tolist())
    pts.append(d_end)
    d = np.clip(np.asarray(pts, dtype=float), d0, ext)
    x = (d - d0) / D
    y = shah_f(d, p) / f0
    x, y = _pad(x, y, nseg)
    return (surface, float(pe) * f0, float(D), x, y)


def tg_curve(rates, tips, h0, ramp, nseg):
    """Tg: ``rates`` [m/d] and root ``tips`` [m] per vegetation type; the
    types whose roots reach the start-of-day head ``h0`` take their rate.

    As the head falls, type v stops over a ramp of width min(ramp,
    h0 - tip_v) above its root tip -- so it starts the day at its full rate
    and is off at its tip. Flat above h0. Returns (SURFACE, RATE, DEPTH,
    PXDP, PETM)."""
    r = np.asarray(rates, dtype=float).ravel()
    z = np.asarray(tips, dtype=float).ravel()
    on = (r > 0.0) & (z < float(h0))
    if not on.any():
        return _flat(h0, nseg)
    r, z = r[on], z[on]
    w = np.minimum(float(ramp), float(h0) - z)
    R = float(r.sum())
    zmin = float(z.min())
    D = float(h0) - zmin

    def T(h):
        return float(np.sum(r * np.clip((h - z) / w, 0.0, 1.0)))

    hs = np.unique(np.r_[z, z + w])
    x = (float(h0) - hs) / D
    y = np.array([T(h) for h in hs]) / R
    x, y = _pad(x, y, nseg)
    return (float(h0), R, D, x, y)
