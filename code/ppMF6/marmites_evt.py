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
    types whose root tips are below the start-of-day head ``h0`` take part.

    Type v takes rate_v * clip((h - tip_v) / ramp, 0, 1): its full rate with
    the head ``ramp`` or more above its root tip, nothing at the tip, linear
    in between -- a ramp ANCHORED AT THE TIP, so the curve is never steeper
    than rate_v / ramp. Flat above h0, at the rate there. Returns (SURFACE,
    RATE, DEPTH, PXDP, PETM).

    Until 2026-10-16 the ramp started at h0, min(ramp, h0 - tip_v) wide:
    the full rate every morning, but a ramp that shrank to nothing as the
    head neared a tip. Where summer water tables sank to the 15 m Q. ilex
    tips that was a step MF6's Newton solve went back and forth across --
    ~200 outer iterations a day, and July-September 2009 of one PEST
    realisation took 2 h 37 min instead of 12 min."""
    r = np.asarray(rates, dtype=float).ravel()
    z = np.asarray(tips, dtype=float).ravel()
    on = (r > 0.0) & (z < float(h0))
    if not on.any():
        return _flat(h0, nseg)
    r, z = r[on], z[on]

    def T(h):
        return float(np.sum(r * np.clip((h - z) / ramp, 0.0, 1.0)))

    R = T(float(h0))
    zmin = float(z.min())
    D = float(h0) - zmin
    hs = np.unique(np.r_[z, np.minimum(z + ramp, float(h0))])
    x = (float(h0) - hs) / D
    y = np.array([T(h) for h in hs]) / R
    x, y = _pad(x, y, nseg)
    return (float(h0), R, D, x, y)


# ---------------------------------------------------------------- #
# ALL CELLS AT ONCE. The coupler builds both curves for every cell every
# day: 4524 cells on La Mata, and in Python that cost as much as the soil
# column itself -- 3.5 s a day under the profiler, 43 % of a day's work
# (2026-10-11), the per-cell np.insert of _pad first. The same arithmetic,
# compiled with numba over the whole cell list; eg_curve / tg_curve above
# stay the reference the tests hold it to. Without numba, curves() loops
# over them.
# ---------------------------------------------------------------- #

try:                                                    # pragma: no cover
    import numba as _numba
    _jit = _numba.njit(cache=True)
except ImportError:                                     # pragma: no cover
    _numba = None
    _jit = None


def _pad_nb(x, y, nseg, ox, oy):
    """_pad into ``ox``/``oy`` (length nseg - 1); returns 0, or the number
    of breakpoints when there are too many for nseg."""
    m = 0
    xs = np.empty(x.size)
    ys = np.empty(x.size)
    for i in range(x.size):
        if x[i] > _EPS and x[i] < 1.0 - _EPS:
            xs[m] = x[i]
            ys[m] = y[i]
            m += 1
    o = np.argsort(xs[:m], kind='mergesort')
    xk = np.empty(m)
    yk = np.empty(m)
    k = 0
    for i in range(m):
        xi = xs[o[i]]
        if k == 0 or xi - xk[k - 1] > _EPS:
            xk[k] = xi
            yk[k] = ys[o[i]]
            k += 1
    if k > nseg - 1:
        return k
    for i in range(k):
        ox[i] = xk[i]
        oy[i] = yk[i]
    n = k
    while n < nseg - 1:
        # the longest segment of (0, 1), (x, y), (1, 0), split in two
        best, kk, prev = -1.0, 0, 0.0
        for i in range(n + 1):
            nxt = ox[i] if i < n else 1.0
            if nxt - prev > best:
                best, kk = nxt - prev, i
            prev = nxt
        x0 = ox[kk - 1] if kk > 0 else 0.0
        y0 = oy[kk - 1] if kk > 0 else 1.0
        x1 = ox[kk] if kk < n else 1.0
        y1 = oy[kk] if kk < n else 0.0
        for i in range(n, kk, -1):
            ox[i] = ox[i - 1]
            oy[i] = oy[i - 1]
        ox[kk] = 0.5 * (x0 + x1)
        oy[kk] = 0.5 * (y0 + y1)
        n += 1
    lo = 1.0
    for i in range(nseg - 1):
        v = min(max(oy[i], 0.0), 1.0)
        lo = min(lo, v)
        oy[i] = lo
    return 0


def _shah_nb(d, y0, b, dll, ext):
    d = d * 100.0
    if d >= ext:
        return 0.0
    if d <= dll:
        return 1.0
    return min(y0 + np.exp(-b * (d - dll)), 1.0)


def _flat_nb(nseg, ox, oy):
    step = 1.0 / nseg
    for i in range(nseg - 1):
        ox[i] = (i + 1) * step
        oy[i] = 0.0


def _eg_nb(pe, land, h0, y0, b, dll_cm, ext_cm, nseg, last, ox, oy, out):
    """eg_curve into out = [SURFACE, RATE, DEPTH]; returns as _pad_nb."""
    d0 = land - h0
    ext = ext_cm / 100.0
    dll = dll_cm / 100.0
    f0 = _shah_nb(max(d0, 0.0), y0, b, dll_cm, ext_cm)
    if pe <= 0.0 or f0 <= 0.0 or d0 >= ext:
        _flat_nb(nseg, ox, oy)
        out[0], out[1], out[2] = h0, 0.0, 1.0
        return 0
    surface = h0 if d0 >= 0.0 else land
    d0 = max(d0, 0.0)
    D = ext - d0
    d_end = ext - last * D
    pts = np.empty(nseg + 1)
    m = 0
    if d0 < dll and dll < d_end:
        pts[m] = dll
        m += 1
    start = max(d0, dll)
    n_in = max(nseg - 2 - m, 0)
    if n_in > 0 and d_end > start:
        bm = b * 100.0
        a = np.exp(-bm * (start - dll) / 2.0)
        z = np.exp(-bm * (d_end - dll) / 2.0)
        step = (z - a) / (n_in + 1)
        for i in range(1, n_in + 1):
            u = i * step + a
            pts[m] = dll - 2.0 * np.log(u) / bm
            m += 1
    pts[m] = d_end
    m += 1
    x = np.empty(m)
    y = np.empty(m)
    for i in range(m):
        d = min(max(pts[i], d0), ext)
        x[i] = (d - d0) / D
        y[i] = _shah_nb(d, y0, b, dll_cm, ext_cm) / f0
    out[0], out[1], out[2] = surface, pe * f0, D
    return _pad_nb(x, y, nseg, ox, oy)


def _tg_nb(rates, tips, h0, ramp, nseg, ox, oy, out):
    """tg_curve into out = [SURFACE, RATE, DEPTH]; returns as _pad_nb."""
    nv = rates.size
    r = np.empty(nv)
    z = np.empty(nv)
    m = 0
    for v in range(nv):
        if rates[v] > 0.0 and tips[v] < h0:
            r[m] = rates[v]
            z[m] = tips[v]
            m += 1
    if m == 0:
        _flat_nb(nseg, ox, oy)
        out[0], out[1], out[2] = h0, 0.0, 1.0
        return 0
    r = r[:m]
    z = z[:m]
    top = np.empty(m)              # each ramp's top, never above h0
    R = 0.0
    zmin = z[0]
    for v in range(m):
        top[v] = min(z[v] + ramp, h0)
        R += r[v] * min(max((h0 - z[v]) / ramp, 0.0), 1.0)
        zmin = min(zmin, z[v])
    D = h0 - zmin
    hs = np.unique(np.concatenate((z, top)))
    x = np.empty(hs.size)
    y = np.empty(hs.size)
    for i in range(hs.size):
        T = 0.0
        for v in range(m):
            T += r[v] * min(max((hs[i] - z[v]) / ramp, 0.0), 1.0)
        x[i] = (h0 - hs[i]) / D
        y[i] = T / R
    out[0], out[1], out[2] = h0, R, D
    return _pad_nb(x, y, nseg, ox, oy)


def _curves_nb(eg_pe, land, h0, shah, tg_rate, tg_tip, ramp, nseg, last,
               S, RT, DP, PX, PM):
    """Both curves for every cell: S, RT, DP (2, n) and PX, PM
    (2, n, nseg - 1), index 0 Eg and 1 Tg. Returns -1, or the first cell
    with too many breakpoints for nseg."""
    out = np.empty(3)
    for k in range(eg_pe.size):
        bad = _eg_nb(eg_pe[k], land[k], h0[k], shah[k, 0], shah[k, 1],
                     shah[k, 2], shah[k, 3], nseg, last, PX[0, k], PM[0, k],
                     out)
        S[0, k], RT[0, k], DP[0, k] = out[0], out[1], out[2]
        if bad:
            return k
        bad = _tg_nb(tg_rate[k], tg_tip[k], h0[k], ramp, nseg, PX[1, k],
                     PM[1, k], out)
        S[1, k], RT[1, k], DP[1, k] = out[0], out[1], out[2]
        if bad:
            return k
    return -1


if _jit is not None:
    _pad_nb = _jit(_pad_nb)
    _shah_nb = _jit(_shah_nb)
    _flat_nb = _jit(_flat_nb)
    _eg_nb = _jit(_eg_nb)
    _tg_nb = _jit(_tg_nb)
    _curves_nb = _jit(_curves_nb)


def curves(eg_pe, land, h0, shah, tg_rate, tg_tip, ramp, nseg, last=0.05,
           compiled=None):
    """The day's Eg and Tg curves for EVERY cell: ``{'eg'|'tg': {SURFACE,
    RATE, DEPTH, PXDP, PETM}}``, as eg_curve and tg_curve give them one cell
    at a time. ``shah`` holds Shah's parameters per cell: a list of dicts,
    or an (n, 4) array of y0, b, dll, ext_d (the compiled soil step's).
    ``compiled``: None uses numba when it is installed; False forces the
    per-cell Python (the reference)."""
    n = int(np.size(eg_pe))
    ns = int(nseg)
    if compiled is None:
        compiled = _jit is not None
    if isinstance(shah, np.ndarray):
        sp = np.ascontiguousarray(shah, dtype=np.float64).reshape(n, 4)
        shah = [{'y0': a, 'b': b, 'dll': c, 'ext_d': d}
                for a, b, c, d in sp.tolist()] if not compiled else None
    else:
        sp = None
    out = {t: {'SURFACE': np.zeros(n), 'RATE': np.zeros(n),
               'DEPTH': np.ones(n), 'PXDP': np.zeros((n, ns - 1)),
               'PETM': np.zeros((n, ns - 1))} for t in ('eg', 'tg')}
    if not compiled:
        for k in range(n):
            for t, (s, r, d, x, y) in (
                    ('eg', eg_curve(eg_pe[k], land[k], h0[k], shah[k], ns,
                                    last=last)),
                    ('tg', tg_curve(tg_rate[k], tg_tip[k], h0[k], ramp,
                                    ns))):
                c = out[t]
                c['SURFACE'][k], c['RATE'][k], c['DEPTH'][k] = s, r, d
                c['PXDP'][k], c['PETM'][k] = x, y
        return out
    if sp is None:
        sp = np.array([[p['y0'], p['b'], p['dll'], p['ext_d']] for p in shah],
                      dtype=np.float64).reshape(n, 4)
    S, RT, DP = np.zeros((2, n)), np.zeros((2, n)), np.ones((2, n))
    PX, PM = np.zeros((2, n, ns - 1)), np.zeros((2, n, ns - 1))
    bad = _curves_nb(np.asarray(eg_pe, dtype=np.float64).ravel(),
                     np.asarray(land, dtype=np.float64).ravel(),
                     np.asarray(h0, dtype=np.float64).ravel(), sp,
                     np.ascontiguousarray(tg_rate, dtype=np.float64),
                     np.ascontiguousarray(tg_tip, dtype=np.float64),
                     float(ramp), ns, float(last), S, RT, DP, PX, PM)
    if bad >= 0:
        # the reference raises the message; let it, for that cell
        p = dict(zip(('y0', 'b', 'dll', 'ext_d'), sp[bad].tolist()))
        eg_curve(eg_pe[bad], land[bad], h0[bad], p, ns, last=last)
        tg_curve(tg_rate[bad], tg_tip[bad], h0[bad], ramp, ns)
        raise ValueError('cell %d: too many breakpoints for et.evt_nseg = %d'
                         % (bad, ns))
    for i, t in enumerate(('eg', 'tg')):
        c = out[t]
        c['SURFACE'], c['RATE'], c['DEPTH'] = S[i], RT[i], DP[i]
        c['PXDP'], c['PETM'] = PX[i], PM[i]
    return out
