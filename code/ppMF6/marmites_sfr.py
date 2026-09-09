# -*- coding: utf-8 -*-
"""Build a MODFLOW 6 SFR network for MARMITES from the PONDw channel map.

MARMITES already carries the stream network: ``inputSTREAMw.asc`` gives a channel
width for every cell the drainage net passes through (La Mata: 244 cells,
1.5-3.0 m wide) and ``inputSTREAMhmax.asc`` the channel depth (1.0-1.5 m). One
SFR reach per stream cell, routed on the sink-filled DEM the flow model already
uses as its top elevation.

Routing (``stream_network``)
---------------------------
Greedy steepest descent restricted to stream cells does NOT work here: the
channel crosses flat stretches where every neighbour has the same DEM value, so
descent terminates in local sinks and whole branches never reach an outlet (on
La Mata it stranded 211 of 244 cells). Instead the drainage tree is grown
*upslope from the known outlets* with a priority flood: the frontier is always
expanded at its lowest cell, and each newly reached cell takes the cell it was
reached from as its receiver. Every stream cell therefore gets exactly one
receiver and a guaranteed path to an outlet, flats included.

Bed elevations
--------------
DEM noise and flats leave a few receivers sitting *higher* than their own cell.
Bed tops are therefore smoothed to be downstream-monotonic (the SFRmaker rule,
Leaf et al. 2021 fig. 2): no reach bed is higher than the lowest bed upstream
of it.

MARMITES coupling
-----------------
``EVAPORATION`` is left at zero: MARMITES computes open-water evaporation
itself (E_ow in Eq. 1), so letting SFR evaporate as well would double-count it.
Runoff is not set here either; the coupler injects the MARMITES runoff of each
stress period into the reaches as INFLOW through the API.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import heapq

import numpy as np

__all__ = ['stream_network', 'build_sfr', 'SFRNetwork']

# defaults (overridable by build_sfr kwargs)
SFR_RBTH = 0.5        # m,   streambed thickness
SFR_RHK = 0.1         # m/d, streambed hydraulic conductivity
SFR_MAN = 0.035       # Manning's n
SFR_MINSLOPE = 1e-4   # minimum reach gradient
SFR_MINLEN = 0.1      # m,   floor on reach length
SFR_MINWID = 0.5      # m,   floor on channel width

_NB8 = ((-1, 0), (1, 0), (0, -1), (0, 1),
        (-1, -1), (-1, 1), (1, -1), (1, 1))


class SFRNetwork(object):
    """Result of :func:`stream_network` / :func:`build_sfr`."""

    def __init__(self, cells, recv, order, acc, outlets):
        self.cells = cells            # [(i, j), ...] in reach order (headwater first)
        self.recv = recv              # (i, j) -> receiving (i, j) or None at an outlet
        self.order = order            # priority-flood pop order (outlet first)
        self.acc = acc                # (i, j) -> number of stream cells draining through
        self.outlets = outlets        # [(i, j), ...]
        self.rno = {c: n for n, c in enumerate(cells)}
        self.packagedata = None
        self.connectiondata = None
        self.reach_len = None
        self.reach_wid = None
        self.reach_top = None
        self.reach_slope = None
        self.nmono = 0

    @property
    def nreaches(self):
        return len(self.cells)

    def __repr__(self):
        return ('<SFRNetwork %d reaches, %d outlet(s), max accumulation %d>'
                % (self.nreaches, len(self.outlets),
                   max(self.acc.values()) if self.acc else 0))


def stream_network(pondw, dem, outlets=None, drn_cells=None):
    """Route the stream cells of ``pondw`` on ``dem``.

    Parameters
    ----------
    pondw : 2-D array
        Channel width per cell; cells with ``pondw > 0`` are stream cells.
    dem : 2-D array
        Land-surface elevation (use the sink-filled DEM).
    outlets : sequence of (i, j), optional
        Cells where water leaves the model. Defaults to ``drn_cells`` that are
        also stream cells, else the single lowest stream cell.
    drn_cells : sequence of (i, j), optional
        Existing outlet DRN cells, used to infer the outlets.

    Returns
    -------
    SFRNetwork
    """
    pondw = np.asarray(pondw, dtype=float)
    dem = np.asarray(dem, dtype=float)
    if pondw.shape != dem.shape:
        raise ValueError('pondw %s and dem %s have different shapes'
                         % (pondw.shape, dem.shape))
    cells = [(int(i), int(j)) for i, j in zip(*np.where(pondw > 0))]
    if not cells:
        raise ValueError('no stream cells: pondw is zero everywhere')
    sset = set(cells)

    if outlets is None:
        cand = [tuple(int(v) for v in c) for c in (drn_cells or [])]
        outlets = [c for c in cand if c in sset]
        if not outlets:
            outlets = [min(cells, key=lambda c: dem[c])]
    else:
        outlets = [tuple(int(v) for v in c) for c in outlets]
        bad = [c for c in outlets if c not in sset]
        if bad:
            raise ValueError('outlet cells are not stream cells: %s' % (bad,))
    # A DRN cell list carries one record per layer, so the same (i, j) arrives
    # several times. Duplicates would seed the flood twice and produce
    # duplicate reaches, so collapse them (order-preserving).
    outlets = list(dict.fromkeys(outlets))

    # ---- priority flood upslope from the outlets ----------------------- #
    recv = {}
    seen = set(outlets)
    heap = []
    for n, c in enumerate(outlets):
        recv[c] = None
        # n breaks ties deterministically without comparing tuples
        heapq.heappush(heap, (float(dem[c]), n, c))
    counter = len(outlets)
    order = []
    while heap:
        z, _, c = heapq.heappop(heap)
        order.append(c)
        for di, dj in _NB8:
            n_ = (c[0] + di, c[1] + dj)
            if n_ in sset and n_ not in seen:
                seen.add(n_)
                recv[n_] = c
                # carry the running maximum so a flat cannot re-order the flood
                heapq.heappush(heap, (max(float(dem[n_]), z), counter, n_))
                counter += 1
    missing = [c for c in cells if c not in seen]
    if missing:
        raise ValueError(
            '%d stream cell(s) do not connect to any outlet, e.g. %s. The '
            'channel map has a disconnected component; either add an outlet '
            'there or fix inputSTREAMw.' % (len(missing), missing[:5]))

    # accumulation: number of stream cells draining through each cell
    acc = {c: 1 for c in cells}
    for c in reversed(order):
        r = recv[c]
        if r is not None:
            acc[r] += acc[c]

    # Reach numbering: reverse of the flood order, so reach 0 is the last cell
    # popped -- necessarily a headwater (nothing can drain into it), which
    # matters because MF6 writes a downstream connection as -rno and -0 has no
    # sign.
    ordered = list(reversed(order))
    return SFRNetwork(ordered, recv, order, acc, outlets)


def build_sfr(net, dem, pondhmax=None, pondw=None, delr=None, delc=None,
              botm=None, idomain=None, incision=None,
              rbth=SFR_RBTH, rhk=SFR_RHK, man=SFR_MAN,
              minslope=SFR_MINSLOPE, cellid=None, verbose=True):
    """Fill ``net.packagedata`` / ``net.connectiondata`` for ModflowGwfsfr.

    ``cellid(k, i, j)`` maps a cell to the grid-appropriate cellid (DIS tuple or
    DISV pair); if omitted, ``(k, i, j)`` is used.
    """
    dem = np.asarray(dem, dtype=float)
    cells = net.cells
    n = len(cells)
    if cellid is None:
        def cellid(k, i, j):
            return (k, i, j)

    dx = float(np.mean(delr)) if delr is not None else 1.0
    dy = float(np.mean(delc)) if delc is not None else dx

    # ---- reach length: distance to the receiving cell ------------------ #
    rlen = []
    for c in cells:
        r = net.recv[c]
        if r is None:                       # outlet: half a cell to the edge
            rlen.append(max(0.5 * (dx + dy) * 0.5, SFR_MINLEN))
            continue
        di, dj = abs(r[0] - c[0]), abs(r[1] - c[1])
        rlen.append(max(float(np.hypot(di * dy, dj * dx)), SFR_MINLEN))

    # ---- width: straight from the PONDw channel map -------------------- #
    if pondw is not None:
        pw = np.asarray(pondw, dtype=float)
        rwid = [max(float(pw[c]), SFR_MINWID) for c in cells]
    else:
        rwid = [1.0] * n

    # ---- bed top: land surface incised by the channel depth ------------ #
    if incision is not None:
        inc = np.asarray(incision, dtype=float)
        depth = [float(inc[c]) if inc.ndim == 2 else float(inc) for c in cells]
    elif pondhmax is not None:
        ph = np.asarray(pondhmax, dtype=float)
        depth = [float(ph[c]) for c in cells]
    else:
        depth = [0.0] * n
    rtp = [float(dem[c]) - d for c, d in zip(cells, depth)]

    # ---- downstream-monotonic smoothing (SFRmaker rule) ---------------- #
    # walk the flood order (outlet first is wrong here -- we need upstream
    # first), so iterate reaches from the headwaters down each path
    top_of = {c: rtp[k] for k, c in enumerate(cells)}
    nmono = 0
    for c in net.order[::-1]:               # headwater -> outlet
        r = net.recv[c]
        if r is not None and top_of[r] > top_of[c]:
            top_of[r] = top_of[c]
            nmono += 1
    rtp = [top_of[c] for c in cells]

    # ---- slope --------------------------------------------------------- #
    rgrd = []
    for k, c in enumerate(cells):
        r = net.recv[c]
        if r is None:
            # outlet: reuse the gradient of the reach flowing into it, else the floor
            ups = [u for u in cells if net.recv[u] == c]
            rgrd.append(max(minslope, min((top_of[u] - top_of[c]) / rlen[net.rno[u]]
                                          for u in ups)) if ups else minslope)
        else:
            rgrd.append(max((top_of[c] - top_of[r]) / rlen[k], minslope))

    # ---- connectivity --------------------------------------------------- #
    ups_of = {c: [] for c in cells}
    for c in cells:
        r = net.recv[c]
        if r is not None:
            ups_of[r].append(c)
    connectiondata, packagedata = [], []
    nlay = 1 if botm is None else int(np.asarray(botm).shape[0])
    for k, c in enumerate(cells):
        conns = [net.rno[u] for u in ups_of[c]]
        r = net.recv[c]
        if r is not None:
            conns.append(-net.rno[r])       # downstream is written negative
        connectiondata.append([k] + conns)
        # place the reach in the topmost active layer deep enough to hold the
        # streambed: MF6 requires rtp - rbth to sit above the cell bottom
        klay = 0
        if botm is not None:
            bed_bot = rtp[k] - rbth
            bo = np.asarray(botm, dtype=float)
            while klay < nlay - 1 and (
                    bo[klay][c] >= bed_bot
                    or (idomain is not None and np.asarray(idomain)[klay][c] <= 0)):
                klay += 1
        packagedata.append([k, cellid(klay, c[0], c[1]), rlen[k], rwid[k],
                            rgrd[k], rtp[k], rbth, rhk, man, len(conns), 1.0, 0])

    net.packagedata = packagedata
    net.connectiondata = connectiondata
    net.reach_len, net.reach_wid = rlen, rwid
    net.reach_top, net.reach_slope = rtp, rgrd
    net.nmono = nmono
    if verbose:
        print('SFR: %d reaches from %d stream cells, %d outlet(s) %s'
              % (n, n, len(net.outlets), net.outlets))
        print('   width %.2f-%.2f m, length %.1f-%.1f m, slope %.2g-%.2g, '
              'bed %.1f-%.1f m'
              % (min(rwid), max(rwid), min(rlen), max(rlen),
                 min(rgrd), max(rgrd), min(rtp), max(rtp)))
        if nmono:
            print('   monotonic bed smoothing: %d reach top(s) lowered '
                  '(DEM flats/noise)' % nmono)
        nflat = sum(1 for g in rgrd if g <= minslope)
        if nflat:
            print('   %d reach(es) sit on the %.2g minimum-slope floor' % (nflat, minslope))
    return net
