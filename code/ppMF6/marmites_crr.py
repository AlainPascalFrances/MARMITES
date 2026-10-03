# -*- coding: utf-8 -*-
"""Cascade routing and reinfiltration (CRR) of surface runoff.  WP5.

Daoud et al. (2022, Hydrogeology Journal 30:899, Eq. 23) send the water a
cell cannot infiltrate to its DOWNSLOPE neighbours, in proportion to the
slope towards each (multiple flow directions, Quinn et al. 1991):

    S_ij     = (elv_i - elv_j) / l_ij          l_ij: centre-to-centre distance
    alpha_ij = beta * S_ij / sum_j(S_ij)       over the neighbours with S_ij > 0

so the fractions of a cell sum to beta, and the rest, 1 - beta, evaporates on
the way (beta "allows for appropriate partitioning between evaporated water
and direct runoff"; Daoud calibrate it at 0.8-1). A cell with no lower
neighbour is a topographic SINK, where all the water evaporates.

What a cell receives depends on what the receiving cell is, in the priority
LAK > SFR > soil:

  * a pond cell (LAK)  -- the water goes to the lake (LAK RUNOFF);
  * a channel cell (SFR) -- to its reach (SFR INFLOW);
  * any other cell     -- onto the MARMITES soil column, where it infiltrates
    by the same law as the rain (Eq. 1b, I = min[Ssurf, D1(phi1 - theta1)])
    and what the soil cannot take runs on downslope.

THE RECEIVER IS THE MARMITES SOIL COLUMN, NOT UZF. CdL does this with MF6's
mover, which can only hand the water to a UZF cell -- beneath the soil
column, skipping its storage, its evaporation and its infiltration law
(cookbook D1). So the cascade runs here, in Python, inside the soil model's
stress period. Not to be confused with the rejected infiltration MF6 hands
back, which enters the column from BELOW and is untouched (cookbook D4).

THE ORDER. Runoff is no longer local: a cell must be solved after every cell
that drains onto it. That is the descending-elevation order, a permutation
computed once (``order``) -- the cell list itself never moves, the coupler
and the DISV mapping both depend on it. A receiver is strictly lower than its
donor, so in that order every receiver comes after all of its donors; this is
asserted when the network is built.

Open-water cells (LAK, SFR) are never donors: their own runoff -- the rain on
the open fraction and what the soil fraction sheds -- goes to their own water
body, as it did before the cascade existed.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

from types import SimpleNamespace

import numpy as np

__all__ = ['CascadeNetwork', 'CRRError', 'KIND_SOIL', 'KIND_SFR', 'KIND_LAK',
           'SINKS', 'network_from_topology', 'structured_topology']

KIND_SOIL, KIND_SFR, KIND_LAK = 0, 1, 2
KIND_NAMES = ('soil', 'stream', 'pond')
SINKS = ('evaporate', 'route')


class CRRError(Exception):
    """The cascade cannot be built, or it lost water."""


class CascadeNetwork:
    """The downslope cascade over the MARMITES cell list.

    Every array is indexed by the position in the cell list (= the cell id).

    Parameters
    ----------
    elev : (n,) land-surface elevation per cell [m]
    neighbours : list of int arrays, the face neighbours of each cell, in
        cell-list positions (inactive cells already left out)
    xy : (n, 2) cell centres [m], for l_ij
    kind : (n,) KIND_SOIL / KIND_SFR / KIND_LAK
    beta : the partitioning factor of Eq. 23, in (0, 1]
    sinks : 'evaporate' (Daoud) -- the water of a sink evaporates; or
        'route' -- it is sent to the nearest stream or pond cell, as if the
        depression spilled into it
    edge : (n,) bool or None, cells on the catchment edge (for the report:
        an edge "sink" may only be a cell whose way down leaves the catchment)

    Attributes
    ----------
    order : (n,) the descending-elevation permutation the cascade visits
    receivers, alpha : per cell, the downslope receivers and their fractions
        (empty for an open-water cell and for a sink)
    sink : (n,) bool, soil cells with no lower neighbour
    route_to : (n,) int, under sinks='route' the open cell a sink spills to
    """

    def __init__(self, elev, neighbours, xy, kind, beta=1.0, sinks='evaporate',
                 edge=None):
        elev = np.asarray(elev, dtype=float).ravel()
        n = elev.size
        kind = np.asarray(kind, dtype=int).ravel()
        xy = np.asarray(xy, dtype=float).reshape(n, 2)
        if kind.size != n or len(neighbours) != n:
            raise CRRError('elevation, kind and neighbours disagree on the '
                           'number of cells (%d, %d, %d)'
                           % (n, kind.size, len(neighbours)))
        if not np.all(np.isfinite(elev)):
            raise CRRError('%d cell(s) have no elevation: the cascade needs '
                           'the land surface everywhere'
                           % int(np.count_nonzero(~np.isfinite(elev))))
        beta = float(beta)
        if not 0.0 < beta <= 1.0:
            raise CRRError('beta must be in (0, 1], got %g' % beta)
        if sinks not in SINKS:
            raise CRRError('sinks must be one of %s, got %r' % (SINKS, sinks))
        self.n, self.elev, self.kind, self.xy = n, elev, kind, xy
        self.beta, self.sinks = beta, sinks
        self.edge = (np.zeros(n, dtype=bool) if edge is None
                     else np.asarray(edge, dtype=bool).ravel())

        self.receivers = [np.zeros(0, dtype=int)] * n
        self.alpha = [np.zeros(0, dtype=float)] * n
        self.sink = np.zeros(n, dtype=bool)
        for i in range(n):
            if kind[i] != KIND_SOIL:
                continue
            # once per neighbour: two cells may share more than one face
            nb = np.unique(np.asarray(neighbours[i], dtype=int))
            nb = nb[nb != i]
            if nb.size:
                d = np.hypot(*(xy[nb] - xy[i]).T)
                s = np.where(d > 0.0, (elev[i] - elev[nb]) / np.where(
                    d > 0.0, d, 1.0), 0.0)
                down = s > 0.0
                nb, s = nb[down], s[down]
            if not nb.size:
                self.sink[i] = True
                continue
            self.receivers[i] = nb
            self.alpha[i] = beta * s / s.sum()

        # where a sink spills under 'route': the nearest open-water cell
        self.route_to = np.full(n, -1, dtype=int)
        if sinks == 'route' and self.sink.any():
            water = np.flatnonzero(kind != KIND_SOIL)
            if water.size:
                for i in np.flatnonzero(self.sink):
                    d = np.hypot(*(xy[water] - xy[i]).T)
                    self.route_to[i] = int(water[int(np.argmin(d))])

        # each donor's receivers split once by what they are, so a stress
        # period only multiplies: soil columns get run-on, open water a
        # delivery (cells are visited ~10^4 times a period)
        self._split = []
        for i in range(n):
            r, a = self.receivers[i], self.alpha[i]
            s = kind[r] == KIND_SOIL
            self._split.append((r[s], a[s], r[~s], a[~s]))

        # descending elevation; a stable sort keeps ties in list order (two
        # cells at the same height never exchange water: S_ij = 0)
        self.order = np.argsort(-elev, kind='stable')
        self._check_order()

    def _check_order(self):
        """The order is a permutation of the cell list, and every receiver is
        visited after its donor (cookbook 5.2)."""
        if not np.array_equal(np.sort(self.order), np.arange(self.n)):
            raise CRRError('the cascade order is not a permutation of the '
                           'cell list')
        rank = np.empty(self.n, dtype=int)
        rank[self.order] = np.arange(self.n)
        for i in range(self.n):
            r = self.receivers[i]
            if r.size and np.any(rank[r] <= rank[i]):
                raise CRRError('cell %d drains onto a cell visited before it'
                               % i)

    # ------------------------------------------------------------------ #

    def route(self, area, cell_fn, conv=1000.0):
        """Run one stress period of the cascade.

        ``cell_fn(k, runon)`` solves cell ``k`` given the run-on it receives,
        ``runon`` in mm/d per unit of CELL area, and returns the runoff that
        leaves it, same units. Cells are visited in ``order``.

        Returns a namespace of per-cell volumes [m3/d]:
          ro      -- runoff leaving each cell (its own and what it passed on)
          runon   -- run-on received by each soil column
          deliver -- what reaches the water body of each open cell (its own
                     runoff, the cascade's, and a routed sink's)
          ecrr    -- evaporated by the cascade, booked at the donor: the
                     1 - beta share, and a sink's water under 'evaporate'
          ecrr_sink -- the sink part of ecrr
          routed  -- sent from a sink to an open cell under 'route' (at the
                     sink; it is also in ``deliver`` at the receiving cell)
        Every cubic metre of ``ro`` ends in exactly one of runon, deliver or
        ecrr; ``check`` asserts it.
        """
        area = np.asarray(area, dtype=float).ravel()
        n = self.n
        ro = np.zeros(n)
        runon = np.zeros(n)
        deliver = np.zeros(n)
        ecrr = np.zeros(n)
        ecrr_sink = np.zeros(n)
        routed = np.zeros(n)
        conv = float(conv)
        kind, split, sink = self.kind, self._split, self.sink
        for k in self.order.tolist():
            r = float(cell_fn(k, runon[k] / area[k] * conv))
            q = max(r, 0.0) * area[k] / conv
            ro[k] = q
            if q == 0.0:
                continue
            if kind[k] != KIND_SOIL:
                deliver[k] += q
                continue
            if not sink[k]:
                rs, a_s, rw, a_w = split[k]
                vs, vw = q * a_s, q * a_w
                runon[rs] += vs              # receivers are distinct cells
                deliver[rw] += vw
                ecrr[k] += q - (vs.sum() + vw.sum())
            elif self.route_to[k] >= 0:
                v = q * self.beta
                deliver[self.route_to[k]] += v
                routed[k] += v
                ecrr[k] += q - v
            else:
                ecrr[k] += q
                ecrr_sink[k] += q
        res = SimpleNamespace(ro=ro, runon=runon, deliver=deliver, ecrr=ecrr,
                              ecrr_sink=ecrr_sink, routed=routed)
        self.check(res)
        return res

    @staticmethod
    def check(res, rtol=1e-9):
        """Mass balance of one pass: runoff out = run-on + delivered +
        evaporated, to machine precision (cookbook 5.7)."""
        out = float(res.ro.sum())
        dest = float(res.runon.sum() + res.deliver.sum() + res.ecrr.sum())
        if abs(out - dest) > rtol * max(abs(out), 1.0):
            raise CRRError('the cascade lost water: %.12g m3/d left the cells, '
                           '%.12g m3/d arrived' % (out, dest))

    # ------------------------------------------------------------------ #

    def summary(self, area=None):
        """One paragraph on the network, for the run log (cookbook 5.3).
        With the cell areas, also where a uniform runoff on the soil cells
        would end if nothing reinfiltrated -- how connected the hillslopes
        are to the streams, before the soil has its say."""
        soil = self.kind == KIND_SOIL
        nrec = np.array([r.size for r in self.receivers])
        donors = soil & ~self.sink
        to_water = np.array([bool(r.size) and bool(np.any(self.kind[r] != KIND_SOIL))
                             for r in self.receivers])
        n_sink = int(self.sink.sum())
        n_edge = int((self.sink & self.edge).sum())
        lines = ['CRR: %d cell(s): %d soil, %d stream, %d pond; beta %g; '
                 'sinks %s' % (self.n, int(soil.sum()),
                               int((self.kind == KIND_SFR).sum()),
                               int((self.kind == KIND_LAK).sum()),
                               self.beta, self.sinks),
                 '     %d soil cell(s) drain downslope, to %.2f neighbour(s) on '
                 'average; %d of them straight into a stream or a pond'
                 % (int(donors.sum()),
                    float(nrec[donors].mean()) if donors.any() else 0.0,
                    int((donors & to_water).sum())),
                 '     %d topographic sink(s), %d of them on the catchment '
                 'edge; their water %s'
                 % (n_sink, n_edge,
                    'evaporates' if self.sinks == 'evaporate' else
                    'goes to the nearest stream or pond cell'
                    + ('' if (self.route_to[self.sink] >= 0).all()
                       else ' (none exists: it evaporates)'))]
        if area is not None and soil.any():
            res = self.route(area, lambda k, runon: runon + (
                1.0 if soil[k] else 0.0))
            gen = float(np.sum(np.asarray(area, float)[soil])) / 1000.0
            sfr = float(res.deliver[self.kind == KIND_SFR].sum())
            lak = float(res.deliver[self.kind == KIND_LAK].sum())
            lines.append('     with no reinfiltration, the soil cells\' runoff '
                         'would end %.1f %% in the streams, %.1f %% in the '
                         'ponds, %.1f %% evaporated (%.1f %% at sinks)'
                         % (100.0 * sfr / gen, 100.0 * lak / gen,
                            100.0 * res.ecrr.sum() / gen,
                            100.0 * res.ecrr_sink.sum() / gen))
        return lines


def network_from_topology(topo, icell, elev, kind, beta=1.0,
                          sinks='evaporate'):
    """Build the cascade on a face topology (``marmites_topology``).

    ``icell[k]`` is the topology cell (icell2d) of MARMITES cell ``k``.
    Topology cells that are not in the cell list are inactive: they receive
    nothing, and a cell next to one is on the catchment edge.
    """
    icell = np.asarray(icell, dtype=int).ravel()
    pos = {int(c): k for k, c in enumerate(icell)}
    if len(pos) != icell.size:
        raise CRRError('two MARMITES cells map onto the same grid cell')
    neighbours, edge = [], np.zeros(icell.size, dtype=bool)
    for k, c in enumerate(icell):
        nb = topo.neighbours[int(c)]
        mine = [pos[int(x)] for x in nb if int(x) in pos]
        edge[k] = bool(topo.boundary[int(c)]) or len(mine) < len(nb)
        neighbours.append(np.asarray(mine, dtype=int))
    return CascadeNetwork(elev, neighbours, topo.xy[icell], kind, beta=beta,
                          sinks=sinks, edge=edge)


def structured_topology(delr, delc):
    """The face topology of a structured grid: the 4 rook neighbours, as on
    a mesh -- MODFLOW connects cells only across faces, and so does the
    cascade. Cell (i, j) is topology cell ``i * ncol + j``."""
    from marmites_grid import disv_from_structured
    from marmites_topology import MeshTopology
    verts, cell2d, ncpl = disv_from_structured(delr, delc)
    return MeshTopology({'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl})
