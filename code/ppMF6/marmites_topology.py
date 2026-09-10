# -*- coding: utf-8 -*-
"""Shared-face neighbour topology for an unstructured (DISV) mesh.  WP1c.5.

On the structured grid "neighbour" is a fixed 8-cell stencil, ``_NB8``, and
``(i-1, j)`` is the cell to the north. On a mesh under the ``(ncpl, 1)``
convention (see ``marmites_mesh``) the row index IS the icell2d, so ``i-1``
is whatever cell the mesh generator happened to number one lower -- somewhere
else entirely. Anything that walks neighbours therefore needs the real
topology, not a stencil:

  * **SFR** (WP3) routes water downslope from cell to cell;
  * **CRR** (WP5) spreads rejected infiltration over the downslope neighbours,
    weighted by ``alpha_ij = beta * S_ij / sum(S_ij)`` (Daoud et al. 2022,
    Eq. 23), which needs the slope to each neighbour;
  * **LAK** (WP4) needs to know which cells touch a lake cell.

Two cells are neighbours when they SHARE A FACE. The faces come from the
``cell2d`` vertex lists: each consecutive pair of vertices is an edge, and an
edge belonging to two cells is an interior face. This works for a quadtree
too, because MODFLOW 6 DISV requires the coarse cell of a refinement step to
list the hanging node, so the split edges match on both sides -- confirmed on
La Mata, where a 25 m cell correctly reports six 12.5 m neighbours.

Vertices are matched by COORDINATE, not by id. flopy's ``VoronoiGrid`` emits
the same point under several different vertex ids -- 127 of them on La Mata's
989-cell mesh -- and matching by id then splits the graph into seven
components, six of which are single boundary cells with no neighbours at all.
One of those slivers was where both of La Mata's west-edge drains landed, so
this is not a cosmetic defect: it produces a drain that nothing can flow to.

A cell touching another only at a corner is NOT a neighbour. That is a real
difference from ``_NB8``, which includes the four diagonals: on a mesh there
is no flow through a point, and MODFLOW 6 itself connects cells only across
shared faces. The routing therefore agrees with the flow model.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import numpy as np

__all__ = ['MeshTopology', 'TopologyError']


class TopologyError(Exception):
    """Raised when a mesh's face topology is not usable."""


def _ring(rec):
    """Vertex ids of one cell2d record, as a closed ring with no repeats.

    flopy's producers differ on whether the first vertex is repeated at the
    end. A repeated vertex would create a zero-length edge (v, v), which then
    looks like a face shared by every cell that closes its ring the same way.
    """
    n = int(rec[3])
    ivs = [int(x) for x in rec[4:4 + n]]
    out = []
    for iv in ivs:
        if not out or iv != out[-1]:
            out.append(iv)
    while len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out


class MeshTopology:
    """Shared-face adjacency of a DISV mesh.

    Parameters
    ----------
    gridprops : dict
        DISV gridprops (``vertices``, ``cell2d``, ``ncpl``).

    Attributes
    ----------
    ncpl : int
    xy : (ncpl, 2)          cell centres, as the mesh reports them
    neighbours : list       ``neighbours[ic]`` -> int array of neighbour cells
    face_length : list      parallel to ``neighbours``, shared face length [m]
    boundary : (ncpl,) bool cells with at least one face on the domain edge
    """

    def __init__(self, gridprops, tol=1e-6):
        cell2d = gridprops['cell2d']
        self.ncpl = int(gridprops.get('ncpl', len(cell2d)))
        self._vxy = {int(v[0]): (float(v[1]), float(v[2]))
                     for v in gridprops['vertices']}
        self.xy = np.array([[float(r[1]), float(r[2])] for r in cell2d],
                           dtype=float)
        self._canon, self.duplicate_vertices = self._canonical_vertices(tol)

        # ---- edge -> cells ------------------------------------------------
        edges = {}
        for ic, rec in enumerate(cell2d):
            ring = [self._canon[v] for v in _ring(rec)]
            # collapse again: two ids that were distinct may now be one
            ring = [v for k, v in enumerate(ring) if v != ring[k - 1]]
            if len(ring) < 3:
                raise TopologyError('cell %d has %d distinct vertices'
                                    % (ic, len(ring)))
            for k in range(len(ring)):
                a, b = ring[k], ring[(k + 1) % len(ring)]
                edges.setdefault((a, b) if a < b else (b, a), []).append(ic)
        self.edges = edges

        bad = {e: cs for e, cs in edges.items() if len(cs) > 2}
        if bad:
            e, cs = next(iter(bad.items()))
            raise TopologyError(
                '%d face(s) are shared by more than two cells, e.g. vertices '
                '%s shared by cells %s. The mesh is not a valid DISV grid: a '
                'face separates exactly two cells, or one cell and the '
                'boundary.' % (len(bad), e, cs))

        # ---- adjacency ----------------------------------------------------
        nb = [[] for _ in range(self.ncpl)]
        fl = [[] for _ in range(self.ncpl)]
        self.boundary = np.zeros(self.ncpl, dtype=bool)
        for (a, b), cs in edges.items():
            if len(cs) == 1:
                self.boundary[cs[0]] = True
                continue
            (x1, y1), (x2, y2) = self._vxy[a], self._vxy[b]
            length = float(np.hypot(x2 - x1, y2 - y1))
            i, j = cs
            nb[i].append(j)
            fl[i].append(length)
            nb[j].append(i)
            fl[j].append(length)
        self.neighbours = [np.asarray(x, dtype=int) for x in nb]
        self.face_length = [np.asarray(x, dtype=float) for x in fl]

    def _canonical_vertices(self, tol):
        """Map every vertex id to one id per DISTINCT LOCATION.

        Coordinates are quantised onto a ``tol``-metre lattice (1 micron by
        default), which is far below any real mesh feature and far above
        floating-point noise on UTM coordinates. Returns ``(map, n_merged)``.
        """
        canon, first = {}, {}
        merged = 0
        for iv, (x, y) in self._vxy.items():
            key = (round(x / tol), round(y / tol))
            if key in first:
                canon[iv] = first[key]
                merged += 1
            else:
                first[key] = iv
                canon[iv] = iv
        return canon, merged

    def vertex_xy(self, iv):
        """Coordinates of a (possibly canonicalised) vertex id."""
        return self._vxy[self._canon.get(int(iv), int(iv))]

    # ------------------------------------------------------------------ #
    @property
    def nfaces(self):
        return len(self.edges)

    @property
    def ninterior(self):
        return int(sum(1 for cs in self.edges.values() if len(cs) == 2))

    def distance(self, ic, jc):
        """Centre-to-centre distance between two cells [m]."""
        d = self.xy[int(jc)] - self.xy[int(ic)]
        return float(np.hypot(d[0], d[1]))

    def slopes(self, ic, z):
        """Downslope gradients from ``ic`` to each of its neighbours.

        ``z`` is a per-cell elevation, shaped ``(ncpl,)`` or ``(ncpl, 1)``.
        Returns ``(neighbours, slope)`` where slope is
        ``(z[ic] - z[jc]) / distance`` -- positive DOWNHILL, zero or negative
        uphill. This is the ``S_ij`` of the CRR multiple-flow-direction
        weights (WP5).
        """
        z = np.asarray(z, dtype=float).reshape(-1)
        nb = self.neighbours[int(ic)]
        if nb.size == 0:
            return nb, np.zeros(0)
        d = np.hypot(*(self.xy[nb] - self.xy[int(ic)]).T)
        with np.errstate(divide='ignore', invalid='ignore'):
            s = np.where(d > 0, (z[int(ic)] - z[nb]) / d, 0.0)
        return nb, s

    def neighbour_cells_ij(self, i, j=0):
        """Neighbours as ``(i, j)`` pairs under the ``(ncpl, 1)`` convention.

        This is the adapter that lets code written against a structured
        stencil -- SFR's priority flood, for one -- run unchanged on a mesh.
        """
        return [(int(k), 0) for k in self.neighbours[int(i)]]

    # ------------------------------------------------------------------ #
    def components(self, mask=None):
        """Connected components over the shared-face graph.

        ``mask`` restricts the walk to a subset of cells (the active ones, or
        the stream cells). Returns a list of lists of icell2d, largest first.
        """
        if mask is None:
            live = np.ones(self.ncpl, dtype=bool)
        else:
            live = np.asarray(mask, dtype=bool).reshape(-1)
        seen = np.zeros(self.ncpl, dtype=bool)
        out = []
        for start in np.flatnonzero(live):
            if seen[start]:
                continue
            stack, comp = [int(start)], []
            seen[start] = True
            while stack:
                c = stack.pop()
                comp.append(c)
                for k in self.neighbours[c]:
                    if live[k] and not seen[k]:
                        seen[k] = True
                        stack.append(int(k))
            out.append(sorted(comp))
        out.sort(key=len, reverse=True)
        return out

    def check(self, mask=None, raise_on_error=True):
        """Validate the topology and report what it found.

        The two properties that matter, and why:

          * **every interior face has exactly two cells** -- otherwise the mesh
            is malformed and any routing over it is meaningless (this one
            raises during construction, and is re-reported here);
          * **the graph is connected** over the cells of interest -- a second
            component cannot drain to the outlet, so SFR would report stream
            cells that reach nothing and CRR would strand water.
        """
        comps = self.components(mask)
        n_live = (self.ncpl if mask is None
                  else int(np.count_nonzero(np.asarray(mask).reshape(-1))))
        rep = {'ncpl': self.ncpl, 'faces': self.nfaces,
               'interior_faces': self.ninterior,
               'merged_vertices': self.duplicate_vertices,
               'isolated_cells': int(sum(1 for n in self.neighbours
                                         if len(n) == 0)),
               'boundary_cells': int(self.boundary.sum()),
               'cells_considered': n_live,
               'components': len(comps),
               'largest_component': (len(comps[0]) if comps else 0),
               'isolated': [c[0] for c in comps[1:][:10] if len(c) == 1],
               'mean_neighbours': (float(np.mean([len(n) for n in self.neighbours]))
                                   if self.ncpl else 0.0)}
        if raise_on_error and len(comps) > 1:
            extra = sum(len(c) for c in comps[1:])
            raise TopologyError(
                'the mesh graph has %d disconnected components: %d cell(s) '
                'outside the largest one (e.g. %s). Water there can never '
                'reach an outlet.'
                % (len(comps), extra, [c[:3] for c in comps[1:4]]))
        return rep
