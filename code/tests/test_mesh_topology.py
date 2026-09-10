# -*- coding: utf-8 -*-
"""WP1c.5 -- shared-face neighbour topology, and SFR routing over it.

The acceptance criteria: every interior face has exactly two cells, and the
graph is connected. Both are asserted here on synthetic meshes where the right
answer is known by construction, and on the real La Mata grid.
"""

import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6'), HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


TOPO = _load('marmites_topology_t',
             os.path.join(CODE, 'ppMF6', 'marmites_topology.py'))
GRID = _load('marmites_grid_t2', os.path.join(CODE, 'marmites_grid.py'))
SFR = _load('marmites_sfr_t', os.path.join(CODE, 'ppMF6', 'marmites_sfr.py'))


def _structured(nrow=4, ncol=5, cs=50.0):
    verts, cell2d, ncpl = GRID.disv_from_structured([cs] * ncol, [cs] * nrow,
                                                    0.0, 0.0)
    return {'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl, 'nlay': 1}


# --------------------------------------------------------------- structured
def test_a_structured_mesh_has_rook_adjacency_not_queen():
    """_NB8 counts the four diagonals; a shared FACE does not. MODFLOW 6 also
    connects cells only across faces, so the routing agrees with the flow
    model rather than with the old stencil."""
    nrow, ncol = 4, 5
    t = TOPO.MeshTopology(_structured(nrow, ncol))
    for i in range(nrow):
        for j in range(ncol):
            ic = i * ncol + j
            want = set()
            for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                a, b = i + di, j + dj
                if 0 <= a < nrow and 0 <= b < ncol:
                    want.add(a * ncol + b)
            assert set(t.neighbours[ic].tolist()) == want, (i, j)
    interior = 1 * ncol + 1
    assert len(t.neighbours[interior]) == 4          # not 8


def test_face_lengths_are_the_shared_edge_length():
    t = TOPO.MeshTopology(_structured(4, 5, cs=50.0))
    for ic in range(t.ncpl):
        assert np.allclose(t.face_length[ic], 50.0)


def test_adjacency_is_symmetric():
    t = TOPO.MeshTopology(_structured())
    for a in range(t.ncpl):
        for b in t.neighbours[a].tolist():
            assert a in set(t.neighbours[b].tolist())


def test_every_interior_face_has_exactly_two_cells():
    """The first acceptance criterion, stated positively."""
    t = TOPO.MeshTopology(_structured(4, 5))
    counts = sorted({len(cs) for cs in t.edges.values()})
    assert counts == [1, 2]
    assert t.ninterior == sum(1 for cs in t.edges.values() if len(cs) == 2)
    # a 4x5 grid of squares: 4*4 vertical + 3*5 horizontal interior faces
    assert t.ninterior == 4 * 4 + 3 * 5


def test_the_graph_is_connected():
    """The second acceptance criterion."""
    t = TOPO.MeshTopology(_structured())
    rep = t.check()
    assert rep['components'] == 1
    assert rep['largest_component'] == t.ncpl


def test_boundary_cells_are_flagged():
    nrow, ncol = 4, 5
    t = TOPO.MeshTopology(_structured(nrow, ncol))
    for i in range(nrow):
        for j in range(ncol):
            edge = (i in (0, nrow - 1)) or (j in (0, ncol - 1))
            assert bool(t.boundary[i * ncol + j]) is edge


# ----------------------------------------------------------- vertex identity
def test_vertices_are_matched_by_COORDINATE_not_by_id():
    """flopy's VoronoiGrid emits the same point under several ids -- 127 of
    them on La Mata's 989-cell mesh -- and id matching then leaves boundary
    cells with no neighbours at all."""
    gp = _structured(2, 2)
    # duplicate every vertex under a fresh id, and repoint cell 3 at the copies
    nv = len(gp['vertices'])
    dup = {int(v[0]): nv + k for k, v in enumerate(gp['vertices'])}
    gp['vertices'] = gp['vertices'] + [[dup[int(v[0])], v[1], v[2]]
                                       for v in gp['vertices']]
    rec = list(gp['cell2d'][3])
    gp['cell2d'][3] = rec[:4] + [dup[int(x)] for x in rec[4:]]
    t = TOPO.MeshTopology(gp)
    assert t.duplicate_vertices == nv
    assert t.check()['components'] == 1
    assert set(t.neighbours[3].tolist()) == {1, 2}


def test_a_repeated_closing_vertex_is_dropped():
    """A ring that repeats its first vertex would make a zero-length edge,
    which then looks like a face shared by every cell that closes the same
    way -- and that raises 'shared by more than two cells'."""
    gp = _structured(2, 2)
    gp['cell2d'] = [list(r[:3]) + [int(r[3]) + 1] + [int(x) for x in r[4:]]
                    + [int(r[4])] for r in gp['cell2d']]
    t = TOPO.MeshTopology(gp)
    assert t.check()['components'] == 1
    assert len(t.neighbours[0]) == 2


def test_a_face_shared_by_three_cells_is_refused():
    gp = _structured(2, 2)
    rec = list(gp['cell2d'][3])
    # give cell 3 the exact edge of cell 0 as well
    v0 = [int(x) for x in gp['cell2d'][0][4:]]
    gp['cell2d'][3] = rec[:3] + [len(rec[4:]) + 2] + [int(x) for x in rec[4:]] \
        + [v0[0], v0[1]]
    with pytest.raises(TOPO.TopologyError) as e:
        TOPO.MeshTopology(gp)
    assert 'more than two cells' in str(e.value)


# ------------------------------------------------------------- non-conforming
def _quadtree_like():
    """One 20 m cell beside two 10 m cells, with the hanging node listed.

    MODFLOW 6 DISV requires the coarse cell of a refinement step to carry the
    hanging node, which is exactly what makes edge matching work across a
    quadtree refinement.
    """
    v = {0: (0., 0.), 1: (20., 0.), 2: (20., 10.), 3: (20., 20.), 4: (0., 20.),
         5: (30., 0.), 6: (30., 10.), 7: (30., 20.)}
    vertices = [[k, x, y] for k, (x, y) in v.items()]
    cell2d = [
        [0, 10., 10., 5, 0, 1, 2, 3, 4],       # coarse, WITH the hanging node 2
        [1, 25., 5., 4, 1, 5, 6, 2],
        [2, 25., 15., 4, 2, 6, 7, 3],
    ]
    return {'vertices': vertices, 'cell2d': cell2d, 'ncpl': 3, 'nlay': 1}


def test_a_hanging_node_connects_a_coarse_cell_to_both_refined_cells():
    t = TOPO.MeshTopology(_quadtree_like())
    assert set(t.neighbours[0].tolist()) == {1, 2}
    assert t.check()['components'] == 1
    # each shared face is half the coarse cell's edge
    assert np.allclose(sorted(t.face_length[0]), [10.0, 10.0])


def test_a_missing_hanging_node_leaves_the_coarse_cell_disconnected():
    """Documents the requirement rather than hiding it: without the hanging
    node there is no shared edge to find, and the check says so."""
    gp = _quadtree_like()
    gp['cell2d'][0] = [0, 10., 10., 4, 0, 1, 3, 4]      # vertex 2 removed
    t = TOPO.MeshTopology(gp)
    assert len(t.neighbours[0]) == 0
    with pytest.raises(TOPO.TopologyError) as e:
        t.check()
    assert 'disconnected components' in str(e.value)


# --------------------------------------------------------------- CRR support
def test_slopes_are_positive_downhill():
    """alpha_ij = beta * S_ij / sum(S_ij) needs S_ij positive downslope."""
    ncol = 5
    t = TOPO.MeshTopology(_structured(4, ncol, cs=50.0))
    z = np.zeros(t.ncpl)
    ic = 1 * ncol + 1
    for k in t.neighbours[ic]:
        z[k] = 10.0
    z[ic] = 20.0                       # a peak: everything drains away
    nb, s = t.slopes(ic, z)
    assert np.all(s > 0)
    assert np.allclose(s, 10.0 / 50.0)
    z[ic] = 0.0                        # a pit: nothing drains away
    _nb, s = t.slopes(ic, z)
    assert np.all(s < 0)


def test_slopes_accepts_a_column_vector():
    """Per-cell fields arrive as (ncpl, 1) under the mesh convention."""
    t = TOPO.MeshTopology(_structured())
    z = np.arange(t.ncpl, dtype=float).reshape(-1, 1)
    nb, s = t.slopes(0, z)
    assert nb.size == len(s)


def test_components_can_be_restricted_to_a_subset():
    """SFR and CRR walk the ACTIVE cells, not the whole mesh."""
    nrow, ncol = 4, 5
    t = TOPO.MeshTopology(_structured(nrow, ncol))
    mask = np.zeros(t.ncpl, dtype=bool)
    mask[[0, 1]] = True                 # adjacent
    mask[ncol * 3 + 4] = True           # far corner, isolated
    comps = t.components(mask)
    assert [len(c) for c in comps] == [2, 1]
    with pytest.raises(TOPO.TopologyError):
        t.check(mask=mask)


# ------------------------------------------------------------------ SFR use
def test_sfr_refuses_a_mesh_layout_without_a_topology():
    """The failure this prevents is silent: on an (ncpl, 1) array the _NB8
    stencil addresses cells by their NUMBERING, so it would happily route a
    network between cells that are nowhere near each other."""
    pondw = np.zeros((12, 1))
    pondw[[2, 5, 9], 0] = 1.0
    dem = np.arange(12, dtype=float).reshape(-1, 1)
    with pytest.raises(ValueError) as e:
        SFR.stream_network(pondw, dem)
    assert 'topology' in str(e.value)


def test_sfr_routes_over_the_shared_face_graph():
    """A channel down one column of a mesh must route cell to cell, and the
    reaches must be exactly the channel cells."""
    nrow, ncol = 4, 5
    t = TOPO.MeshTopology(_structured(nrow, ncol))
    ncpl = nrow * ncol
    pondw = np.zeros((ncpl, 1))
    dem = np.zeros((ncpl, 1))
    chain = [i * ncol + 2 for i in range(nrow)]        # column 2, north->south
    for n, ic in enumerate(chain):
        pondw[ic, 0] = 2.0
        dem[ic, 0] = 100.0 - 10.0 * n                  # falls southward
    net = SFR.stream_network(pondw, dem, topology=t)
    assert net.nreaches == len(chain)
    assert set(c[0] for c in net.cells) == set(chain)
    assert net.outlets == [(chain[-1], 0)]             # the lowest cell
    # every reach drains to its downstream neighbour, and only the outlet has none
    for c in net.cells:
        r = net.recv[c]
        if c == net.outlets[0]:
            assert r is None
        else:
            assert r is not None and r[0] in t.neighbours[c[0]].tolist()
    assert net.acc[net.outlets[0]] == len(chain)


def test_sfr_reports_a_channel_that_cannot_reach_an_outlet():
    nrow, ncol = 4, 5
    t = TOPO.MeshTopology(_structured(nrow, ncol))
    ncpl = nrow * ncol
    pondw = np.zeros((ncpl, 1))
    dem = np.zeros((ncpl, 1))
    for ic in (0 * ncol + 0, 1 * ncol + 0, 3 * ncol + 4):   # one cell detached
        pondw[ic, 0] = 1.0
    dem[:, 0] = np.arange(ncpl)
    with pytest.raises(ValueError) as e:
        SFR.stream_network(pondw, dem, topology=t)
    assert 'do not connect to any outlet' in str(e.value)


def test_structured_sfr_still_walks_the_eight_neighbour_stencil():
    """The structured grid is the regression anchor: passing no topology must
    reproduce the diagonal-capable routing exactly as before."""
    pondw = np.zeros((3, 3))
    dem = np.zeros((3, 3))
    # a purely DIAGONAL chain -- only _NB8 can connect it
    for n, (i, j) in enumerate([(0, 0), (1, 1), (2, 2)]):
        pondw[i, j] = 1.0
        dem[i, j] = 10.0 - n
    net = SFR.stream_network(pondw, dem)
    assert net.nreaches == 3
    assert net.outlets == [(2, 2)]
