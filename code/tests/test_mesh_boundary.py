# -*- coding: utf-8 -*-
"""The voronoi mesh has no slivers on the catchment outline or at the edge
of a stream corridor.

Every vertex of a constraint polygon is a Voronoi GENERATOR. La Mata's
outline has 106 segments under 25 m (one of 3.5 m), so its boundary was a row
of slivers down to 1.86 m2 -- and the stream leaves the catchment through it:
the spin-up's cycle 2 (2026-09-27) failed at outlet reach cells of 225-617 m2.
With stream refinement on, the bands crossing the outline added a vertex at
every crossing: cells of ~0 m2. marmites_meshes.even_ring and band_margins.
"""
import copy
import math
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for p in (CODE, os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'MARMITESutilities')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_meshes as mm        # noqa: E402

DS = os.path.abspath(os.path.join(CODE, '..', 'example', 'LaMata'))


def _segments(ring):
    return [math.hypot(ring[(i + 1) % len(ring)][0] - ring[i][0],
                       ring[(i + 1) % len(ring)][1] - ring[i][1])
            for i in range(len(ring))]


def _area(ring):
    return abs(sum(x1 * y2 - x2 * y1 for (x1, y1), (x2, y2)
                   in zip(ring, ring[1:] + ring[:1]))) / 2.0


def test_the_outline_is_spaced_like_the_cells():
    """A square with a crowded side (3 m) and a bare one (1 km): every
    segment comes out between half the spacing and the spacing."""
    crowded = [(float(x), 0.0) for x in np.arange(0.0, 1000.0, 3.0)]
    ring = crowded + [(1000.0, 0.0), (1000.0, 1000.0), (0.0, 1000.0)]
    out, info = mm.even_ring(ring, 50.0)
    seg = _segments(out)
    assert not info['kept']
    assert min(seg) >= 25.0 - 1e-9 and max(seg) <= 50.0 + 1e-9, \
        (min(seg), max(seg))
    assert abs(_area(out) - 1e6) / 1e6 < 0.01
    assert abs(info['area_change']) < 0.01


def test_the_outline_never_leaves_a_convex_clip():
    """The ring is clipped to the model rectangle first; midpoints and
    points on segments of a convex polygon stay inside it."""
    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0.0, 2.0 * np.pi, 300))
    r = 400.0 + rng.uniform(-3.0, 3.0, t.size)
    ring = [(float(np.clip(500 + ri * np.cos(ti), 200.0, 800.0)),
             float(np.clip(500 + ri * np.sin(ti), 150.0, 850.0)))
            for ri, ti in zip(r, t)]
    out, _info = mm.even_ring(ring, 50.0)
    xs, ys = np.array(out).T
    assert xs.min() >= 200.0 - 1e-9 and xs.max() <= 800.0 + 1e-9
    assert ys.min() >= 150.0 - 1e-9 and ys.max() <= 850.0 + 1e-9


def test_the_band_margins_nest_and_stand_clear():
    """The outermost band stops half a background cell inside the outline,
    each finer band its own cell size further in."""
    m = mm.band_margins([(20.0, 20.0), (60.0, 40.0)], 50.0)
    assert m == [45.0, 25.0]
    m3 = mm.band_margins([(10.0, 10.0), (30.0, 20.0), (70.0, 40.0)], 50.0)
    assert m3 == [55.0, 45.0, 25.0], m3
    assert all(a > b for a, b in zip(m3, m3[1:])), 'finer bands further in'


# ------------------------------------------------ the La Mata mesh itself

def _triangle():
    try:
        import mm_paths
    except Exception:                                   # noqa: BLE001
        return None
    exe = getattr(mm_paths, 'TRIANGLE_EXE', None)
    return exe if exe and os.path.exists(exe) else None


def _lamata(tmp_path, near=None):
    if _triangle() is None:
        pytest.skip('triangle executable not available')
    if not os.path.exists(os.path.join(DS, 'inputWATERSHED.csv')):
        pytest.skip('La Mata dataset not present')
    pytest.importorskip('shapely')
    import marmites_config as mcfg
    cfg = mcfg.RunConfig.from_dict({})
    cfg.grid.kind = 'voronoi'
    v = cfg.grid.voronoi
    v.cell_far = 50.0
    v.refine_ponds = False
    v.stream_refine = near is not None
    if near is not None:
        v.cell_near_stream = float(near)
    v.refresh()
    cmf = SimpleNamespace(nrow=65, ncol=60, nlay=2, delr=np.full(60, 50.0),
                          delc=np.full(65, 50.0), xllcorner=739300.0,
                          yllcorner=4553050.0)
    gp, _info = mm.build_mesh(cfg, cmf, dataset_dir=DS,
                              model_ws=str(tmp_path), warn=lambda m: None,
                              force=True)
    vxy = {int(q[0]): (float(q[1]), float(q[2])) for q in gp['vertices']}
    rings = [[vxy[int(i)] for i in r[4:4 + int(r[3])]] for r in gp['cell2d']]
    areas = np.array([_area(r) for r in rings])
    ctr = np.array([(float(r[1]), float(r[2])) for r in gp['cell2d']])
    return areas, ctr


def _dist_to_outline(ctr):
    from shapely.geometry import Point, Polygon
    edge = Polygon(mm.watershed_ring(os.path.join(DS, 'inputWATERSHED.csv'))).exterior
    return np.array([edge.distance(Point(x, y)) for x, y in ctr])


def test_the_lamata_outline_carries_no_sliver(tmp_path):
    """Was: 51 cells under 250 m2, the smallest 1.9 m2, all on the outline."""
    areas, _ctr = _lamata(tmp_path)
    assert areas.min() >= 100.0, areas.min()
    assert (areas < 250.0).sum() <= 10, (areas < 250.0).sum()


def test_a_20_m_corridor_is_clean_up_to_the_outline(tmp_path):
    """Was: 230 cells under a quarter of the corridor size, 196 of them
    within 5 m of the outline, the smallest ~0 m2."""
    near = 20.0
    areas, ctr = _lamata(tmp_path, near=near)
    assert areas.min() >= 0.04 * near * near, areas.min()
    tiny = areas < 0.25 * near * near
    assert tiny.sum() <= 20, tiny.sum()
    assert not (_dist_to_outline(ctr[tiny]) < 5.0).any(), \
        'a band still crosses the outline'
