# -*- coding: utf-8 -*-
"""WP1d -- wrapping vector layers onto the model grid.

The whole point of the new paradigm is that a cell 37 % covered by alluvium
gets 37, not the class that happened to sit under its centre. So the tests
are about EXACTNESS: clipped areas against hand-computed ones, holes actually
subtracted, and the same answer whichever way round the source file writes
its rings.
"""

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
if CODE not in sys.path:
    sys.path.insert(0, CODE)

import marmites_vector as mv          # noqa: E402

shapefile = pytest.importorskip('shapefile', reason='pyshp not installed')


# =====================================================================
#  geometry primitives
# =====================================================================

SQUARE_CCW = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]


def test_signed_area_sign():
    assert mv._signed_area(SQUARE_CCW) == pytest.approx(1.0)
    assert mv._signed_area(SQUARE_CCW[::-1]) == pytest.approx(-1.0)
    assert mv._signed_area([(0, 0), (1, 1)]) == 0.0


def test_to_ccw_normalises_and_drops_closing_point():
    closed_cw = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
    p = mv._to_ccw(closed_cw)
    assert len(p) == 4
    assert mv._signed_area(p) > 0


def test_clip_square_by_offset_square():
    clip = mv._to_ccw([(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)])
    out = mv._clip_to_convex(SQUARE_CCW, clip)
    assert abs(mv._signed_area(out)) == pytest.approx(0.25)


def test_clip_subject_entirely_inside_is_unchanged_in_area():
    clip = mv._to_ccw([(-1, -1), (2, -1), (2, 2), (-1, 2)])
    out = mv._clip_to_convex(SQUARE_CCW, clip)
    assert abs(mv._signed_area(out)) == pytest.approx(1.0)


def test_clip_disjoint_gives_nothing():
    clip = mv._to_ccw([(5, 5), (6, 5), (6, 6), (5, 6)])
    assert len(mv._clip_to_convex(SQUARE_CCW, clip)) < 3


def test_clip_nonconvex_subject_is_exact():
    """An L: 2x2 with the top-right 1x1 removed, area 3. Clip keeps 0.75."""
    ell = [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]
    clip = mv._to_ccw([(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)])
    out = mv._clip_to_convex(ell, clip)
    assert abs(mv._signed_area(out)) == pytest.approx(0.75)


def test_clip_against_a_triangle_is_exact():
    """Voronoi cells are not rectangles, so the clip window must be general."""
    tri = mv._to_ccw([(0.0, 0.0), (2.0, 0.0), (0.0, 2.0)])
    out = mv._clip_to_convex(SQUARE_CCW, tri)
    assert abs(mv._signed_area(out)) == pytest.approx(1.0)   # unit square is inside


def test_segment_clip_lengths():
    clip = mv._to_ccw(SQUARE_CCW)
    ln, _ = mv._clip_segment_to_convex((-1.0, 0.5), (2.0, 0.5), clip)
    assert ln == pytest.approx(1.0)                      # straight through
    ln, _ = mv._clip_segment_to_convex((0.25, 0.5), (0.75, 0.5), clip)
    assert ln == pytest.approx(0.5)                      # wholly inside
    ln, _ = mv._clip_segment_to_convex((2.0, 2.0), (3.0, 3.0), clip)
    assert ln == 0.0                                     # wholly outside
    ln, _ = mv._clip_segment_to_convex((-1.0, 0.5), (0.5, 0.5), clip)
    assert ln == pytest.approx(0.5)                      # one end inside


def test_point_in_convex():
    clip = mv._to_ccw(SQUARE_CCW)
    assert mv._point_in_convex(0.5, 0.5, clip)
    assert not mv._point_in_convex(1.5, 0.5, clip)


# =====================================================================
#  the target grid
# =====================================================================

def test_structured_grid_geometry():
    g = mv.TargetGrid.structured([10.0] * 3, [10.0] * 2, 0.0, 0.0)
    assert g.shape == (2, 3) and g.ncell == 6
    assert np.allclose(g.area, 100.0)
    # row 0 is the NORTHERNMOST row, as MODFLOW has it
    north = g.polygons[0]
    south = g.polygons[3]
    assert max(y for _, y in north) == pytest.approx(20.0)
    assert max(y for _, y in south) == pytest.approx(10.0)


def test_structured_cell_containing_and_unravel():
    g = mv.TargetGrid.structured([10.0] * 3, [10.0] * 2, 0.0, 0.0)
    ic = g.cell_containing(25.0, 15.0)            # top row, third column
    assert g.unravel(ic) == (0, 2)
    assert g.cell_containing(-5.0, 5.0) == -1


def test_gridprops_grid_uses_the_ncpl_1_convention():
    verts = [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 1.0, 1.0], [3, 0.0, 1.0],
             [4, 2.0, 0.0], [5, 2.0, 1.0]]
    cell2d = [[0, 0.5, 0.5, 4, 0, 1, 2, 3],
              [1, 1.5, 0.5, 4, 1, 4, 5, 2]]
    g = mv.TargetGrid.from_gridprops({'vertices': verts, 'cell2d': cell2d,
                                      'ncpl': 2})
    assert g.shape == (2, 1)
    assert np.allclose(g.area, 1.0)
    assert g.unravel(g.cell_containing(1.5, 0.5)) == (1, 0)


# =====================================================================
#  polygon overlay
# =====================================================================

def _write_poly(path, rings_per_feature, fields, records):
    w = shapefile.Writer(path, shapeType=shapefile.POLYGON)
    for name, typ, size, dec in fields:
        w.field(name, typ, size, dec)
    for rings, rec in zip(rings_per_feature, records):
        w.poly(rings)
        w.record(*rec)
    w.close()
    with open(os.path.splitext(path)[0] + '.prj', 'w') as fh:
        fh.write('PROJCS["test"]')


def _cw(ring):
    """Close a ring and make it clockwise -- the ESRI outer-ring convention."""
    r = list(ring)
    if mv._signed_area(r) > 0:
        r = r[::-1]
    return r + [r[0]]


@pytest.fixture
def grid4():
    """2 x 2 grid of 10 m cells with its origin at (0, 0)."""
    return mv.TargetGrid.structured([10.0, 10.0], [10.0, 10.0], 0.0, 0.0)


def test_overlay_majority_picks_the_dominant_class(tmp_path, grid4):
    """Cell (0,0) spans x 0..10, y 10..20. Zone 1 takes 6 m of it, zone 2 4 m."""
    p = str(tmp_path / 'zones.shp')
    _write_poly(p,
                [[_cw([(0, 10), (6, 10), (6, 20), (0, 20)])],
                 [_cw([(6, 10), (10, 10), (10, 20), (6, 20)])]],
                [('code', 'N', 4, 0)], [[1], [2]])
    arr, rep = mv.overlay_polygons(mv.Layer(p), grid4, field='code',
                                   how='majority', fill=0, dtype=int)
    assert arr.shape == (2, 2)
    assert arr[0, 0] == 1
    assert arr[1, 0] == 0            # nothing reaches the southern row
    assert rep['cells_touched'] == 1


def test_overlay_area_fraction_is_exact(tmp_path, grid4):
    p = str(tmp_path / 'veg.shp')
    # 3 m x 10 m strip inside cell (0,0): 30 of 100 m2 -> 30 %
    _write_poly(p, [[_cw([(0, 10), (3, 10), (3, 20), (0, 20)])]],
                [('sp', 'C', 8, 0)], [['i']])
    arr, rep = mv.overlay_polygons(mv.Layer(p), grid4, how='area_fraction',
                                   fill=0.0, dtype=float)
    assert arr[0, 0] == pytest.approx(30.0)
    assert arr[0, 1] == pytest.approx(0.0)
    assert rep['over_100_cells'] == 0


def test_overlay_area_fraction_never_exceeds_100_for_one_feature(tmp_path, grid4):
    p = str(tmp_path / 'all.shp')
    _write_poly(p, [[_cw([(-5, -5), (25, -5), (25, 25), (-5, 25)])]],
                [('sp', 'C', 8, 0)], [['i']])
    arr, _ = mv.overlay_polygons(mv.Layer(p), grid4, how='area_fraction')
    assert np.allclose(arr, 100.0)


def test_overlay_subtracts_holes(tmp_path, grid4):
    """A 10x10 cell fully covered but with a 4x4 hole -> 84 %."""
    p = str(tmp_path / 'hole.shp')
    outer = _cw([(0, 10), (10, 10), (10, 20), (0, 20)])
    hole = _cw([(3, 13), (7, 13), (7, 17), (3, 17)])[::-1]   # opposite winding
    _write_poly(p, [[outer, hole]], [('sp', 'C', 8, 0)], [['i']])
    arr, _ = mv.overlay_polygons(mv.Layer(p), grid4, how='area_fraction')
    assert arr[0, 0] == pytest.approx(84.0)


def test_overlay_is_indifferent_to_ring_winding(tmp_path, grid4):
    """A file that writes outer rings counter-clockwise must still work."""
    a = str(tmp_path / 'cw.shp')
    b = str(tmp_path / 'ccw.shp')
    ring = [(0, 10), (6, 10), (6, 20), (0, 20)]
    _write_poly(a, [[_cw(ring)]], [('c', 'N', 4, 0)], [[1]])
    _write_poly(b, [[_cw(ring)[::-1]]], [('c', 'N', 4, 0)], [[1]])
    fa, _ = mv.overlay_polygons(mv.Layer(a), grid4, how='area_fraction')
    fb, _ = mv.overlay_polygons(mv.Layer(b), grid4, how='area_fraction')
    assert fa[0, 0] == pytest.approx(60.0)
    assert np.allclose(fa, fb)


def test_overlay_area_mean_weights_by_area(tmp_path, grid4):
    """6 m of thickness 1.0 and 4 m of 2.0 -> 1.4 over that cell."""
    p = str(tmp_path / 'thick.shp')
    _write_poly(p,
                [[_cw([(0, 10), (6, 10), (6, 20), (0, 20)])],
                 [_cw([(6, 10), (10, 10), (10, 20), (6, 20)])]],
                [('thick_m', 'N', 10, 3)], [[1.0], [2.0]])
    arr, _ = mv.overlay_polygons(mv.Layer(p), grid4, field='thick_m',
                                 how='area_mean', fill=np.nan)
    assert arr[0, 0] == pytest.approx(1.4)
    assert np.isnan(arr[1, 0])


def test_overlay_select_callable_filters_features(tmp_path, grid4):
    """'blank Species means grass' is expressed by the caller, not the module."""
    p = str(tmp_path / 'veg.shp')
    _write_poly(p,
                [[_cw([(0, 10), (5, 10), (5, 20), (0, 20)])],
                 [_cw([(5, 10), (10, 10), (10, 20), (5, 20)])]],
                [('Species', 'C', 8, 0)], [['i'], ['']])
    lay = mv.Layer(p)
    ilex, _ = mv.overlay_polygons(lay, grid4, field='Species',
                                  how='area_fraction',
                                  select=lambda v: str(v).strip() == 'i')
    grass, _ = mv.overlay_polygons(lay, grid4, field='Species',
                                   how='area_fraction',
                                   select=lambda v: str(v).strip() == '')
    assert ilex[0, 0] == pytest.approx(50.0)
    assert grass[0, 0] == pytest.approx(50.0)


def test_overlay_presence(tmp_path, grid4):
    p = str(tmp_path / 'p.shp')
    _write_poly(p, [[_cw([(0, 10), (1, 10), (1, 11), (0, 11)])]],
                [('c', 'N', 4, 0)], [[1]])
    arr, _ = mv.overlay_polygons(mv.Layer(p), grid4, how='presence', dtype=int)
    assert arr[0, 0] == 1 and arr.sum() == 1


def test_overlay_rejects_a_line_layer(tmp_path, grid4):
    p = str(tmp_path / 'l.shp')
    w = shapefile.Writer(p, shapeType=shapefile.POLYLINE)
    w.field('c', 'N', 4, 0)
    w.line([[(0, 0), (1, 1)]])
    w.record(1)
    w.close()
    with pytest.raises(mv.VectorError):
        mv.overlay_polygons(mv.Layer(p), grid4, how='presence')


def test_overlay_rejects_an_unknown_mode(tmp_path, grid4):
    p = str(tmp_path / 'z.shp')
    _write_poly(p, [[_cw([(0, 10), (1, 10), (1, 11), (0, 11)])]],
                [('c', 'N', 4, 0)], [[1]])
    with pytest.raises(mv.VectorError):
        mv.overlay_polygons(mv.Layer(p), grid4, how='nearest')


# =====================================================================
#  lines
# =====================================================================

def _write_line(path, parts_per_feature, fields, records):
    w = shapefile.Writer(path, shapeType=shapefile.POLYLINE)
    for name, typ, size, dec in fields:
        w.field(name, typ, size, dec)
    for parts, rec in zip(parts_per_feature, records):
        w.line(parts)
        w.record(*rec)
    w.close()


def test_burn_lines_length_and_longest(tmp_path, grid4):
    """A west-east line at y = 15 crosses the whole northern row: 10 m a cell."""
    p = str(tmp_path / 's.shp')
    _write_line(p, [[[(0, 15), (12, 15)]], [[(12, 15), (20, 15)]]],
                [('w_m', 'N', 10, 3)], [[1.5], [3.0]])
    lay = mv.Layer(p)
    ln, rep = mv.burn_lines(lay, grid4, how='length')
    assert ln[0, 0] == pytest.approx(10.0)
    assert ln[0, 1] == pytest.approx(10.0)
    assert ln[1, 0] == pytest.approx(0.0)
    assert rep['total_length_m'] == pytest.approx(20.0)

    wid, _ = mv.burn_lines(lay, grid4, field='w_m', how='longest', fill=0.0)
    assert wid[0, 0] == pytest.approx(1.5)      # only feature 0 reaches it
    assert wid[0, 1] == pytest.approx(3.0)      # 8 m of feature 1 beats 2 m of 0


def test_burn_lines_presence_and_rejection(tmp_path, grid4):
    p = str(tmp_path / 's2.shp')
    _write_line(p, [[[(1, 11), (2, 12)]]], [('c', 'N', 4, 0)], [[1]])
    arr, _ = mv.burn_lines(mv.Layer(p), grid4, how='presence', dtype=int)
    assert arr[0, 0] == 1 and arr.sum() == 1
    with pytest.raises(mv.VectorError):
        mv.burn_lines(mv.Layer(p), grid4, how='nearest')


# =====================================================================
#  points
# =====================================================================

def test_locate_points_reports_the_ones_outside(tmp_path, grid4):
    p = str(tmp_path / 'obs.shp')
    w = shapefile.Writer(p, shapeType=shapefile.POINT)
    w.field('Name', 'C', 8, 0)
    for nm, x, y in (('A', 5.0, 15.0), ('B', 15.0, 5.0), ('C', -5.0, 5.0)):
        w.point(x, y)
        w.record(nm)
    w.close()
    got = mv.locate_points(mv.Layer(p), grid4, name_field='Name')
    by = {g['name']: g for g in got}
    assert (by['A']['row'], by['A']['col']) == (0, 0)
    assert (by['B']['row'], by['B']['col']) == (1, 1)
    assert by['C']['inside'] is False and by['C']['icell'] == -1


# =====================================================================
#  layer plumbing
# =====================================================================

def test_layer_missing_column_names_what_it_has(tmp_path):
    p = str(tmp_path / 'z.shp')
    _write_poly(p, [[_cw([(0, 0), (1, 0), (1, 1), (0, 1)])]],
                [('SoilCode', 'N', 4, 0)], [[1]])
    lay = mv.Layer(p)
    assert lay.kind == 'polygon'
    assert lay.crs_wkt.startswith('PROJCS')
    with pytest.raises(mv.VectorError) as e:
        lay.require('thick_m')
    assert 'SoilCode' in str(e.value)


def test_layer_missing_file():
    with pytest.raises(mv.VectorError):
        mv.Layer(os.path.join('nowhere', 'nothing.shp'))


def test_coverage_report_mentions_over_100(tmp_path, grid4):
    p = str(tmp_path / 'a.shp')
    ring = _cw([(0, 10), (10, 10), (10, 20), (0, 20)])
    _write_poly(p, [[ring], [ring]], [('c', 'N', 4, 0)], [[1], [2]])
    _, rep = mv.overlay_polygons(mv.Layer(p), grid4, how='area_fraction')
    txt = mv.coverage_report([rep])
    assert 'OVER 100' in txt
