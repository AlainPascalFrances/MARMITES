# -*- coding: utf-8 -*-
"""Phase-4 tests: grid-agnostic ESRI ASCII sampling.

Key acceptance: sampling the real La Mata rasters at DIS cell centres must
reproduce the legacy cell-for-cell reader exactly. That is what allows the
same inputs to feed a DISV (or refined) grid, where the legacy reader cannot
be used at all.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))
sys.path.insert(0, TRUNK)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R = _load('marmites_raster', os.path.join(TRUNK, 'marmites_raster.py'))


def _write_asc(tmp, name, arr, cellsize=10.0, xll=0.0, yll=0.0, nodata=-9999):
    p = os.path.join(str(tmp), name)
    a = np.asarray(arr)
    with open(p, 'w') as f:
        f.write('ncols %d\nnrows %d\nxllcorner %g\nyllcorner %g\ncellsize %g\n'
                'NODATA_value %g\n' % (a.shape[1], a.shape[0], xll, yll, cellsize, nodata))
        for row in a:
            f.write(' '.join('%g' % v for v in row) + '\n')
    return p


# --------------------------------------------------------------------- #

def test_read_header_and_bounds(tmp_path):
    arr = np.arange(6).reshape(2, 3)
    p = _write_asc(tmp_path, 'a.asc', arr, cellsize=10.0, xll=100.0, yll=200.0)
    r = R.read_esri_ascii(p)
    assert (r.nrows, r.ncols) == (2, 3)
    assert r.cellsize == 10.0
    assert r.bounds == (100.0, 200.0, 130.0, 220.0)
    assert np.array_equal(r.array, arr)


def test_sample_nearest_cell_centres(tmp_path):
    # 2x3 grid, 10 m cells, origin (0,0); row 0 is NORTH (y 10..20)
    arr = np.array([[1, 2, 3],
                    [4, 5, 6]])
    p = _write_asc(tmp_path, 'z.asc', arr, cellsize=10.0)
    r = R.read_esri_ascii(p, dtype=int)
    # centres of the north row: y=15 ; south row: y=5
    v = R.sample_points(r, [5.0, 15.0, 25.0, 5.0, 25.0], [15.0, 15.0, 15.0, 5.0, 5.0])
    assert np.array_equal(v, [1, 2, 3, 4, 6])


def test_sample_outside_returns_fill(tmp_path):
    p = _write_asc(tmp_path, 'z.asc', np.ones((2, 2)), cellsize=10.0)
    r = R.read_esri_ascii(p)
    v = R.sample_points(r, [-5.0, 5.0], [5.0, 5.0], outside=-999.0)
    assert v[0] == -999.0 and v[1] == 1.0


def test_nodata_mapped_to_fill(tmp_path):
    arr = np.array([[1.0, -9999.0], [3.0, 4.0]])
    p = _write_asc(tmp_path, 'n.asc', arr, cellsize=10.0)
    xy = R.cell_centres_structured([10.0, 10.0], [10.0, 10.0])
    v = R.sample_to_cells(p, xy, hnoflo=999.0)
    assert v[1] == 999.0                      # NODATA cell
    assert np.allclose([v[0], v[2], v[3]], [1.0, 3.0, 4.0])


def test_zone_raster_stays_integer(tmp_path):
    arr = np.array([[1, 2], [3, 1]])
    p = _write_asc(tmp_path, 'zone.asc', arr, cellsize=10.0)
    xy = R.cell_centres_structured([10.0, 10.0], [10.0, 10.0])
    v = R.sample_to_cells(p, xy, dtype=int)
    assert v.dtype.kind in 'iu'
    assert np.array_equal(v, [1, 2, 3, 1])


def test_cell_centres_structured_orientation():
    xy = R.cell_centres_structured([10.0, 10.0], [10.0, 10.0], xorigin=0.0, yorigin=0.0)
    # cell (0,0) is north-west -> y greater than cell (1,0)
    assert xy[0][1] > xy[2][1]
    assert np.allclose(xy[0], [5.0, 15.0])
    assert np.allclose(xy[3], [15.0, 5.0])


# --------------------------------------------------------------------- #
# real La Mata rasters: sampling == legacy reader
# --------------------------------------------------------------------- #

@pytest.mark.parametrize('fn,dt', [
    ('inputSOILzones.asc', int),
    ('inputMETEOzones.asc', int),
    ('inputSOILthick.asc', float),
    ('inputSTREAMw.asc', float),
])
def test_sampling_matches_legacy_reader_on_lamata(fn, dt):
    path = os.path.join(DS, fn)
    if not os.path.exists(path):
        pytest.skip('La Mata dataset not present')
    ras = R.read_esri_ascii(path, dtype=float)
    # cell centres of the equivalent structured model grid
    delr = [ras.cellsize] * ras.ncols
    delc = [ras.cellsize] * ras.nrows
    xy = R.cell_centres_structured(delr, delc, ras.xllcorner, ras.yllcorner)
    sampled = R.sample_to_cells(path, xy, dtype=dt, hnoflo=np.nan)
    legacy = ras.array.ravel()
    m = ~np.isclose(legacy, float(ras.nodata))
    if dt is int:
        assert np.array_equal(sampled[m], np.rint(legacy[m]).astype(int))
    else:
        assert np.allclose(sampled[m], legacy[m])
    assert m.sum() > 0
