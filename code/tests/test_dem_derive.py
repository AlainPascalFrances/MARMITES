# -*- coding: utf-8 -*-
"""WP1d -- the elevation rasters, derived from the GIS DEM.

``MF_ws/elev_sinkfil.asc`` is the model's land surface and ``[crr] dem`` is
the cascade's copy of it. Both used to be placed by hand, which a new
catchment cannot do. They are now block-averaged from the raster ``[grid]
dem`` names, onto the rectangle the OTHER dataset rasters are on -- because
a DEM on a different rectangle from ibound would misalign against it rather
than fail.
"""

import importlib.util
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))


def _tool():
    spec = importlib.util.spec_from_file_location(
        'gis_to_dataset_probe', os.path.join(CODE, 'tools',
                                             'gis_to_dataset.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


tool = _tool()


def _asc(path, ncols=4, nrows=3, xll=0.0, yll=0.0, cs=10.0, value=1.0):
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('ncols         %d\nnrows         %d\n' % (ncols, nrows))
        fh.write('xllcorner     %g\nyllcorner     %g\n' % (xll, yll))
        fh.write('cellsize      %g\nNODATA_value  -9999\n' % cs)
        for _ in range(nrows):
            fh.write(' '.join(['%g' % value] * ncols) + '\n')
    return path


def test_the_header_is_read_back(tmp_path):
    p = _asc(str(tmp_path / 'a.asc'), ncols=60, nrows=65, xll=739300.0,
             yll=4553050.0, cs=50.0)
    assert tool._read_asc_header(p) == (60, 65, 739300.0, 4553050.0, 50.0)
    assert tool._read_asc_header(str(tmp_path / 'nope.asc')) is None


def test_a_northing_is_not_written_in_scientific_notation(tmp_path):
    """%g turns 4553050 into 4.55305e+06, which the model's reader does not
    parse -- and the file would look right until a run failed on it."""
    out = str(tmp_path / 'dem.asc')
    tool._write_asc(out, np.zeros((2, 3)), 739300.0, 4553050.0, 50.0)
    head = open(out, encoding='utf-8').read().splitlines()[:6]
    assert 'yllcorner     4553050' in head
    assert not any('e+' in ln for ln in head)
    assert tool._read_asc_header(out) == (3, 2, 739300.0, 4553050.0, 50.0)


def test_the_siblings_decide_the_rectangle(tmp_path):
    """Every MF raster shares one header. A DEM on another rectangle would
    misalign against all of them rather than fail outright."""
    mf = tmp_path / 'MF_ws'
    mf.mkdir()
    assert tool._sibling_rectangle(str(tmp_path)) == (None, '')
    _asc(str(mf / 'ibound_l1.asc'), ncols=60, nrows=65, xll=739300.0,
         yll=4553050.0, cs=50.0)
    hdr, named = tool._sibling_rectangle(str(tmp_path))
    assert hdr == (60, 65, 739300.0, 4553050.0, 50.0)
    assert named == 'ibound_l1.asc'


def test_the_dem_is_averaged_not_sampled(tmp_path):
    """The DEM is 5 m and a model cell is 50 m, so one cell covers a hundred
    pixels; taking the one under its centre throws away ninety-nine."""
    rio = pytest.importorskip('rasterio')
    from rasterio.transform import from_origin
    src = str(tmp_path / 'dem.asc')
    # 4 x 4 pixels of 5 m: a 20 m square whose left half is 100 and right 200
    a = np.zeros((4, 4), dtype='float32')
    a[:, :2], a[:, 2:] = 100.0, 200.0
    with rio.open(src, 'w', driver='AAIGrid', height=4, width=4, count=1,
                  dtype='float32', crs='EPSG:23029',
                  transform=from_origin(0.0, 20.0, 5.0, 5.0)) as dst:
        dst.write(a, 1)

    # one 20 m cell over all of it: the MEAN of the two halves
    one, filled, total = tool._resample_dem(src, 1, 1, 0.0, 0.0, 20.0)
    assert (filled, total) == (1, 1)
    assert abs(float(one[0, 0]) - 150.0) < 1e-6

    # two 10 m cells across: each sees one half
    two, filled, _t = tool._resample_dem(src, 2, 1, 0.0, 0.0, 10.0)
    assert filled == 2
    assert abs(float(two[0, 0]) - 100.0) < 1e-6
    assert abs(float(two[0, 1]) - 200.0) < 1e-6


def test_cells_the_dem_does_not_reach_are_nodata(tmp_path):
    """A hole in the model's surface, counted and reported rather than
    filled with something plausible."""
    rio = pytest.importorskip('rasterio')
    from rasterio.transform import from_origin
    src = str(tmp_path / 'dem.asc')
    with rio.open(src, 'w', driver='AAIGrid', height=2, width=2, count=1,
                  dtype='float32', crs='EPSG:23029',
                  transform=from_origin(0.0, 20.0, 10.0, 10.0)) as dst:
        dst.write(np.full((2, 2), 700.0, dtype='float32'), 1)
    arr, filled, total = tool._resample_dem(src, 4, 2, 0.0, 0.0, 10.0)
    assert (filled, total) == (4, 8)
    assert (arr == -9999.0).sum() == 4
    assert (arr[arr != -9999.0] == 700.0).all()
