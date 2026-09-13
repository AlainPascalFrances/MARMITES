# -*- coding: utf-8 -*-
"""WP1d -- the land surface, at its own resolution, onto the model grid.

The elevation used to reach the model as a raster already resampled to the
50 m rectangle, which was then projected AGAIN onto the run's mesh. Two
resamplings, and the first one discarded what the 5 m survey knew -- 0.65 m
rms and 4.4 m at worst on La Mata. So the dataset keeps the DEM fine and
grid-independent and ``marmites_dem`` wraps it onto whichever cells panel 1
produced.

``marmites_dem`` is MODEL-PATH code: numpy only, no rasterio. The converter
does the reading of whatever GDAL format the DEM arrives in.
"""

import importlib.util
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mdem = _load('marmites_dem_probe', os.path.join(CODE, 'marmites_dem.py'))


def _write(tmp_path, arr, xll=0.0, yll=0.0, cs=5.0, name='dem.asc'):
    return mdem.write_asc(str(tmp_path / name), np.asarray(arr, dtype=float),
                          xll, yll, cs)


def test_an_asc_round_trips(tmp_path):
    a = np.arange(12, dtype=float).reshape(3, 4)
    p = _write(tmp_path, a, xll=739300.0, yll=4553050.0, cs=5.0)
    back, head = mdem.read_asc(p)
    assert np.allclose(np.ma.filled(back, -1), a)
    assert (head['ncols'], head['nrows']) == (4, 3)
    assert head['xllcorner'] == 739300.0 and head['yllcorner'] == 4553050.0
    assert head['cellsize'] == 5.0


def test_a_northing_is_not_written_in_scientific_notation(tmp_path):
    """%g turns 4553050 into 4.55305e+06, which readers of this format do not
    parse -- and the file looks right until something fails on it."""
    p = _write(tmp_path, np.zeros((2, 2)), xll=739300.0, yll=4553050.0)
    head = open(p, encoding='utf-8').read().splitlines()[:6]
    assert 'yllcorner     4553050' in head
    assert not any('e+' in ln for ln in head)


def test_nodata_stays_a_hole(tmp_path):
    """A gap in the survey must not average into its neighbours as -9999."""
    a = np.array([[10.0, -9999.0], [30.0, 40.0]])
    back, _h = mdem.read_asc(_write(tmp_path, a))
    assert np.ma.getmaskarray(back)[0, 1]
    assert np.ma.compressed(back).tolist() == [10.0, 30.0, 40.0]


def test_a_truncated_file_is_an_error(tmp_path):
    p = str(tmp_path / 'bad.asc')
    with open(p, 'w', encoding='utf-8') as fh:
        fh.write('ncols 2\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 5\n'
                 'NODATA_value -9999\n1 2\n')
    with pytest.raises(mdem.DEMError):
        mdem.read_asc(p)
    with pytest.raises(mdem.DEMError):
        mdem.read_asc(str(tmp_path / 'absent.asc'))


def _square_grid(n, size, x0=0.0, y0=0.0):
    """A DISV gridprops of n x n squares -- every kind arrives as these."""
    import sys
    sys.path.insert(0, CODE)
    from marmites_grid import disv_from_structured
    verts, cell2d, ncpl = disv_from_structured(
        np.full(n, size), np.full(n, size), x0, y0)
    return {'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl, 'nlay': 1}


def test_the_wrap_is_an_average_not_a_sample(tmp_path):
    """One model cell covers a hundred DEM pixels; taking the one under its
    centre throws away ninety-nine."""
    # 4 x 4 pixels of 5 m over a 20 m square: west half 100, east half 200
    a = np.zeros((4, 4))
    a[:, :2], a[:, 2:] = 100.0, 200.0
    p = _write(tmp_path, a, xll=0.0, yll=0.0, cs=5.0)

    one = _square_grid(1, 20.0)
    arr, info = mdem.wrap_to_grid(p, one)
    assert arr.shape == (1, 1)
    assert abs(float(arr[0, 0]) - 150.0) < 1e-6      # the mean of both halves
    assert info['cellsize'] == 5.0

    two = _square_grid(2, 10.0)
    arr, _i = mdem.wrap_to_grid(p, two)
    got = np.ma.filled(arr, -1).ravel()
    assert sorted(got.tolist()) == [100.0, 100.0, 200.0, 200.0]


def test_a_cell_the_dem_does_not_reach_is_masked_and_reported(tmp_path):
    """A hole in the surface is said out loud, not filled with a plausible
    number."""
    p = _write(tmp_path, np.full((2, 2), 700.0), xll=0.0, yll=0.0, cs=10.0)
    grid = _square_grid(2, 20.0, x0=0.0, y0=0.0)    # twice the DEM's reach
    said = []
    arr, info = mdem.wrap_to_grid(p, grid, warn=said.append)
    assert info['missing'] == 3 and np.ma.getmaskarray(arr).sum() == 3
    assert said and 'no DEM' in said[0]
    assert abs(float(np.ma.compressed(arr)[0]) - 700.0) < 1e-6


def test_the_result_is_cached_on_the_dem_and_the_mesh(tmp_path):
    """Seconds, not milliseconds: 18 s for the 5 m La Mata DEM on a
    5778-cell mesh, so a run must not repeat it."""
    p = _write(tmp_path, np.full((4, 4), 700.0), cs=5.0)
    grid = _square_grid(2, 10.0)
    cache = str(tmp_path / 'cache')
    first, i1 = mdem.wrap_to_grid(p, grid, cache_dir=cache)
    assert i1['cached'] is False
    second, i2 = mdem.wrap_to_grid(p, grid, cache_dir=cache)
    assert i2['cached'] is True
    assert np.allclose(np.ma.filled(first, -1), np.ma.filled(second, -1))

    # A DIFFERENT mesh must not be served the first one's answer.
    other = _square_grid(4, 5.0)
    _third, i3 = mdem.wrap_to_grid(p, other, cache_dir=cache)
    assert i3['cached'] is False and i3['signature'] != i1['signature']


def test_the_model_path_does_not_need_rasterio():
    """Cookbook D3: geopandas and rasterio live in the converter and the app,
    never where the model runs. The converter reads the GIS raster in
    whatever format GDAL opens; this module reads the dataset copy."""
    import ast
    src = open(os.path.join(CODE, 'marmites_dem.py'), encoding='utf-8').read()
    imported = set()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Import):
            imported.update(a.name.split('.')[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split('.')[0])
    assert not ({'rasterio', 'geopandas', 'fiona', 'shapely', 'streamlit'}
                & imported), sorted(imported)
