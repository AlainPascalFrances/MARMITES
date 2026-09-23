# -*- coding: utf-8 -*-
"""Polygon layers onto the model grid: the engine, then La Mata.

Until marmites_overlay nothing could turn a polygon into a grid value, so
the Soil panel's zones and thickness -- asked as a shapefile and a column --
were never read: the run loaded inputSOILzones.asc and inputSOILthick.asc
from filenames written into the driver.
"""

import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(CODE, '..', 'example', 'LaMata'))
for _p in (CODE, os.path.join(CODE, 'ppMF6'),
           os.path.join(CODE, 'MARMITESutilities'),
           os.path.join(CODE, 'ppMF_FloPy')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

pytest.importorskip('shapely')
import marmites_overlay as ov                                  # noqa: E402


def _layer(tmp_path, polys, name='p.geojson'):
    """polys: [(x0, y0, x1, y1, {props})] rectangles."""
    feats = [{'type': 'Feature', 'properties': pr,
              'geometry': {'type': 'Polygon',
                           'coordinates': [[[x0, y0], [x1, y0], [x1, y1],
                                            [x0, y1], [x0, y0]]]}}
             for x0, y0, x1, y1, pr in polys]
    p = tmp_path / name
    p.write_text(json.dumps({'type': 'FeatureCollection', 'features': feats}),
                 encoding='utf-8')
    return str(p)


# ---------------------------------------------------------------- the engine
# One row of two 10 m cells: x 0..10 and 10..20, y 0..10.

def _cells():
    return ov.structured_cells(0.0, 0.0, [10.0, 10.0], [10.0])


def test_the_cells_are_row_major_with_row_zero_at_the_north():
    c = ov.structured_cells(0.0, 0.0, [10.0, 10.0], [5.0, 5.0])
    b = [tuple(round(v, 6) for v in g.bounds) for g in c]
    assert b == [(0, 5, 10, 10), (10, 5, 20, 10),      # row 0: the north
                 (0, 0, 10, 5), (10, 0, 20, 5)]


def test_majority_takes_the_value_covering_most_of_the_cell(tmp_path):
    # cell 0: 70 % code 1, 30 % code 2 -> 1.  cell 1: all code 2 -> 2
    lay = ov.Polygons.from_geojson(_layer(tmp_path, [
        (0, 0, 7, 10, {'c': 1}), (7, 0, 20, 10, {'c': 2})]), ['c'])
    got = ov.majority(_cells(), lay, 'c', fill=-1)
    assert got.tolist() == [1.0, 2.0]


def test_majority_is_by_area_not_by_the_centre(tmp_path):
    # the centre of cell 0 (x = 5) sits in code 2, which covers only 40 %
    lay = ov.Polygons.from_geojson(_layer(tmp_path, [
        (0, 0, 4, 10, {'c': 2}), (4, 0, 10, 10, {'c': 1}),
        (0, 0, 0, 0, {'c': 9})]), ['c'])
    got = ov.majority(ov.structured_cells(0, 0, [10.0], [10.0]), lay, 'c',
                      fill=-1)
    assert got.tolist() == [1.0]


def test_an_uncovered_cell_gets_the_fill(tmp_path):
    lay = ov.Polygons.from_geojson(_layer(tmp_path, [
        (0, 0, 10, 10, {'c': 3})]), ['c'])
    assert ov.majority(_cells(), lay, 'c', fill=-1).tolist() == [3.0, -1.0]


def test_area_mean_weights_by_the_covered_area(tmp_path):
    # cell 0: 1 m over 30 %, 2 m over 70 % -> 1.7 m. Cell 1 half covered,
    # 4 m -> 4 m: the uncovered half is not counted as a zero.
    lay = ov.Polygons.from_geojson(_layer(tmp_path, [
        (0, 0, 3, 10, {'t': 1.0}), (3, 0, 10, 10, {'t': 2.0}),
        (10, 0, 15, 10, {'t': 4.0})]), ['t'])
    got = ov.area_mean(_cells(), lay, 't', fill=-1)
    assert np.allclose(got, [1.7, 4.0])


def test_class_percent_is_against_the_whole_cell(tmp_path):
    # cell 0: 37 % class a; cell 1: 50 % class b, the rest bare
    lay = ov.Polygons.from_geojson(_layer(tmp_path, [
        (0, 0, 3.7, 10, {'s': 'a'}), (10, 0, 15, 10, {'s': 'b'})]), ['s'])
    got = ov.class_percent(_cells(), lay, 's', {'a': 1, 'b': 2}, 2)
    assert np.allclose(got, [[37.0, 0.0], [0.0, 50.0]])


def test_a_class_not_listed_is_ignored(tmp_path):
    lay = ov.Polygons.from_geojson(_layer(tmp_path, [
        (0, 0, 10, 10, {'s': 'z'})]), ['s'])
    got = ov.class_percent(_cells(), lay, 's', {'a': 1}, 1)
    assert np.allclose(got, 0.0)


def test_overlapping_polygons_past_100_percent_are_refused(tmp_path):
    lay = ov.Polygons.from_geojson(_layer(tmp_path, [
        (0, 0, 10, 10, {'s': 'a'}), (0, 0, 10, 10, {'s': 'b'})]), ['s'])
    with pytest.raises(ov.OverlayError) as exc:
        ov.class_percent(_cells(), lay, 's', {'a': 1, 'b': 2}, 2)
    assert 'overlap' in str(exc.value)


def test_a_missing_column_names_the_ones_there(tmp_path):
    with pytest.raises(ov.OverlayError) as exc:
        ov.Polygons.from_geojson(_layer(tmp_path, [
            (0, 0, 10, 10, {'SoilCode': 1})]), ['Soilcode'])
    assert 'Soilcode' in str(exc.value) and 'SoilCode' in str(exc.value)


def test_a_missing_file_points_at_the_converter(tmp_path):
    with pytest.raises(ov.OverlayError) as exc:
        ov.Polygons.from_geojson(str(tmp_path / 'nope.geojson'), ['c'])
    assert 'converter' in str(exc.value)


# ----------------------------------------------------------- La Mata, wired

needs_ds = pytest.mark.skipif(
    not os.path.exists(os.path.join(DS, 'inputSOILZONES.geojson')),
    reason='the La Mata dataset is not present')


@pytest.fixture(scope='module')
def lamata():
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    import marmites_config as mcfg
    import marmites_props as props
    cMF = ppMF.clsMF(MMutils.clsUTILITIES(verbose=0), MM_ws=DS, MM_ws_out=DS,
                     MF_ws=os.path.join(DS, 'MF_ws'),
                     MF_ini_fn='__inputMF_flopy_v3_2s1L.ini',
                     xllcorner=739300.0, yllcorner=4553050.0)
    cfg = mcfg.load_run_config(os.path.join(CODE, 'configs', 'lamata.toml'))
    act = (np.abs(np.asarray(cMF.ibound)) != 0).any(axis=0)
    return cMF, cfg, props, act


@needs_ds
def test_the_thickness_raster_gives_exactly_the_legacy_grid(lamata):
    """The panel names the raster, and a raster wins: nothing changes."""
    cMF, cfg, props, _act = lamata
    assert cfg.soil.thickness.producer() == 'raster'
    got = props.soil_grid(cfg, cMF, DS, 'thickness')
    want = cMF.cPROCESS.inputEsriAscii(grid_fn='inputSOILthick.asc',
                                       datatype=float)
    assert np.array_equal(got, np.asarray(want, dtype=float))


@needs_ds
def test_the_zone_raster_gives_exactly_the_legacy_grid(lamata):
    """Naming the legacy raster keeps the legacy zones, exactly -- the way
    to keep them if the polygon majority is not wanted."""
    import copy
    cMF, cfg, props, _act = lamata
    c = copy.deepcopy(cfg)
    c.soil.zones.raster = 'inputSOILzones.asc'
    got = props.soil_grid(c, cMF, DS, 'zones', kind='int')
    want = cMF.cPROCESS.inputEsriAscii(grid_fn='inputSOILzones.asc',
                                       datatype=int)
    assert np.array_equal(got, np.asarray(want))


@needs_ds
def test_the_zone_polygons_agree_with_the_legacy_grid_but_for_the_edges(
        lamata):
    """The panel names the POLYGONS, by majority of area. On La Mata that
    gives the legacy zone in 1828 of the 1954 active cells; the 126 others
    are cells a zone boundary crosses, where area-majority and the legacy
    rasterisation disagree. Pinned, so a change to either shows up."""
    cMF, cfg, props, act = lamata
    assert cfg.soil.zones.producer() == 'layer'
    got = props.soil_grid(cfg, cMF, DS, 'zones', kind='int')
    old = np.asarray(cMF.cPROCESS.inputEsriAscii(
        grid_fn='inputSOILzones.asc', datatype=int))
    assert set(np.unique(got[act]).tolist()) <= {1, 2, 3}
    assert int((got[act] == old[act]).sum()) == 1828


def test_a_thickness_layer_other_than_the_zone_layer_is_refused():
    import types
    import marmites_config as mcfg
    import marmites_props as props
    cfg = mcfg.RunConfig.from_dict({})
    cfg.soil.thickness = mcfg.VectorSource(layer='other.shp', column='d')
    cMF = types.SimpleNamespace(nrow=1, ncol=1, hnoflo=9999.999)
    with pytest.raises(props.PropertyError) as exc:
        props.soil_grid(cfg, cMF, DS, 'thickness')
    assert 'converter exports only' in str(exc.value)


# ------------------------------------------------------ the vegetation cover
# From the Soil panel's vegetation layer, class column and class table. The
# polygons are the reference: grass ~89 %, where the retired raster said 25.

def _window(r0, r1, c0, c1, nrow=65, xll=739300.0, yll=4553050.0, d=50.0):
    """A sub-grid of La Mata's 50 m grid, as a cMF-shaped namespace."""
    import types
    return types.SimpleNamespace(
        nrow=r1 - r0, ncol=c1 - c0, hnoflo=9999.999,
        xllcorner=xll + d * c0, yllcorner=yll + d * (nrow - r1),
        delr=np.full(c1 - c0, d), delc=np.full(r1 - r0, d))


needs_veg = pytest.mark.skipif(
    not os.path.exists(os.path.join(DS, 'inputVEG.geojson')),
    reason='the La Mata vegetation layer is not present')


@needs_veg
def test_the_trees_come_out_as_they_always_were_and_grass_as_mapped(
        lamata, tmp_path):
    """THE acceptance test, on a 15 x 15 window of the real grid. The tree
    types reproduce the retired rasters to 0.02 points -- which proves the
    overlay; grass is the polygons' own ~90 %, which is the ruling."""
    _cMF, cfg, props, _act = lamata
    r0, r1, c0, c1 = 25, 40, 20, 35
    win = _window(r0, r1, c0, c1)
    got = props.veg_cover(cfg, win, DS, 3, cache_dir=str(tmp_path),
                          verbose=False)
    assert got.shape == (3, 15, 15)
    full = lamata[0]
    for k in (2, 3):
        old = np.asarray(full.cPROCESS.convASCIIraster2array(
            os.path.join(DS, 'inputVEG%darea.asc' % k),
            np.zeros((full.nrow, full.ncol))), dtype=float)[r0:r1, c0:c1]
        assert np.abs(got[k - 1] - old).max() <= 0.02, 'type %d moved' % k
    assert got[0].mean() > 50.0, 'grass is not the polygons\' cover'
    assert got.sum(axis=0).max() <= 100.0 + 1e-4


@needs_veg
def test_the_overlay_is_cached_on_the_content(lamata, tmp_path):
    _cMF, cfg, props, _act = lamata
    win = _window(30, 33, 25, 28)
    a = props.veg_cover(cfg, win, DS, 3, cache_dir=str(tmp_path),
                        verbose=False)
    files = list(tmp_path.glob('veg_cover_*.npz'))
    assert len(files) == 1, 'the overlay was not cached'
    b = props.veg_cover(cfg, win, DS, 3, cache_dir=str(tmp_path),
                        verbose=False)
    assert np.array_equal(a, b)


def _veg_cfg(tmp_path, feats, classes):
    """A configuration and dataset folder holding a toy vegetation layer."""
    import marmites_config as mcfg
    d = tmp_path / 'ds'
    d.mkdir()
    _layer(d, feats, name='inputVEG.geojson')
    cfg = mcfg.RunConfig.from_dict({})
    cfg.soil.veg_column = 'Species'
    cfg.soil.veg_class = [mcfg.VegetationClass(code=c, veg=v)
                          for c, v in classes]
    return cfg, str(d)


def test_a_class_the_table_does_not_map_is_refused(tmp_path):
    """Every class must be defined: an unmapped one used to be skipped, its
    area left unvegetated without a word."""
    import marmites_props as props
    cfg, ds = _veg_cfg(tmp_path, [(0, 0, 10, 10, {'Species': 'g'}),
                                  (10, 0, 20, 10, {'Species': 'x'})],
                       [('g', 1)])
    win = _window(0, 1, 0, 2, nrow=1, xll=0.0, yll=0.0, d=10.0)
    with pytest.raises(props.PropertyError) as exc:
        props.veg_cover(cfg, win, ds, 1, verbose=False)
    assert "'x'" in str(exc.value) and 'defined' in str(exc.value)


def test_a_mapped_type_that_does_not_exist_is_refused(tmp_path):
    import marmites_props as props
    cfg, ds = _veg_cfg(tmp_path, [(0, 0, 10, 10, {'Species': 'g'})],
                       [('g', 4)])
    win = _window(0, 1, 0, 1, nrow=1, xll=0.0, yll=0.0, d=10.0)
    with pytest.raises(props.PropertyError) as exc:
        props.veg_cover(cfg, win, ds, 3, verbose=False)
    assert 'there are 3' in str(exc.value)


def test_the_cover_is_per_type_percent_of_the_cell(tmp_path):
    import marmites_props as props
    cfg, ds = _veg_cfg(tmp_path, [(0, 0, 6, 10, {'Species': 'g'}),
                                  (6, 0, 10, 10, {'Species': 'i'})],
                       [('g', 1), ('i', 2)])
    win = _window(0, 1, 0, 1, nrow=1, xll=0.0, yll=0.0, d=10.0)
    got = props.veg_cover(cfg, win, ds, 2, verbose=False)
    assert np.allclose(got[:, 0, 0], [60.0, 40.0])


def test_the_run_passes_the_cover_in_rather_than_reading_the_rasters():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'),
               encoding='utf-8').read()
    assert 'props.veg_cover(' in src
    assert 'gridVEGarea=_veg' in src, 'the cover is computed and not used'
