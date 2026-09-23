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
