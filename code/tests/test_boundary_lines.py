# -*- coding: utf-8 -*-
"""DRN and GHB placed by a LINE (the modeller's design, 2026-09-24).

The panel names a line shapefile; every cell of the grid the run USES that
the line crosses AND that has a face on the catchment's external boundary
becomes a boundary cell, on every layer listed. A drain sits at the base of
its layer; the conductance is the panel's number per metre of boundary face
(grid-independent) or per cell.
"""

import importlib.util
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'app')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

shapely = pytest.importorskip('shapely')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


props = _load('marmites_props_lines', os.path.join(CODE, 'ppMF6',
                                                   'marmites_props.py'))
import marmites_config as mcfg                                # noqa: E402
import marmites_vector as mv                                  # noqa: E402

# a 4 x 4 grid of 10 m cells, x 0..40, y 0..40; row 0 is the NORTH row
N, D = 4, 10.0


def _cMF(active=None, nlay=2):
    act = np.ones((N, N), dtype=int) if active is None else np.asarray(active)
    ib = np.stack([act] * nlay)
    botm = np.stack([np.full((N, N), 90.0 - 20.0 * l) for l in range(nlay)])
    return SimpleNamespace(nlay=nlay, nrow=N, ncol=N, delr=np.full(N, D),
                           delc=np.full(N, D), xllcorner=0.0, yllcorner=0.0,
                           ibound=ib, outcropL=act.copy(), botm=botm,
                           hnoflo=-9999.0)


def _cells():
    return props.model_cells(_cMF())


def _sel(line, active=None):
    act = np.ones((N, N), dtype=bool) if active is None else active
    return props.boundary_cells_on_line(_cells(), act, line)


# ------------------------------------------------------------ the selection
def test_a_line_along_the_west_edge_takes_the_cells_on_it():
    """x = 0 from y = 20 to 40: the two NORTH-WEST cells, rows 0 and 1."""
    idx, face = _sel(shapely.LineString([(0, 20), (0, 40)]))
    assert idx.tolist() == [0, 4]
    assert np.allclose(face, [20.0, 10.0]), 'the corner cell has 2 faces'


def test_the_cell_past_the_lines_end_is_not_on_it():
    """It touches the line's end vertex only."""
    idx, _ = _sel(shapely.LineString([(0, 30), (0, 40)]))
    assert idx.tolist() == [0]


def test_a_line_across_the_edge_takes_the_boundary_cell_it_crosses():
    idx, face = _sel(shapely.LineString([(-5, 15), (15, 15)]))
    assert idx.tolist() == [8], 'row 2, col 0; not its inner neighbour'
    assert face.tolist() == [10.0]


def test_an_interior_line_takes_nothing():
    idx, _ = _sel(shapely.LineString([(12, 12), (28, 28)]))
    assert idx.size == 0


def test_a_line_a_few_millimetres_off_the_edge_still_counts():
    """The converter rounds to the centimetre."""
    idx, _ = _sel(shapely.LineString([(-0.004, 20), (-0.004, 40)]))
    assert idx.tolist() == [0, 4]


def test_an_inactive_island_is_not_the_catchment_boundary():
    act = np.ones((N, N), dtype=bool)
    act[1, 1] = False                      # a hole inside
    # a line through the island touches cells around it, none on the edge
    idx, _ = _sel(shapely.LineString([(10, 25), (20, 25)]), act)
    assert idx.size == 0


def test_the_boundary_follows_the_active_cells_not_the_grid():
    """Column 0 inactive: the catchment edge is x = 10."""
    act = np.ones((N, N), dtype=bool)
    act[:, 0] = False
    idx, face = _sel(shapely.LineString([(10, 30), (10, 40)]), act)
    assert idx.tolist() == [1] and face.tolist() == [20.0]


# ------------------------------------------------------------ the records
def _cfg(name='drn', **kw):
    cfg = mcfg.RunConfig.from_dict({})
    pkg = getattr(cfg, name)
    pkg.enable = True
    pkg.layers = [1, 2]
    pkg.line = 'outlet.shp'
    pkg.cond = mcfg.VectorSource(value=0.5)
    if name == 'drn':
        pkg.at_layer_base = True
    else:
        pkg.head = mcfg.VectorSource(value=95.0)
    for k, v in kw.items():
        setattr(pkg, k, v)
    return cfg


def _dataset(tmp_path, name='drn', coords=((0, 20), (0, 40))):
    mv.write_geojson(str(tmp_path / props.LINE_TABLES[name]), 'line',
                     [[list(coords)]], [{}])
    return str(tmp_path)


def test_drains_sit_at_the_base_of_each_layer_per_metre_of_face(tmp_path):
    cMF = _cMF()
    done = props.apply_line_boundaries(_cfg(), cMF, _dataset(tmp_path),
                                       verbose=False)
    assert done == [('drn', 2, 4)]
    recs = sorted(cMF.layer_row_column_elevation_cond[0])
    assert [r[:3] for r in recs] == [[0, 0, 0], [0, 1, 0], [1, 0, 0],
                                     [1, 1, 0]]
    assert [round(r[3], 2) for r in recs] == [90.01, 90.01, 70.01, 70.01]
    # 0.5 m2/d per metre x 20 m (the corner) and x 10 m
    assert [r[4] for r in recs] == [10.0, 5.0, 10.0, 5.0]
    assert cMF.drn_cond_array.shape == (2, N, N)
    assert cMF.drn_cond_array[0, 0, 0] == 10.0


def test_per_cell_is_the_same_number_everywhere(tmp_path):
    cMF = _cMF()
    props.apply_line_boundaries(_cfg(cond_per='cell'), cMF,
                                _dataset(tmp_path), verbose=False)
    assert {r[4] for r in cMF.layer_row_column_elevation_cond[0]} == {0.5}


def test_a_layer_absent_under_the_line_gets_no_drain(tmp_path):
    cMF = _cMF()
    cMF.ibound[1, 0, 0] = 0
    props.apply_line_boundaries(_cfg(), cMF, _dataset(tmp_path),
                                verbose=False)
    got = {(r[0], r[1]) for r in cMF.layer_row_column_elevation_cond[0]}
    assert (1, 0) not in got and (0, 0) in got


def test_a_drain_off_the_base_takes_the_elevation_given(tmp_path):
    cMF = _cMF()
    cfg = _cfg(at_layer_base=False)
    cfg.drn.elevation = mcfg.VectorSource(value=88.0)
    props.apply_line_boundaries(cfg, cMF, _dataset(tmp_path), verbose=False)
    assert {r[3] for r in cMF.layer_row_column_elevation_cond[0]} == {88.0}


def test_the_ghb_holds_its_head_on_its_own_line(tmp_path):
    cMF = _cMF()
    ds = _dataset(tmp_path, 'ghb', ((40, 0), (40, 10)))    # the south-east
    done = props.apply_line_boundaries(_cfg('ghb'), cMF, ds, verbose=False)
    assert done == [('ghb', 1, 2)]
    recs = cMF.layer_row_column_head_cond[0]
    assert {(r[1], r[2]) for r in recs} == {(3, 3)}
    assert {r[3] for r in recs} == {95.0}
    assert {r[4] for r in recs} == {10.0}          # 0.5 x (10 + 10) m


def test_a_line_that_misses_the_boundary_is_an_error(tmp_path):
    with pytest.raises(props.PropertyError, match='crosses no cell'):
        props.apply_line_boundaries(
            _cfg(), _cMF(), _dataset(tmp_path, coords=((12, 12), (28, 28))),
            verbose=False)


def test_an_unconverted_line_says_how_to_convert_it(tmp_path):
    with pytest.raises(props.PropertyError, match='Cartography'):
        props.apply_line_boundaries(_cfg(), _cMF(), str(tmp_path),
                                    verbose=False)


def test_a_polygon_is_not_a_line(tmp_path):
    mv.write_geojson(str(tmp_path / 'inputDRN.geojson'), 'polygon',
                     [[[(0, 0), (1, 0), (1, 1), (0, 0)]]], [{}])
    with pytest.raises(props.PropertyError, match='LINE'):
        props.apply_line_boundaries(_cfg(), _cMF(), str(tmp_path),
                                    verbose=False)


def test_with_a_line_the_legacy_rasters_are_not_read():
    """apply_boundaries only announces it: the records come on the final
    grid, and the 50 m rasters would put them on the 50 m cells."""
    cMF = _cMF()
    cMF.cPROCESS = None                    # would fail if a raster were read
    done = props.apply_boundaries(_cfg(), cMF, '/nowhere', verbose=False)
    assert done == [] and cMF.drn_yn == 1
    assert cMF.layer_row_column_elevation_cond == {0: []}


def test_no_line_means_nothing_to_place(tmp_path):
    cfg = _cfg()
    cfg.drn.line = ''
    assert props.apply_line_boundaries(cfg, _cMF(), str(tmp_path)) == []


# ------------------------------------------------------ the configuration
def test_a_drain_on_a_line_at_the_base_needs_no_elevation():
    cfg = _cfg()
    cfg.drn.elevation = mcfg.VectorSource()
    assert not [e for e in cfg.problems() if e.startswith('drn.')]


def test_the_line_must_be_a_shapefile_and_the_unit_known():
    cfg = _cfg(line='outlet.geojson', cond_per='hectare')
    errs = ' '.join(cfg.problems())
    assert 'drn.line must name a shapefile' in errs
    assert "drn.cond_per must be 'length' or 'cell'" in errs


def test_the_default_is_per_metre():
    cfg = mcfg.RunConfig.from_dict({})
    assert cfg.drn.cond_per == 'length' and cfg.ghb.cond_per == 'length'
    assert cfg.drn.line == '' and cfg.ghb.line == ''


def test_the_panel_asks_for_the_line_and_the_unit():
    from lib import schema
    for name in ('drn', 'ghb'):
        rows = schema.DRN_ROWS if name == 'drn' else schema.GHB_ROWS
        flat = [r for row in rows for r in row]
        assert '%s.line' % name in flat and '%s.cond_per' % name in flat
        assert '%s.line' % name in schema.BOUNDARY_FILES
        assert schema.CHOICES['%s.cond_per' % name]() == ['length', 'cell']
    page = open(os.path.join(CODE, 'app', 'pages',
                             '4_4_-_Unsaturated_zone_and_groundwater.py'),
                encoding='utf-8').read()
    assert page.count('files=schema.BOUNDARY_FILES') == 2


def test_the_converter_writes_each_line_where_the_run_reads_it():
    conv = _load('gis_to_dataset_lines',
                 os.path.join(CODE, 'tools', 'gis_to_dataset.py'))
    stems = {role: stem for role, _d, _c, stem in conv.VECTOR_LAYERS}
    assert stems['drn_line'] + '.geojson' == props.LINE_TABLES['drn']
    assert stems['ghb_line'] + '.geojson' == props.LINE_TABLES['ghb']
    from lib import dataset_state
    cfg = mcfg.load_run_config(os.path.join(CODE, 'configs', 'lamata.toml'))
    cfg.drn.line = 'outlet.shp'
    assert ('inputDRN.geojson', 'outlet.shp', []) in dataset_state.tables(cfg)


def test_the_driver_places_them_on_the_final_grid():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert src.index('_apply_dem(cMF, cfg, DS') \
        < src.index('props.apply_line_boundaries(cfg, cMF, DS)') \
        < src.index('mm.build_cell_list(cMF)')


def test_with_a_line_the_conductance_is_the_input_box():
    """The legacy raster holds 0.035 per 50 m CELL; read per metre of face
    it was 50 times the legacy boundary (2026-09-24: DRN out 11.8 mm/yr
    against the legacy 4.4). With a line, only a value is accepted."""
    cfg = _cfg()
    cfg.drn.cond = mcfg.VectorSource(raster='MF_ws/drn_cond_l%d.asc')
    errs = [e for e in cfg.problems() if e.startswith('drn.cond')]
    assert errs and '0.0007 per metre' in errs[0]
    cfg.drn.line = ''                       # the legacy rule reads rasters
    assert not [e for e in cfg.problems() if e.startswith('drn.cond')]
