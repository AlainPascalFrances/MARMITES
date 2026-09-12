# -*- coding: utf-8 -*-
"""WP1d -- the four front-end panels in the configuration schema.

Panel 1 the grid (the catchment polygon first, everything wrapped onto it),
panel 2 the surface (MMsurf), panel 3 the soil (MMsoil) and MODFLOW, panel 4
plotting. The master switches live in one `[run]` block so the shape of a run
is visible in one place.

What these tests are really guarding is the two traps the old ini set: a
parameter that looks right and is not (`kt_s` stored as 1/s), and a parameter
that is read but never reaches the model (`iniMonthHydroYear`).
"""

import importlib.util
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, '..', '..')
REF_CONFIG = os.path.join(REPO, 'code', 'configs', 'lamata.toml')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cfgmod = _load('marmites_config_panels',
               os.path.join(REPO, 'code', 'marmites_config.py'))


def cfg(d):
    return cfgmod.RunConfig.from_dict(d)


# ---------------------------------------------------------------- switches

def test_master_switches_live_in_one_run_block():
    c = cfg({})
    assert (c.run.surface, c.run.model, c.run.plot) == (False, True, True)
    c = cfg({'run': {'surface': True, 'plot': False},
             'surface': {'station': [{'name': 'A'}],
                         'vegetation': [{'name': 'grassMU'}],
                         'soil': [{'name': 'alluvium'}]}})
    assert c.run.surface is True and c.run.plot is False
    assert c.run.model is True          # MMsoil and MF6 are not separable


# ------------------------------------------------------------ panel 1 grid

def test_panel1_domain_defaults():
    c = cfg({})
    assert c.grid.boundary.endswith('.shp')
    assert c.grid.cell_size > 0
    assert c.grid.override.enable is False


@pytest.mark.parametrize('bad, frag', [
    ({'boundary': 'lm_lim.tif'}, 'shapefile'),
    ({'boundary': ''}, 'catchment polygon'),
    ({'cell_size': 0.0}, 'cell_size'),
    ({'buffer': -1.0}, 'buffer'),
    ({'override': {'enable': True}}, 'nrow and ncol'),
])
def test_panel1_guards(bad, frag):
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'grid': bad})
    assert frag in str(e.value)


def test_legacy_grid_override_carries_the_old_rectangle():
    """Deriving the grid from the polygon loses every baseline unless the old
    one can still be rebuilt (WP1d, C.1)."""
    c = cfg({'grid': {'kind': 'structured',
                      'override': {'enable': True, 'xllcorner': 739300.0,
                                   'yllcorner': 4553050.0,
                                   'nrow': 65, 'ncol': 60}}})
    assert c.grid.override.enable is True
    assert (c.grid.override.nrow, c.grid.override.ncol) == (65, 60)


@pytest.mark.parametrize('kind', ['voronoi', 'quadtree'])
def test_the_override_is_refused_on_a_mesh(kind):
    """WP1d panel 1, D3. The override reproduces an EXISTING rectangle, and
    model_rectangle honours it whatever the kind -- so a box left on from a
    structured run would silently override a mesh's domain. The panel hides
    it for these kinds, and hiding is not disabling, so the schema refuses."""
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'grid': {'kind': kind,
                      'override': {'enable': True, 'nrow': 65, 'ncol': 60}}})
    assert 'grid.override.enable' in str(e.value)


def _vor(**kw):
    base = {'cell_far': 100.0, 'cell_near_stream': 40.0,
            'stream_buffer': 1000.0, 'grade_ratio': 1.5}
    base.update(kw)
    return cfg({'grid': {'kind': 'voronoi', 'voronoi': base}}).grid.voronoi


def test_a_band_is_as_wide_as_the_cells_it_carries():
    """The constraint that was missing. Spacing the bands evenly across the
    corridor asked for size changes closer together than the cells were
    wide, and Triangle honoured the band boundaries instead -- which is how
    a corridor meant to hold 15 m cells came out at 3.7 m."""
    v = _vor()
    bands = v.graded_bands()
    assert bands[0] == (40.0, 40.0)             # first band, one cell wide
    prev_d = 0.0
    for d, s in bands:
        assert d - prev_d <= s + 1e-6, 'a band is narrower than its own cells'
        prev_d = d
    sizes = [s for _d, s in bands]
    assert all(b / a <= 1.5 + 1e-9 for a, b in zip(sizes, sizes[1:]))
    assert max(sizes) <= v.cell_far


def test_the_bands_are_derived_whatever_the_file_says():
    v = _vor(trans_levels=[10.0, 20.0, 40.0, 70.0])
    assert v.trans_levels == [d for d, _s in v.graded_bands()]
    assert 10.0 not in v.trans_levels


def test_a_corridor_too_narrow_to_grade_is_clamped_and_reported():
    """The mesh is still built -- the step at the corridor edge is simply
    bigger than the ratio asked for -- and corridor_needed says what it would
    take to finish the grade. The SHIPPED defaults are in this case: 40 to
    100 m at 1.5 needs 190 m and the default corridor is 60."""
    v = _vor(stream_buffer=60.0)
    assert v.corridor_needed() == 190.0
    assert v.trans_levels[-1] == 60.0           # clamped to the corridor
    assert max(s for _d, s in v.graded_bands()) < v.cell_far


def test_a_coarser_grade_ratio_asks_for_fewer_bands():
    """Compared on a corridor wide enough not to clamp either of them."""
    assert len(_vor(grade_ratio=1.2).trans_levels) > \
        len(_vor(grade_ratio=2.5).trans_levels)


@pytest.mark.parametrize('ratio', [1.0, 0.5, 3.5, -1.0])
def test_an_impossible_grade_ratio_is_refused(ratio):
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'grid': {'kind': 'voronoi', 'voronoi': {'grade_ratio': ratio}}})
    assert 'grade_ratio' in str(e.value)


def test_a_ratio_needing_more_bands_than_allowed_is_refused():
    """Capping the bands silently would mean the mesh does not grade the way
    the file says it does."""
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'grid': {'kind': 'voronoi',
                      'voronoi': {'cell_near_stream': 2.0, 'cell_far': 200.0,
                                  'grade_ratio': 1.05}}})
    assert 'transition bands' in str(e.value)


def test_switching_the_refinement_off_clears_the_corridor():
    """Panel 1, D4: with no refinement there is no corridor, so the settings
    that describe one do not stay in the file looking live."""
    v = cfg({'grid': {'kind': 'voronoi',
                      'voronoi': {'stream_refine': False,
                                  'cell_near_stream': 40.0,
                                  'stream_buffer': 60.0}}}).grid.voronoi
    assert v.trans_levels == []
    assert v.cell_near_stream == 0.0 and v.stream_buffer == 0.0


# --------------------------------------------------------- panel 2 surface

def test_surface_parameter_tables_are_arrays_of_tables():
    c = cfg({'surface': {
        'station': [{'name': 'A', 'x': 1.0, 'y': 2.0},
                    {'name': 'B', 'x': 3.0, 'y': 4.0}],
        'vegetation': [{'name': 'grassMU'}],
        'soil': [{'name': 'alluvium', 'porosity': 0.4}]}})
    assert [s.name for s in c.surface.station] == ['A', 'B']
    assert c.surface.station[1].x == 3.0
    assert c.surface.vegetation[0].name == 'grassMU'
    assert c.surface.soil[0].porosity == 0.4


def test_surface_table_rejects_an_unknown_key():
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'surface': {'station': [{'lat': 41.0}]}})
    assert 'unknown key' in str(e.value)


def test_surface_table_rejects_a_bare_table():
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'surface': {'station': {'name': 'A'}}})
    assert 'list of tables' in str(e.value)


def test_kt_s_is_the_slope_not_its_reciprocal():
    """The old ini stored 1/s -- values above 20 -- and the driver inverted
    it. Typing that number here must be refused, not silently used."""
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'surface': {'vegetation': [{'name': 'Qilex',
                                         'kt_s': 21.98151201}]}})
    assert 'not the 1/s' in str(e.value)
    c = cfg({'surface': {'vegetation': [{'name': 'Qilex', 'kt_s': 0.0455}]}})
    assert c.surface.vegetation[0].kt_s == 0.0455


def test_running_the_surface_needs_its_parameter_tables():
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'run': {'surface': True}})
    msg = str(e.value)
    assert 'surface.station is empty' in msg
    assert 'surface.vegetation is empty' in msg


def test_irrigation_needs_fields_and_a_crop():
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'surface': {'irrigation': True}})
    msg = str(e.value)
    assert 'nfield' in msg and 'crop' in msg


# ------------------------------------------------------------ vector sources

def test_vector_source_precedence_is_explicit():
    """raster > layer > value, and 'nothing set' is an error, not a zero."""
    c = cfg({'soil': {'thickness': {'raster': 'inputSOILthick.asc',
                                    'layer': 'Soil_type.shp',
                                    'column': 'SOILthick'}}})
    assert c.soil.thickness.producer() == 'raster'
    c = cfg({'soil': {'thickness': {'layer': 'Soil_type.shp',
                                    'column': 'SOILthick'}}})
    assert c.soil.thickness.producer() == 'layer'
    c = cfg({'soil': {'thickness': 0.75}})
    assert c.soil.thickness.producer() == 'value'
    assert c.soil.thickness.value == 0.75


def test_vector_source_with_nothing_set_is_refused():
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'soil': {'thickness': {}}})
    assert 'nothing to read' in str(e.value)


def test_vector_source_from_a_bare_string():
    assert cfg({'soil': {'zones': 'Soil_type.shp'}}).soil.zones.producer() \
        == 'layer'
    assert cfg({'soil': {'thickness': 'thick.asc'}}).soil.thickness.producer() \
        == 'raster'


def test_vegetation_class_index_must_exist():
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'surface': {'vegetation': [{'name': 'grassMU'}]},
             'soil': {'veg_class': [{'code': 'i', 'veg': 2}]}})
    assert 'not a vegetation index' in str(e.value)


# --------------------------------------------------------- panel 4 plotting

def test_hydro_year_start_is_a_real_field_now():
    """It was in the MM ini but never reached the model: every call site fell
    back to getattr(cMF, 'iniMonthHydroYear', 10), so October was hardcoded."""
    assert cfg({'postproc': {'hydro_year_start': 9}}).postproc.hydro_year_start \
        == 9
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'postproc': {'hydro_year_start': 13}})
    assert 'hydro_year_start' in str(e.value)


def test_wb_unit_is_guarded():
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg({'postproc': {'wb_unit': 'month'}})
    assert 'wb_unit' in str(e.value)


# ------------------------------------------------------- the reference file

def test_reference_config_carries_the_four_panels():
    if not os.path.exists(REF_CONFIG):
        pytest.skip('code/configs/lamata.toml not present')
    c = cfgmod.load_run_config(REF_CONFIG)
    assert c.grid.boundary == 'lm_lim.shp' and c.grid.crs_epsg == 23029
    assert len(c.surface.station) == 1
    assert c.surface.station[0].x == 739508.0      # project CRS, not degrees
    assert [v.name for v in c.surface.vegetation] == ['grassMU', 'Qilex', 'Qpyr']
    assert [s.name for s in c.surface.soil] == ['alluvium', 'regolith', 'outcrop']
    assert c.soil.zones.layer == 'Soil_type.shp'
    assert c.soil.thickness.producer() == 'raster'
    assert {v.code: v.veg for v in c.soil.veg_class} == {'g': 1, 'i': 2, 'p': 3}
    assert c.postproc.hydro_year_start == 10


def test_surface_tables_survive_the_resolved_round_trip(tmp_path):
    """The resolved file written into the run folder is the provenance, so the
    parameter tables must come back byte-for-byte in meaning."""
    if not os.path.exists(REF_CONFIG):
        pytest.skip('code/configs/lamata.toml not present')
    c = cfgmod.load_run_config(REF_CONFIG)
    out = str(tmp_path / 'resolved.toml')
    c.write_toml(out)
    back = cfgmod.load_run_config(out)
    assert back.config_hash() == c.config_hash()
    assert len(back.surface.vegetation) == 3
    assert back.surface.vegetation[1].root_depth == 15.0
    assert back.soil.veg_class[2].code == 'p'
