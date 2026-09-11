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
    c = cfg({'grid': {'override': {'enable': True, 'xllcorner': 739300.0,
                                   'yllcorner': 4553050.0,
                                   'nrow': 65, 'ncol': 60}}})
    assert c.grid.override.enable is True
    assert (c.grid.override.nrow, c.grid.override.ncol) == (65, 60)


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
