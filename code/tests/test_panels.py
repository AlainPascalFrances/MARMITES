# -*- coding: utf-8 -*-
"""WP1d -- the four panels: what they show, and editing what they show.

The pages are thin views over ``app/lib/schema.py`` and ``app/lib/editor.py``,
so what could actually be wrong is tested here rather than found in a browser:
a configuration field no panel shows, a panel naming a section that does not
exist, an edit that would write an invalid file, and a table row that silently
loses a column.
"""

import importlib.util
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
REPO = os.path.abspath(os.path.join(CODE, '..'))
REF = os.path.join(CODE, 'configs', 'lamata.toml')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


cfgmod = _load('marmites_config_p', os.path.join(CODE, 'marmites_config.py'))
schema = _load('app_schema', os.path.join(CODE, 'app', 'lib', 'schema.py'))
editor = _load('app_editor', os.path.join(CODE, 'app', 'lib', 'editor.py'))


@pytest.fixture
def cfg():
    if not os.path.exists(REF):
        pytest.skip('code/configs/lamata.toml not present')
    return cfgmod.load_run_config(REF)


# ---------------------------------------------------------------- panels

def test_the_panels_are_the_modellers_order():
    nums = [p[0] for p in schema.PANELS]
    assert nums == [0, 1, 2, 3, 4]
    titles = [p[1] for p in schema.PANELS]
    assert titles == ['Overview', 'Grid', 'Surface', 'Model', 'Plots']


def test_the_master_switches_are_real_run_keys(cfg):
    """A switch that is not a key the driver reads would be decoration."""
    for _n, _t, _i, switch, _s, _b in schema.PANELS:
        if switch is None:
            continue
        section, key = switch.split('.')
        assert hasattr(getattr(cfg, section), key), switch
        assert isinstance(getattr(getattr(cfg, section), key), bool)


def test_only_surface_model_and_plots_carry_a_switch():
    by = {p[1]: p[3] for p in schema.PANELS}
    assert by['Overview'] is None and by['Grid'] is None
    assert by['Surface'] == 'run.surface'
    assert by['Model'] == 'run.model'
    assert by['Plots'] == 'run.plot'


def test_every_panel_section_exists_in_the_schema(cfg):
    for _n, title, _i, _sw, sections, _b in schema.PANELS:
        for s in sections:
            assert hasattr(cfg, s), '%s names a missing section %r' % (title, s)


def test_every_configuration_section_is_reachable(cfg):
    """A section in no panel is a setting no one can find."""
    shown = {s for p in schema.PANELS for s in p[4]}
    # meta/paths/run live on the Overview panel; pest and ui are not modelling
    # settings and are edited in the TOML.
    exempt = {'meta', 'paths', 'run', 'pest', 'ui'}
    missing = sorted(set(cfgmod._SECTIONS) - shown - exempt)
    assert not missing, 'no panel shows: %s' % ', '.join(missing)


def test_panel_of_finds_the_right_panel():
    assert schema.panel_of('grid') == 1
    assert schema.panel_of('surface') == 2
    assert schema.panel_of('sfr') == 3
    assert schema.panel_of('postproc') == 4
    assert schema.panel_of('nowhere') is None


# ---------------------------------------------------------------- fields

def test_every_field_of_every_panel_has_a_label(cfg):
    """'A friendly explanation' was the brief; this is what enforces it."""
    bare = []
    for _n, _t, _i, _sw, sections, _b in schema.PANELS:
        for s in sections:
            for dotted, _v in schema.fields_of(cfg, s):
                label, _u, _h = schema.describe(dotted)
                if label == dotted.split('.')[-1]:
                    bare.append(dotted)
    assert not bare, 'no label for: %s' % ', '.join(sorted(bare))


def test_describe_falls_back_to_the_table_row_form():
    """surface.vegetation[2].kt_s and vegetation.kt_s are the same field."""
    a = schema.describe('vegetation.kt_s')
    b = schema.describe('surface.vegetation[2].kt_s')
    assert a == b
    assert 'old ini stored 1/s' in a[2]


def test_describe_never_raises_on_an_unknown_key():
    label, units, help_ = schema.describe('nowhere.at.all')
    assert label == 'all' and units == '' and help_ == ''


def test_units_are_given_where_they_exist():
    assert schema.describe('grid.cell_size')[1] == 'm'
    assert schema.describe('seep.cond')[1] == 'm2/d'
    assert schema.describe('postproc.hydro_year_start')[1] == 'month'


def test_fields_of_flattens_nested_blocks(cfg):
    keys = dict(schema.fields_of(cfg, 'grid'))
    assert 'grid.cell_size' in keys
    assert 'grid.voronoi.cell_far' in keys       # nested, flattened
    assert 'grid.override.nrow' in keys


def test_fields_of_skips_the_arrays_of_tables(cfg):
    keys = dict(schema.fields_of(cfg, 'surface'))
    assert 'surface.meteo_ts' in keys
    assert not any(k.startswith('surface.station') for k in keys)


def test_fields_of_rejects_a_missing_section(cfg):
    with pytest.raises(schema.PanelError):
        schema.fields_of(cfg, 'nope')


def test_the_surface_tables_are_declared():
    paths = [t[0] for t in schema.TABLES]
    assert paths == ['surface.station', 'surface.vegetation', 'surface.crop',
                     'surface.soil', 'soil.veg_class']
    assert all(t[1] in (2, 3) for t in schema.TABLES)


# ---------------------------------------------------------------- editing

def test_changes_reports_only_what_moved(cfg):
    assert editor.changes(cfg, {'run.nsp': cfg.run.nsp}) == []
    assert editor.changes(cfg, {'run.nsp': 30}) == ['run.nsp=30']


def test_a_float_that_survived_a_widget_is_not_a_change(cfg):
    assert editor.changes(cfg, {'seep.cond': float(cfg.seep.cond)}) == []


def test_changes_on_a_missing_key_is_an_error(cfg):
    with pytest.raises(editor.EditError):
        editor.changes(cfg, {'run.nope': 1})


def test_apply_changes_leaves_the_original_alone(cfg):
    was = cfg.run.nsp
    new, applied = editor.apply_changes(cfg, {'run.nsp': 42})
    assert applied == ['run.nsp=42']
    assert new.run.nsp == 42
    assert cfg.run.nsp == was


def test_an_edit_that_would_be_invalid_is_refused(cfg):
    with pytest.raises(cfgmod.ConfigError) as e:
        editor.apply_changes(cfg, {'postproc.hydro_year_start': 13})
    assert 'hydro_year_start' in str(e.value)


def test_save_writes_the_resolved_configuration(tmp_path, cfg):
    out = str(tmp_path / 'edited.toml')
    applied, digest = editor.save(cfg, out, {'run.nsp': 30, 'seep.cond': 5000.0})
    assert sorted(applied) == ['run.nsp=30', 'seep.cond=5000.0']
    back = cfgmod.load_run_config(out)
    assert back.run.nsp == 30 and back.seep.cond == 5000.0
    assert back.config_hash() == digest


def test_save_refuses_a_folder_that_does_not_exist(cfg, tmp_path):
    with pytest.raises(editor.EditError):
        editor.save(cfg, str(tmp_path / 'nope' / 'c.toml'))


def test_nothing_is_written_when_the_edit_is_invalid(tmp_path, cfg):
    out = str(tmp_path / 'c.toml')
    with pytest.raises(cfgmod.ConfigError):
        editor.save(cfg, out, {'grid.cell_size': 0.0})
    assert not os.path.exists(out), 'an invalid configuration reached the disk'


# ----------------------------------------------------------- table editing

def test_table_rows_round_trip(cfg):
    rows = editor.table_rows(cfg, 'surface.vegetation')
    assert [r['name'] for r in rows] == ['grassMU', 'Qilex', 'Qpyr']
    new = editor.set_table(cfg, 'surface.vegetation', rows)
    assert new.config_hash() == cfg.config_hash()


def test_editing_a_row_takes_effect(cfg):
    rows = editor.table_rows(cfg, 'surface.vegetation')
    rows[1]['root_depth'] = 12.0
    new = editor.set_table(cfg, 'surface.vegetation', rows)
    assert new.surface.vegetation[1].root_depth == 12.0
    assert cfg.surface.vegetation[1].root_depth == 15.0      # original intact


def test_a_row_missing_a_column_is_refused(cfg):
    """A data editor that drops a column would otherwise reset it to the
    dataclass default -- for a vegetation type, silently to FAO-56 grass."""
    rows = editor.table_rows(cfg, 'surface.vegetation')
    del rows[0]['root_depth']
    with pytest.raises(editor.EditError) as e:
        editor.set_table(cfg, 'surface.vegetation', rows)
    assert 'root_depth' in str(e.value)


def test_an_unknown_column_is_refused(cfg):
    rows = editor.table_rows(cfg, 'surface.vegetation')
    rows[0]['bankfull'] = 1.0
    with pytest.raises(editor.EditError) as e:
        editor.set_table(cfg, 'surface.vegetation', rows)
    assert 'bankfull' in str(e.value)


def test_a_table_edit_is_validated(cfg):
    rows = editor.table_rows(cfg, 'surface.vegetation')
    rows[0]['kt_s'] = 21.98151201            # the old ini's 1/s
    with pytest.raises(cfgmod.ConfigError) as e:
        editor.set_table(cfg, 'surface.vegetation', rows)
    assert 'not the 1/s' in str(e.value)


def test_add_and_drop_a_row(cfg):
    n = len(cfg.surface.station)
    more = editor.add_row(cfg, 'surface.station')
    assert len(more.surface.station) == n + 1
    back = editor.drop_row(more, 'surface.station', n)
    assert len(back.surface.station) == n


def test_dropping_a_row_that_is_not_there_is_an_error(cfg):
    with pytest.raises(editor.EditError):
        editor.drop_row(cfg, 'surface.station', 99)


def test_row_defaults_carry_the_dataclass_defaults(cfg):
    row = editor.row_defaults(cfg, 'surface.vegetation')
    assert row['name'] == 'veg1' and row['kt_s'] == 1.0
    assert set(row) == {f for f in
                        editor.table_rows(cfg, 'surface.vegetation')[0]}


def test_a_table_path_that_is_not_one_is_refused(cfg):
    with pytest.raises(editor.EditError):
        editor.table_rows(cfg, 'surface.meteo_ts')
    with pytest.raises(editor.EditError):
        editor.table_rows(cfg, 'surface.nope')
