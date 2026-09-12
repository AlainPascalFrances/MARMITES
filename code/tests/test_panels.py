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


# ------------------------------------------------- panel 1: the rectangle

meshes = _load('marmites_meshes_p', os.path.join(CODE, 'marmites_meshes.py'))


def test_the_rectangle_is_derived_from_the_polygon(cfg):
    """Panel 1's whole point: the domain comes first and the grid is built
    inside it."""
    cfg.grid.kind = 'structured'                          # 50 m cells
    bbox = (739293.0, 4553110.0, 742223.0, 4556240.0)     # lm_lim.shp
    nrow, ncol, delr, delc, xll, yll = meshes.model_rectangle(cfg, bbox)
    assert (xll, yll) == (739250.0, 4553100.0)            # snapped DOWN
    assert nrow == 63 and ncol == 60
    assert xll + ncol * cfg.grid.cell_size >= bbox[2]     # covers the polygon
    assert yll + nrow * cfg.grid.cell_size >= bbox[3]
    assert len(delr) == ncol and len(delc) == nrow


def test_the_origin_is_snapped_so_the_grid_is_reproducible(cfg):
    """A rectangle that shifted with a re-exported shapefile would invalidate
    every cached mesh for no reason."""
    a = meshes.model_rectangle(cfg, (739293.0, 4553110.0, 742223.0, 4556240.0))
    b = meshes.model_rectangle(cfg, (739299.9, 4553149.9, 742223.0, 4556240.0))
    assert a[4:] == b[4:]


def test_the_buffer_grows_the_rectangle(cfg):
    bbox = (1000.0, 2000.0, 1500.0, 2500.0)
    base = meshes.model_rectangle(cfg, bbox)
    cfg.grid.buffer = 200.0
    wide = meshes.model_rectangle(cfg, bbox)
    assert wide[1] > base[1] and wide[0] > base[0]
    assert wide[4] < base[4] and wide[5] < base[5]


def test_the_override_reproduces_an_existing_grid(cfg):
    cfg.grid.override.enable = True
    cfg.grid.override.nrow, cfg.grid.override.ncol = 65, 60
    cfg.grid.override.xllcorner = 739300.0
    cfg.grid.override.yllcorner = 4553050.0
    nrow, ncol, _dr, _dc, xll, yll = meshes.model_rectangle(cfg)
    assert (nrow, ncol) == (65, 60)
    assert (xll, yll) == (739300.0, 4553050.0)


def test_no_polygon_and_no_override_is_an_error(cfg):
    with pytest.raises(meshes.MeshBuildError) as e:
        meshes.model_rectangle(cfg, None)
    assert 'grid.boundary' in str(e.value)


def test_a_zero_cell_size_is_an_error(cfg):
    cfg.grid.kind = 'structured'
    cfg.grid.cell_size = 0.0
    with pytest.raises(meshes.MeshBuildError) as e:
        meshes.model_rectangle(cfg, (0.0, 0.0, 100.0, 100.0))
    assert 'grid.cell_size' in str(e.value)


def test_a_voronoi_rectangle_is_snapped_to_cell_far(cfg):
    """WP1d panel 1, D1 of 2A.7: the Voronoi producer never uses cell_size --
    cell_far is the size -- so snapping the rectangle to cell_size left one
    number doing nothing and another doing the work."""
    cfg.grid.kind = 'voronoi'
    cfg.grid.cell_size = 50.0
    cfg.grid.voronoi.cell_far = 100.0
    bbox = (739293.0, 4553110.0, 742223.0, 4556240.0)
    _nr, _nc, delr, _dc, xll, yll = meshes.model_rectangle(cfg, bbox)
    assert meshes.rectangle_cell_size(cfg) == 100.0
    assert (xll, yll) == (739200.0, 4553100.0)            # 100 m, not 50 m
    assert float(delr[0]) == 100.0
    cfg.grid.cell_size = 25.0                             # and it does not move
    assert meshes.model_rectangle(cfg, bbox)[4:] == (xll, yll)
    cfg.grid.voronoi.cell_far = 0.0
    with pytest.raises(meshes.MeshBuildError) as e:
        meshes.model_rectangle(cfg, bbox)
    assert 'cell_far' in str(e.value)


def test_grid_stub_carries_what_the_producers_read(cfg):
    """They use nrow, ncol, nlay, delr, delc and the origin -- and nothing
    else, which is what lets a grid be previewed without the MF ini."""
    stub = meshes.grid_stub(cfg, (739293.0, 4553110.0, 742223.0, 4556240.0),
                            nlay=2)
    for attr in ('nrow', 'ncol', 'nlay', 'delr', 'delc', 'xllcorner',
                 'yllcorner'):
        assert hasattr(stub, attr), attr
    assert stub.nlay == 2


# ------------------------------------------------------------- choices

def test_the_enumerated_fields_are_choices():
    kinds = schema.choices_for('grid.kind')
    assert set(kinds) >= {'structured', 'disv', 'voronoi', 'quadtree'}
    assert schema.choices_for('seep.kind') == ['uzf', 'drn']
    assert schema.choices_for('postproc.wb_unit') == ['year', 'day']
    assert schema.choices_for('grid.cell_size') is None


def test_the_choices_come_from_the_schema_not_a_copy():
    """A grid kind added to marmites_config must appear in the panel by
    itself, or the two drift apart."""
    assert schema.choices_for('grid.kind') == list(cfgmod.GRID_KINDS)
    assert schema.choices_for('grid.resample') == list(cfgmod.RESAMPLE_MODES)


def test_the_conditional_blocks_know_what_controls_them():
    assert schema.subpanel_for('grid.voronoi.cell_far') == \
        ('grid.kind', ('voronoi',))
    assert schema.subpanel_for('grid.quadtree.refine_level') == \
        ('grid.kind', ('quadtree',))
    # WP1d: these two used to be permanent, and are not. cell_size is the
    # cell for three kinds and nothing at all for voronoi; the override
    # reproduces an EXISTING rectangle, which only the legacy kinds have.
    assert schema.subpanel_for('grid.cell_size') == \
        ('grid.kind', ('structured', 'disv', 'quadtree'))
    assert schema.subpanel_for('grid.override.nrow') == \
        ('grid.kind', ('structured', 'disv'))


def test_the_conditional_blocks_are_separable_from_the_rest(cfg):
    """``fields_of`` is the full inventory -- that is what the label check
    needs. The PANEL splits it, and this is the rule it splits on."""
    keys = [k for k, _v in schema.fields_of(cfg, 'grid')]
    conditional = [k for k in keys if schema.subpanel_for(k)]
    always = [k for k in keys if not schema.subpanel_for(k)]
    assert all(k.startswith(('grid.voronoi.', 'grid.quadtree.',
                             'grid.override.', 'grid.cell_size'))
               for k in conditional)
    assert 'grid.boundary' in always
    assert 'grid.crs_epsg' in always
    assert 'grid.cell_size' not in always


def test_every_grid_field_is_either_permanent_or_on_one_kinds_subpanel(cfg):
    """The two lists PARTITION the [grid] block: a field in neither is a
    field no panel shows, which is how a live setting goes invisible."""
    keys = {k for k, _v in schema.fields_of(cfg, 'grid')}
    shown = set(schema.GRID_PERMANENT)
    for kind in schema.GRID_SUBPANEL:
        shown.update(schema.grid_fields(kind))
    assert not (keys - shown), 'shown by no panel: %s' % sorted(keys - shown)
    assert not (shown - keys), 'not in the configuration: %s' % sorted(shown - keys)
    assert not (set(schema.GRID_PERMANENT)
                & set(schema.grid_fields('voronoi')))


def test_every_grid_subpanel_field_matches_its_visibility_rule(cfg):
    """A field on the voronoi sub-panel whose rule says 'structured' would be
    drawn and then refused by validate()."""
    for kind in schema.GRID_SUBPANEL:
        for dotted in schema.grid_fields(kind):
            rule = schema.subpanel_for(dotted)
            if rule is not None:
                assert kind in rule[1], \
                    '%s is on the %s sub-panel but applies to %s' % (
                        dotted, kind, ', '.join(rule[1]))


def test_a_gating_switch_is_on_its_own_row_above_what_it_controls(cfg):
    """The nesting IS the layout: a switch shares a row with nothing, and
    everything it gates comes after it."""
    for kind, rows in schema.GRID_SUBPANEL.items():
        flat = schema.grid_fields(kind)
        for switch, dependents in schema.GRID_GATED.items():
            if switch not in flat:
                continue
            row = [r for r in rows if switch in r][0]
            assert len(row) == 1, \
                '%s shares its row with %s' % (switch, row)
            for dep in dependents:
                if dep in flat:
                    assert flat.index(dep) > flat.index(switch), \
                        '%s is drawn before the switch that gates it' % dep


# ------------------------------------------------ panel 1: the mesh viewer

def test_a_mesh_cache_folder_round_trips(tmp_path):
    """The viewer reads the layout the producer writes, wherever it sits --
    panel 1's experiments each get a folder of their own, so the reader
    cannot assume the run's <ws>/MF6_ws_<kind>/_mesh."""
    import json
    from lib import loaders
    d = tmp_path / 'attempt_1'
    d.mkdir()
    (d / 'mesh_voronoi.json').write_text(
        json.dumps({'vertices': [], 'cell2d': [], 'ncpl': 7}), encoding='utf-8')
    (d / 'mesh_voronoi.sig.json').write_text(
        json.dumps({'signature': 'abc123', 'kind': 'voronoi', 'ncpl': 7}),
        encoding='utf-8')
    gp, sig = loaders.read_mesh_at(str(d), 'voronoi')
    assert gp['ncpl'] == 7 and sig['signature'] == 'abc123'
    assert loaders.read_mesh_at(str(tmp_path / 'nothing'), 'voronoi') == (None,
                                                                          None)


def test_promoting_a_mesh_puts_it_where_a_run_looks(tmp_path):
    """Selecting a grid must not leave the run to rebuild it -- and the
    SIGNATURE has to travel with it, or the promoted mesh would be served for
    settings it does not belong to."""
    import json
    from lib import loaders
    src = tmp_path / 'attempt_2'
    src.mkdir()
    (src / 'mesh_voronoi.json').write_text(
        json.dumps({'vertices': [], 'cell2d': [], 'ncpl': 830}),
        encoding='utf-8')
    (src / 'mesh_voronoi.sig.json').write_text(
        json.dumps({'signature': 'deadbeef', 'ncpl': 830}), encoding='utf-8')
    ws = tmp_path / 'ws'
    dst, moved = loaders.promote_mesh(str(src), str(ws), 'voronoi')
    assert len(moved) == 2
    gp, sig = loaders.read_mesh(str(ws), 'voronoi')
    assert gp['ncpl'] == 830 and sig['signature'] == 'deadbeef'
    assert dst == loaders.mesh_cache_paths(str(ws), 'voronoi')[0].parent


def test_the_signature_separates_two_corridors(cfg):
    """The viewer decides whether the mesh on screen is the one the settings
    describe by comparing signatures. Cell counts would not do: two meshes
    can share one and differ."""
    cfg.grid.kind = 'voronoi'
    stub = meshes.grid_stub(cfg, (0.0, 0.0, 3000.0, 3000.0), nlay=1)
    a = meshes.mesh_signature(cfg, stub)
    # Perturbed from whatever the shipped configuration holds, not set to
    # fixed numbers: this test once asserted the values the file already had.
    # And perturbed on the INPUTS -- the corridor is derived from them, so
    # assigning it and validating would simply put it back.
    cfg.grid.voronoi.cell_near_stream = cfg.grid.voronoi.cell_near_stream / 2.0
    cfg.grid.voronoi.refresh()
    b = meshes.mesh_signature(cfg, stub)
    assert b != a, ('a different corridor gives the same signature, so a '
                    'stale mesh is served for it')
    cfg.grid.voronoi.grade_ratio = 1.9
    cfg.grid.voronoi.refresh()
    assert meshes.mesh_signature(cfg, stub) != b


def test_the_quadtree_refinement_is_in_the_signature(cfg):
    """It was not until WP1d: changing refine_level built a different mesh
    and the cache served the old one."""
    cfg.grid.kind = 'quadtree'
    stub = meshes.grid_stub(cfg, (0.0, 0.0, 3000.0, 3000.0), nlay=1)
    a = meshes.mesh_signature(cfg, stub)
    cfg.grid.quadtree.refine_level = 4
    assert meshes.mesh_signature(cfg, stub) != a
