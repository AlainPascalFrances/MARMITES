# -*- coding: utf-8 -*-
"""WP1d -- the four panels: what they show, and editing what they show.

The pages are thin views over ``app/lib/schema.py`` and ``app/lib/editor.py``,
so what could actually be wrong is tested here rather than found in a browser:
a configuration field no panel shows, a panel naming a section that does not
exist, an edit that would write an invalid file, and a table row that silently
loses a column.
"""

import importlib.util
import io
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
REPO = os.path.abspath(os.path.join(CODE, '..'))
REF = os.path.join(CODE, 'configs', 'lamata.toml')

for _p in (CODE, os.path.join(CODE, 'app')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


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
    assert nums == [0, 1, 2, 3, 4, 5, 6]
    titles = [p[1] for p in schema.PANELS]
    assert titles == ['Overview', 'Grid', 'Surface and driving forces',
                      'Soil', 'Unsaturated zone and groundwater',
                      'State variables (calibration)', 'Plots']


def test_the_master_switches_are_real_run_keys(cfg):
    """A switch that is not a key the driver reads would be decoration."""
    for _n, _t, _i, switch, _s, _b in schema.PANELS:
        if switch is None:
            continue
        section, key = switch.split('.')
        assert hasattr(getattr(cfg, section), key), switch
        assert isinstance(getattr(getattr(cfg, section), key), bool)


def test_the_switches_are_where_the_run_is_decided():
    by = {p[1]: p[3] for p in schema.PANELS}
    assert by['Overview'] is None and by['Grid'] is None
    assert by['Surface and driving forces'] == 'run.surface'
    assert by['Plots'] == 'run.plot'
    # Nothing on the calibration panel changes a flux, so it has no switch.
    assert by['State variables (calibration)'] is None
    # ONE switch on TWO panels: MMsoil is stepped from inside the MODFLOW
    # time loop, so there is no running one without the other, and a panel
    # that could turn off half of it would be lying.
    assert by['Soil'] == 'run.model'
    assert by['Unsaturated zone and groundwater'] == 'run.model'


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
    assert schema.panel_of('soil') == 3
    assert schema.panel_of('sfr') == 4
    assert schema.panel_of('obs') == 5
    assert schema.panel_of('postproc') == 6
    assert schema.panel_of('nowhere') is None


def test_the_state_variables_are_the_four_the_model_produces():
    """Groundwater levels, soil moisture, actual evapotranspiration and
    surface runoff. The fourth was not asked for at all before the split,
    although the model computes it."""
    titles = [g[0] for g in schema.OBS_GROUPS]
    assert titles == ['Groundwater levels', 'Soil moisture',
                      'Actual evapotranspiration', 'Surface runoff']


def test_every_state_variable_names_a_field_that_exists(cfg):
    have = dict(schema.fields_of(cfg, 'obs'))
    for _title, dotted, _note in schema.OBS_GROUPS:
        assert dotted in have, '%s names a missing field' % dotted
    for dotted in schema.OBS_COMMON:
        assert dotted in have, '%s names a missing field' % dotted


def test_the_observation_fields_are_all_on_the_panel(cfg):
    """A prefix shown nowhere is a series nobody can point at."""
    laid = set(schema.OBS_COMMON) | {g[1] for g in schema.OBS_GROUPS}
    have = set(dict(schema.fields_of(cfg, 'obs')))
    assert have == laid, 'unplaced: %s' % sorted(have - laid)


def test_actual_evapotranspiration_starts_empty(cfg):
    """La Mata has none, and a blank prefix is an ANSWER -- not measured
    here -- rather than an omission."""
    assert cfg.obs.aet_prefix == ''


def _spoken(path):
    """Every string literal in a file that is NOT a docstring.

    A docstring is written for whoever reads the code and may say whatever is
    clearest; these are the strings that reach the screen. Parsed rather than
    grepped, because a comment and a docstring both look like text to a
    regular expression.
    """
    import ast

    tree = ast.parse(io.open(path, encoding='utf-8').read())
    docs = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            body = getattr(node, 'body', None)
            if body and isinstance(body[0], ast.Expr) \
                    and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                docs.add(id(body[0].value))
    return [(n.lineno, n.value) for n in ast.walk(tree)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)
            and id(n) not in docs]


def test_a_panel_is_named_not_numbered():
    """The modeller sees NAMES in the sidebar, and the numbers move whenever
    a panel is split -- three of them moved when Model became Soil,
    Unsaturated zone and State variables. A number in a sentence is a number
    that will one day point at the wrong panel."""
    import re

    bad = []
    for folder, _dirs, files in os.walk(os.path.join(CODE, 'app')):
        if '__pycache__' in folder:
            continue
        for name in sorted(f for f in files if f.endswith('.py')):
            for line, text in _spoken(os.path.join(folder, name)):
                if re.search(r'\bpanels? [0-9]', text) \
                        and 'panel %d' not in text:
                    bad.append('%s:%d  %s' % (name, line, text[:60]))
    assert not bad, 'a panel is named by its number:\n  ' + '\n  '.join(bad)


def test_panel_name_comes_from_the_panels_themselves():
    assert schema.panel_name(1) == 'Grid'
    assert schema.panel_name(2) == 'Surface and driving forces'
    assert schema.panel_name(5) == 'State variables (calibration)'
    for num, title, _i, _s, _sec, _b in schema.PANELS:
        assert schema.panel_name(num) == title
    # a number that is not a panel says so rather than pretending
    assert 'panel' in schema.panel_name(99)


def test_the_model_path_names_panels_that_exist():
    """marmites_config cannot import the schema -- it is the model path, and
    the panels are the app's -- so its messages carry the names literally.
    That is only safe while the names are real ones."""
    import re

    titles = {p[1] for p in schema.PANELS}
    for line, text in _spoken(os.path.join(CODE, 'marmites_config.py')):
        assert not re.search(r'\bpanels? [0-9]', text), \
            'marmites_config:%d names a panel by number: %s' % (line,
                                                                text[:60])
        for word in re.findall(r'the ([A-Z][A-Za-z ]+?) panel', text):
            assert word in titles, \
                'marmites_config:%d: no panel is called %r' % (line, word)


def test_the_time_discretisation_is_asked_where_the_record_is(cfg):
    """The days are MMsurf's -- it writes one row per day -- so what the run
    does with them belongs beside the record, not in the MODFLOW parameter
    file where the aggregation limit sat under the name `nper`."""
    have = dict(schema.fields_of(cfg, 'run'))
    for row in schema.TIME_ROWS:
        for dotted in row:
            if dotted is None:
                continue
            assert dotted in have, '%s is laid out but does not exist' % dotted
    laid = [d for row in schema.TIME_ROWS for d in row if d]
    assert laid == ['run.daily', 'run.perlen_max', 'run.nsp']


def test_the_aggregation_limit_is_a_length_not_a_count(cfg):
    """`nper` in the MODFLOW ini is NOT the number of periods: it is the
    longest one, and the count is worked out from the rainfall."""
    assert cfg.run.perlen_max >= 1
    cfg.run.perlen_max = 0
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'perlen_max' in str(e.value)


def test_the_record_length_is_read_from_the_dates(tmp_path):
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        'msurf_t', os.path.join(CODE, 'marmites_surface.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['msurf_t'] = mod
    spec.loader.exec_module(mod)

    assert mod.record_days(str(tmp_path)) == 0      # before MMsurf has run
    p = os.path.join(str(tmp_path), mod.FORCING['date'])
    io.open(p, 'w', encoding='utf-8', newline='').write(
        '\n'.join(['# a comment', '01/01/2000 1 1', '02/01/2000 2 2', '',
                   '03/01/2000 3 3', '']))
    assert mod.record_days(str(tmp_path)) == 3


# --------------------------------------------- panel 4: the MF6 packages

def test_the_geometry_asks_only_what_is_not_derivable(cfg):
    """The top is the land surface minus the soil column, so asking for it
    would be asking for an answer that can only disagree with one already
    given.

    IBOUND IS ASKED, and that is a correction: the catchment polygon gives
    the outline, but a layer can pinch out inside it. Measured on La Mata,
    layer 1 is active in 1870 cells and layer 2 in 1954 -- and thick_l1 is
    20-35 m in the 84 cells where layer 1 is absent, so the thickness does
    not express it either. The polygon is the geographic REFERENCE instead:
    every input is checked for sitting inside it."""
    laid = [d for row in schema.GEOMETRY_ROWS for d in row if d]
    assert laid == ['layers.nlay', 'layers.hnoflo',
                    'layers.ibound', 'layers.thickness',
                    'layers.k', 'layers.k33', 'layers.k33_as_ratio',
                    'layers.ss', 'layers.sy', 'layers.convertible']
    # the switch that says what the k33 NUMBER means must sit with it: 2 is
    # an anisotropy ratio or a conductivity depending on that one box
    flat = [d for row in schema.GEOMETRY_ROWS for d in row]
    assert abs(flat.index('layers.k33_as_ratio')
               - flat.index('layers.k33')) <= 2
    have = dict(schema.fields_of(cfg, 'layers'))
    for dotted in laid:
        assert dotted in have, '%s is laid out but does not exist' % dotted
    assert not any('top' in d for d in laid)
    # ... and ibound is a SOURCE, so it can be a raster per layer
    assert schema.is_source(cfg.layers.ibound)


def test_a_per_layer_raster_expands_over_the_layers():
    """`k_%d.asc` with 2 layers is k_1.asc and k_2.asc -- %d is the WHOLE
    placeholder, so no stray `d` and no leading zero."""
    src = cfgmod.VectorSource(raster='k_%d.asc')
    assert src.rasters(2) == ['k_1.asc', 'k_2.asc']
    assert src.rasters(1) == ['k_1.asc']
    # a name without the placeholder is that one raster, whatever nlay says
    assert cfgmod.VectorSource(raster='k.asc').rasters(6) == ['k.asc']
    # and nothing to expand when there is no raster at all
    assert cfgmod.VectorSource(value=1.0).rasters(2) == []


def test_the_layers_help_says_what_the_placeholder_becomes():
    """A pattern is worth nothing if the modeller has to guess whether it
    means k_1.asc or k_01.asc."""
    h = schema.describe('layers.thickness')[2]
    assert 'thick_1.asc' in h and 'thick_2.asc' in h


def test_every_geometry_field_names_its_flopy_argument():
    """The flopy name is the only one that does not drift."""
    for dotted in ('layers.nlay', 'layers.thickness', 'layers.k',
                   'layers.ss', 'layers.sy'):
        help_ = schema.describe(dotted)[2]
        assert 'flopy' in help_.lower(), '%s does not say what it becomes'


def test_hnoflo_may_not_be_zero(cfg):
    """It marks a cell as having nothing to report, and 0 is a perfectly
    good head."""
    cfg.layers.hnoflo = 0.0
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'hnoflo' in str(e.value)


def test_the_layer_properties_take_the_usual_three_producers(cfg):
    """Raster, layer or one value -- the rule every spatial input follows."""
    for name in ('thickness', 'k', 'ss', 'sy'):
        src = getattr(cfg.layers, name)
        assert schema.is_source(src), '%s is not a source' % name
        for producer in ('raster', 'layer', 'value'):
            assert hasattr(src, producer), '%s has no %s' % (name, producer)
        # whatever is answered, exactly one producer wins -- and 'value'
        # is a number rather than the string a hand-edited file may hold
        if src.producer() is not None:
            src.validate('layers.%s' % name)


def test_hnoflo_reaches_the_model():
    """A panel field that the run ignores is decoration. This one is read
    straight after the ini is parsed, before anything uses it."""
    src = io.open(os.path.join(HERE, 'run_lamata_mf6.py'),
                  encoding='utf-8').read()
    assert 'cfg.layers.hnoflo' in src, 'the run never reads layers.hnoflo'
    assert src.index('cfg.layers.hnoflo') < src.index('conv_fact = '), (
        'hnoflo is applied after something has already used it')


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
    # Pinned, not inherited: this reads the live configuration, and a
    # cell size or buffer changed in the front-end would otherwise move the
    # origin this test asserts.
    cfg.grid.kind = 'structured'
    cfg.grid.cell_size, cfg.grid.buffer = 50.0, 0.0
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
    cfg.grid.kind = 'structured'
    cfg.grid.cell_size, cfg.grid.buffer = 50.0, 0.0
    a = meshes.model_rectangle(cfg, (739293.0, 4553110.0, 742223.0, 4556240.0))
    b = meshes.model_rectangle(cfg, (739299.9, 4553149.9, 742223.0, 4556240.0))
    assert a[4:] == b[4:]


def test_the_buffer_grows_the_rectangle(cfg):
    cfg.grid.buffer = 0.0
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
    cfg.grid.cell_size, cfg.grid.buffer = 50.0, 0.0
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


# --------------------------------- panel 1: standing on the model's rasters

def _asc(path, xll, yll, nrow, ncol, cs, value=1.0):
    """A minimal ESRI ASCII grid, written where a dataset would hold one."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with io.open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('ncols %d\nnrows %d\nxllcorner %.10g\nyllcorner %.10g\n'
                 'cellsize %.10g\nnodata_value -9999\n' % (ncol, nrow, xll,
                                                           yll, cs))
        for _ in range(nrow):
            fh.write(' '.join(['%g' % value] * ncol) + '\n')


def _legacy_dataset(root, xll=739300.0, yll=4553050.0, nrow=65, ncol=60,
                    cs=50.0):
    """A dataset shaped like La Mata's: tables at the top, MF rasters below."""
    _asc(os.path.join(root, 'inputSOILzones.asc'), xll, yll, nrow, ncol, cs)
    _asc(os.path.join(root, 'inputVEG1area.asc'), xll, yll, nrow, ncol, cs)
    _asc(os.path.join(root, 'MF_ws', 'elev.asc'), xll, yll, nrow, ncol, cs)
    _asc(os.path.join(root, 'MF_ws', 'ibound_l1.asc'), xll, yll, nrow, ncol, cs)
    return root


def _pin(cfg, kind='structured', cs=50.0):
    cfg.grid.kind = kind
    cfg.grid.cell_size, cfg.grid.buffer = cs, 0.0
    cfg.grid.override.enable = False
    return cfg


BBOX = (739293.0, 4553110.0, 742223.0, 4556240.0)      # lm_lim.shp


def test_the_dataset_rectangle_is_read_off_the_rasters(cfg, tmp_path):
    root = _legacy_dataset(str(tmp_path))
    rect, names, others = meshes.dataset_rectangle(root)
    assert rect == (739300.0, 4553050.0, 65, 60, 50.0)
    assert len(names) == 4 and not others          # the MF_ws ones are included


def test_the_dem_is_not_expected_on_the_model_rectangle(cfg, tmp_path):
    """It is kept at its own resolution and wrapped onto the cells at run
    time, so counting it would report a mismatch that is by design."""
    root = _legacy_dataset(str(tmp_path))
    _asc(os.path.join(root, 'inputDEM.asc'), 739166.1, 4552982.6, 20, 20, 5.0)
    rect, names, others = meshes.dataset_rectangle(root)
    assert rect == (739300.0, 4553050.0, 65, 60, 50.0)
    assert not others and not any('DEM' in n for n in names)


def test_rasters_on_a_third_rectangle_are_reported_not_hidden(cfg, tmp_path):
    """A majority vote that swallowed the odd ones out would be worse than
    saying it: the model reads them as it reads the rest."""
    root = _legacy_dataset(str(tmp_path))
    _asc(os.path.join(root, 'MF_ws', 'vka_l1_old.asc'),
         739325.0, 4553190.0, 69, 72, 40.0)
    _rect, _names, others = meshes.dataset_rectangle(root)
    assert len(others) == 1
    assert others[0][0] == (739325.0, 4553190.0, 69, 72, 40.0)
    assert others[0][1] == ['MF_ws%svka_l1_old.asc' % os.sep]


def test_the_derived_rectangle_overhangs_the_legacy_rasters(cfg, tmp_path):
    """The real La Mata case: the polygon reaches 739293.4 and the snap goes
    DOWN to 739250, while the rasters were frozen at 739300."""
    root = _legacy_dataset(str(tmp_path))
    rc = meshes.rectangle_check(_pin(cfg), BBOX, root)
    assert rc['status'] == 'overhang'
    assert rc['overhang'] == {'west': 50.0}
    assert rc['derived'][0] == 739250.0 and rc['source'][0] == 739300.0
    assert 'west' in rc['detail']


def test_the_overhang_does_not_depend_on_the_grid_kind(cfg, tmp_path):
    """It is a property of the RECTANGLE, so a mesh is in exactly the same
    position as the structured grid -- it simply has no override to escape
    with."""
    root = _legacy_dataset(str(tmp_path))
    cfg.grid.voronoi.cell_far = 50.0
    for kind in ('structured', 'disv', 'voronoi', 'quadtree'):
        rc = meshes.rectangle_check(_pin(cfg, kind), BBOX, root)
        assert rc['status'] == 'overhang', kind
        assert rc['overhang'] == {'west': 50.0}, kind


def test_the_override_puts_the_grid_back_on_the_rasters(cfg, tmp_path):
    """Which is what the panel offers, and the only lever there is until the
    model panel re-derives its own rasters."""
    root = _legacy_dataset(str(tmp_path))
    _pin(cfg)
    cfg.grid.override.enable = True
    cfg.grid.override.nrow, cfg.grid.override.ncol = 65, 60
    cfg.grid.override.xllcorner = 739300.0
    cfg.grid.override.yllcorner = 4553050.0
    rc = meshes.rectangle_check(cfg, BBOX, root)
    assert rc['status'] == 'ok'
    assert rc['derived'] == rc['source']


def test_a_grid_off_the_rasters_lattice_is_flagged(cfg, tmp_path):
    """Covered, so the projection runs -- but every cell is then resampled
    from fractions of four, which blurs the legacy model."""
    root = _legacy_dataset(str(tmp_path))
    _pin(cfg)
    cfg.grid.override.enable = True
    cfg.grid.override.nrow, cfg.grid.override.ncol = 10, 10
    cfg.grid.override.xllcorner = 739325.0            # half a cell east
    cfg.grid.override.yllcorner = 4553050.0
    rc = meshes.rectangle_check(cfg, BBOX, root)
    assert rc['status'] == 'shifted'
    assert '25' in rc['detail']


def test_an_empty_dataset_has_nothing_to_disagree_with(cfg, tmp_path):
    """A brand new catchment: from-scratch must not be told it is wrong."""
    rc = meshes.rectangle_check(_pin(cfg), BBOX, str(tmp_path))
    assert rc['status'] == 'none'
    assert rc['source'] is None and rc['derived'] is not None


def test_a_rectangle_that_cannot_be_derived_reports_error(cfg, tmp_path):
    """The page draws this before the build button, so it must never be the
    thing that breaks the page."""
    root = _legacy_dataset(str(tmp_path))
    rc = meshes.rectangle_check(_pin(cfg), None, root)
    assert rc['status'] == 'error' and 'grid.boundary' in rc['detail']


def test_the_header_is_read_without_the_body(tmp_path):
    dem = _load('marmites_dem_p', os.path.join(CODE, 'marmites_dem.py'))
    path = os.path.join(str(tmp_path), 'x.asc')
    _asc(path, 1000.0, 2000.0, 3, 4, 25.0)
    head = dem.read_asc_header(path)
    assert (head['xllcorner'], head['cellsize']) == (1000.0, 25.0)
    assert int(head['nrows']) == 3 and int(head['ncols']) == 4
    arr, head2 = dem.read_asc(path)
    assert head2 == head and arr.shape == (3, 4)


# ------------------------------------------ panel 1: the ponds in the mesh

def _pond_geojson(path, polys):
    """``inputPONDS.geojson`` as the converter writes it."""
    import json
    feats = [{'type': 'Feature', 'properties': {'id': i + 1},
              'geometry': {'type': 'Polygon',
                           'coordinates': [[list(p) for p in ring]
                                           + [list(ring[0])]]}}
             for i, ring in enumerate(polys)]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    io.open(path, 'w', encoding='utf-8', newline='').write(
        json.dumps({'type': 'FeatureCollection', 'features': feats}))


def _square(cx, cy, half):
    return [(cx - half, cy - half), (cx + half, cy - half),
            (cx + half, cy + half), (cx - half, cy + half)]


def test_the_pond_footprints_are_read_from_the_geojson(tmp_path):
    """inputPONDS.csv carries the centroid and the area, which places a cell
    but cannot draw a footprint or size one -- the outlines are in the
    GeoJSON, and the mesh producer reads it with json alone."""
    root = str(tmp_path)
    _pond_geojson(os.path.join(root, 'inputPONDS.geojson'),
                  [_square(100.0, 200.0, 10.0), _square(500.0, 600.0, 5.0)])
    rings = meshes.pond_rings(root)
    assert len(rings) == 2
    assert len(rings[0]) == 4, 'the closing point should be dropped'
    area, cx, cy = meshes._ring_area_centre(rings[0])
    assert abs(area - 400.0) < 1e-6
    assert (round(cx, 6), round(cy, 6)) == (100.0, 200.0)


def test_a_catchment_with_no_ponds_is_not_an_error(tmp_path):
    assert meshes.pond_rings(str(tmp_path)) == []
    nodes, ponds, gone = meshes.pond_seeds(str(tmp_path))
    assert (nodes, ponds, gone) == ([], [], [])


def test_each_pond_is_seeded_with_a_centre_and_a_ring(tmp_path):
    """The centre gives the cell its position; the ring keeps the cell from
    simply taking the background size (CdL, 2026-07-04)."""
    root = str(tmp_path)
    _pond_geojson(os.path.join(root, 'inputPONDS.geojson'),
                  [_square(100.0, 200.0, 10.0)])
    nodes, ponds, _gone = meshes.pond_seeds(root)
    assert len(ponds) == 1
    assert nodes[0] == (100.0, 200.0), 'the first node is the pond centre'
    assert len(nodes) == 1 + meshes.POND_RING_N
    r = ponds[0][4]
    for x, y in nodes[1:]:
        d = ((x - 100.0) ** 2 + (y - 200.0) ** 2) ** 0.5
        assert abs(d - meshes.POND_RING_FACTOR * r) < 1e-6


def test_a_pond_outside_the_domain_is_refused_and_reported(tmp_path):
    """Triangle discards a node outside the boundary polygon, so without this
    the pond vanished silently -- no cell, no message, and the LAK footprint
    later looking for one."""
    root = str(tmp_path)
    _pond_geojson(os.path.join(root, 'inputPONDS.geojson'),
                  [_square(100.0, 200.0, 10.0), _square(9000.0, 9000.0, 10.0)])
    inside = meshes._inside_ring([(0.0, 0.0), (1000.0, 0.0),
                                  (1000.0, 1000.0), (0.0, 1000.0)])
    nodes, ponds, gone = meshes.pond_seeds(root, inside=inside)
    assert len(ponds) == 1 and len(gone) == 1
    assert all(0.0 <= x <= 1000.0 and 0.0 <= y <= 1000.0 for x, y in nodes)


def test_the_point_in_polygon_test_is_the_crossing_number_rule():
    inside = meshes._inside_ring([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0),
                                  (0.0, 10.0)])
    assert inside(5.0, 5.0)
    assert not inside(15.0, 5.0)
    assert not inside(5.0, -1.0)


def test_a_concave_pond_is_seeded_at_its_own_centroid(tmp_path):
    """The centroid of the POLYGON, not the mean of its vertices: a rim
    mapped with many points down one side would drag the mean that way and
    the cell would not be centred on the water."""
    ring = [(0.0, 0.0), (100.0, 0.0), (100.0, 10.0), (10.0, 10.0),
            (10.0, 100.0), (0.0, 100.0)]
    area, cx, cy = meshes._ring_area_centre(ring)
    assert abs(area - 1900.0) < 1e-6
    mean_x = sum(p[0] for p in ring) / len(ring)
    assert abs(cx - cy) < 1e-9, 'the L is symmetric about the diagonal'
    assert abs(cx - mean_x) > 1.0, 'the vertex mean is NOT the centroid'


def test_a_pond_zone_stands_off_the_pond_it_sizes(cfg):
    """A vertex of a constraint polygon is a GENERATOR, so a zone drawn on
    the rim hands the pond neighbours of its own and its cell can only be a
    fraction of the water. The zone is a ring around it instead."""
    from shapely.geometry import Polygon

    ponds = [(None, 400.0, 100.0, 200.0, 11.28)]
    zone = meshes._pond_zones(ponds, Polygon)[0]
    assert zone.contains(Polygon(_square(100.0, 200.0, 10.0))), \
        'the zone does not even cover the pond'
    # every vertex stands off by the same factor
    for x, y in list(zone.exterior.coords)[:-1]:
        d = ((x - 100.0) ** 2 + (y - 200.0) ** 2) ** 0.5
        assert abs(d - meshes.POND_ZONE_FACTOR * 11.28) < 1e-6


def test_the_zone_ring_is_coarse_enough_to_leave_the_pond_alone(cfg):
    """Its vertices are the generators the pond cell is bounded against, so
    a finely drawn ring would chop the cell up."""
    from shapely.geometry import Polygon

    zone = meshes._pond_zones([(None, 400.0, 0.0, 0.0, 11.28)], Polygon)[0]
    assert len(list(zone.exterior.coords)) - 1 == meshes.POND_ZONE_N
    assert meshes.POND_ZONE_N <= 16, 'too many generators around one pond'


def test_a_pond_size_is_refused_without_the_seeding(cfg):
    cfg.grid.voronoi.refine_ponds = False
    cfg.grid.voronoi.cell_pond = 30.0
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'refine_ponds' in str(e.value)


def test_a_pond_cell_may_not_be_coarser_than_the_background(cfg):
    """Past the background it is not a refinement of anything -- it is a hole
    the cells around it have to grade to reach."""
    cfg.grid.ponds = 'lm_ponds.shp'
    cfg.grid.voronoi.refine_ponds = True
    cfg.grid.voronoi.cell_pond = float(cfg.grid.voronoi.cell_far) + 1.0
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'cell_far' in str(e.value)


def test_the_quadtree_pond_refinement_needs_a_pond_layer(cfg):
    cfg.grid.ponds = ''
    cfg.grid.quadtree.refine_ponds = True
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'grid.ponds' in str(e.value)


def test_the_pond_level_defaults_to_the_stream_level(cfg):
    """0 means "the same as the streams", so the ponds follow them unless
    there is a reason to split further."""
    assert cfg.grid.quadtree.pond_level == 0
    cfg.grid.quadtree.pond_level = -1
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'pond_level' in str(e.value)


# -------------------------------------------------- panel 2: the records

def test_the_surface_rows_name_fields_that_exist(cfg):
    have = dict(schema.fields_of(cfg, 'surface'))
    for row in schema.SURFACE_ROWS:
        for dotted in row:
            if dotted is None:
                continue            # a deliberately empty cell
            assert dotted in have, '%s is laid out but does not exist' % dotted


def test_the_meteorology_and_the_irrigation_are_one_per_column():
    """One column per SUBJECT. Irrigation is the longer story and it is
    optional, so keeping it in one column is what makes the short one
    readable."""
    left = [row[0] for row in schema.SURFACE_ROWS if row[0]]
    right = [row[1] for row in schema.SURFACE_ROWS if len(row) > 1 and row[1]]
    assert left == ['surface.meteo_ts', 'surface.meteo_zones',
                    'surface.out_prefix']
    assert right == ['surface.irrigation', 'surface.irr_zones',
                     'surface.nfield', 'surface.irr_ts',
                     'surface.crop_schedule']
    assert 'surface.plot' not in left + right, \
        'the MMsurf figures belong to panel 4'


def test_every_record_field_offers_the_file_dialog():
    assert set(schema.SURFACE_FILES) == {'surface.meteo_ts', 'surface.irr_ts',
                                         'surface.crop_schedule'}
    for dotted in schema.SURFACE_FILES:
        assert dotted in [d for row in schema.SURFACE_ROWS for d in row]


def test_a_picked_schedule_becomes_a_pattern():
    """One file per field: storing the file that was pointed at would make
    every field read field 1's schedule."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        'panelui_t', os.path.join(CODE, 'app', 'lib', 'panelui.py'))
    if spec is None:                       # pragma: no cover
        pytest.skip('panelui not importable')
    try:
        mod = importlib.util.module_from_spec(spec)
        sys.modules['panelui_t'] = mod
        spec.loader.exec_module(mod)
    except ImportError:
        pytest.skip('streamlit not installed')

    got, note = mod._as_pattern('__inputFIELD1_crop_schedule.txt')
    assert got == '__inputFIELD%d_crop_schedule.txt' and not note
    # already a pattern: left alone
    assert mod._as_pattern(got) == (got, '')
    # no field number at all: said, not silently accepted as a pattern
    kept, why = mod._as_pattern('schedule.txt')
    assert kept == 'schedule.txt' and 'no field number' in why


# --------------------------------- panel 1: the cells the packages will own

def _square_mesh(n=4, side=10.0):
    """An n x n square mesh as DISV gridprops, cells in row-major order."""
    verts, cell2d, vid = [], [], {}
    for j in range(n + 1):
        for i in range(n + 1):
            vid[(i, j)] = len(verts)
            verts.append([len(verts), i * side, j * side])
    for j in range(n):
        for i in range(n):
            ring = [vid[(i, j)], vid[(i + 1, j)], vid[(i + 1, j + 1)],
                    vid[(i, j + 1)]]
            cell2d.append([len(cell2d), (i + 0.5) * side, (j + 0.5) * side,
                           4] + ring)
    return {'vertices': verts, 'cell2d': cell2d, 'ncpl': n * n, 'nlay': 1}


def test_a_pond_owns_the_cells_whose_centre_is_in_the_water():
    """CdL's rule. A Voronoi cell belongs to the generator it surrounds, so
    "the centre is in the water" is the honest test of whether a cell is part
    of the lake."""
    gp = _square_mesh()
    ring = [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)]
    got = meshes.pond_cells(gp, [ring])
    assert len(got) == 1
    # the four cells of the lower-left quarter: centres at 5 and 15
    assert sorted(got[0]) == [0, 1, 4, 5]


def test_a_pond_smaller_than_its_cell_still_owns_one():
    """Without the fallback it would own nothing at all, and LAK would be
    handed an empty footprint."""
    gp = _square_mesh()
    tiny = [(21.0, 21.0), (22.0, 21.0), (22.0, 22.0), (21.0, 22.0)]
    got = meshes.pond_cells(gp, [tiny])
    assert len(got[0]) == 1
    assert got[0][0] == 10, 'the nearest cell is not the one it landed in'


def test_the_stream_cells_are_the_ones_the_line_runs_through():
    gp = _square_mesh()
    line = [[(1.0, 5.0), (39.0, 5.0)]]          # straight along row 0
    got = meshes.stream_cells(gp, line)
    assert got == [0, 1, 2, 3]


def test_no_stream_and_no_pond_is_not_an_error():
    gp = _square_mesh()
    assert meshes.stream_cells(gp, []) == []
    assert meshes.pond_cells(gp, []) == []


def test_touching_is_a_shared_corner_not_an_overlap():
    """The question WP4 asks: a lake that does not touch the network it
    drains into cannot be connected to it by a mover."""
    gp = _square_mesh()
    assert meshes.cells_touch(gp, [0], [1]), 'neighbours do not touch'
    assert meshes.cells_touch(gp, [0], [5]), 'diagonal corner is a touch'
    assert not meshes.cells_touch(gp, [0], [2]), 'cells apart do touch'
    assert not meshes.cells_touch(gp, [0], []), 'nothing cannot be touched'


def test_the_cell_helpers_agree_with_the_gridprops():
    gp = _square_mesh(n=3)
    cen = meshes.cell_centres(gp)
    polys = meshes.cell_polygons(gp)
    assert cen.shape == (9, 2) and len(polys) == 9
    assert len(polys[0]) == 4
    assert abs(cen[0][0] - 5.0) < 1e-9 and abs(cen[0][1] - 5.0) < 1e-9


# ------------------------------------------- panel 1: the level sentences

def test_the_refinement_sentence_is_the_level_on_screen():
    def get(dotted, default=None):
        return {'grid.cell_size': 50.0, 'grid.quadtree.refine_level': 4,
                'grid.quadtree.pond_level': 0}.get(dotted, default)

    streams = schema.GRID_NOTES['grid.quadtree.refine_level'](get)
    assert 'level 4' in streams and '50 m background' in streams
    assert 'gives 3.125 m along the streams' in streams

    ponds = schema.GRID_NOTES['grid.quadtree.pond_level'](get)
    assert '0 follows the streams' in ponds
    assert 'gives 3.125 m over the pond footprints' in ponds


def test_a_pond_level_of_its_own_is_said_as_such():
    def get(dotted, default=None):
        return {'grid.cell_size': 50.0, 'grid.quadtree.refine_level': 2,
                'grid.quadtree.pond_level': 4}.get(dotted, default)

    ponds = schema.GRID_NOTES['grid.quadtree.pond_level'](get)
    assert ponds.startswith('Level 4 ')
    assert '3.125 m over the pond footprints' in ponds
    assert 'follows the streams' not in ponds


def test_every_note_belongs_to_a_field_that_is_shown():
    shown = set(schema.grid_fields('quadtree')) | set(
        schema.grid_fields('voronoi'))
    for dotted in schema.GRID_NOTES:
        assert dotted in shown, '%s has a note but no box' % dotted


# ----------------------------------------------- a key that has been renamed

def _toml_with(tmp_path, old, new):
    """A copy of the reference configuration with one key spelled the old
    way."""
    text = io.open(REF, encoding='utf-8').read().replace(new, old)
    assert old in text, 'the reference file does not carry %r' % new
    p = os.path.join(str(tmp_path), 'old.toml')
    io.open(p, 'w', encoding='utf-8', newline='').write(text)
    return p


def test_a_file_with_the_old_name_still_loads(tmp_path):
    """Unknown keys raise, so a rename has to be declared or every existing
    file stops loading."""
    p = _toml_with(tmp_path, 'seed_ponds = true', 'refine_ponds = true')
    cfg = cfgmod.load_run_config(p)
    assert cfg.grid.voronoi.refine_ponds is True


def test_the_old_name_is_reported_not_silently_accepted(tmp_path):
    """The name in the file is not the name on the panel until something
    saves it back, and a modeller reading the TOML has to know that."""
    p = _toml_with(tmp_path, 'seed_ponds = true', 'refine_ponds = true')
    cfg = cfgmod.load_run_config(p)
    said = ' '.join(getattr(cfg, 'migrated', []))
    assert 'seed_ponds' in said and 'refine_ponds' in said, said
    assert cfgmod.load_run_config(REF).migrated == [], \
        'a file already using the new name should report nothing'


def test_saving_migrates_the_file_for_good(tmp_path):
    p = _toml_with(tmp_path, 'seed_ponds = true', 'refine_ponds = true')
    out = os.path.join(str(tmp_path), 'again.toml')
    cfgmod.load_run_config(p).write_toml(out)
    body = io.open(out, encoding='utf-8').read()
    assert 'refine_ponds' in body and 'seed_ponds' not in body


def test_both_names_at_once_is_refused(tmp_path):
    """One of them would win silently, and it would be the wrong one half
    the time."""
    text = io.open(REF, encoding='utf-8').read().replace(
        'refine_ponds = true', 'refine_ponds = true\nseed_ponds = false')
    p = os.path.join(str(tmp_path), 'both.toml')
    io.open(p, 'w', encoding='utf-8', newline='').write(text)
    with pytest.raises(cfgmod.ConfigError) as e:
        cfgmod.load_run_config(p)
    assert 'seed_ponds' in str(e.value) and 'refine_ponds' in str(e.value)


def test_a_key_that_is_simply_wrong_still_raises(tmp_path):
    """The migration must not become a door for typos."""
    text = io.open(REF, encoding='utf-8').read().replace(
        'refine_ponds = true', 'refine_pnds = true')
    p = os.path.join(str(tmp_path), 'typo.toml')
    io.open(p, 'w', encoding='utf-8', newline='').write(text)
    with pytest.raises(cfgmod.ConfigError) as e:
        cfgmod.load_run_config(p)
    assert 'refine_pnds' in str(e.value)


def test_the_soil_and_the_vegetation_are_separate_subjects(cfg):
    """[soil] carries both because they share a FILE, not because they are
    one question -- and shown together the vegetation read as more soil."""
    have = dict(schema.fields_of(cfg, 'soil'))
    soil = [d for row in schema.SOIL_ROWS for d in row if d]
    veg = [d for row in schema.SOIL_VEG_ROWS for d in row if d]
    assert soil == ['soil.params', 'soil.zones', 'soil.thickness']
    assert veg == ['soil.veg_layer', 'soil.veg_column']
    assert not set(soil) & set(veg), 'a field is on both sub-panels'
    for dotted in soil + veg:
        assert dotted in have, '%s is laid out but does not exist' % dotted
    # nothing of [soil] is left with no home (veg_class is a TABLE)
    laid = set(soil) | set(veg)
    missing = [d for d in have if d not in laid and not d.endswith('veg_class')]
    assert not missing, 'no sub-panel shows %s' % missing


def test_the_vegetation_layer_is_chosen_like_every_other_shapefile():
    assert schema.SOIL_FILES == ('soil.veg_layer',)


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
