# -*- coding: utf-8 -*-
"""WP1d -- every panel renders, headlessly.

``streamlit.testing.v1.AppTest`` runs a page with no browser and reports the
exceptions it raised. That is the only way to catch a NameError on panel 3
without clicking through panels 0 to 2 -- and the failure mode these pages
are most prone to is not logic but plumbing: a duplicate widget key, a field
that is not what the widget expects, a section that no longer exists.

Skipped where streamlit is not installed, which is every machine that only
runs the model (a PEST worker, Spyder). The model must never need it.
"""

import io
import os

import pytest

pytest.importorskip('streamlit', reason='streamlit not installed')
from streamlit.testing.v1 import AppTest            # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
APP = os.path.join(CODE, 'app')

PAGES = ['Home.py'] + [os.path.join('pages', f) for f in (
    '1_Grid.py', '2_Surface.py', '3_Model.py', '4_Plots.py',
    '5_Run.py', '6_Results.py')]


@pytest.mark.parametrize('page', PAGES)
def test_the_page_renders_without_an_exception(page):
    path = os.path.join(APP, page)
    if not os.path.exists(path):
        pytest.skip('%s not present' % page)
    at = AppTest.from_file(path, default_timeout=180)
    at.run()
    problems = ['%s' % str(e.value).splitlines()[0] for e in at.exception]
    assert not problems, '%s raised:\n  %s' % (page, '\n  '.join(problems))


def test_the_panels_are_numbered_in_the_modellers_order():
    """The sidebar order is the file order, so it is the panel order."""
    names = sorted(f for f in os.listdir(os.path.join(APP, 'pages'))
                   if f.endswith('.py') and not f.startswith('_'))
    assert names == ['1_Grid.py', '2_Surface.py', '3_Model.py', '4_Plots.py',
                     '5_Run.py', '6_Results.py']


def test_the_panels_offer_something_to_edit():
    """A panel with no widgets is a panel that cannot do its job."""
    for page, least in (('pages/1_Grid.py', 8), ('pages/2_Surface.py', 6),
                        ('pages/3_Model.py', 20), ('pages/4_Plots.py', 6)):
        at = AppTest.from_file(os.path.join(APP, page), default_timeout=180)
        at.run()
        n = (len(at.text_input) + len(at.number_input) + len(at.checkbox)
             + len(at.selectbox))
        assert n >= least, '%s shows only %d widget(s)' % (page, n)


def test_the_master_switches_are_on_their_panels():
    """Surface, Model and Plots each carry their group's on/off."""
    for page in ('pages/2_Surface.py', 'pages/3_Model.py', 'pages/4_Plots.py'):
        at = AppTest.from_file(os.path.join(APP, page), default_timeout=180)
        at.run()
        keys = [t.key for t in at.toggle]
        assert any(k and k.startswith('sw_run.') for k in keys), \
            '%s has no master switch' % page


def test_the_grid_kind_is_a_choice_not_a_text_box():
    """Typing 'voroni' into a text box and finding out at run time is exactly
    what the panels exist to prevent."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    boxes = {s.key: list(s.options) for s in at.selectbox if s.key}
    kind = boxes.get('grid.kind')
    assert kind is not None, 'grid.kind is not a selectbox'
    assert set(kind) >= {'structured', 'disv', 'voronoi', 'quadtree'}


def _keys(at):
    """Every keyed widget on the page, whatever its type."""
    return {w.key for w in (list(at.number_input) + list(at.checkbox)
                            + list(at.selectbox) + list(at.text_input))
            if w.key}


def test_the_grid_subpanel_follows_the_kind():
    """Each kind shows ITS OWN settings and no one else's -- a setting that
    does nothing for the chosen producer must not sit beside one that does,
    looking equally live. Written without assuming which kind the shipped
    configuration selects."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    at.selectbox(key='grid.kind').select('structured').run()
    keys = _keys(at)
    assert not any(k.startswith('grid.voronoi.') for k in keys), \
        'the voronoi block is shown for a structured grid'
    assert not any(k.startswith('grid.quadtree.') for k in keys)
    # ... and the structured-only settings ARE there. They used to sit in the
    # permanent block, where they looked as though they applied to a mesh.
    assert 'grid.cell_size' in keys and 'grid.override.enable' in keys

    at.selectbox(key='grid.kind').select('voronoi').run()
    keys = _keys(at)
    assert any(k.startswith('grid.voronoi.') for k in keys), \
        'the voronoi block does not appear when voronoi is chosen'
    assert 'grid.override.enable' not in keys, \
        'the override is offered on a mesh, where validate() refuses it'
    assert 'grid.cell_size' not in keys, \
        'cell_size is offered for voronoi, which does not use it'

    at.selectbox(key='grid.kind').select('quadtree').run()
    keys = _keys(at)
    assert any(k.startswith('grid.quadtree.') for k in keys)
    assert not any(k.startswith('grid.voronoi.') for k in keys)
    assert 'grid.cell_size' in keys, 'the quadtree background is not offered'


def test_the_derived_transition_bands_are_read_only():
    """They are computed from the corridor and the grade ratio, so a box that
    accepted a value would be a box that lies."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    at.selectbox(key='grid.kind').select('voronoi').run()
    ro = [t for t in at.text_input if t.key == 'ro_grid.voronoi.trans_levels']
    assert ro, 'the derived bands are not shown'
    assert ro[0].disabled, 'the derived bands are editable'
    assert 'grid.voronoi.trans_levels' not in _keys(at), \
        'the derived bands are ALSO offered as a live field'


def test_the_refinement_settings_are_blocked_when_it_is_off():
    """Panel 1, D4: with the refinement off the corridor does not exist, so
    the settings describing it are greyed and cleared, not left looking live."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    at.selectbox(key='grid.kind').select('voronoi').run()
    assert 'grid.voronoi.cell_near_stream' in _keys(at)
    at.checkbox(key='grid.voronoi.stream_refine').uncheck().run()
    assert 'grid.voronoi.cell_near_stream' not in _keys(at), \
        'the corridor width is still live with the refinement off'
    blocked = [t for t in at.text_input
               if t.key == 'off_grid.voronoi.cell_near_stream']
    assert blocked and blocked[0].disabled and not blocked[0].value


def test_the_grid_panel_offers_both_buttons():
    """Create is the experiment, Select is the commitment -- and nothing else
    on this panel writes, because here the save IS the selection."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    keys = {b.key for b in at.button}
    assert 'mkgrid' in keys, 'no Create grid button'
    assert not any(k and k.startswith('save_') for k in keys), \
        'panel 1 still has a Validate & save'


def test_the_panel_says_when_the_grid_leaves_the_model_rasters():
    """WP1d stopgap. The rectangle panel 1 derives from the polygon and the
    one the model's rasters were frozen on are not the same, and the
    projection that discovers it fails deep inside, talking about tops and
    bottoms. It is stated HERE, where the rectangle is chosen.

    Skipped where the dataset is not on the machine -- the check has nothing
    to compare against and correctly says so.
    """
    import marmites_config as mcfg
    import marmites_meshes as mmesh
    import mm_paths

    cfg = mcfg.load_run_config(os.path.join(CODE, 'configs', 'lamata.toml'))
    ds = str(mm_paths.dataset_dir(cfg.paths.case))
    rect, _names, _others = mmesh.dataset_rectangle(ds)
    if rect is None:
        pytest.skip('no dataset raster on this machine')

    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    said = ' '.join(w.value for w in at.warning) + \
        ' '.join(c.value for c in at.caption)
    assert 'rasters' in said, 'the panel says nothing about the rasters'
    # and it is drawn BEFORE the build, not after it
    assert 'mkgrid' in {b.key for b in at.button}


def test_pinning_the_grid_takes_the_rectangle_from_the_rasters():
    """La Mata's configuration already carries the legacy origin, so a test
    run against it could not tell the button from a coincidence. This one
    starts from an override that is deliberately WRONG and checks the four
    boxes come back holding the rasters' own numbers.
    """
    import shutil

    import marmites_config as mcfg
    import marmites_meshes as mmesh
    import mm_paths

    ref = os.path.join(CODE, 'configs', 'lamata.toml')
    cfg = mcfg.load_run_config(ref)
    rect, _names, _others = mmesh.dataset_rectangle(
        str(mm_paths.dataset_dir(cfg.paths.case)))
    if rect is None:
        pytest.skip('no dataset raster on this machine')

    tmp = os.path.join(CODE, 'configs', '_pintest.toml')
    shutil.copyfile(ref, tmp)
    try:
        text = io.open(tmp, encoding='utf-8').read()
        text = text.replace('xllcorner = 739300.0', 'xllcorner = 111111.0')
        text = text.replace('nrow = 65', 'nrow = 7')
        io.open(tmp, 'w', encoding='utf-8', newline='').write(text)

        at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                               default_timeout=180)
        at.session_state['config_file'] = '_pintest.toml'
        at.run()
        [s for s in at.selectbox if s.key == 'grid.kind'][0] \
            .set_value('structured').run()
        hit = [b for b in at.button if b.key == 'adoptrect']
        assert hit, 'the pin was not offered on a structured grid'
        hit[0].click().run()

        got = dict((n.key, n.value) for n in at.number_input)
        assert got['grid.override.xllcorner'] == rect[0]
        assert got['grid.override.yllcorner'] == rect[1]
        assert got['grid.override.nrow'] == rect[2]
        assert got['grid.override.ncol'] == rect[3]
        assert got['grid.cell_size'] == rect[4]
        assert [c for c in at.checkbox
                if c.key == 'grid.override.enable'][0].value is True
        # and the page now says so, on the same run
        assert any('stands on' in c.value for c in at.caption)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_pin_is_not_offered_on_a_mesh():
    """The override reproduces a legacy DIS grid and means nothing on an
    unstructured mesh, so there the answer is the model panel, not a button
    that would be refused by validate()."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    kind = [s for s in at.selectbox if s.key == 'grid.kind'][0]
    for mesh in ('voronoi', 'quadtree'):
        kind.set_value(mesh).run()
        assert not [b for b in at.button if b.key == 'adoptrect'], mesh


def test_each_layer_picker_keeps_its_title_and_its_help():
    """The four titles are drawn by hand so they align across the columns,
    with the widget's own label collapsed underneath -- which is what took
    the help icon away the first time. It has to ride on the title instead,
    and the bracket has to be in the SAME block or a paragraph gap opens
    between them.
    """
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    want = {'grid.boundary': ('Catchment boundary', '(polygon)'),
            'grid.streams': ('Hydrography layer', '(lines)'),
            'grid.ponds': ('Pond layer', '(polygons)'),
            'grid.dem': ('Elevation raster', '(DEM)')}
    blocks = [(m.value, getattr(m, 'help', None)
               or getattr(getattr(m, 'proto', None), 'help', ''))
              for m in at.markdown]
    for dotted, (title, bracket) in want.items():
        hit = [(v, h) for v, h in blocks if v.startswith('**%s**' % title)]
        assert hit, 'no title block for %s' % dotted
        value, help_ = hit[0]
        assert bracket in value, '%s: the bracket left the title block' % dotted
        assert '\n' in value.split(bracket)[0], \
            '%s: the bracket is not on its own line' % dotted
        assert dotted in (help_ or ''), '%s: the help icon is gone' % dotted


def test_the_quadtree_sentence_follows_the_boxes():
    """It is the sentence the refined size is read off, so quoting a level
    and a background that are no longer on screen is worse than saying
    nothing."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    [s for s in at.selectbox if s.key == 'grid.kind'][0] \
        .set_value('quadtree').run()

    def said():
        got = [c.value for c in at.caption if 'GRIDGEN halves' in c.value]
        assert got, 'the quadtree sentence is not on the page'
        return got[0]

    [n for n in at.number_input if n.key == 'grid.cell_size'][0] \
        .set_value(80.0).run()
    assert 'on a 80 m background' in said(), said()
    lv = [n for n in at.number_input if n.key == 'grid.quadtree.refine_level']
    if lv:                       # drawn only while the refinement is on
        lv[0].set_value(3).run()
        assert 'level 3' in said() and 'gives 10 m' in said(), said()


def test_the_gis_folder_box_says_where_the_folder_comes_from():
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    labels = [x.label for x in at.text_input]
    assert 'Folder with GIS information (defined in panel Home)' in labels, \
        labels


def test_the_map_can_draw_the_pond_footprints():
    """A pond is a polygon the mesh has to cover, and a dot says nothing
    about whether it does -- so the overlay reads the GeoJSON outlines."""
    import marmites_config as mcfg
    import marmites_meshes as mmesh
    import mm_paths

    cfg = mcfg.load_run_config(os.path.join(CODE, 'configs', 'lamata.toml'))
    ds = str(mm_paths.dataset_dir(cfg.paths.case))
    rings = mmesh.pond_rings(ds)
    if not rings:
        pytest.skip('no pond table on this machine')
    assert all(len(r) >= 3 for r in rings), 'a footprint came back as a point'

    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    sel = [m for m in at.multiselect if 'Overlay' in str(m.label)]
    assert sel, 'the overlay chooser is gone'
    sel[0].set_value(['ponds']).run()
    assert not at.exception, [str(e.value) for e in at.exception]


def test_the_pond_settings_are_on_their_sub_panels():
    """The pond cell size belongs to voronoi and the pond refinement to the
    quadtree, each greyed until its own switch is on -- a size for a pond
    nothing is seeding is a box that does nothing."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    def pick(kind):
        # Re-queried every time: an AppTest element is a SNAPSHOT of one run,
        # and reusing a handle across a rerun looks up a widget that no
        # longer exists.
        [s for s in at.selectbox if s.key == 'grid.kind'][0] \
            .set_value(kind).run()

    def seed(on):
        got = [c for c in at.checkbox if c.key == 'grid.voronoi.refine_ponds']
        if got:
            (got[0].check() if on else got[0].uncheck()).run()
        return bool(got)

    pick('voronoi')
    keys = _keys(at)
    assert ('grid.voronoi.cell_pond' in keys
            or 'off_grid.voronoi.cell_pond' in keys), \
        'the pond cell size is not on the voronoi sub-panel'
    if seed(False):
        assert 'grid.voronoi.cell_pond' not in _keys(at), \
            'the pond size is live with nothing seeding a pond'
        seed(True)
        assert 'grid.voronoi.cell_pond' in _keys(at), \
            'the pond size did not come back with the seeding'

    pick('quadtree')
    keys = _keys(at)
    assert ('grid.quadtree.refine_ponds' in keys
            or 'na_grid.quadtree.refine_ponds' in keys), \
        'the quadtree cannot be told to refine at the ponds'
    assert not any(k.startswith('grid.voronoi.') for k in keys)


def test_panel_two_names_its_folder_and_offers_the_dialog():
    """The folder was a caption stating a path that could not be changed --
    fine until the records are somewhere else. And a record is a file on this
    machine, so it is chosen, not typed correctly."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '2_Surface.py'),
                           default_timeout=180)
    at.run()
    assert not at.exception, [str(e.value) for e in at.exception]
    labels = [x.label for x in at.text_input]
    assert 'Folder with SURFACE information' in labels, labels
    keys = {b.key for b in at.button if b.key}
    for dotted in ('surface.meteo_ts', 'surface.irr_ts',
                   'surface.crop_schedule'):
        assert dotted + '.__pick' in keys, '%s has no ... button' % dotted


def test_the_mmsurf_figures_moved_to_the_plots_panel():
    """A plotting choice, asked with the rest of them rather than among the
    records that feed the run."""
    two = AppTest.from_file(os.path.join(APP, 'pages', '2_Surface.py'),
                            default_timeout=180)
    two.run()
    assert 'surface.plot' not in _keys(two), \
        'the MMsurf figures are still on panel 2'

    four = AppTest.from_file(os.path.join(APP, 'pages', '4_Plots.py'),
                             default_timeout=180)
    four.run()
    assert not four.exception, [str(e.value) for e in four.exception]
    assert 'surface.plot' in _keys(four), \
        'the MMsurf figures did not arrive on panel 4'


def test_the_records_are_stated_under_the_table_that_reads_them():
    """A block of green and red dots says whether files are there, and
    nothing about what they are for. Under the table that reads them it is
    one answer -- and it follows the BOXES, since a green light against the
    file named before the last edit is worse than no light at all."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '2_Surface.py'),
                           default_timeout=180)
    at.run()

    def said():
        return ' '.join(c.value for c in at.caption)

    assert '__meteoTB.txt' in said(), said()[:300]
    assert 'Are they there?' not in ' '.join(m.value for m in at.markdown), \
        'the old block of dots is still on the panel'
    [x for x in at.text_input if x.key == 'surface.meteo_ts'][0] \
        .set_value('nowhere_at_all.txt').run()
    assert 'nowhere_at_all.txt' in said(), 'the lines ignored the box'
    assert '🔴' in said(), 'a file that is not there came up green'


def _producer(at, dotted, want):
    """Set a source's producer and rerun. Re-queried, since an AppTest
    element is a snapshot of one run."""
    [s for s in at.selectbox if s.key == dotted + '.__producer'][0] \
        .set_value(want).run()


def test_a_source_set_to_a_layer_offers_the_gis_dialog():
    """A layer is a shapefile in the cartography folder, so it is chosen the
    way every other shapefile on these panels is chosen."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '2_Surface.py'),
                           default_timeout=180)
    at.run()
    for dotted in ('surface.meteo_zones', 'surface.irr_zones'):
        _producer(at, dotted, 'layer')
        keys = {b.key for b in at.button if b.key}
        assert dotted + '.__v.__pick' in keys, \
            '%s as a layer has no ... button' % dotted
        assert any(x.key == dotted + '.__v' for x in at.text_input), \
            '%s as a layer is not a path box' % dotted


def test_a_zone_set_to_a_value_is_a_count():
    """1.0 zones is not a thing, and a box that offers decimals invites
    one."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '2_Surface.py'),
                           default_timeout=180)
    at.run()
    for dotted in ('surface.meteo_zones', 'surface.irr_zones'):
        _producer(at, dotted, 'value')
        box = [n for n in at.number_input if n.key == dotted + '.__v']
        assert box, '%s as a value is not a number box' % dotted
        assert isinstance(box[0].value, int), \
            '%s offers decimals for a zone' % dotted
        assert not [b for b in at.button
                    if b.key == dotted + '.__v.__pick'], \
            'a number has no file to choose'


def test_a_measurement_is_still_a_measurement():
    """soil.thickness is a source too and its value is METRES -- the count
    rule must not reach it, or a 0.8 m soil becomes 1 m."""
    import sys

    if CODE not in sys.path:
        sys.path.insert(0, CODE)
    from lib import schema as sch

    assert 'soil.thickness' not in sch.INTEGER_VALUE
    at = AppTest.from_file(os.path.join(APP, 'pages', '3_Model.py'),
                           default_timeout=180)
    at.run()
    vals = [n.value for n in at.number_input
            if n.key and n.key.endswith('.__v')]
    assert any(isinstance(v, float) and v != int(v) for v in vals), \
        'panel 3 lost its fractional source values'


def test_every_producer_that_names_a_file_offers_the_dialog():
    """By rule, not by listing the fields one at a time: a layer is a
    shapefile in the cartography folder and a raster is a file in the
    dataset, so both are chosen rather than typed correctly. A column names
    an attribute and a value is a number, so neither is."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '3_Model.py'),
                           default_timeout=180)
    at.run()
    keys = {b.key for b in at.button if b.key}
    for dotted, want in (('soil.zones', 'layer'), ('soil.thickness', 'raster')):
        [s for s in at.selectbox
         if s.key == dotted + '.__producer'][0].set_value(want).run()
        keys = {b.key for b in at.button if b.key}
        assert dotted + '.__v.__pick' in keys, (
            '%s as a %s has no ... button' % (dotted, want))
    # A column is not a file. Asked of a ParamSource, since a VectorSource
    # no longer offers column as a producer at all -- there it is an
    # attribute of the layer and is drawn under it.
    [s for s in at.selectbox
     if s.key == 'sfr.width.__producer'][0].set_value('column').run()
    assert 'sfr.width.__v.__pick' not in {b.key for b in at.button}, (
        'a column name is being chosen from the filesystem')
    # ... and a value is not one either
    [s for s in at.selectbox
     if s.key == 'soil.thickness.__producer'][0].set_value('value').run()
    assert 'soil.thickness.__v.__pick' not in {b.key for b in at.button}, (
        'a number is being chosen from the filesystem')


def test_the_soil_parameter_file_is_chosen_too():
    at = AppTest.from_file(os.path.join(APP, 'pages', '3_Model.py'),
                           default_timeout=180)
    at.run()
    assert 'soil.params.__pick' in {b.key for b in at.button if b.key}
    assert any(x.key == 'soil.params' for x in at.text_input)


def test_a_file_below_the_folder_keeps_its_relative_path():
    """The soil parameters live in MF_ws/ under the dataset. Storing an
    absolute path for them would tie the configuration to this machine."""
    import importlib.util
    import sys as _sys

    spec = importlib.util.spec_from_file_location(
        'panelui_rel', os.path.join(APP, 'lib', 'panelui.py'))
    mod = importlib.util.module_from_spec(spec)
    _sys.modules['panelui_rel'] = mod
    spec.loader.exec_module(mod)
    src = io.open(os.path.join(APP, 'lib', 'panelui.py'),
                  encoding='utf-8').read()
    assert 'os.path.relpath(got, os.path.abspath(folder))' in src
    assert "rel.replace(os.sep, '/')" in src, (
        'a stored path must use / so it reads the same on either platform')


def test_irrigation_off_makes_its_fields_read_only():
    """A panel must not be left saying something the run will not do. The
    fields are SHOWN and greyed, not cleared: a switch that is off is a
    statement about the run, not about the answers."""
    gated = ('surface.irr_zones', 'surface.nfield', 'surface.irr_ts',
             'surface.crop_schedule')
    at = AppTest.from_file(os.path.join(APP, 'pages', '2_Surface.py'),
                           default_timeout=180)
    at.run()

    def boxes():
        return dict((w.key, w) for w in
                    (list(at.text_input) + list(at.number_input)
                     + list(at.checkbox) + list(at.selectbox)) if w.key)

    was = boxes()
    assert 'surface.irr_ts' in was, 'the series is not live with irrigation on'
    before = was['surface.irr_ts'].value

    [c for c in at.checkbox if c.key == 'surface.irrigation'][0] \
        .uncheck().run()
    now = boxes()
    for dotted in gated:
        ro = now.get('ro_' + dotted)
        assert ro is not None, '%s vanished instead of greying' % dotted
        assert ro.disabled, '%s is still editable' % dotted
        assert not any(k == dotted or k.startswith(dotted + '.')
                       for k in now), '%s is still live' % dotted
    assert str(before) in str(now['ro_surface.irr_ts'].value), (
        'the greyed field forgot what it held')
    assert not [b for b in at.button
                if b.key == 'surface.irr_ts.__pick'], (
        'a read-only field still offers the file dialog')

    # ... and the meteorology is untouched by any of it
    assert 'surface.meteo_ts' in now and not now['surface.meteo_ts'].disabled


def test_the_irrigation_fields_come_back_with_their_values():
    """Emptying them on an unticked box would mean typing them again on the
    next tick."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '2_Surface.py'),
                           default_timeout=180)
    at.run()
    before = [x.value for x in at.text_input if x.key == 'surface.irr_ts'][0]
    [c for c in at.checkbox if c.key == 'surface.irrigation'][0] \
        .uncheck().run()
    [c for c in at.checkbox if c.key == 'surface.irrigation'][0] \
        .check().run()
    after = [x.value for x in at.text_input if x.key == 'surface.irr_ts']
    assert after and after[0] == before, 'the series was lost on the way'


def test_a_column_is_an_attribute_of_a_layer_not_an_alternative_to_one():
    """VectorSource.producer() never returns `column`, so offering it as a
    choice set the column, cleared the layer, and left a source producing
    nothing. It belongs under the layer it qualifies."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '3_Model.py'),
                           default_timeout=180)
    at.run()
    box = [s for s in at.selectbox if s.key == 'soil.zones.__producer'][0]
    assert 'column' not in list(box.options), \
        'a zone source still offers column as a producer'
    assert set(box.options) == {'raster', 'layer', 'value'}

    # on a layer, the column is there and holds what the file holds. Either
    # widget will do: it is a LIST when the layer can be read and a box when
    # it cannot, and this test is about where the column sits, not which.
    def column(at_):
        return [w for w in (list(at_.text_input) + list(at_.selectbox))
                if w.key == 'soil.zones.__col']

    box.set_value('layer').run()
    col = column(at)
    assert col and col[0].value, 'the attribute column is not shown'

    # on a raster there is no layer, so there is no column either
    [s for s in at.selectbox
     if s.key == 'soil.zones.__producer'][0].set_value('raster').run()
    assert not column(at), 'a column is still offered for a raster'


def test_a_param_source_keeps_column_as_a_producer():
    """There it means a column of the SFR source layer, and producer() does
    return it -- the rule is about what the class means, not about the word."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '3_Model.py'),
                           default_timeout=180)
    at.run()
    box = [s for s in at.selectbox if s.key == 'sfr.width.__producer']
    assert box, 'sfr.width is not a source any more'
    assert 'column' in list(box[0].options), \
        'a ParamSource lost its column producer'
    assert 'layer' not in list(box[0].options), \
        'a ParamSource has no layer of its own to name'


def test_the_attribute_column_is_chosen_from_the_layer():
    """The names are IN the file, so asking someone to remember GRID_CODE
    against GRIDCODE is asking them to make a mistake the panel could have
    prevented."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '3_Model.py'),
                           default_timeout=180)
    at.run()
    box = [s for s in at.selectbox if s.key == 'soil.zones.__col']
    assert box, 'the attribute column is not a list'
    opts = list(box[0].options)
    assert len(opts) > 2, 'the list holds nothing but the none entry: %s' % opts
    assert box[0].value in opts
    assert any('none' in str(o) for o in opts), (
        'no way to say "use the geometry"')


def test_a_layer_with_no_attributes_falls_back_to_a_box():
    """A name typed for a file still to be exported should not be thrown
    away, and a picker with nothing in it is worse than a box."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '3_Model.py'),
                           default_timeout=180)
    at.run()
    [x for x in at.text_input
     if x.key == 'soil.zones.__v'][0].set_value('not_here_at_all.shp').run()
    assert not [s for s in at.selectbox if s.key == 'soil.zones.__col'], (
        'a list is offered for a layer that is not there')
    assert [x for x in at.text_input if x.key == 'soil.zones.__col'], (
        'the column was taken away instead of falling back to a box')


def test_panel_zero_sets_the_machine_paths():
    """They used to be read-only in the sidebar, under a caption telling the
    modeller to go and edit code/mm_paths.py. A path is not source code."""
    at = AppTest.from_file(os.path.join(APP, 'Home.py'), default_timeout=180)
    at.run()
    keys = {w.key for w in at.text_input if w.key}
    for key in ('path_example_root', 'path_data_root', 'path_gis',
                'path_ws_root'):
        assert key in keys, '%s is not settable on panel 0' % key
    assert any(b.key == 'save_paths' for b in at.button), \
        'panel 0 collects the paths but cannot save them'
    said = ' '.join(str(m.value) for m in at.caption) + \
        ' '.join(str(m.value) for m in at.markdown)
    assert 'mm_paths.py`, or set' not in said, \
        'panel 0 still tells the modeller to edit a source file'


def test_panel_zero_offers_the_system_dialog_for_each_path():
    """"Select the folder inside the computer": every path has a ... button
    that opens the operating system's own dialog. The in-page folder browser
    that used to sit beside it is gone -- two ways to do one thing -- and the
    text box is the fallback when no window can be opened."""
    at = AppTest.from_file(os.path.join(APP, 'Home.py'), default_timeout=180)
    at.run()
    keys = {b.key for b in at.button if b.key}
    for k in ('example_root', 'data_root', 'gis', 'ws_root', 'nwt_ref'):
        assert 'path_%s.__native' % k in keys, '%s has no ... button' % k
    assert not any(k.endswith(('.__up', '.__use', '.__usef')) for k in keys), \
        'the in-page folder browser is still on panel 0'


def test_the_derived_bands_follow_the_boxes_without_a_save():
    """They are recomputed from what is ON SCREEN. A derived box that only
    caught up after a save shows the PREVIOUS corridor's bands, and the first
    thing doubted is the number rather than the box."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    at.selectbox(key='grid.kind').select('voronoi').run()
    def shown(a, what='trans_levels'):
        return [t.value for t in a.text_input
                if t.key == 'ro_grid.voronoi.%s' % what][0]

    before, corridor = shown(at), shown(at, 'stream_buffer')
    # A finer cell at the stream takes more bands, and a wider corridor ...
    at.number_input(key='grid.voronoi.cell_near_stream').set_value(10.0).run()
    assert shown(at) != before, 'the bands ignored the size at the stream'
    assert shown(at, 'stream_buffer') != corridor, \
        'the derived corridor ignored the size at the stream'
    # ... and a coarser grade asks for fewer of them. BOTH ratios are set
    # explicitly: which one the configuration starts at is not this test's
    # business, and assuming it is what broke this when the default moved.
    at.number_input(key='grid.voronoi.grade_ratio').set_value(1.2).run()
    n_fine = len(shown(at).split(','))
    at.number_input(key='grid.voronoi.grade_ratio').set_value(3.0).run()
    assert len(shown(at).split(',')) < n_fine, \
        'the bands ignored the grade ratio'
    # The corridor is a READ-OUT now, not a question (panel 1, item 9).
    assert 'grid.voronoi.stream_buffer' not in _keys(at), \
        'the corridor is still offered as a live field'


def _fake_attempt(tmpdir, tag, kind='structured', nrow=3, ncol=4, cell=50.0):
    """A built attempt, cached the way the producer caches one."""
    import json
    import sys
    sys.path.insert(0, CODE)
    import numpy as np
    from marmites_grid import disv_from_structured
    verts, cell2d, ncpl = disv_from_structured(
        np.full(ncol, cell), np.full(nrow, cell), 0.0, 0.0)
    d = os.path.join(str(tmpdir), tag)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, 'mesh_%s.json' % kind), 'w', encoding='utf-8') as f:
        json.dump({'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl,
                   'nlay': 1}, f)
    with open(os.path.join(d, 'mesh_%s.sig.json' % kind), 'w',
              encoding='utf-8') as f:
        json.dump({'signature': 'sig_' + tag, 'kind': kind, 'ncpl': ncpl}, f)
    return d, ncpl


def test_the_mesh_tab_draws_the_attempt_that_is_selected(tmp_path):
    """The combo moved here from the first tab because an attempt is chosen
    by LOOKING at it -- so choosing one has to change the map."""
    import sys
    sys.path.insert(0, CODE)
    import marmites_config as mcfg

    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    attempts = []
    for tag, nrow, ncol in (('structured_1', 3, 4), ('structured_2', 6, 8)):
        d, ncpl = _fake_attempt(tmp_path, tag, nrow=nrow, ncol=ncol)
        trial = mcfg.RunConfig.from_dict({'grid': {'kind': 'structured'}})
        attempts.append({'tag': tag, 'ok': True, 'cache': d, 'cfg': trial,
                         'label': '%d x %d' % (nrow, ncol),
                         'lines': ['built %d cells' % ncpl],
                         'info': {'kind': 'structured', 'ncpl': ncpl,
                                  'area_mean': 2500.0}})
    at.session_state['grid_attempts'] = attempts
    at.run()
    assert not at.exception, [str(e.value) for e in at.exception]

    box = at.selectbox(key='pick_attempt')
    assert box is not None, 'the attempt selector is not on the mesh tab'
    assert any('structured_1' in o for o in box.options)
    assert any('structured_2' in o for o in box.options)

    def cells():
        return [m.value for m in at.metric if 'cells' in str(m.label)][0]

    one = [o for o in box.options if 'structured_1' in o][0]
    at.selectbox(key='pick_attempt').select(one).run()
    assert cells() == '12', 'the map did not follow the selector'
    two = [o for o in box.options if 'structured_2' in o][0]
    at.selectbox(key='pick_attempt').select(two).run()
    assert cells() == '48', 'the map did not follow the selector'


def test_the_attempt_selector_is_not_on_the_first_tab(tmp_path):
    """It was moved, not copied: two selectors disagreeing about which grid
    is on screen is worse than either."""
    src = open(os.path.join(APP, 'pages', '1_Grid.py'), encoding='utf-8').read()
    before, after = src.split('with tab_mesh:', 1)
    assert 'pick_attempt' not in before
    assert 'pick_attempt' in after
    assert 'selgrid' not in before, \
        'Select this grid is still on the Catchment & grid tab'


def _static_map(at):
    """Turn the interactive map off, so the button controls are drawn.

    With plotly installed the map is interactive and carries its own zoom,
    pan and reset; the buttons are the fallback for a machine without it,
    and this is how the test reaches them either way.
    """
    for tg in at.toggle:
        if tg.key == 'interactive_map':
            at.toggle(key='interactive_map').set_value(False).run()
            break
    return at


def test_the_map_has_zoom_and_an_original_extent(tmp_path):
    """A matplotlib figure reaches the browser as a picture, so the view is a
    state the buttons move and the axes are set from."""
    import sys
    sys.path.insert(0, CODE)
    import marmites_config as mcfg

    d, ncpl = _fake_attempt(tmp_path, 'structured_1', nrow=4, ncol=4)
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.session_state['grid_attempts'] = [{
        'tag': 'structured_1', 'ok': True, 'cache': d,
        'cfg': mcfg.RunConfig.from_dict({'grid': {'kind': 'structured'}}),
        'label': '4 x 4', 'lines': ['built %d cells' % ncpl],
        'info': {'kind': 'structured', 'ncpl': ncpl, 'area_mean': 2500.0}}]
    at.run()
    assert not at.exception, [str(e.value) for e in at.exception]
    _static_map(at)
    keys = {b.key for b in at.button if b.key}
    for k in ('zin', 'zout', 'zreset', 'pleft', 'pright', 'pup', 'pdown'):
        assert k in keys, 'the map has no %r control' % k

    def view(a):
        return (a.session_state['view']
                if 'view' in a.session_state else None)

    assert view(at) is None                          # the full extent
    at.button(key='zin').click().run()
    zoomed = view(at)
    assert zoomed is not None, 'zooming in did not take'
    # Compared against the PREVIOUS view rather than against absolute
    # numbers: which mesh the selector lands on is not this test's business.
    at.button(key='zin').click().run()
    closer = view(at)
    assert (closer[2] - closer[0]) < (zoomed[2] - zoomed[0]), \
        'zooming in did not narrow the view'
    at.button(key='zout').click().run()
    assert round(view(at)[2] - view(at)[0], 3) == \
        round(zoomed[2] - zoomed[0], 3), 'zoom out did not undo zoom in'
    at.button(key='pright').click().run()
    panned = view(at)
    assert panned[0] > zoomed[0], 'panning did not move the view'
    assert round(panned[2] - panned[0], 6) == round(zoomed[2] - zoomed[0], 6), \
        'panning changed the zoom'
    at.button(key='zreset').click().run()
    assert view(at) is None, 'Original extent did not reset'


def test_a_different_mesh_resets_the_view(tmp_path):
    """A window from the previous grid would be meaningless on the next."""
    import sys
    sys.path.insert(0, CODE)
    import marmites_config as mcfg

    trial = mcfg.RunConfig.from_dict({'grid': {'kind': 'structured'}})
    attempts = []
    for tag, n in (('structured_1', 4), ('structured_2', 9)):
        d, ncpl = _fake_attempt(tmp_path, tag, nrow=n, ncol=n)
        attempts.append({'tag': tag, 'ok': True, 'cache': d, 'cfg': trial,
                         'label': '%d x %d' % (n, n),
                         'lines': ['built %d cells' % ncpl],
                         'info': {'kind': 'structured', 'ncpl': ncpl,
                                  'area_mean': 2500.0}})
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.session_state['grid_attempts'] = attempts
    at.run()
    _static_map(at)

    def view(a):
        return (a.session_state['view']
                if 'view' in a.session_state else None)
    at.button(key='zin').click().run()
    assert view(at) is not None, 'zooming in did not take'
    box = at.selectbox(key='pick_attempt')
    other = [o for o in box.options if 'structured_1' in o][0]
    at.selectbox(key='pick_attempt').select(other).run()
    assert view(at) is None, \
        'the zoom from the previous mesh was kept'


def test_the_cartography_box_is_on_the_panel_that_uses_it():
    """Panel 1 re-reads the two tables a GRID depends on inside Create grid,
    so it needs no button; the soil, vegetation, observation and pond layers
    belong to panel 3, and so does the button for them."""
    one = open(os.path.join(APP, 'pages', '1_Grid.py'),
               encoding='utf-8').read()
    three = open(os.path.join(APP, 'pages', '3_Model.py'),
                 encoding='utf-8').read()
    assert 'conv_run' not in one, 'panel 1 has a converter button again'
    assert 'conv_run' in three, 'panel 3 has lost the converter button'
    # In a TAB of its own: content outside the tabs is drawn under whichever
    # one is open, so at the foot of the page it read as part of the soil
    # column, which is exactly what it is not.
    assert "'Cartography'" in three and 'with tab_gis:' in three
    # ... and panel 1 still re-reads what a grid needs, by itself.
    assert '_dataset_stale' in one and 'gis_to_dataset.py' in one


def test_the_grid_tabs_are_named_consistently():
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    src = open(os.path.join(APP, 'pages', '1_Grid.py'),
               encoding='utf-8').read()
    assert "'Visualize grid'" in src
    assert 'Visualize mesh' not in src, 'a stale name is still referred to'


def test_panel_one_asks_for_all_three_layers():
    """The working flow is from scratch: the catchment, the hydrography and
    the ponds are all asked for HERE, at the same level, because the
    refinement options below cannot be answered before it is known whether
    there is a network or a pond to refine around."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    keys = _keys(at)
    for k in ('grid.boundary', 'grid.streams', 'grid.ponds'):
        assert k in keys, '%s is not asked for on panel 1' % k
    assert 'gis_folder' in keys, 'no folder to look in'


def test_the_optional_layers_can_be_set_to_none():
    """A catchment with no mapped network has to be able to say so."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    for k in ('grid.streams', 'grid.ponds'):
        box = at.selectbox(key=k)
        assert box is not None, '%s is not a choice' % k
        assert any('none' in str(o).lower() for o in box.options), \
            '%s cannot be left unset' % k
    # ... and the catchment cannot: it is what the grid is built inside.
    assert not any('none' in str(o).lower()
                   for o in at.selectbox(key='grid.boundary').options)


def test_the_refinement_needs_its_layer():
    """Panel 1 offers the switch DISABLED, with the reason, rather than
    offering it and having validate() refuse the save."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    at.selectbox(key='grid.kind').select('voronoi').run()
    assert 'grid.voronoi.stream_refine' in _keys(at), \
        'the refinement is not offered even with a stream layer set'

    # Take the layer away: the switch must go unavailable, not just unticked.
    box = at.selectbox(key='grid.streams')
    none = [o for o in box.options if 'none' in str(o).lower()][0]
    at.selectbox(key='grid.streams').select(none).run()
    assert 'grid.voronoi.stream_refine' not in _keys(at), \
        'the refinement is still live with no hydrography layer'
    blocked = [c for c in at.checkbox
               if c.key == 'na_grid.voronoi.stream_refine']
    assert blocked and blocked[0].disabled and blocked[0].value is False


def test_panel_one_asks_for_the_dem_as_a_raster():
    """It does not shape the grid, but it is where pond rim and bottom come
    from -- and it was hardcoded, which a new catchment cannot use."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_Grid.py'),
                           default_timeout=180)
    at.run()
    assert 'grid.dem' in _keys(at), 'the DEM is not asked for'
    box = at.selectbox(key='grid.dem')
    assert box is not None and any('none' in str(o).lower()
                                   for o in box.options), \
        'the DEM cannot be left unset'
    # It lists RASTERS: an ESRI grid folder is one, and no .shp belongs.
    named = [str(o) for o in box.options]
    assert not any(o.lower().endswith('.shp') for o in named), \
        'the DEM picker is offering shapefiles'
