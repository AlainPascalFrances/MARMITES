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


def _app_dir():
    """The app's folder, as AppTest may open it.

    Streamlit refuses a page on a network path ("Unable to create Page.
    Network paths are not supported") -- which is why the app runs from a
    local mirror (code/tools/launch_app.py). On such a checkout the tests
    open the pages from a mirror too: their OWN, refreshed here in the temp
    folder (never the one the user's app runs from), so what they test is
    the checkout's code byte for byte. A local checkout is used in place.
    (2026-10-09: Home.py and the two panel-0 tests failed on the server.)
    """
    import importlib.util
    import tempfile
    spec = importlib.util.spec_from_file_location(
        '_launch_app_for_tests', os.path.join(CODE, 'tools', 'launch_app.py'))
    la = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(la)
    if not la.on_network_drive(CODE):
        return os.path.join(CODE, 'app')
    mirror = os.path.join(tempfile.gettempdir(), 'MARMITES', 'test_mirror',
                          la.default_mirror(la.REPO).name)
    la.check_mirror_dir(mirror, la.REPO)
    la.sync(la.REPO, mirror)
    la.mark(la.REPO, mirror)
    return os.path.join(mirror, 'code', 'app')


APP = _app_dir()
_MIRROR = (os.path.dirname(os.path.dirname(APP))
           if os.path.dirname(APP) != CODE else None)


@pytest.fixture(autouse=True)
def _the_checkout_stays_first():
    """A page puts its code folders on sys.path (lib/runs.py adds tests/
    too). Run from the mirror, those are the MIRROR's, and a later test
    file then imported its helpers from there: lamata_model looked for the
    dataset beside the mirror and every La Mata test skipped (2026-10-09).
    So each test here leaves sys.path, and the modules it loaded from the
    mirror, as it found them."""
    import sys
    path, before = list(sys.path), set(sys.modules)
    yield
    sys.path[:] = path
    if _MIRROR:
        root = os.path.normcase(os.path.abspath(_MIRROR))
        for name in [n for n in list(sys.modules) if n not in before]:
            f = getattr(sys.modules.get(name), '__file__', None) or ''
            if f and os.path.normcase(os.path.abspath(f)).startswith(root):
                del sys.modules[name]

SURF = 'pages/2_2_-_Surface_and_driving_forces.py'
SOIL = 'pages/3_3_-_Soil.py'
SUB = 'pages/4_4_-_Unsaturated_zone_and_groundwater.py'
OBS = 'pages/5_5_-_State_variables.py'
PLOT = 'pages/6_6_-_Plots.py'

PAGES = ['Home.py'] + [os.path.join('pages', f) for f in (
    '1_1_-_Grid.py', '2_2_-_Surface_and_driving_forces.py', '3_3_-_Soil.py',
    '4_4_-_Unsaturated_zone_and_groundwater.py', '5_5_-_State_variables.py',
    '6_6_-_Plots.py', '7_7_-_Run.py', '8_8_-_Results.py')]


@pytest.mark.parametrize('page', PAGES)
def test_the_page_renders_without_an_exception(page):
    path = os.path.join(APP, page)
    if not os.path.exists(path):
        pytest.skip('%s not present' % page)
    at = AppTest.from_file(path, default_timeout=180)
    at.run()
    problems = ['%s' % str(e.value).splitlines()[0] for e in at.exception]
    assert not problems, '%s raised:\n  %s' % (page, '\n  '.join(problems))


def test_the_results_are_shown_by_tabs_two_to_a_row(tmp_path, monkeypatch):
    """User, 2026-10-09: tabs for input maps, output maps, time series,
    calibration, the total and the ponds' water budget; two figures to a
    row, aligned on top; the title centred below each figure."""
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import mm_paths
    run = tmp_path / 'out_202610091200_try'
    figs = {'_input': ['IN_000_model_map.png', '_sp_plt_IN_007_hk_1_1.png'],
            '_output': ['_sp_plt_GWmap_head_1_1.png', '_0P0_ts.png',
                        'obs_heads.png', 'budget_compartment.png',
                        'lake_budget_years_total.png',
                        'lake_budget_years_by_pond.png',
                        'lake_stage.png']}
    for sub, names in figs.items():
        (run / sub).mkdir(parents=True)
        for n in names:
            fig, ax = plt.subplots(figsize=(1, 1))
            fig.savefig(str(run / sub / n))
            plt.close(fig)
    monkeypatch.setattr(mm_paths, 'WS_ROOT', tmp_path)
    at = AppTest.from_file(os.path.join(APP, 'pages', '8_8_-_Results.py'),
                           default_timeout=180)
    at.run()
    assert not at.exception, [str(e.value) for e in at.exception]
    labels = [t.label for t in at.tabs]
    assert labels == ['Input maps (2)', 'Output maps (1)', 'Time series (1)',
                      'Calibration (state variables) (1)',
                      'Water budgets (1)', 'Ponds water budget (3)'], labels
    said = ' '.join(str(m.value) for m in at.markdown)
    assert 'text-align:center' in said
    assert 'Horizontal hydraulic conductivity' in said
    src = open(os.path.join(APP, 'pages', '8_8_-_Results.py'),
               encoding='utf-8').read()
    assert "st.columns(2, gap='medium', vertical_alignment='top')" in src
    assert "'Columns'" not in src, 'the column count is fixed at two'


def test_the_sidebar_shows_the_number_with_the_name():
    """THE NUMBER IS IN THE FILENAME TWICE, and this says why.

    Streamlit builds a page's sidebar label from its filename with
    ``([0-9]*)[_ -]*(.*)\\.py``: the leading digits ORDER the page and are
    then thrown away, so ``1_Grid.py`` appeared as plain "Grid". The number
    is what lets a modeller match the sidebar against a panel referred to by
    number, so it is written a second time inside the part that survives --
    ``1_1_-_Grid.py`` orders by 1 and shows "1 - Grid".

    Asked of streamlit's own function rather than of the names, so a change
    to that regex fails here instead of quietly renaming every panel in the
    sidebar.
    """
    from pathlib import Path

    from streamlit.source_util import page_icon_and_name

    want = [(1, 'Grid'), (2, 'Surface and driving forces'), (3, 'Soil'),
            (4, 'Unsaturated zone and groundwater'), (5, 'State variables'),
            (6, 'Plots'), (7, 'Run'), (8, 'Results')]
    names = sorted(f for f in os.listdir(os.path.join(APP, 'pages'))
                   if f.endswith('.py') and not f.startswith('_'))
    assert len(names) == len(want), names
    for fn, (num, title) in zip(names, want):
        _icon, raw = page_icon_and_name(Path(fn))
        assert raw.replace('_', ' ') == '%d - %s' % (num, title), (
            '%s would appear in the sidebar as %r'
            % (fn, raw.replace('_', ' ')))
    # Home is the ENTRY SCRIPT, not a page, so it carries no number.
    assert os.path.exists(os.path.join(APP, 'Home.py'))
    assert not any(f.startswith('0') for f in names), \
        'Home has been given a number'


def test_the_panels_are_numbered_in_the_modellers_order():
    """The sidebar order is the file order, so it is the panel order."""
    names = sorted(f for f in os.listdir(os.path.join(APP, 'pages'))
                   if f.endswith('.py') and not f.startswith('_'))
    assert names == ['1_1_-_Grid.py', '2_2_-_Surface_and_driving_forces.py',
                     '3_3_-_Soil.py',
                     '4_4_-_Unsaturated_zone_and_groundwater.py',
                     '5_5_-_State_variables.py', '6_6_-_Plots.py',
                     '7_7_-_Run.py', '8_8_-_Results.py']


def test_the_panels_offer_something_to_edit():
    """A panel with no widgets is a panel that cannot do its job."""
    for page, least in (('pages/1_1_-_Grid.py', 8), (SURF, 6), (SOIL, 6),
                        (SUB, 14), (OBS, 5), (PLOT, 6)):
        at = AppTest.from_file(os.path.join(APP, page), default_timeout=180)
        at.run()
        n = (len(at.text_input) + len(at.number_input) + len(at.checkbox)
             + len(at.selectbox))
        assert n >= least, '%s shows only %d widget(s)' % (page, n)


def test_the_master_switches_are_on_their_panels():
    """Surface, Model and Plots each carry their group's on/off."""
    for page in (SURF, SOIL, SUB, PLOT):
        at = AppTest.from_file(os.path.join(APP, page), default_timeout=180)
        at.run()
        keys = [t.key for t in at.toggle]
        assert any(k and k.startswith('sw_run.') for k in keys), \
            '%s has no master switch' % page


def test_the_grid_kind_is_a_choice_not_a_text_box():
    """Typing 'voroni' into a text box and finding out at run time is exactly
    what the panels exist to prevent."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
                           default_timeout=180)
    at.run()
    at.selectbox(key='grid.kind').select('voronoi').run()
    ro = [t for t in at.text_input if t.key == 'ro_grid.voronoi.trans_levels']
    assert ro, 'the derived bands are not shown'
    assert ro[0].disabled, 'the derived bands are editable'
    assert 'grid.voronoi.trans_levels' not in _keys(at), \
        'the derived bands are ALSO offered as a live field'


def _refine_on(at):
    """Stream refinement ON in this session, whatever the live file says (it
    was switched off there on 2026-09-27). Nothing is saved."""
    cb = at.checkbox(key='grid.voronoi.stream_refine')
    if not cb.value:
        cb.check().run()


def test_the_refinement_settings_are_blocked_when_it_is_off():
    """Panel 1, D4: with the refinement off the corridor does not exist, so
    the settings describing it are greyed and cleared, not left looking live."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
                           default_timeout=180)
    at.run()
    at.selectbox(key='grid.kind').select('voronoi').run()
    _refine_on(at)
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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

    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
                           default_timeout=180)
    at.run()
    said = ' '.join(w.value for w in at.warning) + \
        ' '.join(c.value for c in at.caption) + \
        ' '.join(i.value for i in at.info)
    assert 'rasters' in said, 'the panel says nothing about the rasters'
    # and it is drawn BEFORE the build, not after it
    assert 'mkgrid' in {b.key for b in at.button}


def test_the_grid_stands_on_the_runs_rectangle_whatever_the_override():
    """The panel used to build on the polygon's own rectangle, 50 m west of
    the rasters' on La Mata, and offered a button to pin a structured grid
    back onto them. It builds on the rectangle a RUN uses now
    (marmites_meshes.run_rectangle), so there is nothing to pin -- and an
    override that is deliberately WRONG is said to be unused.
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

    tmp = _scratch_config('_pintest.toml')
    try:
        # the override is refused on a mesh, so the grid is structured
        _force(tmp, 'grid', 'kind', 'kind = "structured"')
        _force(tmp, 'grid.override', 'enable', 'enable = true')
        _force(tmp, 'grid.override', 'xllcorner', 'xllcorner = 111111.0')
        at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
                               default_timeout=180)
        at.session_state['config_file'] = '_pintest.toml'
        at.run()
        assert not at.exception, [str(e.value) for e in at.exception]
        assert not [b for b in at.button if b.key == 'adoptrect'], \
            'a pin is offered for a rectangle the grid already stands on'
        warned = ' '.join(w.value for w in at.warning)
        assert 'not used' in warned, warned
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_pin_is_not_offered_on_a_mesh():
    """The override reproduces a legacy DIS grid and means nothing on an
    unstructured mesh, so there the answer is the model panel, not a button
    that would be refused by validate()."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
                           default_timeout=180)
    at.run()
    labels = [x.label for x in at.text_input]
    assert 'Folder with GIS information (defined in panel Home)' in labels, \
        labels


def _tiny_mesh_cache(ws_root, kind, rings):
    """A 2 x 2 mesh over the ponds, cached where the Grid page looks for the
    run's mesh -- so the map is drawn without a mesh someone built before."""
    import importlib.util
    import numpy as np
    import marmites_meshes as mmesh
    spec = importlib.util.spec_from_file_location(
        '_loaders_for_tests', os.path.join(CODE, 'app', 'lib', 'loaders.py'))
    loaders = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaders)
    xy = np.array([p for r in rings for p in r], dtype=float)
    x0, y0 = xy.min(axis=0) - 50.0
    x1, y1 = xy.max(axis=0) + 50.0
    xs, ys = np.linspace(x0, x1, 3), np.linspace(y0, y1, 3)
    vertices = [[i * 3 + j, float(xs[j]), float(ys[i])]
                for i in range(3) for j in range(3)]
    cell2d = []
    for i in range(2):
        for j in range(2):
            v = [i * 3 + j, i * 3 + j + 1, (i + 1) * 3 + j + 1, (i + 1) * 3 + j]
            cell2d.append([len(cell2d), float(xs[j:j + 2].mean()),
                           float(ys[i:i + 2].mean()), 4] + v)
    gp = {'vertices': vertices, 'cell2d': cell2d, 'ncpl': 4, 'nvert': 9}
    cache = loaders.mesh_cache_paths(ws_root, kind)[0].parent
    mmesh._save_cached(str(cache), kind, 'test-mesh', gp)


def test_the_map_can_draw_the_pond_footprints(tmp_path, monkeypatch):
    """A pond is a polygon the mesh has to cover, and a dot says nothing
    about whether it does -- so the overlay reads the GeoJSON outlines.

    The map needs a mesh: the test caches a tiny one in a workspace of its
    own (2026-10-09) -- it used to need the run's, which a test workspace
    does not have, and failed for that reason alone."""
    import marmites_config as mcfg
    import marmites_meshes as mmesh
    import mm_paths

    cfg = mcfg.load_run_config(os.path.join(CODE, 'configs', 'lamata.toml'))
    ds = str(mm_paths.dataset_dir(cfg.paths.case))
    rings = mmesh.pond_rings(ds)
    if not rings:
        pytest.skip('no pond table on this machine')
    assert all(len(r) >= 3 for r in rings), 'a footprint came back as a point'
    monkeypatch.setattr(mm_paths, 'WS_ROOT', tmp_path)
    _tiny_mesh_cache(str(tmp_path), cfg.grid_kind, rings)

    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    at = AppTest.from_file(os.path.join(APP, SURF),
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
    two = AppTest.from_file(os.path.join(APP, SURF),
                            default_timeout=180)
    two.run()
    assert 'surface.plot' not in _keys(two), \
        'the MMsurf figures are still on panel 2'

    four = AppTest.from_file(os.path.join(APP, PLOT),
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
    at = AppTest.from_file(os.path.join(APP, SURF),
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
    at = AppTest.from_file(os.path.join(APP, SURF),
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
    at = AppTest.from_file(os.path.join(APP, SURF),
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
    # The fractional source values live on the SUBSURFACE panel now -- the
    # channel width, the streambed -- so that is where the rule is checked
    # not to have reached them.
    at = AppTest.from_file(os.path.join(APP, SUB), default_timeout=180)
    at.run()
    vals = [n.value for n in at.number_input
            if n.key and n.key.endswith('.__v')]
    assert any(isinstance(v, float) and v != int(v) for v in vals), \
        'the subsurface panel lost its fractional source values'


def test_every_producer_that_names_a_file_offers_the_dialog():
    """By rule, not by listing the fields one at a time: a layer is a
    shapefile in the cartography folder and a raster is a file in the
    dataset, so both are chosen rather than typed correctly. A column names
    an attribute and a value is a number, so neither is."""
    at = AppTest.from_file(os.path.join(APP, SOIL),
                           default_timeout=180)
    at.run()
    keys = {b.key for b in at.button if b.key}
    for dotted, want in (('soil.zones', 'layer'), ('soil.thickness', 'raster')):
        [s for s in at.selectbox
         if s.key == dotted + '.__producer'][0].set_value(want).run()
        keys = {b.key for b in at.button if b.key}
        assert dotted + '.__v.__pick' in keys, (
            '%s as a %s has no ... button' % (dotted, want))
    # ... and a value is not a file
    [s for s in at.selectbox
     if s.key == 'soil.thickness.__producer'][0].set_value('value').run()
    assert 'soil.thickness.__v.__pick' not in {b.key for b in at.button}, (
        'a number is being chosen from the filesystem')

    # A column is not a file either. Asked of a ParamSource, which lives on
    # the subsurface panel, since a VectorSource no longer offers column as a
    # producer at all -- there it is an attribute of the layer, drawn under
    # it.
    sub = AppTest.from_file(os.path.join(APP, SUB), default_timeout=180)
    sub.run()
    [s for s in sub.selectbox
     if s.key == 'sfr.width.__producer'][0].set_value('column').run()
    assert 'sfr.width.__v.__pick' not in {b.key for b in sub.button}, (
        'a column name is being chosen from the filesystem')


def test_the_soil_column_is_edited_as_tables_not_named_as_a_file():
    """It was a path to inputSOILparam.txt, which the run never read -- the
    driver hard-coded the same path. The column is the two tables now, on
    the Soil column tab, with an import for a case that still has a file."""
    at = AppTest.from_file(os.path.join(APP, SOIL), default_timeout=180)
    at.run()
    assert not at.exception, [str(e.value) for e in at.exception]
    assert not any(x.key == 'soil.params' for x in at.text_input), \
        'the soil column is still asked as a file'
    assert 'soil_import' in {b.key for b in at.button if b.key}, \
        'no way to import an old inputSOILparam.txt'
    src = io.open(os.path.join(APP, SOIL), encoding='utf-8').read()
    tab = src[src.index('with tab_soil:'):src.index('with tab_veg:')]
    assert "'soil.zone', 'soil.horizon'" in tab, \
        'the soil column tables are not on the Soil column tab'
    veg = src[src.index('with tab_veg:'):src.index('panelui.remember(edited)')]
    assert "dotted == 'soil.veg_class'" in veg, \
        'the vegetation tab draws tables that are not its own'


def test_a_file_below_the_folder_keeps_its_relative_path():
    """Dataset files live below the dataset folder (the layer rasters in
    MF_ws/, for one). Storing an absolute path for them would tie the
    configuration to this machine."""
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
    at = AppTest.from_file(os.path.join(APP, SURF),
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
    at = AppTest.from_file(os.path.join(APP, SURF),
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
    at = AppTest.from_file(os.path.join(APP, SOIL),
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
    at = AppTest.from_file(os.path.join(APP, SUB), default_timeout=180)
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
    at = AppTest.from_file(os.path.join(APP, SOIL),
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
    at = AppTest.from_file(os.path.join(APP, SOIL),
                           default_timeout=180)
    at.run()
    [x for x in at.text_input
     if x.key == 'soil.zones.__v'][0].set_value('not_here_at_all.shp').run()
    assert not [s for s in at.selectbox if s.key == 'soil.zones.__col'], (
        'a list is offered for a layer that is not there')
    assert [x for x in at.text_input if x.key == 'soil.zones.__col'], (
        'the column was taken away instead of falling back to a box')


def test_the_overlay_rule_sits_beside_the_column():
    """Both are properties OF the layer named above them, and both are read
    off that layer's own header."""
    at = AppTest.from_file(os.path.join(APP, SOIL),
                           default_timeout=180)
    at.run()
    how = [s for s in at.selectbox if s.key == 'soil.zones.__how']
    assert how, 'the overlay rule is not offered'
    opts = list(how[0].options)
    assert opts[0] == 'auto', 'auto is not the first choice'
    # the zones are a polygon layer, so the polygon rules and not the line
    # ones -- offering `length` for a polygon is offering nothing usable
    assert 'majority' in opts and 'length' not in opts, opts


def test_the_vegetation_class_column_is_a_list_too():
    """It is a plain string field rather than part of a source, but it is
    the same question -- an attribute OF the layer beside it."""
    at = AppTest.from_file(os.path.join(APP, SOIL),
                           default_timeout=180)
    at.run()
    col = [s for s in at.selectbox if s.key == 'soil.veg_column']
    assert col, 'the vegetation class column is still typed from memory'
    assert col[0].value in list(col[0].options)
    assert len(col[0].options) > 2, list(col[0].options)

    # and it follows the LAYER: point that somewhere else and the list goes
    [x for x in at.text_input
     if x.key == 'soil.veg_layer'][0].set_value('not_here.shp').run()
    assert not [s for s in at.selectbox if s.key == 'soil.veg_column'], (
        'a list is offered for a layer that is not there')
    assert [x for x in at.text_input if x.key == 'soil.veg_column'], (
        'the column was taken away instead of falling back to a box')


def test_no_panel_carries_a_save_button_of_its_own():
    """ONE SAVE, IN THE SIDEBAR, FOR THE WHOLE CONFIGURATION.

    There used to be a *Validate & save* on every panel, which meant a panel
    filled and left unsaved was silently discarded, and "is this saved?" was
    a question with six answers. A panel now REMEMBERS and the sidebar
    writes.
    """
    for page in (SURF, SOIL, SUB, OBS, PLOT, 'pages/1_1_-_Grid.py'):
        src = io.open(os.path.join(APP, page), encoding='utf-8').read()
        assert 'panelui.save_button(' not in src, \
            '%s still draws its own save button' % page
        assert 'switch_and_save' not in src, \
            '%s still reserves a slot for one' % page
    src = io.open(os.path.join(APP, 'lib', 'panelui.py'),
                  encoding='utf-8').read()
    assert 'def save_button(' not in src, \
        'panelui still offers a per-panel save button'


def test_every_panel_remembers_last_and_offers_the_sidebar_save():
    """Collected at the END, so every sub-panel has had its say.

    A panel collects tab by tab as the tabs are drawn, so remembering where
    the old button was drawn -- at the top -- would remember an empty dict,
    and remembering inside a tab would miss every tab after it (on the
    driving-forces panel, the whole time discretisation).
    """
    for page in (SURF, SOIL, SUB, OBS, PLOT):
        src = io.open(os.path.join(APP, page), encoding='utf-8').read()
        body = [ln for ln in src.splitlines()
                if ln.strip() and not ln.lstrip().startswith('#')]
        assert body[-1] == 'panelui.sidebar_save(cfg, path)', (
            '%s: the sidebar save is not the last thing on the page: %r'
            % (page, body[-1]))
        assert body[-2] == 'panelui.remember(edited)', (
            '%s: the edits are not remembered just before it: %r'
            % (page, body[-2]))
        assert 'panelui.panel_switch(cfg, panel)' in src, (
            '%s: the master switch is gone' % page)


def test_the_time_discretisation_has_a_sub_panel_of_its_own():
    at = AppTest.from_file(os.path.join(APP, SURF), default_timeout=180)
    at.run()
    keys = {w.key for w in (list(at.checkbox) + list(at.number_input))
            if w.key}
    for dotted in ('run.daily', 'run.perlen_max', 'run.nsp'):
        assert dotted in keys, '%s is not on the panel' % dotted
    said = ' '.join(s.value for s in at.success)
    assert 'day(s) in the record' in said, said[:200]


def _scratch_config(tmp_name='_savetest.toml'):
    """A copy of the reference configuration the app can be pointed at.

    In code/configs, because that is the only folder pick_config lists --
    the app edits the file in place, so a test must not aim it at the real
    one.
    """
    import shutil

    ref = os.path.join(CODE, 'configs', 'lamata.toml')
    tmp = os.path.join(CODE, 'configs', tmp_name)
    shutil.copyfile(ref, tmp)
    return tmp


def test_validate_and_save_actually_writes():
    """It did not. The button was keyed on id(edited) -- the address of a
    dict rebuilt every run -- so it was a DIFFERENT widget on each rerun:
    the click arrived for a key that no longer existed and the new button
    read False. It appeared to work whenever CPython handed the new dict the
    address the old one had just freed, which is most of the time and not
    all of it, and the symptom was a save that silently did nothing.

    So: change something, press the button, and read the FILE back.
    """
    tmp = _scratch_config()
    try:
        def says(key, section='run'):
            # SECTION-AWARE. `model` is a key in two of them -- run.model is
            # the master switch and meta.model is what the model is called
            # -- so a scan of the whole file returns whichever comes first.
            here = None
            for line in io.open(tmp, encoding='utf-8'):
                s = line.strip()
                if s.startswith('[') and s.endswith(']'):
                    here = s[1:-1]
                elif here == section and s.startswith(key + ' '):
                    return s
            return '(missing)'

        assert says('model') == 'model = true', says('model')

        at = AppTest.from_file(os.path.join(APP, SOIL), default_timeout=180)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        [t for t in at.toggle if t.key == 'sw_run.model'][0] \
            .set_value(False).run()
        save = [b for b in at.button if b.key == 'sidebar_save']
        assert save, 'the sidebar save is not keyed predictably'
        save[0].click().run()
        assert says('model') == 'model = false', (
            'Validate & save did not write the switch: %s' % says('model'))
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_save_button_key_does_not_move():
    """A widget key that changes between runs is a widget that cannot be
    clicked. Checked as a rule, not only through the one panel."""
    src = io.open(os.path.join(APP, 'lib', 'panelui.py'),
                  encoding='utf-8').read()
    for line in src.splitlines():
        code = line.split('#')[0]
        if 'key=' in code and ('button(' in code or 'toggle(' in code):
            assert 'id(' not in code, 'unstable widget key: %r' % code.strip()


def test_saving_a_source_does_not_invent_a_value():
    """Clearing the producers wrote value = 0.0 -- which is a legitimate zone
    NUMBER, not an absence -- and the field's default is None, so the editor
    had no type to coerce to and stored the string "0.0". The file then held
    a string where a float belongs, and the count box choked on it.

    Save a panel that carries a layer-source and read the file back.
    """

    tmp = _scratch_config('_sourcetest.toml')
    try:
        at = AppTest.from_file(os.path.join(APP, SURF), default_timeout=180)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        # SOMETHING has to change, or the sidebar save is disabled and this
        # would assert on a file nothing had written.
        at.number_input(key='run.nsp').set_value(7).run()
        save = [b for b in at.button if b.key == 'sidebar_save']
        assert save and not save[0].disabled, 'the sidebar save is not offered'
        save[0].click().run()

        body = io.open(tmp, encoding='utf-8').read()
        head = body.index('[surface.irr_zones]')
        tail = body.find('\n[', head + 1)
        block = body[head:tail if tail > 0 else len(body)]
        for line in block.strip().splitlines():
            key, _, raw = line.partition('=')
            if key.strip() in ('value', 'fill'):
                assert '"' not in raw, (
                    '%s was written as a string: %s' % (key.strip(),
                                                        line.strip()))
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
                           default_timeout=180)
    at.run()
    at.selectbox(key='grid.kind').select('voronoi').run()
    _refine_on(at)
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

    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    src = open(os.path.join(APP, 'pages', '1_1_-_Grid.py'), encoding='utf-8').read()
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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


def test_the_converter_is_on_the_validation_panel_and_launch_runs_it():
    """The converter serves panels 1 to 5, and "is the dataset up to date?"
    is a validation question: its button sits next to the check that lists
    the out-of-date tables -- and a run does not depend on anyone pressing
    it, Launch converts first (lm_veg.shp was edited on 2026-09-23 and every
    run that day read the old vegetation)."""
    one = open(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
               encoding='utf-8').read()
    three = open(os.path.join(APP, SOIL), encoding='utf-8').read()
    seven = open(os.path.join(APP, 'pages', '7_7_-_Run.py'),
                 encoding='utf-8').read()
    run = seven
    assert 'conv_run' in seven and 'dataset_state.stale(cfg' in seven
    assert seven.index('conv_run') < seven.index("'Validate the configuration'"), \
        'the converter belongs above the Validate button'
    for page, name in ((one, 'Grid'), (three, 'Soil')):
        assert 'conv_run' not in page and 'tab_gis' not in page,             'the %s panel has a converter again' % name
    # Create grid still converts the two tables a grid depends on ...
    assert 'only=dataset_state.GRID_TABLES' in one
    # ... and Launch the rest, refusing to run on a failed conversion
    assert 'dataset_state.stale(cfg' in run
    assert 'dataset_state.run_converter(' in run
    assert 'so nothing was launched' in run


def test_the_grid_tabs_are_named_consistently():
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
                           default_timeout=180)
    at.run()
    src = open(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
               encoding='utf-8').read()
    assert "'Visualize grid'" in src
    assert 'Visualize mesh' not in src, 'a stale name is still referred to'


def test_panel_one_asks_for_all_three_layers():
    """The working flow is from scratch: the catchment, the hydrography and
    the ponds are all asked for HERE, at the same level, because the
    refinement options below cannot be answered before it is known whether
    there is a network or a pond to refine around."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
                           default_timeout=180)
    at.run()
    keys = _keys(at)
    for k in ('grid.boundary', 'grid.streams', 'grid.ponds'):
        assert k in keys, '%s is not asked for on panel 1' % k
    assert 'gis_folder' in keys, 'no folder to look in'


def test_the_optional_layers_can_be_set_to_none():
    """A catchment with no mapped network has to be able to say so."""
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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
    at = AppTest.from_file(os.path.join(APP, 'pages', '1_1_-_Grid.py'),
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


# ------------------------------------------------- an unplugged switch stays
# WHAT THE PANELS SHOW MUST BE WHAT RUNS. A master switch is a widget: it
# changes its panel at once and the FILE only when that panel is saved, and
# the Run page launches the FILE. Unplugging MMsurf and launching ran MMsurf,
# because [run] surface was still true. The save works -- these pin down that
# the gap between the two is now visible and blocks the launch.

SWITCH_CFG = '_switchtest.toml'


def _scratch_with_mmsurf_on():
    """A scratch configuration that certainly has MMsurf switched ON.

    NOT just a copy of the reference: the reference is a live working file
    and the first thing anyone does with the fix these tests cover is unplug
    MMsurf in it and save. A test that assumed `surface = true` there passed
    until exactly that happened -- so the starting state is SET, not found.
    """
    tmp = _scratch_config(SWITCH_CFG)
    _force(tmp, 'run', 'surface', 'surface = true')
    return tmp


def _force(path, section, key, line):
    """Put a key in a KNOWN state in a scratch configuration.

    The reference file these copies come from is a LIVE working file: it is
    edited in the browser between test runs, and a test that reads its
    starting state out of it is a test that passes until someone changes
    that setting -- which has now happened twice, once for the MMsurf switch
    and once for paths.libmf6.
    """
    out, here = [], None
    for raw in io.open(path, encoding='utf-8'):
        s = raw.strip()
        if s.startswith('[') and s.endswith(']'):
            here = s[1:-1]
        elif here == section and s.startswith(key + ' '):
            raw = line + '\n'
        out.append(raw)
    io.open(path, 'w', encoding='utf-8', newline='').write(''.join(out))
    assert _says(path, section, key) == line, _says(path, section, key)
    return path


def _says(path, section, key):
    here = None
    for line in io.open(path, encoding='utf-8'):
        s = line.strip()
        if s.startswith('[') and s.endswith(']'):
            here = s[1:-1]
        elif here == section and s.startswith(key + ' '):
            return s
    return '(missing)'


def test_unplugging_mmsurf_and_saving_writes_it():
    """The report was 'I unplug the button and it runs anyway'."""
    tmp = _scratch_with_mmsurf_on()
    try:
        assert _says(tmp, 'run', 'surface') == 'surface = true'
        at = AppTest.from_file(os.path.join(APP, SURF), default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        sw = [t for t in at.toggle if t.key == 'sw_run.surface']
        assert sw, 'panel 2 has no MMsurf switch'
        sw[0].set_value(False).run()
        save = [b for b in at.button if b.key == 'sidebar_save']
        assert save, 'the sidebar save is not keyed predictably'
        save[0].click().run()
        assert _says(tmp, 'run', 'surface') == 'surface = false', (
            'Validate & save did not unplug MMsurf: %s'
            % _says(tmp, 'run', 'surface'))
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_an_unsaved_switch_is_not_silent():
    """Turning it off and NOT saving must say so, not pass for done."""
    tmp = _scratch_with_mmsurf_on()
    try:
        at = AppTest.from_file(os.path.join(APP, SURF), default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        [t for t in at.toggle if t.key == 'sw_run.surface'][0] \
            .set_value(False).run()
        assert _says(tmp, 'run', 'surface') == 'surface = true', \
            'the toggle wrote the file without a save'
        said = ' '.join(str(w.value) for w in at.warning)
        assert 'Not saved' in said, \
            'an unsaved switch says nothing loud: %r' % said
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_run_page_refuses_a_configuration_the_panels_contradict():
    """The launch is the last moment the two can be reconciled."""
    import sys
    if APP not in sys.path:
        sys.path.insert(0, APP)
    from lib import panelui
    tmp = _scratch_with_mmsurf_on()
    try:
        at = AppTest.from_file(os.path.join(APP, 'pages', '7_7_-_Run.py'),
                               default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        # what panel 2 would have left behind, having been unplugged and
        # not saved
        at.session_state['live_run.surface'] = False
        # ... for THIS file, as pick_config records it: a toggle carried
        # over from another file is forgotten, not reported
        at.session_state['__switch_file'] = os.path.basename(tmp)
        at.run()
        launch = [b for b in at.button if b.label == 'Launch']
        assert launch, 'the Run page has no Launch button'
        assert launch[0].disabled, \
            'Launch is offered although the panels contradict the file'
        said = ' '.join(str(e.value) for e in at.error)
        assert 'error' in said.lower(), \
            'the refusal says nothing: %r' % said
        assert 'run.surface' in panelui.SWITCHES
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_switching_configuration_forgets_the_other_ones_switches():
    """A switch belongs to a file; session_state outlives the selectbox."""
    # The file we OPEN has MMsurf on; the session carries an off that
    # belonged to a different file. Opening this one must show what it says.
    tmp = _scratch_with_mmsurf_on()
    try:
        at = AppTest.from_file(os.path.join(APP, SURF), default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.session_state['live_run.surface'] = False
        at.session_state['__switch_file'] = '_some_other_case.toml'
        at.run()
        assert at.session_state['live_run.surface'] is True, \
            'the previous file\'s switch was carried over'
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def _runs_elsewhere(tmp, tmp_path):
    """Point a scratch configuration's run registry at an empty folder, so a
    run going on the real workspace cannot turn Launch off under the test."""
    runs = str(tmp_path / 'runs').replace(chr(92), '/')
    _force(tmp, 'ui', 'runs_dir', 'runs_dir = "%s"' % runs)
    return tmp


def _launch(at):
    b = [x for x in at.button if x.key == 'launch']
    assert b, 'the Run panel has no Launch button'
    return b[0]


def test_the_run_tab_is_frozen_until_the_configuration_is_validated(tmp_path):
    """The model runs only from Launch on the run tab, and only once the
    validation tab has approved the configuration."""
    tmp = _runs_elsewhere(_scratch_config('_frozentest.toml'), tmp_path)
    try:
        at = AppTest.from_file(os.path.join(APP, 'pages', '7_7_-_Run.py'),
                               default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        assert not at.exception, [str(e.value) for e in at.exception]
        assert _launch(at).disabled, 'Launch is live before validation'
        said = ' '.join(str(i.value) for i in at.info)
        assert 'Frozen' in said and 'not been validated' in said, said
        # there is NO other way to start a run on this page
        assert not [b for b in at.button if 'launch' in str(b.label).lower()
                    and b.key != 'launch'], [b.label for b in at.button]
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_validate_opens_the_run_tab_and_any_edit_freezes_it_again(tmp_path):
    """Validate saves what the panels hold and approves it; Launch is live
    for exactly that, and the next edit anywhere freezes it again. Launch is
    NEVER pressed here -- it would start a model."""
    tmp = _runs_elsewhere(_scratch_config('_validatetest.toml'), tmp_path)
    _force(tmp, 'run', 'allow_bad_budget', 'allow_bad_budget = false')
    try:
        at = AppTest.from_file(os.path.join(APP, 'pages', '7_7_-_Run.py'),
                               default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        # an edit, not saved: Validate must save it before approving
        at.checkbox(key='run.allow_bad_budget').check().run()
        [b for b in at.button if b.key == 'validate'][0].click().run()
        assert not at.exception, [str(e.value) for e in at.exception]
        assert _says(tmp, 'run', 'allow_bad_budget') == \
            'allow_bad_budget = true', 'Validate did not save the pending edit'
        errors = [e for e in at.error if 'cannot run' in str(e.value)]
        if errors:
            pytest.skip('the scratch configuration has errors here: %s'
                        % errors[0].value)
        assert not _launch(at).disabled, 'Launch stays frozen after Validate'
        ok = ' '.join(str(s.value) for s in at.success)
        assert 'Validated' in ok, ok
        # ... and the next edit freezes it again
        at.checkbox(key='run.allow_bad_budget').uncheck().run()
        assert _launch(at).disabled, 'an edit after validation left it live'
        said = ' '.join(str(i.value) for i in at.info)
        assert 'Frozen' in said and 'since it was validated' in said, said
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_library_can_be_set_from_the_run_page():
    """It was described in the schema and drawn by no panel, so the one
    field between a build and a coupled run needed a text editor."""
    tmp = _force(_scratch_config('_libtest.toml'), 'paths', 'libmf6',
                 'libmf6 = ""')
    try:
        at = AppTest.from_file(os.path.join(APP, 'pages', '7_7_-_Run.py'),
                               default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        assert not at.exception, [str(e.value) for e in at.exception]
        box = [w for w in at.text_input if w.key == 'paths.libmf6']
        assert box, 'paths.libmf6 is on no panel: %s' % [
            w.key for w in at.text_input]
        box[0].set_value('auto').run()
        save = [b for b in at.button if b.key == 'sidebar_save']
        assert save, 'the Run page does not offer the sidebar save'
        save[0].click().run()
        assert _says(tmp, 'paths', 'libmf6') == 'libmf6 = "auto"', \
            _says(tmp, 'paths', 'libmf6')
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


# ------------------------------------------- one save, for every panel at once
# THE POINT OF THE SIDEBAR SAVE. A panel used to carry its own Validate &
# save, so a panel filled and left was silently discarded and "is this saved?"
# had six answers. Edits now accumulate across pages and one button writes
# them.

def test_edits_made_on_one_panel_are_saved_from_another():
    tmp = _force(_scratch_config('_crosstest.toml'), 'run', 'nsp', 'nsp = 0')
    try:
        at = AppTest.from_file(os.path.join(APP, SURF), default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        at.number_input(key='run.nsp').set_value(42).run()
        # NOT saved by leaving the panel: carried.
        assert _says(tmp, 'run', 'nsp') == 'nsp = 0', _says(tmp, 'run', 'nsp')
        carried = at.session_state['__pending_edits']
        assert carried.get('run.nsp') == 42, carried

        # ... and written by the sidebar save on a DIFFERENT panel.
        at2 = AppTest.from_file(os.path.join(APP, PLOT), default_timeout=300)
        at2.session_state['config_file'] = os.path.basename(tmp)
        at2.session_state['__switch_file'] = os.path.basename(tmp)
        at2.session_state['__pending_edits'] = dict(carried)
        at2.run()
        save = [b for b in at2.button if b.key == 'sidebar_save']
        assert save and not save[0].disabled, \
            'the sidebar save does not offer another panel\'s edits'
        save[0].click().run()
        assert _says(tmp, 'run', 'nsp') == 'nsp = 42', _says(tmp, 'run', 'nsp')
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_remembered_edits_are_dropped_on_changing_configuration():
    """They were entered against a different file."""
    tmp = _scratch_config('_crosstest2.toml')
    try:
        at = AppTest.from_file(os.path.join(APP, PLOT), default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.session_state['__switch_file'] = '_some_other_case.toml'
        at.session_state['__pending_edits'] = {'run.nsp': 999}
        at.run()
        left = at.session_state['__pending_edits']
        assert left.get('run.nsp') != 999, \
            'another file\'s edits survived the change: %r' % left
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_validation_panel_lists_what_it_finds():
    tmp = _force(_scratch_config('_valtest.toml'), 'run', 'nsp', 'nsp = 365')
    try:
        at = AppTest.from_file(
            os.path.join(APP, 'pages', '7_7_-_Run.py'),
            default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        assert not at.exception, [str(e.value) for e in at.exception]
        said = ' '.join(str(m.value) for m in at.markdown)
        assert 'run.nsp' in said, 'the cap is not listed: %r' % said[:400]
        warned = ' '.join(str(w.value) for w in at.warning)
        assert 'warning' in warned.lower(), warned[:200]
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_launch_is_blocked_while_the_configuration_has_an_error():
    """Errors are a wall, and the button says so before it is pressed."""
    tmp = _scratch_with_mmsurf_on()
    try:
        at = AppTest.from_file(os.path.join(APP, 'pages', '7_7_-_Run.py'),
                               default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.session_state['live_run.surface'] = False     # panel says off
        at.session_state['__switch_file'] = os.path.basename(tmp)
        at.run()
        launch = [b for b in at.button if b.label == 'Launch']
        assert launch and launch[0].disabled, \
            'Launch is offered on an invalid configuration'
        # ... and it cannot be validated either
        [b for b in at.button if b.key == 'validate'][0].click().run()
        assert [b for b in at.button if b.label == 'Launch'][0].disabled
        assert any('not validated' in str(e.value) for e in at.error), \
            [str(e.value) for e in at.error]
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_a_check_does_not_print_its_key_twice():
    """Most checks name the key they are about as their first words, so
    printing the key underneath the sentence said the same thing twice:

        spinup.steady_means = 'hi_spinup' has no scope sidecar
        spinup.steady_means
    """
    src = io.open(
        os.path.join(APP, 'pages', '7_7_-_Run.py'),
        encoding='utf-8').read()
    assert 'if c.key and c.key not in c.title:' in src, \
        'the key is added to the heading unconditionally'


def test_the_coupling_is_asked_on_the_run_panel_and_saved():
    """How the run executes was TOML-only: the audit behind the cookbook's
    Appendix B found it on no panel at all. The lagged/iterative choice it
    was about is gone (2026-10-07) -- the coupling is said, not asked -- and
    the rest stays on the Run panel."""
    tmp = _force(_scratch_config('_couplingtest.toml'), 'run',
                 'allow_bad_budget', 'allow_bad_budget = false')
    try:
        at = AppTest.from_file(os.path.join(APP, 'pages', '7_7_-_Run.py'),
                               default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        assert not at.exception, [str(e.value) for e in at.exception]
        keys = {w.key for w in list(at.selectbox) + list(at.number_input)
                + list(at.checkbox) if w.key}
        for k in ('run.ats', 'run.ats_dtmin', 'run.max_discrepancy',
                  'run.allow_bad_budget', 'run.build_only'):
            assert k in keys, '%s is not on the Run panel' % k
        for gone in ('run.mode', 'run.relax'):
            assert gone not in keys, '%s is back on the Run panel' % gone
        said = ' '.join(str(c.value) for c in at.caption)
        assert 'MMsoil runs first' in said, 'the coupling is not said'
        at.checkbox(key='run.allow_bad_budget').check().run()
        save = [b for b in at.button if b.key == 'sidebar_save']
        assert save and not save[0].disabled
        save[0].click().run()
        assert _says(tmp, 'run', 'allow_bad_budget') == \
            'allow_bad_budget = true', _says(tmp, 'run', 'allow_bad_budget')
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_solver_is_asked_on_the_run_panel_and_saved():
    tmp = _force(_scratch_config('_solvertest.toml'), 'solver',
                 'outer_dvclose', 'outer_dvclose = 0.001')
    try:
        at = AppTest.from_file(os.path.join(APP, 'pages', '7_7_-_Run.py'),
                               default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        assert not at.exception, [str(e.value) for e in at.exception]
        keys = {w.key for w in list(at.selectbox) + list(at.number_input)
                if w.key}
        for k in ('solver.complexity', 'solver.outer_dvclose',
                  'solver.outer_maximum', 'solver.inner_dvclose',
                  'solver.inner_rclose', 'solver.cell_averaging'):
            assert k in keys, '%s is not on the Run panel' % k
        at.number_input(key='solver.outer_dvclose').set_value(0.002).run()
        save = [b for b in at.button if b.key == 'sidebar_save']
        assert save and not save[0].disabled
        save[0].click().run()
        assert _says(tmp, 'solver', 'outer_dvclose') == \
            'outer_dvclose = 0.002', _says(tmp, 'solver', 'outer_dvclose')
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_run_tag_is_saved_with_the_configuration():
    """The tag was a widget on the run tab that nothing wrote (2026-10-06):
    typed for one run, back to the file's value at the next. It is asked on
    the Validation tab and saved as meta.name, like every other answer."""
    tmp = _force(_scratch_config('_runtagtest.toml'), 'meta', 'name',
                 'name = "old_tag"')
    try:
        at = AppTest.from_file(os.path.join(APP, 'pages', '7_7_-_Run.py'),
                               default_timeout=300)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        assert not at.exception, [str(e.value) for e in at.exception]
        assert not [w for w in at.text_input if w.key == 'run_tag'], \
            'an unsaved run-tag box is still on the run tab'
        box = [w for w in at.text_input if w.key == 'meta.name']
        assert box, 'the run tag is not asked on the Run panel'
        box[0].set_value('evt_fix').run()
        save = [b for b in at.button if b.key == 'sidebar_save']
        assert save and not save[0].disabled
        save[0].click().run()
        assert _says(tmp, 'meta', 'name') == 'name = "evt_fix"', \
            _says(tmp, 'meta', 'name')
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_plots_panel_asks_for_the_input_maps_first():
    """One "Input maps" switch, on the top; the rest is the OUTPUT, and
    says so -- "Post-process" said how, not what."""
    at = AppTest.from_file(os.path.join(APP, PLOT), default_timeout=300)
    at.run()
    assert not at.exception, [str(e.value) for e in at.exception]
    keys = [c.key for c in at.checkbox if str(c.key).startswith('postproc.')]
    assert keys[0] == 'postproc.input_maps', keys
    labels = [c.label for c in at.checkbox]
    assert sum('Input maps' in lb for lb in labels) == 1, labels
    assert any('Output maps and plots' in lb for lb in labels), labels
    assert not any('Post-process' in lb for lb in labels), labels
