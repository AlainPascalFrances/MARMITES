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


def test_panel_zero_offers_a_browser_for_each_path():
    """"Select the folder inside the computer": typed OR browsed to, without
    a native dialog -- the server process is only the user's own machine
    while the app runs locally, which ui.execution says is not a given."""
    at = AppTest.from_file(os.path.join(APP, 'Home.py'), default_timeout=180)
    at.run()
    keys = {b.key for b in at.button if b.key}
    assert 'path_gis.__up' in keys and 'path_gis.__use' in keys
    # ... and the two that are FILES, not folders, are picked as files.
    assert 'path_nwt_ref.__usef' in keys or 'path_nwt_ref.__up' in keys


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
    # ... and a coarser grade asks for fewer of them.
    n_before = len(shown(at).split(','))
    at.number_input(key='grid.voronoi.grade_ratio').set_value(2.5).run()
    assert len(shown(at).split(',')) < n_before, \
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
