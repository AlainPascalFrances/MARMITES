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
