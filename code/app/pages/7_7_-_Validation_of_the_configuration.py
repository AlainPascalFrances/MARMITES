# -*- coding: utf-8 -*-
"""Panel 7 -- VALIDATION OF THE CONFIGURATION.  WP1d.

Everything that can be said about this configuration before a run, in one
place, BEFORE the Run panel. The order is the point: fill the panels, save,
look here, then launch.

It answers one question -- *is this configuration ready?* -- and it answers
it in the same words the Launch button uses, because both ask ``lib.checks``.
A validation page that disagreed with the thing guarding the run would be
worse than no page at all.

  ERROR    the run cannot proceed. Launch is blocked and sends the modeller
           here.
  WARNING  the run will proceed and may not mean what was intended. Launch
           sends the modeller here to SEE it, and it can then be taken
           anyway from this page.
  INFO     worth knowing, nothing to fix.

Nothing on this page edits anything. It reads the SAVED file, because that
is what a run reads.
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

from lib import checks as chk                # noqa: E402
from lib import panelui, schema              # noqa: E402

_TITLES = dict((p[0], p[1]) for p in schema.PANELS)
_TITLES.setdefault(7, 'Validation of the configuration')
_TITLES.setdefault(8, 'Run')


def _panel_title(num):
    return _TITLES.get(num, '')

st.set_page_config(page_title='7 Validation', page_icon='🔎', layout='wide')

cfg, path = panelui.pick_config()
panelui.dataset_banner(cfg)

st.title('🔎  7 — Validation of the configuration')
st.caption('Everything that can be said about this configuration before a '
           'run. The Launch button on the Run panel asks exactly these '
           'questions and sends you back here if any of them has an '
           'answer.')

# THE ONE THING THE SIDEBAR SAVE DOES NOT COMMIT, said before the list so it
# is not mistaken for an oversight. Everything else on every panel is written
# by *Validate & save*; a grid is not, because committing one also promotes
# the mesh it produced to where the driver looks, and only the Grid panel
# knows which mesh that is.
st.info('**The grid is committed on the Grid panel, not by the sidebar '
        'save.** *Select this grid for the model* writes the `[grid]` '
        'settings AND puts the mesh they produced where a run will find '
        'it — the sidebar save cannot do the second half, so a grid built '
        'and not selected is not the model\'s grid. Which one a run would '
        'use is listed below, and again on the Run panel.')


# UNSAVED EDITS FIRST, and on their own. They are not a property of the
# configuration -- they are a property of this browser session -- so they are
# not one of the checks; but validating a file while the panels hold
# something else is exactly the trap the Run panel now refuses to walk into.
todo = panelui.unsaved_changes(cfg)
if todo:
    st.warning('**%d unsaved change(s).** What follows is the FILE, which is '
               'what a run reads — not what the panels currently show. Press '
               '*Validate & save* in the sidebar first.' % len(todo))
    with st.expander('What is unsaved'):
        for line in todo:
            st.code(line, language='ini')

found = chk.collect_default(cfg, panelui.unsaved_switches(cfg))
counts = chk.count_by_level(found)

c1, c2, c3, c4 = st.columns(4)
c1.metric('Errors', counts[chk.ERROR])
c2.metric('Warnings', counts[chk.WARNING])
c3.metric('Notes', counts[chk.INFO])
c4.metric('Hash', cfg.config_hash())

if counts[chk.ERROR]:
    st.error('**This configuration cannot run.** Fix the errors below; the '
             'Launch button stays disabled until they are gone.')
elif counts[chk.WARNING]:
    st.warning('**It will run.** Nothing here stops it, but each warning is '
               'something that has misled a run before — read them, then '
               'launch from the Run panel or from the button at the foot '
               'of this page.')
else:
    st.success('**Ready.** Nothing to report; the run can be launched.')

st.markdown('---')

ICON = {chk.ERROR: '⛔', chk.WARNING: '⚠️', chk.INFO: 'ℹ️'}
LABEL = {chk.ERROR: 'Errors', chk.WARNING: 'Warnings', chk.INFO: 'Notes'}

for level in (chk.ERROR, chk.WARNING, chk.INFO):
    here = [c for c in found if c.level == level]
    if not here:
        continue
    st.markdown('### %s %s' % (ICON[level], LABEL[level]))
    for c in here:
        where = ('panel %d — %s' % (c.panel, _panel_title(c.panel))
                 if c.panel is not None and _panel_title(c.panel) else '')
        # The key is shown ONLY when the sentence does not already carry it.
        # Most checks name the key they are about as their first words --
        # "spinup.steady_means = 'hi_spinup' has no scope sidecar" -- so
        # adding it underneath printed the same thing twice.
        head = '**%s**' % c.title
        if c.key and c.key not in c.title:
            head += '  \n`%s`' % c.key
        box = st.container(border=True)
        with box:
            st.markdown(head)
            if c.detail:
                st.caption(c.detail)
            if where:
                st.caption('Answered on %s.' % where)

st.markdown('---')

# LAUNCH FROM HERE TOO. Panel 8 sends the modeller here on a warning; making
# them navigate back to act on what they have just read would only teach
# them to skip this page.
if counts[chk.ERROR]:
    st.button('Launch the run', disabled=True,
              help='There are errors. Fix them first.')
    st.caption('Blocked by %d error(s).' % counts[chk.ERROR])
else:
    if st.button('Launch the run', type='primary'):
        st.session_state['validated'] = cfg.config_hash()
        if not panelui.go_to(panelui.RUN_PAGE):
            st.session_state.pop('validated', None)
            st.error('Could not open the Run panel — go there and press '
                     'Launch.')
    st.caption('Goes to the Run panel and starts it there, so the log is '
               'where runs are followed.')
