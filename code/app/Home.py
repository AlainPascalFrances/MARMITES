# -*- coding: utf-8 -*-
"""Panel 0 -- MARMITES / MODFLOW 6 front-end.  WP1d.

    streamlit run code/app/Home.py

A VIEWER, AN EDITOR AND A LAUNCHER. It reads the inputs, edits the
configuration, starts a run and shows the output. It is never part of the
model path: ``code/app/`` may import the model, the model may never import
``code/app/`` or streamlit, and ``code/tests/test_repo_hygiene.py`` asserts
it -- so the model keeps running headless from Spyder and from a PEST worker
where Streamlit is not installed.

The panels are in the order the work is done, and the grid comes first
because every other input is wrapped onto it.
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

import mm_paths                                  # noqa: E402
from lib import panelui, runs as runlib          # noqa: E402

st.set_page_config(page_title='MARMITES / MF6', page_icon='💧', layout='wide')


def case_selector():
    cases = []
    ex = mm_paths.REPO / 'example'
    if ex.is_dir():
        cases = sorted(d.name for d in ex.iterdir()
                       if d.is_dir() and not d.name.startswith('.'))
    cases = cases or ['LaMata']
    default = st.session_state.get('case', 'LaMata')
    idx = cases.index(default) if default in cases else 0
    return st.sidebar.selectbox('Case study', cases, index=idx, key='case')


def runs_dir(cfg=None):
    return (cfg.ui.runs_dir if cfg and cfg.ui.runs_dir else
            os.path.join(str(mm_paths.WS_ROOT), 'runs'))


def main():
    case = case_selector()
    cfg, path = panelui.pick_config()
    panelui.dataset_banner(cfg)
    st.sidebar.caption('Edit `code/mm_paths.py`, or set the `MM_*` '
                       'environment variables.')

    st.title('💧  MARMITES / MODFLOW 6')
    st.caption('A soil water balance coupled to MODFLOW 6 through the API. '
               'Case **%s**, configuration `%s`.'
               % (case, os.path.basename(path)))

    st.markdown('### What this configuration will run')
    cols = st.columns(3)
    for c, (sw, what) in zip(cols, [
            ('run.surface', 'MMsurf — the daily forcing'),
            ('run.model', 'MMsoil + MODFLOW 6'),
            ('run.plot', 'Figures')]):
        section, key = sw.split('.')
        on = bool(getattr(getattr(cfg, section), key))
        c.metric(what, 'ON' if on else 'off')
    if not cfg.run.surface:
        st.caption('MMsurf is off, so the daily forcing must already exist. It '
                   'is checked for presence AND shape before the run starts — '
                   'a stale or truncated file stops the run rather than being '
                   'used.')

    st.markdown('### The panels, in the order to fill them')
    st.markdown("""
| Panel | What it settles | Switch |
|---|---|---|
| **1 Grid** | the catchment polygon, and the grid built inside it | always |
| **2 Surface** | the meteorological record → the daily forcing | `run.surface` |
| **3 Model** | the soil column, the aquifer, the stream and the ponds | `run.model` |
| **4 Plots** | what to draw afterwards | `run.plot` |
| **5 Run** | launch it, and follow the log | — |
| **6 Results** | the figures a run wrote | — |
| **7 Inputs** | every file the model reads, and the source cartography | — |

**The grid comes first**, because every other input is *wrapped onto it*: the
soil zones, the vegetation cover, the stream network and the observation
points are vector layers, projected onto whichever grid panel 1 produced.
Change the grid and they follow — they do not have to be re-made.

The switches are not decoration. They are the same `[run]` keys the driver
reads, so turning one off means that half of the model does not execute.
""")

    ds = mm_paths.dataset_dir(case)
    c1, c2, c3 = st.columns(3)
    c1.metric('Grid', '%s @ %g m' % (cfg.grid.kind, cfg.grid.cell_size))
    c2.metric('Vegetation types', len(cfg.surface.vegetation))
    c3.metric('Runs on record', len(runlib.list_runs(runs_dir(cfg))))

    if not ds.is_dir():
        st.error('The dataset folder does not exist: %s' % ds)

    recent = runlib.list_runs(runs_dir(cfg), limit=5)
    if recent:
        st.markdown('### Recent runs')
        for r in recent:
            state = r.get('state')
            icon = ('⏳' if state == 'running'
                    else '✅' if r.get('outcome') == 'completed'
                    else '⛔' if r.get('outcome') == 'stopped' else '❌')
            st.write('%s `%s` — %s, started %s'
                     % (icon, r['run_id'], state, r.get('started')))

    st.markdown('---')
    st.caption('A run is launched as a DETACHED process: it survives closing '
               'this tab, and its log keeps being written. Nothing here runs a '
               'model inline — a coupled run is minutes and a calibration is '
               'hours.')


main()
