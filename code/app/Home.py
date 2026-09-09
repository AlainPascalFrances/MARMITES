# -*- coding: utf-8 -*-
"""MARMITES / MODFLOW 6 -- Streamlit front-end.  WP1b, stage 1.

    streamlit run code/app/Home.py

A VIEWER AND A LAUNCHER. It reads the inputs, edits the configuration, starts a
run and shows the output. It is never part of the model path: `code/app/` may
import the model, the model may never import `code/app/` or streamlit, and
`code/tests/test_repo_hygiene.py` asserts it -- so the model keeps running
headless from Spyder and from a PEST worker where Streamlit is not installed.
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
from lib import runs as runlib                   # noqa: E402

st.set_page_config(page_title='MARMITES / MF6', page_icon='💧', layout='wide')


def case_selector():
    """The case study every page works against, kept in session state."""
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
    d = (cfg.ui.runs_dir if cfg and cfg.ui.runs_dir else
         os.path.join(str(mm_paths.WS_ROOT), 'runs'))
    return d


def sidebar_paths(case):
    """The machine banner: where this session is actually reading and writing."""
    st.sidebar.markdown('### Paths')
    rows = [('dataset', mm_paths.dataset_dir(case)),
            ('GIS (Tier B)', mm_paths.GIS),
            ('workspace', mm_paths.WS_ROOT),
            ('libmf6', mm_paths.LIBMF6)]
    for label, p in rows:
        ok = os.path.exists(str(p))
        st.sidebar.markdown('%s **%s**  \n`%s`'
                            % ('🟢' if ok else '🔴', label, p))
    st.sidebar.caption('Edit `code/mm_paths.py`, or set the `MM_*` environment '
                       'variables.')


def main():
    case = case_selector()
    sidebar_paths(case)

    st.title('MARMITES / MODFLOW 6')
    st.caption('Soil-water balance coupled to MODFLOW 6 — inputs, configuration, '
               'runs and results for **%s**.' % case)

    c1, c2, c3 = st.columns(3)
    ds = mm_paths.dataset_dir(case)
    n_inputs = sum(len(fs) for _, fs in
                   __import__('lib.loaders', fromlist=['x']).inventory(ds)) \
        if ds.is_dir() else 0
    with c1:
        st.metric('Tier-A inputs', n_inputs)
    with c2:
        rd = runs_dir()
        st.metric('Runs on record', len(runlib.list_runs(rd)))
    with c3:
        st.metric('Case', case)

    st.markdown("""
### Where to go

| Page | What it does |
|---|---|
| **Inputs** | every file MM and MF read, plus the source cartography on a map |
| **Configuration** | edit and validate a run configuration; the MMsurf parameters, already filled |
| **Run** | launch on this machine or the server, and follow the log |
| **Results** | the figures a run wrote, with a run picker |
| **Calibration** | PEST++-IES output *(stage 3, WP7)* |

The model is launched as a **detached** process: a run survives closing this
browser tab, and its log keeps being written. Nothing here runs a model inline —
a coupled run is minutes and an IES run is hours.
""")

    recent = runlib.list_runs(runs_dir(), limit=5)
    if recent:
        st.markdown('### Recent runs')
        for r in recent:
            state = r.get('state')
            icon = ('⏳' if state == 'running'
                    else '✅' if r.get('outcome') == 'completed'
                    else '⛔' if r.get('outcome') == 'stopped' else '❌')
            st.write('%s `%s` — %s, started %s'
                     % (icon, r['run_id'], state, r.get('started')))


if __name__ == '__main__':
    main()
else:
    main()
