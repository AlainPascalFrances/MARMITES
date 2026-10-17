# -*- coding: utf-8 -*-
"""Panel 9 -- PEST++ PROGRESS: follow a calibration while it runs, and see
what each iteration gave.  WP7 front-end, phase 1
(docs/MM-MF6/WP7_PEST_frontend_plan.md).

READ-ONLY. Everything shown is read from the PEST workspace's master folder
by lib/pest_progress.py; nothing here starts, stops or changes a run (that
is panel 8's, phase 4). The PEST workspace is set on Home and only shown
here (user, 2026-10-17).

REFRESHED AS A FRAGMENT, like the Run panel's log, but every REFRESH_S
seconds rather than every few: one run of a calibration takes one to three
hours, and the record files grow to megabytes.
"""

import html
import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

import mm_paths                       # noqa: E402
from lib import pest_progress as pp   # noqa: E402

REFRESH_S = 30

st.set_page_config(page_title='PEST++ progress', page_icon='📈',
                   layout='wide')
st.title('📈  9 — PEST++ progress')
WS = str(mm_paths.PEST_WS) if mm_paths.PEST_WS else ''
st.caption('PEST workspace (set on Home): `%s`' % (WS or 'not set'))
if not WS:
    st.info('No PEST workspace is set yet: choose one on **Home**.')
    st.stop()

auto = st.toggle('Refresh every %d s' % REFRESH_S, value=True)

# what the report figures are, by file name
_TITLES = {'phi_groups.png': 'Phi per observation group'}


def _title(fn):
    if fn in _TITLES:
        return _TITLES[fn]
    if fn.startswith('ens_'):
        return 'Ensemble against the observations: %s' % fn[4:-4]
    return fn


def _figures(folder):
    """Two to a row, aligned on top, the title centred below each -- as on
    the Results panel."""
    pngs = sorted(f for f in os.listdir(folder) if f.endswith('.png'))
    pngs.sort(key=lambda f: (f != 'phi_groups.png', f))
    for i in range(0, len(pngs), 2):
        cols = st.columns(2, gap='medium', vertical_alignment='top')
        for col, fn in zip(cols, pngs[i:i + 2]):
            col.image(os.path.join(folder, fn), width='stretch')
            col.markdown('<div style="text-align:center">%s</div>'
                         % html.escape(_title(fn)), unsafe_allow_html=True)


@st.fragment(run_every=REFRESH_S if auto else None)
def _show():
    import pandas as pd
    p = pp.progress(WS)
    if p is None:
        st.info('No calibration in `%s` yet: the master folder appears when '
                'a run is launched.' % pp.master_dir(WS))
        return
    c = st.columns(6)
    c[0].metric('State', p['state'])
    c[1].metric('Iteration', '%s of %s' % (p['iteration'],
                                           p['noptmax'] if p['noptmax']
                                           is not None else '?'))
    c[2].metric('Stage', p['stage'])
    c[3].metric('Runs going', len(p['agents']))
    c[4].metric('Elapsed', '%.1f h' % p['elapsed_h']
                if p['elapsed_h'] is not None else '—')
    c[5].metric('Estimated end', p['eta'].strftime('%a %d %b %H:%M')
                if p['eta'] else '—')
    st.caption('%s realisation(s), %d agent(s); mean run %s; master `%s`.'
               % (p['nreal'] or '?', p['n_agents'],
                  ('%.0f min' % p['mean_run_min']) if p['mean_run_min']
                  else '—', p['master'])
               + (' The end is a rough estimate: what is left of this '
                  'batch, then the iterations to come, at the mean run '
                  'time.' if p['eta'] else ''))
    if p['state'] == 'stopped':
        st.warning('No pestpp-ies process is running and the run did not '
                   'finish: it was stopped or interrupted.')

    tabs = st.tabs(['Runs', 'Phi', 'Failures (%d)' % len(p['failures']),
                    'Reports (%d)' % len(p['reports'])])
    with tabs[0]:
        st.markdown('**Per batch** — an iteration tests its lambdas on a few '
                    'realisations, then runs the rest with the best one.')
        st.dataframe(pd.DataFrame(p['batches']), hide_index=True,
                     width='stretch',
                     column_config={
                         'mean run (min)': st.column_config.NumberColumn(
                             format='%.0f'),
                         'longest run (min)': st.column_config.NumberColumn(
                             format='%.0f')})
        if p['agents']:
            st.markdown('**Going now**')
            st.dataframe(pd.DataFrame(p['agents']), hide_index=True,
                         width='stretch',
                         column_config={'elapsed (min)':
                                        st.column_config.NumberColumn(
                                            format='%.0f')})
    with tabs[1]:
        if p['phi'] is None:
            st.info('No phi yet: the first is written when the initial '
                    'ensemble has run.')
        else:
            st.markdown('**Phi per iteration** (without the measurement '
                        'noise)')
            st.dataframe(p['phi'], hide_index=True,
                         width='stretch')
            st.scatter_chart(p['phi_long'], x='iteration', y='phi',
                             height=300)
            if p['phi_groups'] is not None:
                st.markdown('**Mean phi per observation group**')
                st.dataframe(p['phi_groups'], hide_index=True,
                             width='stretch')
        if p['rec']['lambda_tests']:
            st.markdown('**Lambda tests** (the subset of realisations; '
                        '"best" is the lambda kept)')
            st.dataframe(pd.DataFrame(p['rec']['lambda_tests']),
                         hide_index=True, width='stretch')
    with tabs[2]:
        if not p['failures']:
            st.success('No run has failed or been killed.')
        else:
            st.markdown('A failed realisation is dropped from the ensemble. '
                        'The cause is read from its agent\'s model log while '
                        'that agent has had no run since.')
            st.dataframe(pd.DataFrame(p['failures']), hide_index=True,
                         width='stretch')
    with tabs[3]:
        if not p['reports']:
            st.info('No report yet: one is written at the end of each '
                    'iteration (tools/pest_lamata_report.py).')
        for k, folder in reversed(p['reports']):
            with st.expander('Iteration %d' % k, expanded=(k == max(
                    r[0] for r in p['reports']))):
                summ = os.path.join(folder, 'summary.txt')
                if os.path.exists(summ):
                    st.code(open(summ, encoding='utf-8',
                                 errors='replace').read(), language=None)
                _figures(str(folder))


_show()
