# -*- coding: utf-8 -*-
"""Panel 5 -- STATE VARIABLES: what the run is compared against.  WP1d.

The observations used to be one tab of the old Model panel, three prefixes in
a row with no heading -- which said nothing about WHAT was being compared,
and left out actual evapotranspiration entirely although the model produces
it. Here each state variable the model computes gets its own block, and each
says whether the series are actually on disk.

Nothing on this panel changes a flux. It decides what the run is measured
against, which is why it sits after the two panels that decide what the run
IS.
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP, os.path.join(CODE, 'ppMF6')):
    if p not in sys.path:
        sys.path.insert(0, p)

from lib import panelui, schema             # noqa: E402

st.set_page_config(page_title='5 State variables', page_icon='📏',
                   layout='wide')
case = st.session_state.get('case', 'LaMata')
cfg, path = panelui.pick_config()
ds = panelui.dataset_banner(cfg)
panel = panelui.header(5)

edited = {}

# ------------------------------------------------------------ the points
st.markdown('#### The observation points')
st.caption('One table for every state variable: the points are the same '
           'places, and a piezometer that is also a soil-moisture site '
           'should not be two rows that can disagree.')
edited.update(panelui.rows_form(
    cfg, (schema.OBS_COMMON,), 'obs', columns=len(schema.OBS_COMMON)))

tbl = os.path.join(str(ds), panelui.live('obs.table', cfg.obs.table) or '')
if os.path.exists(tbl):
    try:
        from marmites_postprocess import obs_points
        pts = obs_points(str(ds))
        st.success('%d active point(s)' % len(pts))
        st.dataframe([{'name': p['name'], 'x': p['x'], 'y': p['y'],
                       'layer': p.get('lay')} for p in pts],
                     width='stretch', hide_index=True)
    except Exception as exc:
        st.info('Could not read the points: %r' % exc)
        pts = []
    with st.expander('Show the file'):
        st.code(open(tbl, encoding='utf-8', errors='replace').read())
    st.caption('`##` before a name means the point is not drawn on the maps; '
               'a single `#` comments it out of the run entirely.')
else:
    st.error('Missing: %s' % tbl)
    pts = []

# --------------------------------------------------- the state variables
st.markdown('---')
st.markdown('#### The measured series')
st.caption('One block per state variable the model produces. A prefix left '
           'blank means that variable is not measured here — which is an '
           'answer, not an omission.')

tabs = st.tabs([title for title, _f, _n in schema.OBS_GROUPS])
for tab, (title, dotted, note) in zip(tabs, schema.OBS_GROUPS):
    with tab:
        st.caption(note)
        edited.update(panelui.rows_form(cfg, ((dotted, None),), 'obs',
                                        columns=2))
        prefix = panelui.live(dotted, getattr(cfg.obs,
                                              dotted.split('.')[-1], '')) or ''
        if not prefix:
            st.info('Not measured in this catchment, so nothing is compared '
                    'against %s.' % title.lower())
            continue
        # WHICH files that prefix actually finds. A prefix is only useful if
        # it names something, and "<prefix>_<point>.txt" is easy to get one
        # character wrong -- the kind of thing that shows up as an empty
        # calibration plot hours later.
        found = sorted(f for f in (os.listdir(str(ds))
                                   if os.path.isdir(str(ds)) else [])
                       if f.startswith(prefix + '_') and f.endswith('.txt'))
        if found:
            st.success('%d series: %s' % (len(found), ', '.join(found)))
            named = {f[len(prefix) + 1:-4] for f in found}
            loose = [p['name'] for p in pts if p['name'] not in named]
            if loose:
                st.caption('No series for %s — those points are not compared '
                           'against this variable.' % ', '.join(loose))
        else:
            st.warning('`%s_*.txt` matches nothing in the dataset, so this '
                       'variable has no measurements to be compared against.'
                       % prefix)

panelui.save_button(cfg, path, edited)
