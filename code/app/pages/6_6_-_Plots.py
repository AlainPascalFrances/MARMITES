# -*- coding: utf-8 -*-
"""Panel 6 -- PLOTS.  WP1d.

What to draw once the run finishes. Nothing here changes a flux, which is
what makes it the last panel -- and what makes its master switch cheap to
turn off while a configuration is being tried out.

Laid out by what the figures ARE (user, 2026-10-11: the generic form put
checkboxes and numbers side by side in field order): the maps, then the
time series and the water balance, then MMsurf's own figures; then the
time axis every series shares; the notes last. An option that only refines
another (the map days, the Sankey's variants) is asked under it, and only
while it is on -- off, its value is kept, not lost.

One field here is not cosmetic and used to do nothing at all:
``hydro_year_start``. It was in the MM ini, it reached ``MMConfig``, and
nothing ever set it on the model -- every call site fell back to
``getattr(cMF, 'iniMonthHydroYear', 10)``, so October was hardcoded. It is a
real field now, and it drives the x-axis of roughly twenty time-series
figures and the Sankey year index.
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

from lib import panelui, schema             # noqa: E402

st.set_page_config(page_title='6 Plots', page_icon='📊', layout='wide')
cfg, path = panelui.pick_config()
panelui.dataset_banner(cfg)
panel = panelui.header(6)

edited = panelui.panel_switch(cfg, panel)


def ask(dotted):
    """One field, with its label, units and help; its edit is remembered.
    Returns the value the widget shows."""
    sec, leaf = dotted.split('.', 1)
    value = getattr(getattr(cfg, sec), leaf)
    got = panelui._widget(dotted, value)
    if got is not None:
        edited[dotted] = got
        return got
    return value


def under():
    """An indented column for the options that refine the switch above."""
    return st.columns([1, 14])[1]


# ---------------------------------------------------------------- draw
st.markdown('#### What to draw')
c_maps, c_wb, c_surf = st.columns(3)
with c_maps:
    st.markdown('**Maps**')
    ask('postproc.input_maps')
    if ask('postproc.result_maps'):
        with under():
            ask('postproc.map_days')
            st.caption('The head and soil-storage maps on that many days, '
                       'evenly spread over the run, on one colour ramp; 0 '
                       'draws the time means only.')
with c_wb:
    st.markdown('**Time series and water balance**')
    ask('postproc.obs_series')
    if ask('postproc.sankey'):
        with under():
            ask('postproc.sankey_full')
            ask('postproc.sankey_obs_years')
            ask('postproc.sankey_min_flux')
with c_surf:
    st.markdown('**MMsurf**')
    # an int on file (0 PNG, 1 on screen), asked in words
    _how = ['PNG files in the workspace', 'on screen']
    _cur = int(cfg.surface.plot)
    pick = st.selectbox('MMsurf figures', _how,
                    index=min(max(_cur, 0), 1), key='surface.plot',
                    help='`surface.plot`. %s' % schema.describe(
                        'surface.plot')[2])
    edited['surface.plot'] = _how.index(pick)
    st.caption('Drawn by MMsurf itself while it turns the hourly record into '
               'the daily forcing, so only on a run with **%s** switched '
               'on. Everything else here is drawn after the model.'
               % schema.panel_name(2))

# ----------------------------------------------------------- time axis
st.markdown('#### Time axis and units')
MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
          'August', 'September', 'October', 'November', 'December']
c1, c2, c3, c4 = st.columns(4)
cur = int(cfg.postproc.hydro_year_start)
pick = c1.selectbox('Hydrological year starts in', MONTHS, index=cur - 1,
                    key='postproc.hydro_year_start',
                    help='`postproc.hydro_year_start`. Drives the x-axis of '
                         'every time series, the Sankey year index and the '
                         'ponds\' yearly budgets.')
edited['postproc.hydro_year_start'] = MONTHS.index(pick) + 1
with c2:
    ask('postproc.wb_unit')
with c3:
    ask('postproc.tick_trimester_years')
with c4:
    ask('postproc.tick_semester_years')
st.caption('A run that does not span a whole hydrological year still gets a '
           'whole-period panel, scaled to mm/y and labelled as such — the '
           'figures say so rather than quietly annualising a fortnight. The '
           'ticks: a quarter per tick on a run up to the first number of '
           'years, a half-year up to the second, a year beyond.')

# --------------------------------------------------------------- notes
with st.expander('Notes: where the figures go, and the open-water '
                 'evaporation in the water balance'):
    st.markdown('Figures are written into the run folder under the '
                'workspace, and the **Results** panel shows them with a run '
                'picker. Turning this group off leaves the run itself '
                'untouched — a finished run can always be re-drawn later: '
                'turn **Run MMsoil + MODFLOW 6** off with the figures on, '
                'and launch.')
    st.info('**Open-water evaporation is drawn as a SPLIT OF RUNOFF**, not as '
            'a loss from a surface store — there is no longer such a store. '
            'The surface box therefore closes exactly:\n\n'
            '`Pe + Exf_1 = I + E_ow + Ro_net`\n\n'
            'A reach or a pond also evaporates water that reached it from '
            'upstream (and from the aquifer), so in its cell `E_ow` can '
            'exceed the runoff that cell generated — a pond cell evaporates '
            'close to `Eo` all year. The split is clamped and the run SAYS SO '
            'rather than drawing a negative flow — expect that note on the '
            'panels of points in stream or pond cells, not on the catchment '
            'one.')

# REMEMBERED, not written: the one save is in the sidebar (see panelui).
panelui.remember(edited)
panelui.sidebar_save(cfg, path)
