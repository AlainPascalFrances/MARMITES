# -*- coding: utf-8 -*-
"""Panel 4 -- PLOTS.  WP1d.

What to draw once the run finishes. Nothing here changes a flux, which is
what makes it the last panel -- and what makes its master switch cheap to
turn off while a configuration is being tried out.

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

from lib import panelui                     # noqa: E402

st.set_page_config(page_title='4 Plots', page_icon='📊', layout='wide')
cfg, path = panelui.pick_config()
panelui.dataset_banner(cfg)
panel = panelui.header(4)

edited = {}
edited.update(panelui.master_switch(cfg, panel[3]) or {})

MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
          'August', 'September', 'October', 'November', 'December']

st.markdown('#### The hydrological year')
c1, c2 = st.columns([1, 2])
cur = int(cfg.postproc.hydro_year_start)
pick = c1.selectbox('Starts in', MONTHS, index=cur - 1,
                    help='`postproc.hydro_year_start`. Drives the x-axis of '
                         'every time series and the Sankey year index.')
edited['postproc.hydro_year_start'] = MONTHS.index(pick) + 1
c2.caption('A run that does not span a whole hydrological year still gets a '
           'whole-period panel, scaled to mm/y and labelled as such — the '
           'figures say so rather than quietly annualising a fortnight.')

st.markdown('#### What to draw')
edited.update(panelui.section_form(cfg, 'postproc', columns=3))

st.markdown('#### The water balance')
st.info('**Open-water evaporation is drawn as a SPLIT OF RUNOFF**, not as a '
        'loss from a surface store — there is no longer such a store. The '
        'surface box therefore closes exactly:\n\n'
        '`Pe + Exf_1 = I + E_ow + Ro_net`\n\n'
        'The stream is also fed by groundwater, so over a dry window `E_ow` '
        'can exceed the runoff generated in that cell. The split is clamped '
        'and the run SAYS SO rather than drawing a negative flow — expect '
        'that note on per-point panels, not on the catchment one.')

panelui.save_button(cfg, path, edited)

st.markdown('---')
st.caption('Figures are written into the run folder under the workspace, and '
           'the **Results** panel shows them with a run picker. Turning this '
           'group off leaves the run itself untouched — a finished run can '
           'always be re-drawn later with `postproc.only`.')
