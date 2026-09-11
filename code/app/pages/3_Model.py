# -*- coding: utf-8 -*-
"""Panel 3 -- MODEL: the soil column and MODFLOW 6.  WP1d.

One switch, because they are no longer separable: the legacy "run MMsoil
alone" mode belonged to the Picard loop that Phase 1 removed, and MMsoil is
now stepped from inside the MODFLOW time loop through the API.

The spatial inputs are VECTOR LAYERS, wrapped onto whichever grid panel 1
produced. Each declares where its value comes from, with an explicit
precedence -- raster beats layer beats a single value -- rather than being
decided by which file happens to exist.
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP, os.path.join(CODE, 'ppMF6')):
    if p not in sys.path:
        sys.path.insert(0, p)

import mm_paths                             # noqa: E402
from lib import panelui, schema             # noqa: E402

st.set_page_config(page_title='3 Model', page_icon='🌍', layout='wide')
case = st.session_state.get('case', 'LaMata')
cfg, path = panelui.pick_config()
ds = panelui.dataset_banner(cfg)
panel = panelui.header(3)

edited = {}
edited.update(panelui.master_switch(cfg, panel[3]) or {})

tab_soil, tab_aq, tab_water, tab_obs = st.tabs(
    ['Soil column', 'Aquifer & solver', 'Streams, ponds & runoff',
     'Observations'])

# ------------------------------------------------------------------ soil
with tab_soil:
    st.markdown('#### Where each cell gets its soil from')
    edited.update(panelui.section_form(cfg, 'soil', columns=2))

    st.info('**Precedence is explicit, not decided by which file exists:** a '
            'raster beats a polygon attribute, and a polygon attribute beats '
            'a single value. Nothing set is an error, not a silent zero.')

    st.markdown('#### Vegetation classes')
    st.caption('What the layer\'s class column holds, and which vegetation '
               'type of panel 2 it means. The share of each cell covered is '
               'an exact AREA OVERLAY — a cell 37 % covered gets 37, not the '
               'class that happened to sit under its centre.')
    for dotted, num, singular, _c in schema.TABLES:
        if num == 3:
            panelui.table_form(cfg, dotted, singular, path)

    names = [v.name for v in cfg.surface.vegetation]
    if names:
        st.caption('Vegetation types available from panel 2: ' +
                   ', '.join('**%d** %s' % (k + 1, n)
                             for k, n in enumerate(names)))

    par = os.path.join(str(ds), cfg.soil.params.replace('/', os.sep))
    st.markdown('#### Soil column parameters')
    if os.path.exists(par):
        st.caption('`%s` — the zone ORDER in this file is what the zone codes '
                   'above refer to.' % par)
        with st.expander('Show the file'):
            st.code(open(par, encoding='utf-8', errors='replace').read())
    else:
        st.error('Missing: %s' % par)

# --------------------------------------------------------------- aquifer
with tab_aq:
    for section, title in [('layers', 'Layers'), ('uzf', 'Unsaturated zone'),
                           ('seep', 'Seepage face'), ('et', 'Evapotranspiration'),
                           ('spinup', 'Spin-up & initial state')]:
        st.markdown('#### %s' % title)
        edited.update(panelui.section_form(cfg, section, columns=3))
        st.markdown('')
    if not cfg.spinup.strt_heads:
        st.warning('No saved initial heads. A cold start puts the water table '
                   'above ground over much of the catchment, and the first '
                   'weeks measure how the grid relaxes that rather than the '
                   'hydrology — runoff and exfiltration reach tens of times '
                   'precipitation. Fine for a smoke test, misleading for '
                   'anything else.')

# ------------------------------------------------------- surface water
with tab_water:
    st.markdown('#### Stream routing (SFR)')
    st.caption('The network is the hydrography you MAPPED, burned onto the '
               'grid at run time: a cell is a stream cell when a line '
               'actually crosses it. Width and incision are resolved after '
               'routing, because contributing area is only known once the '
               'reaches are ordered.')
    edited.update(panelui.section_form(cfg, 'sfr', columns=2))

    stream = os.path.join(str(ds), cfg.sfr.source)
    if os.path.exists(stream):
        n = sum(1 for ln in open(stream, encoding='utf-8')
                if ln.strip() and not ln.startswith('#'))
        st.caption('`%s` — %d vertex row(s)' % (cfg.sfr.source, n - 1))
    else:
        st.error('Missing: %s — run the converter from panel 1.' % stream)

    st.markdown('#### Ponds (LAK)')
    st.caption('One EMBEDDEDV lake per pond: every La Mata pond is smaller '
               'than a cell, so there is nothing to excavate. The builder '
               'needs the POLYGONS, not the centroid table.')
    edited.update(panelui.section_form(cfg, 'lak', columns=2))

    st.markdown('#### Runoff cascade (CRR)')
    edited.update(panelui.section_form(cfg, 'crr', columns=2))

    st.info('Open-water evaporation is **MODFLOW\'s** now. MMsoil used to '
            'evaporate from its own surface store; that store is gone, so the '
            'coupler writes the rate into SFR and LAK every stress period '
            'from the `Eo` forcing, and reads back what was actually removed '
            'so it still appears in the water balance.')

# ---------------------------------------------------------- observations
with tab_obs:
    edited.update(panelui.section_form(cfg, 'obs', columns=2))
    tbl = os.path.join(str(ds), cfg.obs.table)
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
        with st.expander('Show the file'):
            st.code(open(tbl, encoding='utf-8', errors='replace').read())
        st.caption('`##` before a name means the point is not drawn on the '
                   'maps; a single `#` comments it out of the run entirely.')
    else:
        st.error('Missing: %s' % tbl)

panelui.save_button(cfg, path, edited)
