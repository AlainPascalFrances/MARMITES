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

tab_soil, tab_veg, tab_aq, tab_water, tab_obs, tab_gis = st.tabs(
    ['Soil column', 'Vegetation characteristics', 'Aquifer & solver',
     'Streams, ponds & runoff', 'Observations', 'Cartography'])

# ------------------------------------------------------------------ soil
with tab_soil:
    edited.update(panelui.rows_form(cfg, schema.SOIL_ROWS, 'soil',
                                    columns=2, folder=str(ds),
                                    files=schema.SOIL_DATASET_FILES))

    st.info('**Precedence is explicit, not decided by which file exists:** a '
            'raster beats a polygon attribute, and a polygon attribute beats '
            'a single value. Nothing set is an error, not a silent zero.')

    # The BOX, not the saved file: a green light against the file named
    # before the last edit is worse than no light at all.
    _par = panelui.live('soil.params', cfg.soil.params) or ''
    par = (_par if os.path.isabs(_par)
           else os.path.join(str(ds), _par.replace('/', os.sep)))
    st.markdown('#### Soil column parameters')
    if os.path.exists(par):
        st.caption('`%s` — the zone ORDER in this file is what the zone codes '
                   'above refer to.' % par)
        with st.expander('Show the file'):
            st.code(open(par, encoding='utf-8', errors='replace').read())
    else:
        st.error('Missing: %s' % par)

# ------------------------------------------------------------ vegetation
# [soil] carries the vegetation COVER as well, because they share a file --
# not because they are one question. Shown beside the soil settings the
# second was read as more of the first.
with tab_veg:
    edited.update(panelui.rows_form(cfg, schema.SOIL_VEG_ROWS, 'soil',
                                    columns=2, folder=str(mm_paths.GIS),
                                    files=schema.SOIL_FILES))

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

# ------------------------------------------------------- the cartography
# HERE and not on panel 1, because these are the layers PANEL 3 wraps: the
# soil zones, the vegetation cover, the observation points, the pond
# outlines. Panel 1 needs no button of its own -- the two tables a GRID
# depends on, the catchment ring and the stream network, are checked and
# re-read by Create grid itself.
#
# In a TAB of its own, not at the foot of the page: content outside the tabs
# is drawn under whichever one is open, so this read as part of the soil
# column -- which is exactly what it is not.
with tab_gis:
    st.caption(
        'A run never opens a shapefile. The converter does, once: it reads '
        'the GIS folder and writes GRID-INDEPENDENT tables into `%s`, and '
        'those are what a run reads and what this panel wraps onto the grid. '
        'So press this when you have EDITED OR REPLACED a shapefile — '
        'changing the grid does not need it, which is the whole point of the '
        'two tiers.' % ds)

    def _converter(dry):
        """The converter as a SUBPROCESS: it keeps geopandas out of this
        process's import graph, and it is the command a user would type."""
        import subprocess
        cmd = [sys.executable,
               os.path.join(CODE, 'tools', 'gis_to_dataset.py'),
               '--case', case, '--config', path]
        if dry:
            cmd.append('--dry-run')
        r = subprocess.run(cmd, capture_output=True, text=True,
                           cwd=str(mm_paths.REPO), timeout=900)
        return (r.stdout or '') + (('\n' + r.stderr) if r.stderr else '')

    c1, c2 = st.columns(2)
    if c1.button('Preview (dry run)', key='conv_dry'):
        st.session_state['conv'] = _converter(True)
    if c2.button('Update dataset', type='primary', key='conv_run'):
        st.session_state['conv'] = _converter(False)
    if st.session_state.get('conv'):
        st.code(st.session_state['conv'], language='text')
