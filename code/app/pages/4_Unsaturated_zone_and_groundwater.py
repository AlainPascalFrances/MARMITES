# -*- coding: utf-8 -*-
"""Panel 4 -- UNSATURATED ZONE AND GROUNDWATER: what MODFLOW 6 reads.  WP1d.

The layers, the unsaturated zone, the seepage face,
evapotranspiration, the surface-water packages and the state a run
starts from.

The master switch is the SAME field as panel 3's. MMsoil is stepped
from inside the MODFLOW time loop through the API, so there is no
"run the soil balance alone" any more -- the legacy mode belonged to
the Picard loop Phase 1 removed. Turning it off on either panel turns
it off on both.
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP, os.path.join(CODE, 'ppMF6')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg             # noqa: E402
import mm_paths                            # noqa: E402
from lib import panelui, schema             # noqa: E402

st.set_page_config(page_title='4 Subsurface', page_icon='🌍', layout='wide')
case = st.session_state.get('case', 'LaMata')
cfg, path = panelui.pick_config()
ds = panelui.dataset_banner(cfg)
panel = panelui.header(4)

edited, save_slot = panelui.switch_and_save(cfg, panel)

tab_geom, tab_aq, tab_water = st.tabs(
    ['Geometry', 'Aquifer & solver', 'Streams, ponds & runoff'])

# ---------------------------------------------------------------- geometry
# The first of one sub-panel per MF6 package. Every field says which flopy
# argument it becomes, because that is the only name that does not drift.
with tab_geom:
    st.caption('The layer stack, and what each layer is made of. Every field '
               'names the flopy argument it becomes.')
    edited.update(panelui.rows_form(cfg, schema.GEOMETRY_ROWS, 'layers',
                                    columns=2))

    st.info('**The top is not asked.** The aquifer top is the land surface '
            'minus the soil column — elevation from the DEM on %s, thickness '
            'from %s — which is the bottom of the MARMITES soil column and '
            'the surface groundwater discharges at. A third answer could only '
            'disagree with the other two.\n\n'
            '**Nor is ibound.** A cell is active when it is inside the '
            'catchment polygon of %s: the same polygon the grid was built '
            'inside, so a separate map could only contradict it.'
            % (schema.panel_name(1), schema.panel_name(3),
               schema.panel_name(1)))

    # WHAT THE RUN ACTUALLY READS, said plainly. hnoflo and nlay reach the
    # model; the four properties are still taken from the MODFLOW parameter
    # file until the converter resolves them, and a panel that implied
    # otherwise would be the decorative switch all over again.
    _wired = [d for d in ('layers.nlay', 'layers.hnoflo')]
    _pending = [d for d in ('layers.thickness', 'layers.k', 'layers.ss',
                            'layers.sy')
                if getattr(cfg.layers, d.split('.')[-1]).producer() is None]
    if _pending:
        st.warning('Not read by a run yet: %s. Those still come from the '
                   'MODFLOW parameter file `%s`; naming a raster or a layer '
                   'here records the intent but does not change the run until '
                   'the converter resolves them.'
                   % (', '.join('`%s`' % d for d in _pending),
                      'MF_ws/__inputMF_flopy_v3_*.ini'))

# --------------------------------------------------------------- aquifer
with tab_aq:
    for section, title in [('uzf', 'Unsaturated zone'),
                           ('seep', 'Seepage face'),
                           ('et', 'Evapotranspiration'),
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
    else:
        # The SAME question the run asks at launch, asked here instead: saved
        # state belongs to the grid and layer set that produced it, and the
        # run refuses to reuse it otherwise. Hearing that after pressing Run
        # is hearing it too late.
        why = mcfg.state_problem(cfg, mcfg.state_workspace(cfg,
                                                           mm_paths.WS_ROOT))
        if why:
            st.error('**This run will not start.** %s' % why)
        else:
            st.success('The saved state belongs to this grid and layer set.')

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
        st.error('Missing: %s — build the grid on the %s panel, which\n'
                 'reads the cartography first.'
                 % (stream, schema.panel_name(1)))

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

panelui.save_button(cfg, path, edited, slot=save_slot)
