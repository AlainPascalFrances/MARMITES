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
    ['MODFLOW aquifer layers', 'Aquifer & solver',
     'Streams, ponds & runoff'])

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

    # WHAT %d RESOLVES TO, spelled out. A pattern is worth nothing if the
    # modeller has to guess whether `k_%d.asc` means k_1.asc or k_01.asc --
    # so the panel opens the expansion and says whether the files are there.
    _props = ('thickness', 'k', 'k33', 'ss', 'sy')
    _pat = [(n, getattr(cfg.layers, n)) for n in _props
            if '%d' in (getattr(cfg.layers, n).raster or '')]
    if _pat:
        lines, wrong_case = [], []
        for name, src in _pat:
            marks = []
            for fn in src.rasters(cfg.layers.nlay):
                # CASE-STRICT on purpose. Windows opens Ss_l2.asc when the
                # file is called ss_l2.asc and Linux does not, so a check
                # that trusted the filesystem would hide the one mistake a
                # pattern makes easy.
                folder = os.path.dirname(os.path.join(str(ds or ''), fn))
                base = os.path.basename(fn)
                try:
                    here = os.listdir(folder)
                except OSError:
                    here = []
                if base in here:
                    marks.append('🟢 `%s`' % fn)
                elif base.lower() in [h.lower() for h in here]:
                    real = [h for h in here if h.lower() == base.lower()][0]
                    marks.append('🟠 `%s` (on disk: `%s`)' % (fn, real))
                    wrong_case.append((fn, real))
                else:
                    marks.append('🔴 `%s`' % fn)
            lines.append('- **%s** — %s' % (name, ', '.join(marks)))
        st.markdown('`%%d` over %d layer(s):\n\n%s'
                    % (cfg.layers.nlay, '\n'.join(lines)))
        if wrong_case:
            st.warning('Same name, different case: %s. Windows opens these '
                       'and Linux does not — rename the file or the '
                       'pattern so one spelling covers every layer.'
                       % '; '.join('`%s` vs `%s`' % (a, b)
                                   for a, b in wrong_case))

    # WHAT THE RUN ACTUALLY READS, said plainly -- and it now reads all of
    # it. Each answer replaces what the MODFLOW parameter file parsed,
    # before the cell list, the soil model or any MF6 package has seen it;
    # what is left blank still comes from that file, which is the only
    # reason it is still opened at all.
    _blank = [n for n in _props if getattr(cfg.layers, n).producer() is None]
    if _blank:
        st.warning('Still read from `MF_ws/__inputMF_flopy_v3_*.ini`: %s. '
                   'Answer it here and the run uses the answer instead.'
                   % ', '.join('`layers.%s`' % n for n in _blank))
    else:
        st.success('Every field on this tab reaches the run: `thickness` '
                   'becomes `botm`, `k` and `k33` become `ModflowGwfnpf`, '
                   '`ss` and `sy` become `ModflowGwfsto`, and none of them '
                   'is taken from the MODFLOW parameter file any more.')

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
