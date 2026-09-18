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

tab_geom, tab_ghb, tab_drn, tab_uzf, tab_water, tab_init = st.tabs(
    ['MODFLOW aquifer layers', 'GHB', 'DRN', 'UZF', 'SFR, LAK and CRR',
     'Initial heads & spin-up'])

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

# -------------------------------------------------------------------- ghb
# One sub-panel per MF6 package. The switch is the ini's ghb_yn: with it
# off nothing below is read, and the fields are drawn read-only rather than
# live so the panel cannot be left saying something the run will not do.
with tab_ghb:
    st.caption('A head held outside the model, and the conductance of what '
               'lies between it and the boundary cells — `ModflowGwfghb`.')
    edited.update(panelui.rows_form(
        cfg, schema.GHB_ROWS, 'ghb', columns=2,
        gated={'ghb.enable': ('ghb.layers', 'ghb.head', 'ghb.cond')}))
    panelui.boundary_note(cfg, 'ghb', 'head', ds)

# -------------------------------------------------------------------- drn
with tab_drn:
    st.caption('The outflow boundary of the catchment — `ModflowGwfdrn`. '
               'Water leaves a cell once its head rises above the drain.')
    edited.update(panelui.rows_form(
        cfg, schema.DRN_ROWS, 'drn', columns=2,
        gated={'drn.enable': ('drn.layers', 'drn.elevation', 'drn.cond',
                              'drn.at_layer_base')}))
    panelui.boundary_note(cfg, 'drn', 'elevation', ds)

    # THE SECOND DRAIN PACKAGE. MARMITES builds two, and they were on
    # different tabs, which is how someone hunting the seepage face ends up
    # switching off the catchment outlet. Both are `ModflowGwfdrn`, so both
    # belong on the DRN tab -- with the difference stated rather than
    # cross-referenced.
    st.markdown('---')
    st.markdown('#### Seepage face — the second drain package')
    st.caption('Groundwater leaving at the LAND SURFACE, over the whole '
               'catchment, as `drn_seep`. The boundary drain above is six '
               'cells at the outlet; this one is every surface cell, and '
               'MARMITES takes its discharge back into the soil column as '
               'exfiltration rather than routing it away.')
    edited.update(panelui.section_form(cfg, 'seep', columns=3))
    if cfg.seep.kind != 'drn':
        st.warning('`seep.kind = %r` builds the seepage face inside UZF '
                   '(`SIMULATE_GWSEEP`) instead, so no `drn_seep` package '
                   'exists and the conductance above is not read. That '
                   'option is deprecated in MODFLOW 6 and switches discharge '
                   'on and off discontinuously; `drn` is the validated '
                   'choice and what La Mata runs.' % cfg.seep.kind)

# -------------------------------------------------------------------- uzf
with tab_uzf:
    st.caption('A vertical column of UZF objects per active cell — '
               '`ModflowGwfuzf`. The soil column above it is MARMITES\'s; '
               'this is what happens between the bottom of that column and '
               'the water table.')
    edited.update(panelui.rows_form(cfg, schema.UZF_ROWS, 'uzf', columns=2))
    # `gated` takes a SWITCH, and this one is a choice of source, so the
    # dependency is said in words rather than drawn as a greyed box.
    if cfg.uzf.vks_from != 'raster':
        st.caption('`uzf.vks` above is not read: the unsaturated vertical K '
                   'comes from each layer\'s own `k33`, so the column and '
                   'the aquifer cannot disagree.')

    st.markdown('#### Evapotranspiration inside MODFLOW')
    st.info('**Groundwater ET is not asked, here or anywhere.** MARMITES '
            'computes ETg and applies it through the WEL package, so a '
            'second answer inside MODFLOW could only remove the same water '
            'twice. The old switch for it was forced off and read by '
            'nothing, which is worse than absent.')
    edited.update(panelui.rows_form(
        cfg, schema.UZF_ET_ROWS, 'et', columns=2,
        gated={'et.uzf_et': ('et.extdp', 'et.extwc_source',
                             'et.unsat_form')}))
    if cfg.et.uzf_et:
        st.warning('`simulate_et` is currently hard-coded to False when the '
                   'model is built, so this switch and the extinction depth '
                   'are not read yet. Wiring it means UZF drying the '
                   'unsaturated zone alongside MMsoil, which is a modelling '
                   'decision — ask before relying on it.')

    with st.expander('Where every UZF1 name went'):
        st.caption('The parameter file carried twenty-two of these and '
                   'MODFLOW 6 keeps eight. "Gone" is only a useful answer '
                   'with the reason attached.')
        st.markdown('\n'.join('- `%s` — %s' % (a, b)
                               for a, b in schema.UZF_LEGACY))

# ------------------------------------------------------------ initial heads
# Last, because it is the only tab that is not about WHAT the model is: it
# is where the model starts from, which is the last thing answered and the
# first thing a run reads.
with tab_init:
    st.caption('The heads a run starts from, and the spin-up that produces '
               'them — `ModflowGwfic(strt=)`. Nothing here changes the '
               'model; it changes where the model begins.')
    edited.update(panelui.section_form(cfg, 'spinup', columns=3))
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
