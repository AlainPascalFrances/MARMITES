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
from lib import loaders, panelui, schema    # noqa: E402

st.set_page_config(page_title='4 Subsurface', page_icon='🌍', layout='wide')
case = st.session_state.get('case', 'LaMata')
cfg, path = panelui.pick_config()
ds = panelui.dataset_banner(cfg)
panel = panelui.header(4)

edited = panelui.panel_switch(cfg, panel)

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
            '**Active cells ARE asked**, one map per layer. The catchment '
            'polygon on %s gives the outline, but a layer can pinch out '
            'inside it — here layer 1 is absent in 84 cells where layer 2 '
            'is present, with a thickness still written there — and a '
            'polygon cannot see that. What the polygon is instead is the '
            'GEOGRAPHIC REFERENCE: the run checks that the active cells sit '
            'inside it, and stops if they do not.'
            % (schema.panel_name(1), schema.panel_name(3),
               schema.panel_name(1)))

    # WHAT %d RESOLVES TO, spelled out. A pattern is worth nothing if the
    # modeller has to guess whether `k_%d.asc` means k_1.asc or k_01.asc --
    # so the panel opens the expansion and says whether the files are there.
    _props = ('ibound', 'thickness', 'k', 'k33', 'ss', 'sy')
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
        st.success('Every field on this tab reaches the run: `ibound` '
                   'becomes `idomain`, `thickness` becomes `botm`, `k` and '
                   '`k33` become `ModflowGwfnpf`, `ss` and `sy` become '
                   '`ModflowGwfsto`, and none of them is taken from the '
                   'MODFLOW parameter file any more.')

# -------------------------------------------------------------------- ghb
# One sub-panel per MF6 package. The switch is the ini's ghb_yn: with it
# off nothing below is read, and the fields are drawn read-only rather than
# live so the panel cannot be left saying something the run will not do.
with tab_ghb:
    st.caption('A head held outside the model, and the conductance of what '
               'lies between it and the boundary cells — `ModflowGwfghb`.')
    edited.update(panelui.rows_form(
        cfg, schema.GHB_ROWS, 'ghb', columns=2, folder=str(mm_paths.GIS),
        files=schema.BOUNDARY_FILES,
        gated={'ghb.enable': ('ghb.line', 'ghb.layers', 'ghb.head',
                              'ghb.cond', 'ghb.cond_per')}))
    panelui.boundary_note(cfg, 'ghb', 'head', ds)

# -------------------------------------------------------------------- drn
with tab_drn:
    st.caption('The outflow boundary of the catchment — `ModflowGwfdrn`. '
               'Water leaves a cell once its head rises above the drain.')
    edited.update(panelui.rows_form(
        cfg, schema.DRN_ROWS, 'drn', columns=2, folder=str(mm_paths.GIS),
        files=schema.BOUNDARY_FILES,
        gated={'drn.enable': ('drn.line', 'drn.layers', 'drn.elevation',
                              'drn.cond', 'drn.at_layer_base',
                              'drn.cond_per')}))
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
    _sd = float(cfg.uzf.surfdep)
    st.caption('To reproduce UZF\'s own seepage (SIMULATE_GWSEEP, as in '
               'Daoud et al. 2022) with this drain: base %g m (SURFDEP/2), '
               'smoothing depth %g m (SURFDEP), conductance from `uzf` -- '
               'with the UZF tab\'s SURFDEP of %g m. The two use the same '
               'cubic ramp in MODFLOW 6, so this is exact; the stream cells '
               'stay with SFR.' % (_sd / 2.0, _sd, _sd))
    if cfg.seep.kind != 'drn':
        st.warning('`seep.kind = %r` builds the seepage face inside UZF '
                   '(`SIMULATE_GWSEEP`) instead, so no `drn_seep` package '
                   'exists and the conductance above is not read. That '
                   'option is deprecated since MODFLOW 6.5, which '
                   'recommends the drain; a coupled run refuses it. The '
                   'drain reproduces it exactly (see the help of *Seepage '
                   'mechanism*).' % cfg.seep.kind)

# -------------------------------------------------------------------- uzf
with tab_uzf:
    st.caption('A vertical column of UZF objects per active cell — '
               '`ModflowGwfuzf`. The soil column above it is MARMITES\'s; '
               'this is what happens between the bottom of that column and '
               'the water table.')
    edited.update(panelui.rows_form(cfg, schema.UZF_ROWS, 'uzf', columns=2))
    # `gated` takes a SWITCH, and this one is a choice of source, so the
    # dependency is said in words rather than drawn as a greyed box.
    if cfg.uzf.thtr_from == 'sy':
        st.caption('`uzf.thtr` above is not read: the residual water content '
                   'is `thts - Sy` in every cell, so the unsaturated zone '
                   'drains exactly what the aquifer stores. An initial '
                   'content below it is raised to it at the build.')
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
    st.caption('**UZF always simulates unsaturated-zone ET — there is no '
               'switch for it.** Total ET has three sources and the deep '
               'unsaturated zone is one of them, so an off position could '
               'only ever produce a model that evaporates nothing from it. '
               'What is asked here is *how*: the formulation, the depth it '
               'stops at, and the water content it stops at.')
    edited.update(panelui.rows_form(cfg, schema.UZF_ET_ROWS, 'et',
                                    columns=2))
    if cfg.et.extdp_from == 'vegetation':
        st.caption('**Per vegetation zone.** Each cell takes the rooting '
                   'depth of its cover — %s — weighted by the cover; the '
                   'fraction nothing covers takes the extinction depth '
                   'above, and an irrigated field its crops\' root depth.'
                   % ', '.join('%s %g m' % (v.name, float(v.root_depth))
                               for v in cfg.surface.vegetation))
    st.info('**Wired end to end (WP2).** UZF simulates unsaturated-zone ET '
            'with this extinction depth, and groundwater ET is deliberately '
            'NOT asked of it — `linear_gwet` and `square_gwet` are both '
            'omitted, so ETg stays with MARMITES through the WEL package.\n\n'
            'The demand chain: what the soil column leaves of PE and PT is '
            'written to UZF each stress period (PETMAX); UZF\'s **actual** '
            'ET is read back from its budget and only the residual is left '
            'to ETg — never the demand, which a dry deep zone cannot meet. '
            'In lagged mode the actual is the previous stress period\'s.')

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
    # The saved-state names are PICKED, below, from what is on disk: they
    # were free text boxes, blank until a name only a log line gave was
    # typed in, with nothing saying where the run looks or what fits.
    edited.update(panelui.section_form(
        cfg, 'spinup', columns=3,
        exclude=('spinup.strt_heads', 'spinup.steady_means',
                 'spinup.save_strt', 'spinup.save_means')))
    state_dir = loaders.state_dir_for(cfg, mm_paths.WS_ROOT)
    grid_now = loaders.run_grid(cfg, state_dir)
    st.markdown('#### Where the run starts')
    st.caption('Saved states are looked for in `%s`, which is where a run '
               'saves them. %s'
               % (state_dir,
                  ('This run\'s grid: %s cells (the mesh cached there), so '
                   'each entry says whether it fits it.'
                   % mcfg.cells_text(grid_now['shape'])) if grid_now else
                  'Each entry says whether its recorded grid and layer set '
                  'match these settings.'))
    c1, c2 = st.columns(2)
    with c1:
        strt = panelui.state_picker(
            cfg, 'spinup.strt_heads', state_dir, grid_now,
            'none - start from the DEM, through a steady period')
    with c2:
        means = panelui.state_picker(
            cfg, 'spinup.steady_means', state_dir, grid_now,
            'none - one uniform recharge, no groundwater ET',
            disabled=bool(strt))
        if strt:
            st.caption('Not used: a run from saved heads has no steady '
                       'period to drive.')
    edited['spinup.strt_heads'] = strt
    edited['spinup.steady_means'] = means

    st.markdown('#### Where the run ends')
    # ONE name for the end state. The heads (+ _state.npz) and the steady
    # means were two fields always given the same name -- the spin-up button
    # set both -- and a single run that named only the heads saved no means.
    got = panelui.section_form(cfg, 'spinup', columns=2,
                               only='spinup.save_strt', skip_subpanels=False)
    save = (got.get('spinup.save_strt') or '').strip()
    edited.update(got)
    edited['spinup.save_means'] = save
    cycles = int(edited.get('spinup.cycles', cfg.spinup.cycles) or 1)
    _there = {d['name']: d for d in mcfg.saved_states(
        state_dir, int(cfg.layers.nlay), 'spinup.strt_heads')}
    if save and save in _there:
        _d = _there[save]
        st.caption('⚠ `%s` is already saved there (%s cells): this run '
                   'overwrites it when it ends.'
                   % (save, mcfg.cells_text(_d['shape'])
                      if _d['shape'] else '?'))
    elif not save:
        st.caption('Blank: a single run saves nothing; a spin-up (cycles > 1) '
                   'saves as `hi_spinup` anyway.')
    st.markdown('')
    if not strt and cycles > 1:
        st.info('**This run is a spin-up.** It starts cold, from the DEM '
                '(`spinup.strt_dem`, elevation - 2 m by default) through a '
                'steady period, repeats the forcing up to %d times until the '
                'water table moves less than %g m between cycles, and saves '
                'where it ends as `%s`. Then pick `%s` under *Initial heads* '
                'and set cycles to 1 for the runs that follow.'
                % (cycles, float(edited.get('spinup.tol', cfg.spinup.tol)),
                   save or 'hi_spinup', save or 'hi_spinup'))
    elif not strt:
        st.warning('No saved initial heads. A cold start puts the water table '
                   'above ground over much of the catchment, and the first '
                   'weeks measure how the grid relaxes that rather than the '
                   'hydrology — runoff and exfiltration reach tens of times '
                   'precipitation. Fine for a smoke test, misleading for '
                   'anything else.')
    else:
        # The SAME question the run asks at launch, asked here instead: saved
        # state belongs to the grid and layer set that produced it. Hearing
        # that after pressing Run is hearing it too late.
        #
        # IT IS NOT A REFUSAL. This said "This run will not start", which was
        # true when the guard raised a CONFIG ERROR and false from the moment
        # resolve_initial_heads started falling back to the land surface and
        # saying so. A panel that threatens a refusal the run does not make
        # is the same defect as a switch the run does not read.
        why = mcfg.prefix_problem(cfg, state_dir, 'spinup.strt_heads',
                                  strt, grid_now)
        if why:
            st.warning('**The saved heads will not be used.** %s\n\nThe run '
                       'does not refuse: it starts the water table from the '
                       'land surface (`spinup.strt_dem`) and says so in the '
                       'log. That is a cold start, with the first weeks '
                       'measuring how the grid relaxes rather than the '
                       'hydrology.' % why)
            # HOW TO REGENERATE IT, which the sentence above asks for and
            # nothing else on this panel explained. A spin-up is not a
            # button: it is a RUN -- minutes to hours -- so it goes through
            # the Run panel like any other. What this does is set the fields
            # that make the next run a spin-up, and say what to do next.
            st.markdown('**Regenerate it on this grid**')
            c1, c2 = st.columns([1, 3])
            _name = '%s_%s' % (cfg.grid_kind, cfg.meta.model_name(
                cfg.paths.case) or 'spinup')
            if c1.button('Set up a spin-up run', key='mkspinup'):
                # Parked rather than assigned: these widgets already exist
                # this run, and session_state cannot be written to once they
                # do -- see panelui.park.
                panelui.park('spinup.cycles', max(int(cfg.spinup.cycles), 6))
                panelui.park('spinup.strt_heads', '')
                panelui.park('spinup.steady_means', '')
                panelui.park('spinup.save_strt', _name)
                st.rerun()
            c2.caption('Sets `cycles` to 6, clears the state that does not '
                       'fit, and names the new one `%s`. Then **Validate & '
                       'save** in the sidebar and launch from the %s panel: '
                       'the run repeats the forcing until the water table '
                       'stops moving between cycles (`tol`), writes the '
                       'equilibrated heads, and prints the line to paste '
                       'back into `strt_heads`.'
                       % (_name, schema.panel_name(8)))
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

# REMEMBERED, not written: the one save is in the sidebar (see panelui).
# Collected HERE because a panel collects tab by tab as the tabs are drawn,
# so remembering earlier would remember an empty dict.
panelui.remember(edited)
panelui.sidebar_save(cfg, path)
