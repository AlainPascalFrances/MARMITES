# -*- coding: utf-8 -*-
"""Panel 7 -- RUN: validate the configuration, then launch it and follow it.

Two tabs, in the order they are used:

  Validation of the configuration   COMPULSORY. Everything that can be said
      about the configuration before a run -- the run's own settings, the
      one-off overrides, the checks (``lib.checks``), whether the dataset
      matches the cartography -- and ONE button, *Validate the configuration*: it saves what
      the panels hold, checks the result and, when nothing stops it,
      approves it. Errors keep it from approving; warnings are read here and
      approved with it.
  Run   FROZEN until then, and again the moment anything changes -- an edit
      on any panel, a different override, another configuration. Its
      *Launch* is the only thing in the app that starts a model, and it
      never sends the modeller anywhere else.

It used to be two panels, Validation (7) and Run (8): Launch on the Run
panel diverted to Validation on any warning, and Launch on Validation
started the SAVED file -- so an edit not yet saved could be left out of the
run. A one-year run went without the stream that way on 2026-09-26.

The run is DETACHED. Streamlit reruns this script on every interaction, so a
model call in the page body would block the app for the duration of a
coupled run and die on a browser refresh. Instead the run writes a log and a
status file, and this page polls them -- so a run survives closing the tab.

Server mode is not a different mechanism: run Streamlit ON the server
(`streamlit run code/app/Home.py --server.address 0.0.0.0 --server.port 8501`)
and it launches there, next to the workspace.
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg        # noqa: E402
import mm_paths                       # noqa: E402
from lib import checks as chk          # noqa: E402
from lib import dataset_state          # noqa: E402
from lib import panelui               # noqa: E402
from lib import runs as runlib        # noqa: E402
from lib import schema                # noqa: E402

st.set_page_config(page_title='7 Run', page_icon='▶️', layout='wide')

cfg, path = panelui.pick_config()
chosen = os.path.basename(path)
panelui.dataset_banner(cfg)

st.title('▶️  7 — Run')
panelui.saved_note()
st.caption('Validate the configuration on the first tab; the run tab stays '
           'frozen until you have, and freezes again as soon as anything '
           'changes. The model runs only from **Launch** on the run tab.')

RUNS = cfg.ui.runs_dir or os.path.join(str(mm_paths.WS_ROOT), 'runs')


def _panel_title(num):
    return schema.panel_name(num) if num is not None else ''


def _effective(overrides):
    """(the configuration that would RUN -- the saved file with the
    overrides applied -- or None, and the error that stopped it)"""
    if not overrides:
        return cfg, None
    try:
        probe = mcfg.load_run_config(path)
        probe.apply_overrides(overrides, echo=False)
        return probe, None
    except mcfg.ConfigError as exc:
        return None, str(exc)


tab_val, tab_run = st.tabs(['🔎 Validation of the configuration', '▶️ Run'])

# ================================================================ validation
with tab_val:
    _vm = st.session_state.pop('__validate_msg', None)
    if _vm:
        (st.success if _vm[0] else st.error)(_vm[1])
    # WHICH GRID THIS RUN WILL USE, stated rather than inferred. The grid is
    # not committed by the sidebar save: the Grid panel's *Select this grid
    # for the model* writes the [grid] settings AND promotes the mesh they
    # produced to where the driver looks. So it is the one setting whose
    # current value cannot be read off the panel you are standing on.
    st.markdown('#### The grid this run will use')
    _g = chk.describe_grid(cfg, str(mm_paths.WS_ROOT))
    _gc1, _gc2, _gc3 = st.columns(3)
    _gc1.metric('Producer', _g['kind'])
    # NOT 'cell size' on a voronoi mesh: that number is grid.voronoi.cell_far,
    # the far-field target, and the cells themselves vary a long way from it.
    _gc2.metric(_g['size_what'].capitalize(),
                '%g m' % _g['size'] if _g['size'] else '—')
    _gc3.metric('Cells', '{:,}'.format(_g['ncpl']) if _g['ncpl']
                else ('rectangle' if _g['kind'] in ('structured', 'dis')
                      else 'not built'))
    if _g['kind'] in ('structured', 'dis'):
        st.caption('A structured run builds its rectangle from the catchment '
                   'boundary; there is no mesh to cache.')
    elif _g['cached']:
        st.caption('Read from `%s`. It got there from the Grid panel — '
                   '*Select this grid for the model* is what commits a grid, '
                   'not the sidebar save, because it also promotes the mesh '
                   'it built. A run rebuilds it if `[grid]` has changed '
                   'since.' % _g['where'])
    else:
        st.warning('No mesh is cached for this producer, so **the run will '
                   'build it first** — minutes, on La Mata. Build it on the '
                   'Grid panel and press *Select this grid for the model*.')

    # THE RUN TAG, saved as meta.name. It sat on the run tab as a widget that
    # was never written (2026-10-06): a tag typed there named one run and was
    # gone at the next, back to whatever the file said. It is part of the
    # configuration (and of its hash), so it is asked HERE and saved by
    # Validate like every other answer on this tab.
    st.markdown('#### Run tag')
    _tc1, _tc2 = st.columns([2, 3])
    with _tc1:
        _tag_edit = panelui.rows_form(cfg, (('meta.name', None),), 'meta',
                                      columns=1)
    with _tc2:
        st.markdown('')
        _tag = str(_tag_edit.get('meta.name', cfg.meta.name) or '').strip()
        st.caption('Optional. The run writes `out_<timestamp>_%s`. Saved with '
                   'the configuration, so it stays until you change it.'
                   % (_tag or cfg.run_tag))
    panelui.remember(_tag_edit)

    # HOW THE RUN IS COUPLED, and what it must satisfy to count as a result.
    # They are part of what is validated, so they are asked HERE, before the
    # button. The lagged/iterative choice is gone (2026-10-07): there is one
    # coupling, said in words rather than asked.
    st.markdown('#### How the run is coupled')
    st.caption('Each stress period MMsoil runs first, from the heads the '
               'previous one ended with; then MODFLOW 6 solves the period '
               'with its own time step, retried shorter (ATS) when it fails. '
               'Groundwater ET is taken by EVT at the head MODFLOW solves '
               'for, so it follows the water table within the day.')
    _coupling = panelui.rows_form(cfg, schema.RUN_COUPLING_ROWS, 'run',
                                  columns=2)
    panelui.remember(_coupling)

    # The IMS solver. Its tolerances came from the legacy NWT ini and no
    # panel asked them; HEADTOL 0.05 m left a 4 % mass-balance discrepancy
    # once the UZF/Sy runaway was fixed and the flows became what they are.
    st.markdown('#### MODFLOW 6 solver')
    _solver = panelui.rows_form(cfg, schema.SOLVER_ROWS, 'solver', columns=2)
    panelui.remember(_solver)
    st.caption('Tighter tolerances close the water budget and cost '
               'iterations. A run whose cumulative discrepancy exceeds the '
               'limit above stops rather than being reported as a result.')

    # ASKED HERE, not edited into the file by hand: the one field standing
    # between a build and a coupled run. Blank builds and stops, set runs.
    st.markdown('#### MODFLOW 6 library')
    _lc1, _lc2 = st.columns([2, 3])
    with _lc1:
        _lib_edit = panelui.rows_form(cfg, (('paths.libmf6', None),),
                                      'paths', columns=1)
    libmf6 = str(_lib_edit.get('paths.libmf6', cfg.paths.libmf6) or '').strip()
    with _lc2:
        st.markdown('')
        if not libmf6:
            st.warning('Blank: the run will stop after writing the MF6 '
                       'files. Type `auto` for a coupled run.')
        else:
            # SETTLED HERE, where it can still be corrected, rather than
            # after the build -- the same resolver the driver uses.
            _rl = runlib.import_runner()
            try:
                _lib = _rl.check_libmf6(_rl.resolve_libmf6(libmf6))
            except _rl.LibMF6Error as exc:
                st.error(str(exc))
            else:
                st.success('`%s`' % _lib)
        st.caption('`auto` finds it next to the other MODFLOW binaries and '
                   'keeps the file portable. It is the shared LIBRARY — '
                   '`libmf6.dll` — not `mf6.exe`: the coupler steps MODFLOW '
                   'one stress period at a time through the API, which the '
                   'executable cannot do.')
    panelui.remember(_lib_edit)

    st.markdown('#### One-off overrides')
    st.caption('`section.key=value`, one per line, applied on top of the file '
               'and echoed by the driver. They are part of what is '
               'validated: change them and the run tab freezes again. There '
               'is deliberately no free-text command box: this page starts a '
               'process, so it only ever runs the validated entry point.')
    raw = st.text_area('Overrides', height=90, key='run_overrides',
                       placeholder='run.nsp=365\nsfr.enable=true')
    overrides = [x.strip() for x in raw.splitlines() if x.strip()]
    eff, bad = _effective(overrides)
    if bad:
        st.error(bad)
    elif overrides:
        st.success('Overrides valid — resulting hash `%s`'
                   % eff.config_hash())

    st.markdown('---')
    st.markdown('#### The checks')
    # UNSAVED EDITS FIRST, and on their own. They are not a property of the
    # configuration but of this browser session -- and they are exactly
    # what a run would silently leave out. Validating saves them.
    todo = panelui.unsaved_changes(cfg)
    if todo:
        st.warning('**%d unsaved change(s).** What follows is the FILE, '
                   'which is what a run reads. *Validate the configuration* saves '
                   'them first.'
                   % len(todo))
        with st.expander('What is unsaved'):
            for line in todo:
                st.code(line, language='ini')
    found = chk.collect_default(eff if eff is not None else cfg,
                                panelui.unsaved_switches(cfg))
    counts = chk.count_by_level(found)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric('Errors', counts[chk.ERROR])
    c2.metric('Warnings', counts[chk.WARNING])
    c3.metric('Notes', counts[chk.INFO])
    c4.metric('Hash', (eff if eff is not None else cfg).config_hash())
    if counts[chk.ERROR]:
        st.error('**This configuration cannot run.** Fix the errors below; '
                 'it cannot be validated until they are gone.')
    elif counts[chk.WARNING]:
        st.warning('**It will run.** Nothing here stops it, but each warning '
                   'is something that has misled a run before — read them, '
                   'then validate.')
    else:
        st.success('**Ready.** Nothing to report.')

    ICON = {chk.ERROR: '⛔', chk.WARNING: '⚠️', chk.INFO: 'ℹ️'}
    LABEL = {chk.ERROR: 'Errors', chk.WARNING: 'Warnings', chk.INFO: 'Notes'}
    for level in (chk.ERROR, chk.WARNING, chk.INFO):
        here = [c for c in found if c.level == level]
        if not here:
            continue
        st.markdown('##### %s %s' % (ICON[level], LABEL[level]))
        for c in here:
            where = ('panel %d — %s' % (c.panel, _panel_title(c.panel))
                     if c.panel is not None and _panel_title(c.panel) else '')
            # The key is shown ONLY when the sentence does not already carry
            # it: most checks name the key they are about as their first
            # words, and adding it underneath printed the same thing twice.
            head = '**%s**' % c.title
            if c.key and c.key not in c.title:
                head += '  \n`%s`' % c.key
            with st.container(border=True):
                st.markdown(head)
                if c.detail:
                    st.caption(c.detail)
                if where:
                    st.caption('Answered on %s.' % where)
    st.info('**The grid is committed on the Grid panel, not by the sidebar '
            'save.** *Select this grid for the model* writes the `[grid]` '
            'settings AND puts the mesh they produced where a run will find '
            'it — so a grid built and not selected is not the model\'s grid. '
            'Which one a run would use is at the top of this tab.')

    # THE CARTOGRAPHY -> DATASET CONVERTER, next to the list that reports
    # it. A run reads the converted tables, never a shapefile; lm_veg.shp was
    # edited on 2026-09-23 and every run that day read the vegetation of the
    # 13th. Launch converts by itself; this is the same conversion on demand.
    st.markdown('#### Cartography → dataset')
    _ds = mm_paths.dataset_dir(cfg.paths.case)
    st.caption('A run never opens a shapefile. The converter reads the GIS '
               'folder (`%s`) and writes grid-independent tables into the '
               'dataset (`%s`); those are what a run reads. **Launch '
               'converts first whenever a table is out of date**, and '
               '*Create grid* does it for the catchment ring and the '
               'streams.' % (mm_paths.GIS, _ds))
    try:
        _why = dataset_state.stale(cfg, _ds, mm_paths.GIS)
    except Exception as exc:                            # noqa: BLE001
        _why = None
        st.warning('The dataset state could not be read: %s' % exc)
    if _why:
        st.warning('**Out of date with the cartography:**\n\n- ' +
                   '\n- '.join(_why))
    elif _why is not None:
        st.success('Every converted table matches its shapefile.')
    _k1, _k2 = st.columns(2)
    if _k1.button('Preview (dry run)', key='conv_dry'):
        st.session_state['conv'] = dataset_state.run_converter(
            cfg.paths.case, path, dry=True)[1]
    if _k2.button('Update dataset', key='conv_run'):
        _ok, _out = dataset_state.run_converter(cfg.paths.case, path)
        st.session_state['conv'] = _out
        if _ok:
            st.rerun()                  # the list above is now out of date
    if st.session_state.get('conv'):
        with st.expander('Converter output', expanded=True):
            st.code(st.session_state['conv'], language='text',
                    wrap_lines=True)

    # THE ONE BUTTON THAT APPROVES. It saves first, because what is
    # validated must be what runs: the file, as the panels have it now.
    st.markdown('---')
    if st.button('Validate the configuration', type='primary',
                 key='validate',
                 disabled=bool(bad)):
        ok, msg = True, ''
        if todo:
            _applied, _whys, _saved = panelui.save_now(cfg, path)
            if not _saved:
                ok, msg = False, ('NOT saved, so not validated — the '
                                  'configuration would be invalid:\n\n%s'
                                  % _whys)
            else:
                msg = 'Saved %d change(s). ' % len(_applied)
        if ok:
            cfg = mcfg.load_run_config(path)
            eff, bad = _effective(overrides)
            _n = chk.count_by_level(chk.collect_default(
                eff, panelui.unsaved_switches(cfg)))
            if _n[chk.ERROR]:
                ok, msg = False, msg + ('%d error(s): not validated.'
                                        % _n[chk.ERROR])
                st.session_state.pop(panelui.VALIDATED, None)
            else:
                st.session_state[panelui.VALIDATED] = \
                    panelui.validation_record(chosen, eff.config_hash())
                msg += ('Validated (%d warning(s) read). Open the **Run** '
                        'tab to launch.' % _n[chk.WARNING])
        if not ok:
            st.session_state.pop(panelui.VALIDATED, None)
        st.session_state['__validate_msg'] = (ok, msg)
        st.rerun()
    st.caption('Saves what the panels hold, checks it again and approves it. '
               'The run tab opens only for exactly what was validated.')

# The one save, as on every other panel: this page edits run, solver and
# paths.libmf6. After the tab that remembers them, before the run tab's
# refresh loop.
panelui.sidebar_save(cfg, path)

# ======================================================================= run
with tab_run:
    todo_now = panelui.unsaved_changes(cfg)
    record = st.session_state.get(panelui.VALIDATED)
    eff_hash = eff.config_hash() if eff is not None else None
    unlocked = panelui.run_unlocked(record, chosen, eff_hash, todo_now)
    if not unlocked:
        st.info('🔒 **Frozen — %s** The model runs only from here, and only '
                'what the *Validation of the configuration* tab approved.'
                % panelui.frozen_reason(record, chosen, eff_hash, todo_now))
    else:
        st.success('✅ **Validated** — `%s`, hash `%s`%s.'
                   % (chosen, eff_hash,
                      (' with %d override(s)' % len(overrides))
                      if overrides else ''))
    # Read from what was validated, not typed here: an edit on this tab would
    # name the run without being saved (see the Run tag on the first tab).
    run_tag = (eff.meta.name if eff is not None else cfg.meta.name) or ''
    st.markdown('**Run tag** `%s` — set on the *Validation of the '
                'configuration* tab.' % (run_tag or '(none)'))

    def _start():
        """Actually start the detached run -- converting the cartography
        first when a table is out of date with its shapefile."""
        _busy_now = runlib.active(RUNS)
        if _busy_now:
            # checked BEFORE converting: the dataset is not rewritten under
            # a run that is still reading it
            st.error('Not launched: %s' % runlib.RunBusy(_busy_now[0]))
            return
        _stale = dataset_state.stale(cfg, mm_paths.dataset_dir(cfg.paths.case),
                                     mm_paths.GIS)
        if _stale:
            with st.spinner('Updating the dataset from the cartography…'):
                _ok, _out = dataset_state.run_converter(cfg.paths.case, path)
            if not _ok:
                st.error('The dataset could not be updated from the '
                         'cartography, so nothing was launched: a run would '
                         'read the tables as they were. The converter said:')
                st.code(_out[-4000:], language='text', wrap_lines=True)
                return
            st.info('Dataset updated from the cartography before the run: %s'
                    % '; '.join(_stale))
            with st.expander('Converter output'):
                st.code(_out, language='text', wrap_lines=True)
        try:
            run_id, payload = runlib.launch(
                path, RUNS, overrides=overrides, run_tag=run_tag or None,
                python_exe=(mm_paths.PYTHON_EXE
                            if os.path.exists(mm_paths.PYTHON_EXE)
                            else sys.executable))
        except runlib.RunBusy as exc:
            st.error('Not launched: %s' % exc)
            st.session_state['watch'] = exc.run.get('run_id')
        except Exception as exc:                        # noqa: BLE001
            st.error('Launch failed: %s' % exc)
        else:
            st.session_state['watch'] = run_id
            st.success('Launched `%s` (pid %s)' % (run_id, payload['pid']))

    # One run at a time: a second writes the same MODFLOW workspace
    # (runs.launch refuses it too -- this only says so before the press).
    _busy = runlib.active(RUNS)
    if _busy:
        st.warning('Run `%s` is still going (started %s). Launch is off until '
                   'it finishes or is stopped below.'
                   % (_busy[0].get('run_id'), _busy[0].get('started')))
    if st.button('Launch', type='primary', key='launch',
                 disabled=not unlocked or bool(_busy)):
        # asked again AT THE PRESS: a page drawn a minute ago may be stale
        if panelui.run_unlocked(st.session_state.get(panelui.VALIDATED),
                                chosen, eff_hash,
                                panelui.unsaved_changes(cfg)):
            _start()
        else:
            st.error('Not launched: the configuration changed after it was '
                     'validated. Validate it again.')

    st.markdown('---')
    st.markdown('#### Follow a run')
    known = runlib.list_runs(RUNS)
    if not known:
        st.info('No runs yet. `%s`' % RUNS)
    else:
        ids = [r['run_id'] for r in known]
        watch = st.session_state.get('watch')
        sel = st.selectbox('Follow', ids,
                           index=ids.index(watch) if watch in ids else 0)
        _first = runlib.status(RUNS, sel)
        auto = st.checkbox('Auto-refresh every %d s' % cfg.ui.poll_secs,
                           value=(_first.get('state') == 'running'))

        # REFRESHED AS A FRAGMENT, not by sleeping and rerunning the whole
        # page: that froze every widget above for the poll interval, and it
        # never let the script finish while a run was going (an AppTest of
        # this page timed out whenever a real run was live).
        @st.fragment(run_every=(max(1, int(cfg.ui.poll_secs))
                                if auto else None))
        def _follow(sel=sel):
            info = runlib.status(RUNS, sel)
            f1, f2, f3, f4 = st.columns(4)
            f1.metric('State', info.get('state', '?'))
            f2.metric('Outcome', info.get('outcome', '—'))
            # WHEN IT ENDED, under the outcome: the last write to the run's
            # log, and how long the run took
            if info.get('state') == 'finished':
                _end = runlib.finished_at(RUNS, sel)
                _took = runlib.duration(info.get('started'), _end)
                f2.caption('at %s%s' % (_end, (' — took %s' % _took)
                                        if _took else ''))
            f3.metric('Started', (info.get('started') or '')[-8:])
            f4.metric('PID', info.get('pid', '—'))
            st.code(' '.join(str(x) for x in info.get('cmd', [])),
                    language='bash', wrap_lines=True)

            # AN EMPTY LOG MEANS THE RUN WAS KILLED, not that it went well.
            # MF6 aborts the process from inside the library when it refuses
            # an input file, which takes python with it before anything is
            # written.
            if info.get('outcome') == 'unknown':
                _ws = mcfg.state_workspace(cfg, str(mm_paths.WS_ROOT))
                st.error('**This run left no log at all**, which means it was '
                         'killed rather than that it finished. MODFLOW 6 ends '
                         'the process itself when it refuses an input file, '
                         'which takes Python with it. The reason is in the '
                         'model workspace, not here:\n\n'
                         '- `%s` — the last lines name the file and the '
                         'records it rejected\n'
                         '- `%s` — the model listing, if it got that far'
                         % (os.path.join(str(_ws), 'mfsim.lst'),
                            os.path.join(str(_ws), '%s.lst'
                                         % cfg.meta.model_name(
                                             cfg.paths.case))))

            # st.code, NOT st.text_area: a text_area given both `value` and
            # `key` takes its content from session state after the first
            # render, so the log froze at whatever it was when first drawn.
            # Wrapped: the progress and warning lines run past the box, and
            # scrolling sideways to read each one is what nobody does.
            st.markdown('**run.log** (tail)')
            st.code(runlib.log_tail(RUNS, sel, 400) or '(empty)',
                    language='text', height=460, wrap_lines=True)
            if info.get('state') == 'running' and st.button('Stop this run'):
                st.warning('Stopped.' if runlib.stop(RUNS, sel)
                           else 'Could not stop it.')

        _follow()
