# -*- coding: utf-8 -*-
"""Run — launch a model run and follow it.  WP1b, steps 1b.6 / 1b.7.

The run is DETACHED. Streamlit reruns this script on every interaction, so a
model call in the page body would block the app for the ~11.5 minutes of a
coupled run and die on a browser refresh. Instead the run writes a log and a
status file, and this page polls them — so a run survives closing the tab.

Server mode is not a different mechanism: run Streamlit ON the server
(`streamlit run code/app/Home.py --server.address 0.0.0.0 --server.port 8501`)
and it launches there, next to the workspace. That is one setting, not an SSH
credential path.
"""

import os
import sys
import time

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg        # noqa: E402
import mm_paths                       # noqa: E402
from lib import checks as chk          # noqa: E402
from lib import panelui               # noqa: E402
from lib import runs as runlib        # noqa: E402

st.set_page_config(page_title='Run', page_icon='▶️', layout='wide')
CONFIG_DIR = os.path.join(CODE, 'configs')

st.title('▶️  8 — Run')
panelui.saved_note()

files = sorted(f for f in os.listdir(CONFIG_DIR) if f.endswith('.toml')) \
    if os.path.isdir(CONFIG_DIR) else []
if not files:
    st.error('No configuration in %s' % CONFIG_DIR)
    st.stop()

# THE ONE THE PANELS ARE EDITING, not a fixed favourite: this page used to
# open on lamata.toml whatever the panels had been pointed at, so the model
# that was launched need not have been the model that was filled in.
_panels = st.session_state.get('config_file', 'lamata.toml')
_default = _panels if _panels in files else (
    'lamata.toml' if 'lamata.toml' in files else files[0])
c1, c2 = st.columns([2, 3])
chosen = c1.selectbox('Configuration', files, index=files.index(_default))
cfg_path = os.path.join(CONFIG_DIR, chosen)
try:
    cfg = mcfg.load_run_config(cfg_path)
except mcfg.ConfigError as exc:
    st.error(str(exc))
    st.stop()

RUNS = cfg.ui.runs_dir or os.path.join(str(mm_paths.WS_ROOT), 'runs')
run_tag = c2.text_input('Run tag (optional)', value=cfg.meta.name or '')

# WHICH GRID THIS RUN WILL USE, stated rather than inferred. The grid is not
# committed by the sidebar save: the Grid panel's *Select this grid for the
# model* writes the [grid] settings AND promotes the mesh they produced to
# where the driver looks. So it is the one setting whose current value
# cannot be read off the panel you are standing on, and this page is the
# last place to notice that it is not the one you meant.
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
               'not the sidebar save, because it also promotes the mesh it '
               'built. A run rebuilds it if `[grid]` has changed since.'
               % _g['where'])
else:
    st.warning('No mesh is cached for this producer, so **the run will '
               'build it first** — minutes, on La Mata. Build it on the '
               'Grid panel and press *Select this grid for the model*.')

st.markdown('#### One-off overrides')
st.caption('`section.key=value`, one per line. Applied on top of the file and '
           'echoed by the driver — a setting that changes a run without '
           'appearing anywhere is how a multi-hour run gets wasted. '
           'There is deliberately no free-text command box: this page starts a '
           'process, so it only ever runs the validated entry point.')
raw = st.text_area('Overrides', value='', height=90,
                   placeholder='run.nsp=365\nsfr.enable=true')
overrides = [x.strip() for x in raw.splitlines() if x.strip()]

# Validate the overrides BEFORE launching, so a typo fails here in a second
# rather than in the run's first minute.
preview = None
if overrides:
    try:
        probe = mcfg.load_run_config(cfg_path)
        probe.apply_overrides(overrides, echo=False)
        preview = probe
        st.success('Overrides valid — resulting hash `%s`' % probe.config_hash())
    except mcfg.ConfigError as exc:
        st.error(str(exc))

# ASKED HERE, not edited into the file by hand. It was described in the
# schema and drawn by no panel, so the one field standing between a build and
# a coupled run could only be set in a text editor -- which is the thing the
# front-end exists to stop. It belongs on this page because this is where its
# effect is: blank builds and stops, set runs coupled.
st.markdown('#### MODFLOW 6 library')
_lc1, _lc2 = st.columns([2, 3])
with _lc1:
    _lib_edit = panelui.rows_form(cfg, (('paths.libmf6', None),), 'paths',
                                  columns=1)
libmf6 = str(_lib_edit.get('paths.libmf6', cfg.paths.libmf6) or '').strip()
with _lc2:
    st.markdown('')
    if not libmf6:
        st.warning('Blank: the run will stop after writing the MF6 files. '
                   'Type `auto` for a coupled run.')
    else:
        # SETTLED HERE, where it can still be corrected, rather than after
        # the build. The same resolver the driver uses, so what this page
        # says is what the run will do -- a folder is completed to the
        # library inside it, and mf6.exe is named as the wrong one of the two.
        _rl = runlib.import_runner()
        try:
            _lib = _rl.check_libmf6(_rl.resolve_libmf6(libmf6))
        except _rl.LibMF6Error as exc:
            st.error(str(exc))
        else:
            st.success('`%s`' % _lib)
    st.caption('`auto` finds it next to the other MODFLOW binaries and keeps '
               'the file portable to another machine. It is the shared '
               'LIBRARY — `libmf6.dll` — not `mf6.exe`: the coupler steps '
               'MODFLOW one stress period at a time through the API, which '
               'the executable cannot do. A folder is completed to the '
               'library inside it.')
panelui.remember(_lib_edit)

# ------------------------------------------------------------------ launch
# LAUNCH SAVES, THEN VALIDATES, THEN RUNS -- in that order, and it is the
# whole point of the button. The run reads the FILE, so anything the panels
# hold and have not written would be silently left out; and a configuration
# nobody has looked at is how a multi-hour run gets spent on a stress-period
# cap left over from a trial.
#
#   errors   -> nothing is launched, and the Validation panel is opened
#               on them
#   warnings -> nothing is launched YET, the Validation panel is opened so
#               they are SEEN,
#               and the run can be started from there
#   neither  -> it runs
#
# Panel 7 asks lib.checks exactly as this does, so the two cannot disagree.
st.markdown('---')
st.markdown('#### Launch')
_todo = panelui.unsaved_changes(cfg)
if _todo:
    st.info('%d unsaved change(s). Launch saves them first.' % len(_todo))

# The checks are run HERE as well as on the press, so the state of the button
# tells the truth before it is pressed: an error is a wall, not a surprise.
_now = chk.count_by_level(chk.collect_default(
    cfg, panelui.unsaved_switches(cfg)))
if _now[chk.ERROR]:
    st.error('**%d error(s).** This configuration cannot run — the '
             'Validation panel says what they are.' % _now[chk.ERROR])
elif _now[chk.WARNING]:
    st.warning('**%d warning(s).** Launch will take you to the Validation '
               'panel to read them; the run can be started from there.'
               % _now[chk.WARNING])
panelui.panel_link(panelui.VALIDATION_PAGE,
                   'Validation of the configuration', icon='🔎')

def _start():
    """Actually start the detached run."""
    try:
        run_id, payload = runlib.launch(
            cfg_path, RUNS, overrides=overrides, run_tag=run_tag or None,
            python_exe=(mm_paths.PYTHON_EXE
                        if os.path.exists(mm_paths.PYTHON_EXE)
                        else sys.executable))
    except Exception as exc:                            # noqa: BLE001
        st.error('Launch failed: %s' % exc)
    else:
        st.session_state['watch'] = run_id
        st.success('Launched `%s` (pid %s)' % (run_id, payload['pid']))


# Arriving from the Validation panel's own Launch: the checks were just run
# and looked at there, so they are not asked again -- the hash says it is the
# same configuration that was approved.
_approved = st.session_state.pop('validated', None) == cfg.config_hash()

can_launch = ((not overrides) or (preview is not None)) \
    and not _now[chk.ERROR]
if _approved:
    st.info('Validated — starting.')
    _start()
elif st.button('Launch', type='primary', disabled=not can_launch):
    if _todo:
        _applied, _why, _ok = panelui.save_now(cfg, cfg_path)
        if not _ok:
            st.error('NOT saved, so nothing was launched — the configuration '
                     'would be invalid:\n\n%s' % _why)
            st.stop()
        cfg = mcfg.load_run_config(cfg_path)
        st.session_state['__saved_note'] = (
            'Saved %d change(s) before launching.' % len(_applied))
    _counts = chk.count_by_level(
        chk.collect_default(cfg, panelui.unsaved_switches(cfg)))
    if _counts[chk.ERROR] or _counts[chk.WARNING]:
        # NOT launched. Panel 7 says what is wrong, in full, and offers the
        # launch again for the cases that are only warnings.
        if not panelui.go_to(panelui.VALIDATION_PAGE):
            st.error('This configuration has %d error(s) and %d warning(s) '
                     'and was NOT launched. Open the Validation of the '
                     'configuration panel to read them.'
                     % (_counts[chk.ERROR], _counts[chk.WARNING]))
            st.stop()
    _start()

st.markdown('---')
known = runlib.list_runs(RUNS)
if not known:
    st.info('No runs yet. `%s`' % RUNS)
    st.stop()

ids = [r['run_id'] for r in known]
watch = st.session_state.get('watch')
sel = st.selectbox('Follow', ids, index=ids.index(watch) if watch in ids else 0)
info = runlib.status(RUNS, sel)

c1, c2, c3, c4 = st.columns(4)
c1.metric('State', info.get('state', '?'))
c2.metric('Outcome', info.get('outcome', '—'))
c3.metric('Started', (info.get('started') or '')[-8:])
c4.metric('PID', info.get('pid', '—'))
st.code(' '.join(str(x) for x in info.get('cmd', [])), language='bash')

auto = st.checkbox('Auto-refresh every %d s' % cfg.ui.poll_secs,
                   value=(info.get('state') == 'running'))
# st.code, NOT st.text_area: a text_area given both `value` and `key` takes its
# content from session state after the first render, so the log froze at
# whatever it was when the page first drew it -- empty, for a run just started.
st.markdown('**run.log** (tail)')
st.code(runlib.log_tail(RUNS, sel, 400) or '(empty)', language='text',
        height=460)

if info.get('state') == 'running' and st.button('Stop this run'):
    st.warning('Stopped.' if runlib.stop(RUNS, sel) else 'Could not stop it.')

if auto and info.get('state') == 'running':
    time.sleep(max(1, int(cfg.ui.poll_secs)))
    st.rerun()

# The one save, as on every other panel: this page edits paths.libmf6.
panelui.sidebar_save(cfg, cfg_path)
