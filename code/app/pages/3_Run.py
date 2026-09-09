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
from lib import runs as runlib        # noqa: E402

st.set_page_config(page_title='Run', page_icon='▶️', layout='wide')
CONFIG_DIR = os.path.join(CODE, 'configs')

st.title('Run')

files = sorted(f for f in os.listdir(CONFIG_DIR) if f.endswith('.toml')) \
    if os.path.isdir(CONFIG_DIR) else []
if not files:
    st.error('No configuration in %s' % CONFIG_DIR)
    st.stop()

c1, c2 = st.columns([2, 3])
chosen = c1.selectbox('Configuration', files,
                      index=files.index('lamata.toml') if 'lamata.toml' in files else 0)
cfg_path = os.path.join(CONFIG_DIR, chosen)
try:
    cfg = mcfg.load_run_config(cfg_path)
except mcfg.ConfigError as exc:
    st.error(str(exc))
    st.stop()

RUNS = cfg.ui.runs_dir or os.path.join(str(mm_paths.WS_ROOT), 'runs')
run_tag = c2.text_input('Run tag (optional)', value=cfg.meta.name or '')

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

libmf6 = (cfg.paths.libmf6 or '').strip()
if not libmf6:
    st.warning('`paths.libmf6` is blank, so the run will stop after writing the '
               'MF6 files. Set it to `auto` (use `mm_paths.LIBMF6`) or to a path '
               'for a coupled run.')

can_launch = (not overrides) or (preview is not None)
if st.button('Launch', type='primary', disabled=not can_launch):
    try:
        run_id, st_payload = runlib.launch(
            cfg_path, RUNS, overrides=overrides, run_tag=run_tag or None,
            python_exe=(mm_paths.PYTHON_EXE if os.path.exists(mm_paths.PYTHON_EXE)
                        else sys.executable))
    except Exception as exc:
        st.error('Launch failed: %s' % exc)
    else:
        st.session_state['watch'] = run_id
        st.success('Launched `%s` (pid %s)' % (run_id, st_payload['pid']))

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
st.text_area('run.log (tail)', value=runlib.log_tail(RUNS, sel, 400),
             height=460, key='log_' + sel)

if info.get('state') == 'running' and st.button('Stop this run'):
    st.warning('Stopped.' if runlib.stop(RUNS, sel) else 'Could not stop it.')

if auto and info.get('state') == 'running':
    time.sleep(max(1, int(cfg.ui.poll_secs)))
    st.rerun()
