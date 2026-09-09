# -*- coding: utf-8 -*-
"""The run registry: launch a model run and follow it.  WP1b, step 1b.6.

STREAMLIT-FREE, like loaders.py, so it can be tested in the model environment.

THE RULE THIS MODULE EXISTS TO ENFORCE: a model run is NEVER executed inside a
Streamlit callback. Streamlit reruns the whole script on every widget
interaction, and a coupled La Mata run is ~11.5 minutes (an IES run is hours).
An in-line call would block the app for the duration and die on a browser
refresh. So a run is launched as a DETACHED subprocess that writes

    <runs_dir>/<run_id>/run.log       the driver's stdout+stderr, live
    <runs_dir>/<run_id>/status.json   pid, command, state, timestamps

and the page polls those. A run therefore survives closing the browser, and its
log can be read afterwards by anyone.
"""

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

__all__ = ['launch', 'status', 'log_tail', 'list_runs', 'stop', 'run_dir']

_STATUS = 'status.json'
_LOG = 'run.log'


def run_dir(runs_dir, run_id):
    return Path(runs_dir) / run_id


def _write_status(d, payload):
    tmp = Path(d) / (_STATUS + '.tmp')
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, default=str)
    os.replace(tmp, Path(d) / _STATUS)     # atomic: a reader never sees half a file


def launch(config_path, runs_dir, overrides=None, run_tag=None,
           python_exe=None, driver=None, cwd=None, env=None):
    """Start a detached model run. Returns (run_id, status dict).

    The command is built from a VALIDATED configuration path and a list of
    ``section.key=value`` overrides -- never from free text. The page must not
    offer a command box: this function executes a process.
    """
    config_path = str(config_path)
    if not os.path.exists(config_path):
        raise FileNotFoundError('configuration not found: %s' % config_path)
    for item in overrides or ():
        if '=' not in item or item.strip().startswith('-'):
            raise ValueError('override must be section.key=value, got %r' % item)

    here = Path(__file__).resolve()
    code_dir = here.parents[2]                       # .../code
    driver = driver or str(code_dir / 'tests' / 'run_lamata_mf6.py')
    python_exe = python_exe or sys.executable
    cwd = str(cwd or code_dir.parent)                # the repo root

    run_id = time.strftime('%Y%m%d%H%M%S')
    if run_tag:
        run_id += '_' + str(run_tag).strip().replace(' ', '_')
    d = run_dir(runs_dir, run_id)
    d.mkdir(parents=True, exist_ok=True)

    cmd = [python_exe, driver, '--config', config_path]
    for item in overrides or ():
        cmd += ['--set', item]
    if run_tag:
        cmd += ['--run-tag', str(run_tag)]

    log_path = d / _LOG
    logfh = open(log_path, 'w', encoding='utf-8', errors='replace')
    # Detach so the run outlives the Streamlit session that started it.
    kwargs = {}
    if os.name == 'nt':
        kwargs['creationflags'] = (subprocess.CREATE_NEW_PROCESS_GROUP
                                   | getattr(subprocess, 'DETACHED_PROCESS', 0))
    else:
        kwargs['start_new_session'] = True
    proc = subprocess.Popen(cmd, stdout=logfh, stderr=subprocess.STDOUT,
                            cwd=cwd, env=env, **kwargs)
    payload = {'run_id': run_id, 'pid': proc.pid, 'state': 'running',
               'cmd': cmd, 'cwd': cwd, 'config': config_path,
               'overrides': list(overrides or ()), 'run_tag': run_tag,
               'started': time.strftime('%Y-%m-%d %H:%M:%S'),
               'finished': None, 'returncode': None}
    _write_status(d, payload)
    return run_id, payload


def _pid_alive(pid):
    if pid is None:
        return False
    try:
        if os.name == 'nt':
            out = subprocess.run(
                ['tasklist', '/FI', 'PID eq %d' % int(pid), '/NH'],
                capture_output=True, text=True, timeout=10).stdout
            return str(pid) in out
        os.kill(int(pid), 0)
        return True
    except Exception:
        return False


def status(runs_dir, run_id, refresh=True):
    """Read a run's status, refreshing 'running' against the real process."""
    p = run_dir(runs_dir, run_id) / _STATUS
    if not p.exists():
        return None
    with open(p, encoding='utf-8') as fh:
        st = json.load(fh)
    if refresh and st.get('state') == 'running' and not _pid_alive(st.get('pid')):
        # The process is gone. We cannot recover its exit code from another
        # process on Windows, so record what we know and say so plainly.
        st['state'] = 'finished'
        st['finished'] = time.strftime('%Y-%m-%d %H:%M:%S')
        tail = log_tail(runs_dir, run_id, 40)
        st['outcome'] = ('failed' if _looks_failed(tail) else 'completed')
        _write_status(run_dir(runs_dir, run_id), st)
    return st


_FAIL_MARKERS = ('Traceback (most recent call last)', 'SystemExit',
                 'CONFIG ERROR', 'Error:', 'ERROR:', 'did not converge')


def _looks_failed(text):
    return any(m in text for m in _FAIL_MARKERS)


def log_tail(runs_dir, run_id, lines=200):
    """Last ``lines`` of a run's log. Reads only the tail of the file."""
    p = run_dir(runs_dir, run_id) / _LOG
    if not p.exists():
        return ''
    size = p.stat().st_size
    want = min(size, max(4096, lines * 200))
    with open(p, 'rb') as fh:
        fh.seek(size - want)
        raw = fh.read().decode('utf-8', errors='replace')
    return '\n'.join(raw.splitlines()[-lines:])


def list_runs(runs_dir, limit=50):
    """Known runs, newest first."""
    root = Path(runs_dir)
    if not root.is_dir():
        return []
    out = []
    for d in sorted((x for x in root.iterdir() if x.is_dir()), reverse=True):
        st = status(runs_dir, d.name, refresh=False)
        if st:
            out.append(st)
        if len(out) >= limit:
            break
    return out


def stop(runs_dir, run_id):
    """Ask a running model to stop. Returns True if a signal was delivered."""
    st = status(runs_dir, run_id, refresh=False)
    if not st or st.get('state') != 'running':
        return False
    pid = st.get('pid')
    try:
        if os.name == 'nt':
            subprocess.run(['taskkill', '/PID', str(pid), '/T', '/F'],
                           capture_output=True, timeout=15)
        else:
            os.killpg(os.getpgid(int(pid)), signal.SIGTERM)
    except Exception:
        return False
    st['state'] = 'finished'
    st['outcome'] = 'stopped'
    st['finished'] = time.strftime('%Y-%m-%d %H:%M:%S')
    _write_status(run_dir(runs_dir, run_id), st)
    return True
