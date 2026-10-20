# -*- coding: utf-8 -*-
"""PEST++ jobs started from the PEST++ set-up page (WP7 front-end, phase
4): the SET-UP (code/tools/pest_lamata_setup.py) and the LAUNCH
(pest_lamata_run.py).

STREAMLIT-FREE, like runs.py, and for the same rule: a job is never run
inside a Streamlit callback. Each is a DETACHED process -- closing the
browser stops nothing -- started through code/tools/pest_job.py, which
records the tool's exit code, in

    <pest_ws>/jobs/<YYYYmmddHHMMSS>_<kind>/
        job.log       the tool's output, live
        status.json   written here at start: kind, pid, command, started
        result.json   written by the job when the tool returns: its exit
                      code (Windows cannot give it to another process)

ONE PEST JOB AT A TIME: a set-up rewrites the template every agent copies
and a launch runs 20 agents on it, so a second job is refused while one
runs, or while pestpp-ies runs at all. The launch is offered only after a
set-up that SUCCEEDED -- its exit code 0 and a template written by it.
The tools run from the checkout (mm_paths.REPO), as they always have:
the set-up snapshots that code for the workers.
"""

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

__all__ = ['KINDS', 'JobBusy', 'jobs_dir', 'launch', 'status', 'jobs',
           'latest', 'active', 'setup_ready', 'log_tail', 'stop']

KINDS = {'setup': 'pest_lamata_setup.py', 'launch': 'pest_lamata_run.py'}
_STATUS, _RESULT, _LOG = 'status.json', 'result.json', 'job.log'


class JobBusy(RuntimeError):
    """A PEST job (or pestpp-ies) is running: one at a time."""


def jobs_dir(pest_ws):
    return Path(pest_ws) / 'jobs'


def _read(p):
    try:
        with open(p, encoding='utf-8') as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _write(p, payload):
    tmp = Path(str(p) + '.tmp')
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=1, default=str)
    os.replace(tmp, p)


def _pid_alive(pid):
    from lib import runs
    return runs._pid_alive(pid)


def status(pest_ws, job_id):
    """A job's status: what was started, and -- once it ended -- its exit
    code and outcome ('ok', 'failed', or 'unknown' for a job that died
    without writing one, killed or stopped)."""
    d = jobs_dir(pest_ws) / job_id
    st = _read(d / _STATUS)
    if st is None:
        return None
    res = _read(d / _RESULT)
    if res is not None:
        st.update(res, state='finished',
                  outcome='ok' if int(res.get('returncode', 1)) == 0
                  else 'failed')
    elif st.get('state') == 'running' and not _pid_alive(st.get('pid')):
        st.update(state='finished', outcome=st.get('outcome') or 'unknown')
    return st


def jobs(pest_ws, kind=None):
    """The jobs, newest first."""
    root = jobs_dir(pest_ws)
    if not root.is_dir():
        return []
    out = []
    for d in sorted((x for x in root.iterdir() if x.is_dir()), reverse=True):
        st = status(pest_ws, d.name)
        if st and (kind is None or st.get('kind') == kind):
            out.append(st)
    return out


def latest(pest_ws, kind=None):
    j = jobs(pest_ws, kind)
    return j[0] if j else None


def active(pest_ws):
    return [j for j in jobs(pest_ws) if j.get('state') == 'running']


def setup_ready(pest_ws):
    """``(ready, why)``: may a calibration be launched? Only after the
    LATEST set-up succeeded and left the template it writes."""
    s = latest(pest_ws, 'setup')
    if s is None:
        return False, 'no set-up has been run from this page yet'
    if s.get('state') == 'running':
        return False, 'the set-up is still running'
    if s.get('outcome') != 'ok':
        return False, ('the last set-up (%s) did not succeed -- see its log'
                       % s.get('job_id'))
    pst = Path(pest_ws) / 'template' / 'lamata.pst'
    if not pst.exists():
        return False, 'the set-up left no template (%s)' % pst
    try:
        t0 = time.mktime(time.strptime(s.get('started', ''),
                                       '%Y-%m-%d %H:%M:%S'))
    except ValueError:
        t0 = 0.0
    if pst.stat().st_mtime < t0 - 1.0:
        return False, 'the template is older than the last set-up'
    return True, 'set-up %s succeeded' % s.get('job_id')


def launch(kind, config_path, pest_ws, python_exe=None, tools_dir=None,
           running=None):
    """Start a detached PEST job. Returns ``(job_id, status)``. Raises
    JobBusy while another job or pestpp-ies runs (``running``: the
    pestpp-ies check, injectable for tests)."""
    if kind not in KINDS:
        raise ValueError('kind must be one of %s' % ', '.join(KINDS))
    busy = active(pest_ws)
    if busy:
        raise JobBusy('job %s is still running' % busy[0].get('job_id'))
    if running is None:
        from lib import pest_progress
        running = pest_progress.pestpp_running
    if running():
        raise JobBusy('pestpp-ies is running on this machine')
    if not os.path.exists(str(config_path)):
        raise FileNotFoundError('configuration not found: %s' % config_path)
    if tools_dir is None:
        import mm_paths
        tools_dir = Path(mm_paths.REPO) / 'code' / 'tools'
    tools_dir = Path(tools_dir)
    job_id = time.strftime('%Y%m%d%H%M%S') + '_' + kind
    d = jobs_dir(pest_ws) / job_id
    d.mkdir(parents=True, exist_ok=False)
    python_exe = python_exe or sys.executable
    cmd = [str(python_exe), '-u', str(tools_dir / 'pest_job.py'), str(d),
           str(tools_dir / KINDS[kind]), '--config', str(config_path)]
    kw = {}
    if os.name == 'nt':
        kw['creationflags'] = (subprocess.CREATE_NEW_PROCESS_GROUP
                               | getattr(subprocess, 'DETACHED_PROCESS', 0))
    else:
        kw['start_new_session'] = True
    log = open(d / _LOG, 'w', encoding='utf-8', errors='replace')
    proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                            cwd=str(pest_ws), **kw)
    st = {'job_id': job_id, 'kind': kind, 'pid': proc.pid, 'cmd': cmd,
          'config': str(config_path), 'state': 'running',
          'started': time.strftime('%Y-%m-%d %H:%M:%S')}
    _write(d / _STATUS, st)
    return job_id, st


def log_tail(pest_ws, job_id, lines=40):
    p = jobs_dir(pest_ws) / job_id / _LOG
    if not p.exists():
        return ''
    size = p.stat().st_size
    with open(p, 'rb') as fh:
        fh.seek(max(0, size - max(65536, lines * 400)))
        raw = fh.read().decode('utf-8', errors='replace')
    # a line rewritten in place (PANTHER's progress line, every second, by
    # carriage returns) is shown as it stands now, once -- not as the 40
    # last rewrites
    out = [ln.rstrip('\r').split('\r')[-1]
           for ln in raw.replace('\r\n', '\n').split('\n')]
    return '\n'.join(out[-lines:])


def stop(pest_ws, job_id):
    """Stop a running job and everything it started (the agents of a
    launch with it). Returns True if the stop was sent."""
    st = status(pest_ws, job_id)
    if not st or st.get('state') != 'running':
        return False
    pid = st.get('pid')
    try:
        if os.name == 'nt':
            subprocess.run(['taskkill', '/PID', str(pid), '/T', '/F'],
                           capture_output=True, timeout=30)
        else:
            os.killpg(os.getpgid(int(pid)), signal.SIGTERM)
    except Exception:                                   # noqa: BLE001
        return False
    st.update(state='finished', outcome='stopped',
              finished=time.strftime('%Y-%m-%d %H:%M:%S'))
    _write(jobs_dir(pest_ws) / job_id / _STATUS, st)
    return True
