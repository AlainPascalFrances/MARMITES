# -*- coding: utf-8 -*-
"""PEST++ jobs from the PEST++ set-up page (WP7 front-end, phase 4; user,
2026-10-20): the set-up and the launch as DETACHED jobs with their exit
code recorded, one at a time, the launch offered only after a set-up that
succeeded. Fake tools stand in for pest_lamata_setup / _run."""

import json
import os
import shutil
import sys
import time

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for p in (CODE, os.path.join(CODE, 'app')):
    if p not in sys.path:
        sys.path.insert(0, p)

from lib import pest_jobs                                       # noqa: E402

NOT_RUNNING = staticmethod(lambda: False)


def _tools(tmp_path, setup_body, run_body='print("run")'):
    """A tools folder: the real job wrapper, fake tools."""
    d = tmp_path / 'tools'
    d.mkdir()
    shutil.copy(os.path.join(CODE, 'tools', 'pest_job.py'), d)
    (d / 'pest_lamata_setup.py').write_text(setup_body, encoding='utf-8')
    (d / 'pest_lamata_run.py').write_text(run_body, encoding='utf-8')
    return d


def _wait(ws, job_id, timeout=300.0):
    t0 = time.time()
    while time.time() - t0 < timeout:
        st = pest_jobs.status(ws, job_id)
        if st and st['state'] != 'running':
            return st
        time.sleep(0.2)
    raise AssertionError('job %s still running' % job_id)


SETUP_OK = ('import os, sys\n'
            'print("setting up", sys.argv[1:], flush=True)\n'
            'os.makedirs("template", exist_ok=True)\n'
            'open(os.path.join("template", "lamata.pst"), "w").write("pst")\n')
SETUP_FAIL = 'import sys\nprint("no base")\nsys.exit(3)\n'


def test_a_set_up_that_succeeds_opens_the_launch(tmp_path):
    ws = tmp_path / 'pest'
    ws.mkdir()
    cfg = tmp_path / 'c.toml'
    cfg.write_text('[meta]\n')
    tools = _tools(tmp_path, SETUP_OK)
    assert pest_jobs.setup_ready(str(ws))[0] is False
    job_id, st = pest_jobs.launch('setup', str(cfg), str(ws),
                                  python_exe=sys.executable,
                                  tools_dir=str(tools),
                                  running=lambda: False)
    assert st['state'] == 'running' and job_id.endswith('_setup')
    done = _wait(str(ws), job_id)
    assert done['outcome'] == 'ok' and done['returncode'] == 0
    assert 'setting up' in pest_jobs.log_tail(str(ws), job_id)
    assert '--config' in pest_jobs.log_tail(str(ws), job_id)
    ready, why = pest_jobs.setup_ready(str(ws))
    assert ready, why


def test_a_set_up_that_fails_keeps_the_launch_closed(tmp_path):
    ws = tmp_path / 'pest'
    ws.mkdir()
    cfg = tmp_path / 'c.toml'
    cfg.write_text('[meta]\n')
    tools = _tools(tmp_path, SETUP_FAIL)
    job_id, _st = pest_jobs.launch('setup', str(cfg), str(ws),
                                   python_exe=sys.executable,
                                   tools_dir=str(tools),
                                   running=lambda: False)
    done = _wait(str(ws), job_id)
    assert done['outcome'] == 'failed' and done['returncode'] == 3
    ready, why = pest_jobs.setup_ready(str(ws))
    assert not ready and 'did not succeed' in why


def test_one_job_at_a_time_and_a_stop(tmp_path):
    ws = tmp_path / 'pest'
    ws.mkdir()
    cfg = tmp_path / 'c.toml'
    cfg.write_text('[meta]\n')
    tools = _tools(tmp_path, 'import time\ntime.sleep(60)\n')
    job_id, _st = pest_jobs.launch('setup', str(cfg), str(ws),
                                   python_exe=sys.executable,
                                   tools_dir=str(tools),
                                   running=lambda: False)
    time.sleep(1.0)
    with pytest.raises(pest_jobs.JobBusy):
        pest_jobs.launch('launch', str(cfg), str(ws), tools_dir=str(tools),
                         running=lambda: False)
    assert pest_jobs.stop(str(ws), job_id)
    st = pest_jobs.status(str(ws), job_id)
    assert st['state'] == 'finished' and st['outcome'] == 'stopped'
    assert not pest_jobs._pid_alive(st['pid'])
    assert not pest_jobs.setup_ready(str(ws))[0]
    # pestpp-ies running elsewhere refuses a job too
    with pytest.raises(pest_jobs.JobBusy, match='pestpp-ies'):
        pest_jobs.launch('setup', str(cfg), str(ws), tools_dir=str(tools),
                         running=lambda: True)


def test_a_template_older_than_the_set_up_is_not_ready(tmp_path):
    ws = tmp_path / 'pest'
    (ws / 'template').mkdir(parents=True)
    pst = ws / 'template' / 'lamata.pst'
    pst.write_text('old')
    old = time.time() - 3600
    os.utime(pst, (old, old))
    cfg = tmp_path / 'c.toml'
    cfg.write_text('[meta]\n')
    tools = _tools(tmp_path, 'print("did nothing")\n')
    job_id, _st = pest_jobs.launch('setup', str(cfg), str(ws),
                                   python_exe=sys.executable,
                                   tools_dir=str(tools),
                                   running=lambda: False)
    _wait(str(ws), job_id)
    ready, why = pest_jobs.setup_ready(str(ws))
    assert not ready and 'older' in why


def test_a_line_rewritten_in_place_is_shown_once(tmp_path):
    """PANTHER rewrites its progress line every second with carriage
    returns: the page shows it as it stands, not its last 40 versions."""
    ws = tmp_path / 'pest'
    d = ws / 'jobs' / '20261020081102_launch'
    d.mkdir(parents=True)
    (d / 'job.log').write_bytes(
        b'pestpp-ies, noptmax 3\nPANTHER progress\n'
        + b''.join(b'10/20 08:11:%02d runs(C%d)\r' % (s, s)
                   for s in range(40, 60)) + b'10/20 08:12:00 runs(C20)')
    tail = pest_jobs.log_tail(str(ws), '20261020081102_launch', 40)
    assert tail.splitlines() == ['pestpp-ies, noptmax 3', 'PANTHER progress',
                                 '10/20 08:12:00 runs(C20)']


def test_a_stop_that_does_not_kill_leaves_the_job_running(tmp_path,
                                                          monkeypatch):
    """2026-10-20: taskkill was refused ('Access is denied') and the job
    was marked stopped all the same, its PEST++ processes running on. Now
    the stop says so, and the job stays 'running' for another try."""
    ws = tmp_path / 'pest'
    d = ws / 'jobs' / 'j1'
    d.mkdir(parents=True)
    st0 = {'job_id': 'j1', 'kind': 'launch', 'pid': 4242,
           'state': 'running', 'started': '2026-10-20 08:11:02'}
    (d / pest_jobs._STATUS).write_text(json.dumps(st0))

    class _R:
        stdout = ''
        stderr = ('ERROR: The process with PID 4242 (child process of PID '
                  '11) could not be terminated.\nReason: Access is denied.')

    monkeypatch.setattr(pest_jobs.subprocess, 'run', lambda *a, **k: _R())
    monkeypatch.setattr(pest_jobs.os, 'killpg', lambda *a: None,
                        raising=False)
    monkeypatch.setattr(pest_jobs.os, 'getpgid', lambda p: p, raising=False)
    monkeypatch.setattr(pest_jobs, '_pid_alive', lambda pid: True)
    with pytest.raises(pest_jobs.StopFailed, match='still running') as e:
        pest_jobs.stop(str(ws), 'j1', wait=0.2)
    if os.name == 'nt':
        assert 'Access is denied' in str(e.value)
    assert pest_jobs.status(str(ws), 'j1')['state'] == 'running'
