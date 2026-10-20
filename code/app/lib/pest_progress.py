# -*- coding: utf-8 -*-
"""A PEST++ run's progress, read from what pestpp-ies's master writes.
WP7 front-end, phase 1 (docs/MM-MF6/WP7_PEST_frontend_plan.md).

STREAMLIT-FREE, like runs.py, so it is tested in the model environment.

READ, NEVER WRITTEN: ``<pest_ws>/master/``

    <case>.rec        the options, where the run is (stage markers), the
                      lambda tests' phi, "analysis complete" at the end
    <case>.rmr        the run manager: one batch per "running model N
                      times"; each run sent to an agent, then received,
                      failed, or killed as overdue
    <case>.phi.actual.csv, <case>.phi.group.csv   phi per iteration
    report_<k>/       the iteration reports (tools/pest_lamata_report.py)

and, for a failed run's cause, the tail of the agent's model log
(``<pest_ws>/workers/<agent>/model.log``) while that agent has had no run
since.

Run ids start again at 0 in every batch, so a run is (batch, run id). A
run can be sent twice (PANTHER sends an overdue run to a second agent);
it is completed as soon as one copy comes back. A lambda test's runs say
their lambda (``lambda:0.01 scale_fac:1``); that is what names the batch.
"""

import datetime as dt
import math
import os
import re
import subprocess
from pathlib import Path

__all__ = ['parse_rmr', 'parse_rec', 'read_phi', 'failure_cause',
           'progress', 'pestpp_running', 'master_dir', 'reports']

_TS = re.compile(r'^(\d\d/\d\d/\d\d \d\d:\d\d:\d\d)->')
_BATCH = re.compile(r'running model (\d+) times')
_SEND = re.compile(r'Sending run_id:(\d+) to:(\S+) group_id:(\d+)')
_ITER = re.compile(r'\biteration:(-?\d+)')
_REAL = re.compile(r'\brealization:(\S+)')
_LAMBDA = re.compile(r'\blambda:(\S+)')
_RECV = re.compile(r'run_id:(\d+) received from:(\S+)\s+run_time:([\d.]+)')
_FAIL = re.compile(r'run_id:(\d+) failed on agent:(\S+)')
_KILL = re.compile(r'sending kill request\. reason:([^.]+)\. '
                   r'duration:([\d.]+), avg:([\d.]+) run_id:(\d+)')
_NUM = r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'


def master_dir(pest_ws):
    return Path(pest_ws) / 'master'


def _when(s):
    return dt.datetime.strptime(s, '%m/%d/%y %H:%M:%S')


def _agent(s):
    """'HOST$X:\\...\\workers\\worker_3' -> 'worker_3'."""
    return re.split(r'[\\/$]', s)[-1]


def parse_rmr(text):
    """The run manager's record as a list of batches, in order:
    {planned, iteration, stage, lambdas, first, last, runs: {run_id: {real,
    lambda, agent, sent, end, state, time}}}; state is completed, failed,
    killed or running (sent and not back -- 'unfinished' once the master
    has stopped)."""
    batches, cur, ts = [], None, None
    for line in text.splitlines():
        m = _TS.match(line)
        if m:
            ts = _when(m.group(1))
        b = _BATCH.search(line)
        if b:
            cur = {'planned': int(b.group(1)), 'iteration': None,
                   'lambdas': [], 'first': None, 'last': ts, 'runs': {}}
            batches.append(cur)
            continue
        if cur is None:
            continue
        s = _SEND.search(line)
        if s:
            rid = int(s.group(1))
            it, real = _ITER.search(line), _REAL.search(line)
            lam = _LAMBDA.search(line)
            r = cur['runs'].setdefault(rid, {
                'real': real.group(1) if real else '', 'sent': ts,
                'lambda': float(lam.group(1)) if lam else None,
                'end': None, 'time': None, 'state': 'running'})
            r['agent'] = _agent(s.group(2))
            if cur['iteration'] is None and it:
                cur['iteration'] = int(it.group(1))
            if lam and float(lam.group(1)) not in cur['lambdas']:
                cur['lambdas'].append(float(lam.group(1)))
            cur['first'] = cur['first'] or ts
            cur['last'] = ts
            continue
        s = _RECV.search(line)
        if s:
            r = cur['runs'].get(int(s.group(1)))
            if r is not None:
                r.update(state='completed', time=float(s.group(3)), end=ts,
                         agent=_agent(s.group(2)))
                cur['last'] = ts
            continue
        s = _FAIL.search(line)
        if s:
            r = cur['runs'].get(int(s.group(1)))
            if r is not None and r['state'] != 'completed':
                r.update(state='failed', end=ts, agent=_agent(s.group(2)))
                cur['last'] = ts
            continue
        s = _KILL.search(line)
        if s:
            r = cur['runs'].get(int(s.group(4)))
            if r is not None and r['state'] == 'running':
                r.update(state='killed', end=ts, reason=s.group(1).strip())
                cur['last'] = ts
    # what each batch is: a lambda test says its lambdas; otherwise the
    # initial ensemble (iteration 0) or an iteration's remaining
    # realisations
    for b in batches:
        b['stage'] = ('lambda test' if b['lambdas'] else
                      'initial ensemble' if b['iteration'] in (None, 0)
                      else 'remaining realisations')
    return batches


def parse_rec(text):
    """The record file: options, where the run is, the lambda tests."""
    def opt(pat, cast=float):
        m = re.search(pat, text, re.M)
        return cast(m.group(1)) if m else None
    out = {'noptmax': opt(r'^\s*noptmax\s*=\s*(-?\d+)', int),
           'nreal': opt(r'^ies_num_reals:\s*(\d+)', int),
           'overdue_giveup_fac': opt(r'^overdue_giveup_fac:\s*' + _NUM),
           'lambda_mults': [float(x) for x in re.findall(
               _NUM, (re.search(r'^ies_lambda_mults:\s*(.*)$', text, re.M)
                      or [None, ''])[1])],
           'finished': 'pestpp-ies analysis complete' in text}
    its = re.findall(r'starting solve for iteration:\s*(\d+)', text)
    out['iteration'] = int(its[-1]) if its else 0
    marks = [('initialization complete', 'initial ensemble done'),
             ('starting solve for iteration', 'computing the upgrades'),
             ('running upgrade ensembles', 'lambda test'),
             ('evaluating upgrade ensembles', 'evaluating the lambda test'),
             ('running remaining realizations', 'remaining realisations'),
             ('pestpp-ies analysis complete', 'finished')]
    pos = [(text.rfind(k), v) for k, v in marks if text.rfind(k) >= 0]
    out['stage'] = max(pos)[1] if pos else 'initial ensemble'
    # the lambda tests: each 'phi summary for lambda, scale fac:L , 1' is
    # followed by its measured and actual phi lines
    tests, it = [], 0
    lines = text.splitlines()
    for i, line in enumerate(lines):
        m = re.search(r'starting solve for iteration:\s*(\d+)', line)
        if m:
            it = int(m.group(1))
            continue
        m = re.search(r'phi summary for (best )?lambda, scale fac:\s*'
                      + _NUM, line)
        if not m:
            continue
        row = {'iteration': it, 'lambda': float(m.group(2)),
               'best': bool(m.group(1))}
        for nxt in lines[i + 1:i + 6]:
            a = re.match(r'\s*actual\s+' + r'\s+'.join([_NUM] * 4), nxt)
            if a:
                row.update(mean=float(a.group(1)), sd=float(a.group(2)),
                           min=float(a.group(3)), max=float(a.group(4)))
                break
        tests.append(row)
    out['lambda_tests'] = tests
    for k, pat in (('started_at', r'started at (.+)$'),
                   ('finished_at', r'finished at (.+)$')):
        m = re.search(pat, text, re.M)
        out[k] = m.group(1).strip() if m else None
    return out


def read_phi(master, case):
    """(per iteration: total_runs, mean, sd, min, max, n; the realisations'
    phi long: iteration, realization, phi; the group means per iteration)
    as DataFrames, or None where the file is not there yet."""
    import pandas as pd
    out = [None, None, None]
    fa = Path(master) / ('%s.phi.actual.csv' % case)
    if fa.exists():
        a = pd.read_csv(fa)
        reals = [c for c in a.columns if c not in (
            'iteration', 'total_runs', 'mean', 'standard_deviation', 'min',
            'max')]
        t = a[['iteration', 'total_runs', 'mean', 'standard_deviation',
               'min', 'max']].copy()
        t['realisations'] = a[reals].notna().sum(axis=1).to_numpy()
        out[0] = t
        long = a.melt(id_vars=['iteration'], value_vars=reals,
                      var_name='realization', value_name='phi').dropna()
        out[1] = long
    fg = Path(master) / ('%s.phi.group.csv' % case)
    if fg.exists():
        g = pd.read_csv(fg)
        groups = [c for c in g.columns if c not in (
            'iteration', 'total_runs', 'obs_realization', 'par_realization')]
        out[2] = g.groupby('iteration')[groups].mean().reset_index()
    return tuple(out)


def failure_cause(pest_ws, agent, lines=80):
    """Why a run failed, from the end of its model log: the newest one
    the forward run KEPT for that agent (failed_runs/, 2026-10-20), else
    the agent's own model.log -- which its next run overwrites within
    seconds. The error and the lines that explain it; '' when there is
    nothing to read."""
    kept = sorted((Path(pest_ws) / 'failed_runs').glob(
        '*_%s_model.log' % agent))
    fn = kept[-1] if kept else Path(pest_ws) / 'workers' / agent / 'model.log'
    try:
        with open(fn, 'rb') as fh:
            fh.seek(0, 2)
            fh.seek(max(fh.tell() - 20000, 0))
            tail = fh.read().decode('utf-8', 'replace').splitlines()[-lines:]
    except OSError:
        return ''
    for i in range(len(tail) - 1, -1, -1):
        line = tail[i].strip()
        if re.match(r'^[\w.]*(Error|Exception|NonConvergence)\b.*:', line):
            more = [x.strip() for x in tail[i + 1:i + 6]
                    if x.strip().startswith('-')]
            return ' '.join([line] + more)
    return ''


def pestpp_running():
    """Whether a pestpp-ies process is alive on this machine."""
    try:
        if os.name == 'nt':
            out = subprocess.run(['tasklist', '/FI',
                                  'IMAGENAME eq pestpp-ies.exe', '/NH'],
                                 capture_output=True, text=True,
                                 timeout=10).stdout
            return 'pestpp-ies' in out.lower()
        out = subprocess.run(['pgrep', '-f', 'pestpp-ies'],
                             capture_output=True, text=True, timeout=10)
        return bool(out.stdout.strip())
    except Exception:
        return False


def reports(master):
    """The iteration reports: [(k, folder)], in order."""
    out = []
    for d in Path(master).glob('report_*'):
        try:
            out.append((int(d.name.split('_')[1]), d))
        except (IndexError, ValueError):
            continue
    return sorted(out)


def _case(master):
    pst = sorted(Path(master).glob('*.pst'))
    return pst[0].stem if pst else None


def _read(fn):
    try:
        return Path(fn).read_text(encoding='utf-8', errors='replace')
    except OSError:
        return ''


def progress(pest_ws, now=None, running=None):
    """Everything the progress page shows, as plain data; None when
    ``<pest_ws>/master`` holds no run. ``running``: whether pestpp-ies is
    alive (asked of the machine when None)."""
    master = master_dir(pest_ws)
    case = _case(master) if master.is_dir() else None
    if case is None:
        return None
    now = now or dt.datetime.now()
    rec = parse_rec(_read(master / ('%s.rec' % case)))
    batches = parse_rmr(_read(master / ('%s.rmr' % case)))
    alive = pestpp_running() if running is None else running
    state = ('finished' if rec['finished']
             else 'running' if alive else 'stopped')
    rows, times, agents, failures = [], [], [], []
    last_of = {}                     # agent -> (batch, run id) of its last run
    for bi, b in enumerate(batches):
        for rid, r in b['runs'].items():
            last_of[r['agent']] = max(last_of.get(r['agent'], (-1, -1)),
                                      (bi, rid))
    for bi, b in enumerate(batches):
        runs = b['runs'].values()
        n = {k: sum(r['state'] == k for r in runs)
             for k in ('completed', 'failed', 'killed', 'running')}
        t = [r['time'] for r in runs if r['time'] is not None]
        times += t
        live = state == 'running' and bi == len(batches) - 1
        # sent and not back: running while the master runs, unfinished once
        # it has stopped
        rows.append({'iteration': b['iteration'], 'stage': b['stage'],
                     'planned': b['planned'], 'completed': n['completed'],
                     'failed': n['failed'], 'killed (overdue)': n['killed'],
                     ('running' if state == 'running' else 'unfinished'):
                     n['running'],
                     'mean run (min)': (sum(t) / len(t)) if t else None,
                     'longest run (min)': max(t) if t else None,
                     'started': b['first'], 'last event': b['last']})
        for rid, r in b['runs'].items():
            if r['state'] == 'running' and live:
                agents.append({'agent': r['agent'], 'run': rid,
                               'realisation': r['real'], 'started': r['sent'],
                               'elapsed (min)': (now - r['sent'])
                               .total_seconds() / 60.0})
            if r['state'] in ('failed', 'killed'):
                cause = (r.get('reason', '') if r['state'] == 'killed'
                         else (failure_cause(pest_ws, r['agent'])
                               if last_of.get(r['agent']) == (bi, rid)
                               else ''))
                failures.append({'iteration': b['iteration'],
                                 'stage': b['stage'], 'run': rid,
                                 'realisation': r['real'],
                                 'agent': r['agent'], 'when': r['end'],
                                 'outcome': r['state'], 'cause': cause})
    mean_t = (sum(times) / len(times)) if times else None
    phi, phi_long, phi_groups = read_phi(master, case)
    started = batches[0]['first'] if batches else None
    return {'case': case, 'master': str(master), 'state': state,
            'iteration': rec['iteration'], 'noptmax': rec['noptmax'],
            'nreal': rec['nreal'], 'stage': rec['stage'], 'rec': rec,
            'batches': rows, 'agents': sorted(agents,
                                              key=lambda a: a['agent']),
            'n_agents': len({r['agent'] for b in batches
                             for r in b['runs'].values()}),
            'failures': failures, 'mean_run_min': mean_t,
            'started': started,
            'elapsed_h': ((now - started).total_seconds() / 3600.0
                          if started else None),
            'eta': _eta(rec, batches, state, mean_t, now),
            'phi': phi, 'phi_long': phi_long, 'phi_groups': phi_groups,
            'reports': reports(master)}


def _eta(rec, batches, state, mean_t, now):
    """A rough end: what is left of the current batch, then the iterations
    still to come at this run's mean run time. None when it cannot be
    said."""
    if state != 'running' or not batches or not mean_t:
        return None
    b = batches[-1]
    agents = max(len({r['agent'] for r in b['runs'].values()}), 1)
    left = b['planned'] - sum(r['state'] != 'running'
                              for r in b['runs'].values())
    hours = math.ceil(max(left, 0) / agents) * mean_t / 60.0
    nreal = rec['nreal'] or b['planned']
    nlam = max(len(rec['lambda_mults']), 1)
    subset = max(int(math.ceil(0.1 * nreal)), 1)
    per_it = (math.ceil(nlam * subset / agents)
              + math.ceil(max(nreal - subset, 0) / agents)) * mean_t / 60.0
    todo = max((rec['noptmax'] or 0) - max(rec['iteration'], 0), 0)
    if b['stage'] == 'lambda test':
        hours += math.ceil(max(nreal - subset, 0) / agents) * mean_t / 60.0
    return now + dt.timedelta(hours=hours + todo * per_it)
