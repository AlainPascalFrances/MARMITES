# -*- coding: utf-8 -*-
"""What a finished run leaves to build on (WP7 front-end, phase 3): the
configuration it RAN -- the file with its overrides -- and the state it
saved where it ended. A PEST++ calibration starts from both: the
configuration is the model it calibrates, the end state is where every
realisation starts (docs/MM-MF6/WP7_PEST_frontend_plan.md).

A run the front-end launches (app/lib/runs.py) is told its registry folder,
``<runs_dir>/<run_id>/``, in the environment (RUN_DIR_ENV), and writes there

    resolved_config.toml   the configuration as the run resolved it
    run_info.json          its results folder, its state workspace and, once
                           saved, its end state

An OLDER run has neither. Both are recovered from what it did leave: its
log names the state it saved ("final heads saved: <prefix>_l1.asc, ...")
and its results folder ``out_<stamp>_<tag>`` holds
``_input/resolved_config.toml``, written by every run since WP0.5.

A saved state carries a NAME, and the next run asked to save under the same
name overwrites it. So an end state is only this run's while its files are
no newer than the run's end -- a later one is refused, never used with a
configuration that did not produce it.
"""

import glob
import json
import os
import re
import time

__all__ = ['RUN_DIR_ENV', 'INFO', 'CONFIG', 'run_dir_from_env', 'runs_dir',
           'record', 'read', 'base_of']

RUN_DIR_ENV = 'MM_RUN_DIR'
INFO = 'run_info.json'
CONFIG = 'resolved_config.toml'
# a state file written this long after the run's last log line is not its own
_SLACK_S = 300.0


def runs_dir(cfg, ws_root):
    """The run registry: ``ui.runs_dir``, blank -> WS_ROOT/runs."""
    return (cfg.ui.runs_dir if cfg.ui.runs_dir
            else os.path.join(str(ws_root), 'runs'))


def run_dir_from_env():
    """The registry folder of the run this process is, or ''."""
    return os.environ.get(RUN_DIR_ENV, '') or ''


def record(run_dir, **items):
    """Merge ``items`` into the run's run_info.json. Atomic: a page reading
    it never sees half a file. Nothing happens without a folder."""
    if not run_dir:
        return None
    os.makedirs(run_dir, exist_ok=True)
    info = read(run_dir)
    info.update(items)
    path = os.path.join(run_dir, INFO)
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(info, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)
    return path


def read(run_dir):
    """The run's run_info.json as a dict ({} when there is none)."""
    try:
        with open(os.path.join(run_dir, INFO), encoding='utf-8') as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def _log_lines(run_dir):
    try:
        with open(os.path.join(run_dir, 'run.log'), encoding='utf-8',
                  errors='replace') as fh:
            return fh.read().splitlines()
    except OSError:
        return []


_SAVED = re.compile(r'^(final|equilibrated|NOT-converged spin-up) heads '
                    r'saved: (.+?)_l1\.asc')
_WRITTEN = re.compile(r'(?:results|file\(s\)) written to (.+?)\s*$')


def _from_log(run_dir, ws_root, status):
    """An older run's configuration copy and end state, from its log and its
    results folder."""
    out = {}
    lines = _log_lines(run_dir)
    for ln in lines:
        m = _SAVED.match(ln.strip())
        if m:
            out['state'] = os.path.basename(m.group(2))
            out['converged'] = m.group(1) != 'NOT-converged spin-up'
    out_dir = ''
    for ln in lines:
        m = _WRITTEN.search(ln)
        if m:
            d = m.group(1).strip()
            d = os.path.dirname(d) if os.path.basename(d) in (
                '_input', '_output') else d
            if os.path.basename(d).startswith('out_'):
                out_dir = d
    if not out_dir and ws_root and status.get('started'):
        # out_<YYYYMMDDHHMM>_<tag>, stamped when the driver started: the
        # registry's minute or the next
        try:
            t0 = time.mktime(time.strptime(status['started'],
                                           '%Y-%m-%d %H:%M:%S'))
        except ValueError:
            t0 = None
        for dt_ in ((0, 60) if t0 is not None else ()):
            stamp = time.strftime('%Y%m%d%H%M', time.localtime(t0 + dt_))
            hits = sorted(glob.glob(os.path.join(str(ws_root),
                                                 'out_%s_*' % stamp)))
            if hits:
                out_dir = hits[0]
                break
    if out_dir:
        out['out_dir'] = out_dir
        cfg = os.path.join(out_dir, '_input', CONFIG)
        if os.path.exists(cfg):
            out['config'] = cfg
    return out


def _status(run_dir):
    try:
        with open(os.path.join(run_dir, 'status.json'), encoding='utf-8') as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def _ended(run_dir, status):
    """When the run ended: its log's last write (what the registry records
    too), else the status' finished stamp."""
    p = os.path.join(run_dir, 'run.log')
    if os.path.exists(p):
        return os.path.getmtime(p)
    try:
        return time.mktime(time.strptime(status.get('finished') or '',
                                         '%Y-%m-%d %H:%M:%S'))
    except ValueError:
        return None


def base_of(run_dir, ws_root=None):
    """What a calibration can start from in this run, as a dict:

        config     the configuration copy (a path), or None
        state      the end state's name, or None
        state_dir  the workspace it is in
        out_dir    the run's results folder, when known
        converged  False for a spin-up saved before it converged
        source     'run_info' (a run that recorded it) or 'log'
        why        '' when usable, else the reason in a sentence

    ``ws_root``: WS_ROOT, where an older run's results folder and, by
    default, its state workspace are.
    """
    import marmites_config as mcfg
    status = _status(run_dir)
    info = read(run_dir)
    if info:
        res = {'source': 'run_info', 'out_dir': info.get('out_dir'),
               'state': info.get('state'), 'state_dir': info.get('state_dir'),
               'converged': info.get('converged', True)}
        c = os.path.join(run_dir, info.get('config') or CONFIG)
        res['config'] = c if os.path.exists(c) else None
    else:
        res = dict({'source': 'log', 'state': None, 'config': None,
                    'converged': True, 'state_dir': None},
                   **_from_log(run_dir, ws_root, status))
    res['why'] = ''
    if status.get('state') == 'running':
        res['why'] = 'it is still running'
        return res
    if status.get('outcome') not in (None, 'completed'):
        res['why'] = 'it did not complete (%s)' % status.get('outcome')
        return res
    if not res.get('config'):
        res['why'] = ('no copy of its configuration was found (neither in its '
                      'registry folder nor in its results folder)')
        return res
    if not res.get('state'):
        res['why'] = ('it saved no end state (spinup.save_strt was blank and '
                      'it was a single run)')
        return res
    try:
        cfg = mcfg.load_run_config(res['config'])
    except Exception as exc:                            # noqa: BLE001
        res['why'] = 'its configuration copy no longer loads: %s' % exc
        return res
    res['nlay'] = int(cfg.layers.nlay)
    if not res.get('state_dir'):
        res['state_dir'] = (cfg.paths.ws or
                            mcfg.state_workspace(cfg, ws_root or ''))
    files = ['%s_l%d.asc' % (res['state'], k + 1)
             for k in range(max(res['nlay'], 1))]
    paths = [os.path.join(str(res['state_dir']), f) for f in files]
    missing = [f for f, p in zip(files, paths) if not os.path.exists(p)]
    if missing:
        res['why'] = ('its end state %s is not in %s any more (%s missing)'
                      % (res['state'], res['state_dir'], ', '.join(missing)))
        return res
    newest = max(os.path.getmtime(p) for p in paths)
    res['state_time'] = newest
    end = _ended(run_dir, status)
    if end is not None and newest > end + _SLACK_S:
        res['why'] = ('its end state %s was saved again by a later run (%s): '
                      'it is not this run\'s state any more'
                      % (res['state'], time.strftime(
                          '%d %b %Y %H:%M', time.localtime(newest))))
    return res
