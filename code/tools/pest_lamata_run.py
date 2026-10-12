# -*- coding: utf-8 -*-
"""Run pestpp-ies on the La Mata template (WP7): the master and its local
PANTHER agents, in the PEST workspace (Home, pest_ws).

    python code/tools/pest_lamata_run.py --workers 10 --noptmax 0
    python code/tools/pest_lamata_run.py --join 4      # agents for a
                                                        # running master

``--noptmax 0`` is the base run through PEST++ itself (one forward run in
an agent: the check that the whole chain works), ``-1`` the prior Monte
Carlo, ``2``/``3`` the IES iterations. The template is copied into
``<pest_ws>/workers/worker_<k>`` (each worker its own dataset copy and
workspace) and the master runs in ``<pest_ws>/master``.

``--join N`` adds N agents to a master that is ALREADY running -- PANTHER
takes agents at any time -- in new worker folders after the existing ones,
never touching those (pyemu.start_workers would recreate them).

At most MAX_WORKERS agents in all: the user allows 10 of the server's 16
logical cores (2026-10-12; was 6, the physical cores minus two).
"""

import argparse
import os
import shutil
import socket
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
if CODE not in sys.path:
    sys.path.insert(0, CODE)

MAX_WORKERS = 10


def _running_agents(workers):
    """Worker folders whose agent still holds the folder (a pestpp-ies
    process started there): the ones --join must leave alone."""
    names = sorted(d for d in os.listdir(workers)
                   if d.startswith('worker_')) if os.path.isdir(workers) \
        else []
    return names


def join(pest_ws, n, port):
    """Start ``n`` more agents for the master listening on ``port``, in
    worker folders numbered after the existing ones. Returns the Popen
    handles (the agents end with the master)."""
    import mm_paths
    tpl = os.path.join(pest_ws, 'template')
    workers = os.path.join(pest_ws, 'workers')
    have = _running_agents(workers)
    if len(have) + n > MAX_WORKERS:
        raise SystemExit('%d agents exist, %d more would pass the limit of '
                         '%d' % (len(have), n, MAX_WORKERS))
    first = 1 + max([int(d.split('_')[1]) for d in have] or [-1])
    exe = str(mm_paths.PESTPP_IES)
    host = socket.gethostname()
    procs = []
    for k in range(first, first + n):
        d = os.path.join(workers, 'worker_%d' % k)
        if os.path.exists(d):
            raise SystemExit('%s exists: not overwriting a worker' % d)
        shutil.copytree(tpl, d)
        procs.append(subprocess.Popen(
            [exe, 'lamata.pst', '/h', '%s:%d' % (host, port)], cwd=d,
            stdout=open(os.path.join(d, 'agent_stdout.txt'), 'w'),
            stderr=subprocess.STDOUT))
        print('agent %d started in %s' % (k, d))
        time.sleep(2)
    return procs


def main(argv=None):
    import pyemu
    import mm_paths
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--workers', type=int, default=MAX_WORKERS)
    ap.add_argument('--noptmax', type=int, default=None)
    ap.add_argument('--nreal', type=int, default=None,
                    help='ies_num_reals (default: as set up)')
    ap.add_argument('--port', type=int, default=4004)
    ap.add_argument('--join', type=int, default=0,
                    help='add this many agents to the RUNNING master')
    a = ap.parse_args(argv)
    pest_ws = str(mm_paths.PEST_WS or '')
    if a.join:
        procs = join(pest_ws, a.join, a.port)
        for p in procs:
            p.wait()
        return
    if a.noptmax is None:
        raise SystemExit('--noptmax is required to start a run')
    if a.workers > MAX_WORKERS:
        raise SystemExit('at most %d workers (the user allows %d cores)'
                         % (MAX_WORKERS, MAX_WORKERS))
    tpl = os.path.join(pest_ws, 'template')
    if not os.path.exists(os.path.join(tpl, 'lamata.pst')):
        raise SystemExit('no template in %s: run pest_lamata_setup.py' % tpl)
    pst = pyemu.Pst(os.path.join(tpl, 'lamata.pst'))
    pst.control_data.noptmax = a.noptmax
    if a.nreal:
        pst.pestpp_options['ies_num_reals'] = a.nreal
    pst.write(os.path.join(tpl, 'lamata.pst'), version=2)
    master = os.path.join(pest_ws, 'master')
    if os.path.isdir(master):
        shutil.rmtree(master)
    workers = os.path.join(pest_ws, 'workers')
    os.makedirs(workers, exist_ok=True)
    exe = str(mm_paths.PESTPP_IES)
    print('pestpp-ies %s, noptmax %d, %d worker(s), master %s'
          % (exe, a.noptmax, a.workers, master))
    pyemu.os_utils.start_workers(
        tpl, exe, 'lamata.pst', num_workers=a.workers,
        worker_root=workers, master_dir=master,
        port=a.port, cleanup=False)


if __name__ == '__main__':
    main()
