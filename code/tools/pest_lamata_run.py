# -*- coding: utf-8 -*-
"""Run pestpp-ies on the La Mata template (WP7): the master and its local
PANTHER agents, in the PEST workspace (Home, pest_ws).

    python code/tools/pest_lamata_run.py --config code/configs/lamata.toml
    python code/tools/pest_lamata_run.py --join 4      # agents for a
                                                        # running master

A launch takes its settings from the configuration's [pest] section (the
PEST++ set-up page), read at launch, so they can change after the set-up:
``agents``, ``noptmax`` (``-1`` the prior only, ``0`` one run at the
starting values -- the check that the whole chain works --, ``n`` the IES
iterations), the lambda multipliers, the bad-phi sigma and the overdue
limit, and the START:

  prior    the prior ensemble drawn at set-up (prior_pe.jcb)
  warm     a parameter ensemble (.jcb, ``start_from``) after a MODEL change:
           the realisations are run first with the current model, then
           updated (the warm start of 2026-10-16, after the Tg curve was
           fixed)
  restart  an interrupted run with the SAME model: ``start_from`` is its
           master folder, whose completed runs' outputs are kept

The template is copied into ``<pest_ws>/workers/worker_<k>`` (each worker
its own dataset copy and workspace) and the master runs in
``<pest_ws>/master``.

``--join N`` adds N agents to a master that is ALREADY running -- PANTHER
takes agents at any time -- in new worker folders after the existing ones,
never touching those (pyemu.start_workers would recreate them).

At most MAX_WORKERS agents in all: the user allows 20 (2026-10-16; the
server has 32 cores / 64 threads since 2026-10-15; was 10 of 16 threads).
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

MAX_WORKERS = 20


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


def pestpp_options(p):
    """The pestpp-ies options of a configuration's [pest] (``p``,
    marmites_config.Pest) that a launch may change after the set-up."""
    lm = [float(x) for x in p.lambda_mults.replace(',', ' ').split()]
    return {'ies_lambda_mults': ','.join(repr(x) for x in lm),
            'ies_bad_phi_sigma': float(p.bad_phi_sigma),
            # a run is given up at this x the running mean: at 2x the
            # production lost 12 of 43 runs of iteration 1 (2026-10-16) --
            # early in a batch the mean is the first, fastest runs'
            'overdue_giveup_fac': float(p.overdue_giveup_fac)}


def prior_files(tpl, pst, nreal):
    """Start from the prior drawn at set-up: ``nreal`` of its
    realisations, no outputs of an earlier run. Returns the pestpp
    options."""
    import pyemu
    n = len(pyemu.ParameterEnsemble.from_binary(
        pst, os.path.join(tpl, 'prior_pe.jcb'))._df)
    if nreal > n:
        raise SystemExit('the prior has %d realisation(s), pest.nreal asks '
                         'for %d: set up again to draw more' % (n, nreal))
    for k in ('ies_restart_observation_ensemble', 'ies_observation_ensemble'):
        pst.pestpp_options.pop(k, None)
    return {'ies_parameter_ensemble': 'prior_pe.jcb', 'ies_num_reals': nreal}


def restart_files(src_master, tpl):
    """Restart an interrupted run from what its master completed: the
    parameters and outputs of the COMPLETED runs of the prior, read back
    from the run storage (lamata.rns), and their noise realisations --
    written into the template for pestpp-ies to start at iteration 1
    (2026-10-15: a server reboot stopped the production prior at 88 of
    100). Returns the number of realisations kept and the pestpp options."""
    import pyemu
    pst = pyemu.Pst(os.path.join(src_master, 'lamata.pst'))
    par, obs, meta = pyemu.helpers.read_pestpp_runstorage(
        os.path.join(src_master, 'lamata.rns'), irun='all',
        with_metadata=True)
    done = meta['status'] == 'completed'
    names = [t.split('realization:')[-1].strip()
             for t in meta.loc[done, 'info_txt']]
    if not names:
        raise SystemExit('no completed run in %s' % src_master)
    p = par.loc[done, pst.par_names].copy()
    o = obs.loc[done, pst.obs_names].copy()
    p.index = o.index = names
    noise = pyemu.ObservationEnsemble.from_binary(
        pst, os.path.join(src_master, 'lamata.obs+noise.jcb'))._df
    noise.index = noise.index.astype(str)
    keep = [n for n in names if n in noise.index]
    pyemu.ParameterEnsemble(pst, p.loc[keep]).to_binary(
        os.path.join(tpl, 'restart_par.jcb'))
    pyemu.ObservationEnsemble(pst, o.loc[keep]).to_binary(
        os.path.join(tpl, 'restart_obs.jcb'))
    pyemu.ObservationEnsemble(pst, noise.loc[keep]).to_binary(
        os.path.join(tpl, 'restart_noise.jcb'))
    opts = {'ies_parameter_ensemble': 'restart_par.jcb',
            'ies_restart_observation_ensemble': 'restart_obs.jcb',
            'ies_observation_ensemble': 'restart_noise.jcb',
            'ies_num_reals': len(keep)}
    print('restart from %s: %d completed realisation(s) of %d (%s)'
          % (src_master, len(keep), len(meta),
             meta['status'].value_counts().to_dict()))
    return len(keep), opts


def start_files(ens, tpl, pst):
    """A warm start: the parameter ensemble ``ens`` (.jcb) as the starting
    one, copied into the template; no observation ensemble -- pestpp-ies
    runs the realisations with the current model first. Returns the number
    of realisations and the pestpp options."""
    import numpy as np
    import pyemu
    pe = pyemu.ParameterEnsemble.from_binary(pst, ens)._df
    missing = sorted(set(pst.par_names) - set(pe.columns))
    extra = sorted(set(pe.columns) - set(pst.par_names))
    if missing or extra:
        raise SystemExit('%s does not carry the template\'s parameters: '
                         '%d missing (%s), %d extra (%s)'
                         % (ens, len(missing), missing[:3], len(extra),
                            extra[:3]))
    pe = pe[pst.par_names]
    pd_ = pst.parameter_data.loc[pst.par_names]
    lo, hi = pd_['parlbnd'].values, pd_['parubnd'].values
    # a value ON a bound comes back a rounding off it (25 of iteration 1's,
    # 2026-10-16): within 1e-6 it is set on the bound, beyond it refused
    tol = 1e-6 * np.maximum(np.abs(np.c_[lo, hi]), 1.0)
    out = (pe < lo - tol[:, 0]) | (pe > hi + tol[:, 1])
    if out.values.any():
        raise SystemExit('%s: %d value(s) outside the template\'s bounds'
                         % (ens, int(out.values.sum())))
    pe = pe.clip(lower=lo, upper=hi, axis=1)
    pyemu.ParameterEnsemble(pst, pe).to_binary(
        os.path.join(tpl, 'start_par.jcb'))
    print('warm start from %s: %d realisation(s)' % (ens, len(pe)))
    return len(pe), {'ies_parameter_ensemble': 'start_par.jcb',
                     'ies_num_reals': len(pe)}


def main(argv=None):
    import pyemu
    import marmites_config as mcfg
    import mm_paths
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--config', default='',
                    help='the configuration whose [pest] section sets the '
                         'launch (the PEST++ set-up page)')
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
    if not a.config:
        raise SystemExit('--config is required to start a run')
    cfg = mcfg.load_run_config(os.path.abspath(a.config))
    p = cfg.pest
    if p.agents > MAX_WORKERS:
        raise SystemExit('at most %d agents (the user allows %d cores), '
                         'pest.agents is %d' % (MAX_WORKERS, MAX_WORKERS,
                                                p.agents))
    tpl = os.path.join(pest_ws, 'template')
    if not os.path.exists(os.path.join(tpl, 'lamata.pst')):
        raise SystemExit('no template in %s: run pest_lamata_setup.py' % tpl)
    pst = pyemu.Pst(os.path.join(tpl, 'lamata.pst'))
    pst.control_data.noptmax = int(p.noptmax)
    pst.pestpp_options.update(pestpp_options(p))
    master = os.path.join(pest_ws, 'master')
    if p.start == 'prior':
        pst.pestpp_options.update(prior_files(tpl, pst, int(p.nreal)))
    elif p.start == 'warm':
        # the outputs and noise of an earlier run belong to the old model
        for k in ('ies_restart_observation_ensemble',
                  'ies_observation_ensemble'):
            pst.pestpp_options.pop(k, None)
        _n, opts = start_files(os.path.abspath(p.start_from), tpl, pst)
        pst.pestpp_options.update(opts)
    else:                                                    # restart
        src = os.path.abspath(p.start_from)
        if os.path.normcase(src) == os.path.normcase(master):
            raise SystemExit('move the interrupted master aside first: '
                             'this run recreates %s' % src)
        _n, opts = restart_files(src, tpl)
        pst.pestpp_options.update(opts)
    pst.write(os.path.join(tpl, 'lamata.pst'), version=2)
    if os.path.isdir(master):
        # an earlier calibration's master is KEPT, moved aside under its
        # last date (2026-10-20: it was deleted -- a launch from the PEST++
        # set-up page must never lose a calibration's results)
        aside = '%s_%s' % (master, time.strftime(
            '%Y%m%d_%H%M%S', time.localtime(os.path.getmtime(master))))
        os.rename(master, aside)
        print('the previous master moved to %s' % aside)
    workers = os.path.join(pest_ws, 'workers')
    os.makedirs(workers, exist_ok=True)
    exe = str(mm_paths.PESTPP_IES)
    print('pestpp-ies %s, noptmax %d, start %s, %d agent(s), master %s'
          % (exe, p.noptmax, p.start, p.agents, master))
    pyemu.os_utils.start_workers(
        tpl, exe, 'lamata.pst', num_workers=int(p.agents),
        worker_root=workers, master_dir=master,
        port=a.port, cleanup=False)


if __name__ == '__main__':
    main()
