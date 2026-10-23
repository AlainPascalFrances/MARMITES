# -*- coding: utf-8 -*-
"""The calibrated model as a configuration of its own (WP7 front-end,
phase 5): <model>_PEST.toml beside the input configuration (lamata.toml ->
lamata_PEST.toml), whose parameters are one realisation of a finished
calibration -- by default the 'base' realisation of the last iteration,
the one PEST++ carried from the starting values (user, 2026-10-23).

    python code/tools/pest_lamata_export.py [--pest-ws WS] [--iter N]
                                            [--real base|best|NAME]
                                            [--out PATH] [--force]

1. The realisation's values go through the calibration's own templates in
   a copy of its template folder (<pest ws>/results/it<N>_<real>/): the
   pilot points kriged into the per-cell multipliers (mm_mult/), the
   scalars written to pars/*.csv -- what a worker gets before its run.
2. The configuration is the BASE RUN's (setup.json: model_config), with
   what a calibration run differs by: one spin-up cycle from the base
   state, a cold start's UZF empty, [pest] multipliers -> that mm_mult
   folder (absolute), and the scalars applied (pest_lamata_forward). The
   end state is saved as <state>_PEST, never over the base state.
3. A pilot-point property is a multiplier of its reference raster: one
   changed since the set-up (layer 1's Sy, 0.01 -> 0.025 on 2026-10-22)
   would multiply into another model. Its multipliers are rescaled when
   both rasters are uniform; otherwise the export refuses.

Nothing is run. Provenance in <pest ws>/results/it<N>_<real>/export.json.
"""

import argparse
import datetime as dt
import filecmp
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np                                              # noqa: E402

# the template folder without the model's own copies (dataset, state, code)
SKIP = shutil.ignore_patterns('ws', 'example', 'code', 'base_exports',
                              'model.log', 'pest_obs.csv', '__pycache__')
# the pilot-point families: mm_mult/<stem>_l<L>.txt multiplies layers.<field>
PILOT = (('k', 'k'), ('sy', 'sy'), ('ss', 'ss'))


def _phi_row(master, it):
    import pandas as pd
    phi = pd.read_csv(os.path.join(master, 'lamata.phi.actual.csv'))
    row = phi[phi['iteration'] == it]
    if row.empty:
        raise SystemExit('iteration %d is not in %s' % (it, master))
    return row.iloc[0].drop(['iteration', 'total_runs', 'mean',
                             'standard_deviation', 'min', 'max'],
                            errors='ignore').dropna().astype(float)


def realisation(master, it, real):
    """(pst, parameter values, name, phi) of realisation ``real`` of
    iteration ``it``: 'base', 'best' (lowest phi) or a name."""
    import pyemu
    pst = pyemu.Pst(os.path.join(master, 'lamata.pst'))
    pe = pyemu.ParameterEnsemble.from_binary(
        pst, os.path.join(master, 'lamata.%d.par.jcb' % it))._df
    pe.index = pe.index.astype(str)
    phi = _phi_row(master, it)
    phi.index = phi.index.astype(str)
    if real == 'best':
        real = phi[phi.index.isin(pe.index)].idxmin()
    for key in (real, real.lower(), real.upper()):
        if key in pe.index:
            return pst, pe.loc[key], key, float(phi.get(key, np.nan))
    raise SystemExit('no realisation %r in iteration %d (%d there)'
                     % (real, it, len(pe.index)))


def write_inputs(pst, values, work):
    """The values into the template files in ``work``, then the pilot
    points kriged into mm_mult/ -- forward_run.py's first step."""
    import pyemu
    par = pst.parameter_data
    par.loc[values.index, 'parval1'] = values.values.astype(float)
    pst.write_input_files(pst_path=work)
    cwd = os.getcwd()
    os.chdir(work)
    try:
        pyemu.helpers.apply_list_and_array_pars(
            arr_par_file='mult2model_info.csv', chunk_len=50)
    finally:
        os.chdir(cwd)


def _uniform(path):
    a = np.loadtxt(path, skiprows=6)
    v = np.unique(a[(a > -9999.0) & (a != 0.0)])
    return float(v[0]) if v.size == 1 else None


def rescale(cfg, mult, tpl_ds, user_ds):
    """Rescale the multipliers of a reference raster that changed since
    the set-up; a list of what was done. SystemExit when it cannot be."""
    notes = []
    nlay = int(cfg.layers.nlay)
    for stem, field in PILOT:
        src = getattr(cfg.layers, field)
        if not getattr(src, 'raster', ''):
            continue
        for L in range(nlay):
            f = os.path.join(mult, '%s_l%d.txt' % (stem, L + 1))
            if not os.path.exists(f):
                continue
            name = src.raster % (L + 1) if '%d' in src.raster else src.raster
            old, new = os.path.join(tpl_ds, name), os.path.join(user_ds, name)
            if filecmp.cmp(old, new, shallow=False):
                continue
            a, b = _uniform(old), _uniform(new)
            if a is None or b is None:
                raise SystemExit(
                    '%s changed since the set-up and is not uniform: its '
                    'multipliers (%s) cannot be carried over. Put the set-up\'s '
                    'raster back (%s) or set up the calibration again.'
                    % (new, f, old))
            m = np.loadtxt(f)
            np.savetxt(f, m * (a / b), fmt='%14.6E')
            notes.append('%s: %s was %g at the set-up and is %g now -> '
                         'multipliers x%g' % (os.path.basename(f), name, a, b,
                                              a / b))
    return notes


def export(pest_ws, it=None, real='base', out=None, force=False):
    import marmites_config as mcfg
    import mm_paths
    import pest_lamata_forward as plf
    master = os.path.join(pest_ws, 'master')
    tpl = os.path.join(pest_ws, 'template')
    with open(os.path.join(tpl, 'setup.json')) as fh:
        setup = json.load(fh)
    if it is None:
        it = int(_phi_row_last(master))
    pst, values, name, phi = realisation(master, it, real)
    tag = 'it%d_%s' % (it, name)
    work = os.path.join(pest_ws, 'results', tag)
    if os.path.isdir(work):
        if not force:
            raise SystemExit('%s exists (--force to make it again)' % work)
        shutil.rmtree(work)
    shutil.copytree(tpl, work, ignore=SKIP)
    import pyemu
    pst = pyemu.Pst(os.path.join(work, 'lamata.pst'))
    write_inputs(pst, values, work)
    # the configuration: the base run's, as a calibration run has it
    cfg = mcfg.load_run_config(setup['model_config'])
    state = setup['state']
    cfg.spinup.cycles = 1
    cfg.spinup.strt_heads = state
    cfg.spinup.save_strt = state + '_PEST'
    cfg.spinup.save_means = ''
    cfg.uzf.thti_from = 'thtr'
    cfg.meta.name = (cfg.meta.name or 'model') + '_PEST'
    mult = os.path.join(work, 'mm_mult')
    plf.apply_scalars(cfg, plf.read_scalars(work), mult_dir=mult)
    cfg.pest.multipliers = mult
    case = cfg.paths.case
    notes = rescale(cfg, mult, os.path.join(tpl, 'example', case),
                    str(mm_paths.dataset_dir(case)))
    cfg.validate()
    if out is None:
        stem = os.path.splitext(os.path.basename(setup['config']))[0]
        out = os.path.join(str(mm_paths.CONFIG_DIR), stem + '_PEST.toml')
    if os.path.exists(out) and not force:
        raise SystemExit('%s exists (--force to replace it)' % out)
    cfg.write_toml(out)
    info = {'made': dt.datetime.now().isoformat(timespec='seconds'),
            'master': master, 'iteration': it, 'realisation': name,
            'phi': phi, 'base': setup.get('base'), 'state': state,
            'model_config': setup['model_config'], 'config': out,
            'multipliers': mult, 'rescaled': notes}
    with open(os.path.join(work, 'export.json'), 'w') as fh:
        json.dump(info, fh, indent=1)
    return info


def _phi_row_last(master):
    import pandas as pd
    return pd.read_csv(os.path.join(master,
                                    'lamata.phi.actual.csv'))['iteration'].max()


def main(argv=None):
    import mm_paths
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--pest-ws', default=getattr(mm_paths, 'PEST_WS', None))
    ap.add_argument('--iter', type=int, default=None,
                    help='the iteration (default: the last)')
    ap.add_argument('--real', default='base',
                    help="base (default), best, or a realisation's name")
    ap.add_argument('--out', default=None)
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args(argv)
    if not a.pest_ws:
        raise SystemExit('no PEST workspace (--pest-ws, or set it on Home)')
    info = export(a.pest_ws, a.iter, a.real, a.out, a.force)
    print('realisation %s of iteration %d (phi %.1f) -> %s'
          % (info['realisation'], info['iteration'], info['phi'],
             info['config']))
    print('multipliers: %s' % info['multipliers'])
    for n in info['rescaled']:
        print('rescaled: %s' % n)


if __name__ == '__main__':
    main()
