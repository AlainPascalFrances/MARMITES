# -*- coding: utf-8 -*-
"""The forward run of the La Mata PEST++ calibration (WP7), in a worker.

forward_run.py (written by pest_lamata_setup into the template, so in every
worker) first lets pyemu apply the parameters -- the multiplier files under
``mm_mult/`` (k, sy, ss per cell from the pilot points; drn, sr per layer)
and the one-value parameters ``pars/<name>.csv`` -- then calls
:func:`forward`, which

  1. writes the scalars into the worker's configuration (pest_lamata.toml):
     the soil horizons' field capacity (as its position in [sr, smax]) and
     Ks, the streambed K, the pond-bed leakance, UZF's eps and vks scale,
     CRR's beta -- values --set cannot address (list items);
  2. runs the COUPLED model (code/tests/run_lamata_mf6.py) in the worker:
     its own dataset copy (a run rewrites the forcing files there) and its
     own workspace, MMsurf off, no figures;
  3. reads the observation exports and writes ``pest_obs.csv`` (obsnme,
     value) in the targets' order -- what the instruction file reads;
  4. removes the run's results folder, so a worker does not fill the disk.

A run that fails still writes pest_obs.csv, every value the 1e10 sentinel:
PEST++ then sees a failed realisation, not a missing file.
"""

import glob
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

CONFIG = 'pest_lamata.toml'
# the model's multiplier files (PstFrom reserves 'mult' for its own)
MULT_DIR = 'mm_mult'
# one file per scalar parameter, pars/<name>.csv (name,value): PstFrom
# parameterises a file once in 'direct' style
PARS = 'pars'
TARGETS = 'pest_targets.csv'
OBS_OUT = 'pest_obs.csv'


def apply_scalars(cfg, table, mult_dir=None):
    """The scalar parameters into the configuration. ``table``: DataFrame
    name, value. Names: sfc_f_h<k> and ks_h<k> (horizon k, 1-based, in
    table order), sfr_rhk, lak_bedleak, uzf_eps, uzf_vks, crr_beta; and the
    per-layer drn_l<L> (a multiplier) and sr_l<L> (the specific retention,
    a water content: thts = Sy + Sr), written as the one-value files
    ``mult_dir/<name>.txt`` the model reads (marmites_pest). Ss is on pilot
    points now (2026-10-18), a per-cell file of its own."""
    import marmites_pest as mp
    for name, value in zip(table['name'], table['value']):
        v = float(value)
        if name.startswith(('drn_l', 'sr_l')):
            if mult_dir is None:
                raise ValueError('%s needs the multiplier folder' % name)
            with open(os.path.join(mult_dir, '%s.txt' % name), 'w') as fh:
                fh.write('%.10g\n' % v)
        elif name.startswith('sfc_f_h'):
            h = cfg.soil.horizon[int(name[7:]) - 1]
            f = min(max(v, 0.0), 1.0)
            h.sfc = h.sr + f * (h.smax - h.sr)
        elif name.startswith('ks_h'):
            cfg.soil.horizon[int(name[4:]) - 1].ks = v
        elif name == 'sfr_rhk':
            mp.set_param(cfg, 'sfr.rhk.value', v)
        elif name == 'lak_bedleak':
            mp.set_param(cfg, 'lak.bedleak', v)
        elif name == 'uzf_eps':
            mp.set_param(cfg, 'uzf.eps.value', min(max(v, 3.5), 14.0))
        elif name == 'uzf_vks':
            mp.set_param(cfg, 'uzf.vks_scale', v)
        elif name == 'crr_beta':
            mp.set_param(cfg, 'crr.beta', min(max(v, 1e-3), 1.0))

        else:
            raise ValueError('unknown scalar parameter %r' % name)
    return cfg


def read_scalars(wd):
    """The scalar parameters of ``wd/pars/*.csv`` as one table."""
    import pandas as pd
    files = sorted(glob.glob(os.path.join(wd, PARS, '*.csv')))
    if not files:
        raise FileNotFoundError('no scalar parameter in %s'
                                % os.path.join(wd, PARS))
    t = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    return t.rename(columns={'scalar': 'name'})


def write_scalars(wd, table):
    """One file per row of ``table`` (name, value), pars/<name>.csv. The
    column is called 'scalar' in the file: PstFrom copies the index column
    into the parameter table, where 'name' is an alias of parnme."""
    d = os.path.join(wd, PARS)
    os.makedirs(d, exist_ok=True)
    for name, value in zip(table['name'], table['value']):
        with open(os.path.join(d, '%s.csv' % name), 'w') as fh:
            fh.write('scalar,value\n%s,%.10g\n' % (name, float(value)))


def _sentinel(targets, path):
    import pandas as pd
    pd.DataFrame({'obsnme': targets['obsnme'], 'value': 1.0e10}).to_csv(
        path, index=False, float_format='%.10g')


EXPORTS = os.path.join('ws', 'MF6_ws_voronoi', 'obs_exports')


def run_model(wd, python=None):
    """Steps 1-2: the scalars into the configuration, then the coupled
    model in ``wd``. Returns the exports folder, or None if the run failed."""
    import pandas as pd
    import marmites_config as mcfg
    wd = os.path.abspath(wd)
    cfg = mcfg.load_run_config(os.path.join(wd, CONFIG))
    apply_scalars(cfg, read_scalars(wd), mult_dir=os.path.join(wd, MULT_DIR))
    cfg.validate()
    run_cfg = os.path.join(wd, 'pest_lamata_run.toml')
    cfg.write_toml(run_cfg)
    env = dict(os.environ)
    env['MM_WS_ROOT'] = os.path.join(wd, 'ws')
    env['MM_EXAMPLE_ROOT'] = os.path.join(wd, 'example')
    env['PYTHONUNBUFFERED'] = '1'
    exe = python or sys.executable
    driver = os.path.join(CODE, 'tests', 'run_lamata_mf6.py')
    exports = os.path.join(wd, EXPORTS)
    if os.path.isdir(exports):
        shutil.rmtree(exports, ignore_errors=True)
    t0 = time.time()
    with open(os.path.join(wd, 'model.log'), 'w') as log:
        rc = subprocess.call([exe, '-u', driver, '--config', run_cfg],
                             cwd=wd, env=env, stdout=log,
                             stderr=subprocess.STDOUT)
    print('model run: exit %s, %.1f min' % (rc, (time.time() - t0) / 60.0))
    if rc != 0 or not os.path.exists(os.path.join(exports, 'obs_et.csv')):
        print('the model failed, see %s' % os.path.join(wd, 'model.log'))
        return None
    return exports


def forward(wd='.', python=None, keep=False):
    """One forward run in the worker directory ``wd``. Returns 0 on
    success."""
    import pandas as pd
    import marmites_config as mcfg
    import pest_lamata_obs as po
    wd = os.path.abspath(wd)
    targets = pd.read_csv(os.path.join(wd, TARGETS))
    out = os.path.join(wd, OBS_OUT)
    if os.path.exists(out):
        os.remove(out)
    try:
        exports = run_model(wd, python=python)
        if exports is None:
            _sentinel(targets, out)
            return 1
        # the targets' own settings (q_eps): the [pest] of the copy
        pest = mcfg.load_run_config(os.path.join(wd, CONFIG)).pest
        po.collect(exports, targets, pest=pest).to_csv(
            out, index=False, float_format='%.10g')
        return 0
    except Exception as exc:                          # pragma: no cover
        print('forward run failed: %r' % exc)
        _sentinel(targets, out)
        return 1
    finally:
        if not keep:
            for d in glob.glob(os.path.join(wd, 'ws', 'out_*')):
                shutil.rmtree(d, ignore_errors=True)
