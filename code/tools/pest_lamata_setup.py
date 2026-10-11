# -*- coding: utf-8 -*-
"""Set up the La Mata PEST++-IES calibration (WP7) in the PEST workspace.

    python code/tools/pest_lamata_setup.py --config code/configs/lamata.toml \\
        --state hi_voronoi_lamata_evt_drn1 [--nreal 30] [--step all]

The PEST workspace is the one set on Home (``pest_ws``; mm_paths.PEST_WS),
never WS_ROOT and never the repository. Steps (``--step``):

  org     the model as every worker will run it, in ``<pest_ws>/org``:
            code/         a SNAPSHOT of the repository's code, so a
                          calibration of days is not changed under it
            pest_lamata.toml  the configuration copy: MMsurf off, no
                          figures, one cycle from the saved state, nothing
                          saved, [pest] multipliers = mult
            example/LaMata the dataset copy (a run rewrites the forcing
                          files there), the observation table with C4, C5
                          and W1 switched on
            ws/MF6_ws_voronoi the mesh and overlay caches and the saved
                          state every realisation starts from
            mm_mult/         k_l<L>, sy_l<L> per cell (the pilot points)
            pars/scalars.csv the scalar parameters (soil horizons, SFR, LAK,
                          UZF, CRR, Ss and DRN per layer)
            modis_et_catchment.csv, ec_et.csv  the ET records
  base    the base run in org (the test run of the setup, ~45 min), then
          the targets (pest_targets.csv) and their weights, balanced so each
          group carries its share of the initial phi
  pst     PstFrom -> ``<pest_ws>/template``: pilot points for K and Sy per
          layer on the mesh, the scalars, the observations, the prior
          ensemble drawn from the geostatistical structure, pestpp-ies
          options; noptmax = 0 until the pilot is launched

Model runs happen only here, in the PEST workspace (CLAUDE.md rule 1, the
PEST++ exception of 2026-10-11).
"""

import argparse
import datetime as dt
import json
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

CALIB_START = dt.date(2009, 10, 1)      # after the warm-up hydrological year
PP_SPACE = 400.0                        # m between pilot points
VARIO_A = 1500.0                        # m, exponential variogram range
from pest_lamata_forward import MULT_DIR  # noqa: E402
EC_XLSX = os.path.join('MMsurf_ws', 'Eddy_MM_202110.xlsx')   # in DATA_ROOT


# ------------------------------------------------------------------ org

def _check_ws(pest_ws, mm_paths):
    if not pest_ws:
        raise SystemExit('the PEST workspace is not set: set it on Home '
                         '(pest_ws)')
    real = os.path.normcase(os.path.abspath(str(pest_ws)))
    for name, root in (('WS_ROOT', mm_paths.WS_ROOT),
                       ('the repository', mm_paths.REPO)):
        r = os.path.normcase(os.path.abspath(str(root)))
        if real == r or real.startswith(r + os.sep):
            raise SystemExit('the PEST workspace %s is inside %s' % (pest_ws,
                                                                    name))


def _snapshot_code(dst):
    """The repository's code, without the app, caches or machine files --
    and with the machine paths file, which the snapshot's mm_paths reads."""
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    shutil.copytree(CODE, dst, ignore=shutil.ignore_patterns(
        '__pycache__', '*.pyc', 'app', '.pytest_cache', 'legacy',
        '*.nbi', '*.nbc'))
    cfgd = os.path.join(dst, 'configs')
    for f in os.listdir(cfgd):
        if f.endswith('.toml') and f != 'paths.local.toml':
            os.remove(os.path.join(cfgd, f))


def _obs_table_with(path, names=('C4', 'C5', 'W1')):
    """Switch the commented-out points ``names`` on in an inputObs.txt."""
    lines = open(path, encoding='utf-8', errors='replace').read().splitlines()
    out, on = [], []
    for ln in lines:
        s = ln.lstrip('#').strip()
        tok = s.split()
        if ln.startswith('#') and tok and tok[0] in names and len(tok) >= 5:
            out.append(s)
            on.append(tok[0])
        else:
            out.append(ln)
    with open(path, 'w', encoding='utf-8', newline='\r\n') as fh:
        fh.write('\n'.join(out) + '\n')
    return on


def scalar_table(cfg):
    """The scalar parameters and their starting values, as the forward run
    reads them (pest_lamata_forward.apply_scalars)."""
    import pandas as pd
    rows = []
    for k, h in enumerate(cfg.soil.horizon):
        rows.append(('sfc_f_h%d' % (k + 1),
                     (h.sfc - h.sr) / (h.smax - h.sr)))
    for k, h in enumerate(cfg.soil.horizon):
        rows.append(('ks_h%d' % (k + 1), h.ks))
    rows += [('sfr_rhk', cfg.sfr.rhk.value), ('lak_bedleak', cfg.lak.bedleak),
             ('uzf_eps', cfg.uzf.eps.value), ('uzf_vks', cfg.uzf.vks_scale),
             ('crr_beta', cfg.crr.beta)]
    for L in range(cfg.layers.nlay):
        rows.append(('ss_l%d' % (L + 1), 1.0))
    for L in range(cfg.layers.nlay):
        rows.append(('drn_l%d' % (L + 1), 1.0))
    return pd.DataFrame(rows, columns=['name', 'value'])


def scalar_bounds(name, value):
    """(transform, lower, upper) of a scalar parameter -- the draft's
    table, decisions of 2026-10-11."""
    if name.startswith('sfc_f_h'):
        return 'none', max(value - 0.10, 0.01), min(value + 0.10, 0.99)
    if name.startswith('ks_h'):
        return 'log', value * 0.1, value * 10.0
    return {'sfr_rhk': ('log', 0.001, 1.0),
            'lak_bedleak': ('log', 1.0e-4, 1.0e-2),
            'uzf_eps': ('none', 3.5, 8.0),
            'uzf_vks': ('log', 0.2, 5.0),
            'crr_beta': ('none', 0.8, 1.0)}.get(
        name, ('log', 0.1, 100.0) if name.startswith('drn_l')
        else ('log', 0.1, 10.0))                       # ss_l<L>


def step_org(cfg_path, state, pest_ws):
    import pandas as pd
    import marmites_config as mcfg
    import mm_paths
    import pest_modis_et as pme
    org = os.path.join(pest_ws, 'org')
    os.makedirs(org, exist_ok=True)
    print('PEST workspace %s' % pest_ws)
    _snapshot_code(os.path.join(pest_ws, 'code'))
    print('   code snapshot: %s' % os.path.join(pest_ws, 'code'))
    cfg = mcfg.load_run_config(cfg_path)
    case = cfg.paths.case
    cfg.meta.name = 'pest'
    cfg.run.surface = False
    cfg.run.model = True
    cfg.run.plot = False
    cfg.postproc.input_maps = False
    cfg.spinup.cycles = 1
    cfg.spinup.strt_heads = state
    cfg.spinup.save_strt = ''
    cfg.spinup.save_means = ''
    cfg.paths.ws = ''
    cfg.pest.multipliers = MULT_DIR
    cfg.validate()
    cfg.write_toml(os.path.join(org, 'pest_lamata.toml'))
    # the dataset copy, with the points the calibration adds
    ds = os.path.join(org, 'example', case)
    if os.path.isdir(ds):
        shutil.rmtree(ds)
    shutil.copytree(str(mm_paths.dataset_dir(case)), ds)
    on = _obs_table_with(os.path.join(ds, cfg.obs.table))
    print('   dataset copy, observation points switched on: %s'
          % (', '.join(on) or 'none'))
    # the caches and the saved state
    src = os.path.join(str(mm_paths.WS_ROOT), 'MF6_ws_voronoi')
    dst = os.path.join(org, 'ws', 'MF6_ws_voronoi')
    os.makedirs(dst, exist_ok=True)
    for d in ('_mesh', '_overlay'):
        if os.path.isdir(os.path.join(src, d)):
            shutil.copytree(os.path.join(src, d), os.path.join(dst, d),
                            dirs_exist_ok=True)
    n = 0
    for f in os.listdir(src):
        if f.startswith(state + '_') or f == state + '.scope.json':
            shutil.copy2(os.path.join(src, f), os.path.join(dst, f))
            n += 1
    if not n:
        raise SystemExit('no saved state %s in %s' % (state, src))
    print('   saved state %s: %d file(s)' % (state, n))
    # the parameters' files
    md = os.path.join(org, MULT_DIR)
    os.makedirs(md, exist_ok=True)
    import pest_lamata_forward as pfw
    pd_ = os.path.join(org, pfw.PARS)
    if os.path.isdir(pd_):
        shutil.rmtree(pd_)
    pfw.write_scalars(org, scalar_table(cfg))
    # the ET records
    poly = pme.catchment_polygon(os.path.join(ds, 'inputBOUNDARY.geojson'))
    dates = [l.split(',')[0].strip()[:10]
             for l in open(os.path.join(ds, 'inputDATE.txt'))
             if l.strip() and not l.startswith('#')]
    modis_dir = os.path.join(str(mm_paths.RS), 'MODIS4%s' % case)
    pme.catchment_series(modis_dir, poly, cfg.grid.crs_epsg,
                         dt.date.fromisoformat(dates[0]),
                         dt.date.fromisoformat(dates[-1])).to_csv(
        os.path.join(org, 'modis_et_catchment.csv'), index=False,
        float_format='%.6g')
    ec = pd.read_excel(os.path.join(str(mm_paths.DATA_ROOT), EC_XLSX),
                       sheet_name='compil4python')
    ec = pd.DataFrame({'date': pd.to_datetime(ec['Date']).dt.date,
                       'et_mmd': pd.to_numeric(ec['ET_EC'], errors='coerce')})
    ec.dropna().to_csv(os.path.join(org, 'ec_et.csv'), index=False)
    json.dump({'config': os.path.abspath(cfg_path), 'state': state,
               'made': dt.datetime.now().isoformat(timespec='seconds')},
              open(os.path.join(org, 'setup.json'), 'w'), indent=1)
    return org


GRB = 'lamata.disv.grb'


def _keep_only_inputs(org, state):
    """After the base run: the grid file copied to org, then org's model
    workspace back to what a run needs -- the mesh and overlay caches and the
    saved state -- so the template, copied into every worker, does not carry
    the base run's GB of MF6 output."""
    ws = os.path.join(org, 'ws', 'MF6_ws_voronoi')
    shutil.copy2(os.path.join(ws, GRB), os.path.join(org, GRB))
    for f in os.listdir(ws):
        if f in ('_mesh', '_overlay') or f.startswith(state + '_') \
                or f == state + '.scope.json':
            continue
        p = os.path.join(ws, f)
        if os.path.isdir(p):
            shutil.rmtree(p, ignore_errors=True)
        else:
            os.remove(p)
    for f in os.listdir(os.path.join(org, 'ws')):
        if f.startswith('out_'):
            shutil.rmtree(os.path.join(org, 'ws', f), ignore_errors=True)


def _write_unit_mults(org, nlay):
    """k_l<L>, sy_l<L>: one 1.0 per mesh cell (the grid file gives ncpl)."""
    import flopy
    ncpl = flopy.mf6.utils.MfGrdFile(os.path.join(org, GRB)).modelgrid.ncpl
    for stem in ('k', 'sy'):
        for L in range(nlay):
            np.savetxt(os.path.join(org, MULT_DIR, '%s_l%d.txt' % (stem, L + 1)),
                       np.ones(ncpl), fmt='%15.6E')   # pyemu: fixed width
    return ncpl


# ----------------------------------------------------------------- base

def step_base(org, python=None, calib_start=CALIB_START):
    import pandas as pd
    import pest_lamata_forward as pf
    import pest_lamata_obs as po
    import marmites_config as mcfg
    cfg = mcfg.load_run_config(os.path.join(org, 'pest_lamata.toml'))
    exports = pf.run_model(org, python=python)
    if exports is None:
        raise SystemExit('the base run failed: see %s'
                         % os.path.join(org, 'model.log'))
    modis = pd.read_csv(os.path.join(org, 'modis_et_catchment.csv'))
    ec = pd.read_csv(os.path.join(org, 'ec_et.csv'))
    tg = po.build_targets(exports, calib_start, modis=modis, ec=ec)
    sim = po.collect(exports, tg)
    # the base run's exports, kept for the record, then org trimmed
    keep = os.path.join(org, 'base_exports')
    if os.path.isdir(keep):
        shutil.rmtree(keep)
    shutil.copytree(exports, keep)
    _keep_only_inputs(org, cfg.spinup.strt_heads)
    _write_unit_mults(org, cfg.layers.nlay)
    tg['weight'] = po.balance_weights(tg, sim)
    tg.to_csv(os.path.join(org, pf.TARGETS), index=False,
              float_format='%.10g')
    sim.to_csv(os.path.join(org, pf.OBS_OUT), index=False,
               float_format='%.10g')
    rep = tg.assign(sim=sim['value'].to_numpy())
    rep['wres2'] = ((rep['obsval'] - rep['sim']) * rep['weight']) ** 2
    print('targets: %d, fitted %d' % (len(tg), int((tg['weight'] > 0).sum())))
    print(rep.groupby('obgnme').agg(n=('obsnme', 'size'),
                                    fitted=('weight', lambda w: (w > 0).sum()),
                                    phi=('wres2', 'sum')).to_string())
    return tg


# ------------------------------------------------------------------ pst

def step_pst(org, pest_ws, nreal, python=None):
    import warnings
    warnings.filterwarnings('ignore', message='.*Stringified.*')
    import flopy
    import pandas as pd
    import pyemu
    import marmites_config as mcfg
    import mm_paths
    import pest_lamata_forward as pf
    cfg = mcfg.load_run_config(os.path.join(org, 'pest_lamata.toml'))
    nlay = cfg.layers.nlay
    tpl = os.path.join(pest_ws, 'template')
    grb = flopy.mf6.utils.MfGrdFile(os.path.join(org, GRB))
    mg = grb.modelgrid
    idomain = np.asarray(grb.idomain).reshape(nlay, -1)
    pf_ = pyemu.utils.PstFrom(original_d=org, new_d=tpl, remove_existing=True,
                              spatial_reference=mg, longnames=True,
                              zero_based=False, echo=False)
    gs = pyemu.geostats.GeoStruct(
        variograms=pyemu.geostats.ExpVario(contribution=1.0, a=VARIO_A),
        transform='log', name='pp_gs')
    for stem, lo, hi in (('k', 0.1, 10.0), ('sy', 0.1, 10.0)):
        for L in range(nlay):
            pf_.add_parameters(
                filenames=os.path.join(MULT_DIR, '%s_l%d.txt' % (stem, L + 1)),
                par_type='pilotpoints', par_name_base='%s%d' % (stem, L + 1),
                pargp='%s_l%d' % (stem, L + 1), lower_bound=lo,
                upper_bound=hi, geostruct=gs,
                zone_array=(idomain[L] > 0).astype(int),
                pp_options={'pp_space': PP_SPACE, 'try_use_ppu': False})
    sc = pf.read_scalars(org)
    for name, value in zip(sc['name'], sc['value']):
        trans, lo, hi = scalar_bounds(name, float(value))
        pf_.add_parameters(
            filenames=os.path.join(pf.PARS, '%s.csv' % name), par_type='grid',
            par_style='d', index_cols=['scalar'], use_cols=['value'],
            par_name_base=name, pargp=name.rstrip('0123456789')
            .rstrip('_lh') or name, transform=trans, lower_bound=lo,
            upper_bound=hi)
    # the observations: our own instruction file, our own names
    tg = pd.read_csv(os.path.join(org, pf.TARGETS))
    ins = os.path.join(tpl, pf.OBS_OUT + '.ins')
    with open(ins, 'w') as fh:
        fh.write('pif ~\nl1\n')
        for n in tg['obsnme']:
            fh.write('l1 ~,~ !%s!\n' % n)
    pf_.add_observations_from_ins(ins, out_file=os.path.join(tpl,
                                                             pf.OBS_OUT),
                                  pst_path='.')
    exe = python or str(mm_paths.PYTHON_EXE) or sys.executable
    pf_.mod_sys_cmds.append('%s forward_lamata.py' % exe)
    pst = pf_.build_pst(version=2)
    od = pst.observation_data
    t = tg.set_index('obsnme')
    od.loc[t.index, 'obsval'] = t['obsval']
    od.loc[t.index, 'weight'] = t['weight']
    od.loc[t.index, 'obgnme'] = t['obgnme']
    od.loc[t.index, 'standard_deviation'] = t['standard_deviation']
    # the forward run the workers call (PstFrom's own applies the
    # parameters, then runs this)
    with open(os.path.join(tpl, 'forward_lamata.py'), 'w') as fh:
        fh.write('import os, sys\n'
                 'sys.path.insert(0, os.path.join(os.path.dirname('
                 'os.path.abspath(__file__)), "code", "tools"))\n'
                 'sys.path.insert(0, %r)\n'
                 'from pest_lamata_forward import forward\n'
                 'if __name__ == "__main__":\n'
                 '    import multiprocessing\n'
                 '    multiprocessing.freeze_support()\n'
                 '    sys.exit(forward("."))\n'
                 % os.path.join(pest_ws, 'code', 'tools'))
    pst.control_data.noptmax = 0
    pst.pestpp_options.update({
        'ies_num_reals': int(nreal), 'ies_parameter_ensemble': 'prior_pe.jcb',
        'ies_autoadaloc': True, 'ies_bad_phi_sigma': 2.5,
        'ies_lambda_mults': '0.1,1.0,10.0', 'lambda_scale_fac': '0.75,1.0',
        'ies_save_rescov': False, 'save_binary': True,
        'overdue_giveup_fac': 1.5, 'max_run_fail': 1,
        'panther_agent_freeze_on_fail': False})
    pe = pf_.draw(num_reals=int(nreal), use_specsim=False)
    pe.enforce()
    pe.to_binary(os.path.join(tpl, 'prior_pe.jcb'))
    # PstFrom's metadata brings a 'name' column beside parnme/obsnme -- an
    # alias the version-2 reader refuses to see twice
    for df in (pst.parameter_data, pst.observation_data):
        if 'name' in df.columns:
            df.drop(columns='name', inplace=True)
    pst.write(os.path.join(tpl, 'lamata.pst'), version=2)
    # the code snapshot travels with the template, so every worker has it
    shutil.copytree(os.path.join(pest_ws, 'code'), os.path.join(tpl, 'code'),
                    dirs_exist_ok=True)
    print('template %s: %d parameter(s) (%s), %d observation(s), %d with '
          'weight; prior ensemble of %d'
          % (tpl, pst.npar, ', '.join('%s %d' % (g, n) for g, n in
                                      pst.parameter_data.groupby('pargp')
                                      .size().items()),
             pst.nobs, int((pst.observation_data['weight'] > 0).sum()),
             int(nreal)))
    return tpl


def main(argv=None):
    import mm_paths
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--config', required=True)
    ap.add_argument('--state', required=True,
                    help='the saved state every realisation starts from')
    ap.add_argument('--nreal', type=int, default=30)
    ap.add_argument('--step', choices=('org', 'base', 'pst', 'all'),
                    default='all')
    ap.add_argument('--calib-start', default=CALIB_START.isoformat(),
                    help='the first day fitted (default: after the warm-up '
                         'hydrological year)')
    a = ap.parse_args(argv)
    pest_ws = str(mm_paths.PEST_WS) if mm_paths.PEST_WS else ''
    _check_ws(pest_ws, mm_paths)
    a.config = os.path.abspath(a.config)
    os.chdir(pest_ws)            # PstFrom writes its log where it runs
    org = os.path.join(pest_ws, 'org')
    if a.step in ('org', 'all'):
        step_org(a.config, a.state, pest_ws)
    if a.step in ('base', 'all'):
        step_base(org, calib_start=dt.date.fromisoformat(a.calib_start))
    if a.step in ('pst', 'all'):
        step_pst(org, pest_ws, a.nreal)


if __name__ == '__main__':
    main()
