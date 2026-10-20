# -*- coding: utf-8 -*-
"""Set up the La Mata PEST++-IES calibration (WP7) in the PEST workspace.

    python code/tools/pest_lamata_setup.py --config code/configs/lamata.toml \\
        [--step all]

The calibration's choices are the configuration's ``[pest]`` section,
edited on the PEST++ set-up page (docs/MM-MF6/WP7_PEST_frontend_plan.md):
the saved state every realisation starts from (``base_state``), the first
day fitted (``calib_start``), the parameter table (``par``: empty means
:func:`default_pars`), the pilot points, the observation groups and the
ensemble. The PEST workspace is the one set on Home (``pest_ws``;
mm_paths.PEST_WS), never WS_ROOT and never the repository. Steps
(``--step``):

  org     the model as every worker will run it, in ``<pest_ws>/org``:
            code/         a SNAPSHOT of the repository's code, so a
                          calibration of days is not changed under it
            pest_lamata.toml  the configuration copy: MMsurf off, no
                          figures, one cycle from the saved state, nothing
                          saved, [pest] multipliers = mult
            example/LaMata the dataset copy (a run rewrites the forcing
                          files there), the observation table with C4, C5
                          and W1 switched on; the measured ET at the points
                          (<obs.aet_prefix>_<point>.txt) is read by the
                          runs themselves, into their obs_et.csv
            ws/MF6_ws_voronoi the mesh and overlay caches and the saved
                          state every realisation starts from
            mm_mult/         k_l<L>, sy_l<L>, ss_l<L> per cell (the pilot
                          points)
            pars/<name>.csv  the one-value parameters (soil horizons, SFR,
                          LAK, UZF, CRR, Sr and DRN per layer)
            modis_et_catchment.csv  MODIS's catchment ET per composite
  base    the base run in org (the test run of the setup, ~45 min), then
          the targets (pest_targets.csv) and their weights, balanced so each
          group carries its share of the initial phi
  pst     PstFrom -> ``<pest_ws>/template``: pilot points for K, Sy and Ss
          per layer on the mesh, the one-value parameters, the
          observations, the prior ensemble drawn from the geostatistical
          structure, pestpp-ies options; noptmax = 0 until a launch
          (pest_lamata_run.py) sets it

Model runs happen only here, in the PEST workspace (CLAUDE.md rule 1, the
PEST++ exception of 2026-10-11).
"""

import argparse
import dataclasses
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

import marmites_config as mcfg            # noqa: E402
from pest_lamata_forward import MULT_DIR  # noqa: E402


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


def points_on(cfg):
    """The commented-out observation points the calibration switches on
    ([pest] points_on, comma- or space-separated)."""
    return tuple(n for n in cfg.pest.points_on.replace(',', ' ').split() if n)


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


# THE DEFAULT PARAMETER TABLE -- what an empty [pest] par means
# (default_pars): the calibration of 2026-10-18.
#
# starting values that differ from the configuration's: inside their
# ranges, so the prior does not pile on a bound (64 % and 57 % of the first
# prior's draws sat on MF6's eps minimum and on beta = 1; user, 2026-10-12)
PRIOR_START = {'uzf_eps': 4.5, 'crr_beta': 0.95}

# THE HYDRAULIC PARAMETERS, from the user's table
# (lamata_sardon_hydraulic_parameters.xlsx, 2026-10-18) -- PEST only, the
# reference model keeps its own (Sy 0.01, thts 0.45, K 0.05 m/d).
#
# Sy and the specific retention Sr per layer, as water contents: (start,
# lower, upper). UZF's thts = Sy + Sr in every cell and thtr = Sr. The
# maxima widened so the table's medians sit inside (user): saprolite Sr 25
# -> 30 % (n up to 35 %), fissured Sy 1 -> 2 %.
SY_SR = {1: {'sy': (0.025, 0.005, 0.05), 'sr': (0.275, 0.045, 0.30)},
         2: {'sy': (0.010, 0.001, 0.02), 'sr': (0.040, 0.009, 0.09)}}
# K [m/d]: the reference's K is layer 1's start; the bounds span every K
# value gathered -- 2.3e-5 (the table's saprolite minimum) to 11.2
# (Legchenko et al. 2004,
# docs/MM-MF6/WP7_granite_porosity_specific_retention_review.md) -- split
# so layer 1's lie [pest] k_order_decades below layer 2's, and layer 2
# starts that far above layer 1. The model holds layer 1's K at most
# 10**-k_order_decades x layer 2's, cell by cell (marmites_pest).
K_RANGE = (2.3e-5, 11.2)
# The anisotropy Kh/Kv per layer, (start, lower, upper): the saprolite
# about 2 (weakly anisotropic, 1-3), the fissured layer about 10 (2-30,
# subhorizontal weathering fractures) -- user, 2026-10-19, from
# docs/MM-MF6/WP7_granite_anisotropy_review.md. A deeper layer takes the
# last. Calibrated only where layers.k33 is a value given as a ratio.
K33_RANGE = {1: (mcfg.K33_RATIO_DEFAULT[0], 1.0, 3.0),
             2: (mcfg.K33_RATIO_DEFAULT[1], 2.0, 30.0)}
# one-value parameters of a package that may be off
PACKAGE_OF = {'sfr_rhk': 'sfr', 'lak_bedleak': 'lak', 'crr_beta': 'crr'}


def _raster_median(path):
    """The median of an ASCII raster's valid cells."""
    a = np.loadtxt(path, skiprows=6, dtype=float)
    a = a[np.isfinite(a) & (a > -9999.0) & (a > 0.0)]
    if not a.size:
        raise SystemExit('%s has no valid cell' % path)
    return float(np.median(a))


def reference_median(cfg, what, layer, dataset_dir):
    """The reference value of a layer property (``layers.<what>``) that the
    calibration's values are relative to: the median of its raster's valid
    cells (uniform on La Mata, so exact), or its one value."""
    src = getattr(cfg.layers, what)
    names = src.rasters(int(cfg.layers.nlay))
    if names:
        return _raster_median(os.path.join(dataset_dir, names[layer - 1]))
    if src.value is not None:
        return float(src.values(int(cfg.layers.nlay))[layer - 1])
    raise SystemExit('layers.%s has neither a raster nor a value: the '
                     'calibration needs a reference' % what)


def _row(family, layer, start, lower, upper, log):
    """A parameter row (marmites_config.PestPar's fields), its start held
    inside its bounds."""
    start = min(max(float(start), float(lower)), float(upper))
    return {'family': family, 'layer': int(layer), 'enable': True,
            'start': start, 'lower': float(lower), 'upper': float(upper),
            'log': bool(log)}


def default_pars(cfg, dataset_dir):
    """The parameter table of the calibration of 2026-10-18, for this
    configuration -- what an EMPTY [pest] par means. K and Sy (pilot points)
    and Sr from the user's table (SY_SR, K_RANGE), Ss (pilot points) and the
    outlet drains as multipliers, each soil horizon's field capacity (its
    place between the wilting point and the porosity, +-0.2, never outside
    them) and Ks (x0.1 to x10), and the package scalars; the anisotropy Kh/Kv per layer (K33_RANGE)
    where layers.k33 is a ratio. A list of dicts."""
    nlay = int(cfg.layers.nlay)
    dec = float(cfg.pest.k_order_decades)
    rows = []
    for L in range(1, nlay + 1):
        k0 = reference_median(cfg, 'k', L, dataset_dir)
        lo, hi = K_RANGE
        if nlay >= 2:
            if L == 1:
                hi = hi / 10.0 ** dec
            else:
                lo = lo * 10.0 ** dec
        rows.append(_row('k', L, k0 * (10.0 ** dec if L >= 2 else 1.0),
                         lo, hi, True))
    for L in range(1, nlay + 1):
        if L in SY_SR:
            rows.append(_row('sy', L, *SY_SR[L]['sy'], True))
        else:                          # a layer the table does not cover
            sy0 = reference_median(cfg, 'sy', L, dataset_dir)
            rows.append(_row('sy', L, sy0, 0.03 * sy0, min(30.0 * sy0, 0.5),
                             True))
    for L in range(1, nlay + 1):
        rows.append(_row('ss', L, 1.0, 0.1, 10.0, True))
    for L in range(1, nlay + 1):
        if L in SY_SR:                 # the table's layers; others keep thts
            rows.append(_row('sr', L, *SY_SR[L]['sr'], True))
    for L in range(1, nlay + 1):
        rows.append(_row('drn', L, 1.0, 0.01, 100.0, True))
    if cfg.layers.k33_as_ratio and cfg.layers.k33.producer() == 'value':
        for L in range(1, nlay + 1):
            rows.append(_row('k33', L, *K33_RANGE[min(L, max(K33_RANGE))],
                             True))
    for k, h in enumerate(cfg.soil.horizon, 1):
        # +-0.2 of its place (was 0.1: 20 % of the pilot on the lower
        # bound, 2026-10-13)
        w = h.smax - h.sr
        f = (h.sfc - h.sr) / w
        rows.append(_row('sfc', k, h.sfc, h.sr + max(f - 0.20, 0.01) * w,
                         h.sr + min(f + 0.20, 0.99) * w, False))
    for k, h in enumerate(cfg.soil.horizon, 1):
        rows.append(_row('ks', k, h.ks, 0.1 * h.ks, 10.0 * h.ks, True))
    # streambed K: at 1 m/d the model crawls (64 min alone, 63 linear
    # iterations a day) and no better a fit -- capped at 0.5 (2026-10-12)
    if cfg.sfr.enable and cfg.sfr.rhk.value is not None:
        rows.append(_row('sfr_rhk', 0, cfg.sfr.rhk.value, 0.001, 0.5, True))
    if cfg.lak.enable:
        rows.append(_row('lak_bedleak', 0, cfg.lak.bedleak, 1.0e-4, 1.0e-2,
                         True))
    rows.append(_row('uzf_eps', 0, PRIOR_START['uzf_eps'], 3.5, 8.0,
                     False))                         # 3.5: MF6's minimum
    rows.append(_row('uzf_vks', 0, cfg.uzf.vks_scale, 0.2, 5.0, True))
    if cfg.crr.enable:
        rows.append(_row('crr_beta', 0, PRIOR_START['crr_beta'], 0.8, 1.0,
                         False))
    return rows


def pars(cfg, dataset_dir):
    """The calibration's parameter table: [pest] par, or -- empty -- the
    defaults for this configuration. A list of dicts."""
    if cfg.pest.par:
        return [dataclasses.asdict(r) for r in cfg.pest.par]
    return default_pars(cfg, dataset_dir)


def pilot_specs(cfg, dataset_dir, rows):
    """(stem, layer, start, lower, upper, log) of the pilot-point families
    switched on, as multipliers of the reference: K in m/d and Sy as a water
    content, divided by the reference's median; Ss is a multiplier already
    (user, 2026-10-17: Ss with K and Sy, it was one value per layer)."""
    out = []
    for r in rows:
        if r['family'] not in mcfg.PEST_PILOT or not r['enable']:
            continue
        L = int(r['layer'])
        ref = (reference_median(cfg, r['family'], L, dataset_dir)
               if r['family'] in ('k', 'sy') else 1.0)
        out.append((r['family'], L, r['start'] / ref, r['lower'] / ref,
                    r['upper'] / ref, bool(r['log'])))
    return out


def scalar_specs(cfg, rows, out=print):
    """(name, start, transform, lower, upper) of the one-value parameters
    switched on, as the forward run reads them
    (pest_lamata_forward.apply_scalars): sfc_f_h<k> and ks_h<k> (horizon k),
    sr_l<L> and drn_l<L> (layer L), the package scalars by their family's
    name. A field capacity is calibrated as its place between the wilting
    point and the porosity, so no draw can leave them (user, 2026-10-17). A
    package that is off has nothing to calibrate."""
    res = []
    for r in rows:
        fam = r['family']
        if fam in mcfg.PEST_PILOT or not r['enable']:
            continue
        pkg = PACKAGE_OF.get(fam)
        if pkg and not getattr(cfg, pkg).enable:
            out('   %s: %s is off, not calibrated' % (fam, pkg.upper()))
            continue
        L = int(r['layer'])
        trans = 'log' if r['log'] else 'none'
        s, lo, hi = float(r['start']), float(r['lower']), float(r['upper'])
        if fam == 'sfc':
            h = cfg.soil.horizon[L - 1]
            w = h.smax - h.sr
            res.append(('sfc_f_h%d' % L, (s - h.sr) / w, trans,
                        (lo - h.sr) / w, (hi - h.sr) / w))
        elif fam == 'ks':
            res.append(('ks_h%d' % L, s, trans, lo, hi))
        elif fam in ('sr', 'drn', 'k33'):
            res.append(('%s_l%d' % (fam, L), s, trans, lo, hi))
        else:
            res.append((fam, s, trans, lo, hi))
    return res


def scalar_table(specs):
    """The one-value parameters' starting values (name, value), as the
    forward run reads them."""
    import pandas as pd
    return pd.DataFrame([(sp[0], sp[1]) for sp in specs],
                        columns=['name', 'value'])


def calib_start(cfg, dataset_dir):
    """The first day fitted: [pest] calib_start, or -- blank -- the first
    hydrological year that starts after one whole year of record, the
    warm-up (La Mata: 2009-10-01, the record starting on 2008-05-31)."""
    s = (cfg.pest.calib_start or '').strip()
    if s:
        return dt.date.fromisoformat(s)
    first = None
    with open(os.path.join(dataset_dir, 'inputDATE.txt')) as fh:
        for ln in fh:
            if ln.strip() and not ln.lstrip().startswith('#'):
                first = dt.date.fromisoformat(ln.split(',')[0].strip()[:10])
                break
    if first is None:
        raise SystemExit('no date in %s'
                         % os.path.join(dataset_dir, 'inputDATE.txt'))
    try:
        after = first.replace(year=first.year + 1)
    except ValueError:                                  # 29 February
        after = first.replace(year=first.year + 1, day=28)
    m = int(cfg.postproc.hydro_year_start)
    d = dt.date(after.year, m, 1)
    return d if d >= after else dt.date(after.year + 1, m, 1)


def k_order_prior(df, pdata, decades=1.0):
    """The prior's layer-1 K ``decades`` under layer 2's at each pair of
    pilot points at the same place -- the model holds it cell by cell, and
    this keeps the prior from spending realisations on draws it would clip:
    where log10 k1 > log10 k2 - decades, both move halfway towards each
    other. ``df`` (realisations x parameters, multipliers of rasters with
    the same median) is changed in place. Returns how many draws moved."""
    import re

    def xy(name):
        m = re.search(r'x:([-\d.]+)_y:([-\d.]+)', name)
        return ((round(float(m.group(1)), 2), round(float(m.group(2)), 2))
                if m else None)
    grp = pdata['pargp']
    k1 = {xy(n): n for n in df.columns if grp.get(n) == 'k_l1' and xy(n)}
    k2 = {xy(n): n for n in df.columns if grp.get(n) == 'k_l2' and xy(n)}
    moved = 0
    for place, n1 in k1.items():
        n2 = k2.get(place)
        if n2 is None:
            continue
        a = np.log10(df[n1].to_numpy(dtype=float))
        b = np.log10(df[n2].to_numpy(dtype=float))
        d = a - (b - decades)
        bad = d > 0.0
        if bad.any():
            df.loc[bad, n1] = 10.0 ** (a[bad] - d[bad] / 2.0 - 1e-9)
            df.loc[bad, n2] = 10.0 ** (b[bad] + d[bad] / 2.0 + 1e-9)
            moved += int(bad.sum())
    return moved


def _print_pars(rows):
    print('   parameters (%d row(s), %d on):' % (len(rows), sum(
        1 for r in rows if r['enable'])))
    for r in rows:
        print('      %-3s %-12s %2d  start %-10.4g [%.4g, %.4g]%s'
              % ('on' if r['enable'] else 'off', r['family'], r['layer'],
                 r['start'], r['lower'], r['upper'],
                 '  log' if r['log'] else ''))


def resolve_base(cfg, ws_root):
    """The model calibrated and the state every realisation starts from:
    ``(model configuration, state name, state folder, what they are)``.

    [pest] base_run, a finished run of the registry: its configuration copy
    is the model and its end state the start (marmites_runinfo.base_of),
    with THIS configuration's [pest] -- the calibration's settings are not
    the run's. Blank: this configuration as it is and base_state, a saved
    state chosen by hand, in its state workspace."""
    import copy
    import marmites_runinfo as ri
    run = (cfg.pest.base_run or '').strip()
    if not run:
        state = (cfg.pest.base_state or '').strip()
        if not state:
            raise SystemExit('no base for the realisations: choose a run '
                             'solution (pest.base_run) or a saved state '
                             '(pest.base_state) on the PEST++ set-up page')
        return (cfg, state, mcfg.state_workspace(cfg, ws_root),
                'this configuration and the saved state %s, chosen by hand'
                % state)
    b = ri.base_of(os.path.join(ri.runs_dir(cfg, ws_root), run), ws_root)
    if b['why']:
        raise SystemExit('run %s cannot be the base: %s' % (run, b['why']))
    model = mcfg.load_run_config(b['config'])
    model.pest = copy.deepcopy(cfg.pest)
    return (model, b['state'], b['state_dir'],
            'run %s: its configuration (%s) and its end state %s%s'
            % (run, b['config'], b['state'],
               '' if b.get('converged', True) else
               ' -- a spin-up saved BEFORE it converged'))


def step_org(cfg_path, pest_ws):
    import mm_paths
    import pest_modis_et as pme
    cfg, state, src, what = resolve_base(mcfg.load_run_config(cfg_path),
                                         mm_paths.WS_ROOT)
    org = os.path.join(pest_ws, 'org')
    os.makedirs(org, exist_ok=True)
    print('PEST workspace %s' % pest_ws)
    print('   base: %s' % what)
    _snapshot_code(os.path.join(pest_ws, 'code'))
    print('   code snapshot: %s' % os.path.join(pest_ws, 'code'))
    case = cfg.paths.case
    # the parameter table, explicit in the copy: every later step reads
    # the same rows, whatever the defaults would be by then
    rows = pars(cfg, str(mm_paths.dataset_dir(case)))
    cfg.pest.par = [mcfg.PestPar(**r) for r in rows]
    _print_pars(rows)
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
    # a realisation that stalls stops as non-convergence at once and
    # counts as failed (user, 2026-10-20)
    cfg.run.max_solver_iterations = int(cfg.pest.max_solver_iterations)
    # where the saved state carries no UZF water content, the unsaturated
    # zone starts EMPTY, as the reference model's does (its thti 0.15 is
    # under its thtr 0.44, so the build raises it to thtr); a realisation's
    # thts can fall below 0.15 -- the fissured layer's 0.05 -- which would
    # lower it to thts and start that zone full (2026-10-18)
    cfg.uzf.thti.value = 1.0e-3
    cfg.validate()
    cfg.write_toml(os.path.join(org, 'pest_lamata.toml'))
    # the dataset copy, with the points the calibration adds
    ds = os.path.join(org, 'example', case)
    if os.path.isdir(ds):
        shutil.rmtree(ds)
    shutil.copytree(str(mm_paths.dataset_dir(case)), ds)
    on = _obs_table_with(os.path.join(ds, cfg.obs.table), points_on(cfg))
    print('   dataset copy, observation points switched on: %s'
          % (', '.join(on) or 'none'))
    print('   first day fitted: %s' % calib_start(cfg, ds))
    # the caches and the saved state, from the base's state workspace
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
    pfw.write_scalars(org, scalar_table(scalar_specs(cfg, rows)))
    # the order of the layers' K, held by the model (marmites_pest)
    dec = float(cfg.pest.k_order_decades)
    if cfg.layers.nlay >= 2 and dec > 0.0:
        with open(os.path.join(md, 'k_l1_max_ratio.txt'), 'w') as fh:
            fh.write('%.10g\n' % 10.0 ** -dec)
    # the base run at the calibration's starting values: K, Sy and Ss, one
    # value per layer (step_base then writes the per-cell ones the pilot
    # points multiply)
    for stem, L, start, _lo, _hi, _log in pilot_specs(cfg, ds, rows):
        with open(os.path.join(md, '%s_l%d.txt' % (stem, L)), 'w') as fh:
            fh.write('%.10g\n' % start)
    # MODIS's catchment ET
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
    json.dump({'config': os.path.abspath(cfg_path), 'state': state,
               'state_dir': str(src), 'base': what,
               'base_run': cfg.pest.base_run,
               'model_config': cfg.source_path,
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


def _write_unit_mults(org, specs):
    """The pilot-point families' files (k_l<L>, sy_l<L>, ss_l<L> of
    ``specs``): one 1.0 per mesh cell (the grid file gives ncpl)."""
    import flopy
    ncpl = flopy.mf6.utils.MfGrdFile(os.path.join(org, GRB)).modelgrid.ncpl
    for stem, L, *_rest in specs:
        np.savetxt(os.path.join(org, MULT_DIR, '%s_l%d.txt' % (stem, L)),
                   np.ones(ncpl), fmt='%15.6E')   # pyemu: fixed width
    return ncpl


def _org_dataset(org, cfg):
    return os.path.join(org, 'example', cfg.paths.case)


# ----------------------------------------------------------------- base

def _targets(org, cfg, exports):
    """The targets, the model's values of them and the balanced weights,
    from a run's exports -- the [pest] observation settings of the copy."""
    import pandas as pd
    import pest_lamata_obs as po
    modis = pd.read_csv(os.path.join(org, 'modis_et_catchment.csv'))
    tg = po.build_targets(exports, calib_start(cfg, _org_dataset(org, cfg)),
                          modis=modis,
                          hy_month=int(cfg.postproc.hydro_year_start),
                          pest=cfg.pest)
    sim = po.collect(exports, tg, pest=cfg.pest)
    tg['weight'] = po.balance_weights(tg, sim, pest=cfg.pest)
    return tg, sim


def step_targets(org):
    """The targets and their weights again, from the base run's exports
    kept in org/base_exports -- no model run; the scalar files take the
    copy's starting values."""
    import pest_lamata_forward as pf
    ex = os.path.join(org, 'base_exports')
    if not os.path.isdir(ex):
        raise SystemExit('no base exports in %s: run --step base' % org)
    cfg = mcfg.load_run_config(os.path.join(org, 'pest_lamata.toml'))
    shutil.rmtree(os.path.join(org, pf.PARS), ignore_errors=True)
    rows = pars(cfg, _org_dataset(org, cfg))
    pf.write_scalars(org, scalar_table(scalar_specs(cfg, rows)))
    tg, sim = _targets(org, cfg, ex)
    tg.to_csv(os.path.join(org, pf.TARGETS), index=False,
              float_format='%.10g')
    sim.to_csv(os.path.join(org, pf.OBS_OUT), index=False,
               float_format='%.10g')
    _report_targets(tg, sim)
    return tg


def _report_targets(tg, sim):
    rep = tg.assign(sim=sim['value'].to_numpy())
    rep['wres2'] = ((rep['obsval'] - rep['sim']) * rep['weight']) ** 2
    print('targets: %d, fitted %d' % (len(tg), int((tg['weight'] > 0).sum())))
    print(rep.groupby('obgnme').agg(
        n=('obsnme', 'size'), fitted=('weight', lambda w: (w > 0).sum()),
        weight=('weight', 'max'), phi=('wres2', 'sum')).to_string())


def step_base(org, python=None):
    import pest_lamata_forward as pf
    cfg = mcfg.load_run_config(os.path.join(org, 'pest_lamata.toml'))
    exports = pf.run_model(org, python=python)
    if exports is None:
        raise SystemExit('the base run failed: see %s'
                         % os.path.join(org, 'model.log'))
    tg, sim = _targets(org, cfg, exports)
    # the base run's exports, kept for the record, then org trimmed
    keep = os.path.join(org, 'base_exports')
    if os.path.isdir(keep):
        shutil.rmtree(keep)
    shutil.copytree(exports, keep)
    _keep_only_inputs(org, cfg.spinup.strt_heads)
    _write_unit_mults(org, pilot_specs(cfg, _org_dataset(org, cfg),
                                       pars(cfg, _org_dataset(org, cfg))))
    tg.to_csv(os.path.join(org, pf.TARGETS), index=False,
              float_format='%.10g')
    sim.to_csv(os.path.join(org, pf.OBS_OUT), index=False,
               float_format='%.10g')
    _report_targets(tg, sim)
    return tg


# ------------------------------------------------------------------ pst

def step_pst(org, pest_ws, python=None):
    import warnings
    warnings.filterwarnings('ignore', message='.*Stringified.*')
    import flopy
    import pandas as pd
    import pyemu
    import mm_paths
    import pest_lamata_forward as pf
    from pest_lamata_run import pestpp_options
    cfg = mcfg.load_run_config(os.path.join(org, 'pest_lamata.toml'))
    p = cfg.pest
    nlay = cfg.layers.nlay
    nreal = int(p.nreal)
    tpl = os.path.join(pest_ws, 'template')
    grb = flopy.mf6.utils.MfGrdFile(os.path.join(org, GRB))
    mg = grb.modelgrid
    idomain = np.asarray(grb.idomain).reshape(nlay, -1)
    pf_ = pyemu.utils.PstFrom(original_d=org, new_d=tpl, remove_existing=True,
                              spatial_reference=mg, longnames=True,
                              zero_based=False, echo=False)
    gs = pyemu.geostats.GeoStruct(
        variograms=pyemu.geostats.ExpVario(contribution=1.0,
                                           a=float(p.vario_range)),
        transform='log', name='pp_gs')
    # K, Sy and Ss per layer on pilot points (the parameter table)
    rows = pars(cfg, _org_dataset(org, cfg))
    for stem, L, start, lo, hi, log in pilot_specs(cfg, _org_dataset(org, cfg),
                                                   rows):
        print('   %s layer %d: start x%.4g, bounds x%.4g .. x%.4g'
              % (stem, L, start, lo, hi))
        pf_.add_parameters(
            filenames=os.path.join(MULT_DIR, '%s_l%d.txt' % (stem, L)),
            par_type='pilotpoints', par_name_base='%s%d' % (stem, L),
            pargp='%s_l%d' % (stem, L), lower_bound=lo, upper_bound=hi,
            initial_value=start, geostruct=gs,
            transform='log' if log else 'none',
            zone_array=(idomain[L - 1] > 0).astype(int),
            pp_options={'pp_space': float(p.pp_space), 'try_use_ppu': False})
    for name, _start, trans, lo, hi in scalar_specs(cfg, rows):
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
        'ies_num_reals': nreal, 'ies_parameter_ensemble': 'prior_pe.jcb',
        'ies_autoadaloc': True, 'lambda_scale_fac': '1.0',
        'ies_save_rescov': False, 'save_binary': True, 'max_run_fail': 1,
        'panther_agent_freeze_on_fail': False})
    # lambdas, bad phi, overdue limit: [pest], as a launch sets them again
    pst.pestpp_options.update(pestpp_options(p))
    pe = pf_.draw(num_reals=nreal, use_specsim=False)
    pe.enforce()
    dec = float(p.k_order_decades)
    if dec > 0.0:
        n_ord = k_order_prior(pe._df, pst.parameter_data, dec)
        if n_ord:
            print('   prior: layer 1 K brought %g decade(s) under layer 2 at '
                  '%d pilot-point draw(s)' % (dec, n_ord))
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
             nreal))
    return tpl


def main(argv=None):
    import mm_paths
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--config', required=True,
                    help='the configuration; its [pest] section is the '
                         'calibration (the PEST++ set-up page)')
    ap.add_argument('--step', choices=('org', 'base', 'targets', 'pst',
                                       'all'), default='all')
    a = ap.parse_args(argv)
    pest_ws = str(mm_paths.PEST_WS) if mm_paths.PEST_WS else ''
    _check_ws(pest_ws, mm_paths)
    a.config = os.path.abspath(a.config)
    os.chdir(pest_ws)            # PstFrom writes its log where it runs
    org = os.path.join(pest_ws, 'org')
    if a.step in ('org', 'all'):
        step_org(a.config, pest_ws)
    if a.step in ('base', 'all'):
        step_base(org)
    if a.step == 'targets':
        step_targets(org)
    if a.step in ('pst', 'all'):
        step_pst(org, pest_ws)


if __name__ == '__main__':
    main()
