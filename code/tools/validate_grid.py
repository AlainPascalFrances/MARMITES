# -*- coding: utf-8 -*-
"""The WP1c.8 validation ladder: does a mesh run mean the same thing?

Three rungs, in increasing distance from the validated model:

  (a) DIS  vs  DISV-from-DIS      MUST be identical. The vertex grid is the
                                  same geometry re-expressed, so any
                                  difference is a bug in the projection, not
                                  discretisation.
  (b) DISV-from-DIS  vs  Voronoi  Different discretisation, so the water
                                  balance may move. The question is whether it
                                  moves by an amount a modeller would accept,
                                  and the tool states the number rather than
                                  hiding it behind a pass mark.
  (c) production Voronoi  vs  the validated multi-year DIS run
                                  Not run here: it needs the real 45-year run,
                                  which the modeller launches. This tool does
                                  the COMPARING -- point it at two finished
                                  workspaces.

Rung (c) is the one that decides whether Voronoi becomes the default, and it
cannot be automated away: every difference has to be attributed to
discretisation rather than to a fault, and that is a judgement.

Usage
-----
    # run rungs (a) and (b) end to end, short
    python code/tools/validate_grid.py --ladder --nsp 30

    # compare two runs that already exist (rung c, or any pair)
    python code/tools/validate_grid.py --compare <WS_A> <WS_B> \\
        --labels structured voronoi --report out.md
"""

import argparse
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'tests'), os.path.join(CODE, 'ppMF6'),
           os.path.join(CODE, 'MARMITESutilities')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# (short label, INDEX_MM key) -- the canonical catchment water-balance terms,
# the same list plot_water_budget draws, so the ladder and the figures cannot
# disagree about what the budget IS.
FLUX_KEYS = [
    ('P', 'iP'), ('Pe', 'iPe'), ('Ei', 'iEi'), ('Ro', 'iRo'), ('I', 'iI'),
    ('Eow', 'iEow'), ('ETsoil', 'iETsoil'), ('Eg', 'iEg'), ('Tg', 'iTg'),
    ('ETg', 'iETg'), ('Rp', 'iperc'), ('EXFg', 'iEXFg'),
    ('dSsoil', 'idSsoil'), ('dSsurf', 'idSsurf'),
]


class LadderError(Exception):
    pass


# --------------------------------------------------------------------- #
# reading a finished run
# --------------------------------------------------------------------- #

def read_run(ws, mode='lagged'):
    """Everything the ladder compares, from one finished workspace.

    Reads the coupled HDF5 only -- no model rebuild, no .ini -- so it works on
    any run that finished, including one launched from Spyder months ago.
    """
    import h5py
    from marmites_indices import INDEX_MM
    fn = os.path.join(ws, '_coupled_%s.h5' % mode)
    if not os.path.exists(fn):
        raise LadderError('no coupled results in %s (looked for %s)'
                          % (ws, os.path.basename(fn)))
    out = {'ws': ws, 'file': fn}
    with h5py.File(fn, 'r') as h:
        wb = np.asarray(h['wb_ts'][:])
        out['nper'] = int(wb.shape[0])
        out['ncell'] = int(np.asarray(h['perc']).shape[1])
        out['grid_shape'] = (tuple(np.asarray(h['grid_shape'][:]).tolist())
                             if 'grid_shape' in h else None)
        out['cell_ij'] = (np.asarray(h['cell_ij'][:]) if 'cell_ij' in h
                          else None)
        out['obs_names'] = ([n.decode() if isinstance(n, bytes) else str(n)
                             for n in h['obs_names'][:]]
                            if 'obs_names' in h else [])
        out['obs_ij'] = (np.asarray(h['obs_ij'][:]) if 'obs_ij' in h else None)
        # Heads at the observation points are not stored separately: `heads`
        # is (nper, ncell) over the MM cell list and `obs_idx` indexes into
        # it. Gathering here keeps the ladder reading only the results file.
        out['obs_heads'] = None
        if 'heads' in h and 'obs_idx' in h:
            hd = np.asarray(h['heads'][:])
            idx = np.asarray(h['obs_idx'][:], dtype=int)
            if hd.ndim == 2 and idx.size:
                out['obs_heads'] = hd[:, idx]        # (nper, nobs)
        # mm/day over the run -> annual-equivalent mm/y, exactly as
        # plot_water_budget's totals panel does it
        out['flux'] = {}
        for short, key in FLUX_KEYS:
            if key in INDEX_MM:
                out['flux'][short] = float(wb[:, INDEX_MM[key]].sum()
                                           / out['nper'] * 365.25)
    return out


def mass_balance(flux):
    """Closure of the soil water balance, mm/y.

    in  = P + EXFg
    out = Ei + Ro + Eow + ETsoil + ETg + Rp
    residual = in - out - dS
    """
    g = flux.get
    inflow = g('P', 0.0) + g('EXFg', 0.0)
    outflow = (g('Ei', 0.0) + g('Ro', 0.0) + g('Eow', 0.0) + g('ETsoil', 0.0)
               + g('ETg', 0.0) + g('Rp', 0.0))
    ds = g('dSsoil', 0.0) + g('dSsurf', 0.0)
    return {'in': inflow, 'out': outflow, 'dS': ds,
            'residual': inflow - outflow - ds}


# --------------------------------------------------------------------- #
# comparing
# --------------------------------------------------------------------- #

def spinup_dominated(flux, factor=3.0):
    """Is this run still dominated by its initial condition?

    On a cold start the water table sits above ground over much of the
    catchment, so the seepage face discharges enormously and that water runs
    off and re-infiltrates. La Mata at 10 stress periods from a flat start
    shows Ro and EXFg near 8000 mm/y against 237 mm/y of precipitation -- a
    recirculation, not a hydrology.

    Comparing two such runs measures how each GRID relaxes its initial
    condition, which is not what the ladder is asking. Returns the offending
    terms, so the tool can refuse to call it a pass or a fail.
    """
    p = abs(flux.get('P', 0.0))
    if p <= 0:
        return []
    bad = []
    for term in ('Ro', 'EXFg'):
        v = abs(flux.get(term, 0.0))
        if v > factor * p:
            bad.append((term, v, v / p))
    return bad


def compare(a, b, label_a='A', label_b='B'):
    """Term-by-term comparison of two runs. Returns rows + a summary."""
    rows = []
    for short, _key in FLUX_KEYS:
        va, vb = a['flux'].get(short), b['flux'].get(short)
        if va is None or vb is None:
            continue
        d = vb - va
        rel = (abs(d) / abs(va) * 100.0) if abs(va) > 1e-12 else float('nan')
        rows.append({'term': short, 'a': va, 'b': vb, 'diff': d, 'rel': rel})
    mba, mbb = mass_balance(a['flux']), mass_balance(b['flux'])
    # The comparison that matters is on the terms that actually carry water:
    # a 0.001 mm/y term differing by 300 % says nothing.
    big = [r for r in rows if abs(r['a']) >= 1.0]
    worst = max(big, key=lambda r: r['rel']) if big else None
    summary = {
        'label_a': label_a, 'label_b': label_b,
        'ncell_a': a['ncell'], 'ncell_b': b['ncell'],
        'nper_a': a['nper'], 'nper_b': b['nper'],
        'residual_a': mba['residual'], 'residual_b': mbb['residual'],
        'max_abs_diff': max((abs(r['diff']) for r in rows), default=0.0),
        'worst_term': (worst['term'] if worst else None),
        'worst_rel': (worst['rel'] if worst else 0.0),
        'identical': all(abs(r['diff']) < 1e-9 for r in rows),
        'spinup_a': spinup_dominated(a['flux']),
        'spinup_b': spinup_dominated(b['flux']),
    }
    return rows, summary


def compare_heads(a, b):
    """Mean modelled head at each observation POINT, both runs.

    Points are compared by NAME, not by cell: on different grids the same
    piezometer lives in different cells, which is the whole point.
    """
    if a.get('obs_heads') is None or b.get('obs_heads') is None:
        return []
    HA, HB = np.asarray(a['obs_heads']), np.asarray(b['obs_heads'])
    out = []
    for n, nm in enumerate(a['obs_names']):
        if nm not in b['obs_names'] or n >= HA.shape[1]:
            continue
        m = b['obs_names'].index(nm)
        if m >= HB.shape[1]:
            continue
        va, vb = float(np.nanmean(HA[:, n])), float(np.nanmean(HB[:, m]))
        out.append({'name': nm, 'a': va, 'b': vb, 'diff': vb - va})
    return out


def format_report(rows, summary, heads=None, rung=None, tol=None):
    L = []
    ttl = 'Rung %s' % rung if rung else 'Run comparison'
    L.append('## %s: %s vs %s' % (ttl, summary['label_a'], summary['label_b']))
    L.append('')
    L.append('| | %s | %s |' % (summary['label_a'], summary['label_b']))
    L.append('|---|---:|---:|')
    L.append('| active cells | %d | %d |' % (summary['ncell_a'],
                                             summary['ncell_b']))
    L.append('| stress periods | %d | %d |' % (summary['nper_a'],
                                               summary['nper_b']))
    L.append('| balance residual [mm/y] | %.4g | %.4g |'
             % (summary['residual_a'], summary['residual_b']))
    L.append('')
    L.append('| term | %s [mm/y] | %s [mm/y] | diff | diff %% |'
             % (summary['label_a'], summary['label_b']))
    L.append('|---|---:|---:|---:|---:|')
    for r in rows:
        L.append('| %s | %.3f | %.3f | %+.3f | %s |'
                 % (r['term'], r['a'], r['b'], r['diff'],
                    ('%.2f' % r['rel']) if np.isfinite(r['rel']) else '-'))
    L.append('')
    if heads:
        L.append('| obs point | head %s [m] | head %s [m] | diff [m] |'
                 % (summary['label_a'], summary['label_b']))
        L.append('|---|---:|---:|---:|')
        for h in heads:
            L.append('| %s | %.3f | %.3f | %+.3f |'
                     % (h['name'], h['a'], h['b'], h['diff']))
        L.append('')
    if summary['identical']:
        L.append('**Identical** to 1e-9 mm/y on every term.')
    else:
        L.append('Largest relative difference on a term carrying >= 1 mm/y: '
                 '**%s, %.2f %%**; largest absolute difference %.4g mm/y.'
                 % (summary['worst_term'], summary['worst_rel'],
                    summary['max_abs_diff']))
    spun = summary['spinup_a'] + summary['spinup_b']
    if spun:
        L.append('')
        L.append('> **Not a verdict: still spinning up.** '
                 + '; '.join('%s is %.0f mm/y, %.0fx precipitation'
                             % (t, v, r) for t, v, r in spun)
                 + '. On a cold start the water table sits above ground, the '
                   'seepage face discharges hugely, and that water runs off '
                   'and re-infiltrates. Comparing two such runs measures how '
                   'each GRID relaxes its initial condition, not its '
                   'hydrology. Re-run from an equilibrated state before '
                   'reading these numbers.')
    if tol is not None:
        L.append('')
        if spun and tol != 0:
            L.append('Tolerance %.2f %%: **INCONCLUSIVE** (see above).' % tol)
        else:
            ok = (summary['identical'] if tol == 0
                  else summary['worst_rel'] <= tol)
            L.append('Tolerance %s: **%s**'
                     % ('exact' if tol == 0 else '%.2f %%' % tol,
                        'PASS' if ok else 'FAIL'))
    return '\n'.join(L)


# --------------------------------------------------------------------- #
# running a rung
# --------------------------------------------------------------------- #

def run_case(overrides, tag, nsp, config=None, python_exe=None, quiet=True):
    """Run the driver once and return its MF6 workspace.

    A SUBPROCESS, not an import: each case builds a different grid, and the
    driver caches meshes and model state on module-level state that is not
    designed to be set up twice in one interpreter.
    """
    import mm_paths
    config = config or os.path.join(CODE, 'configs', 'lamata.toml')
    exe = python_exe or sys.executable
    cmd = [exe, '-u', os.path.join(CODE, 'tests', 'run_lamata_mf6.py'),
           '--config', config, '--run-tag', tag]
    sets = dict(overrides)
    sets.setdefault('run.nsp', str(nsp))
    sets.setdefault('run.build_only', 'false')
    sets.setdefault('paths.libmf6', 'auto')
    sets.setdefault('spinup.strt_heads', '')
    sets.setdefault('spinup.steady_means', '')
    sets.setdefault('spinup.cycles', '1')
    sets.setdefault('postproc.enable', 'false')
    sets.setdefault('postproc.preproc', 'false')
    for k, v in sets.items():
        cmd += ['--set', '%s=%s' % (k, v)]
    env = dict(os.environ, MPLBACKEND='Agg')
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(mm_paths.REPO),
                       env=env)
    if r.returncode != 0:
        raise LadderError('case %r failed (exit %d):\n%s'
                          % (tag, r.returncode, (r.stdout or '')[-3000:]
                             + (r.stderr or '')[-2000:]))
    ws = None
    for line in (r.stdout or '').splitlines():
        if 'simulation written to' in line:
            ws = line.split('written to')[1].strip().split('  (')[0].strip()
    if ws is None:
        raise LadderError('could not find the workspace in the output of %r' % tag)
    if not quiet:
        print(r.stdout)
    return ws


RUNGS = {
    'a': dict(label_a='structured (DIS)', label_b='DISV-from-DIS',
              a={'grid.kind': 'structured'}, b={'grid.kind': 'disv'},
              tol=0.0,
              why='The vertex grid is the SAME geometry re-expressed, so any '
                  'difference is a projection bug, not discretisation.'),
    'b': dict(label_a='DISV-from-DIS', label_b='voronoi (uniform)',
              a={'grid.kind': 'disv'},
              b={'grid.kind': 'voronoi',
                 'grid.voronoi.stream_refine': 'false'},
              tol=15.0,
              why='A genuinely different discretisation. The water balance is '
                  'expected to move; the question is by how much.'),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ladder', action='store_true',
                    help='run rungs (a) and (b) end to end')
    ap.add_argument('--rungs', nargs='*', default=['a', 'b'], choices=['a', 'b'])
    ap.add_argument('--nsp', type=int, default=30)
    ap.add_argument('--config', default=None)
    ap.add_argument('--mode', default='lagged')
    ap.add_argument('--compare', nargs=2, metavar=('WS_A', 'WS_B'),
                    help='compare two finished workspaces (rung c)')
    ap.add_argument('--labels', nargs=2, default=['A', 'B'])
    ap.add_argument('--tol', type=float, default=None,
                    help='PASS if the worst relative difference is below this')
    ap.add_argument('--report', default=None, help='write markdown here')
    a = ap.parse_args()

    chunks = []
    if a.compare:
        ra = read_run(a.compare[0], a.mode)
        rb = read_run(a.compare[1], a.mode)
        rows, summary = compare(ra, rb, a.labels[0], a.labels[1])
        txt = format_report(rows, summary, compare_heads(ra, rb), tol=a.tol)
        print(txt)
        chunks.append(txt)
    elif a.ladder:
        for r in a.rungs:
            spec = RUNGS[r]
            print('=== rung %s: %s vs %s (nsp=%d) ==='
                  % (r, spec['label_a'], spec['label_b'], a.nsp))
            wsa = run_case(spec['a'], 'ladder_%s_a' % r, a.nsp, a.config)
            wsb = run_case(spec['b'], 'ladder_%s_b' % r, a.nsp, a.config)
            ra, rb = read_run(wsa, a.mode), read_run(wsb, a.mode)
            rows, summary = compare(ra, rb, spec['label_a'], spec['label_b'])
            txt = format_report(rows, summary, compare_heads(ra, rb),
                                rung=r, tol=spec['tol'])
            print(txt)
            chunks.append('*%s*\n\n' % spec['why'] + txt)
    else:
        ap.error('give --ladder or --compare')

    if a.report:
        with open(a.report, 'w', encoding='utf-8') as fh:
            fh.write('# WP1c.8 grid validation\n\n' + '\n\n---\n\n'.join(chunks)
                     + '\n')
        print('\nreport written to %s' % a.report)


if __name__ == '__main__':
    main()
