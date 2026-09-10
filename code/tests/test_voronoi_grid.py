# -*- coding: utf-8 -*-
"""WP1c.8 -- the grid validation ladder.

The ladder itself needs real runs, so what is tested here is the part that
decides what a comparison MEANS: the water-balance reduction, the term-by-term
comparison, and above all the guard that refuses to call a cold-start
comparison a pass or a fail.

That guard is the point of the whole file. Rung (b) on two 10-stress-period
runs reports a 73 % difference and would happily print FAIL -- but Ro and EXFg
are then 30x precipitation, so what it measured was two initial conditions
relaxing, not two discretisations.
"""

import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'tools'), os.path.join(CODE, 'ppMF6'), HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


VG = _load('validate_grid_t', os.path.join(CODE, 'tools', 'validate_grid.py'))


def _flux(**kw):
    """A plausible equilibrated La Mata balance, overridable per test."""
    base = {'P': 237.0, 'Pe': 217.0, 'Ei': 20.5, 'Ro': 5.0, 'I': 139.0,
            'Eow': 0.8, 'ETsoil': 120.0, 'Eg': 3.0, 'Tg': 2.0, 'ETg': 5.0,
            'Rp': 60.0, 'EXFg': 4.0, 'dSsoil': 45.0, 'dSsurf': 1.0}
    base.update(kw)
    return base


def _run(flux, ncell=1954, nper=365, names=None, heads=None):
    return {'ws': 'x', 'file': 'x', 'flux': flux, 'ncell': ncell, 'nper': nper,
            'obs_names': names or [], 'obs_heads': heads, 'obs_ij': None,
            'cell_ij': None, 'grid_shape': None}


# ------------------------------------------------------------- reduction
def test_mass_balance_closes_on_a_consistent_budget():
    f = {'P': 100.0, 'EXFg': 10.0, 'Ei': 5.0, 'Ro': 5.0, 'Eow': 0.0,
         'ETsoil': 40.0, 'ETg': 10.0, 'Rp': 30.0,
         'dSsoil': 18.0, 'dSsurf': 2.0}
    mb = VG.mass_balance(f)
    assert mb['in'] == 110.0
    assert mb['out'] == 90.0
    assert mb['dS'] == 20.0
    assert abs(mb['residual']) < 1e-12


def test_mass_balance_reports_a_gap_rather_than_hiding_it():
    f = {'P': 100.0, 'EXFg': 0.0, 'Ei': 0.0, 'Ro': 0.0, 'Eow': 0.0,
         'ETsoil': 0.0, 'ETg': 0.0, 'Rp': 0.0, 'dSsoil': 0.0, 'dSsurf': 0.0}
    assert VG.mass_balance(f)['residual'] == 100.0


# ------------------------------------------------------------ comparison
def test_identical_runs_are_reported_as_identical():
    """Rung (a)'s criterion: DISV-from-DIS is the same geometry re-expressed,
    so every term must match exactly."""
    f = _flux()
    rows, summary = VG.compare(_run(f), _run(dict(f)), 'dis', 'disv')
    assert summary['identical'] is True
    assert summary['max_abs_diff'] == 0.0
    assert all(r['diff'] == 0.0 for r in rows)


def test_a_single_changed_term_breaks_identity():
    """Guard against a vacuous pass."""
    a, b = _flux(), _flux(Rp=60.001)
    _rows, summary = VG.compare(_run(a), _run(b))
    assert summary['identical'] is False


def test_the_worst_term_ignores_fluxes_that_carry_no_water():
    """A 0.001 mm/y term differing by 300 % says nothing about the grid, and
    would otherwise always be 'the worst'."""
    a = _flux(Eow=0.001, ETsoil=120.0)
    b = _flux(Eow=0.004, ETsoil=126.0)          # +300 % vs +5 %
    _rows, summary = VG.compare(_run(a), _run(b))
    assert summary['worst_term'] == 'ETsoil'
    assert 4.0 < summary['worst_rel'] < 6.0


def test_comparison_records_both_grids_shape():
    _rows, s = VG.compare(_run(_flux(), ncell=1954), _run(_flux(), ncell=469),
                          'dis', 'voronoi')
    assert (s['ncell_a'], s['ncell_b']) == (1954, 469)


# --------------------------------------------------------- the spin-up guard
def test_a_cold_start_is_detected():
    """Observed on La Mata at nsp=10: Ro 7309 mm/y and EXFg 8027 mm/y against
    237 mm/y of precipitation."""
    bad = VG.spinup_dominated(_flux(Ro=7309.0, EXFg=8027.0))
    terms = [t for t, _v, _r in bad]
    assert terms == ['Ro', 'EXFg']
    assert all(r > 30 for _t, _v, r in bad)


def test_an_equilibrated_balance_is_not_flagged():
    assert VG.spinup_dominated(_flux()) == []


def test_a_spinning_up_comparison_is_inconclusive_not_a_failure():
    """The whole point: the tool must not print FAIL for a number that says
    nothing about the grid."""
    a = _flux(Ro=7309.0, EXFg=8027.0)
    b = _flux(Ro=5543.0, EXFg=6189.0)
    rows, summary = VG.compare(_run(a), _run(b), 'disv', 'voronoi')
    txt = VG.format_report(rows, summary, tol=15.0)
    assert 'INCONCLUSIVE' in txt
    assert 'FAIL' not in txt
    assert 'still spinning up' in txt


def test_an_equilibrated_comparison_does_give_a_verdict():
    a, b = _flux(), _flux(ETsoil=126.0)
    rows, summary = VG.compare(_run(a), _run(b))
    assert 'PASS' in VG.format_report(rows, summary, tol=15.0)
    assert 'FAIL' in VG.format_report(rows, summary, tol=1.0)


def test_rung_a_demands_exactness_even_while_spinning_up():
    """Identity does not depend on equilibrium: the same geometry must give
    the same answer whatever state the run is in, so rung (a) still gets a
    verdict."""
    f = _flux(Ro=7309.0, EXFg=8027.0)
    rows, summary = VG.compare(_run(f), _run(dict(f)))
    txt = VG.format_report(rows, summary, rung='a', tol=0.0)
    assert 'PASS' in txt and 'INCONCLUSIVE' not in txt


# ------------------------------------------------------------------ heads
def test_heads_are_matched_by_NAME_not_by_position():
    """On different grids a piezometer lives in a different cell, and the obs
    lists can differ in order or length -- matching by index would silently
    compare two different boreholes."""
    ha = np.array([[700.0, 800.0]])              # (nper, nobs)
    hb = np.array([[801.0, 701.0]])
    a = _run(_flux(), names=['A', 'B'], heads=ha)
    b = _run(_flux(), names=['B', 'A'], heads=hb)
    got = {d['name']: d['diff'] for d in VG.compare_heads(a, b)}
    assert abs(got['A'] - 1.0) < 1e-9
    assert abs(got['B'] - 1.0) < 1e-9


def test_heads_are_skipped_when_a_run_has_none():
    a = _run(_flux(), names=['A'], heads=None)
    assert VG.compare_heads(a, a) == []


def test_a_point_present_in_only_one_run_is_dropped():
    a = _run(_flux(), names=['A', 'B'], heads=np.array([[1.0, 2.0]]))
    b = _run(_flux(), names=['A'], heads=np.array([[1.5]]))
    got = VG.compare_heads(a, b)
    assert [d['name'] for d in got] == ['A']


# ------------------------------------------------------------------- I/O
def test_read_run_reduces_a_results_file(tmp_path):
    h5py = pytest.importorskip('h5py')
    from marmites_indices import INDEX_MM
    nper, nidx = 10, max(INDEX_MM.values()) + 1
    wb = np.zeros((nper, nidx))
    wb[:, INDEX_MM['iP']] = 2.0                  # 2 mm/d every day
    fn = tmp_path / '_coupled_lagged.h5'
    with h5py.File(fn, 'w') as h:
        h.create_dataset('wb_ts', data=wb)
        h.create_dataset('perc', data=np.zeros((nper, 7)))
        h.create_dataset('heads', data=np.full((nper, 7), 700.0))
        h.create_dataset('obs_idx', data=np.array([1, 3]))
        h.create_dataset('obs_names', data=np.array([b'A', b'B']))
    r = VG.read_run(str(tmp_path))
    assert r['nper'] == nper and r['ncell'] == 7
    # 2 mm/d -> annual-equivalent
    assert abs(r['flux']['P'] - 2.0 * 365.25) < 1e-6
    assert r['obs_names'] == ['A', 'B']
    assert r['obs_heads'].shape == (nper, 2)


def test_read_run_says_which_file_was_missing(tmp_path):
    with pytest.raises(VG.LadderError) as e:
        VG.read_run(str(tmp_path))
    assert '_coupled_lagged.h5' in str(e.value)


def test_the_rungs_are_declared_with_their_tolerances():
    assert VG.RUNGS['a']['tol'] == 0.0, 'rung (a) must demand exactness'
    assert VG.RUNGS['b']['tol'] > 0.0
    for spec in VG.RUNGS.values():
        assert spec['why'] and spec['label_a'] and spec['label_b']
