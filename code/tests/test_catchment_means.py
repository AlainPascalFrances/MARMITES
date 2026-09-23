# -*- coding: utf-8 -*-
"""Catchment means weight each cell by its area; a spin-up says if it failed.

The run of 2026-09-23 on the Voronoi mesh: half the 15,750 cells cover 5 % of
La Mata, refined along the streams, and every "catchment mean" was a plain
mean over cells -- runoff 339 mm/yr against 61 area-weighted (more than the
318 mm/yr of rain), exfiltration 310 against 16, where MF6's own seepage
drain said 15.9. And the spin-up ran out of its 6 cycles with |dWT| still
0.38 m against a 0.05 m tolerance, and saved its heads as 'equilibrated'.
"""

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'MARMITESutilities', 'MARMITESplot'),
           os.path.join(CODE, 'MARMITESutilities')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _src(*parts):
    return open(os.path.join(CODE, *parts), encoding='utf-8').read()


def test_the_coupler_weights_the_catchment_series_by_area():
    src = _src('marmites_coupler.py')
    assert 'self.wb_ts[n] = mm_cells.mean(axis=0)' not in src
    assert ('self.wb_ts[n] = np.average(mm_cells, axis=0,\n'
            '                                           weights=self.area)') in src
    assert 'self.wb_ts_soil[n] = np.average(mms_cells, axis=0,' in src


def test_an_area_weighted_mean_is_what_a_small_corridor_cannot_skew():
    """Ten 1 m2 corridor cells with 1000 mm, one 10000 m2 hillslope cell
    with 10 mm: the catchment gets ~11, the plain mean 910."""
    v = np.array([1000.0] * 10 + [10.0])
    a = np.array([1.0] * 10 + [10000.0])
    assert v.mean() == pytest.approx(910.0)
    assert np.average(v, weights=a) == pytest.approx(10.99, abs=0.01)


def test_the_results_file_keeps_the_cell_areas_and_the_run_line_uses_them():
    src = _src('tests', 'run_lamata_mf6.py')
    assert "f.create_dataset('cell_area'" in src
    assert "np.average(res['perc'], axis=1, weights=_w)" in src
    assert "res['perc'].mean(), res['etg'].mean()" not in src


def test_the_budget_summary_and_the_coupling_figures_weight_by_area():
    pwb = _src('tests', 'plot_water_budget.py')
    assert "'cell_area'):" in pwb
    assert "np.average(new['rejinf'], axis=1, weights=w)" in pwb
    assert pwb.count("area=new.get('cell_area')") == 2


def test_the_plot_helper_weights_when_given_areas():
    mmplot = pytest.importorskip('MARMITESplot_v3')
    v = np.array([[1000.0] * 10 + [10.0]])
    a = np.array([1.0] * 10 + [10000.0])
    assert mmplot._cell_mean(v, a)[0] == pytest.approx(10.99, abs=0.01)
    assert mmplot._cell_mean(v)[0] == pytest.approx(910.0)
    assert mmplot._cell_mean(v, np.ones(3))[0] == pytest.approx(910.0), \
        'areas of the wrong length must not be applied'


def test_a_spin_up_out_of_cycles_says_so_and_is_not_called_equilibrated():
    src = _src('tests', 'run_lamata_mf6.py')
    assert 'spin_converged = True' in src
    assert 'the spin-up did NOT converge' in src
    assert "'NOT-converged spin-up'" in src
    assert "print('equilibrated heads saved:" not in src
