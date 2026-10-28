# -*- coding: utf-8 -*-
"""The soil and aquifer maps also on a PERCENTILE colour ramp (user,
2026-10-28: "To avoid to have a map that shows only outlier, please produce
also the MMmaps and GWmaps using a percentil color ramp").

On the linear ramp the springs' Rn of -107 mm/d took the whole range and
the rest of the catchment came out one colour. Every time-mean MMmap and
GWmap now has a twin, <map>_pct, whose classes each hold an equal share of
the cells; zero is a class of its own, and a field of both signs keeps zero
at the neutral middle of its diverging ramp.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in ('', 'MARMITESutilities', 'MARMITESutilities/MARMITESplot'):
    sys.path.insert(0, os.path.join(TRUNK, p))

import matplotlib  # noqa: E402
matplotlib.use('agg')
import matplotlib.pyplot as plt  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P = _load('mmplot_v3_pct', os.path.join(TRUNK, 'MARMITESutilities',
                                        'MARMITESplot', 'MARMITESplot_v3.py'))
R = _load('_app_results_pct', os.path.join(TRUNK, 'app', 'lib', 'results.py'))
BLUES = matplotlib.colormaps['Blues']


def _share(v, bounds):
    """How many values fall in each class [b_i, b_i+1)."""
    return np.histogram(v, bins=bounds)[0]


# ---- the classes ----------------------------------------------------------

def test_an_outlier_takes_one_class_not_the_ramp():
    v = np.r_[np.linspace(0.1, 1.0, 99), 100.0]
    b, c = P.percentile_classes(v, BLUES)
    assert len(b) == 11 and len(c) == 10
    assert np.all(np.diff(b) > 0) and b[0] == v.min() and b[-1] > v.max()
    assert _share(v, b).tolist() == [10] * 10, 'ten cells in each class'
    # 99 of the 100 cells in the first 9 classes, below 1.0; the outlier
    # alone at the top
    assert b[-2] < 1.0


def test_zero_is_a_class_and_the_percentiles_are_over_the_others():
    v = np.r_[np.zeros(80), np.arange(1.0, 21.0)]   # a stream map: 0 away
    b, c = P.percentile_classes(v, BLUES)
    assert b[0] == 0.0 and b[1] == 1.0, 'the zero class: [0, smallest > 0)'
    s = _share(v, b)
    assert s[0] == 80 and s[1:].tolist() == [2] * 10
    assert c[0] == BLUES(0.0), 'zero at the foot of the ramp'


def test_a_rounding_zero_is_zero():
    v = np.r_[np.full(10, -1e-12), np.arange(1.0, 11.0)]
    b, c = P.percentile_classes(v, BLUES)
    assert b[0] == 0.0, 'one-signed: no diverging ramp for a -1e-12'


def test_both_signs_keep_zero_at_the_neutral_middle():
    rng = np.random.default_rng(1)
    v = np.r_[-rng.random(30) * 100.0, np.zeros(5), rng.random(70) * 10.0]
    b, c = P.percentile_classes(v, BLUES)
    ramp = matplotlib.colormaps['coolwarm_r']
    z = int(np.flatnonzero(b == 0.0)[0])
    assert c[z] == ramp(0.5), 'the zero class neutral'
    assert len(b[:z]) == 5 and len(c) - z - 1 == 5, 'five classes a side'
    # the negative side red (coolwarm_r's low end), the positive blue
    assert all(cc[0] > cc[2] for cc in c[:z])
    assert all(cc[2] > cc[0] for cc in c[z + 1:])
    s = _share(v, b)
    assert s[:z].tolist() == [6] * 5 and s[z] == 5
    assert s[z + 1:].tolist() == [14] * 5


def test_a_field_of_one_value_has_no_classes():
    assert P.percentile_classes(np.full(7, 3.0), BLUES) is None


def test_the_bar_labels():
    assert [P._pct_label(x) for x in (0.0, 733.25, 13.64, -106.66, 0.01234,
                                      2.5e-9)] == \
        ['0', '733', '13.6', '-107', '0.0123', '2.50e-09']


# ---- plotLAYER ------------------------------------------------------------

def test_each_layer_has_its_own_classes(tmp_path, monkeypatch):
    """The legacy percentile mode took every layer's classes from the
    first layer's values."""
    nrow, ncol = 6, 8
    rng = np.random.default_rng(0)
    V = np.stack([rng.random((nrow, ncol)), 50.0 + rng.random((nrow, ncol))])
    V = V[None]
    mask = np.zeros((2, nrow, ncol), bool)
    seen = []
    real = plt.savefig

    def grab(*a, **k):
        f = plt.gcf()
        bars = [ax for ax in f.axes if hasattr(ax, '_colorbar')]
        seen.append([ax.collections[0].norm.boundaries for ax in f.axes
                     if ax not in bars])
        seen.append(bars)
        return real(*a, **k)
    monkeypatch.setattr(plt, 'savefig', grab)
    P.plotLAYER(days=[0], str_per=[0], Date='NA', JD='NA', ncol=ncol,
                nrow=nrow, nlay=2, nplot=2, V=V, cmap=BLUES, CBlabel='x',
                msg='', plt_title='T_pct', MM_ws=str(tmp_path),
                interval_type='percentile', Vmax=[V.max()], Vmin=[V.min()],
                fmt='%5.2f', mask=mask, hnoflo=-999.9)
    (b1, b2), bars = seen
    assert b1[0] == V[0, 0].min() and b1[-1] > V[0, 0].max() and b1[-1] < 1.01
    assert b2[0] == V[0, 1].min() and b2[-1] > V[0, 1].max()
    assert len(bars) == 2, 'one colour bar per layer'
    assert os.path.exists(str(tmp_path / '_sp_plt_T_pct_1_1.png'))


# ---- the post-processing draws the twins ----------------------------------

def test_every_mean_map_has_its_percentile_twin():
    src = open(os.path.join(TRUNK, 'ppMF6', 'marmites_postprocess.py'),
               encoding='utf-8').read()
    body = src[src.index('def _native_result_maps('):]
    body = body[:body.index('\ndef ')]
    assert ("twin = (days is None and prefix in ('MMmap', 'GWmap')" in body)
    assert "(('percentile', '_pct'),)" in body
    assert "plt_title='%s_%s%s' % (prefix, stem, suffix)" in body


# ---- the Results panel ----------------------------------------------------

def test_the_twin_on_the_results_panel():
    f = '_sp_plt_GWmap_Rn_pct_1_1.png'
    assert R.classify('_output', f) == 'output'
    assert R.title('_output', f) == ('Aquifer: net recharge (Rg + Rsw - Exf '
                                     '- ETg) (percentile colour ramp)')
    assert R.title('_output', '_sp_plt_MMmap_P_pct_1_1.png') == \
        'Soil column: rainfall (percentile colour ramp)'
    names = ['_sp_plt_MMmap_P_1_1.png', '_sp_plt_MMmap_P_pct_1_1.png',
             '_sp_plt_MMmap_Pe_1_1.png', '_sp_plt_MMmap_Pe_pct_1_1.png',
             '_sp_plt_GWmap_head_1_1.png', '_sp_plt_GWmap_head_pct_1_1.png',
             '_sp_plt_GWmap_head_series_day00001_1_1.png',
             '_sp_plt_GWmap_Rg_1_1.png', '_sp_plt_GWmap_Rg_pct_1_1.png']
    got = [it[1] for it in R.arrange([('_output', n, n)
                                      for n in names[::-1]])['output']]
    assert got == names, got


def test_a_map_and_its_twin_share_a_row():
    names = ['_sp_plt_MMmap_A_1_1.png',                  # no twin
             '_sp_plt_MMmap_B_1_1.png', '_sp_plt_MMmap_B_pct_1_1.png',
             '_sp_plt_MMmap_C_1_1.png', '_sp_plt_MMmap_C_pct_1_1.png']
    items = [('_output', n, n, n) for n in names]
    rows = [[it[1] for it in r] for r in R.rows(items)]
    assert rows == [names[:1], names[1:3], names[3:]]
