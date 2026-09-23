# -*- coding: utf-8 -*-
"""The Plots panel's settings reach the figures.

The audit behind the cookbook's Appendix B found eight of them read by
nothing. The hydrological year was asked of cMF by the figures -- and never
set, so every one fell back to October. The Sankey, the per-point series,
the result maps and the input maps were drawn whatever the switches said
(the Sankey's was a literal True). The water-balance unit had no figure left
to govern: the Sankey was hard-wired to mm per year.
"""

import os
import sys
import types

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6'),
           os.path.join(CODE, 'MARMITESutilities'),
           os.path.join(CODE, 'MARMITESutilities', 'MARMITESplot')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

pp = pytest.importorskip('marmites_postprocess')
import marmites_config as mcfg                                 # noqa: E402
import marmites_props as props                                 # noqa: E402


def test_the_figure_settings_are_set_where_the_figures_read_them():
    cfg = mcfg.RunConfig.from_dict({'postproc': {
        'hydro_year_start': 9, 'wb_unit': 'day',
        'tick_trimester_years': 3, 'tick_semester_years': 7}})
    cMF = types.SimpleNamespace()
    props.apply_plot_settings(cfg, cMF)
    assert cMF.iniMonthHydroYear == 9, 'the hydrological year stays October'
    assert cMF.plt_WB_unit == 'day'
    assert pp._tick_years(cMF) == {'maxYearsTickTrimester': 3,
                                   'maxYearsTickSemester': 7}


def test_the_tick_defaults_are_the_ones_the_routines_had():
    assert pp._tick_years(types.SimpleNamespace()) == {
        'maxYearsTickTrimester': 5, 'maxYearsTickSemester': 10}


@pytest.fixture
def drawn(monkeypatch):
    """native_suite with every figure family replaced by a recorder."""
    calls = []
    monkeypatch.setattr(pp, '_mmplot', lambda trunk=None: object())
    for fn in ('_native_obs_timeseries', '_native_result_maps',
               '_native_sankey', '_native_sankey_obs'):
        monkeypatch.setattr(pp, fn,
                            lambda *a, _n=fn, **k: calls.append(_n) or [])
    monkeypatch.setattr(pp, '_aquifer_pass', lambda *a, **k: None)
    res = {'perc': np.zeros((2, 1))}
    cMF = types.SimpleNamespace(modelname='m')
    ctx = types.SimpleNamespace(cells=[(0, 0, 0)])

    def run(**kw):
        del calls[:]
        pp.native_suite('out', cMF, ctx, res, sim_ws='ws', verbose=False,
                        **kw)
        return set(calls)
    return run


def test_every_switch_on_draws_every_family(drawn):
    got = drawn(sankey=True, obs_series=True, result_maps=True)
    assert got == {'_native_obs_timeseries', '_native_result_maps',
                   '_native_sankey', '_native_sankey_obs'}


@pytest.mark.parametrize('off, gone', [
    ({'sankey': False}, {'_native_sankey', '_native_sankey_obs'}),
    ({'obs_series': False}, {'_native_obs_timeseries'}),
    ({'result_maps': False}, {'_native_result_maps'}),
])
def test_a_switch_off_draws_nothing_of_its_family(drawn, off, gone):
    kw = dict(sankey=True, obs_series=True, result_maps=True)
    kw.update(off)
    got = drawn(**kw)
    assert not (got & gone), '%s drawn although switched off' % (got & gone)


def test_the_sankey_can_be_drawn_per_day():
    """It was hard-wired to per-year; the unit is postproc.wb_unit now, and
    the threshold stays a per-year quantity, converted."""
    src = open(os.path.join(CODE, 'MARMITESutilities', 'MARMITESplot',
                            'MARMITESplot_v3.py'), encoding='utf-8').read()
    body = src[src.index('def plotWBsankey'):]
    body = body[:body.index('\ndef ')]
    assert 'per_day=False' in body
    assert 'treshold = treshold / 365.0' in body
    assert "(1.0 if per_day else 365.0)" in body
    assert "'d' if per_day else 'y'" in body
    post = open(os.path.join(CODE, 'ppMF6', 'marmites_postprocess.py'),
                encoding='utf-8').read()
    assert "getattr(smf, 'plt_WB_unit', 'year'))" in post


def test_the_run_passes_every_switch():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'),
               encoding='utf-8').read()
    for needle in ('props.apply_plot_settings(cfg, cMF)',
                   'cfg.postproc.input_maps', 'cfg.postproc.sankey',
                   'cfg.postproc.obs_series', 'cfg.postproc.result_maps'):
        assert needle in src, '%s is not passed' % needle
    assert 'sankey=True,' not in src, 'the Sankey is a literal True again'
