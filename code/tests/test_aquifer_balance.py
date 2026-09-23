# -*- coding: utf-8 -*-
"""The aquifer balance line printed at the end of a run.

The run of 2026-09-23 printed "recharge to WT 2263793.8 mm/yr vs discharge
2065812.5 mm/yr". Two faults: the area was (active cells) x 50 m x 50 m,
wrong on the Voronoi mesh whose cells run from ~1 to ~4000 m2; and the rate
was a plain mean of per-step rates. Corrected, the line said 7210 mm/yr --
true, and 100 % of it one pulse on day 2, which the line now says.
"""

import importlib.util
import os
import sys
import types

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, HERE, os.path.join(CODE, 'ppMF6')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

pd = pytest.importorskip('pandas')


def _driver():
    name = '_mm_runner'
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, 'run_lamata_mf6.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


rl = _driver()


# ------------------------------------------------------------------ the area
def test_the_area_is_the_true_cells_on_the_structured_grid():
    cMF = types.SimpleNamespace(delr=np.array([10.0, 20.0]),
                                delc=np.array([5.0]))
    b = types.SimpleNamespace(nrow=1, ncol=2, surf_cells=[(0, 1, 0)])
    assert rl.active_area(cMF, b) == pytest.approx(100.0)   # 20 x 5 only


def test_the_area_is_the_true_cells_on_a_mesh():
    """A 1 m2 cell and a 16 m2 one: the old formula made both 2500 m2."""
    gp = {'ncpl': 2,
          'vertices': [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 1.0, 1.0],
                       [3, 0.0, 1.0], [4, 5.0, 0.0], [5, 5.0, 4.0],
                       [6, 1.0, 4.0]],
          'cell2d': [[0, 0.5, 0.5, 4, 0, 3, 2, 1],      # 1 x 1
                     [1, 3.0, 2.0, 4, 1, 6, 5, 4]]}     # 4 x 4
    cMF = types.SimpleNamespace(mesh_proj=object(), mesh_gridprops=gp)
    both = types.SimpleNamespace(nrow=2, ncol=1,
                                 surf_cells=[(0, 0, 0), (1, 0, 0)])
    assert rl.active_area(cMF, both) == pytest.approx(1.0 + 16.0)
    one = types.SimpleNamespace(nrow=2, ncol=1, surf_cells=[(1, 0, 1)])
    assert rl.active_area(cMF, one) == pytest.approx(16.0)


def test_an_area_array_of_the_wrong_size_is_refused():
    cMF = types.SimpleNamespace(delr=np.ones(3), delc=np.ones(2))
    b = types.SimpleNamespace(nrow=4, ncol=1, surf_cells=[])
    with pytest.raises(ValueError):
        rl.active_area(cMF, b)


# --------------------------------------------------------------- the balance
def _cum(times, **rates):
    """Cumulative volumes from per-period rates [m3/d] held over each period."""
    t = np.asarray(times, float)
    dt = np.diff(np.concatenate([[0.0], t]))
    return pd.DataFrame({k: np.cumsum(np.asarray(v, float) * dt)
                         for k, v in rates.items()}), t


def test_rates_are_volumes_over_time_not_a_mean_of_rates():
    """A 1-day step at 1000 m3/d and a 29-day one at 0: 1000 m3 over 30 d
    after the first period -- not the 500 m3/d a plain mean gives."""
    cum, t = _cum([1.0, 2.0, 31.0], **{'UZF-GWRCH_IN': [7.0, 1000.0, 0.0],
                                        'DRN_OUT': [0.0, 0.0, 0.0]})
    bal = rl.aquifer_balance(cum, t, area=1000.0)
    assert bal['days'] == pytest.approx(30.0)
    assert bal['recharge'] == pytest.approx(1000.0 / 30 / 1000 * 1000 * 365)
    assert bal['peak_share'] == pytest.approx(1.0)
    assert bal['peak_time'] == pytest.approx(2.0)
    assert bal['peak_mm'] == pytest.approx(1000.0)     # 1000 m3 on 1000 m2


def test_every_drain_well_and_ghb_outflow_is_discharge_nothing_else():
    cum, t = _cum([1.0, 11.0], **{
        'UZF-GWRCH_IN': [0.0, 10.0], 'DRN_OUT': [0.0, 1.0],
        'DRN2_OUT': [0.0, 2.0], 'WEL_OUT': [0.0, 3.0], 'GHB_OUT': [0.0, 4.0],
        'STO-SY_OUT': [0.0, 50.0], 'UZF-GWRCH_OUT': [0.0, 60.0],
        'DRN_IN': [0.0, 70.0]})
    bal = rl.aquifer_balance(cum, t, area=365.0)
    assert sorted(bal['discharge_terms']) == ['DRN2_OUT', 'DRN_OUT',
                                              'GHB_OUT', 'WEL_OUT']
    assert bal['discharge'] == pytest.approx(10.0 * 1000.0)   # 1+2+3+4 m3/d


def test_the_lines_warn_of_a_pulse_and_only_then():
    base = {'days': 60.0, 'area_km2': 4.8, 'recharge': 7209.8,
            'discharge': 6579.3, 'discharge_terms': {'DRN2_OUT': 6501.7},
            'peak_time': 2.0, 'peak_mm': 1185.0}
    pulse = rl.balance_lines(dict(base, peak_share=1.0))
    assert '7209.8' in pulse[0] and '4.800 km2' in pulse[0]
    assert any('ONE step' in l and 't = 2 d' in l for l in pulse)
    steady = rl.balance_lines(dict(base, peak_share=0.05))
    assert not any('ONE step' in l for l in steady)


def test_the_run_uses_the_new_balance_and_says_when_it_cannot():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'),
               encoding='utf-8').read()
    body = src[src.index('def main('):]
    assert 'aquifer_balance(' in body and 'active_area(cMF, b)' in body
    assert "b.idomain[0] > 0)) * float(np.mean(cMF.delr))" not in src
    assert 'aquifer balance: not computed' in body


LAMATA_LST = r'E:\00code_ws\LaMata_MM-MF6\MF6_ws_voronoi\lamata.lst'


@pytest.mark.skipif(not os.path.exists(LAMATA_LST),
                    reason='no La Mata run in the workspace')
def test_a_real_list_file_reads():
    """Whatever run is there: the columns the balance needs are found."""
    flopy = pytest.importorskip('flopy')
    lst = flopy.utils.Mf6ListBudget(LAMATA_LST)
    bal = rl.aquifer_balance(lst.get_dataframes(diff=False)[1],
                             lst.get_times(), area=4.8e6)
    assert bal['days'] > 0 and bal['discharge_terms']
    assert np.isfinite(bal['recharge']) and np.isfinite(bal['discharge'])
