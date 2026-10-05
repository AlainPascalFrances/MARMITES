# -*- coding: utf-8 -*-
"""Binary outputs read PER STRESS PERIOD (cookbook WP6.4, the CdL trap).

MF6 saves a head and a budget record at every TIME STEP; adaptive time
stepping cuts a period into as many as it needs (1354 records for 365
periods on 2026-10-05). The post-processing read "the last nper records" as
the nper periods -- the raw heads at the observation points, the aquifer
arms of the Sankeys, the per-layer aquifer series and maps were all shifted
in time. Now a state is the period's LAST record and a rate its steps'
time-weighted mean.
"""
import os
import shutil
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for p in (CODE, os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'MARMITESutilities'),
          os.path.join(CODE, 'MARMITESutilities', 'MARMITESplot')):
    if p not in sys.path:
        sys.path.insert(0, p)

PP = pytest.importorskip('marmites_postprocess')


def test_records_are_grouped_into_their_periods():
    # a steady period, then two periods of 1 d in 3 and 2 unequal steps
    times = [1.0, 1.2, 1.5, 2.0, 2.25, 3.0]
    kk = [(0, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2)]
    g = PP.period_steps(times, kk, steady=True)
    assert [[r for r, _w in grp] for grp in g] == [[1, 2, 3], [4, 5]]
    assert [w for _r, w in g[0]] == pytest.approx([0.2, 0.3, 0.5])
    assert [w for _r, w in g[1]] == pytest.approx([0.25, 0.75])
    # without a steady period the first one is a period like the others
    assert len(PP.period_steps(times, kk, steady=False)) == 3

    class _H:
        def get_times(self):
            return times

        def get_kstpkper(self):
            return kk
    assert PP.period_end_kk(_H(), steady=True) == [(2, 1), (1, 2)]


def _mf6_exe():
    exe = shutil.which('mf6')
    if exe:
        return exe
    try:
        import mm_paths
        cand = os.path.join(os.path.dirname(mm_paths.LIBMF6), 'mf6.exe')
        return cand if os.path.exists(cand) else None
    except Exception:                                  # noqa: BLE001
        return None


@pytest.fixture(scope='module')
def stepped_run(tmp_path_factory):
    """A steady period, then 3 periods of 4 time steps each, pumping more
    every period: heads and storage change step by step."""
    flopy = pytest.importorskip('flopy')
    exe = _mf6_exe()
    if exe is None:
        pytest.skip('no mf6 executable')
    ws = str(tmp_path_factory.mktemp('stepped'))
    name = 'st'
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=4, perioddata=[(1.0, 1, 1.0)]
                          + [(1.0, 4, 1.5)] * 3, time_units='DAYS')
    flopy.mf6.ModflowIms(sim, complexity='SIMPLE')
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, delr=50.0, delc=50.0,
                            top=100.0, botm=0.0)
    flopy.mf6.ModflowGwfic(gwf, strt=90.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0, save_flows=True)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=1e-4, sy=0.1,
                            steady_state={0: True}, transient={1: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 90.0]])
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={
        0: [[(0, 2, 2), 0.0]], 1: [[(0, 2, 2), -20.0]],
        2: [[(0, 2, 2), -60.0]], 3: [[(0, 2, 2), -5.0]]}, save_flows=True)
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord=f'{name}.hds',
                           budget_filerecord=f'{name}.cbc',
                           saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    if not ok:
        pytest.skip('the stepped mf6 run failed here')
    return ws, name


def test_heads_are_the_end_of_each_period(stepped_run):
    import flopy
    ws, name = stepped_run
    hds = flopy.utils.HeadFile(os.path.join(ws, '%s.hds' % name))
    assert len(hds.get_kstpkper()) == 13                  # 1 + 3 x 4 steps
    assert PP.period_end_kk(hds, steady=True) == [(3, 1), (3, 2), (3, 3)]
    got = PP._obs_head_series(ws, name, 1, [(2, 2)], 3)[0, 0]
    want = [hds.get_data(kstpkper=(3, p))[0, 2, 2] for p in (1, 2, 3)]
    assert got == pytest.approx(want)
    # the old reading (the last 3 saved steps) is not that
    old = [hds.get_data(kstpkper=k)[0, 2, 2] for k in hds.get_kstpkper()[-3:]]
    assert not np.allclose(got, old)


def test_rates_are_each_periods_time_weighted_mean(stepped_run):
    import flopy
    ws, name = stepped_run
    cbc = flopy.utils.CellBudgetFile(os.path.join(ws, '%s.cbc' % name))
    cMF = SimpleNamespace(nlay=1, nrow=3, ncol=3)
    agg = PP._aquifer_pass(ws, name, cMF, [[(2, 2)]], 3, cache_dir=None,
                           verbose=False)
    t = np.asarray(cbc.get_times())
    dt = np.diff(np.r_[0.0, t])
    kk = cbc.get_kstpkper()
    for p, per in enumerate((1, 2, 3)):
        rows = [i for i, k in enumerate(kk) if k[1] == per]
        sto = [np.asarray(cbc.get_data(kstpkper=kk[i], text='STO-SS',
                                       full3D=True)[0]).reshape(1, -1)[0, 8]
               for i in rows]
        w = dt[rows] / dt[rows].sum()
        assert agg['STO-SS'][p, 0, 0] == pytest.approx(float(np.dot(w, sto)))
        # the well is constant within a period
        assert agg['WEL'][p, 0, 0] == pytest.approx((-20.0, -60.0, -5.0)[p])


def test_package_means_weight_each_step_by_its_length(stepped_run):
    ws, name = stepped_run
    b = PP.package_budget(ws, '%s.cbc' % name, max_samples=None)
    # equal 1-day periods pumping 20, 60, 5: the time mean is 85 / 3
    assert b['WEL'] == pytest.approx(-85.0 / 3.0)
