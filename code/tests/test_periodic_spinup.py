# -*- coding: utf-8 -*-
"""The periodic spin-up (2026-09-25).

Each spin-up cycle used to start with a STEADY period driven by the mean
forcing. MF6 ignores the initial heads in a steady period, so the previous
cycle's last heads never reached the next one: every cycle restarted from
the same mean-forcing equilibrium, ended 0.6 m below it, and the cycles
agreed with each other rather than with themselves. Now every cycle after
the first starts where the last one ended -- heads, the water the
unsaturated zone holds, the soil, and the previous day's exchange terms --
with no steady period.
"""

import importlib.util
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(CODE, '..', 'example', 'LaMata'))
for _p in (CODE, HERE, os.path.join(CODE, 'ppMF6'),
           os.path.join(CODE, 'MARMITESutilities'),
           os.path.join(CODE, 'MARMITESsoil'), os.path.join(CODE, 'ppMF_FloPy')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load('t_coupler_mock_spin', os.path.join(HERE, 'test_coupler_mock.py'))


# ---------------------------------------------------------- the coupler
def _periodic(nper=3, carry=None):
    cpl, api, ctx = M._setup(nper=nper)
    cpl.mf6b.steady_first = False
    cpl.carry_in = carry
    seen = []
    real = cpl.mm.step

    def spy(*a, **k):
        seen.append((np.array(a[4]), np.array(k['rejinf_cell']),
                     np.array(k['etuzf_cell'])))
        return real(*a, **k)
    cpl.mm.step = spy
    cpl.run(api)
    return cpl, api, ctx, seen


def test_without_a_steady_period_every_step_is_the_run():
    cpl, api, ctx, _seen = _periodic()
    assert len(api.finf_at_advance) == ctx.cMF.nper, 'a steady step was run'
    cpl2, api2, ctx2 = M._setup(nper=3)
    cpl2.run(api2)
    assert len(api2.finf_at_advance) == ctx2.cMF.nper + 1


def test_the_first_day_reads_what_the_last_cycle_left():
    n = M._setup(nper=1)[2].ncell
    carry = {'exf': np.full(n, 7.0), 'rej': np.full(n, 3.0),
             'etuzf': np.full(n, 0.5)}
    _cpl, _api, _ctx, seen = _periodic(carry=carry)
    exf0, rej0, et0 = seen[0]
    assert np.allclose(exf0, 7.0) and np.allclose(rej0, 3.0)
    assert np.allclose(et0, 0.5)
    exf1, _r, _e = seen[1]
    assert not np.allclose(exf1, 7.0), 'only the first day is carried'


def test_a_cycle_hands_its_last_day_on():
    cpl, _api, ctx, _seen = _periodic()
    assert set(cpl.carry_out) == {'exf', 'rej', 'etuzf'}
    assert cpl.carry_out['exf'].shape == (ctx.ncell,)
    assert np.allclose(cpl.carry_out['etuzf'], cpl.etuzf_prev)


# ------------------------------------------------------------- the soil
def test_a_carried_soil_starts_from_its_own_state():
    """The first period ignored the state and read the panel's initial
    moisture -- a carried state has to be what it starts from."""
    cpl, api, ctx = M._setup(nper=1)
    mm = cpl.mm
    heads = np.full(ctx.ncell, 699.0)
    zero = np.zeros(ctx.ncell)
    fresh = mm.init_state(ctx)
    a = mm.step(ctx, 0, 0, heads, zero, fresh)
    wet = mm.init_state(ctx)
    wet.Ssoil_ini[:] = 1e3                 # mm, far wetter than the panel's
    b = mm.step(ctx, 0, 0, heads, zero, wet)
    assert np.allclose(a['MM'], b['MM']), 'an uncarried state is ignored'
    wet2 = mm.init_state(ctx)
    wet2.Ssoil_ini[:] = 1e3
    wet2.carried = True
    c = mm.step(ctx, 0, 0, heads, zero, wet2)
    assert not np.allclose(a['MM'], c['MM'])
    assert mm.init_state(ctx).carried is False


# ------------------------------------------------------------ the build
@pytest.fixture(scope='module')
def cmf():
    pytest.importorskip('flopy')
    import lamata_model
    # La Mata's model description as the run builds it -- no
    # parameter file (lamata_model derives the outcrop layer too)
    c = lamata_model.lamata_cmf()
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]
    return c


def _build(cmf, tmp, **kw):
    mf6mod = _load('marmites_mf6_spin', os.path.join(CODE, 'ppMF6',
                                                     'marmites_mf6.py'))
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp))
    for k, v in kw.items():
        setattr(b, k, v)
    b.build()
    b.write()
    return b


def test_a_periodic_cycle_is_transient_from_the_start(cmf, tmp_path):
    pp = _load('marmites_postprocess_spin', os.path.join(
        CODE, 'ppMF6', 'marmites_postprocess.py'))
    b = _build(cmf, tmp_path / 'p', steady_first=False)
    assert b.nper == 3
    assert pp.steady_first(str(tmp_path / 'p')) is False
    b = _build(cmf, tmp_path / 's')
    assert b.nper == 4
    assert pp.steady_first(str(tmp_path / 's')) is True


def test_the_unsaturated_zone_keeps_its_water(cmf, tmp_path):
    """MF6's WCNEW per object -> the next cycle's thti, clipped to [thtr,
    thts]; an object under the water table keeps the panel's."""
    b0 = _build(cmf, tmp_path / 'a')
    wc = np.zeros(b0.nuzfcells)
    wc[0], wc[1] = 0.4455, 0.9            # 0.9 is above thts: clipped
    carry = b0.thti_from_wc(wc)
    k, i, j = b0.uzf_obj_kij[0]
    assert carry[k, i, j] == 0.4455
    assert np.isnan(carry[b0.uzf_obj_kij[2]]), '0 = under the water table'
    b1 = _build(cmf, tmp_path / 'b', uzf_thti_carry=carry)
    thti = {int(r[0]): float(r[8]) for r in b1.uzf_packagedata}
    thts = {int(r[0]): float(r[7]) for r in b1.uzf_packagedata}
    assert thti[0] == pytest.approx(0.4455)
    assert thti[1] == pytest.approx(thts[1])
    assert thti[2] == {int(r[0]): float(r[8]) for r in b0.uzf_packagedata}[2]


# ---------------------------------------------------- the driver's pieces
def test_the_aquifer_balance_counts_from_zero_without_a_steady_period():
    pd = pytest.importorskip('pandas')
    rl = _load('_mm_runner_spin', os.path.join(HERE, 'run_lamata_mf6.py'))
    cum = pd.DataFrame({'UZF-GWRCH_IN': [10.0, 20.0, 30.0],
                        'DRN2_OUT': [1.0, 2.0, 3.0]})
    t = [1.0, 2.0, 3.0]
    s = rl.aquifer_balance(cum, t, area=365.0)
    p = rl.aquifer_balance(cum, t, area=365.0, steady_first=False)
    assert s['days'] == 2.0 and p['days'] == 3.0
    assert p['recharge'] == pytest.approx(30.0 / 3.0 / 365.0 * 1000 * 365)


def test_the_spin_up_loop_carries_all_three_states():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    loop = src[src.index('for cyc in range(ncyc):'):]
    loop = loop[:loop.index('MMsoil.report_tallies()')]
    for piece in ('b.strt_array = prev_heads', 'b.steady_first = False',
                  'b.thti_from_wc(cpl.uzf_wc_final)', 'st.carried = True',
                  'cpl.carry_in = _carry'):
        assert piece in loop, piece


# ------------------------------------------- a run started from saved state
def _rl():
    return _load('_mm_runner_saved', os.path.join(HERE, 'run_lamata_mf6.py'))


def _fake(nlay=2, nrow=3, ncol=1, ncell=3, nsl=2):
    b = SimpleNamespace(nlay=nlay, nrow=nrow, ncol=ncol,
                        thti_from_wc=lambda wc: np.full((nlay, nrow, ncol), 0.2))
    cpl = SimpleNamespace(uzf_wc_final=np.zeros(4),
                          carry_out={'exf': np.full(ncell, 1.5),
                                     'rej': np.full(ncell, 0.5),
                                     'etuzf': np.full(ncell, 0.1)})
    st = SimpleNamespace(Ssoil_ini=np.arange(ncell * nsl, dtype=float)
                         .reshape(ncell, nsl))
    ctx = SimpleNamespace(ncell=ncell, _nslmax=nsl)
    return b, cpl, st, ctx


def test_the_saved_state_round_trips(tmp_path):
    rl = _rl()
    b, cpl, st, ctx = _fake()
    fn = rl.save_run_state(str(tmp_path / 'hi'), b, cpl, st)
    assert fn.endswith('hi_state.npz')
    got = rl.load_run_state(str(tmp_path / 'hi'), b, ctx)
    assert np.allclose(got['uzf_wc'], 0.2)
    assert np.allclose(got['soil'], st.Ssoil_ini)
    assert np.allclose(got['carry']['exf'], 1.5)
    assert np.allclose(got['carry']['etuzf'], 0.1)


def test_a_state_that_does_not_fit_is_not_used(tmp_path, capsys):
    rl = _rl()
    b, cpl, st, ctx = _fake()
    rl.save_run_state(str(tmp_path / 'hi'), b, cpl, st)
    other = SimpleNamespace(ncell=5, _nslmax=2)
    assert rl.load_run_state(str(tmp_path / 'hi'), b, other) is None
    assert 'does not fit' in capsys.readouterr().out
    assert rl.load_run_state(str(tmp_path / 'none'), b, ctx) is None
    assert 'panel\'s initial values' in capsys.readouterr().out


def test_a_run_from_saved_heads_has_no_steady_period():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    block = src[src.index("if _kind == 'saved':"):src.index('b.build()\n    b.write()')]
    for piece in ('b.steady_first = False', 'load_run_state(pref, b, ctx)',
                  "b.uzf_thti_carry = saved_state['uzf_wc']"):
        assert piece in block, piece
    loop = src[src.index('for cyc in range(ncyc):'):]
    loop = loop[:loop.index('MMsoil.report_tallies()')]
    assert "st.Ssoil_ini[:] = saved_state['soil']" in loop
    assert "saved_state['carry']" in loop
    assert 'save_run_state(pref, b, cpl, st)' in src
    assert '(skips the spin-up)' not in src


# ------------------------- a cycle that passed survives a later failure
def test_a_passed_cycle_is_saved_whole_and_dropped_at_the_end(tmp_path):
    """2026-10-04: cycle 2 failed one sub-step of SP74 and the converged
    cycle 1 was lost with it. Each passed cycle is now on disk at once as
    <save>_lastcycle -- a full state the panel lists -- and removed when the
    spin-up's own state is saved."""
    import types
    import marmites_config as mcfg
    from marmites_mf6 import clsMF6
    rl = _rl()
    b, cpl, st, ctx = _fake()
    cMF = SimpleNamespace(xllcorner=0.0, yllcorner=0.0, delr=[50.0],
                          nrow=b.nrow, ncol=b.ncol)
    b.cMF, b.idomain = cMF, np.ones((b.nlay, b.nrow, b.ncol), int)
    b.save_heads_asc = types.MethodType(clsMF6.save_heads_asc, b)
    ctx.cells = [(k, k, 0, k) for k in range(ctx.ncell)]
    a = SimpleNamespace(state_dir=str(tmp_path),
                        config=mcfg.RunConfig.from_dict({}))
    res = {'perc': np.full((4, ctx.ncell), 1e-4),
           'etg': np.full((4, ctx.ncell), 2e-5)}
    heads = np.full((b.nlay, b.nrow, b.ncol), 700.0)
    name = rl.save_cycle_state(a, b, cMF, ctx, 'hi_x', heads, cpl, st, res,
                               grid={'shape': [3, 1], 'signature': 'sig'})
    assert name == 'hi_x_lastcycle'
    got = rl.load_run_state(str(tmp_path / name), b, ctx)
    assert np.allclose(got['soil'], st.Ssoil_ini)
    listed = mcfg.saved_states(str(tmp_path), b.nlay)
    assert [s['name'] for s in listed] == [name] and listed[0]['full']
    assert [s['name'] for s in mcfg.saved_states(
        str(tmp_path), b.nlay, what='spinup.steady_means')] == [name]
    rl.drop_cycle_state(a, 'hi_x', b.nlay)
    assert os.listdir(str(tmp_path)) == []


def test_the_loop_keeps_the_last_passed_cycle():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    loop = src[src.index('for cyc in range(ncyc):'):src.index('check = None')]
    # a failure -- or an interruption -- names the cycle that is on disk
    body = loop[loop.index('try:'):loop.index('prev_heads = _final_heads()')]
    assert 'res = cpl.run(api)' in body and 'cpl.check_solution(' in body
    assert 'except BaseException:' in body and 'raise' in body
    # saved as it passes, but not the last cycle: the final save writes it
    assert 'if cycle_base and cyc + 1 < ncyc:' in loop
    assert 'save_cycle_state(' in loop
    tail = src[src.index('check = None'):src.index('def _run_postproc')]
    assert tail.index('drop_cycle_state(') > tail.index("mp + '_etg.asc'")
