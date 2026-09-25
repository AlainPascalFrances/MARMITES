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
    cpl, api, ctx = M._setup(nper=nper, mode='lagged')
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
    cpl2, api2, ctx2 = M._setup(nper=3, mode='lagged')
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
    cpl, api, ctx = M._setup(nper=1, mode='lagged')
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
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s1L.ini')):
        pytest.skip('La Mata dataset not present')
    pytest.importorskip('flopy')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=0), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s1L.ini',
                   xllcorner=739300.0, yllcorner=4553050.0)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
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
