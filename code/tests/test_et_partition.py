# -*- coding: utf-8 -*-
"""WP2: total ET from five sources, PET spent once (cookbook 2b, 2.3-2.8).

    Ei -> Eow -> ETsoil (MMsoil) -> ETuzf (UZF, MF6) -> ETg = Eg + Tg (MMsoil)

What the soil leaves of PE and PT is UZF's demand (PETuzf), written to UZF's
PETMAX; groundwater ET then sees what remains after UZF's ACTUAL uptake --
read back from UZF's UZET, one stress period late in lagged mode -- never
the demand, which a dry deep zone cannot meet.
"""

import importlib.util
import os
import sys
import types

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for _p in (TRUNK, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T = _load('t_runmmsoil_et', os.path.join(HERE, 'test_runmmsoil.py'))
F = _load('t_soilflux_et', os.path.join(HERE, 'test_soil_flux.py'))
M = _load('t_coupler_mock_et', os.path.join(HERE, 'test_coupler_mock.py'))

IESOIL, ITSOIL, IEG, ITG, IPETUZF = 4, 5, 8, 9, 15
PT = [1.0, 2.0]                 # per vegetation type, mm/d
VEG = [40.0, 30.0]              # % of the cell
PE = 4.0                        # bare-soil potential evaporation, mm/d


def _flux(etuzf_prev=0.0, soil_frac=0.05, dgwt=100.0):
    """One loam column, one dry day (Pe = 0), water table 0.1 m below the
    soil: Eg takes all the evaporation demand left to it."""
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    cMF = T._FakeMF(nper=1, perlen=[1])
    Sm, Sfc, Sr, Ks, Tl = F._column(2)
    top = 700_000.0
    toplay, botlay = np.zeros(2), np.zeros(2)
    for l in range(2):
        toplay[l] = top if l == 0 else botlay[l - 1]
        botlay[l] = toplay[l] - Tl[l]
    Zr = [top - 600.0, top - 1500.0]
    ssoil = [soil_frac * t for t in Tl]
    return mm.flux(cMF, 1.0, 0.0, np.array(PT), PE, Zr, np.array(VEG),
                   botlay[-1] - dgwt, toplay, botlay, Tl, 2, Sm, Sfc, Sr, Ks,
                   ssoil, 0.0, dgwt, 'loam', 0, 0, 0,
                   [0.9, 0.7], [0.9, 0.7], [0.5, 0.5], [0.1, 0.1], 2,
                   np.array([2.0, 1.5]), ETUZF_prev=etuzf_prev)


# ------------------------------------------------------------- the chain
def test_uzf_demand_is_what_the_soil_left():
    """A dry soil takes nothing: UZF's demand is all of PE and PT."""
    out = _flux()
    assert np.allclose(out[IESOIL], 0.0) and np.allclose(out[ITSOIL], 0.0)
    want = PE + sum(p * v * 0.01 for p, v in zip(PT, VEG))
    assert out[IPETUZF] == pytest.approx(want)


def test_the_soil_is_served_first():
    """Moist soil: UZF gets exactly what the soil did not take."""
    dry, wet = _flux(soil_frac=0.05), _flux(soil_frac=0.35)
    et_soil = float(np.sum(wet[IESOIL]) + np.sum(wet[ITSOIL]))
    assert et_soil > 0.0
    assert wet[IPETUZF] + et_soil == pytest.approx(dry[IPETUZF], rel=1e-9)


def test_groundwater_et_sees_what_uzf_actually_took():
    base = _flux(0.0)
    f = 1.0 - 1.0 / base[IPETUZF]              # UZF took 1 mm/d
    got = _flux(1.0)
    assert got[IEG] == pytest.approx(base[IEG] * f, rel=1e-9)
    assert got[ITG] == pytest.approx(base[ITG] * f, rel=1e-6)


def test_a_dry_deep_zone_does_not_starve_groundwater_et():
    """The cookbook's case: UZF's actual << its demand. Subtracting the
    DEMAND would leave Eg nothing; the actual leaves almost all of it."""
    base = _flux(0.0)
    got = _flux(0.02 * base[IPETUZF])           # actual = 2 % of demand
    assert got[IEG] == pytest.approx(0.98 * base[IEG], rel=1e-9)
    assert got[IEG] > 0.9 * base[IEG]


def test_uzf_can_never_take_more_than_was_left():
    got = _flux(1e3)
    assert got[IEG] == pytest.approx(0.0, abs=1e-12)
    assert got[ITG] == pytest.approx(0.0, abs=1e-12)


def test_the_new_columns_are_appended_after_the_legacy_24():
    """The MODFLOW-NWT reference file keeps its 24 columns: new indices
    only ever go at the end."""
    import marmites_indices as mi
    legacy = [k for k, v in mi.INDEX_MM.items() if v < 24]
    assert len(legacy) == 24
    assert [mi.INDEX_MM[k] for k in ('iPETuzf', 'iETuzf', 'iETtot',
                                     'iRejInf')] == [24, 25, 26, 27]


# ---------------------------------------------------------- the coupler
def _run(nper=4):
    cpl, api, ctx = M._setup(nper=nper, mode='lagged')
    seen = []
    real = cpl.mm.step

    def spy(*a, **k):
        seen.append(np.array(k.get('etuzf_cell')))
        out = real(*a, **k)
        seen[-1] = (seen[-1], np.array(out['petuzf']))
        return out
    cpl.mm.step = spy
    res = cpl.run(api)
    return cpl, api, ctx, res, seen


def test_the_demand_is_written_to_petmax_after_prepare_solve():
    """PETMAX, not PET: MF6 resets PET from PETMAX every solve iteration.
    The mock's prepare_solve resets both to -1 like uzf_ad, so a demand
    written before it would not survive."""
    cpl, api, ctx, res, seen = _run()
    for n, (_prev, petuzf) in enumerate(seen):
        assert np.allclose(api.petmax_at_advance[1 + n][:ctx.ncell], petuzf)
        assert np.all(api.pet_used[1 + n][:ctx.ncell] >= 0.0), \
            'the solve saw the period input, not the demand'


def test_the_actual_uzf_et_is_read_back_and_used_one_period_late():
    cpl, api, ctx, res, seen = _run()
    want = 2.0 / cpl.area * 1000.0             # UZET = -2 m3/d per object
    for n in range(ctx.cMF.nper):
        assert np.allclose(res['etuzf'][n], want)
    assert np.allclose(seen[0][0], 0.0), 'nothing is known before SP 1'
    for n in range(1, ctx.cMF.nper):
        assert np.allclose(seen[n][0], want), 'SP %d used the wrong ETuzf' % n


def test_the_water_balance_carries_the_uzf_arm():
    cpl, api, ctx, res, seen = _run()
    ix = ctx.index
    ts = np.asarray(res['wb_ts'])
    w = cpl.area / cpl.area.sum()
    assert np.allclose(ts[:, ix['iETuzf']],
                       [np.sum(e * w) for e in res['etuzf']])
    parts = (ts[:, ix['iEi']] + ts[:, ix['iEow']] + ts[:, ix['iETsoil']]
             + ts[:, ix['iETuzf']] + ts[:, ix['iETg']])
    assert np.allclose(ts[:, ix['iETtot']], parts)
    assert len(res['pet_unmet']) == ctx.cMF.nper


# ------------------------------------------------------------ the guards
def _lst(tmp_path, uzf_disc=0.01, gwet=None):
    body = ['  UZF BUDGET FOR ENTIRE MODEL AT END OF TIME STEP    1, STRESS PERIOD   2',
            '           INFILTRATION =      100.0000          INFILTRATION =   1.0',
            ' PERCENT DISCREPANCY =  %.2f     PERCENT DISCREPANCY =  0.00' % uzf_disc,
            '  VOLUME BUDGET FOR ENTIRE MODEL AT END OF TIME STEP    1, STRESS PERIOD   2']
    if gwet is not None:
        body.append('            UZF-GWET =  %.4f           UZF-GWET =  0.1' % gwet)
    body.append(' PERCENT DISCREPANCY =  0.01     PERCENT DISCREPANCY =  0.00')
    (tmp_path / 'toy.lst').write_text('\n'.join(body) + '\n', encoding='utf-8')
    return str(tmp_path)


def test_a_uzf_budget_out_of_balance_fails_the_run(tmp_path):
    d, g = M.coup.MF6Coupler._uzf_budget_checks(_lst(tmp_path, uzf_disc=92.49))
    assert d == pytest.approx(92.49) and g == 0.0


def test_groundwater_et_in_modflow_is_caught(tmp_path):
    """2.6: a deliberately mis-set linear_gwet shows up as UZF-GWET."""
    d, g = M.coup.MF6Coupler._uzf_budget_checks(_lst(tmp_path, gwet=12.5))
    assert g == pytest.approx(12.5)
    d, g = M.coup.MF6Coupler._uzf_budget_checks(_lst(tmp_path))
    assert g == 0.0 and d == pytest.approx(0.01)


def test_check_solution_refuses_both(tmp_path):
    cpl, api, ctx = M._setup(nper=2, mode='lagged')
    cpl.run(api)
    cpl.sim_ws = _lst(tmp_path, uzf_disc=92.49, gwet=12.5)
    rep = cpl.check_solution(max_discrepancy=1.0, raise_on_fail=False)
    text = ' '.join(rep['messages'])
    assert not rep['ok']
    assert 'UZF package budget is 92.49%' in text
    assert 'UZF-GWET' in text


def test_the_builder_lists_each_columns_uzf_objects():
    src = open(os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'),
               encoding='utf-8').read()
    assert 'self.uzf_columns = [[n] + [no for no, _kk in col_children[n]]' in src
