# -*- coding: utf-8 -*-
"""Sprinkler irrigation down the whole soil column (soil.irr_infiltration).

2026-10-04, the CRR run: on 11 Aug 2008 La Mata's field got 25 mm of
irrigation in one daily step. Eq. 1b lets a day's water into the TOP horizon
only (~19 mm of room there), so the rest ran off, and the cascade carried it
into the streams as a pulse that a dry summer never sees. Sprinklers run for
hours at a rate the soil absorbs and the water soaks down the column. With
'column', what of the IRRIGATION the top horizon cannot take fills the
horizons below, top-down, each up to saturation, before any runs off; rain
keeps Eq. 1b.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


T = _load('t_runmmsoil_irrcol', os.path.join(HERE, 'test_runmmsoil.py'))
IX, IXS = T.INDEX_MM, T.INDEX_MM_S

# flux() return tuple
IRO, IRP, ISSOIL, II, IREINF = 2, 3, 6, 14, 16

SM = np.array([0.41, 0.40, 0.38])
TL = np.array([200.0, 300.0, 500.0])          # mm


def _flux(Pe, irr, frac, runon=0.0, nsl=3):
    """One day, no ET demand, no groundwater: only the surface split."""
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    cMF = T._FakeMF(nper=1, perlen=[1])
    Sm, Tl = SM[:nsl], TL[:nsl]
    top = 700_000.0
    tops = top - np.concatenate([[0.0], np.cumsum(Tl)[:-1]])
    bots = tops - Tl
    return mm.flux(cMF, 1.0, Pe, np.zeros(1), 0.0, [top - 50.0], np.zeros(1),
                   top - 50_000.0, tops, bots, Tl, nsl, Sm, Sm * 0.75,
                   Sm * 0.2, np.full(nsl, 1e-9), np.asarray(frac)[:nsl] * Tl,
                   0.0, 50_000.0, 'loam', 0, 0, 0, [0.0], [0.0], [0.5], [0.1],
                   1, np.zeros(1), RUNON=runon, IRR=irr)


def test_top_is_eq_1b_and_column_fills_the_horizons_below():
    frac = [0.38, 0.39, 0.30]              # room: 6, 3, 40 mm
    top = _flux(30.0, 0.0, frac)
    col = _flux(30.0, 30.0, frac)
    assert top[IRO] == pytest.approx(24.0)                 # 30 - 6
    assert top[II] == pytest.approx(6.0)
    # top-down: 6 into the top horizon, 3 fills the second, 21 the third
    assert col[IRO] == pytest.approx(0.0, abs=1e-12)
    assert col[II] == pytest.approx(30.0)
    s0 = np.asarray(frac) * TL
    assert (col[ISSOIL] - s0).tolist() == pytest.approx([6.0, 3.0, 21.0])
    # booked as percolation through the horizons it crossed
    assert np.asarray(col[IRP])[:2].tolist() == pytest.approx([24.0, 21.0])
    # nothing above saturation
    assert np.all(col[ISSOIL] <= SM * TL + 1e-9)


def test_only_the_irrigation_goes_down():
    """Rain keeps Eq. 1b: of a 30 mm day with 5 mm of irrigation, only those
    5 mm can go below the top horizon."""
    frac = [0.38, 0.30, 0.30]
    top = _flux(30.0, 0.0, frac)
    col = _flux(30.0, 5.0, frac)
    assert top[IRO] - col[IRO] == pytest.approx(5.0)
    assert col[II] - top[II] == pytest.approx(5.0)


def test_a_full_column_still_runs_off():
    frac = [0.41, 0.40, 0.375]             # 2.5 mm of room, at the bottom
    col = _flux(20.0, 20.0, frac)
    assert col[IRO] == pytest.approx(17.5)
    assert np.allclose(col[ISSOIL], SM * TL)


def test_one_horizon_has_nowhere_to_send_it():
    col = _flux(10.0, 10.0, [0.38], nsl=1)
    assert col[IRO] == pytest.approx(10.0 - 0.03 * 200.0)


def test_zero_irrigation_is_the_old_law_exactly():
    frac = [0.38, 0.30, 0.30]
    a = _flux(30.0, 0.0, frac)
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    assert 'IRR' in mm.flux.__code__.co_varnames
    for k in (IRO, II):
        assert a[k] == pytest.approx(30.0 - 6.0 if k == IRO else 6.0)


def test_run_on_reinfiltration_counts_what_it_displaced_down():
    """With CRR, run-on fills the top horizon too, so more of the cell's own
    irrigation goes deeper. The reinfiltration is what the soil took BEYOND
    the cell's own water, so the catchment identity still holds:
    I(with run-on) - I(without) = REinf."""
    frac = [0.30, 0.30, 0.30]               # top room 22 mm
    own = _flux(20.0, 20.0, frac)
    both = _flux(20.0, 20.0, frac, runon=15.0)
    assert both[IREINF] == pytest.approx(both[II] - own[II])
    assert both[IRO] == pytest.approx(own[IRO] + 15.0 - both[IREINF])


# ------------------------------------------------------- an irrigated cell
def _irrigated(P_rain, irr, frac_top, irr_column):
    cMF = T._FakeMF(nrow=1, ncol=2, nper=1, perlen=[1])
    cMF.modelname = 'toy'
    cMF.outcropL[:] = 1
    inp = T._build_inputs(cMF)
    inp['P_veg_zoneSP'][:] = P_rain
    inp['Pe_veg_zonesSP'][:] = 0.9 * P_rain
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    cells = mm.build_cell_list(cMF)
    one = np.ones((1, 1, 1))
    ctx = mm.build_context(
        cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'], inp['_Sm'],
        inp['_Sfc'], inp['_Sr'], inp['_slprop'], inp['_Ssoil_ini'],
        inp['botm_l0'], inp['_Ks'], inp['gridSOIL'], inp['gridSOILthick'],
        inp['TopSoil'], inp['gridMETEO'], T.INDEX_MM, T.INDEX_MM_S,
        inp['P_veg_zoneSP'], inp['Eo_zonesSP'], inp['PT_veg_zonesSP'],
        inp['Pe_veg_zonesSP'], inp['PE_zonesSP'], inp['gridVEGarea'],
        inp['LAI_veg_zonesSP'], inp['Zr'], inp['kTg_min'], inp['kTg_max'],
        inp['kT_f'], inp['kT_s'], inp['NVEG'], 1000.0, 1,
        one * (P_rain + irr), one * 3.0, one * 0.95 * (P_rain + irr),
        np.ones((1, 1), int), np.array([[1, 0]]), [0.5], [0.0], [0.1],
        [0.5], [0.1])
    ctx.irr_column = irr_column
    st = mm.init_state(ctx)
    st.carried = True
    sm = np.asarray(inp['_Sm'][0])
    tl = 2000.0 * np.asarray(inp['_slprop'][0])
    st.Ssoil_ini[:, :2] = np.array([frac_top, 0.15]) * tl
    out = mm.step(ctx, 0, 0, np.full(2, 690.0), np.zeros(2), st)
    return out['MM'], out['MM_S'], sm * tl


@pytest.mark.parametrize('rain,irr', [(0.0, 25.0), (5.0, 25.0), (30.0, 2.0)])
def test_an_irrigated_cell_keeps_its_sprinkler_water(rain, irr):
    """Cell 0 is irrigated, cell 1 is not: only cell 0 changes, by exactly
    the irrigation the top horizon could not take, and every balance of the
    soil -- surface, column, each layer -- still closes."""
    top_mm, top_s, cap = _irrigated(rain, irr, 0.38, False)
    col_mm, col_s, _ = _irrigated(rain, irr, 0.38, True)
    assert np.array_equal(top_mm[1], col_mm[1])            # not irrigated
    pe = 0.95 * (rain + irr)
    room = cap[0] - 0.38 * 2000.0 * 0.3
    excess = max(pe - room, 0.0)
    to_depth = min(excess, pe * irr / (rain + irr))
    assert top_mm[0, IX['iRo']] == pytest.approx(excess)
    assert col_mm[0, IX['iRo']] == pytest.approx(excess - to_depth)
    assert col_mm[0, IX['iI']] - top_mm[0, IX['iI']] == pytest.approx(to_depth)
    for MM, MMS in ((top_mm, top_s), (col_mm, col_s)):
        assert abs(MM[0, IX['iMBsurf']]) < 1e-9
        assert abs(MM[0, IX['iMB']]) < 1e-9
        assert np.abs(MMS[0, :, IXS['iMB_s']]).max() < 1e-9


def test_the_panel_and_the_driver_carry_the_switch():
    import marmites_config as mcfg
    assert mcfg.RunConfig.from_dict({}).soil.irr_infiltration == 'top'
    with pytest.raises(mcfg.ConfigError) as e:
        mcfg.RunConfig.from_dict({'soil': {'irr_infiltration': 'flood'}})
    assert 'soil.irr_infiltration' in str(e.value)
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert "irr_column=(cfg.soil.irr_infiltration == 'column')" in src
    assert "ctx.irr_column = bool(getattr(a, 'irr_column', False))" in src
