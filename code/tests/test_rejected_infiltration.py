# -*- coding: utf-8 -*-
"""UZF's rejected infiltration, and the cap on the soil's percolation.

MODFLOW's UZF can refuse part of the percolation MARMITES delivers (its
unsaturated zone is full, or the rate exceeds VKS). The first coupled La
Mata run showed ~17 % refused, and that water vanished from the balance;
from 2026-09 it went back into the BOTTOM soil layer the next day, like
groundwater exfiltration (Frances & Lubczynski 2023, Eq. 1 and sec. 2.3).

THE FLIP-FLOP (2026-10-17). The bottom layer percolates only on a day when
nothing enters it from below, so the water returned one day stopped the
next day's percolation, and the day after it percolated everything again,
which UZF refused again: every rejecting cell alternated day by day. 25 %
of the base run's percolation went round that loop, 90-95 % in calibrated
parameter sets, and in one realisation it drove UZF's ET to 184 mm/d.

THE RULE SINCE 2026-10-18 (analysis doc 8.26):

  * the bottom layer percolates no more than the unsaturated zone below can
    take that day, PERC_CAP (the coupler's estimate: UZF's VKS and the room
    above the water table); the soil keeps the rest;
  * what UZF still refuses -- the small remainder -- RUNS OFF (with the CRR
    cascade, downslope); it no longer enters the soil from below.
"""
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for _p in (TRUNK, HERE, os.path.join(TRUNK, 'MARMITESsoil')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T = _load('t_runmmsoil_rej', os.path.join(HERE, 'test_runmmsoil.py'))
F = _load('t_soilflux_rej', os.path.join(HERE, 'test_soil_flux.py'))

# indices into the flux() return tuple
IEOW, ISSURF, IRO, IRP, IESOIL, ITSOIL, ISSOIL = 0, 1, 2, 3, 4, 5, 6
ISAT, IREXF, II = 12, 13, 14


def _call(rejinf=0.0, Pe=0.0, exf=0.0, soil_frac=(0.20, 0.22), nsl=2,
          cap=np.inf):
    """One cell, one day, with explicit exfiltration, rejected infiltration
    and percolation cap."""
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    cMF = T._FakeMF(nper=1, perlen=[1])
    Sm, Sfc, Sr, Ks, Tl = F._column(nsl)
    TopSoil = 700_000.0
    TopSoilLay = np.zeros(nsl); BotSoilLay = np.zeros(nsl)
    for l in range(nsl):
        TopSoilLay[l] = TopSoil if l == 0 else BotSoilLay[l - 1]
        BotSoilLay[l] = TopSoilLay[l] - Tl[l]
    NVEG = 2
    Zr_elev = [TopSoil - 600.0, TopSoil - 1500.0]
    Ssoil_ini = [soil_frac[l] * Tl[l] for l in range(nsl)]
    return mm.flux(cMF, 1.0, Pe, np.zeros(NVEG), 0.0, Zr_elev,
                   np.array([0.0, 0.0]), TopSoil - 5000.0, TopSoilLay, BotSoilLay,
                   Tl, nsl, Sm, Sfc, Sr, Ks,
                   Ssoil_ini, exf, 5000.0, 'loam',
                   0, 0, 0, [0.0] * NVEG, [0.9, 0.7], [0.5] * NVEG, [0.1] * NVEG,
                   NVEG, np.array([2.0, 1.5]), REJINF_ini=rejinf,
                   PERC_CAP=cap)


def _f(x):
    return np.asarray(x, dtype=float)


# --------------------------------------------------------------------- #
# the refused water runs off
# --------------------------------------------------------------------- #

def test_rejected_water_runs_off_and_the_soil_does_not_see_it():
    base = _call(rejinf=0.0, Pe=3.0)
    rej = _call(rejinf=50.0, Pe=3.0)
    assert float(rej[IRO]) - float(base[IRO]) == np.float64(50.0)
    for k in (ISSOIL, IRP, IREXF, IESOIL, ITSOIL, II):
        assert np.allclose(_f(rej[k]), _f(base[k])), k


def test_returned_water_no_longer_stops_the_bottom_percolation():
    """The flip-flop's cause: a wet column refused yesterday must still
    percolate today. The old rule gave Rp[-1] = 0 here."""
    wet = (0.35, 0.38)
    base = _call(rejinf=0.0, soil_frac=wet)
    rej = _call(rejinf=30.0, soil_frac=wet)
    assert float(_f(base[IRP])[-1]) > 0.0
    assert float(_f(rej[IRP])[-1]) == float(_f(base[IRP])[-1])


def test_negative_rejinf_is_treated_as_magnitude():
    """MF6 reports rejection with a negative sign."""
    assert np.isclose(float(_call(rejinf=40.0)[IRO]),
                      float(_call(rejinf=-40.0)[IRO]))


# --------------------------------------------------------------------- #
# the cap: the soil keeps what UZF cannot take
# --------------------------------------------------------------------- #

def test_the_bottom_percolation_is_capped_and_the_soil_keeps_the_rest():
    wet = (0.35, 0.38)
    free = _call(soil_frac=wet)
    rp = float(_f(free[IRP])[-1])
    assert rp > 1.0
    capped = _call(soil_frac=wet, cap=0.5)
    assert float(_f(capped[IRP])[-1]) == 0.5
    kept = float(np.sum(_f(capped[ISSOIL]))) - float(np.sum(_f(free[ISSOIL])))
    assert np.isclose(kept, rp - 0.5)
    assert float(capped[IRO]) == float(free[IRO])         # nothing ran off


def test_a_cap_above_the_percolation_changes_nothing():
    wet = (0.35, 0.38)
    free = _call(soil_frac=wet, Pe=5.0)
    for cap in (np.inf, 1.0e6):
        same = _call(soil_frac=wet, Pe=5.0, cap=cap)
        for x, y in zip(free, same):
            if x is None and y is None:      # the EVT route's extra, off here
                continue
            assert np.array_equal(_f(x), _f(y))


def test_a_zero_cap_stops_the_percolation():
    out = _call(soil_frac=(0.35, 0.38), cap=0.0)
    assert float(_f(out[IRP])[-1]) == 0.0


# --------------------------------------------------------------------- #
# groundwater exfiltration keeps its route
# --------------------------------------------------------------------- #

def test_groundwater_exfiltration_still_enters_from_below():
    base = _call()
    exf = _call(exf=50.0)
    stored = float(np.sum(_f(exf[ISSOIL]))) - float(np.sum(_f(base[ISSOIL])))
    assert stored > 0.0
    assert float(_f(exf[IRP])[-1]) == 0.0       # no percolation that day


# --------------------------------------------------------------------- #
# conservation
# --------------------------------------------------------------------- #

def test_surface_mass_balance_closes_with_the_refused_water():
    """dSsurf/dt = Pe + Exf_g1 + Rrej - I - Eow - Ro."""
    for rej in (0.0, 20.0, 400.0):
        for exf in (0.0, 900.0):
            out = _call(rejinf=rej, Pe=4.0, exf=exf)
            Eow, Ssurf, Ro = (float(out[IEOW]), float(out[ISSURF]),
                              float(out[IRO]))
            I, Exf_g1 = float(out[II]), float(_f(out[IREXF])[0])
            mb = (4.0 + Exf_g1 + rej) - (Eow + Ro + I + Ssurf)
            assert abs(mb) < 1e-9, (rej, exf, mb)


def test_soil_column_mass_balance_closes():
    """In (I) minus out equals the storage change -- the refused water is
    not an input of the soil any more, the cap only keeps water in it."""
    Sm, _Sfc, _Sr, _Ks, Tl = F._column(2)
    for rej, cap in ((0.0, np.inf), (25.0, np.inf), (25.0, 0.3)):
        out = _call(rejinf=rej, Pe=3.0, soil_frac=(0.30, 0.33), cap=cap)
        Ssoil, Rexf, Rp = _f(out[ISSOIL]), _f(out[IREXF]), _f(out[IRP])
        Es, Ts = _f(out[IESOIL]), _f(out[ITSOIL])
        ini = np.array([0.30 * Tl[0], 0.33 * Tl[1]])
        dS = float(np.sum(Ssoil - ini))
        mb = float(out[II]) - (Rexf[0] + Es.sum() + Ts.sum() + dS + Rp[-1])
        assert abs(mb) < 1e-6, (rej, cap, mb)


# --------------------------------------------------------------------- #
# through step(): per cell, compiled or not
# --------------------------------------------------------------------- #

def _ctx_state(mm, cMF):
    inp = T._build_inputs(cMF)
    cells = mm.build_cell_list(cMF)
    c = mm.build_context(cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'],
                         inp['_Sm'], inp['_Sfc'], inp['_Sr'], inp['_slprop'],
                         inp['_Ssoil_ini'], inp['botm_l0'], inp['_Ks'],
                         inp['gridSOIL'], inp['gridSOILthick'], inp['TopSoil'],
                         inp['gridMETEO'], T.INDEX_MM, T.INDEX_MM_S,
                         inp['P_veg_zoneSP'], inp['Eo_zonesSP'],
                         inp['PT_veg_zonesSP'], inp['Pe_veg_zonesSP'],
                         inp['PE_zonesSP'], inp['gridVEGarea'],
                         inp['LAI_veg_zonesSP'], inp['Zr'], inp['kTg_min'],
                         inp['kTg_max'], inp['kT_f'], inp['kT_s'], inp['NVEG'],
                         1000.0, 0, None, None, None, None, None,
                         None, None, None, None, None)
    return c, mm.init_state(c)


def test_step_takes_the_rejection_and_the_cap_per_cell():
    cMF = T._FakeMF(nper=1, perlen=[1])
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    ro, pi = T.INDEX_MM['iRo'], T.INDEX_MM['iperc']
    runs = {}
    for name, kw in (('base', {}),
                     ('rej', {'rejinf_cell': 7.0}),
                     ('cap0', {'perc_cap': 0.0})):
        c, s = _ctx_state(mm, cMF)
        h = np.full(c.ncell, 699.0)
        z = np.zeros(c.ncell)
        kw = {k: np.full(c.ncell, v) for k, v in kw.items()}
        runs[name] = (mm.step(c, 0, 0, h, z, s, **kw), s, c)
    base, rej, cap0 = (np.asarray(runs[k][0]['MM']) for k in ('base', 'rej',
                                                              'cap0'))
    c = runs['base'][2]
    fs = 1.0 - np.asarray(getattr(c, 'f_open', np.zeros(c.ncell)), float)
    # the refused water (per cell area) is the cell's runoff wherever the
    # cell has a soil share to have refused it from
    assert np.allclose(rej[:, ro] - base[:, ro], np.where(fs > 0, 7.0, 0.0))
    assert np.allclose(rej[:, pi], base[:, pi])
    # a zero cap reaches every cell: nothing percolates (what the soil
    # keeps is tested on flux() above; this fixture's first day is dry)
    assert np.all(cap0[:, pi] == 0.0)
    assert np.all(np.asarray(runs['cap0'][0]['perc']) == 0.0)
