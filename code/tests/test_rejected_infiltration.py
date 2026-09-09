# -*- coding: utf-8 -*-
"""Rejected-infiltration routing (Frances & Lubczynski 2023, Eq. 1 / sec. 2.3).

MODFLOW's UZF can refuse part of the percolation MARMITES delivers (the cell
is saturated, or the rate exceeds VKS). The first coupled La Mata run showed
~17% rejected, and that water previously vanished from the water balance.

Routing follows the published formulation, NOT a direct surface injection:

    Eq. 1   dSsurf/dt = Pe + Exf_g1 - I - Eow - Ro

the only subsurface->surface input being Exf_g1, i.e. exfiltration arriving
*from the topmost soil layer*. Per section 2.3, water the subsurface cannot
accept enters the soil column and "eventually creat[es] saturation-excess
overland flow (Dunnian flow) if all the soil layers turn saturated".

So rejected infiltration enters the BOTTOM soil layer, fills the column from
below, spills into the surface store, and only the excess above Ssurf_max
becomes runoff -- identical to the groundwater-exfiltration pathway.
"""
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, TRUNK)
sys.path.insert(0, HERE)


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


def _call(rejinf, Ssurf_max=10.0, Pe=0.0, Ssurf_ini=0.0, exf=0.0,
          soil_frac=(0.20, 0.22), nsl=2):
    """One cell, one day, with explicit exfiltration / rejected infiltration."""
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
    return mm.flux(cMF, 1.0, Pe, np.zeros(NVEG), 0.0, 4.0, Zr_elev,
                   np.array([0.0, 0.0]), TopSoil - 5000.0, TopSoilLay, BotSoilLay,
                   Tl, nsl, Sm, Sfc, Sr, Ks, Ssurf_max,
                   Ssoil_ini, Ssurf_ini, exf, 5000.0, 'loam',
                   0, 0, 0, [0.0] * NVEG, [0.9, 0.7], [0.5] * NVEG, [0.1] * NVEG,
                   NVEG, np.array([2.0, 1.5]), REJINF_ini=rejinf)


def _capacity(soil_frac=(0.20, 0.22), nsl=2):
    """Free pore space in the soil column [mm]."""
    Sm, _Sfc, _Sr, _Ks, Tl = F._column(nsl)
    return sum(Sm[l] * Tl[l] - soil_frac[l] * Tl[l] for l in range(nsl))


# --------------------------------------------------------------------- #
# the pathway: bottom layer -> up the column -> surface -> runoff
# --------------------------------------------------------------------- #

def test_rejected_water_enters_the_soil_column_from_below():
    """Below the column's free capacity it is stored in the soil, and
    reaches neither the surface nor runoff (no Dunnian flow yet)."""
    base = _call(rejinf=0.0)
    rej = _call(rejinf=50.0)
    stored = float(np.sum(rej[ISSOIL])) - float(np.sum(base[ISSOIL]))
    assert stored > 0.0, 'rejected water must be stored in the soil column'
    assert np.isclose(float(rej[IREXF][0]), 0.0), 'no Exf_g1 while unsaturated'
    assert np.isclose(float(rej[IRO]), 0.0), 'no runoff while unsaturated'


def test_saturating_the_column_produces_exf_g1_at_the_surface():
    """Once the column saturates, the excess appears as Exf_g1 = Rexf[0],
    the Eq. 1 surface input, and the SAT flags are raised."""
    cap = _capacity()
    out = _call(rejinf=cap + 200.0, Ssurf_max=1000.0)
    assert float(out[IREXF][0]) > 0.0, 'excess must reach the surface as Exf_g1'
    assert bool(np.asarray(out[ISAT])[-1]), 'bottom layer must be flagged saturated'
    assert float(out[ISSURF]) > 0.0, 'surface store must receive it'


def test_dunnian_runoff_once_soil_and_surface_are_full():
    """Soil saturated AND surface store full -> saturation-excess runoff."""
    cap = _capacity()
    out = _call(rejinf=cap + 500.0, Ssurf_max=1.0)
    assert float(out[IRO]) > 0.0, 'Dunnian runoff expected'
    assert float(out[ISSURF]) <= 1.0 + 1e-9, 'store cannot exceed Ssurf_max'


def test_more_rejection_gives_more_runoff():
    cap = _capacity()
    ro = [float(_call(rejinf=r, Ssurf_max=1.0)[IRO])
          for r in (0.0, cap + 100.0, cap + 600.0)]
    assert ro[0] == 0.0 and ro[1] < ro[2], ro


def test_rejected_infiltration_behaves_like_exfiltration():
    """Same quantity of water, delivered as exfiltration or as rejected
    infiltration, must give identical results -- both enter the column base."""
    a = _call(rejinf=0.0, exf=80.0, Ssurf_max=1.0)
    b = _call(rejinf=80.0, exf=0.0, Ssurf_max=1.0)
    for x, y in zip(a, b):
        assert np.allclose(np.asarray(x, dtype=float),
                           np.asarray(y, dtype=float)), 'pathways must coincide'


def test_bottom_layer_percolation_stops_when_subsurface_refuses():
    """If the subsurface cannot accept water, the soil must not keep
    percolating into it that period (Rp of the last layer is suppressed)."""
    out = _call(rejinf=30.0, soil_frac=(0.35, 0.38))
    assert np.isclose(float(np.asarray(out[IRP])[-1]), 0.0)


# --------------------------------------------------------------------- #
# conservation and compatibility
# --------------------------------------------------------------------- #

def test_surface_mass_balance_closes_eq1():
    """Eq. 1: dSsurf/dt = Pe + Exf_g1 - I - Eow - Ro."""
    cap = _capacity()
    for rej in (0.0, 20.0, cap + 50.0, cap + 800.0):
        out = _call(rejinf=rej, Ssurf_max=2.0, Pe=4.0)
        Eow, Ssurf, Ro = float(out[IEOW]), float(out[ISSURF]), float(out[IRO])
        I, Exf_g1 = float(out[II]), float(out[IREXF][0])
        mb = (4.0 + Exf_g1) - (Eow + Ro + I + (Ssurf - 0.0))
        assert abs(mb) < 1e-9, (rej, mb)


def test_soil_column_mass_balance_closes():
    """Water in (I + returned water) minus out equals the storage change."""
    cap = _capacity()
    for rej in (0.0, 25.0, cap + 300.0):
        out = _call(rejinf=rej, Ssurf_max=2.0, Pe=3.0)
        Ssoil = np.asarray(out[ISSOIL], dtype=float)
        Rexf = np.asarray(out[IREXF], dtype=float)
        Rp = np.asarray(out[IRP], dtype=float)
        Es = np.asarray(out[IESOIL], dtype=float)
        Ts = np.asarray(out[ITSOIL], dtype=float)
        Sm, _Sfc, _Sr, _Ks, Tl = F._column(2)
        ini = np.array([0.20 * Tl[0], 0.22 * Tl[1]])
        dS = float(np.sum(Ssoil - ini))
        mb = (float(out[II]) + rej) - (Rexf[0] + Es.sum() + Ts.sum() + dS + Rp[-1])
        assert abs(mb) < 1e-6, (rej, mb)


def test_zero_rejection_is_unchanged():
    """Default 0.0 must leave the uncoupled path bit-identical."""
    a = _call(rejinf=0.0, Pe=6.0)
    b = _call(rejinf=0.0, Pe=6.0)
    for x, y in zip(a, b):
        assert np.array_equal(np.asarray(x, dtype=float), np.asarray(y, dtype=float))


def test_negative_rejinf_is_treated_as_magnitude():
    """MF6 reports rejection with a negative sign."""
    cap = _capacity()
    pos = _call(rejinf=cap + 300.0, Ssurf_max=1.0)
    neg = _call(rejinf=-(cap + 300.0), Ssurf_max=1.0)
    assert np.isclose(float(pos[IRO]), float(neg[IRO]))


def test_step_accepts_rejinf_per_cell():
    """The coupler passes an (ncell,) array through step(); a large rejection
    must raise runoff once the columns saturate."""
    cMF = T._FakeMF(nper=1, perlen=[1])
    inp = T._build_inputs(cMF)
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    cells = mm.build_cell_list(cMF)

    def ctx_state():
        c = mm.build_context(cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'],
                             inp['_Sm'], inp['_Sfc'], inp['_Sr'], inp['_slprop'],
                             inp['_Ssoil_ini'], inp['botm_l0'], inp['_Ks'],
                             inp['gridSOIL'], inp['gridSOILthick'], inp['TopSoil'],
                             inp['gridMETEO'], T.INDEX_MM, T.INDEX_MM_S,
                             inp['gridSsurfhmax'], inp['gridSsurfw'],
                             inp['P_veg_zoneSP'], inp['Eo_zonesSP'],
                             inp['PT_veg_zonesSP'], inp['Pe_veg_zonesSP'],
                             inp['PE_zonesSP'], inp['gridVEGarea'],
                             inp['LAI_veg_zonesSP'], inp['Zr'], inp['kTg_min'],
                             inp['kTg_max'], inp['kT_f'], inp['kT_s'], inp['NVEG'],
                             1000.0, 0, None, None, None, None, None,
                             None, None, None, None, None)
        return c, mm.init_state(c)

    c0, s0 = ctx_state()
    h = np.full(c0.ncell, 699.0)
    z = np.zeros(c0.ncell)
    base = mm.step(c0, 0, 0, h, z, s0)
    c1, s1 = ctx_state()
    big = mm.step(c1, 0, 0, h, z, s1, rejinf_cell=np.full(c1.ncell, 5000.0))
    ro = T.INDEX_MM['iRo']
    assert np.all(np.asarray(big['MM'])[:, ro] >= np.asarray(base['MM'])[:, ro])
    assert np.any(np.asarray(big['MM'])[:, ro] > np.asarray(base['MM'])[:, ro])
