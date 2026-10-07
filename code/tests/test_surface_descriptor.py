# -*- coding: utf-8 -*-
"""WP4.6 / WP4.4 -- open water bypasses the MM soil column (cookbook §4a).

A cell is part soil, part open water: f_lake (a pond), f_stream (the
channel's share). The MM soil column runs on the soil fraction only, its
fluxes per unit of CELL area are the column's times f_soil, and over the open
fraction the rain goes straight to the water body as runoff -- which the
coupler hands to SFR (INFLOW) and LAK (RUNOFF). Evaporation there is MF6's.
"""
import copy
import importlib.util
import os
import sys
from types import SimpleNamespace

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


T = _load('t_runmmsoil_sd', os.path.join(HERE, 'test_runmmsoil.py'))
coup = _load('marmites_coupler_sd', os.path.join(TRUNK, 'marmites_coupler.py'))
mf6mod = _load('marmites_mf6_sd', os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'))

IX = T.INDEX_MM
EXTENSIVE = ('iPT', 'iPE', 'iEi', 'iETsoil', 'iEg', 'iTg', 'iETg', 'iperc',
             'iI', 'iPETuzf', 'idSsoil')


def _context(nper=3):
    cMF = T._FakeMF(nper=nper, perlen=[1] * nper)
    cMF.modelname = 'toy'
    inp = T._build_inputs(cMF)
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    cells = mm.build_cell_list(cMF)
    ctx = mm.build_context(cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'],
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
    return mm, ctx, cMF


def _run(mm, ctx, f_open, nper=3, exf=0.0):
    ctx.f_open = f_open
    state = mm.init_state(ctx)
    heads = np.asarray([float(ctx.TopSoil[c[1], c[2]]) * 1e-3 - 3.0
                        for c in ctx.cells])
    outs = []
    for n in range(nper):
        outs.append(mm.step(ctx, n, n, heads, np.full(ctx.ncell, exf), state))
    return outs, state


def test_a_pond_cell_has_no_soil_column_and_its_rain_goes_to_the_pond():
    mm, ctx, _ = _context()
    base, s0 = _run(mm, ctx, None)
    fo = np.zeros(ctx.ncell)
    fo[1] = 1.0
    outs, s1 = _run(mm, ctx, fo)
    for o, b in zip(outs, base):
        MM, MB = o['MM'], b['MM']
        for k in EXTENSIVE:
            assert MM[1, IX[k]] == 0.0, k
        assert o['perc'][1] == 0.0 and o['etg'][1] == 0.0
        assert o['petuzf'][1] == 0.0
        assert MM[1, IX['iRo']] == pytest.approx(MM[1, IX['iP']])
        # every other cell is exactly as it was
        others = [c for c in range(ctx.ncell) if c != 1]
        assert np.array_equal(MM[others], MB[others])
    # the column's own state goes on as it would: it is only not counted
    assert np.array_equal(s1.Ssoil_ini, s0.Ssoil_ini)


def test_a_channel_cell_keeps_its_soil_column_on_the_land_share():
    mm, ctx, _ = _context()
    base, _ = _run(mm, ctx, None)
    fo = np.zeros(ctx.ncell)
    fo[2] = 0.3
    outs, _ = _run(mm, ctx, fo)
    for o, b in zip(outs, base):
        MM, MB = o['MM'], b['MM']
        for k in EXTENSIVE:
            assert MM[2, IX[k]] == pytest.approx(0.7 * MB[2, IX[k]]), k
        assert o['perc'][2] == pytest.approx(0.7 * b['perc'][2])
        assert o['etg'][2] == pytest.approx(0.7 * b['etg'][2])
        assert MM[2, IX['iRo']] == pytest.approx(
            0.7 * MB[2, IX['iRo']] + 0.3 * MB[2, IX['iP']])
        # states are the column's, not scaled
        for k in ('iSsoil_pc', 'idgwt'):
            assert MM[2, IX[k]] == pytest.approx(MB[2, IX[k]])


def test_rain_reaches_the_ground_or_the_water_whole():
    """P = Ei + Pe per cell: interception only over the soil fraction, and
    the rain on the open fraction reaches the water body undiminished."""
    mm, ctx, _ = _context()
    fo = np.linspace(0.0, 1.0, ctx.ncell)
    outs, _ = _run(mm, ctx, fo)
    for o in outs:
        MM = o['MM']
        assert np.allclose(MM[:, IX['iP']],
                           MM[:, IX['iEi']] + MM[:, IX['iPe']], atol=1e-9)


def test_seepage_under_open_water_goes_to_the_water_body():
    mm, ctx, cMF = _context()
    fo = np.zeros(ctx.ncell)
    fo[0] = 1.0
    outs, _ = _run(mm, ctx, fo, exf=2.0)
    for n, o in enumerate(outs):
        MM = o['MM']
        assert MM[0, IX['iRo']] == pytest.approx(
            MM[0, IX['iP']] + 2.0 / cMF.perlen[n])


# ------------------------------ the builder ------------------------------ #

def _builder(ponds=(), reaches=None):
    b = SimpleNamespace(ponds=list(ponds), sfr=reaches is not None,
                        sfr_net=None)
    if reaches is not None:
        cells = [c for c, _w, _l in reaches]
        b.sfr_net = SimpleNamespace(
            rno={c: k for k, c in enumerate(cells)},
            reach_wid=[w for _c, w, _l in reaches],
            reach_len=[ln for _c, _w, ln in reaches])
    return b


def test_the_descriptor_spreads_a_pond_over_its_footprint():
    pond = SimpleNamespace(cells=[(0, 0), (1, 0)], area=150.0)
    b = _builder([pond], reaches=[((2, 0), 2.0, 10.0), ((3, 0), 30.0, 10.0)])
    cells = [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)]
    fl, fs = mf6mod.clsMF6.surface_fractions(b, cells, np.full(5, 100.0))
    assert fl.tolist() == pytest.approx([0.75, 0.75, 0.0, 0.0, 0.0])
    assert fs.tolist() == pytest.approx([0.0, 0.0, 0.2, 1.0, 0.0])  # capped


def test_a_resolved_pond_is_all_water_and_a_sub_grid_one_is_its_share():
    big = SimpleNamespace(cells=[(0, 0), (1, 0)], area=260.0)   # > footprint
    small = SimpleNamespace(cells=[(2, 0)], area=40.0)          # sub-grid
    b = _builder([big, small])
    fl, fs = mf6mod.clsMF6.surface_fractions(
        b, [(0, 0), (1, 0), (2, 0)], np.full(3, 100.0))
    assert fl.tolist() == pytest.approx([1.0, 1.0, 0.4])
    assert not fs.any()


# ------------------------------ the coupler ------------------------------ #

class _Cpl(coup.MF6Coupler):
    """The surface-water writes without the rest of the coupler."""

    def __init__(self, area, lak_cells=(), reach_idx=None):
        self.ncell = len(area)
        self.area = np.asarray(area, dtype=float)
        self.conv_fact = 1000.0
        self.nlakes = len(lak_cells)
        self.lak_cells = [(np.asarray(k, int), np.asarray(w, float))
                          for k, w in lak_cells]
        self.lak_cell_idx = np.full(max(self.nlakes, 1), -1, dtype=int)
        self.nreaches = 0 if reach_idx is None else int(max(reach_idx) + 1)
        self.sfr_reach_idx = (np.full(self.ncell, -1, dtype=int)
                              if reach_idx is None else np.asarray(reach_idx))
        self.p_sfr_inflow = (None if reach_idx is None
                             else np.zeros(self.nreaches))
        self.p_lak_runoff = np.zeros(self.nlakes) if self.nlakes else None
        self.p_sfr_simevap = None
        self.p_lak_simevap = None
        self._iRo = 0


def test_a_pond_takes_the_runoff_of_its_whole_footprint():
    """WP4.4: pond inflow follows the MM runoff -- here the rain on its
    cells, which MMsoil hands over as runoff."""
    c = _Cpl([100.0, 100.0, 50.0, 100.0],
             lak_cells=[([0, 1], [0.5, 0.5]), ([2], [1.0])],
             reach_idx=[-1, -1, -1, 0])
    ro = np.array([2.0, 4.0, 10.0, 3.0])              # mm/d
    c._write_runoff(ro)
    assert c.p_lak_runoff.tolist() == pytest.approx([0.6, 0.5])   # m3/d
    assert c.p_sfr_inflow.tolist() == pytest.approx([0.3])
    c._write_runoff(2 * ro)
    assert c.p_lak_runoff.tolist() == pytest.approx([1.2, 1.0])


def test_a_ponds_evaporation_is_booked_over_its_footprint():
    c = _Cpl([100.0, 300.0, 100.0], lak_cells=[([0, 1], [0.25, 0.75])])
    c.p_lak_simevap = np.array([-0.4])                 # m3/d, MF6's sign
    _sfr, lak = c._read_openwater_evap_split()
    assert lak.tolist() == pytest.approx([1.0, 1.0, 0.0])          # mm/d
    assert float(np.sum(lak * c.area)) / 1000.0 == pytest.approx(0.4)


# --------------------------- La Mata, 50 m grid ------------------------- #

def test_lamata_structured_ponds_are_sub_grid_and_streams_skirt_them(tmp_path):
    """On the 50 m grid a pond is its host cell's share (0.14-0.81), and no
    cell is both pond and channel: the reaches in a footprint were cut."""
    lk = _load('test_lak_ponds_sd', os.path.join(HERE, 'test_lak_ponds.py'))
    if not os.path.exists(lk.SHP):
        pytest.skip('pond shapefile not present')
    pytest.importorskip('flopy')
    import matplotlib
    matplotlib.use('agg')
    import lamata_model
    import marmites_channel as mch
    import marmites_config as cfgmod
    import marmites_vector as mv
    DS = lk.DS
    # La Mata's model description as the run builds it -- no
    # parameter file (lamata_model derives the outcrop layer too)
    c = lamata_model.lamata_cmf()
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]
    b = mf6mod.clsMF6(c, top=np.asarray(c.elev, float),
                      botm=np.asarray(c.botm, float), sim_ws=str(tmp_path),
                      daily=True)
    b.verbose = False
    lines, seg = mch.read_stream_lines(os.path.join(DS, 'inputSTREAM.csv'),
                                       os.path.join(DS, 'inputSTREAM_param.csv'))
    vg = mv.TargetGrid.structured(c.delr, c.delc, c.xllcorner, c.yllcorner)
    present, seg_of_cell, ch_len = mch.burn_channel(lines, vg, (c.nrow, c.ncol))
    act = np.asarray(c.outcropL) > 0
    b.sfr_pondw = np.where(act, present, 0.0)
    b.sfr_pondhmax = np.zeros_like(b.sfr_pondw)
    b.sfr_seg_of_cell, b.sfr_seg_params = seg_of_cell, seg
    b.sfr_cell_length = ch_len
    b.cell_area = 2500.0
    b.sfr_width_source = cfgmod.ParamSource(
        drainage={'w_min': 1.5, 'w_max': 3.0, 'power': 2.0})
    b.sfr_depth_source = cfgmod.ParamSource(value=1.0)
    b.lak_shapefile = lk.SHP
    b.lak_depth = 1.5
    b.build()
    cells = [(int(i), int(j)) for i, j in zip(*np.where(act))]
    fl, fs = b.surface_fractions(cells, np.full(len(cells), 2500.0))
    pos = {cc: k for k, cc in enumerate(cells)}
    for p in b.ponds:
        assert fl[pos[p.cell]] == pytest.approx(p.area / 2500.0)
    assert 0.1 < fl[fl > 0].min() and fl.max() < 0.85
    assert not np.any((fl > 0) & (fs > 0))
    assert np.all(fl + fs <= 1.0 + 1e-12)
    assert fs.max() <= 1.0 and (fs > 0).sum() == b.sfr_net.nreaches
