# -*- coding: utf-8 -*-
"""Phase-2 integration tests for the grid-agnostic soil driver.

Exercise runMMsoil() and step() on a small synthetic 2x3 grid with an
in-memory HDF5, verifying:
  * the run completes and fills the structured output;
  * the soil-column mass balance closes at every active cell / SP;
  * step() (the MF6-API interface) reproduces runMMsoil's percolation/ETg
    exactly (same physics, different entry point);
  * inactive cells are masked with hnoflo.
"""
import importlib.util
import os

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.join(HERE, '..')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


new = _load('mmsoil_new2', os.path.join(TRUNK, 'MARMITESsoil', 'MARMITESsoil_v3.py'))
proc = _load('mmproc_new2', os.path.join(TRUNK, 'MARMITESutilities', 'MARMITESprocess_v3.py'))
_idx = _load('marmites_indices', os.path.join(TRUNK, 'marmites_indices.py'))

HNOFLO = -999.9
HDRY = 1.0e30

# shared index maps (single source of truth)
INDEX_MM = _idx.INDEX_MM
INDEX_MM_S = _idx.INDEX_MM_SOIL


class _FakeProcess:
    def __init__(self, sy, nlay, nrow, ncol):
        self._sy = sy
        self.nlay, self.nrow, self.ncol = nlay, nrow, ncol

    def float2array(self, _):
        return np.full((self.nlay, self.nrow, self.ncol), self._sy)


class _FakeMF:
    def __init__(self, nrow=2, ncol=3, nper=4, sy=0.05, uzf_yn=1, wel_yn=1, perlen=None):
        self.nrow, self.ncol, self.nlay = nrow, ncol, 1
        self.nper = nper
        # default mixes multi-day SPs (exercises carry-over); daily=[1]*nper
        self.perlen = perlen if perlen is not None else [1, 3, 1, 2][:nper]
        self.nstp = [1] * nper
        self.hnoflo = HNOFLO
        self.hdry = HDRY
        self.uzf_yn = uzf_yn
        self.wel_yn = wel_yn
        self.delr = [100.0] * ncol
        self.delc = [100.0] * nrow
        # first cell of the grid is inactive, the rest outcrop in layer 1
        self.outcropL = np.ones((nrow, ncol), dtype=int)
        self.outcropL[0, 0] = 0
        self.sy_actual = sy
        self.cPROCESS = _FakeProcess(sy, self.nlay, nrow, ncol)


def _build_inputs(cMF, nmeteo=1, nveg=2, nsoil=1):
    nper = cMF.nper
    rng = np.random.default_rng(42)
    # one soil zone, two layers (loam), everywhere
    Sm = [0.41, 0.41]; Sfc = [0.30, 0.30]; Sr = [0.05, 0.05]; Ks = [50.0, 20.0]
    slprop = [np.array([0.3, 0.7])]
    _nsl = [2]
    _st = ['loam']
    _Sm = [Sm]; _Sfc = [Sfc]; _Sr = [Sr]; _Ks = [Ks]
    _Ssoil_ini = [[0.20, 0.22]]
    _nslmax = 2
    grid1 = np.ones((cMF.nrow, cMF.ncol), dtype=int)
    gridSOIL = grid1.copy()
    gridMETEO = grid1.copy()
    gridSOILthick = np.full((cMF.nrow, cMF.ncol), 2.0)   # 2 m soil
    TopSoil = np.full((cMF.nrow, cMF.ncol), 700_000.0)   # mm elevation
    botm_l0 = np.full((cMF.nrow, cMF.ncol), 650.0)
    gridSsurfhmax = np.full((cMF.nrow, cMF.ncol), 0.1)
    gridSsurfw = np.full((cMF.nrow, cMF.ncol), 1.0)
    gridVEGarea = np.zeros((nveg, cMF.nrow, cMF.ncol))
    gridVEGarea[0] = 40.0
    gridVEGarea[1] = 30.0
    # meteo/veg time series over stress periods (index by tstart_MF = SP index here)
    P_veg = rng.uniform(0, 20, size=(nmeteo, nper))
    Eo = np.full((nmeteo, nper), 4.0)
    PT_veg = np.full((nmeteo, nveg, nper), 3.0)
    Pe_veg = P_veg[:, None, :] * 0.9 * np.ones((1, nveg, 1))
    LAI_veg = np.full((nveg, nper), 2.0)
    PE = np.full((nmeteo, nsoil, nper), 3.5)
    Zr = [0.6, 1.5]
    kTg_min = [0.0, 0.0]; kTg_max = [0.9, 0.7]; kT_f = [0.5, 0.5]; kT_s = [0.1, 0.1]
    return dict(_nsl=_nsl, _nslmax=_nslmax, _st=_st, _Sm=_Sm, _Sfc=_Sfc, _Sr=_Sr, _slprop=slprop,
                _Ssoil_ini=_Ssoil_ini, botm_l0=botm_l0, _Ks=_Ks, gridSOIL=gridSOIL,
                gridSOILthick=gridSOILthick, TopSoil=TopSoil, gridMETEO=gridMETEO,
                gridSsurfhmax=gridSsurfhmax, gridSsurfw=gridSsurfw, P_veg_zoneSP=P_veg,
                Eo_zonesSP=Eo, PT_veg_zonesSP=PT_veg, Pe_veg_zonesSP=Pe_veg, PE_zonesSP=PE,
                gridVEGarea=gridVEGarea, LAI_veg_zonesSP=LAI_veg, Zr=Zr, kTg_min=kTg_min,
                kTg_max=kTg_max, kT_f=kT_f, kT_s=kT_s, NVEG=nveg)


def _make_h5_MF(path, cMF, heads_val=699_000.0 / 1000.0, exf_val=0.0):
    """heads in m; exf4MM stored as MODFLOW volumetric (converted back inside)."""
    ndays = sum(cMF.perlen)
    h5 = h5py.File(path, 'w')
    heads = np.full((ndays, cMF.nrow, cMF.ncol), heads_val, dtype=np.float32)
    h5.create_dataset('heads4MM', data=heads)
    # exf4MM: volumetric; runMMsoil multiplies by conv_fact/(delr*delc)
    exf = np.full((ndays, cMF.nrow, cMF.ncol), exf_val, dtype=np.float32)
    h5.create_dataset('exf4MM', data=exf)
    return h5


def _make_h5_MM(path, cMF, nslmax):
    h5 = h5py.File(path, 'w')
    nd = sum(cMF.perlen)
    h5.create_dataset('MM', shape=(nd, cMF.nrow, cMF.ncol, len(INDEX_MM)), dtype=np.float32)
    h5.create_dataset('MM_S', shape=(nd, cMF.nrow, cMF.ncol, nslmax, len(INDEX_MM_S)), dtype=np.float32)
    h5.create_dataset('perc', shape=(cMF.nper, cMF.nrow, cMF.ncol), dtype=np.float32)
    h5.create_dataset('ETg', shape=(cMF.nper, cMF.nrow, cMF.ncol), dtype=np.float32)
    return h5


def _run(tmp, cMF, exf_val=0.0):
    inp = _build_inputs(cMF)
    conv_fact = 1000.0
    h5_MF = _make_h5_MF(os.path.join(tmp, '_h5_MF.h5'), cMF, exf_val=exf_val)
    h5_MM = _make_h5_MM(os.path.join(tmp, '_h5_MM.h5'), cMF, inp['_nslmax'])
    mm = new.clsMMsoil(hnoflo=HNOFLO)
    mm.runMMsoil(inp['_nsl'], inp['_nslmax'], inp['_st'], inp['_Sm'], inp['_Sfc'], inp['_Sr'],
                 inp['_slprop'], inp['_Ssoil_ini'], inp['botm_l0'], inp['_Ks'],
                 inp['gridSOIL'], inp['gridSOILthick'], inp['TopSoil'] , inp['gridMETEO'],
                 INDEX_MM, INDEX_MM_S, inp['gridSsurfhmax'], inp['gridSsurfw'],
                 inp['P_veg_zoneSP'], inp['Eo_zonesSP'], inp['PT_veg_zonesSP'],
                 inp['Pe_veg_zonesSP'], inp['PE_zonesSP'], inp['gridVEGarea'],
                 inp['LAI_veg_zonesSP'], inp['Zr'], inp['kTg_min'], inp['kTg_max'],
                 inp['kT_f'], inp['kT_s'], inp['NVEG'], cMF, conv_fact, h5_MF, h5_MM, irr_yn=0,
                 verbose=1)
    h5_MF.close()
    return os.path.join(tmp, '_h5_MM.h5'), inp


def test_runmmsoil_completes_and_masks_inactive(tmp_path):
    cMF = _FakeMF()
    mm_fn, inp = _run(str(tmp_path), cMF)
    with h5py.File(mm_fn, 'r') as h5:
        MM = h5['MM'][:]
        perc = h5['perc'][:]
    # inactive cell (0,0) masked with hnoflo
    assert np.allclose(MM[:, 0, 0, :], HNOFLO)
    # active cell (1,1) has finite, non-hnoflo P
    assert MM[0, 1, 1, INDEX_MM['iP']] != HNOFLO
    # percolation is finite everywhere
    assert np.all(np.isfinite(perc))
    assert perc[0, 0, 0] == 0.0  # inactive


def test_runmmsoil_soil_mass_balance_closes_multiday(tmp_path):
    # soil-column MB (iMB) closes for any SP length
    cMF = _FakeMF()  # mixed [1,3,1,2]
    mm_fn, inp = _run(str(tmp_path), cMF)
    with h5py.File(mm_fn, 'r') as h5:
        MM = h5['MM'][:]
    for i in range(cMF.nrow):
        for j in range(cMF.ncol):
            if cMF.outcropL[i, j] <= 0:
                continue
            assert np.max(np.abs(MM[:, i, j, INDEX_MM['iMB']])) < 1e-3, (i, j, 'soil')


def test_runmmsoil_full_mass_balance_closes_daily(tmp_path):
    # with daily stress periods (the Phase-3 target regime), BOTH the soil
    # and surface mass balances close. Surface MB does not close for
    # perlen>1 SPs because flux() runs one representative day then
    # rate-averages -- this is the SP-averaging approximation that the
    # daily-SP default (review section 4.3) removes.
    cMF = _FakeMF(perlen=[1, 1, 1, 1])
    mm_fn, inp = _run(str(tmp_path), cMF)
    with h5py.File(mm_fn, 'r') as h5:
        MM = h5['MM'][:]
    for i in range(cMF.nrow):
        for j in range(cMF.ncol):
            if cMF.outcropL[i, j] <= 0:
                continue
            assert np.max(np.abs(MM[:, i, j, INDEX_MM['iMB']])) < 1e-3, (i, j, 'soil')
            assert np.max(np.abs(MM[:, i, j, INDEX_MM['iMBsurf']])) < 1e-3, (i, j, 'surf')


def test_runmmsoil_mass_balance_with_exfiltration(tmp_path):
    # modest exfiltration (a few mm/d up into the soil), daily SPs
    # exf4MM volumetric; runMMsoil applies conv_fact/(delr*delc)=0.1,
    # then negates -> exf_cell = +3 mm/d
    cMF = _FakeMF(perlen=[1, 1, 1, 1])
    mm_fn, inp = _run(str(tmp_path), cMF, exf_val=-30.0)
    with h5py.File(mm_fn, 'r') as h5:
        MM = h5['MM'][:]
    for i in range(cMF.nrow):
        for j in range(cMF.ncol):
            if cMF.outcropL[i, j] <= 0:
                continue
            assert np.max(np.abs(MM[:, i, j, INDEX_MM['iMB']])) < 1e-3, (i, j)


def test_step_matches_runmmsoil(tmp_path):
    """step() (MF6-API entry) must reproduce runMMsoil's perc/ETg exactly."""
    cMF = _FakeMF()
    mm_fn, inp = _run(str(tmp_path), cMF)
    with h5py.File(mm_fn, 'r') as h5:
        perc_ref = h5['perc'][:]
        etg_ref = h5['ETg'][:]

    # replay through step() directly
    conv_fact = 1000.0
    mm = new.clsMMsoil(hnoflo=HNOFLO)
    cells = mm.build_cell_list(cMF)
    ctx = mm.build_context(cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'], inp['_Sm'],
                           inp['_Sfc'], inp['_Sr'], inp['_slprop'], inp['_Ssoil_ini'],
                           inp['botm_l0'], inp['_Ks'], inp['gridSOIL'], inp['gridSOILthick'],
                           inp['TopSoil'], inp['gridMETEO'], INDEX_MM, INDEX_MM_S,
                           inp['gridSsurfhmax'], inp['gridSsurfw'], inp['P_veg_zoneSP'],
                           inp['Eo_zonesSP'], inp['PT_veg_zonesSP'], inp['Pe_veg_zonesSP'],
                           inp['PE_zonesSP'], inp['gridVEGarea'], inp['LAI_veg_zonesSP'],
                           inp['Zr'], inp['kTg_min'], inp['kTg_max'], inp['kT_f'], inp['kT_s'],
                           inp['NVEG'], conv_fact, 0, None, None, None, None, None,
                           None, None, None, None, None)
    state = mm.init_state(ctx)
    i_arr = np.array([c[1] for c in cells])
    j_arr = np.array([c[2] for c in cells])
    heads_val = 699_000.0 / 1000.0
    tstart_MF = 0
    for n in range(cMF.nper):
        heads_cell = np.full(ctx.ncell, heads_val)
        exf_cell = np.zeros(ctx.ncell)
        out = mm.step(ctx, n, tstart_MF, heads_cell, exf_cell, state)
        # scatter and compare to the reference structured perc/ETg
        perc_grid = np.zeros((cMF.nrow, cMF.ncol), dtype=np.float32)
        etg_grid = np.zeros((cMF.nrow, cMF.ncol), dtype=np.float32)
        perc_grid[i_arr, j_arr] = out['perc']
        etg_grid[i_arr, j_arr] = out['etg']
        assert np.allclose(perc_grid, perc_ref[n], atol=1e-6), f'perc SP{n}'
        assert np.allclose(etg_grid, etg_ref[n], atol=1e-6), f'etg SP{n}'
        tstart_MF += cMF.nstp[n]
