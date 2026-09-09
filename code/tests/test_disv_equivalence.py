# -*- coding: utf-8 -*-
"""Phase-4 primary acceptance test: DISV-from-DIS equivalence.

A DISV grid built from a (square-celled) DIS grid is geometrically the same
grid. Running the soil model through the DISV geometry must therefore
reproduce the DIS results exactly -- every flux, every stress period.

This is the correctness baseline that lets the unstructured path be trusted
before any true quadtree refinement exists: if DISV-from-DIS diverges, the
geometry plumbing is wrong, not the grid.
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


G = _load('marmites_grid_eq', os.path.join(TRUNK, 'marmites_grid.py'))
T = _load('t_runmmsoil_eq', os.path.join(HERE, 'test_runmmsoil.py'))


def _make_ctx(grid, cMF, inp, mm):
    cells = mm.build_cell_list(cMF)
    geom = G.geometry_for(cMF, cells, grid=grid)
    ctx = mm.build_context(cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'], inp['_Sm'],
                           inp['_Sfc'], inp['_Sr'], inp['_slprop'], inp['_Ssoil_ini'],
                           inp['botm_l0'], inp['_Ks'], inp['gridSOIL'], inp['gridSOILthick'],
                           inp['TopSoil'], inp['gridMETEO'], T.INDEX_MM, T.INDEX_MM_S,
                           inp['gridSsurfhmax'], inp['gridSsurfw'], inp['P_veg_zoneSP'],
                           inp['Eo_zonesSP'], inp['PT_veg_zonesSP'], inp['Pe_veg_zonesSP'],
                           inp['PE_zonesSP'], inp['gridVEGarea'], inp['LAI_veg_zonesSP'],
                           inp['Zr'], inp['kTg_min'], inp['kTg_max'], inp['kT_f'], inp['kT_s'],
                           inp['NVEG'], 1000.0, 0, None, None, None, None, None,
                           None, None, None, None, None, geom=geom)
    return ctx, mm.init_state(ctx), geom


def _march(grid, nper=4, heads=699.0, exf=0.0, seed_shift=0.0):
    """Run nper stress periods through step() on the given grid kind."""
    cMF = T._FakeMF(nper=nper, perlen=[1] * nper)      # square 100 x 100 cells
    inp = T._build_inputs(cMF)
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    ctx, state, geom = _make_ctx(grid, cMF, inp, mm)
    outs = []
    for n in range(nper):
        h = np.full(ctx.ncell, heads + seed_shift)
        e = np.full(ctx.ncell, exf)
        outs.append(mm.step(ctx, n, n, h, e, state))
    return outs, geom, state


def test_disv_geometry_matches_dis_on_square_grid():
    _, g_dis, _ = _march('dis', nper=1)
    _, g_disv, _ = _march('disv', nper=1)
    assert g_dis.kind == 'structured' and g_disv.kind == 'vertex'
    assert np.allclose(g_dis.area, g_disv.area)
    assert np.allclose(g_dis.width, g_disv.width)


def test_disv_reproduces_dis_fluxes_exactly():
    """Every output column, every SP: DISV == DIS on the equivalent grid."""
    out_dis, _, st_dis = _march('dis', nper=4)
    out_disv, _, st_disv = _march('disv', nper=4)
    assert len(out_dis) == len(out_disv) == 4
    for n, (a, b) in enumerate(zip(out_dis, out_disv)):
        for key in ('MM', 'MM_S', 'perc', 'etg'):
            assert np.array_equal(np.asarray(a[key]), np.asarray(b[key])), \
                f'SP{n} differs in {key}'
    # carry-over state must also match, otherwise later SPs would diverge
    assert np.array_equal(st_dis.Ssoil_ini, st_disv.Ssoil_ini)
    assert np.array_equal(st_dis.Ssurf_ini, st_disv.Ssurf_ini)


def test_disv_reproduces_dis_with_exfiltration():
    """The exfiltration path divides by cell area -- the most geometry-
    sensitive exchange; it must agree too."""
    out_dis, _, _ = _march('dis', nper=3, exf=2.5, heads=699.2)
    out_disv, _, _ = _march('disv', nper=3, exf=2.5, heads=699.2)
    for n, (a, b) in enumerate(zip(out_dis, out_disv)):
        assert np.array_equal(a['MM'], b['MM']), f'SP{n} MM differs'
        assert np.array_equal(a['perc'], b['perc']), f'SP{n} perc differs'
    # exfiltration actually reached the soil (guard against a vacuous test)
    iexf = T.INDEX_MM['iEXFg']
    assert np.any(np.asarray(out_dis[0]['MM'])[:, iexf] > 0)


def test_disv_mass_balance_closes():
    out, _, _ = _march('disv', nper=3, exf=1.0)
    imb, imbs = T.INDEX_MM['iMB'], T.INDEX_MM['iMBsurf']
    for n, o in enumerate(out):
        mm = np.asarray(o['MM'])
        assert np.max(np.abs(mm[:, imb])) < 1e-6, f'SP{n} soil MB'
        assert np.max(np.abs(mm[:, imbs])) < 1e-6, f'SP{n} surface MB'


def test_cell_order_contract_preserved():
    """The coupler asserts MM cell order == MF6 surf_cells order; the DISV
    node (icell2d) must stay i*ncol+j so both agree."""
    cMF = T._FakeMF(nper=1, perlen=[1])
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    cells = mm.build_cell_list(cMF)
    for cid, i, j, node in cells:
        assert node == i * cMF.ncol + j
        assert cells[cid][0] == cid
