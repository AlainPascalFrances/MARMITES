# -*- coding: utf-8 -*-
"""Phase-4 tests: the refined-grid (quadtree) input layout.

A genuinely refined grid has no (row, column), so MARMITES uses the
degenerate layout agreed in the design: cells are ``(cid, cid, 0, icell2d)``
and every spatial input is an ``(ncell, 1)`` column vector, which makes the
kernel's legacy ``grid[i, j]`` lookups resolve to ``grid[cid, 0]``.

The test that matters: feeding the SAME physical model through the 2-D
layout and through the degenerate layout must give identical results. That
validates the refined-grid data path without needing the gridgen binary
(refinement then only changes cell areas and count, both already covered).
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


G = _load('marmites_grid_ref', os.path.join(TRUNK, 'marmites_grid.py'))
GG = _load('marmites_gridgen_ref', os.path.join(TRUNK, 'marmites_gridgen.py'))
T = _load('t_runmmsoil_ref', os.path.join(HERE, 'test_runmmsoil.py'))


class _DegenProcess:
    """cPROCESS stand-in returning (nlay, ncell, 1) for the degenerate layout."""

    def __init__(self, sy, nlay, ncell):
        self._sy, self.nlay, self.ncell = sy, nlay, ncell

    def float2array(self, _):
        return np.full((self.nlay, self.ncell, 1), self._sy)


class _DegenMF:
    """Degenerate-layout twin of tests._FakeMF (nrow := ncell, ncol := 1)."""

    def __init__(self, src, ncell):
        self.nper, self.perlen, self.nstp = src.nper, src.perlen, src.nstp
        self.hnoflo, self.hdry = src.hnoflo, src.hdry
        self.uzf_yn, self.wel_yn = src.uzf_yn, src.wel_yn
        self.nlay = src.nlay
        self.nrow, self.ncol = ncell, 1
        self.delr, self.delc = [1.0], [1.0] * ncell   # unused: geom supplies width
        self.outcropL = np.ones((ncell, 1), dtype=int)
        self.sy_actual = src.sy_actual
        self.cPROCESS = _DegenProcess(src.sy_actual, src.nlay, ncell)


def _ctx_2d(cMF, inp, mm, grid='disv'):
    cells = mm.build_cell_list(cMF)
    geom = G.geometry_for(cMF, cells, grid=grid)
    ctx = mm.build_context(cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'], inp['_Sm'],
                           inp['_Sfc'], inp['_Sr'], inp['_slprop'], inp['_Ssoil_ini'],
                           inp['botm_l0'], inp['_Ks'], inp['gridSOIL'], inp['gridSOILthick'],
                           inp['TopSoil'], inp['gridMETEO'], T.INDEX_MM, T.INDEX_MM_S,
                           inp['P_veg_zoneSP'],
                           inp['Eo_zonesSP'], inp['PT_veg_zonesSP'], inp['Pe_veg_zonesSP'],
                           inp['PE_zonesSP'], inp['gridVEGarea'], inp['LAI_veg_zonesSP'],
                           inp['Zr'], inp['kTg_min'], inp['kTg_max'], inp['kT_f'], inp['kT_s'],
                           inp['NVEG'], 1000.0, 0, None, None, None, None, None,
                           None, None, None, None, None, geom=geom)
    return ctx, mm.init_state(ctx), cells


def _ctx_degenerate(cMF, inp, mm, cells2d):
    """Same model, column-vector inputs, cells (cid, cid, 0, icell2d)."""
    nodes = np.array([c[3] for c in cells2d], dtype=int)      # keep icell2d
    ncell = len(cells2d)
    cells = GG.refined_cell_list(nodes)
    assert cells[0] == (0, 0, 0, int(nodes[0]))
    dmf = _DegenMF(cMF, ncell)
    # column-vector inputs sampled at the same physical cells
    ii = np.array([c[1] for c in cells2d], dtype=int)
    jj = np.array([c[2] for c in cells2d], dtype=int)
    col = lambda a: np.asarray(a)[ii, jj].reshape(-1, 1)      # noqa: E731
    veg = np.stack([inp['gridVEGarea'][v][ii, jj].reshape(-1, 1)
                    for v in range(inp['NVEG'])])
    # geometry over the SAME polygons -> identical areas/widths
    verts, cell2d, _ = G.disv_from_structured(cMF.delr, cMF.delc,
                                              getattr(cMF, 'xllcorner', 0.0),
                                              getattr(cMF, 'yllcorner', 0.0))
    geom = G.VertexGeometry.from_vertices(verts, cell2d, nodes, nlay=cMF.nlay)
    ctx = mm.build_context(dmf, cells, inp['_nsl'], inp['_nslmax'], inp['_st'], inp['_Sm'],
                           inp['_Sfc'], inp['_Sr'], inp['_slprop'], inp['_Ssoil_ini'],
                           col(inp['botm_l0']), inp['_Ks'], col(inp['gridSOIL']),
                           col(inp['gridSOILthick']), col(inp['TopSoil']),
                           col(inp['gridMETEO']), T.INDEX_MM, T.INDEX_MM_S,

                           inp['P_veg_zoneSP'], inp['Eo_zonesSP'], inp['PT_veg_zonesSP'],
                           inp['Pe_veg_zonesSP'], inp['PE_zonesSP'], veg,
                           inp['LAI_veg_zonesSP'], inp['Zr'], inp['kTg_min'], inp['kTg_max'],
                           inp['kT_f'], inp['kT_s'], inp['NVEG'], 1000.0, 0,
                           None, None, None, None, None, None, None, None, None, None,
                           geom=geom)
    return ctx, mm.init_state(ctx), cells


def _march(ctx, state, mm, nper, heads=699.0, exf=0.0):
    outs = []
    for n in range(nper):
        h = np.full(ctx.ncell, heads)
        e = np.full(ctx.ncell, exf)
        outs.append(mm.step(ctx, n, n, h, e, state))
    return outs


def _pair(nper=3, exf=0.0):
    cMF = T._FakeMF(nper=nper, perlen=[1] * nper)
    cMF.xllcorner = cMF.yllcorner = 0.0
    inp = T._build_inputs(cMF)
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    c2, s2, cells2d = _ctx_2d(cMF, inp, mm)
    cd, sd, _ = _ctx_degenerate(cMF, inp, mm, cells2d)
    return (_march(c2, s2, mm, nper, exf=exf),
            _march(cd, sd, mm, nper, exf=exf), s2, sd)


# --------------------------------------------------------------------- #

def test_refined_cell_list_shape():
    cells = GG.refined_cell_list([7, 3, 11])
    assert cells == [(0, 0, 0, 7), (1, 1, 0, 3), (2, 2, 0, 11)]
    # i == cid and j == 0 -> grid[i, j] resolves to a column vector
    for cid, i, j, _node in cells:
        assert i == cid and j == 0


def test_degenerate_layout_matches_2d_exactly():
    out2, outd, s2, sd = _pair(nper=3)
    for n, (a, b) in enumerate(zip(out2, outd)):
        for key in ('MM', 'MM_S', 'perc', 'etg'):
            assert np.array_equal(np.asarray(a[key]), np.asarray(b[key])), \
                f'SP{n} differs in {key}'
    assert np.array_equal(s2.Ssoil_ini, sd.Ssoil_ini)


def test_degenerate_layout_matches_with_exfiltration():
    out2, outd, _, _ = _pair(nper=3, exf=2.0)
    for n, (a, b) in enumerate(zip(out2, outd)):
        assert np.array_equal(a['MM'], b['MM']), f'SP{n} MM'
        assert np.array_equal(a['perc'], b['perc']), f'SP{n} perc'
    iexf = T.INDEX_MM['iEXFg']
    assert np.any(np.asarray(out2[0]['MM'])[:, iexf] > 0)


def test_column_vector_inputs_are_actually_used():
    """Guard against a vacuous test: perturbing one cell's soil thickness in
    the degenerate inputs must change only that cell's result."""
    nper = 1
    cMF = T._FakeMF(nper=nper, perlen=[1])
    cMF.xllcorner = cMF.yllcorner = 0.0
    inp = T._build_inputs(cMF)
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    _c2, _s2, cells2d = _ctx_2d(cMF, inp, mm)
    cd, sd, _ = _ctx_degenerate(cMF, inp, mm, cells2d)
    base = mm.step(cd, 0, 0, np.full(cd.ncell, 699.0), np.zeros(cd.ncell), sd)
    # perturb cell 0 only
    cd2, sd2, _ = _ctx_degenerate(cMF, inp, mm, cells2d)
    cd2.gridSOILthick[0, 0] *= 1.5
    pert = mm.step(cd2, 0, 0, np.full(cd2.ncell, 699.0), np.zeros(cd2.ncell), sd2)
    b = np.asarray(base['MM'])
    p = np.asarray(pert['MM'])
    assert not np.array_equal(b[0], p[0]), 'perturbed cell must change'
    assert np.array_equal(b[1:], p[1:]), 'other cells must be unaffected'


def test_gridgen_module_requires_executable():
    """build_quadtree fails clearly when the gridgen binary is absent."""
    class _C:
        nlay, delr, delc = 1, [10.0], [10.0]
        top, botm = np.zeros((1, 1)), np.zeros((1, 1, 1))
        xllcorner = yllcorner = 0.0
    try:
        GG.build_quadtree(_C(), '/definitely/not/here/gridgen', '/tmp/gg')
    except FileNotFoundError as e:
        assert 'gridgen' in str(e)
    else:
        raise AssertionError('missing gridgen executable must raise')
