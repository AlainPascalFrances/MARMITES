# -*- coding: utf-8 -*-
"""WP6.1 -- the streams and the ponds in the water-balance Sankey.

SFR and LAK exchange with the aquifer directly, both ways. The aquifer pass
reads them from the cell budget (records 'SFR' and 'LAK', + into the
aquifer, as MF6 writes them -- checked on a toy: the lake's seepage +8.5, a
gaining stream -6.8), every layer block of the Sankey draws them as signed
arms, and they enter the layer's closure. Missing, a losing stream's
seepage and a gaining stream's baseflow read as aquifer storage change.
"""
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


PP = _load('marmites_postprocess_sw',
           os.path.join(TRUNK, 'ppMF6', 'marmites_postprocess.py'))
IDX = _load('marmites_indices_sw', os.path.join(TRUNK, 'marmites_indices.py'))

NPER, NLAY, NROW, NCOL = 60, 2, 1, 4


def _model():
    cMF = SimpleNamespace(nlay=NLAY, nrow=NROW, ncol=NCOL,
                          delr=np.full(NCOL, 50.0), delc=np.full(NROW, 50.0),
                          ibound=np.ones((NLAY, NROW, NCOL), dtype=int),
                          lenuni=2, modelname='toy')
    ctx = SimpleNamespace(cells=[(k, 0, j, j) for j in range(NCOL)
                                 for k in [j]],
                          index=dict(IDX.INDEX_MM),
                          index_S=dict(IDX.INDEX_MM_SOIL), geom=None)
    wb = np.zeros((NPER, len(IDX.INDEX_MM)))
    ix = IDX.INDEX_MM
    wb[:, ix['iP']] = 2.0
    wb[:, ix['iEi']] = 0.2
    wb[:, ix['iPe']] = 1.8
    wb[:, ix['iI']] = 1.5
    wb[:, ix['iRo']] = 0.3
    wb[:, ix['iperc']] = 0.6
    wb[:, ix['iETsoil']] = 0.9
    wb[:, ix['iEg']] = 0.05
    wb[:, ix['iTg']] = 0.05
    wb[:, ix['iETg']] = 0.1
    res = {'perc': np.zeros((NPER, NCOL)), 'wb_ts': wb}
    return cMF, ctx, res


def _agg():
    """m3/d per (nper, target, layer), as _aquifer_pass returns it."""
    z = np.zeros((NPER, 1, NLAY))
    agg = {k: z.copy() for k, _t, _p in PP._AQ_RECORDS}
    agg['FLF'] = z.copy()
    agg['UZF-GWRCH'][:, 0, 0] = 5.0
    agg['SFR'][:, 0, 0] = -4.0         # a gaining stream: baseflow OUT
    agg['SFR'][:, 0, 1] = 1.0          # ... losing in layer 2
    agg['LAK'][:, 0, 0] = 2.0          # a pond seeping IN
    agg['GHB'][:, 0, 1] = -0.5
    return agg


def test_the_aquifer_pass_reads_the_streams_ponds_and_ghb():
    keys = [k for k, _t, _p in PP._AQ_RECORDS]
    for k in ('SFR', 'LAK', 'GHB'):
        assert k in keys


def test_the_aquifer_digest_is_keyed_on_the_records_read():
    """A digest cached before SFR/LAK were read must not be reused."""
    src = open(os.path.join(TRUNK, 'ppMF6', 'marmites_postprocess.py'),
               encoding='utf-8').read()
    body = src[src.index('def _aquifer_pass'):]
    body = body[:body.index('\ndef ')]
    assert "'+'.join(k for k, _t, _p in _AQ_RECORDS)" in body


def test_the_layer_fluxes_carry_the_streams_and_ponds_signed():
    cMF, ctx, res = _model()
    out, _n, _d = PP._aquifer_layer_fluxes('', 'toy', cMF, ctx, res,
                                           agg=_agg(), target=0)
    area = NROW * NCOL * 2500.0
    to_mm = 1000.0 / area
    assert out['iSFR_1'] == pytest.approx(np.full(NPER, -4.0 * to_mm))
    assert out['iSFR_2'] == pytest.approx(np.full(NPER, 1.0 * to_mm))
    assert out['iLAK_1'] == pytest.approx(np.full(NPER, 2.0 * to_mm))
    assert out['iGHB_2'] == pytest.approx(np.full(NPER, -0.5 * to_mm))
    assert PP._ghb_cells(out, NLAY) == [0, 1]


def test_the_sankey_draws_them_and_closes_the_layers(tmp_path):
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('agg')
    MMplot = PP._mmplot()
    if MMplot is None:
        pytest.skip('MARMITESplot not importable')
    cMF, ctx, res = _model()
    aq, ncell, drn = PP._aquifer_layer_fluxes('', 'toy', cMF, ctx, res,
                                              agg=_agg(), target=0)
    mmsv = np.zeros((NPER, 2, len(IDX.INDEX_MM_SOIL)))
    mmsv[:, 0, IDX.INDEX_MM_SOIL['iEsoil']] = 0.4
    mmsv[:, 0, IDX.INDEX_MM_SOIL['iTsoil']] = 0.5
    flx, fi = PP._assemble_flx(ctx.index, ctx.index_S, res['wb_ts'], mmsv,
                               aq, NPER)
    assert 'iSFR_1' in fi and 'iLAK_1' in fi
    import pandas as pd
    DATE = np.asarray(matplotlib.dates.date2num(
        pd.date_range('2008-06-01', periods=NPER, freq='D')), float)
    HY, years = PP._hydro_year_index(DATE, 10)
    smf = PP._SankeyMF(cMF, ncell, drn, DATE,
                       ghbcells=PP._ghb_cells(aq, NLAY))
    written = PP._render_sankey_inner(MMplot, str(tmp_path), DATE, flx, fi,
                                      HY, years, smf, ncell, [1, 1], 'toy',
                                      'toy', 0.0, False)
    assert written, 'no Sankey page was written'
