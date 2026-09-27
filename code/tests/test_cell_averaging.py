# -*- coding: utf-8 -*-
"""[solver] cell_averaging -- NPF's alternative_cell_averaging, asked on the
Run panel. harmonic is MF6's own default and writes nothing; amt-hmk
(arithmetic-mean saturated thickness x harmonic-mean K) is what the CdL
model runs with: it keeps the conductance up as a cell dewaters, and
removed a float-overflow crash there."""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

import marmites_config as mcfg  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_default_is_mf6s_own_harmonic_mean():
    cfg = mcfg.RunConfig.from_dict({})
    assert cfg.solver.cell_averaging == 'harmonic'
    assert 'amt-hmk' in mcfg.CELL_AVERAGING


def test_an_unknown_averaging_is_refused():
    cfg = mcfg.RunConfig.from_dict({})
    cfg.solver.cell_averaging = 'geometric'
    with pytest.raises(mcfg.ConfigError) as e:
        cfg.validate()
    assert 'cell_averaging' in str(e.value)


@pytest.mark.parametrize('avg, written', [('harmonic', None),
                                          ('amt-hmk', 'AMT-HMK')])
def test_the_npf_carries_the_averaging(tmp_path, avg, written):
    flopy = pytest.importorskip('flopy')
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s1L.ini')):
        pytest.skip('La Mata dataset not present')
    import matplotlib
    matplotlib.use('agg')
    import dataclasses
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    mf6mod = _load('marmites_mf6_avg', os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'))
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s1L.ini',
                   xllcorner=739300.0, yllcorner=4553050.0)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]
    b = mf6mod.clsMF6(c, top=np.asarray(c.elev, float),
                      botm=np.asarray(c.botm, float), sim_ws=str(tmp_path),
                      daily=True)
    b.verbose = False
    sv = dataclasses.asdict(mcfg.Solver())
    sv['cell_averaging'] = avg
    b.solver = sv
    b.build()
    b.write()
    npf = open(os.path.join(tmp_path, '%s.npf' % c.modelname.lower())).read()
    if written is None:
        assert 'ALTERNATIVE_CELL_AVERAGING' not in npf.upper()
    else:
        assert 'ALTERNATIVE_CELL_AVERAGING  %s' % written in npf.upper() \
            or 'ALTERNATIVE_CELL_AVERAGING %s' % written in npf.upper(), npf[:300]
