# -*- coding: utf-8 -*-
"""The MODFLOW 6 solver is asked on the Run panel, not taken from NWT.

The run of 2026-09-23 stopped on a 4.09 % cumulative discrepancy: 35-65 m3/d
unaccounted in every step, from the steady period on. The IMS tolerances
were the legacy NWT ini's -- HEADTOL 0.05 m as outer_dvclose, the inner ones
MF6's COMPLEX defaults -- and a 5 cm head error with Sy 0.01 is ~0.15 m3/d
of storage per average cell. It had been there all along, hidden while a
5.7e6 m3 UZF pulse dominated the budget.
"""

import dataclasses
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(CODE, '..', 'example', 'LaMata'))
for _p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6',
           'app'):
    _d = os.path.join(CODE, _p)
    if _d not in sys.path:
        sys.path.insert(0, _d)

import marmites_config as mcfg                                 # noqa: E402

APPROVED = {'complexity': 'complex', 'outer_dvclose': 0.001,
            'outer_maximum': 500, 'inner_dvclose': 0.0001,
            'inner_rclose': 0.01}


# ------------------------------------------------------ the configuration
def test_the_defaults_are_the_approved_values():
    assert dataclasses.asdict(mcfg.RunConfig.from_dict({}).solver) == APPROVED


def test_lamata_asks_for_them():
    """The section is in the file -- not its current values: lamata.toml is
    a live working file, edited on the panels between test runs."""
    path = os.path.join(CODE, 'configs', 'lamata.toml')
    assert '[solver]' in open(path, encoding='utf-8').read()
    cfg = mcfg.load_run_config(path)
    assert set(dataclasses.asdict(cfg.solver)) == set(APPROVED)


@pytest.mark.parametrize('bad, needle', [
    ({'complexity': 'hard'}, 'solver.complexity'),
    ({'outer_dvclose': 0.0}, 'solver.outer_dvclose'),
    ({'inner_rclose': -1.0}, 'solver.inner_rclose'),
    ({'outer_maximum': 0}, 'solver.outer_maximum'),
    ({'inner_dvclose': 0.01, 'outer_dvclose': 0.001}, 'inner_dvclose'),
])
def test_an_unusable_solver_is_refused(bad, needle):
    with pytest.raises(mcfg.ConfigError) as exc:
        mcfg.RunConfig.from_dict({'solver': bad})
    assert needle in str(exc.value)


# --------------------------------------------------------------- the build
def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope='module')
def cmf():
    pytest.importorskip('flopy')
    if not os.path.exists(os.path.join(DS, 'MF_ws',
                                       '__inputMF_flopy_v3_2s1L.ini')):
        pytest.skip('La Mata dataset not present')
    import matplotlib
    matplotlib.use('agg')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s1L.ini',
                   xllcorner=739300.0, yllcorner=4553050.0)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 2, [1, 1], [1, 1]
    return c


def _ims(cmf, tmp, solver=None):
    mf6mod = _load('marmites_mf6_solver',
                   os.path.join(CODE, 'ppMF6', 'marmites_mf6.py'))
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp), grid='dis')
    if solver is not None:
        b.solver = solver
    b.build()
    b.sim.write_simulation(silent=True)
    ims = [f for f in os.listdir(str(tmp)) if f.endswith('.ims')]
    txt = open(os.path.join(str(tmp), ims[0])).read().upper()
    # KEY -> VALUE, whatever flopy's column spacing
    vals = {}
    for line in txt.splitlines():
        p = line.split()
        if len(p) == 2 and not line.startswith('#'):
            vals[p[0]] = p[1]
    return b, txt, vals


def test_the_build_writes_the_panel_values_not_the_nwt_headtol(cmf, tmp_path):
    assert float(cmf.headtol) == pytest.approx(0.05), 'the ini changed'
    b, txt, v = _ims(cmf, tmp_path)
    assert float(v['OUTER_DVCLOSE']) == pytest.approx(0.001), txt
    assert float(v['INNER_DVCLOSE']) == pytest.approx(1e-4), txt
    assert float(v['INNER_RCLOSE']) == pytest.approx(0.01), txt
    assert 'STRICT' not in txt and 'L2NORM' not in txt, 'an option crept in'
    assert v['COMPLEXITY'] == 'COMPLEX'
    assert b.outer_maximum == 500


def test_a_configured_solver_reaches_the_ims_file(cmf, tmp_path):
    sv = dict(APPROVED, outer_dvclose=0.002, outer_maximum=77,
              complexity='moderate')
    b, txt, v = _ims(cmf, tmp_path, solver=sv)
    assert float(v['OUTER_DVCLOSE']) == pytest.approx(0.002), txt
    assert v['OUTER_MAXIMUM'] == '77' and v['COMPLEXITY'] == 'MODERATE'
    assert b.outer_maximum == 77, 'the coupler caps its loop on this'


def test_the_run_passes_the_whole_section():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'),
               encoding='utf-8').read()
    assert 'b.solver = dataclasses.asdict(cfg.solver)' in src
