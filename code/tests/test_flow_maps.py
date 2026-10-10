# -*- coding: utf-8 -*-
"""The maps added on 2026-10-10 (user): the streamflow along the network,
the exchange with the streams and the ponds, and the flow across the UPPER
face of each layer -- the companion of FLF, so the valley's upward flow can
be followed out of the aquifer."""
import os
import shutil
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import marmites_postprocess as PP                    # noqa: E402


def test_the_upper_face_of_layer_1_is_the_net_exchange_at_the_top():
    """+ into the aquifer: recharge in, seepage / ET / a gaining stream out;
    layer 2's upper face is layer 1's lower face, plus what is attached to
    layer 2 itself (an outcrop cell, where layer 1 is absent)."""
    nlay, n = 2, 3
    z = np.zeros((nlay, n))
    maps = {'UZF-GWRCH': z + [[2.0, 0.0, 1.0], [0.0, 0.0, 0.5]],
            'DRN_SEEP': z + [[-0.5, 0.0, 0.0], [0.0, 0.0, 0.0]],
            'EVT_EG': z + [[-0.2, -0.1, 0.0], [0.0, 0.0, -0.1]],
            'SFR': z + [[0.0, -3.0, 0.0], [0.0, 0.0, 0.0]],
            'DRN': z + [[-9.0, -9.0, -9.0], [-9.0, -9.0, -9.0]],   # lateral
            'WEL': z + [[-7.0, 0.0, 0.0], [0.0, 0.0, 0.0]],       # pumping
            'FLF': z + [[1.0, -2.5, 0.0], [0.0, 0.0, 0.0]]}
    fuf = PP.upper_face_flow(maps, nlay)
    assert np.allclose(fuf[0], [1.3, -3.1, 1.0])
    assert np.allclose(fuf[1], [1.0, -2.5, 0.4])
    assert PP.upper_face_flow({'DRN': z}, nlay) is None


def _mf6_exe():
    import mm_paths
    if os.path.isfile(mm_paths.MF6_EXE):
        return mm_paths.MF6_EXE
    return shutil.which('mf6')


@pytest.fixture(scope='module')
def stream_run(tmp_path_factory):
    """Three reaches in a row over a two-layer aquifer, fed 100 m3/d, with a
    streambed so tight the exchange is negligible: each reach hands on the
    100 m3/d, and the last one hands it out of the model."""
    flopy = pytest.importorskip('flopy')
    exe = _mf6_exe()
    if not exe:
        pytest.skip('no mf6')
    ws = str(tmp_path_factory.mktemp('stream'))
    name = 'strm'
    nlay, nrow, ncol = 2, 3, 4
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 2, 1.0)] * 2,
                          time_units='DAYS')
    flopy.mf6.ModflowIms(sim, complexity='SIMPLE',
                         linear_acceleration='BICGSTAB')
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True,
                               newtonoptions='NEWTON')
    flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=20.0,
                            delc=20.0, top=100.0, botm=[90.0, 50.0])
    flopy.mf6.ModflowGwfic(gwf, strt=95.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0, k33=0.1)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.1,
                            transient={0: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(1, 0, 0), 95.0]]})
    pkg = [(r, (0, 1, r + 1), 20.0, 1.0, 0.01, 99.0 - r, 0.5, 1e-9, 0.035,
            1 if r in (0, 2) else 2, 1.0, 0) for r in range(3)]
    con = [(0, -1), (1, 0, -2), (2, 1)]
    flopy.mf6.ModflowGwfsfr(gwf, nreaches=3, packagedata=pkg,
                            connectiondata=con,
                            perioddata={0: [(0, 'INFLOW', 100.0)]},
                            save_flows=True,
                            budget_filerecord='%s.sfr.cbc' % name)
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord='%s.cbc' % name,
                           saverecord=[('BUDGET', 'ALL')])
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    if not ok:
        pytest.skip('the stream model did not run')
    return ws, name, nlay, nrow, ncol


def test_each_reach_hands_on_its_flow(stream_run):
    ws, name, nlay, nrow, ncol = stream_run
    got = PP.sfr_reach_flows(ws, name)
    assert np.allclose(got['q'], 100.0, rtol=1e-3), got['q']
    # the GWF user node of each reach: layer 1, row 2, columns 2..4
    assert list(got['node']) == [1 * ncol + c + 1 for c in (1, 2, 3)]


def test_no_sfr_budget_no_streamflow(tmp_path):
    assert PP.sfr_reach_flows(str(tmp_path), 'none') is None


def test_layer_2s_upper_face_is_layer_1s_lower_face(stream_run):
    """On a real MF6 budget: nothing attached to layer 2 but a CHD, so what
    crosses its top is exactly what leaves layer 1 at the bottom."""
    ws, name, nlay, nrow, ncol = stream_run
    cMF = SimpleNamespace(nlay=nlay, nrow=nrow, ncol=ncol)
    maps = PP._aquifer_map_pass(ws, name, cMF, 2, cache_dir=None,
                                verbose=False)
    assert 'FLF' in maps and 'SFR' in maps
    fuf = PP.upper_face_flow(maps, nlay)
    flf = np.asarray(maps['FLF']).reshape(nlay, -1)
    assert np.allclose(fuf[1], flf[0])
    assert np.allclose(fuf[0], np.asarray(maps['SFR']).reshape(nlay, -1)[0])
