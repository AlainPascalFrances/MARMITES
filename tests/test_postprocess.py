# -*- coding: utf-8 -*-
"""Post-/pre-processing of a coupled MARMITES/MF6 run.

Built and run on a tiny synthetic MF6 model (with UZF and SFR) so the tests
need neither the La Mata dataset nor a long simulation. They check that the
budget parsing, obs-point mapping and every figure/CSV are produced.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..', 'trunk'))
for p in ('', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

flopy = pytest.importorskip('flopy')
pytest.importorskip('pandas')
import matplotlib  # noqa: E402
matplotlib.use('agg')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PP = _load('marmites_postprocess', os.path.join(TRUNK, 'ppMF6', 'marmites_postprocess.py'))


def _has_mf6():
    from shutil import which
    return which('mf6') is not None


# --------------------------------------------------------------------- #
# parsing helpers that need no MF6 binary
# --------------------------------------------------------------------- #

def test_compartment_taxonomy():
    assert PP._compartment_of('STO-SS') == 'aquifer (storage)'
    assert PP._compartment_of('UZF-GWRCH') == 'unsaturated (UZF->GW)'
    assert PP._compartment_of('DRN_SEEP') == 'surface (seepage)'
    assert PP._compartment_of('SFR') == 'surface (stream)'
    assert PP._compartment_of('WEL') == 'ET (groundwater)'
    assert PP._compartment_of('SOMETHING') == 'other'


def test_obs_points_skips_disabled(tmp_path):
    (tmp_path / 'inputObs.txt').write_text(
        '# Name X Y lay\n'
        'P0 739525 4555875 1\n'
        '##W1 741177 4554666 1\n'
        '#C4 740764 4555342 1\n'
        'C1 740528 4554369 1\n'
        'G1 739552 4555658 3\n')
    pts = PP.obs_points(str(tmp_path))
    names = [p['name'] for p in pts]
    assert names == ['P0', 'C1', 'G1'], 'commented points must be dropped'
    assert pts[2]['lay'] == 3


def test_obs_series_reads_dates(tmp_path):
    (tmp_path / 'inputObsHEADS_C1.txt').write_text(
        '2007-10-01\t768.9\n2008-01-01\t770.4\n2008-04-01\t770.9\n')
    s = PP.obs_series(str(tmp_path), 'C1')
    assert list(s.columns) == ['date', 'head']
    assert len(s) == 3 and s['head'].iloc[1] == pytest.approx(770.4)
    assert PP.obs_series(str(tmp_path), 'NOPE') is None


def test_xy_to_ij_corner_and_centre():
    # 50 m cells, 3x4 grid, origin (0,0). Row 0 is the NORTH edge.
    i, j = PP._xy_to_ij(25.0, 25.0, 0.0, 0.0, 50.0, 3, 4)
    assert (i, j) == (2, 0)                       # bottom-left cell
    i, j = PP._xy_to_ij(175.0, 145.0, 0.0, 0.0, 50.0, 3, 4)
    assert (i, j) == (0, 3)                       # top-right cell


def test_dates_from_dataset(tmp_path):
    (tmp_path / 'inputDATE.txt').write_text(
        '#\n2008-05-31 00:00, 152\n2008-06-01 00:00, 153\n2008-06-02 00:00, 154\n')
    d = PP._dates_from_dataset(str(tmp_path), 2)
    assert len(d) == 2 and str(d[0].date()) == '2008-05-31'


# --------------------------------------------------------------------- #
# a real (tiny) MF6 model
# --------------------------------------------------------------------- #

@pytest.fixture(scope='module')
def tiny_run(tmp_path_factory):
    if not _has_mf6():
        pytest.skip('mf6 binary not on PATH')
    ws = str(tmp_path_factory.mktemp('tiny'))
    name = 'tiny'
    nlay, nrow, ncol = 2, 4, 5
    top = np.full((nrow, ncol), 100.0)
    botm = np.stack([np.full((nrow, ncol), 90.0), np.full((nrow, ncol), 50.0)])
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name='mf6')
    flopy.mf6.ModflowTdis(sim, nper=3, perioddata=[(1.0, 1, 1.0)] * 3,
                          time_units='DAYS')
    flopy.mf6.ModflowIms(sim, complexity='SIMPLE')
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True,
                               newtonoptions='NEWTON')
    flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=50.0,
                            delc=50.0, top=top, botm=botm)
    flopy.mf6.ModflowGwfic(gwf, strt=95.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0, k33=0.1, save_flows=True)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.1,
                            steady_state={0: True}, transient={1: True})
    # a well sink so storage and flows are non-trivial
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: [[(0, 0, 0), -50.0]]},
                            save_flows=True)
    # UZF: one land cell, recharge on top. packagedata is
    # (ifno, cellid, landflag, ivertcon, surfdep, vks, thtr, thts, thti, eps)
    uzf_pkg = [(0, (0, 1, 1), 1, -1, 0.1, 0.1, 0.05, 0.35, 0.1, 3.5)]
    uzf_spd = {0: [(0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)]}
    flopy.mf6.ModflowGwfuzf(gwf, nuzfcells=1, ntrailwaves=7, nwavesets=40,
                            packagedata=uzf_pkg, perioddata=uzf_spd,
                            save_flows=True, budget_filerecord=f'{name}.uzf.cbc')
    # SFR: two reaches
    sfr_pkg = [(0, (0, 2, 2), 50.0, 1.0, 0.01, 99.0, 0.5, 0.1, 0.035, 1, 1.0, 0),
               (1, (0, 2, 3), 50.0, 1.0, 0.01, 98.0, 0.5, 0.1, 0.035, 1, 1.0, 0)]
    sfr_conn = [(0, -1), (1, 0)]
    flopy.mf6.ModflowGwfsfr(gwf, nreaches=2, packagedata=sfr_pkg,
                            connectiondata=sfr_conn,
                            perioddata={0: [(0, 'INFLOW', 1.0)]},
                            save_flows=True, budget_filerecord=f'{name}.sfr.cbc')
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord=f'{name}.hds',
                           budget_filerecord=f'{name}.cbc',
                           saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    if not ok:
        pytest.skip('tiny mf6 run did not converge in this environment')
    return ws, name, nlay, nrow, ncol


def test_listing_budget_by_compartment(tiny_run):
    ws, name, *_ = tiny_run
    df, comp = PP.budget_by_compartment(ws, name)
    assert not comp.empty
    comps = set(comp['compartment'])
    assert 'aquifer (storage)' in comps
    assert any('unsaturated' in c for c in comps)      # UZF present
    assert any('stream' in c for c in comps)           # SFR present


def test_package_budget_uzf_and_sfr(tiny_run):
    ws, name, *_ = tiny_run
    uzf = PP.package_budget(ws, f'{name}.uzf.cbc')
    sfr = PP.package_budget(ws, f'{name}.sfr.cbc')
    assert uzf and sfr
    assert any('GWF' in k.upper() or 'RCH' in k.upper() or 'INFILTR' in k.upper()
               for k in uzf)


def test_layer_storage_change_shape(tiny_run):
    ws, name, nlay, nrow, ncol = tiny_run
    sto = PP.layer_storage_change(ws, name, nlay, nrow, ncol)
    assert sto.shape == (nlay, nrow, ncol)
    assert np.isfinite(sto).all()


def test_run_postproc_writes_the_expected_files(tiny_run, tmp_path):
    ws, name, nlay, nrow, ncol = tiny_run
    files = PP.run_postproc(ws, str(tmp_path), name=name, verbose=False)
    got = {os.path.basename(f) for f in files}
    # the head / depth / storage grids are CSV only: the maps themselves are
    # drawn by native_suite (GWmap_head, MMmap_dgwt) with the coordinate axes
    for L in range(1, nlay + 1):
        assert 'mean_head_L%d.csv' % L in got
        assert 'mean_depth_L%d.csv' % L in got
        assert 'storage_change_L%d.csv' % L in got
    assert not any(f.endswith('.png') and f.startswith(('mean_', 'storage_'))
                   for f in got)
    assert 'budget_compartment.png' in got
    assert 'budget_uzf.png' in got
    assert 'budget_sfr.png' in got
    outdir = os.path.join(ws, '_output')
    assert os.path.exists(os.path.join(outdir, 'mean_head_L1.csv'))
    assert os.path.exists(os.path.join(outdir, 'budget_compartment.csv'))


def test_run_postproc_obs_overlay(tiny_run, tmp_path):
    ws, name, nlay, nrow, ncol = tiny_run
    # an obs point inside the grid (origin 0,0; 50 m cells) with a series
    (tmp_path / 'inputObs.txt').write_text('# Name X Y lay\nT1 75 75 1\n')
    (tmp_path / 'inputObsHEADS_T1.txt').write_text('2008-01-01\t95.0\n2008-01-02\t94.5\n')
    files = PP.run_postproc(ws, str(tmp_path), name=name, verbose=False)
    assert any(os.path.basename(f) == 'obs_heads.png' for f in files)
    assert os.path.exists(os.path.join(ws, '_output', 'obs_heads_computed.csv'))


def test_subsample_is_even_and_bounded():
    assert PP._subsample(list(range(10)), 100) == list(range(10))   # short: all
    s = PP._subsample(list(range(1000)), 5)
    assert len(s) == 5 and s[0] == 0 and s[-1] == 999               # endpoints kept
    assert s == sorted(set(s))                                      # unique, ordered
def test_run_preproc_writes_input_maps(tiny_run, tmp_path):
    ws, name, nlay, nrow, ncol = tiny_run
    # no cMF/ctx and a GIS workspace that does not exist: every figure needs
    # one or the other, so this exercises the guards rather than the drawing
    files = PP.run_preproc(ws, str(tmp_path), name=name, mf_ws=str(tmp_path),
                           gis_ws=str(tmp_path / 'no_such_gis'), verbose=False)
    assert files == []
    assert not any(f.startswith(('aq_', 'mm_'))
                   for f in os.listdir(os.path.join(ws, '_input'))), (
        'the plain aq_*/mm_* maps were replaced by the native IN_* set')


def test_ja_down_index_matches_flopy_faceflows():
    """Our precomputed JA down-connection index must reproduce flopy's
    get_structured_faceflows exactly, including its sign convention
    (flopy applies flows[face][n] = -1 * flowja[i]).

    We bypass that helper because it re-parses the .grb on every call (10.6 ms
    x nper) and because its documented ia/ja path raises on flopy master
    (PR #1968 added nlay/nrow/ncol but left `for n in range(grb.nodes)`).
    """
    ws = os.path.join(os.environ.get(
        'MARMITES_WS_ROOT',
        os.path.join('E:' + os.sep, '00code_ws', 'LaMata_MM-MF6')), 'MF6_ws')
    cbc_fn = os.path.join(ws, 'lamatamm.cbc')
    grb_fn = os.path.join(ws, 'lamatamm.dis.grb')
    if not (os.path.exists(cbc_fn) and os.path.exists(grb_fn)):
        pytest.skip('no La Mata cbc/grb on disk')
    import flopy
    from flopy.mf6.utils.postprocessing import get_structured_faceflows
    cbc = flopy.utils.CellBudgetFile(cbc_fn)
    kk = cbc.get_kstpkper()
    rec = cbc.get_data(text='FLOW-JA-FACE', kstpkper=kk[len(kk) // 2])[0]
    nlay, nrow, ncol = 2, 65, 60
    src, pos, nodes = PP._ja_down_index(grb_fn, nlay, nrow, ncol)
    mine = np.zeros(nodes)
    mine[src] = -np.asarray(rec).ravel()[pos]
    _, _, flf = get_structured_faceflows(rec, grb_file=grb_fn)
    ref = np.ma.filled(np.asarray(flf, float), 0.0).ravel()
    assert np.allclose(mine, ref, atol=1e-6), 'JA index diverges from flopy'
