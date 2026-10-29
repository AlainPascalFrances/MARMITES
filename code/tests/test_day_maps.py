# -*- coding: utf-8 -*-
"""Maps of ONE day from the Results panel (user, 2026-10-29: "Make a panel
into results where you ask for a day (calendar) and a flux (list) to
produce this type of maps"): the run's map of a flux -- Rn first -- on the
day asked, from its own model outputs while they are still in its
workspace (ppMF6/marmites_daymap)."""
import os
import shutil
import sys
import time
from datetime import date, timedelta

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'app')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import marmites_config as mcfg                       # noqa: E402
import marmites_daymap as DM                         # noqa: E402
from lib import results as R                         # noqa: E402

NAME = 'dm'
NLAY, NROW, NCOL = 2, 2, 3
START = date(2010, 6, 1)
NPER = 3


def _mf6_exe():
    import mm_paths
    if os.path.isfile(mm_paths.MF6_EXE):
        return mm_paths.MF6_EXE
    return shutil.which('mf6')


@pytest.fixture(scope='module')
def run(tmp_path_factory):
    """Two layers of 2 x 3 cells, three days: UZF on layer 1 fed a
    different infiltration each day, the seepage drains and groundwater
    evaporation of a MARMITES run (their package names), a fixed head in
    layer 2 -- and a results folder of the run, with its configuration
    copy and its per-period dates."""
    flopy = pytest.importorskip('flopy')
    exe = _mf6_exe()
    if not exe:
        pytest.skip('no mf6')
    root = tmp_path_factory.mktemp('dayws')
    ws = str(root / 'MF6_ws')
    sim = flopy.mf6.MFSimulation(sim_name=NAME, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=NPER, perioddata=[(1.0, 1, 1.0)] * NPER,
                          time_units='DAYS')
    flopy.mf6.ModflowIms(sim, complexity='MODERATE')
    gwf = flopy.mf6.ModflowGwf(sim, modelname=NAME, save_flows=True,
                               newtonoptions='NEWTON')
    flopy.mf6.ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0,
                            delc=10.0, top=100.0, botm=[95.0, 80.0])
    flopy.mf6.ModflowGwfic(gwf, strt=99.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0, k33=0.5)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.1,
                            transient={0: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={
        0: [[(1, 1, 2), 98.0]]})
    cells = [(0, i, j) for i in range(NROW) for j in range(NCOL)]
    pkg = [(n, c, 1, -1, 0.1, 5.0, 0.05, 0.35, 0.1, 3.5)
           for n, c in enumerate(cells)]
    finf = (0.002, 0.03, 0.0)
    per = {k: [(n, finf[k] * (1 + n), 0.0, 1.0, 0.05, 0.0, 0.0, 0.0)
               for n in range(len(cells))] for k in range(NPER)}
    flopy.mf6.ModflowGwfuzf(gwf, nuzfcells=len(cells), packagedata=pkg,
                            perioddata=per, ntrailwaves=7, nwavesets=40)
    flopy.mf6.ModflowGwfdrn(gwf, pname='drn_seep', stress_period_data={
        0: [[(0, 0, 0), 99.2, 50.0]]})
    flopy.mf6.ModflowGwfevt(gwf, pname='evt_eg', stress_period_data={
        0: [[(0, 0, 1), 100.0, 0.002, 2.0]]})
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord='%s.cbc' % NAME,
                           head_filerecord='%s.hds' % NAME,
                           saverecord=[('BUDGET', 'ALL'), ('HEAD', 'ALL')])
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    if not ok:
        pytest.skip('the day-map model did not run')
    out = _results(root, ws, 'out_%s_try')
    return {'ws': ws, 'out': out, 'root': root}


def _results(root, ws, pattern, shift=-120.0):
    """A results folder stamped ``shift`` s from the model's outputs, its
    figures written 60 s after them."""
    t = os.path.getmtime(os.path.join(ws, '%s.cbc' % NAME))
    out = root / (pattern % time.strftime('%Y%m%d%H%M',
                                           time.localtime(t + shift)))
    (out / '_input').mkdir(parents=True)
    (out / '_output').mkdir()
    cfg = mcfg.RunConfig.from_dict({})
    cfg.paths.ws = ws
    cfg.write_toml(str(out / '_input' / 'resolved_config.toml'))
    days = [START + timedelta(days=k) for k in range(NPER)]
    p = out / '_output' / 'sp_run_timing.csv'
    p.write_text('date,period,seconds,iterations\n' + ''.join(
        '%s,%d,1.0,10\n' % (d.isoformat(), k + 1) for k, d in enumerate(days)))
    for f in (p, out / '_input' / 'resolved_config.toml'):
        os.utime(str(f), (t + shift + 180.0, t + shift + 180.0))
    return str(out)


def _cbc(ws, text, pak2=None, k=0):
    import flopy
    cbc = flopy.utils.CellBudgetFile(os.path.join(ws, '%s.cbc' % NAME))
    try:
        kk = cbc.get_kstpkper()
        d = cbc.get_data(text=text, paknam2=pak2, kstpkper=kk[k], full3D=True)
        return np.ma.filled(np.asarray(d[0], float), 0.0).reshape(NLAY, -1)
    finally:
        cbc.close()


# ---- which run, which days -------------------------------------------------

def test_the_outputs_are_this_runs_while_no_later_run_wrote_them(run):
    info = DM.outputs_of(run['out'])
    assert info['why'] == '', info['why']
    assert info['name'] == NAME and os.path.samefile(info['sim_ws'], run['ws'])
    assert info['dates'] == [START + timedelta(days=k) for k in range(NPER)]
    # a run whose figures were done 10 minutes before the outputs: a later
    # run wrote them
    old = _results(run['root'], run['ws'], 'out_%s_older', shift=-900.0)
    assert 'by a later run' in DM.outputs_of(old)['why']
    # a results folder made after them: it did not run the model
    new = _results(run['root'], run['ws'], 'out_%s_redraw', shift=600.0)
    assert 'older than this run' in DM.outputs_of(new)['why']


def test_why_a_day_cannot_be_drawn(tmp_path, run):
    assert 'no copy of its configuration' in DM.outputs_of(
        str(tmp_path / 'out_202610291200_x'))['why']
    assert DM.model_name(run['ws']) == NAME
    assert DM.model_name(str(tmp_path)) is None
    with pytest.raises(ValueError, match='not a day of this run'):
        DM.draw(run['out'], START - timedelta(days=1), 'Rn')


# ---- the values ------------------------------------------------------------

def test_a_day_is_that_period_of_the_cell_budget(run):
    """Each term on day k is MF6's own record of period k, per cell area in
    mm/d -- the run's GWmap terms on one day, not their time mean."""
    ws = run['ws']
    area = 100.0
    for k in range(NPER):
        rg = DM.day_values(ws, NAME, k, NPER, 'Rg')
        assert np.allclose(rg, _cbc(ws, 'UZF-GWRCH', k=k) / area * 1000.0)
        ds = DM.day_values(ws, NAME, k, NPER, 'dSg')
        assert np.allclose(ds, (_cbc(ws, 'STO-SS', k=k) + _cbc(ws, 'STO-SY', k=k))
                           / area * 1000.0)
    days = [DM.day_values(ws, NAME, k, NPER, 'Rg') for k in range(NPER)]
    assert not np.allclose(days[0], days[1]), 'each day its own'
    k = 1
    v = {key: DM.day_values(ws, NAME, k, NPER, key)
         for key in ('Rn', 'Re', 'Rg', 'Rsw', 'EXFg', 'ETg', 'Eg')}
    assert np.allclose(v['EXFg'],
                       -_cbc(ws, 'DRN', 'DRN_SEEP', k=k) / area * 1000.0)
    assert np.allclose(v['Eg'], -_cbc(ws, 'EVT', 'EVT_EG', k=k) / area * 1000.0)
    assert np.allclose(v['Re'], v['Rg'] + v['Rsw'] - v['EXFg'])
    assert np.allclose(v['Rn'], v['Re'] - v['ETg'])
    assert v['EXFg'].max() > 0 and v['Eg'].max() > 0
    # the head: the period's end
    import flopy
    h = flopy.utils.HeadFile(os.path.join(ws, '%s.hds' % NAME))
    want = h.get_data(kstpkper=h.get_kstpkper()[k]).reshape(NLAY, -1)
    h.close()
    assert np.allclose(DM.day_values(ws, NAME, k, NPER, 'head'), want)
    # the flow to the layer below: none below the bottom layer, so that
    # layer is not drawn
    assert DM.day_values(ws, NAME, k, NPER, 'FLF').shape == (NLAY - 1, 6)


def _h5(ws, tmp_path):
    """A coupled h5 of 4 soil columns: cells 0, 1, 3, 5 of the 2 x 3 grid,
    two horizons of 100 and 300 mm."""
    import h5py
    from marmites_indices import INDEX_MM_SOIL as IXS
    p = str(tmp_path / 'c.h5')
    ij = np.array([[0, 0], [0, 1], [1, 0], [1, 2]])
    ms = np.zeros((4, 2, 9))
    ms[:, 0, IXS['iSsoil']], ms[:, 0, IXS['iSsoil_pc_s']] = 20.0, 0.2
    ms[:, 1, IXS['iSsoil']], ms[:, 1, IXS['iSsoil_pc_s']] = 30.0, 0.1
    ss = np.array([[40.0, 80.0, 120.0, 160.0]] * NPER)
    ss[2] *= 0.5
    with h5py.File(p, 'w') as h:
        h['cell_ij'] = ij
        h['ssoil'] = ss
        h['perc'] = np.full((NPER, 4), 0.002)            # m/d
        h['etuzf'] = np.ones((NPER, 4))
        h['runoff'] = np.zeros((NPER, 4))
        h['wb_map_soil'] = ms
        h['obs_names'] = np.array([b'P0'])
        h['obs_ij'] = np.array([[0, 1]])
    return p


def test_the_soil_columns_moisture_is_its_storage_over_its_thickness(
        run, tmp_path):
    """The h5 keeps each column's storage; its thickness comes from the
    run's means of each horizon's storage and moisture (100 + 300 mm)."""
    h5 = _h5(run['ws'], tmp_path)
    th = DM.day_values(run['ws'], NAME, 0, NPER, 'theta', h5=h5)
    assert th.shape == (1, NROW * NCOL)
    cells = [0, 1, 3, 5]
    assert np.allclose(th[0, cells], [10.0, 20.0, 30.0, 40.0])
    assert np.isnan(th[0, [2, 4]]).all(), 'no soil column there'
    assert np.allclose(DM.day_values(run['ws'], NAME, 2, NPER, 'theta',
                                     h5=h5)[0, cells], [5.0, 10.0, 15.0, 20.0])
    assert np.allclose(DM.day_values(run['ws'], NAME, 0, NPER, 'Rp',
                                     h5=h5)[0, cells], 2.0), 'm/d -> mm/d'
    with pytest.raises(ValueError):
        DM.day_values(run['ws'], NAME, 0, NPER, 'Ssoil', h5=None)


# ---- the figures and the panel ---------------------------------------------

def test_the_day_is_drawn_into_the_runs_folder(run, tmp_path):
    out = run['out']
    day = START + timedelta(days=1)
    got = DM.draw(out, day, 'Rn')
    names = [os.path.basename(p) for p in got]
    assert names == ['daymap_GW_Rn_%s.png' % day, 'daymap_GW_Rn_%s_pct.png'
                     % day], names
    assert all(os.path.dirname(p) == os.path.join(out, '_output', DM.DAY_DIR)
               for p in got)
    info = dict(DM.outputs_of(out), h5=_h5(run['ws'], tmp_path))
    DM.draw(out, START, 'theta', info=info)
    # nothing left under plotLAYER's own names
    left = os.listdir(os.path.join(out, '_output', DM.DAY_DIR))
    assert all(f.startswith('daymap_') for f in left), left
    # the latest first, each map before its twin
    order = [f for _s, f, _p in DM.drawn(out)]
    assert order[:4] == ['daymap_MM_theta_%s.png' % START,
                         'daymap_MM_theta_%s_pct.png' % START,
                         'daymap_GW_Rn_%s.png' % day,
                         'daymap_GW_Rn_%s_pct.png' % day], order
    from PIL import Image
    with Image.open(got[0]) as im:
        assert im.info.get('MM-panels') == '2', 'one panel per layer'
    # the run's own tabs never list them
    assert not any(f.startswith('daymap_')
                   for f in os.listdir(os.path.join(out, '_output')))


def test_the_panel_titles_and_rows_of_the_day_maps():
    assert R.title('_output/day_maps', 'daymap_GW_Rn_2013-06-15.png') == \
        'Aquifer: net recharge (Rg + Rsw - Exf - ETg) on 2013-06-15'
    assert R.title('_output/day_maps', 'daymap_GW_Rn_2013-06-15_pct.png') == \
        ('Aquifer: net recharge (Rg + Rsw - Exf - ETg) on 2013-06-15 '
         '(percentile colour ramp)')
    assert R.title('_output/day_maps', 'daymap_MM_theta_2013-06-15.png') == \
        'Soil column: soil moisture on 2013-06-15'
    assert R.title('_output/day_maps', 'daymap_GW_head_2013-06-15.png') == \
        'Aquifer: head on 2013-06-15'
    # every flux of the list has its words on the panel
    for key, fam, *_ in DM.FLUXES:
        t = R.title('x', 'daymap_%s_%s_2013-06-15.png' % (fam, key))
        assert key == 'head' or ': %s on' % key not in t, t
    names = ['daymap_MM_theta_2013-06-15.png',
             'daymap_MM_theta_2013-06-15_pct.png',
             'daymap_MM_Rp_2013-06-14.png', 'daymap_MM_Rp_2013-06-14_pct.png']
    rows = [[it[1] for it in r] for r in R.rows([('s', n, n, n)
                                                for n in names])]
    assert rows == [names[:2], names[2:]], 'a map and its twin share a row'
    assert DM.label('Rn').startswith('Aquifer: Rn, net recharge')
    assert DM.label('theta').startswith('Soil column: theta, soil moisture')
