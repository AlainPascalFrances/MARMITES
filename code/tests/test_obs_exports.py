# -*- coding: utf-8 -*-
"""The observation exports, WP7's contract (cookbook WP6.4).

obs_heads.csv, obs_sm.csv, obs_sfr.csv and obs_et.csv: long format, ONE ROW
PER STRESS PERIOD for each point (and layer) whatever adaptive time stepping
did, the observed value of the period beside the simulated one, a stable
obsnme -- written on every run, figures or not.
"""
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for p in (CODE, os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'MARMITESutilities'),
          os.path.join(CODE, 'MARMITESutilities', 'MARMITESplot')):
    if p not in sys.path:
        sys.path.insert(0, p)

pd = pytest.importorskip('pandas')
PP = pytest.importorskip('marmites_postprocess')
from marmites_indices import INDEX_MM, INDEX_MM_SOIL         # noqa: E402

NPER = 6
DAYS = pd.date_range('2008-06-01', periods=NPER, freq='D')


def _dataset(ds):
    w = lambda fn, s: open(os.path.join(ds, fn), 'w').write(s)   # noqa: E731
    w('inputDATE.txt', ''.join('%s\n' % d.strftime('%Y-%m-%d') for d in DAYS))
    w('inputObs.txt', '# Name X Y lay\nP1 0 0 2 2.0\nP2 1 1 1 2.0\n')
    w('inputObsHEADS_P1.txt', '2008-06-02\t770.5\n2008-06-05\t770.1\n'
                              '2009-01-01\t771.0\n')
    w('inputObsSM_P1.txt', ''.join('%s\t%.3f\t%.3f\n'
                                   % (d.strftime('%Y-%m-%d'), 0.10 + k / 100,
                                      0.20 + k / 100)
                                   for k, d in enumerate(DAYS[:3])))
    w('inputObsRo_catchment.txt', '2008-06-03\t0.5\n2008-06-04\t9999.999\n')


def _run():
    nidx, nidxs = len(INDEX_MM), len(INDEX_MM_SOIL)
    mm_obs = np.zeros((NPER, 2, nidx))
    mm_obs[:, :, INDEX_MM['ihcorr']] = 768.0
    mm_obs[:, 0, INDEX_MM['iETtot']] = np.arange(NPER)
    mms = np.zeros((NPER, 2, 2, nidxs))
    mms[:, 0, 1, INDEX_MM_SOIL['iSsoil_pc_s']] = 0.33
    wb = np.zeros((NPER, nidx))
    wb[:, INDEX_MM['iETtot']] = 2.5
    res = {'mm_obs': mm_obs, 'mms_obs': mms,
           'obs_ij': np.array([[0, 0], [1, 0]]),
           'obs_names': np.array([b'P1', b'P2']), 'wb_ts': wb}
    ctx = SimpleNamespace(index=INDEX_MM, index_S=INDEX_MM_SOIL, _nsl=[2],
                          gridSOIL=np.ones((2, 1), int),
                          geom=SimpleNamespace(area=np.array([5e5, 5e5])))
    cMF = SimpleNamespace(nlay=2, perlen=[1.0] * NPER)
    return cMF, ctx, res


def test_four_files_one_row_per_period(tmp_path, monkeypatch):
    ds = str(tmp_path / 'ds')
    os.makedirs(ds)
    _dataset(ds)
    PP.OBS.update(table='inputObs.txt', heads='inputObsHEADS',
                  sm='inputObsSM', ro='inputObsRo', aet='')
    heads = np.full((2, 2, NPER), 765.0)
    heads[0, 1] = 766.0 + np.arange(NPER)            # P1 sits in layer 2
    monkeypatch.setattr(PP, '_obs_head_series', lambda *a, **k: heads)
    monkeypatch.setattr(PP, 'sfr_observations', lambda *a, **k: pd.DataFrame(
        {'outflow': -np.full(NPER, 1000.0), 'perlen': np.ones(NPER)}))
    cMF, ctx, res = _run()
    outs = [str(tmp_path / 'a'), str(tmp_path / 'b')]
    written = PP.export_observations('ws', 'x', ds, cMF, ctx, res, outs,
                                     verbose=False)
    assert len(written) == 8                          # 4 files x 2 places
    h = pd.read_csv(os.path.join(outs[0], 'obs_heads.csv'))
    sm = pd.read_csv(os.path.join(outs[0], 'obs_sm.csv'))
    q = pd.read_csv(os.path.join(outs[0], 'obs_sfr.csv'))
    et = pd.read_csv(os.path.join(outs[1], 'obs_et.csv'))
    # deterministic: points x layers x periods, whatever was observed
    assert len(h) == 2 * NPER and len(sm) == 2 * 2 * NPER
    assert len(q) == NPER and len(et) == 3 * NPER     # + the catchment
    for df in (h, sm, q, et):
        assert df['obsnme'].is_unique
    # heads: the point's own layer, the observation in its period only
    p1 = h[h.point == 'P1']
    assert (p1['layer'] == 2).all()
    assert p1['sim_head_m'].tolist() == pytest.approx(766.0 + np.arange(NPER))
    assert p1['obs_head_m'].iloc[1] == pytest.approx(770.5)
    assert p1['obs_head_m'].iloc[4] == pytest.approx(770.1)
    assert p1['obs_head_m'].isna().sum() == NPER - 2
    assert (p1['sim_hcorr_m'] == 768.0).all()
    # soil moisture per horizon
    s2 = sm[(sm.point == 'P1') & (sm.soil_layer == 2)]
    assert s2['sim_theta'].tolist() == pytest.approx([0.33] * NPER)
    assert s2['obs_theta'].iloc[:3].tolist() == pytest.approx([0.20, 0.21,
                                                               0.22])
    assert s2['obs_theta'].iloc[3:].isna().all()
    # streamflow: m3/d and mm/d over the catchment (1 km2); a gap is a gap
    assert q['sim_q_mmd'].tolist() == pytest.approx([1.0] * NPER)
    assert q['obs_q_mmd'].iloc[2] == pytest.approx(0.5)
    assert q['obs_q_mmd'].isna().sum() == NPER - 1
    # ET: every source, per point and the catchment
    assert et[et.site == 'P1']['sim_et_total_mmd'].tolist() == \
        pytest.approx(list(range(NPER)))
    assert (et[et.site == 'catchment']['sim_et_total_mmd'] == 2.5).all()
    assert {'sim_e_g_mmd', 'sim_t_g_mmd', 'sim_et_uzf_mmd',
            'obs_et_mmd'} <= set(et.columns)


def test_a_run_without_point_series_writes_nothing(tmp_path):
    assert PP.export_observations('ws', 'x', str(tmp_path), None, None, {},
                                  [str(tmp_path)], verbose=False) == []


def test_the_driver_exports_on_every_run():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    body = src[src.index('def _run_postproc'):]
    assert body.index('export_observations(') < \
        body.index('if a.postproc or a.preproc:')
    assert "os.path.join(a.ws, 'obs_exports')" in body
    main = src[src.index('def main('):]
    blk = main[main.index('ON EVERY RUN, figures or not'):]
    assert blk.index('resolve_obs_cells(cMF, ctx, DS') < blk.index('ncyc =')
    # the panel's actual-ET prefix is read
    pp = open(os.path.join(CODE, 'ppMF6', 'marmites_postprocess.py'),
              encoding='utf-8').read()
    assert "aet=(o.aet_prefix or '').strip()" in pp


def test_soil_moisture_at_depth_one_page_per_point(tmp_path):
    """WP6.3: one page per observation point, one curve per horizon, with
    or without measurements."""
    import matplotlib
    matplotlib.use('agg')
    ds = str(tmp_path / 'ds')
    os.makedirs(ds)
    _dataset(ds)
    PP.OBS.update(sm='inputObsSM')
    cMF, ctx, res = _run()
    ctx._Sm, ctx._Sfc, ctx._Sr = [[0.4, 0.4]], [[0.3, 0.3]], [[0.05, 0.05]]
    ctx._slprop = [[0.3, 0.7]]
    ctx.gridSOILthick = np.full((2, 1), 0.5)
    out = str(tmp_path / 'fig')
    os.makedirs(out)
    got = PP._fig_sm_depth(out, cMF, ctx, res, ds, verbose=False)
    assert [os.path.basename(f) for f in got] == ['sm_depth_P1.png',
                                                  'sm_depth_P2.png']
    assert all(os.path.getsize(f) > 0 for f in got)
