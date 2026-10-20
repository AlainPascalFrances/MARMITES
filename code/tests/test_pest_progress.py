# -*- coding: utf-8 -*-
"""The PEST++ progress panel's reading of a calibration (WP7 front-end,
phase 1): pestpp-ies's record files turned into batches, runs, phi and
failures. On a small synthetic master folder in the formats pestpp-ies
5.2 writes -- no PEST++ needed."""

import datetime as dt
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.abspath(os.path.join(HERE, '..', 'app'))
if APP not in sys.path:
    sys.path.insert(0, APP)

from lib import pest_progress as pp                       # noqa: E402

W = r'HOST$X:\P\workers\worker_%d'

RMR = '\n'.join([
    'running model 4 times',
    '10/16/26 17:51:46->    3 agents ready',
    '10/16/26 17:51:46->Sending run_id:0 to:' + W % 0 + ' group_id:2  '
    'da_cycle:-9999   iteration:0  realization:0 concurrent:1)',
    '10/16/26 17:51:47->Sending run_id:1 to:' + W % 1 + ' group_id:2  '
    'da_cycle:-9999   iteration:0  realization:5 concurrent:1)',
    '10/16/26 17:51:48->Sending run_id:2 to:' + W % 2 + ' group_id:2  '
    'da_cycle:-9999   iteration:0  realization:7 concurrent:1)',
    '10/16/26 19:40:00-> run_id:0 received from:' + W % 0 + '  '
    'run_time:108.2  avg_run_time:108.2 group_id:2  da_cycle:-9999',
    '10/16/26 19:40:01->Sending run_id:3 to:' + W % 0 + ' group_id:2  '
    'da_cycle:-9999   iteration:0  realization:BASE concurrent:1)',
    '10/16/26 20:10:00-> run_id:2 failed on agent:' + W % 2 + ' group_id:2 '
    'group_id:2 run_id:2  da_cycle:-9999',
    '10/16/26 21:00:00->sending kill request. reason:overdue. '
    'duration:432.0, avg:108.0 run_id:1,  num_previous_fails:0 agent:'
    + W % 1,
    '10/16/26 21:00:01-> run_id:1 killed on agent:' + W % 1,
    '10/16/26 21:40:00-> run_id:3 received from:' + W % 0 + '  '
    'run_time:120.0  avg_run_time:114.1 group_id:2  da_cycle:-9999',
    '   2 runs complete :  2 runs failed',
    'running model 2 times',
    '10/17/26 02:19:09->Sending run_id:0 to:' + W % 0 + ' group_id:3  '
    'da_cycle:-9999   lambda:0.01 scale_fac:1 iteration:1  realization:0 '
    'concurrent:1)',
    '10/17/26 02:19:09->Sending run_id:1 to:' + W % 1 + ' group_id:3  '
    'da_cycle:-9999   lambda:0.1 scale_fac:1 iteration:1  realization:0 '
    'concurrent:1)',
    '10/17/26 03:49:09-> run_id:0 received from:' + W % 0 + '  '
    'run_time:90.0  avg_run_time:90.0 group_id:3',
    # the overdue run sent to a second agent, which brings it back first
    '10/17/26 03:50:00->Sending run_id:1 to:' + W % 0 + ' group_id:3  '
    'da_cycle:-9999   lambda:0.1 scale_fac:1 iteration:1  realization:0 '
    'concurrent:2)',
    '10/17/26 04:10:00-> run_id:1 received from:' + W % 0 + '  '
    'run_time:20.0  avg_run_time:55.0 group_id:3',
    'running model 3 times',
    '10/17/26 05:00:00->Sending run_id:0 to:' + W % 0 + ' group_id:4  '
    'da_cycle:-9999   iteration:1  realization:5 concurrent:1)',
    '10/17/26 05:00:00->Sending run_id:1 to:' + W % 1 + ' group_id:4  '
    'da_cycle:-9999   iteration:1  realization:BASE concurrent:1)',
]) + '\n'

REC = '\n'.join([
    '    noptmax = 2',
    'da_noptmax_schedule:',
    'ies_num_reals: 4',
    'overdue_giveup_fac: 4',
    'ies_lambda_mults: 0.1,1,10,',
    '  ---  initialization complete  ---  ',
    '  ---  starting solve for iteration: 1  ---  ',
    '  ---  running upgrade ensembles  ---  ',
    '  ---  evaluating upgrade ensembles  ---  ',
    '  ---  phi summary for lambda, scale fac:0.01 , 1 ,   ---  ',
    '       phi type           mean            std            min            '
    'max',
    '       measured        7216.27        689.985        5923.34        '
    '8059.82',
    '         actual        6538.98        483.049        5923.34        '
    '7403.08',
    '  ---  phi summary for lambda, scale fac:0.1 , 1 ,   ---  ',
    '       phi type           mean            std            min            '
    'max',
    '         actual        6554.89        490.769        5912.97        '
    '7417.42',
    '  ---  phi summary for best lambda, scale fac: 0.01 , 1 ,   ---  ',
    '       phi type           mean            std            min            '
    'max',
    '         actual        6538.98        483.049        5923.34        '
    '7403.08',
    '  ---  running remaining realizations for best lambda, scale:0.01 , '
    '1 ,   ---',
]) + '\n'

MODEL_LOG = '\n'.join([
    '   stress period 1949/1949 (100%): ...',
    'Traceback (most recent call last):',
    '  File "run_lamata_mf6.py", line 1683, in main',
    'marmites_coupler.CouplingError: MF6 solution is not valid:',
    '  - ETsoil + ETuzf + ETg exceeded PE + PT in 9 cell-period(s), by up '
    'to 184 mm/d. Total ET can never be above PET.',
]) + '\n'


def fake_pest_ws(root):
    """A PEST workspace with a master folder caught in iteration 1."""
    ws = os.path.join(str(root), 'PEST')
    m = os.path.join(ws, 'master')
    os.makedirs(os.path.join(m, 'report_1'))
    for k in range(3):
        os.makedirs(os.path.join(ws, 'workers', 'worker_%d' % k))
    w = lambda fn, s: open(fn, 'w').write(s)                 # noqa: E731
    w(os.path.join(m, 'lamata.pst'), 'pcf\n')
    w(os.path.join(m, 'lamata.rec'), REC)
    w(os.path.join(m, 'lamata.rmr'), RMR)
    w(os.path.join(m, 'lamata.phi.actual.csv'),
      'iteration,total_runs,mean,standard_deviation,min,max,0,3,BASE\n'
      '0,4,6900.0,700.0,6100.0,7600.0,6100.0,7000.0,7600.0\n'
      '1,6,6500.0,450.0,5900.0,7100.0,5900.0,6500.0,7100.0\n')
    w(os.path.join(m, 'lamata.phi.group.csv'),
      'iteration,total_runs,obs_realization,par_realization,hds_p0,sm_p0\n'
      '0,4,0,0,900.0,1700.0\n0,4,3,3,1000.0,1800.0\n'
      '1,6,0,0,800.0,1600.0\n')
    w(os.path.join(m, 'report_1', 'summary.txt'), 'iteration 1: ...\n')
    w(os.path.join(ws, 'workers', 'worker_2', 'model.log'), MODEL_LOG)
    return ws


def test_the_batches_are_counted_run_by_run():
    b = pp.parse_rmr(RMR)
    assert [x['stage'] for x in b] == ['initial ensemble', 'lambda test',
                                       'remaining realisations']
    assert [x['iteration'] for x in b] == [0, 1, 1]
    assert [x['planned'] for x in b] == [4, 2, 3]
    st = {k: r['state'] for k, r in b[0]['runs'].items()}
    assert st == {0: 'completed', 1: 'killed', 2: 'failed', 3: 'completed'}
    assert b[0]['runs'][1]['reason'] == 'overdue'
    assert b[0]['runs'][3]['real'] == 'BASE'
    assert b[0]['runs'][3]['time'] == pytest.approx(120.0)
    # a run sent twice is completed when one copy comes back
    assert b[1]['lambdas'] == [0.01, 0.1]
    assert b[1]['runs'][1]['state'] == 'completed'
    assert b[1]['runs'][1]['agent'] == 'worker_0'
    assert all(r['state'] == 'running' for r in b[2]['runs'].values())


def test_the_record_says_where_the_run_is():
    r = pp.parse_rec(REC)
    assert (r['noptmax'], r['nreal'], r['overdue_giveup_fac']) == (2, 4, 4.0)
    assert r['lambda_mults'] == [0.1, 1.0, 10.0]
    assert r['iteration'] == 1 and r['stage'] == 'remaining realisations'
    assert not r['finished']
    t = r['lambda_tests']
    assert [(x['lambda'], x['best']) for x in t] == [(0.01, False),
                                                     (0.1, False),
                                                     (0.01, True)]
    assert t[1]['mean'] == pytest.approx(6554.89)
    r = pp.parse_rec(REC + 'pestpp-ies analysis complete...\n'
                     'started at 10/12/26 10:37:58\n'
                     'finished at 10/13/26 06:28:38\n')
    assert r['finished'] and r['stage'] == 'finished'
    assert r['finished_at'] == '10/13/26 06:28:38'


def test_the_progress_of_a_running_calibration(tmp_path):
    ws = fake_pest_ws(tmp_path)
    now = dt.datetime(2026, 10, 17, 6, 0, 0)
    p = pp.progress(ws, now=now, running=True)
    assert p['case'] == 'lamata' and p['state'] == 'running'
    assert (p['iteration'], p['noptmax'], p['nreal']) == (1, 2, 4)
    rows = p['batches']
    assert [r['completed'] for r in rows] == [2, 2, 0]
    assert [r['failed'] for r in rows] == [1, 0, 0]
    assert [r['killed (overdue)'] for r in rows] == [1, 0, 0]
    assert [r['running'] for r in rows] == [0, 0, 2]
    assert rows[0]['mean run (min)'] == pytest.approx(114.1)
    assert rows[0]['shortest run (min)'] == pytest.approx(108.2)
    assert rows[0]['longest run (min)'] == pytest.approx(120.0)
    assert [(a['agent'], a['realisation']) for a in p['agents']] == [
        ('worker_0', '5'), ('worker_1', 'BASE')]
    # each agent's runs so far, in place of PEST++'s run id (user,
    # 2026-10-20): worker_0 brought back four, worker_1's one was killed
    assert [(a['runs done'], a['of which failed']) for a in p['agents']] \
        == [(4, 0), (1, 1)]
    assert 'run' not in p['agents'][0]
    assert p['agents'][0]['elapsed (min)'] == pytest.approx(60.0)
    f = {x['run']: x for x in p['failures']}
    assert f[1]['outcome'] == 'killed' and f[1]['cause'] == 'overdue'
    # worker_2 has had no run since: its log says why
    assert f[2]['outcome'] == 'failed'
    assert 'CouplingError' in f[2]['cause']
    assert 'exceeded PE + PT' in f[2]['cause']
    assert list(p['phi']['realisations']) == [3, 3]
    assert list(p['phi_long']['phi'])[:2] == [6100.0, 5900.0]  # real 0
    assert p['phi_groups'].set_index('iteration').loc[0, 'hds_p0'] == \
        pytest.approx(950.0)
    assert [k for k, _ in p['reports']] == [1]
    assert p['eta'] is not None and p['eta'] > now


def test_a_stopped_calibration_has_unfinished_runs_and_no_agents(tmp_path):
    p = pp.progress(fake_pest_ws(tmp_path), running=False)
    assert p['state'] == 'stopped' and p['eta'] is None
    assert not p['agents']
    assert [r['unfinished'] for r in p['batches']] == [0, 0, 2]


def test_no_master_folder_is_no_calibration(tmp_path):
    assert pp.progress(str(tmp_path)) is None


def test_each_failure_reads_its_own_kept_log(tmp_path):
    """2026-10-20: the forward run keeps a failed run's log in
    failed_runs/; each failure reads the one kept at its time, not only
    the agent's latest."""
    ws = fake_pest_ws(tmp_path)
    kept = os.path.join(ws, 'failed_runs')
    os.makedirs(kept)
    open(os.path.join(kept, '20261016_200950_worker_2_model.log'), 'w').write(
        'marmites_coupler.NonConvergence: NON-CONVERGENCE: stress period 247'
        ' needed more than 1500 linear-solver iterations\n')
    open(os.path.join(kept, '20261017_120000_worker_2_model.log'), 'w').write(
        'RuntimeError: a later failure\n')
    p = pp.progress(ws, now=dt.datetime(2026, 10, 17, 6, 0, 0), running=True)
    f = {x['run']: x for x in p['failures']}
    assert f[2]['cause'].startswith('marmites_coupler.NonConvergence')
