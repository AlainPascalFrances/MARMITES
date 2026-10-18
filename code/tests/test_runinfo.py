# -*- coding: utf-8 -*-
"""What a finished run leaves to build on (marmites_runinfo; WP7 front-end,
phase 3): the configuration it ran and the state it saved where it ended
-- recorded by a run the front-end launches, recovered from an older run's
log and results folder, and refused once a later run has saved over it."""

import json
import os
import sys
import time
from types import SimpleNamespace

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in (TRUNK, os.path.join(TRUNK, 'tools'), os.path.join(TRUNK, 'app')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg                                  # noqa: E402
import marmites_runinfo as ri                                   # noqa: E402


def _cfg(nlay=2):
    cfg = mcfg.RunConfig.from_dict({})
    cfg.layers.nlay = nlay
    return cfg


def _status(run_dir, outcome='completed', started='2026-10-11 21:49:04'):
    os.makedirs(run_dir, exist_ok=True)
    with open(os.path.join(run_dir, 'status.json'), 'w') as fh:
        json.dump({'run_id': os.path.basename(run_dir), 'state': 'finished',
                   'outcome': outcome, 'started': started,
                   'finished': '2026-10-11 22:45:38'}, fh)


def _state(state_dir, prefix, nlay=2, when=None):
    os.makedirs(state_dir, exist_ok=True)
    for k in range(nlay):
        p = os.path.join(state_dir, '%s_l%d.asc' % (prefix, k + 1))
        open(p, 'w').write('x\n')
        if when is not None:
            os.utime(p, (when, when))


def _log(run_dir, lines, when=None):
    p = os.path.join(run_dir, 'run.log')
    open(p, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
    if when is not None:
        os.utime(p, (when, when))


def test_record_merges_into_one_file(tmp_path):
    d = str(tmp_path / 'run')
    ri.record(d, config=ri.CONFIG, out_dir='x')
    ri.record(d, state='hi_a')
    assert ri.read(d) == {'config': ri.CONFIG, 'out_dir': 'x',
                          'state': 'hi_a'}
    assert ri.record('', state='x') is None          # no folder, nothing
    assert ri.read(str(tmp_path / 'none')) == {}


def test_a_recorded_run_is_a_base(tmp_path):
    rd, ws = str(tmp_path / 'runs' / 'r1'), str(tmp_path / 'ws')
    _status(rd)
    cfg = _cfg()
    cfg.write_toml(os.path.join(rd, ri.CONFIG))
    now = time.time()
    _state(ws, 'hi_end', when=now - 100)
    ri.record(rd, config=ri.CONFIG, out_dir='o', state_dir=ws,
              state='hi_end', converged=True)
    _log(rd, ['run'], when=now)
    b = ri.base_of(rd, str(tmp_path))
    assert b['why'] == '' and b['source'] == 'run_info'
    assert b['state'] == 'hi_end' and b['state_dir'] == ws
    assert b['config'] == os.path.join(rd, ri.CONFIG)


def test_an_older_run_is_recovered_from_its_log(tmp_path):
    """Run 20261011214903 left no record, but its log names the state it
    saved and its results folder holds its resolved configuration."""
    root = tmp_path / 'WS_ROOT'
    rd = str(root / 'runs' / '20261011214903_2lay')
    _status(rd)
    out = root / 'out_202610112149_2lay'
    (out / '_input').mkdir(parents=True)
    cfg = _cfg()
    cfg.write_toml(str(out / '_input' / ri.CONFIG))
    now = time.time()
    _state(mcfg.state_workspace(cfg, str(root)), 'hi_drn1', when=now - 50)
    _log(rd, ['config: x.toml  (hash 81b349e8a030b63f)',
              'final heads saved: hi_drn1_l1.asc, hi_drn1_l2.asc, '
              'hi_drn1_state.npz',
              'results written to %s' % out], when=now)
    b = ri.base_of(rd, str(root))
    assert b['why'] == '', b['why']
    assert b['source'] == 'log' and b['state'] == 'hi_drn1'
    assert b['config'] == str(out / '_input' / ri.CONFIG)
    # no line naming the folder: found by the minute the run started
    _log(rd, ['final heads saved: hi_drn1_l1.asc, hi_drn1_l2.asc'],
         when=now)
    assert ri.base_of(rd, str(root))['config'] == str(out / '_input' /
                                                       ri.CONFIG)
    # a spin-up saved before it converged: usable, and said
    _log(rd, ['NOT-converged spin-up heads saved: hi_drn1_l1.asc',
              'results written to %s' % out], when=now)
    b = ri.base_of(rd, str(root))
    assert b['why'] == '' and b['converged'] is False


def test_a_state_saved_again_by_a_later_run_is_refused(tmp_path):
    root = tmp_path / 'WS_ROOT'
    rd = str(root / 'runs' / 'r1')
    _status(rd)
    cfg = _cfg()
    cfg.write_toml(os.path.join(rd, ri.CONFIG))
    now = time.time()
    _state(str(tmp_path / 'ws'), 'hi_evt', when=now)
    ri.record(rd, config=ri.CONFIG, state='hi_evt',
              state_dir=str(tmp_path / 'ws'))
    _log(rd, ['run'], when=now - 3600)          # ended an hour before
    b = ri.base_of(rd, str(root))
    assert 'saved again by a later run' in b['why']


@pytest.mark.parametrize('case, frag', [
    ('failed', 'did not complete'),
    ('no_state', 'saved no end state'),
    ('no_config', 'no copy of its configuration'),
    ('gone', 'not in'),
])
def test_what_cannot_be_a_base_says_why(tmp_path, case, frag):
    root = tmp_path / 'WS_ROOT'
    rd = str(root / 'runs' / 'r1')
    _status(rd, outcome='failed' if case == 'failed' else 'completed')
    cfg = _cfg()
    if case != 'no_config':
        cfg.write_toml(os.path.join(rd, ri.CONFIG))
    info = {'config': ri.CONFIG, 'state_dir': str(tmp_path / 'ws')}
    if case != 'no_state':
        info['state'] = 'hi_x'
    if case not in ('gone', 'no_state'):
        _state(str(tmp_path / 'ws'), 'hi_x', when=time.time() - 60)
    ri.record(rd, **info)
    _log(rd, ['run'])
    assert frag in ri.base_of(rd, str(root))['why']


def test_the_runs_folder_is_the_configurations_or_ws_root(tmp_path):
    cfg = _cfg()
    assert ri.runs_dir(cfg, str(tmp_path)) == os.path.join(str(tmp_path),
                                                           'runs')
    cfg.ui.runs_dir = str(tmp_path / 'elsewhere')
    assert ri.runs_dir(cfg, 'x') == str(tmp_path / 'elsewhere')


def test_the_registry_and_the_record_agree_on_the_name():
    from lib import runs as runlib
    assert runlib._RUN_DIR_ENV == ri.RUN_DIR_ENV


def test_a_launched_run_is_told_its_folder(tmp_path):
    """The registry hands the run its folder, so the run can record its
    configuration and end state beside its log."""
    from lib import runs as runlib
    cfgp = tmp_path / 'c.toml'
    cfgp.write_text('[meta]\nconfig_version = 1\n', encoding='utf-8')
    driver = tmp_path / 'fake_driver.py'
    driver.write_text('import os\nprint("dir=" + os.environ.get("%s", ""), '
                      'flush=True)\n' % ri.RUN_DIR_ENV, encoding='utf-8')
    runs = tmp_path / 'runs'
    rid, st = runlib.launch(str(cfgp), str(runs), driver=str(driver),
                            python_exe=sys.executable, cwd=str(tmp_path))
    for _ in range(100):
        if not runlib._pid_alive(st['pid']):
            break
        time.sleep(0.1)
    assert ('dir=%s' % (runs / rid)) in runlib.log_tail(str(runs), rid)


def test_the_driver_records_the_configuration_it_ran(tmp_path):
    """run_lamata_mf6.record_run: the resolved configuration and where the
    run writes, into its registry folder -- nothing without one."""
    from lib import runs as runlib
    runner = runlib.import_runner()
    cfg = _cfg()
    cfg.meta.name = 'tagged'
    a = SimpleNamespace(out_dir=str(tmp_path / 'out_x'),
                        ws=str(tmp_path / 'ws'))
    rd = str(tmp_path / 'runs' / 'r1')
    runner.record_run(cfg, a, run_dir=rd)
    back = mcfg.load_run_config(os.path.join(rd, ri.CONFIG))
    assert back.meta.name == 'tagged'
    assert ri.read(rd) == {'config': ri.CONFIG, 'out_dir': a.out_dir,
                           'state_dir': a.ws}
    os.environ.pop(ri.RUN_DIR_ENV, None)
    assert runner.record_run(cfg, a) is None


def test_the_setup_takes_the_run_as_its_base(tmp_path, monkeypatch):
    """pest.base_run: the run's configuration is the model, with THIS
    configuration's [pest]; its end state the start. Blank: this
    configuration and the state chosen by hand."""
    import pest_lamata_setup as pls
    root = tmp_path / 'WS_ROOT'
    rd = str(root / 'runs' / 'r1')
    _status(rd)
    run_cfg = _cfg()
    run_cfg.run.nsp = 365
    run_cfg.write_toml(os.path.join(rd, ri.CONFIG))
    _state(str(tmp_path / 'ws'), 'hi_end', when=time.time() - 60)
    ri.record(rd, config=ri.CONFIG, state='hi_end',
              state_dir=str(tmp_path / 'ws'))
    _log(rd, ['run'])
    live = _cfg()
    live.pest.nreal = 42
    with pytest.raises(SystemExit, match='no base'):
        pls.resolve_base(live, str(root))
    live.pest.base_state = 'hi_hand'
    model, state, sdir, what = pls.resolve_base(live, str(root))
    assert model is live and state == 'hi_hand'
    assert sdir == mcfg.state_workspace(live, str(root))
    assert 'chosen by hand' in what
    live.pest.base_run = 'r1'
    model, state, sdir, what = pls.resolve_base(live, str(root))
    assert model.run.nsp == 365 and model.pest.nreal == 42
    assert state == 'hi_end' and sdir == str(tmp_path / 'ws')
    assert 'run r1' in what
    live.pest.base_run = 'nope'
    with pytest.raises(SystemExit, match='cannot be the base'):
        pls.resolve_base(live, str(root))
