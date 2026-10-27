# -*- coding: utf-8 -*-
"""The calibrated model as a configuration of its own (WP7 phase 5,
tools/pest_lamata_export.py): a pilot-point property is a multiplier of its
reference raster, so a raster changed since the set-up (layer 1's Sy, 0.01
-> 0.025 on 2026-10-22) must rescale its multipliers -- or the exported run
would carry 2.5 times the calibrated Sy."""

import os
import sys
import types

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'tools')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import pest_lamata_export as plx                                # noqa: E402


def _asc(path, value, nodata_cell=True):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write('ncols 3\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 10\n'
                 'NODATA_value -9999\n')
        fh.write('%s %s %s\n%s %s %s\n' % ((-9999 if nodata_cell else value,)
                                           + (value,) * 5))


def _cfg(raster='MF_ws/Sy_l%d.asc'):
    src = types.SimpleNamespace(raster=raster)
    blank = types.SimpleNamespace(raster='')
    return types.SimpleNamespace(layers=types.SimpleNamespace(
        nlay=2, sy=src, k=blank, ss=blank))


def test_a_changed_uniform_reference_rescales_its_multipliers(tmp_path):
    old, new, mult = tmp_path / 'tpl', tmp_path / 'user', tmp_path / 'mm'
    mult.mkdir()
    _asc(str(old / 'MF_ws' / 'Sy_l1.asc'), 0.01)
    _asc(str(new / 'MF_ws' / 'Sy_l1.asc'), 0.025)
    _asc(str(old / 'MF_ws' / 'Sy_l2.asc'), 0.01)      # layer 2 unchanged
    _asc(str(new / 'MF_ws' / 'Sy_l2.asc'), 0.01)
    np.savetxt(str(mult / 'sy_l1.txt'), [2.5, 3.0])
    np.savetxt(str(mult / 'sy_l2.txt'), [1.5, 0.5])
    notes = plx.rescale(_cfg(), str(mult), str(old), str(new))
    assert np.allclose(np.loadtxt(str(mult / 'sy_l1.txt')), [1.0, 1.2])
    assert np.allclose(np.loadtxt(str(mult / 'sy_l2.txt')), [1.5, 0.5])
    assert len(notes) == 1 and 'x0.4' in notes[0]
    # the calibrated Sy is what the run gets: 0.025 x 1.0 = 0.01 x 2.5
    assert 0.025 * 1.0 == pytest.approx(0.01 * 2.5)


def test_a_changed_raster_that_is_not_uniform_is_refused(tmp_path):
    old, new, mult = tmp_path / 'tpl', tmp_path / 'user', tmp_path / 'mm'
    mult.mkdir()
    _asc(str(old / 'MF_ws' / 'Sy_l1.asc'), 0.01)
    p = new / 'MF_ws' / 'Sy_l1.asc'
    _asc(str(p), 0.02)
    txt = open(str(p)).read().replace('0.02 0.02 0.02\n', '0.03 0.02 0.02\n',
                                      1)
    open(str(p), 'w').write(txt)
    np.savetxt(str(mult / 'sy_l1.txt'), [2.0])
    with pytest.raises(SystemExit, match='not uniform'):
        plx.rescale(_cfg(), str(mult), str(old), str(new))


def test_the_export_starts_from_the_base_run_as_a_calibration_run():
    src = open(plx.__file__, encoding='utf-8').read()
    for piece in ("cfg.spinup.cycles = 1", "cfg.spinup.strt_heads = state",
                  "cfg.spinup.save_strt = state + '_PEST'",
                  "cfg.uzf.thti_from = 'thtr'",
                  "plf.apply_scalars(cfg, plf.read_scalars(work)",
                  "cfg.pest.multipliers = mult"):
        assert piece in src, piece


def _master(tmp_path, made='2026-10-25T15:47:53'):
    """A PEST workspace: a master that ran two iterations, a template."""
    import json
    ws = tmp_path / 'pest'
    m = ws / 'master'
    m.mkdir(parents=True)
    (m / 'lamata.phi.actual.csv').write_text(
        'iteration,total_runs,mean,standard_deviation,min,max,0,1,2,base\n'
        '0,4,10,1,9,11,11,9,,10\n1,4,8,1,7,9,9,7,6.5,8\n')
    for k in (0, 1):
        (m / ('lamata.%d.par.jcb' % k)).write_bytes(b'')
    (m / 'lamata.rec').write_text('pestpp-ies\nstarted at 10/25/26 17:05:19\n')
    t = ws / 'template'
    (t / 'ws' / 'MF6_ws').mkdir(parents=True)
    sdir = tmp_path / 'WS_ROOT' / 'MF6_ws'
    sdir.mkdir(parents=True)
    for f in ('s_l1.asc', 's_state.npz'):
        (t / 'ws' / 'MF6_ws' / f).write_text('h ' + f)
        (sdir / f).write_text('h ' + f)
    (sdir / 's_PEST_l1.asc').write_text('another state')
    (t / 'setup.json').write_text(json.dumps(
        {'config': str(tmp_path / 'lamata.toml'), 'state': 's',
         'state_dir': str(sdir), 'made': made, 'model_config': 'x.toml'}))
    return ws, sdir


def test_the_calibration_and_its_realisations(tmp_path):
    ws, _sdir = _master(tmp_path)
    m = str(ws / 'master')
    assert plx.calibration_start(m).isoformat() == '2026-10-25T17:05:19'
    assert plx.calibration_start(str(tmp_path)) is None
    assert plx.iterations(m) == [0, 1]
    # base first, then from the lowest phi; a failed one (NaN) left out
    assert list(plx.choices(m, 1).index) == ['base', '2', '1', '0']
    assert list(plx.choices(m, 0).index) == ['base', '1', '0']


def test_a_base_state_saved_over_since_the_set_up_is_refused(tmp_path):
    """Every run of the calibration started from the set-up's copy of the
    base state; the exported run starts from the live one."""
    ws, sdir = _master(tmp_path)
    tpl = str(ws / 'template')
    setup = plx.read_setup(str(ws))
    assert plx.state_changes(tpl, setup) == []
    (sdir / 's_l1.asc').write_text('a later run')
    assert plx.state_changes(tpl, setup) == [str(sdir / 's_l1.asc')]
    with pytest.raises(SystemExit, match='changed since the set-up'):
        plx.export(str(ws), 1, 'base', out=str(tmp_path / 'o.toml'))


def test_a_template_set_up_after_the_calibration_is_refused(tmp_path):
    ws, _sdir = _master(tmp_path, made='2026-10-26T09:00:00')
    with pytest.raises(SystemExit, match='set up again'):
        plx.export(str(ws), 1, 'base', out=str(tmp_path / 'o.toml'))


def test_an_existing_configuration_is_kept_without_force(tmp_path):
    ws, _sdir = _master(tmp_path)
    out = tmp_path / 'o.toml'
    out.write_text('[meta]\n')
    with pytest.raises(SystemExit, match='--force to replace'):
        plx.export(str(ws), 1, 'base', out=str(out))


def test_the_provenance_of_an_exported_configuration(tmp_path):
    import json
    work = tmp_path / 'results' / '202610251705_it3_base'
    (work / 'mm_mult').mkdir(parents=True)
    (work / 'export.json').write_text(json.dumps(
        {'made': '2026-10-27T18:00:00', 'iteration': 3,
         'realisation': 'base', 'phi': 6051.0,
         'calibration': '202610251705'}))
    cfg = tmp_path / 'lamata_PEST.toml'
    cfg.write_text('[pest]\nmultipliers = %r\n' % str(work / 'mm_mult'))
    import datetime as dt
    t = dt.datetime(2026, 10, 27, 18, 0, 30).timestamp()
    os.utime(str(cfg), (t, t))
    info = plx.provenance(str(cfg))
    assert info['realisation'] == 'base' and info['edited'] is False
    t = dt.datetime(2026, 10, 28, 9, 0, 0).timestamp()
    os.utime(str(cfg), (t, t))
    assert plx.provenance(str(cfg))['edited'] is True
    # not written by an export
    other = tmp_path / 'lamata.toml'
    other.write_text('[pest]\nmultipliers = ""\n')
    assert plx.provenance(str(other)) is None
    assert plx.provenance(str(tmp_path / 'none.toml')) is None
