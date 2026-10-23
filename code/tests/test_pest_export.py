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
