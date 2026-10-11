# -*- coding: utf-8 -*-
"""The model side of the PEST++ calibration (WP7, marmites_pest): the
multiplier files a forward run writes, applied on the grid the model runs
on, and configuration values set by name, list items included."""

import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
if TRUNK not in sys.path:
    sys.path.insert(0, TRUNK)

import marmites_pest as mp                                      # noqa: E402
import marmites_config as mcfg                                  # noqa: E402


def _cmf():
    """A mesh-shaped model (nrow 1, ncol = ncpl), K as a list of arrays and
    Ss as a list of floats -- both forms the projection leaves."""
    return SimpleNamespace(
        nlay=2, nrow=1, ncol=4,
        hk_actual=[np.array([[0.05, 0.1, 0.5, 0.05]]),
                   np.array([[0.05, 0.05, 0.1, 0.5]])],
        sy_actual=np.full((2, 1, 4), 0.01),
        ss_actual=[1e-4, 4e-7],
        layer_row_column_elevation_cond={0: [(0, 0, 1, 700.0, 0.004),
                                             (1, 0, 1, 650.0, 0.004)]})


def _write(folder, name, values):
    np.savetxt(os.path.join(folder, name), np.atleast_1d(values))


def test_per_cell_and_per_layer_multipliers(tmp_path):
    c = _cmf()
    _write(tmp_path, 'k_l1.txt', [2.0, 0.5, 10.0, 1.0])
    _write(tmp_path, 'sy_l2.txt', [3.0])
    _write(tmp_path, 'ss_l1.txt', [0.1])
    _write(tmp_path, 'drn_l2.txt', [50.0])
    done = mp.apply_multipliers(c, str(tmp_path), out=None)
    assert np.allclose(c.hk_actual[0], [[0.1, 0.05, 5.0, 0.05]])
    assert np.allclose(c.hk_actual[1], [[0.05, 0.05, 0.1, 0.5]])  # no file
    assert np.allclose(c.sy_actual[1], 0.03) and np.allclose(c.sy_actual[0],
                                                             0.01)
    assert np.allclose(c.ss_actual[0], 1e-5) and np.allclose(c.ss_actual[1],
                                                             4e-7)
    recs = c.layer_row_column_elevation_cond[0]
    assert recs[0][4] == pytest.approx(0.004)            # layer 1 untouched
    assert recs[1][4] == pytest.approx(0.2)              # layer 2 x 50
    assert set(done) == {'k_l1', 'sy_l2', 'ss_l1', 'drn_l2'}


def test_a_bad_multiplier_file_is_refused(tmp_path):
    c = _cmf()
    _write(tmp_path, 'k_l1.txt', [1.0, 2.0])             # 2 of 4 cells
    with pytest.raises(mp.PestError):
        mp.apply_multipliers(c, str(tmp_path), out=None)
    _write(tmp_path, 'k_l1.txt', [1.0, -2.0, 1.0, 1.0])  # negative
    with pytest.raises(mp.PestError):
        mp.apply_multipliers(c, str(tmp_path), out=None)
    with pytest.raises(mp.PestError):
        mp.apply_multipliers(c, str(tmp_path / 'absent'), out=None)


def test_no_folder_changes_nothing():
    c = _cmf()
    assert mp.apply_multipliers(c, '', out=None) == {}
    assert np.allclose(c.hk_actual[0], [[0.05, 0.1, 0.5, 0.05]])


def test_parameters_are_set_by_name_list_items_included():
    cfg = mcfg.RunConfig.from_dict({})
    cfg.soil.horizon = [mcfg.SoilHorizon(zone=1, ks=1.5),
                        mcfg.SoilHorizon(zone=2, ks=5.0)]
    assert mp.set_param(cfg, 'soil.horizon[1].ks', 7.5) == 5.0
    assert cfg.soil.horizon[1].ks == 7.5
    mp.set_param(cfg, 'drn.cond.value', 2e-4)
    assert cfg.drn.cond.value == 2e-4
    mp.set_param(cfg, 'crr.beta', 0.9)
    assert mp.get_param(cfg, 'crr.beta') == 0.9
    with pytest.raises(mp.PestError):
        mp.set_param(cfg, 'soil.horizon[9].ks', 1.0)
    with pytest.raises(mp.PestError):
        mp.set_param(cfg, 'sfr.enable', 1.0)              # a switch
    with pytest.raises(mp.PestError):
        mp.set_param(cfg, 'nope.x', 1.0)


def test_the_driver_applies_the_multipliers_before_building():
    src = open(os.path.join(TRUNK, 'tests', 'run_lamata_mf6.py'),
               encoding='utf-8').read()
    i = src.index('marmites_pest.apply_multipliers(cMF, a.pest_multipliers)')
    assert i < src.index('b = clsMF6(cMF, top=top')
    assert 'pest_multipliers=_pest_folder(cfg)' in src
