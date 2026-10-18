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


def test_specific_retention_sets_thts_and_keeps_the_reference(tmp_path):
    """sr_l<L> is a water content: thts = Sy + Sr in every cell (so the
    builder's thtr = thts - Sy = Sr), after the Sy multipliers; the
    reference (thtr, thts) is kept for a carried UZF state (2026-10-18)."""
    c = _cmf()
    c.thts_actual = [0.45, 0.45]
    _write(tmp_path, 'sy_l1.txt', [2.5])
    _write(tmp_path, 'sr_l1.txt', [0.275])
    _write(tmp_path, 'sr_l2.txt', [0.04])
    done = mp.apply_multipliers(c, str(tmp_path), out=None)
    thts = np.asarray(c.thts)
    assert thts.shape == (2, 1, 4)
    assert np.allclose(thts[0], 0.025 + 0.275)
    assert np.allclose(thts[1], 0.01 + 0.04)
    assert np.allclose(np.asarray(c.thts_actual), thts)
    thtr0, thts0 = c.uzf_carry_ref
    assert np.allclose(thts0, 0.45) and np.allclose(thtr0, 0.44)
    assert {'sr_l1', 'sr_l2', 'sy_l1'} <= set(done)


def test_layer_1_k_is_held_under_layer_2_k(tmp_path):
    c = _cmf()
    _write(tmp_path, 'k_l1_max_ratio.txt', [0.1])
    said = []
    mp.apply_multipliers(c, str(tmp_path), out=said.append)
    # layer 2 [0.05, 0.05, 0.1, 0.5] -> layer 1 at most [5e-3, 5e-3, 0.01,
    # 0.05]; it was [0.05, 0.1, 0.5, 0.05]
    assert np.allclose(c.hk_actual[0], [[0.005, 0.005, 0.01, 0.05]])
    assert np.allclose(c.hk_actual[1], [[0.05, 0.05, 0.1, 0.5]])
    assert any('3 cell(s)' in s for s in said)


def test_sy_and_sr_above_one_are_refused(tmp_path):
    c = _cmf()
    c.thts_actual = [0.45, 0.45]
    _write(tmp_path, 'sr_l1.txt', [0.995])
    _write(tmp_path, 'sy_l1.txt', [2.0])
    with pytest.raises(mp.PestError, match='above 1'):
        mp.apply_multipliers(c, str(tmp_path), out=None)
