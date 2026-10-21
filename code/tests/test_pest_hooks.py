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


def test_the_log_gives_the_change_from_the_reference_in_percent(tmp_path):
    """user, 2026-10-21: 'drn_l1 -87.1 %', not 'drn_l1 x0.1294'."""
    c = _cmf()
    _write(tmp_path, 'k_l1.txt', [0.006111, 1.0, 1.0, 21.68])
    _write(tmp_path, 'sy_l2.txt', [2.5])
    _write(tmp_path, 'drn_l1.txt', [0.1294])
    said = []
    mp.apply_multipliers(c, str(tmp_path), out=said.append)
    line = [s for s in said if s.startswith('PEST multipliers')][0]
    assert 'change from the reference model' in line
    assert 'drn_l1 -87.1 %' in line
    assert 'k_l1 -99.4 % .. +2068 %' in line
    assert 'sy_l2 +150 %' in line
    assert ' x' not in line.split(': ', 1)[1]
    assert [mp.pct_change(m) for m in (1.0, 1.0004, 0.9, 1.25, 224.0)] == \
        ['0 %', '0 %', '-10.0 %', '+25.0 %', '+22300 %']


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


def test_a_sy_multiplier_keeps_the_uzf_reference_and_sr_files_are_gone(
        tmp_path):
    """The specific retention is layers.sr since 2026-10-20 (the forward
    run writes it into the configuration): an sr_l<L>.txt is not read. A
    Sy multiplier keeps the model's (thtr, thts) = (Sr, Sy + Sr) from
    before it, for a carried UZF state saved without its own."""
    c = _cmf()
    c.sr_actual = [0.275, 0.04]
    _write(tmp_path, 'sy_l1.txt', [2.5])
    _write(tmp_path, 'sr_l1.txt', [0.9])                 # not read
    done = mp.apply_multipliers(c, str(tmp_path), out=None)
    thtr0, thts0 = c.uzf_carry_ref
    assert np.allclose(thtr0[0], 0.275) and np.allclose(thtr0[1], 0.04)
    assert np.allclose(thts0[0], 0.01 + 0.275)           # Sy before x2.5
    assert np.allclose(thts0[1], 0.01 + 0.04)
    assert 'sr_l1' not in done and 'sy_l1' in done
    assert not hasattr(c, 'thts') or not np.allclose(
        np.asarray(c.thts, dtype=float), 0.9)


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


# (Sy + Sr above 1 is refused by the build since 2026-10-20, wherever Sr
# comes from: marmites_props.uzf_water_contents, test_uzf_sy)
