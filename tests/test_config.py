# -*- coding: utf-8 -*-
"""Phase-2 tests for the TOML config module and the legacy-ini converter."""
import importlib.util
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.join(HERE, '..', 'trunk')
DATASET_INI = os.path.join(HERE, '..', 'DataSet_LaMata', '__inputMM_v3.ini')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cfgmod = _load('marmites_config', os.path.join(TRUNK, 'marmites_config.py'))


def test_convert_real_dataset_ini(tmp_path):
    if not os.path.exists(DATASET_INI):
        pytest.skip('DataSet_LaMata ini not present')
    toml_path = str(tmp_path / 'mm.toml')
    cfg = cfgmod.convert_ini_file(DATASET_INI, toml_path)
    # known values from the LaMata annotated ini (includes maxYearsTick fields)
    assert cfg.run_name.startswith('2s3L')
    assert cfg.verbose == 1
    assert cfg.iniMonthHydroYear == 10
    assert cfg.maxYearsTickTrimester == 10
    assert cfg.maxYearsTickSemester == 20
    assert cfg.plt_WB_unit == 'year'
    assert cfg.irr_yn is True
    assert cfg.gridIRR_fn == 'inputIRRzones.asc'
    assert cfg.xllcorner == 739300.0
    assert cfg.yllcorner == 4553050.0
    # a TOML file was written and reloads to an equivalent config
    assert os.path.exists(toml_path)
    cfg2 = cfgmod.load_mm_config(toml_path)
    assert cfg2.run_name == cfg.run_name
    assert cfg2.maxYearsTickSemester == cfg.maxYearsTickSemester
    assert cfg2.xllcorner == cfg.xllcorner


def test_roundtrip_toml(tmp_path):
    cfg = cfgmod.MMConfig(run_name='demo', iniMonthHydroYear=1, plt_WB_unit='day',
                          MMsoil_yn=-1, chunks=1)
    p = str(tmp_path / 'x.toml')
    cfgmod.dump_toml(cfg, p)
    back = cfgmod.load_mm_config(p)
    assert back.run_name == 'demo'
    assert back.iniMonthHydroYear == 1
    assert back.plt_WB_unit == 'day'
    assert back.MMsoil_yn == -1
    assert back.chunks == 1


def test_validation_rejects_bad_month():
    with pytest.raises(cfgmod.ConfigError):
        cfgmod.MMConfig(run_name='x', iniMonthHydroYear=13).validate()


def test_validation_rejects_bad_wb_unit():
    with pytest.raises(cfgmod.ConfigError):
        cfgmod.MMConfig(run_name='x', plt_WB_unit='week').validate()


def test_validation_rejects_empty_run_name():
    with pytest.raises(cfgmod.ConfigError):
        cfgmod.MMConfig(run_name='').validate()


def test_deprecated_fields_captured(tmp_path):
    if not os.path.exists(DATASET_INI):
        pytest.skip('DataSet_LaMata ini not present')
    cfg = cfgmod.convert_ini_file(DATASET_INI, str(tmp_path / 'm.toml'))
    # Picard-loop fields are parsed but quarantined as deprecated
    assert 'convcrit' in cfg.deprecated
    assert 'ccnum' in cfg.deprecated
