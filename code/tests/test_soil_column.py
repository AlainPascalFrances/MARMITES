# -*- coding: utf-8 -*-
"""The soil column comes from the Soil panel, not from inputSOILparam.txt.

The audit behind the cookbook's Appendix B found that soil.params was never
read: MMsoil took its soil column from MF_ws/inputSOILparam.txt, at a path
hard-coded in the driver. The column is now two tables in the TOML --
[[soil.zone]] and [[soil.horizon]] -- and the acceptance test is the one the
drains had: what the run hands MMsoil is what the legacy reader handed it,
value for value.
"""

import importlib.util
import io
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(CODE, '..', 'example', 'LaMata'))
PARAM = os.path.join(DS, 'MF_ws', 'inputSOILparam.txt')
for _p in (CODE, os.path.join(CODE, 'ppMF6'),
           os.path.join(CODE, 'MARMITESutilities'),
           os.path.join(CODE, 'ppMF_FloPy')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


cfgmod = _load('marmites_config_sc', os.path.join(CODE, 'marmites_config.py'))
props = _load('marmites_props_sc', os.path.join(CODE, 'ppMF6',
                                                'marmites_props.py'))

needs_ds = pytest.mark.skipif(not os.path.exists(PARAM),
                              reason='the La Mata dataset is not present')


def _cfg_from_file():
    zones, horizons = cfgmod.soil_tables_from_param_file(PARAM)
    return cfgmod.RunConfig.from_dict({'soil': {'zone': zones,
                                                'horizon': horizons}})


def _legacy():
    import lamata_model
    # La Mata's model description as the run builds it -- no
    # parameter file (lamata_model derives the outcrop layer too)
    cMF = lamata_model.lamata_cmf()
    return cMF.cPROCESS.inputSoilParam(
        SOILparam_fn=os.path.join('MF_ws', 'inputSOILparam.txt'), NSOIL=3)


@needs_ds
def test_the_tables_give_mmsoil_what_the_file_gave_it():
    """THE acceptance test. Every array, every zone, every horizon."""
    got = props.soil_parameters(_cfg_from_file(), nsoil=3)
    want = _legacy()
    names = ('nsl', 'name', 'st', 'slprop', 'Sm', 'Sfc', 'Sr', 'Si', 'Ks')
    for n, a, b in zip(names, got, want):
        assert a == b, '%s differs: tables %r, file %r' % (n, a, b)


@needs_ds
def test_la_mata_is_three_zones_of_two_two_and_one_horizons():
    nsl, name, st, *_ = props.soil_parameters(_cfg_from_file())
    assert nsl == [2, 2, 1]
    assert name == ['SOIL1 - alluvium', 'SOIL2 - regolith', 'SOIL3 - outcrop']
    assert st == ['sandy loam field', 'sandy loam field', 'sand']


@needs_ds
def test_an_old_file_is_imported_not_defaulted(tmp_path):
    """A config that still names inputSOILparam.txt and has no tables must
    not fall back to the one-zone default -- it would run, silently, on a
    soil that is not La Mata's."""
    cfg = cfgmod.RunConfig.from_dict(
        {'paths': {'case': 'LaMata'},
         'soil': {'params': 'MF_ws/inputSOILparam.txt'}})
    assert len(cfg.soil.zone) == 3, 'the default replaced the real column'
    assert len(cfg.soil.horizon) == 5
    assert any('soil.params' in m and 'imported' in m for m in cfg.migrated)
    assert not hasattr(cfg.soil, 'params')


def test_an_old_file_that_is_not_there_is_refused():
    with pytest.raises(cfgmod.ConfigError) as exc:
        cfgmod.RunConfig.from_dict(
            {'soil': {'params': 'C:/nowhere/inputSOILparam.txt'}})
    assert 'not there' in str(exc.value)


def test_the_tables_are_not_overridden_by_a_stale_params():
    """Once migrated, the tables are the truth; a leftover params is dropped."""
    cfg = cfgmod.RunConfig.from_dict(
        {'soil': {'params': 'C:/nowhere/x.txt',
                  'zone': [{'name': 'a', 'type': 't'}],
                  'horizon': [{'zone': 1, 'slprop': 1.0}]}})
    assert [z.name for z in cfg.soil.zone] == ['a']


def test_the_default_column_is_valid_and_one_zone():
    cfg = cfgmod.RunConfig.from_dict({})
    assert len(cfg.soil.zone) == 1 and len(cfg.soil.horizon) == 1
    assert cfg.problems() == []


@pytest.mark.parametrize('change, needle', [
    ({'smax': 0.1, 'sfc': 0.2}, 'smax > sfc > sr'),
    ({'si': 0.9}, 'smax >= si >= sr'),
    ({'ks': 0.0}, 'ks must be > 0'),
    ({'slprop': 0.4}, 'sum to'),
    ({'zone': 2}, 'is not a soil zone'),
])
def test_the_rules_the_legacy_reader_stopped_the_run_on(change, needle):
    """inputSoilParam stopped the RUN on these, after the model was built.
    They are refused at the save now, on the panel."""
    h = {'zone': 1, 'slprop': 1.0, 'smax': 0.3, 'sfc': 0.2, 'sr': 0.05,
         'si': 0.2, 'ks': 1.0}
    h.update(change)
    # Refused as the configuration is BUILT: from_dict validates.
    with pytest.raises(cfgmod.ConfigError) as exc:
        cfgmod.RunConfig.from_dict(
            {'soil': {'zone': [{'name': 'a'}], 'horizon': [h]}})
    assert needle in str(exc.value), str(exc.value)


def test_a_zone_with_no_horizon_is_refused():
    with pytest.raises(cfgmod.ConfigError) as exc:
        cfgmod.RunConfig.from_dict(
            {'soil': {'zone': [{'name': 'a'}, {'name': 'b'}],
                      'horizon': [{'zone': 1, 'slprop': 1.0}]}})
    assert 'has no horizon' in str(exc.value)


def test_more_zones_than_pe_series_is_refused():
    cfg = cfgmod.RunConfig.from_dict(
        {'soil': {'zone': [{'name': 'a'}, {'name': 'b'}],
                  'horizon': [{'zone': 1}, {'zone': 2}]}})
    with pytest.raises(props.PropertyError) as exc:
        props.soil_parameters(cfg, nsoil=1)
    assert 'no evaporation' in str(exc.value)


def test_the_run_no_longer_hard_codes_the_file():
    src = io.open(os.path.join(HERE, 'run_lamata_mf6.py'),
                  encoding='utf-8').read()
    assert 'props.soil_parameters(cfg' in src
    body = src[src.index('props.soil_parameters(cfg'):]
    # The legacy call survives only in the no-configuration branch.
    assert body.index('else:') < body.index('inputSoilParam(')
