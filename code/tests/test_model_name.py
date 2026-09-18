# -*- coding: utf-8 -*-
"""The model has a name, asked once on Overview.

It used to live in the MODFLOW parameter file as "lamataMM": the model was
named in a file the modeller never opened, and MODFLOW 6 names every file
it writes after it.
"""

import importlib.util
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
REF = os.path.join(CODE, 'configs', 'lamata.toml')
HOME = os.path.join(CODE, 'app', 'Home.py')

for _p in (CODE, os.path.join(CODE, 'app')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


cfgmod = _load('marmites_config_mn', os.path.join(CODE, 'marmites_config.py'))
from lib import schema                                        # noqa: E402


@pytest.fixture
def cfg():
    return cfgmod.load_run_config(REF)


def test_the_model_is_named_in_the_configuration(cfg):
    assert cfg.meta.model == 'lamata'
    assert cfg.meta.model_name(cfg.paths.case) == 'lamata'


def test_blank_falls_back_to_the_case(cfg):
    """A new configuration should not need this answered before anything
    else works; the case name is a reasonable stand-in."""
    cfg.meta.model = ''
    assert cfg.meta.model_name(cfg.paths.case) == 'lamata'
    fresh = cfgmod.RunConfig.from_dict({})
    assert fresh.meta.model_name(fresh.paths.case) == 'lamata'


@pytest.mark.parametrize('bad,frag', [
    ('la mata', 'letters'),            # a space
    ('la-mata', 'letters'),            # a hyphen
    ('2lamata', 'letters'),            # starts with a digit
    ('averyverylongmodelname', '16'),  # 22 characters
])
def test_a_name_modflow_cannot_use_is_refused_here(cfg, bad, frag):
    """MODFLOW 6 names its files after the model, so an unusable name is
    otherwise discovered only once the input has been written."""
    cfg.meta.model = bad
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert frag in str(e.value)


def test_the_run_takes_the_name_from_the_panel():
    """A field the run ignores is decoration -- and this one decides what
    every output file is called."""
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert 'cfg.meta.model_name(' in src
    assert 'cMF.modelname = _model' in src


def test_overview_asks_for_it_before_anything_else():
    page = open(HOME, encoding='utf-8').read()
    assert "'meta.model'" in page
    assert (page.index('What this model is called')
            < page.index('What this configuration will run')), (
        'the name is asked after the things named after it')


def test_the_configuration_file_is_not_renamed_behind_the_modeller(cfg):
    """A file may be open in an editor, named in a launch command, or read
    by a run still going. The panel says the two disagree and stops there."""
    page = open(HOME, encoding='utf-8').read()
    assert 'nothing is renamed automatically' in page
    for verb in ('os.rename(', 'shutil.move(', 'os.remove('):
        assert verb not in page, '%s on the Overview panel' % verb


def test_the_field_says_what_it_becomes():
    help_ = schema.describe('meta.model')[2]
    assert '.hds' in help_ and '16' in help_
