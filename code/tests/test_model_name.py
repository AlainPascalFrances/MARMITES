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
    """It is the first field on the page, before the things named after
    it. No heading: the field's own label says what it is."""
    page = open(HOME, encoding='utf-8').read()
    assert "'meta.model'" in page
    assert (page.index("'meta.model'")
            < page.index('What this configuration will run')), (
        'the name is asked after the things named after it')


def test_the_page_says_the_file_follows_the_name():
    """Automatic, so the panel has to say it is about to happen rather
    than leave a rename to be discovered."""
    page = open(HOME, encoding='utf-8').read()
    assert 'renames this file to match' in page


def test_the_field_says_what_it_becomes():
    help_ = schema.describe('meta.model')[2]
    assert '.hds' in help_ and '16' in help_


# ------------------------------------------- the file follows the model name

def _panelui():
    """panelui imports streamlit; skip where it is not installed."""
    pytest.importorskip('streamlit')
    from lib import panelui
    return panelui


def test_saving_renames_the_file_after_the_model(tmp_path, monkeypatch):
    """The model's name IS its identity: MODFLOW writes <model>.hds, and the
    configuration describing it should not be called something else."""
    panelui = _panelui()
    src = tmp_path / 'old_name.toml'
    src.write_text(open(REF, encoding='utf-8').read(), encoding='utf-8')
    cfg = cfgmod.load_run_config(str(src))
    cfg.meta.model = 'newname'

    monkeypatch.setattr(panelui.st, 'session_state', {}, raising=False)
    new, why = panelui.rename_to_model(cfg, str(src))

    assert new and os.path.basename(new) == 'newname.toml'
    assert os.path.exists(new) and not os.path.exists(str(src)), (
        'it must MOVE, not leave both'
    )
    assert 'Renamed' in why
    # the sidebar picks the file from this, so it has to follow
    assert panelui.st.session_state['config_file'] == 'newname.toml'


def test_a_name_already_taken_is_refused_rather_than_overwritten(tmp_path,
                                                                 monkeypatch):
    """That file describes a DIFFERENT model. Silently overwriting it would
    lose it, and automatic is not the same as careless."""
    panelui = _panelui()
    body = open(REF, encoding='utf-8').read()
    src = tmp_path / 'one.toml'
    src.write_text(body, encoding='utf-8')
    other = tmp_path / 'taken.toml'
    other.write_text('# another model entirely\n', encoding='utf-8')

    cfg = cfgmod.load_run_config(str(src))
    cfg.meta.model = 'taken'
    monkeypatch.setattr(panelui.st, 'session_state', {}, raising=False)
    new, why = panelui.rename_to_model(cfg, str(src))

    assert not new
    assert 'another configuration' in why
    assert os.path.exists(str(src)), 'the file being saved was lost'
    assert other.read_text(encoding='utf-8').startswith('# another model')


def test_a_matching_name_is_left_alone(tmp_path, monkeypatch):
    panelui = _panelui()
    src = tmp_path / 'lamata.toml'
    src.write_text(open(REF, encoding='utf-8').read(), encoding='utf-8')
    cfg = cfgmod.load_run_config(str(src))
    monkeypatch.setattr(panelui.st, 'session_state', {}, raising=False)
    assert panelui.rename_to_model(cfg, str(src)) == ('', '')
