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
    # set at construction, in props.model_from_config (2026-10-07)
    assert 'props.model_from_config(' in src
    prp = open(os.path.join(HERE, '..', 'ppMF6', 'marmites_props.py'),
               encoding='utf-8').read()
    body = prp[prp.index('def model_from_config('):]
    body = body[:body.index('\ndef ')]
    assert 'modelname=cfg.meta.model_name(cfg.paths.case)' in body


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
    than leave it to be discovered: a new name starts a new file."""
    page = open(HOME, encoding='utf-8').read()
    assert 'A new model name starts a NEW configuration' in page
    assert 'The configuration file is named after the model' in page


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
    # (on the next run: a drawn widget's own key cannot be assigned)
    assert panelui.st.session_state['__config_next'] == 'newname.toml'


# ------------------------------- a new model name starts a new configuration

def _pending(panelui, monkeypatch, edits):
    monkeypatch.setattr(panelui.st, 'session_state',
                        {'__pending_edits': dict(edits)}, raising=False)


def test_a_new_model_name_writes_a_new_file_and_keeps_the_old(tmp_path,
                                                              monkeypatch):
    """CdL is started from lamata.toml (2026-11-12): typing its name on the
    Home panel and saving gives cdl.toml, and La Mata's file is untouched."""
    panelui = _panelui()
    src = tmp_path / 'lamata.toml'
    body = open(REF, encoding='utf-8').read()
    src.write_text(body, encoding='utf-8')
    cfg = cfgmod.load_run_config(str(src))
    _pending(panelui, monkeypatch, {'meta.model': 'CdL', 'paths.case': 'CdL'})

    applied, _digest, ok = panelui.save_now(cfg, str(src))

    assert ok, _digest
    new = tmp_path / 'cdl.toml'
    assert new.exists()
    assert src.read_text(encoding='utf-8') == body, 'the old file changed'
    got = cfgmod.load_run_config(str(new))
    assert got.meta.model == 'CdL' and got.paths.case == 'CdL'
    state = panelui.st.session_state
    assert state['__config_next'] == 'cdl.toml'
    assert 'cdl.toml' in state['__saved_as'] and 'unchanged' in state['__saved_as']
    assert not state.get('__pending_edits')


def test_a_new_name_another_configuration_has_is_refused(tmp_path,
                                                         monkeypatch):
    panelui = _panelui()
    src = tmp_path / 'lamata.toml'
    body = open(REF, encoding='utf-8').read()
    src.write_text(body, encoding='utf-8')
    other = tmp_path / 'cdl.toml'
    other.write_text('# another model entirely\n', encoding='utf-8')
    cfg = cfgmod.load_run_config(str(src))
    _pending(panelui, monkeypatch, {'meta.model': 'cdl'})

    applied, why, ok = panelui.save_now(cfg, str(src))

    assert not ok and 'another' in why and 'cdl.toml' in why
    assert src.read_text(encoding='utf-8') == body
    assert other.read_text(encoding='utf-8').startswith('# another model')
    # nothing is lost: the edits are still remembered
    assert panelui.st.session_state['__pending_edits'] == {'meta.model': 'cdl'}


def test_a_save_without_a_new_name_goes_to_its_own_file(tmp_path,
                                                        monkeypatch):
    panelui = _panelui()
    src = tmp_path / 'lamata.toml'
    src.write_text(open(REF, encoding='utf-8').read(), encoding='utf-8')
    cfg = cfgmod.load_run_config(str(src))
    for edits in ({'run.nsp': 7}, {'meta.model': cfg.meta.model},
                  {'meta.model': 'LAMATA'}):
        assert panelui.save_path(cfg, str(src), edits) == (str(src), '')


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
