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
    ('la.mata', 'letters'),            # a dot; the underscore IS allowed
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
    # only the NAME is given: the case study is the model name, not asked
    _pending(panelui, monkeypatch, {'meta.model': 'CdL'})

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
    assert panelui.st.session_state['__pending_edits']['meta.model'] == 'cdl'


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


# ------------------------------------- the folders go through the same save
# (user, 2026-11-13: no button of their own, no second way)

def _own_settings(panelui, monkeypatch, tmp_path):
    """mm_paths on a settings file of its own, restored afterwards."""
    mp = panelui.mm_paths
    for env in ('MM_RS_WS', 'MM_PEST_WS', 'MM_SM_WS'):
        monkeypatch.delenv(env, raising=False)
    monkeypatch.setattr(mp, 'SETTINGS', tmp_path / 'paths.local.toml')
    monkeypatch.setattr(mp, 'CASE', '')
    mp.reload_paths()
    return mp


def test_the_folders_are_written_by_the_one_save_for_the_case(tmp_path,
                                                              monkeypatch):
    panelui = _panelui()
    mp = _own_settings(panelui, monkeypatch, tmp_path)
    try:
        src = tmp_path / 'lamata.toml'
        body = open(REF, encoding='utf-8').read()
        src.write_text(body, encoding='utf-8')
        cfg = cfgmod.load_run_config(str(src))
        boxes = {k: panelui._path_now(k) for k in mp.SETTABLE}
        monkeypatch.setattr(panelui.st, 'session_state', {}, raising=False)
        panelui.remember_paths(boxes)
        assert panelui.unsaved_paths() == [], 'untouched boxes are no change'
        assert panelui.unsaved_changes(cfg) == []

        boxes['rs'] = str(tmp_path / 'cdl_rs')
        boxes['pest_ws'] = str(tmp_path / 'cdl_pest')
        panelui.remember_paths(boxes)
        panelui.remember({'meta.model': 'CdL'})
        todo = panelui.unsaved_changes(cfg)
        assert any(t.startswith('folder rs = ') for t in todo), todo
        assert any(t.startswith('folder pest_ws = ') for t in todo), todo

        applied, why, ok = panelui.save_now(cfg, str(src))
        assert ok, why
        assert any(a.startswith('folder rs') for a in applied)
        # for the case of the configuration that was written, and in force
        assert mp.read_case_settings('CdL')['rs'] == str(tmp_path / 'cdl_rs')
        # (under the MODEL's name, whatever its letter case)
        assert mp.CASE == 'cdl' and mp.RS == tmp_path / 'cdl_rs'
        assert mp.read_case_settings('cdl') == mp.read_case_settings('CdL')
        assert mp.PEST_WS == tmp_path / 'cdl_pest'
        assert (tmp_path / 'cdl.toml').exists()
        assert src.read_text(encoding='utf-8') == body
        assert not panelui.pending_paths() and not panelui.pending()
    finally:
        monkeypatch.undo()
        mp.reload_paths()


def test_a_refused_save_writes_no_folder(tmp_path, monkeypatch):
    panelui = _panelui()
    mp = _own_settings(panelui, monkeypatch, tmp_path)
    try:
        src = tmp_path / 'lamata.toml'
        src.write_text(open(REF, encoding='utf-8').read(), encoding='utf-8')
        # panelui's own marmites_config, so the refusal is ITS ConfigError
        cfg = panelui.mcfg.load_run_config(str(src))
        boxes = {k: panelui._path_now(k) for k in mp.SETTABLE}
        boxes['rs'] = str(tmp_path / 'other_rs')
        monkeypatch.setattr(panelui.st, 'session_state', {}, raising=False)
        panelui.remember_paths(boxes)
        panelui.remember({'run.nsp': -5})                  # invalid
        applied, why, ok = panelui.save_now(cfg, str(src))
        assert not ok and 'run.nsp' in why
        assert not mp.SETTINGS.exists()
        assert panelui.pending_paths()['rs'] == str(tmp_path / 'other_rs')
    finally:
        monkeypatch.undo()
        mp.reload_paths()


def test_the_case_study_is_the_model_name(tmp_path, monkeypatch):
    """User, 2026-11-13: not asked. A new name is a new case; an existing
    configuration keeps its spelling when only the letter case differs
    (model lamata, dataset LaMata)."""
    panelui = _panelui()
    cfg = cfgmod.load_run_config(REF)
    assert cfg.meta.model == 'lamata' and cfg.paths.case == 'LaMata'
    assert panelui.case_of(cfg, {}) == 'LaMata'
    assert panelui.case_of(cfg, {'meta.model': 'LAMATA'}) == 'LaMata'
    assert panelui.case_of(cfg, {'meta.model': 'CdL'}) == 'CdL'
    assert panelui.case_of(cfg, {'meta.model': ''}) == 'LaMata'
    page = open(HOME, encoding='utf-8').read()
    assert "'paths.case'" not in page, 'the case study is still asked on Home'
    help_ = schema.describe('meta.model')[2]
    assert 'underscores' in help_ and 'folders' in help_
    # an underscore is a valid name (user, 2026-11-13)
    cfg.meta.model = 'cdl_v2'
    assert not [m for m in cfg.problems() if 'meta.model' in m]
    assert panelui.model_of(cfg, {}) == 'cdl_v2'
    assert panelui.model_of(cfg, {'meta.model': 'CdL'}) == 'cdl'


def test_a_new_model_takes_the_folders_on_screen_as_its_own(tmp_path,
                                                            monkeypatch):
    """Unchanged boxes included: a model without a set of its own would
    read whatever was saved before models had theirs."""
    panelui = _panelui()
    mp = _own_settings(panelui, monkeypatch, tmp_path)
    try:
        src = tmp_path / 'lamata.toml'
        src.write_text(open(REF, encoding='utf-8').read(), encoding='utf-8')
        cfg = panelui.mcfg.load_run_config(str(src))
        boxes = {k: panelui._path_now(k) for k in mp.SETTABLE}
        monkeypatch.setattr(panelui.st, 'session_state', {}, raising=False)
        panelui.remember_paths(boxes)
        panelui.remember({'meta.model': 'cdl_v2'})
        assert panelui.unsaved_paths() == []
        applied, why, ok = panelui.save_now(cfg, str(src))
        assert ok, why
        assert (tmp_path / 'cdl_v2.toml').exists()
        assert mp.CASE == 'cdl_v2'
        table = mp.read_case_settings('cdl_v2')
        # the environment pins some roots in the tests: whatever is not
        # pinned and not blank was stored
        free = [k for k in mp.CASE_KEYS if boxes.get(k)
                and 'environment' not in mp.source_of(k)]
        assert free and all(table.get(k) == boxes[k] for k in free), table
    finally:
        monkeypatch.undo()
        mp.reload_paths()


def test_the_dataset_is_one_of_the_models_folders(tmp_path, monkeypatch):
    """User, 2026-11-13: the folder with the data need not be called
    after the model -- it is named on Home with the others."""
    import importlib.util as _u
    spec = _u.spec_from_file_location('mm_paths_ds',
                                      os.path.join(CODE, 'mm_paths.py'))
    for k in list(os.environ):
        if k.startswith(('MM_', 'MARMITES_')):
            monkeypatch.delenv(k, raising=False)
    mp = _u.module_from_spec(spec)
    spec.loader.exec_module(mp)
    mp.SETTINGS = tmp_path / 'p.toml'
    mp.reload_paths()
    mp.save_settings({'dataset': str(tmp_path / 'anywhere'),
                      'example_root': str(tmp_path / 'ex')}, case='cdl_v2')
    mp.use_case('cdl_v2', 'CdL')
    assert mp.dataset_dir('CdL') == tmp_path / 'anywhere'
    assert mp.active_dataset() == str(tmp_path / 'anywhere')
    # another case is not this model's: the examples root, as before
    assert mp.dataset_dir('LaMata') == tmp_path / 'ex' / 'LaMata'
    # a model with no dataset of its own: <examples root>/<case>
    mp.use_case('lamata', 'LaMata')
    assert mp.dataset_dir('LaMata') == tmp_path / 'ex' / 'LaMata'
    # a process that moves the examples root (a PEST worker, the tests)
    # is not reached around
    mp.use_case('cdl_v2', 'CdL')
    monkeypatch.setenv('MM_EXAMPLE_ROOT', str(tmp_path / 'worker'))
    mp.reload_paths()
    assert mp.dataset_dir('CdL') == tmp_path / 'worker' / 'CdL'


def test_home_has_one_way_to_save():
    page = open(HOME, encoding='utf-8').read()
    assert 'save_paths' not in page and 'Save these paths' not in page
    assert 'selectbox(' not in page, 'a selectbox of its own on Home'
    assert 'remember_paths' in page
