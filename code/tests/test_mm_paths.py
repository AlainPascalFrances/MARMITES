# -*- coding: utf-8 -*-
"""WP1d -- the machine-specific paths, and where they come from.

Three sources with a fixed precedence: an ``MM_*`` environment variable, then
``code/configs/paths.local.toml`` (what panel 0 writes), then the built-in
default. What these tests guard is the precedence itself -- a settings file
that quietly beat the environment would break a PEST worker, which sets the
environment precisely to override one run.

The module is imported FRESH in each test: its roots are module-level
constants, so a test that re-resolved them in place would leak into the next.
"""

import importlib.util
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
MODULE = os.path.join(CODE, 'mm_paths.py')


def _fresh(monkeypatch, settings_path=None, **env):
    """Import mm_paths with a given environment and settings file."""
    for k in list(os.environ):
        if k.startswith(('MM_', 'MARMITES_')):
            monkeypatch.delenv(k, raising=False)
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    spec = importlib.util.spec_from_file_location('mm_paths_probe', MODULE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if settings_path is not None:
        mod.SETTINGS = settings_path
        mod.reload_paths()
    return mod


def test_the_defaults_resolve_without_anything_set(monkeypatch, tmp_path):
    mp = _fresh(monkeypatch, settings_path=tmp_path / 'absent.toml')
    assert str(mp.DATA_ROOT)
    assert mp.EXAMPLE_ROOT == mp.REPO / 'example'
    assert mp.GIS == mp.DATA_ROOT / 'GIS'
    assert mp.source_of('data_root') == 'the built-in default'


def test_a_missing_or_broken_settings_file_falls_back(monkeypatch, tmp_path):
    """It is read on every import, so it must never be able to stop one."""
    bad = tmp_path / 'broken.toml'
    bad.write_text('this is not = = toml [[[', encoding='utf-8')
    mp = _fresh(monkeypatch, settings_path=bad)
    assert mp.read_settings() == {}
    assert str(mp.DATA_ROOT)


def test_panel_zero_settings_are_read(monkeypatch, tmp_path):
    f = tmp_path / 'paths.local.toml'
    f.write_text('[paths]\ndata_root = "%s"\nws_root = "%s"\n'
                 % (str(tmp_path / 'data').replace('\\', '\\\\'),
                    str(tmp_path / 'ws').replace('\\', '\\\\')),
                 encoding='utf-8')
    mp = _fresh(monkeypatch, settings_path=f)
    assert mp.DATA_ROOT == tmp_path / 'data'
    assert mp.WS_ROOT == tmp_path / 'ws'
    assert mp.GIS == tmp_path / 'data' / 'GIS'     # still derived
    assert mp.source_of('ws_root') == 'code/configs/paths.local.toml'


def test_the_environment_beats_the_settings_file(monkeypatch, tmp_path):
    """A PEST worker sets MM_WS_ROOT to override ONE run. If the file it
    happens to find on disk won, that override would silently do nothing."""
    f = tmp_path / 'paths.local.toml'
    f.write_text('[paths]\nws_root = "%s"\n'
                 % str(tmp_path / 'from_file').replace('\\', '\\\\'),
                 encoding='utf-8')
    mp = _fresh(monkeypatch, settings_path=f,
                MM_WS_ROOT=str(tmp_path / 'from_env'))
    assert mp.WS_ROOT == tmp_path / 'from_env'
    assert 'environment variable' in mp.source_of('ws_root')


def test_saving_writes_only_what_was_given_and_applies_it(monkeypatch, tmp_path):
    mp = _fresh(monkeypatch, settings_path=tmp_path / 'p.toml')
    mp.save_settings({'ws_root': str(tmp_path / 'out'), 'data_root': '',
                      'gis': str(tmp_path / 'shp')})
    text = (tmp_path / 'p.toml').read_text(encoding='utf-8')
    assert 'ws_root' in text and 'gis' in text
    assert 'data_root' not in text          # blank -> back to the default
    assert mp.WS_ROOT == tmp_path / 'out'   # applied without a restart
    assert mp.GIS == tmp_path / 'shp'


def test_a_windows_path_survives_the_round_trip(monkeypatch, tmp_path):
    r"""Backslashes are why panel 0 does not write Python source: 'E:\next'
    read back as an escape would be a different folder, or a syntax error."""
    mp = _fresh(monkeypatch, settings_path=tmp_path / 'p.toml')
    win = 'E:\\00code_ws\\next\\tab\\new'
    mp.save_settings({'data_root': win})
    assert mp.read_settings()['data_root'] == win
    assert str(mp.DATA_ROOT) == win


def test_the_dataset_follows_the_example_root(monkeypatch, tmp_path):
    mp = _fresh(monkeypatch, settings_path=tmp_path / 'p.toml')
    assert mp.dataset_dir('LaMata') == mp.REPO / 'example' / 'LaMata'
    mp.save_settings({'example_root': str(tmp_path / 'cases')})
    assert mp.dataset_dir('LaMata') == tmp_path / 'cases' / 'LaMata'


def test_the_report_points_at_panel_zero_not_at_a_file_to_edit(monkeypatch,
                                                               tmp_path):
    """The whole point of putting these on panel 0: nobody should be told to
    go and edit a source file."""
    import io
    mp = _fresh(monkeypatch, settings_path=tmp_path / 'p.toml',
                MM_DATA_ROOT=str(tmp_path / 'nope'))
    buf = io.StringIO()
    mp.report_paths('LaMata', stream=buf)
    out = buf.getvalue()
    assert 'MISSING' in out
    assert 'edit code/mm_paths.py' not in out
    assert 'PANEL 0' in out


@pytest.mark.parametrize('key', ['example_root', 'data_root', 'gis', 'ws_root'])
def test_every_folder_panel_zero_offers_is_documented(monkeypatch, tmp_path,
                                                      key):
    mp = _fresh(monkeypatch, settings_path=tmp_path / 'p.toml')
    env, _default, doc = mp.SETTABLE[key]
    assert env.startswith('MM_')
    assert len(doc) > 30, '%s has no explanation for the panel to show' % key
