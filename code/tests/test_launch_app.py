# -*- coding: utf-8 -*-
"""The front-end runs from a local mirror of the code; the checkout keeps
everything else (user, 2026-10-09).

Streamlit does not run from a mapped network drive, so the code used to be
copied by hand to C: -- and the configuration the panels saved stayed THERE,
out of the repository. code/tools/launch_app.py now mirrors the code at
every launch and marks the mirror, and mm_paths sends every configuration,
setting and dataset path back to the checkout.
"""
import importlib.util
import os
import shutil

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


la = _load('_launch_app', os.path.join(CODE, 'tools', 'launch_app.py'))


def _write(path, text='x'):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding='utf-8')


@pytest.fixture
def repo(tmp_path):
    r = tmp_path / 'checkout'
    _write(r / 'code' / 'a.py', 'a = 1\n')
    _write(r / 'code' / 'sub' / 'b.py', 'b = 2\n')
    _write(r / 'code' / 'configs' / 'lamata.toml', '[meta]\n')
    _write(r / 'code' / 'configs' / 'paths.local.toml', '[paths]\n')
    _write(r / 'code' / '__pycache__' / 'a.cpython-312.pyc', 'bytecode')
    _write(r / '.streamlit' / 'config.toml', '[server]\n')
    _write(r / 'example' / 'LaMata' / 'inputDEM.asc', 'grid')
    _write(r / '.git' / 'HEAD', 'ref')
    return r


def _tree(root):
    return sorted(str(p.relative_to(root)).replace(os.sep, '/')
                  for p in root.rglob('*') if p.is_file())


def test_the_mirror_receives_the_code_and_nothing_else(repo, tmp_path):
    m = tmp_path / 'mirror'
    copied, removed = la.sync(repo, m)
    assert _tree(m) == ['.streamlit/config.toml', 'code/a.py', 'code/sub/b.py']
    assert (copied, removed) == (3, 0)
    assert la.sync(repo, m) == (0, 0), 'an unchanged checkout copies nothing'


def test_a_change_and_a_deletion_reach_the_mirror(repo, tmp_path):
    m = tmp_path / 'mirror'
    la.sync(repo, m)
    before = _tree(repo)
    _write(repo / 'code' / 'a.py', 'a = 10  # changed\n')
    (repo / 'code' / 'sub' / 'b.py').unlink()
    _write(m / 'code' / 'configs' / 'lamata.toml', 'a stale copy')
    _write(m / 'code' / '__pycache__' / 'a.cpython-312.pyc', 'cache')
    copied, removed = la.sync(repo, m)
    assert (m / 'code' / 'a.py').read_text() == 'a = 10  # changed\n'
    assert not (m / 'code' / 'sub').exists(), 'the emptied folder goes too'
    assert not (m / 'code' / 'configs').exists(), \
        'no configuration may live in the mirror'
    assert (m / 'code' / '__pycache__').exists(), 'caches are left alone'
    assert copied == 1 and removed == 2
    assert _tree(repo) == [p for p in before if p != 'code/sub/b.py'], \
        'the checkout is never written'


def test_the_mark_names_the_checkout_and_survives_a_refresh(repo, tmp_path):
    import mm_paths
    m = tmp_path / 'mirror'
    la.sync(repo, m)
    fn = la.mark(repo, m)
    assert fn.name == mm_paths.MIRROR_MARK
    assert fn.read_text(encoding='utf-8') == str(repo)
    la.sync(repo, m)
    assert fn.exists()


def test_a_folder_the_refresh_could_damage_is_refused(repo, tmp_path):
    with pytest.raises(SystemExit, match='overlaps'):
        la.check_mirror_dir(repo / 'code', repo)
    with pytest.raises(SystemExit, match='overlaps'):
        la.check_mirror_dir(repo, repo)
    other = tmp_path / 'other_checkout'
    _write(other / '.git' / 'HEAD', 'ref')
    with pytest.raises(SystemExit, match='git checkout'):
        la.check_mirror_dir(other, repo)
    la.check_mirror_dir(tmp_path / 'mirror', repo)          # a folder of its own


def test_a_unc_path_is_a_network_drive(tmp_path):
    assert la.on_network_drive(r'\\server\share\MARMITES')
    assert la.on_network_drive('//server/share/MARMITES')


def test_mm_paths_in_a_mirror_points_at_the_checkout(repo, tmp_path,
                                                     monkeypatch):
    """Configurations, machine settings and the default dataset are the
    CHECKOUT's; only the code runs from the mirror."""
    for env in ('MM_EXAMPLE_ROOT',):
        monkeypatch.delenv(env, raising=False)
    m = tmp_path / 'mirror'
    (m / 'code').mkdir(parents=True)
    shutil.copy2(os.path.join(CODE, 'mm_paths.py'), m / 'code' / 'mm_paths.py')
    alone = _load('_mmp_alone', str(m / 'code' / 'mm_paths.py'))
    assert not alone.MIRROR and alone.CONFIG_DIR == m / 'code' / 'configs'

    la.mark(repo, m)
    mp = _load('_mmp_mirror', str(m / 'code' / 'mm_paths.py'))
    assert mp.MIRROR and mp.REPO == repo
    assert mp.CODE_DIR == (m / 'code').resolve()
    assert mp.CONFIG_DIR == repo / 'code' / 'configs'
    assert mp.SETTINGS == repo / 'code' / 'configs' / 'paths.local.toml'
    assert mp.EXAMPLE_ROOT == repo / 'example'
    _write(repo / 'code' / 'configs' / 'paths.local.toml',
           '[paths]\nws_root = "%s"\n'
           % str(tmp_path / 'ws').replace('\\', '\\\\'))
    monkeypatch.delenv('MM_WS_ROOT', raising=False)
    monkeypatch.delenv('MARMITES_WS_ROOT', raising=False)
    mp = _load('_mmp_mirror2', str(m / 'code' / 'mm_paths.py'))
    assert mp.WS_ROOT == tmp_path / 'ws', 'the settings are read in the checkout'


def test_the_app_keeps_its_configurations_in_the_checkout():
    """panelui lists and saves them in mm_paths.CONFIG_DIR, and the machine
    settings beside them are not offered as a run configuration."""
    src = open(os.path.join(CODE, 'app', 'lib', 'panelui.py'),
               encoding='utf-8').read()
    assert 'CONFIG_DIR = str(mm_paths.CONFIG_DIR)' in src
    assert "os.path.join(CODE, 'configs')" not in src
    assert 'f != mm_paths.SETTINGS.name' in src
