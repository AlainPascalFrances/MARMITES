# -*- coding: utf-8 -*-
"""Saved state belongs to its MESH, not only to its grid kind -- and the
panel picks it from what is on disk.

2026-09-27: the workspace held ``voronoi_lamata`` saved on the 15 915-cell
stream-refined mesh while the model had moved to the 2126-cell one. Its
scope sidecar named the grid KIND, the layer count and the case -- all the
same -- so the panel would have shown "The saved state belongs to this grid
and layer set" and the run would have handed MF6 an array of the wrong
length. The sidecar now records the mesh (signature and shape), and every
state's own ESRI header is compared too, which catches states saved before.
"""
import importlib.util
import io
import json
import os
import sys
import time
import types

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for p in (CODE, os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'app')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg        # noqa: E402
import marmites_props as props        # noqa: E402

GRID = {'shape': (2126, 1), 'signature': 'fca2d233303fc153'}


def _asc(path, nrow, ncol=1):
    with io.open(path, 'w', encoding='utf-8') as fh:
        fh.write('ncols %d\nnrows %d\nxllcorner 0\nyllcorner 0\n'
                 'cellsize 1\nnodata_value -9999\n' % (ncol, nrow))
        fh.write('\n'.join('1.0' for _ in range(nrow * ncol)) + '\n')


def _cfg(**spinup):
    cfg = mcfg.RunConfig.from_dict({})
    cfg.grid.kind = 'voronoi'
    cfg.layers.nlay = 2
    for k, v in spinup.items():
        setattr(cfg.spinup, k, v)
    return cfg


def _state(d, name, cfg, nrow, grid=None, full=True, means=False, nlay=2):
    """A state as the run writes it: heads per layer, the rest, the
    sidecar -- ``grid`` absent for one saved before 2026-09-27."""
    for k in range(nlay):
        _asc(os.path.join(d, '%s_l%d.asc' % (name, k + 1)), nrow)
    if full:
        open(os.path.join(d, name + '_state.npz'), 'wb').close()
    if means:
        _asc(os.path.join(d, name + '_perc.asc'), nrow)
        _asc(os.path.join(d, name + '_etg.asc'), nrow)
    side = {'state_hash': cfg.state_hash(), 'scope': cfg.state_scope()}
    if grid:
        side['grid'] = {'shape': list(grid['shape']),
                        'signature': grid['signature']}
    with io.open(mcfg.state_sidecar(d, name), 'w', encoding='utf-8') as fh:
        json.dump(side, fh, default=str)


# ---------------------------------------------------------- the question

def test_a_state_from_another_mesh_of_the_same_kind_is_caught(tmp_path):
    """The case that passed: same kind, same layers, same case -- another
    mesh. The sidecar predates the grid record, so the HEADER says it."""
    d = str(tmp_path)
    cfg = _cfg(strt_heads='voronoi_lamata')
    _state(d, 'voronoi_lamata', cfg, 15915)
    assert mcfg.state_problem(cfg, d) == '', 'the old check (no grid)'
    why = mcfg.state_problem(cfg, d, grid=GRID)
    assert '15915' in why and '2126' in why, why


def test_the_recorded_mesh_signature_is_compared(tmp_path):
    """Two meshes can share a cell count; the signature cannot."""
    d = str(tmp_path)
    cfg = _cfg(strt_heads='other')
    _state(d, 'other', cfg, 2126,
           grid={'shape': (2126, 1), 'signature': '0123456789abcdef'})
    why = mcfg.state_problem(cfg, d, grid=GRID)
    assert 'another voronoi mesh' in why and GRID['signature'] in why, why


def test_a_state_from_this_mesh_fits(tmp_path):
    d = str(tmp_path)
    cfg = _cfg(strt_heads='mine')
    _state(d, 'mine', cfg, 2126, grid=GRID)
    assert mcfg.state_problem(cfg, d, grid=GRID) == ''


def test_the_question_can_be_asked_of_one_field(tmp_path):
    """The heads must not fall back because the MEANS do not fit -- a run
    from saved heads does not even use the means."""
    d = str(tmp_path)
    cfg = _cfg(strt_heads='mine', steady_means='old')
    _state(d, 'mine', cfg, 2126, grid=GRID)
    _state(d, 'old', cfg, 15915, means=True)
    assert mcfg.state_problem(cfg, d, grid=GRID,
                              keys=('spinup.strt_heads',)) == ''
    assert '15915' in mcfg.state_problem(cfg, d, grid=GRID,
                                         keys=('spinup.steady_means',))


def test_the_run_starts_from_the_land_surface_on_another_mesh(tmp_path):
    d = str(tmp_path)
    cfg = _cfg(strt_heads='voronoi_lamata')
    _state(d, 'voronoi_lamata', cfg, 15915)
    kind, _p, why = props.resolve_initial_heads(cfg, d, verbose=False,
                                                grid=GRID)
    assert kind == 'dem' and '15915' in why
    _state(d, 'mine', cfg, 2126, grid=GRID)
    cfg.spinup.strt_heads = 'mine'
    kind, payload, _why = props.resolve_initial_heads(cfg, d, verbose=False,
                                                      grid=GRID)
    assert (kind, payload) == ('saved', 'mine')


# ------------------------------------------------------ what is on disk

def test_the_saved_states_are_listed_newest_first(tmp_path):
    d = str(tmp_path)
    cfg = _cfg()
    _state(d, 'older', cfg, 2126, grid=GRID, full=False)
    past = time.time() - 3600
    for fn in os.listdir(d):
        os.utime(os.path.join(d, fn), (past, past))
    _state(d, 'newer', cfg, 2126, grid=GRID, means=True)
    _asc(os.path.join(d, 'broken_l1.asc'), 2126)      # layer 2 missing
    got = mcfg.saved_states(d, 2)
    assert [s['name'] for s in got] == ['newer', 'older']
    assert got[0]['full'] and not got[1]['full'], 'heads-only not flagged'
    assert got[0]['shape'] == (2126, 1)
    means = mcfg.saved_states(d, 2, 'spinup.steady_means')
    assert [s['name'] for s in means] == ['newer']


def test_no_workspace_lists_nothing(tmp_path):
    assert mcfg.saved_states(str(tmp_path / 'absent'), 2) == []


# ------------------------------------------------------------ the run

def _runner_mod():
    name = '_runner_state_grid'
    if name in sys.modules:
        return sys.modules[name]
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, 'run_lamata_mf6.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_the_run_records_the_mesh_beside_the_state(tmp_path):
    r = _runner_mod()
    cfg = _cfg()
    a = types.SimpleNamespace(state_dir=str(tmp_path))
    cmf = types.SimpleNamespace(nrow=2126, ncol=1,
                                mesh_signature=GRID['signature'])
    r._write_state_scope(a, cfg, 'mine', grid=r._run_grid(cmf))
    side = json.load(open(mcfg.state_sidecar(str(tmp_path), 'mine')))
    assert side['grid'] == {'shape': [2126, 1],
                            'signature': GRID['signature']}
    assert side['state_hash'] == cfg.state_hash()
    # a structured model has no mesh signature, and says so
    assert r._run_grid(types.SimpleNamespace(nrow=65, ncol=60)) == {
        'shape': (65, 60), 'signature': None}


def test_the_driver_passes_the_grid_everywhere_it_decides():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert "cMF.mesh_signature = info['signature']" in src
    assert src.count('grid=_grid') >= 4, 'heads, means, and both saves'
    assert 'resolve_initial_heads(cfg, a.state_dir,' in src


# ------------------------------------------------ the validation panel

def _checks():
    spec = importlib.util.spec_from_file_location(
        '_checks_state_grid', os.path.join(CODE, 'app', 'lib', 'checks.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _mesh_cache(d, ncpl, signature):
    m = os.path.join(d, '_mesh')
    os.makedirs(m, exist_ok=True)
    json.dump({'ncpl': ncpl}, open(os.path.join(m, 'mesh_voronoi.json'), 'w'))
    json.dump({'signature': signature, 'ncpl': ncpl},
              open(os.path.join(m, 'mesh_voronoi.sig.json'), 'w'))


def test_validation_warns_about_a_state_from_another_mesh(tmp_path):
    d = str(tmp_path)
    cfg = _cfg(strt_heads='voronoi_lamata')
    _state(d, 'voronoi_lamata', cfg, 15915)
    _mesh_cache(d, 2126, GRID['signature'])
    got = list(_checks().check_state_scope(cfg, workspace=d))
    assert len(got) == 1 and 'will not be used' in got[0].title
    assert '15915' in got[0].detail and 'land surface' in got[0].detail


def test_validation_is_silent_on_a_state_that_fits(tmp_path):
    d = str(tmp_path)
    cfg = _cfg(strt_heads='mine', steady_means='old')
    _state(d, 'mine', cfg, 2126, grid=GRID)
    _state(d, 'old', cfg, 15915, means=True)     # unused: heads are named
    _mesh_cache(d, 2126, GRID['signature'])
    assert list(_checks().check_state_scope(cfg, workspace=d)) == []


# -------------------------------------------------------------- the panel

def test_the_panel_picks_the_state_and_asks_one_name_to_save_it(tmp_path):
    pytest.importorskip('streamlit', reason='streamlit not installed')
    from streamlit.testing.v1 import AppTest
    import shutil

    d = str(tmp_path).replace('\\', '/')
    ref = os.path.join(CODE, 'configs', 'lamata.toml')
    tmp = os.path.join(CODE, 'configs', '_stategridtest.toml')
    shutil.copyfile(ref, tmp)
    text = io.open(tmp, encoding='utf-8').read()
    out, here = [], None
    for raw in text.splitlines(True):
        s = raw.strip()
        if s.startswith('[') and s.endswith(']'):
            here = s[1:-1]
        elif here == 'paths' and s.startswith('ws '):
            raw = 'ws = "%s"\n' % d
        elif here == 'spinup' and s.split(' ')[0] in (
                'strt_heads', 'steady_means', 'save_strt', 'save_means'):
            raw = '%s = ""\n' % s.split(' ')[0]
        elif here == 'spinup' and s.startswith('cycles '):
            raw = 'cycles = 6\n'
        out.append(raw)
    io.open(tmp, 'w', encoding='utf-8', newline='').write(''.join(out))
    try:
        cfg = mcfg.load_run_config(tmp)
        assert cfg.paths.ws == d, 'the scratch file did not take paths.ws'
        _state(d, 'voronoi_lamata', cfg, 15915, nlay=int(cfg.layers.nlay))
        _state(d, 'coarse', cfg, 2126, grid=GRID, means=True,
               nlay=int(cfg.layers.nlay))
        _mesh_cache(d, 2126, GRID['signature'])
        at = AppTest.from_file(os.path.join(
            CODE, 'app', 'pages', '4_4_-_Unsaturated_zone_and_groundwater.py'),
            default_timeout=180)
        at.session_state['config_file'] = os.path.basename(tmp)
        at.run()
        assert not at.exception, [str(e.value) for e in at.exception]
        boxes = {s.key: s for s in at.selectbox if s.key}
        assert 'spinup.strt_heads' in boxes, 'still a text box'
        assert 'spinup.steady_means' in boxes
        labels = list(boxes['spinup.strt_heads'].options)
        assert any(l.startswith('coarse') and '2126 cells' in l
                   and l.endswith('✓') for l in labels), labels
        assert any(l.startswith('voronoi_lamata') and '15915 cells' in l
                   and 'does not fit' in l for l in labels), labels
        assert labels[0].startswith('none'), labels
        texts = {t.key for t in at.text_input if t.key}
        assert 'spinup.save_strt' in texts
        assert 'spinup.save_means' not in texts, 'two names asked again'
        # a cold start with cycles > 1 is a spin-up, and is said to be one
        assert any('This run is a spin-up' in str(i.value) for i in at.info)
        # picking the state from the other mesh says so, with both counts
        at.selectbox(key='spinup.strt_heads').select('voronoi_lamata').run()
        warned = ' '.join(str(w.value) for w in at.warning)
        assert 'will not be used' in warned and '15915' in warned, warned
        # ... and the one from this mesh is accepted
        at.selectbox(key='spinup.strt_heads').select('coarse').run()
        assert any('belongs to this grid' in str(s.value) for s in at.success)
        # the means are not asked of a run from saved heads
        assert at.selectbox(key='spinup.steady_means').disabled
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_the_saved_state_keeps_the_ponds_stages(tmp_path):
    """A run from a saved state starts each pond where the spin-up left it;
    a state saved with no LAK carries none."""
    r = _runner_mod()
    import numpy as np
    b = types.SimpleNamespace(nlay=2, nrow=5, ncol=1,
                              thti_from_wc=lambda wc: wc)
    st = types.SimpleNamespace(Ssoil_ini=np.zeros((5, 3)))
    ctx = types.SimpleNamespace(ncell=5, _nslmax=3)
    pref = str(tmp_path / 'mine')
    cpl = types.SimpleNamespace(uzf_wc_final=None, carry_out=None,
                                lak_stage_final=np.array([731.2, 745.9]))
    r.save_run_state(pref, b, cpl, st)
    got = r.load_run_state(pref, b, ctx)
    assert np.allclose(got['lak_stage'], [731.2, 745.9])
    cpl.lak_stage_final = None
    r.save_run_state(pref, b, cpl, st)
    assert r.load_run_state(pref, b, ctx)['lak_stage'] is None


def test_the_driver_hands_the_stages_on():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert "b.lak_strt_carry = saved_state.get('lak_stage')" in src
    assert "b.lak_strt_carry = getattr(cpl, 'lak_stage_final', None)" in src
    cpl = open(os.path.join(CODE, 'marmites_coupler.py'), encoding='utf-8').read()
    assert "('XNEWPAK', f'{name}/LAK')" in cpl
