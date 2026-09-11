# -*- coding: utf-8 -*-
"""WP1d -- MMsurf driven by the configuration, and the handover file removed.

``__inputMMsurf4MMsoil.txt`` was written by MMsurf and read back by the
driver, positionally, and it WON: with MMsurf not running, editing Zr or kT*
in the ini changed nothing. These tests pin the replacement -- one source of
truth, a shape check that refuses a stale or truncated forcing file, and the
kt_s inversion happening in exactly one place.
"""

import importlib.util
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
REPO = os.path.abspath(os.path.join(CODE, '..'))
DS = os.path.join(REPO, 'example', 'LaMata')
REF_CONFIG = os.path.join(CODE, 'configs', 'lamata.toml')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


msurf = _load('marmites_surface_t', os.path.join(CODE, 'marmites_surface.py'))
cfgmod = _load('marmites_config_t', os.path.join(CODE, 'marmites_config.py'))


@pytest.fixture
def cfg():
    if not os.path.exists(REF_CONFIG):
        pytest.skip('reference configuration not present')
    return cfgmod.load_run_config(REF_CONFIG)


# ------------------------------------------------------- the handover is gone

def test_the_handover_file_is_not_in_the_dataset():
    assert not os.path.exists(os.path.join(DS, '__inputMMsurf4MMsoil.txt')), \
        'the handover file is back -- WP1d removed it deliberately'


def test_forcing_spec_carries_what_the_handover_used_to(cfg):
    spec = msurf.forcing_spec(cfg, DS)
    assert (spec.nmeteo, spec.nveg, spec.nsoil) == (1, 3, 3)
    assert (spec.ncrop, spec.nfield) == (1, 1)
    assert spec.veg_names == ['grassMU', 'Qilex', 'Qpyr']
    assert spec.Zr == [0.4, 15.0, 10.0]
    # The old driver did kT_s = 1.0 / <file value>. The configuration holds
    # the slope, so the value arrives ready to use.
    assert all(0.0 < s <= 1.0 for s in spec.kT_s)


def test_forcing_file_names_are_the_committed_ones(cfg):
    spec = msurf.forcing_spec(cfg, DS)
    assert spec.files['rf_veg'] == 'inputZONRF_veg_d.txt'
    assert spec.files['tf_veg'] == 'inputZONTF_veg_d.txt'
    assert os.path.isabs(spec.path('rf_veg'))


# ------------------------------------------------------------- the shape check

def test_block_counts_match_the_committed_forcing(cfg):
    """Nothing in the files says how many blocks they hold; this does."""
    spec = msurf.forcing_spec(cfg, DS)
    assert spec.blocks('rf_veg') == 1          # per meteo zone
    assert spec.blocks('tf_veg') == 3          # meteo x vegetation
    assert spec.blocks('pt_veg') == 3
    assert spec.blocks('lai_veg') == 3         # vegetation only
    assert spec.blocks('pe') == 3              # meteo x surface soil
    assert spec.blocks('eo') == 1


def test_the_committed_forcing_passes_the_shape_check(cfg):
    spec = msurf.forcing_spec(cfg, DS)
    ndays = msurf.check_forcing(spec, must_exist=True)
    assert ndays == 1949


def test_a_truncated_forcing_file_is_refused(cfg, tmp_path):
    """The committed inputZONRFe_veg_d.txt held 4869 values against a 1949-day
    record -- not a whole number of blocks -- and went unnoticed for years."""
    ws = str(tmp_path)
    spec = msurf.forcing_spec(cfg, ws)
    for role in spec.roles():
        n = 3 if role == 'date' else spec.blocks(role) * 3
        body = ('2008-01-01 00:00, 1\n' * 3 if role == 'date'
                else '0.0\n' * n)
        with open(spec.path(role), 'w') as fh:
            fh.write('#\n' + body)
    assert msurf.check_forcing(spec) == 3
    with open(spec.path('tf_veg'), 'w') as fh:
        fh.write('#\n' + '0.0\n' * 7)          # not a whole number of blocks
    with pytest.raises(msurf.MMsurfError) as e:
        msurf.check_forcing(spec)
    assert 'expected' in str(e.value)


def test_a_missing_forcing_file_says_what_to_do(cfg, tmp_path):
    spec = msurf.forcing_spec(cfg, str(tmp_path))
    with pytest.raises(msurf.MMsurfError) as e:
        msurf.check_forcing(spec, must_exist=True)
    msg = str(e.value)
    assert 'run.surface' in msg and 'inputDATE.txt' in msg


def test_a_forcing_shorter_than_the_run_is_refused(cfg):
    spec = msurf.forcing_spec(cfg, DS)
    with pytest.raises(msurf.MMsurfError) as e:
        msurf.check_forcing(spec, must_exist=True, nper=100000)
    assert 'stress period' in str(e.value)


# ------------------------------------------------- the generated parameter file

def test_generated_par_file_round_trips_through_the_app_parser(cfg, tmp_path):
    """The file MMsurf parses must be readable by the front-end's own reader,
    which is the independent check that the layout is right."""
    msc = _load('mmsurf_cfg_t', os.path.join(CODE, 'mmsurf_config.py'))
    out = str(tmp_path / 'gen.ini')
    msurf.write_par_file(cfg, out, config_hash='deadbeef')
    ms = msc.load_mmsurf_ini(out)
    assert ms.counts == {'NMETEO': 1, 'NVEG': 3, 'NCRP': 1, 'NFIELD': 1,
                         'NSOIL': 3}
    assert [n for n, _v in ms.veg] == ['grassMU', 'Qilex', 'Qpyr']
    assert [v['Zr'] for _n, v in ms.veg] == [0.4, 15.0, 10.0]
    assert [n for n, _v in ms.soil] == ['alluvium', 'regolith', 'outcrop']


def test_kt_s_is_inverted_exactly_once(cfg, tmp_path):
    """The configuration holds s; the FILE holds 1/s. If both ends agreed the
    model would silently run with the reciprocal of the intended slope."""
    out = str(tmp_path / 'gen.ini')
    msurf.write_par_file(cfg, out)
    rows = [ln.split() for ln in open(out) if ln.startswith('Qilex')]
    assert rows, 'Qilex row not written'
    written = float(rows[0][-1])
    s = cfg.surface.vegetation[1].kt_s
    assert written == pytest.approx(1.0 / s)
    assert written > 1.0 and s < 1.0


def test_the_generated_file_names_its_provenance(cfg, tmp_path):
    out = str(tmp_path / 'gen.ini')
    msurf.write_par_file(cfg, out, config_hash='abc123')
    head = open(out).read(400)
    assert 'GENERATED' in head and 'abc123' in head


# ------------------------------------------------------ the station conversion

def test_station_is_converted_from_project_crs_to_geographic(cfg):
    pytest.importorskip('pyproj')
    phi, lm = msurf.geographic(cfg, cfg.surface.station[0])
    # La Mata sits near 41.12 N, 6.15 W; the old ini said 41.045 / 6.16,
    # which is 6.8 km south of the catchment.
    assert 41.0 < phi < 41.3
    assert 6.0 < lm < 6.3
    assert abs(phi - 41.045) > 0.05, 'this is the old ini value, not the real one'


def test_a_missing_crs_is_an_error_not_a_guess(cfg):
    cfg.grid.crs_epsg = 0
    with pytest.raises(msurf.MMsurfError) as e:
        msurf.geographic(cfg, cfg.surface.station[0])
    assert 'crs_epsg' in str(e.value)


def test_surface_ws_is_outside_the_repository(cfg):
    ws = msurf.surface_ws(cfg, r'E:\00code_ws', 'LaMata')
    assert 'LaMata_MMsurf' in ws
    assert os.path.abspath(REPO) not in os.path.abspath(ws)
