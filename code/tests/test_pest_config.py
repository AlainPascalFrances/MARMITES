# -*- coding: utf-8 -*-
"""The configuration's [pest] section (WP7 front-end, phase 2): what the
PEST++ set-up page edits, and the rules that keep a calibration PEST++ can
run -- refused at save, not hours into a set-up."""

import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
if TRUNK not in sys.path:
    sys.path.insert(0, TRUNK)

import marmites_config as mcfg                                  # noqa: E402


def _cfg(nlay=2):
    cfg = mcfg.RunConfig.from_dict({})
    cfg.layers.nlay = nlay
    cfg.soil.horizon = [mcfg.SoilHorizon(zone=1, smax=0.30, sfc=0.20,
                                         sr=0.05, si=0.1, ks=1.5)]
    return cfg


def _pest_problems(cfg):
    return [e for e in cfg.problems() if e.startswith('pest')]


def _with(cfg, *rows):
    cfg.pest.par = [mcfg.PestPar(**r) for r in rows]
    return _pest_problems(cfg)


def test_an_empty_section_is_the_calibration_of_2026_10_18():
    p = mcfg.Pest()
    assert (p.pp_space, p.vario_range, p.k_order_decades) == (400.0, 1500.0,
                                                              1.0)
    assert (p.nreal, p.noptmax, p.agents, p.start) == (100, 3, 20, 'prior')
    assert (p.bad_phi_sigma, p.overdue_giveup_fac) == (2.5, 4.0)
    assert p.points_on == 'C4, C5, W1' and p.soft_fraction == 0.25
    assert p.par == []                       # the defaults, at set-up
    share = {g.name: g.share for g in p.group}
    assert share == {'hds_p0': 0.35, 'hds_c': 0.0, 'sm_p0': 0.175,
                     'sm_sm': 0.175, 'q_log': 0.10, 'q_dry': 0.05,
                     'q_year': 0.05, 'et_rs': 0.10, 'et_pt': 0.0}
    assert _pest_problems(_cfg()) == []


def test_the_section_survives_the_round_trip(tmp_path):
    cfg = _cfg()
    cfg.pest.par = [mcfg.PestPar('k', 1, True, 0.05, 2.3e-5, 1.12, True),
                    mcfg.PestPar('uzf_eps', 0, False, 4.5, 3.5, 8.0, False)]
    cfg.pest.base_state = 'hi_voronoi_lamata_evt_drn1'
    cfg.pest.base_run = '20261011214903_2lay_fullperiod_numba'
    cfg.pest.start, cfg.pest.start_from = 'warm', 'x/lamata.1.par.jcb'
    f = str(tmp_path / 'c.toml')
    cfg.write_toml(f)
    back = mcfg.load_run_config(f)
    assert back.pest.par == cfg.pest.par
    assert back.pest.group == cfg.pest.group
    assert back.pest.base_state == 'hi_voronoi_lamata_evt_drn1'
    assert back.pest.base_run == '20261011214903_2lay_fullperiod_numba'
    assert (back.pest.start, back.pest.start_from) == ('warm',
                                                       'x/lamata.1.par.jcb')


def test_the_retired_keys_are_reported_and_dropped(tmp_path):
    f = tmp_path / 'c.toml'
    f.write_text('[pest]\nuse_pest_params = false\n'
                 'params_npz = "pest_optimised.npz"\nmultipliers = ""\n')
    cfg = mcfg.load_run_config(str(f))
    said = ' '.join(cfg.migrated)
    assert 'pest.use_pest_params is gone' in said
    assert 'pest.params_npz is gone' in said
    assert not hasattr(cfg.pest, 'use_pest_params')


@pytest.mark.parametrize('row, frag', [
    (dict(family='kk', layer=1), 'no such family'),
    (dict(family='k', layer=3), 'no such layer (1..2)'),
    (dict(family='ks', layer=2), 'no such horizon (1..1)'),
    (dict(family='uzf_eps', layer=1, start=4.5, lower=3.5, upper=8.0),
     'layer must be 0'),
    (dict(family='k', layer=1, start=1.0, lower=2.0, upper=10.0),
     'start 1 outside'),
    (dict(family='k', layer=1, start=1.0, lower=1.0, upper=1.0),
     'must be below upper'),
    (dict(family='drn', layer=1, start=0.0, lower=0.0, upper=1.0),
     'needs lower > 0'),
    (dict(family='sy', layer=1, start=0.5, lower=0.1, upper=1.0),
     'a water content'),
    (dict(family='uzf_eps', layer=0, start=3.0, lower=2.0, upper=8.0,
          log=False), '[3.5, 14]'),
    (dict(family='crr_beta', layer=0, start=1.0, lower=0.8, upper=1.2,
          log=False), 'a share'),
])
def test_a_range_pest_cannot_sample_is_refused(row, frag):
    errs = _with(_cfg(), row)
    assert any(frag in e for e in errs), errs


def test_a_row_given_twice_is_refused():
    r = dict(family='ss', layer=1)
    assert any('given twice' in e for e in _with(_cfg(), r, r))


def test_the_field_capacity_stays_between_wilting_point_and_porosity():
    """User, 2026-10-17: the soil field capacity is bounded by the
    porosity and the wilting point (residual) -- the range cannot leave
    them, whatever the Soil panel's values become."""
    ok = dict(family='sfc', layer=1, start=0.2, lower=0.1, upper=0.25,
              log=False)
    assert _with(_cfg(), ok) == []
    for lo, hi in ((0.05, 0.25), (0.1, 0.30), (0.04, 0.35)):
        errs = _with(_cfg(), dict(ok, lower=lo, upper=hi))
        assert any('between the horizon\'s wilting point sr 0.05 and '
                   'porosity smax 0.3' in e for e in errs), errs


def test_uzf_porosity_sy_plus_sr_stays_below_one():
    errs = _with(_cfg(), dict(family='sy', layer=1, start=0.3, lower=0.1,
                              upper=0.5),
                 dict(family='sr', layer=1, start=0.3, lower=0.1,
                      upper=0.5))
    assert any('Sy + Sr' in e for e in errs), errs


def test_a_k_order_no_pair_can_hold_is_refused():
    k1 = dict(family='k', layer=1, start=1.0, lower=1.0, upper=10.0)
    k2 = dict(family='k', layer=2, start=5.0, lower=1.0, upper=5.0)
    assert any('decade' in e for e in _with(_cfg(), k1, k2))
    cfg = _cfg()
    cfg.pest.k_order_decades = 0.0                # no order asked
    assert _with(cfg, k1, k2) == []


def test_the_shares_make_phi_and_soft_data_take_none():
    cfg = _cfg()
    cfg.pest.group[0].share = 0.30                         # hds_p0
    assert any('sum to 0.95' in e for e in _pest_problems(cfg))
    cfg = _cfg()
    cfg.pest.group[1].share = 0.1                          # hds_c, soft
    assert any('soft data' in e for e in _pest_problems(cfg))
    cfg = _cfg()
    cfg.pest.group = cfg.pest.group[:-1]                   # no et_pt
    assert any('no row for et_pt' in e for e in _pest_problems(cfg))
    cfg = _cfg()
    cfg.pest.group[2].noise = 0.0
    assert any('noise > 0' in e for e in _pest_problems(cfg))


@pytest.mark.parametrize('key, value, frag', [
    ('calib_start', '1 Oct 2009', 'YYYY-MM-DD'),
    ('pp_space', 0.0, 'pp_space must be > 0'),
    ('k_order_decades', -1.0, 'k_order_decades must be >= 0'),
    ('nreal', 1, 'at least 2'),
    ('noptmax', -2, 'noptmax must be >= -1'),
    ('lambda_mults', '0.1, x', 'lambda_mults'),
    ('lambda_mults', '0.1, -1', 'lambda_mults'),
    ('overdue_giveup_fac', 1.0, 'overdue_giveup_fac must be > 1'),
    ('agents', 0, 'agents must be >= 1'),
    ('start', 'resume', 'pest.start must be one of'),
    ('start', 'warm', 'needs pest.start_from'),
    ('start', 'restart', 'master folder'),
])
def test_the_ensemble_settings_are_guarded(key, value, frag):
    cfg = _cfg()
    setattr(cfg.pest, key, value)
    errs = _pest_problems(cfg)
    assert any(frag in e for e in errs), errs


def test_the_set_up_page_edits_go_through_the_one_save(tmp_path):
    """A field and a table of [pest], applied as the panels apply them
    (lib.editor, the sidebar's save): validated, then written."""
    sys.path.insert(0, os.path.join(TRUNK, 'app'))
    try:
        from lib import editor
    finally:
        sys.path.pop(0)
    cfg = _cfg()
    f = str(tmp_path / 'c.toml')
    rows = [dict(family='sr', layer=1, enable=True, start=0.275,
                 lower=0.045, upper=0.30, log=True)]
    applied, _h = editor.save(cfg, f, {'pest.nreal': 50,
                                       'pest.lambda_mults': '0.5, 1.0'},
                              {'pest.par': rows})
    assert 'pest.par = 1 row(s)' in applied
    back = mcfg.load_run_config(f)
    assert back.pest.nreal == 50 and back.pest.lambda_mults == '0.5, 1.0'
    assert back.pest.par[0].family == 'sr'
    bad = [dict(rows[0], upper=0.02)]
    with pytest.raises(mcfg.ConfigError, match='pest.par'):
        editor.save(cfg, f, {}, {'pest.par': bad})
