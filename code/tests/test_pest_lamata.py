# -*- coding: utf-8 -*-
"""The La Mata PEST++ tools (WP7): MODIS composites, the targets and their
collection, the scalars a forward run writes, the setup's own helpers. On
small synthetic files -- no model run, no MODIS download needed."""

import datetime as dt
import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in (TRUNK, os.path.join(TRUNK, 'tools')):
    if p not in sys.path:
        sys.path.insert(0, p)

import pest_lamata_obs as po                                    # noqa: E402
import pest_lamata_forward as pfw                               # noqa: E402
import pest_lamata_setup as pst                                 # noqa: E402
import pest_modis_et as pme                                     # noqa: E402
import marmites_config as mcfg                                  # noqa: E402


# ------------------------------------------------------------------ MODIS

def test_a_composite_covers_eight_days_and_the_last_ends_the_year():
    assert pme.composite_period(2009, 1) == (dt.date(2009, 1, 1), 8)
    assert pme.composite_period(2010, 361) == (dt.date(2010, 12, 27), 5)
    assert pme.composite_period(2008, 361) == (dt.date(2008, 12, 26), 6)


def test_a_composite_downloaded_twice_is_read_once(tmp_path):
    for d in ('a', 'b'):
        (tmp_path / d).mkdir()
    for d, doy in (('a', 321), ('a', 329), ('b', 329), ('b', 337)):
        (tmp_path / d / ('MOD16A2GF.A2010%03d.h17v04.061.2021.hdf' % doy)
         ).write_bytes(b'')
    got = pme.list_composites(str(tmp_path))
    assert [c[1] for c in got] == [321, 329, 337]


def test_pixels_weigh_by_the_area_they_share_with_the_catchment():
    from rasterio.transform import Affine
    from shapely.geometry import box
    tr = Affine(10.0, 0.0, 0.0, 0.0, -10.0, 30.0)       # 3 x 3 pixels of 10
    r, c, a = pme.pixel_weights(tr, (3, 3), box(5.0, 5.0, 15.0, 25.0))
    got = {(int(i), int(j)): float(s) for i, j, s in zip(r, c, a)}
    assert got == {(0, 0): 25.0, (0, 1): 25.0, (1, 0): 50.0, (1, 1): 50.0,
                   (2, 0): 25.0, (2, 1): 25.0}


def test_fills_are_left_out_of_the_mean():
    et, frac = pme.weighted_mean([100, 200, 32767], np.array([1.0, 1.0, 2.0]))
    assert et == pytest.approx(15.0) and frac == pytest.approx(0.5)
    et, frac = pme.weighted_mean([32761], np.array([1.0]))
    assert np.isnan(et) and frac == 0.0


# --------------------------------------------------------------- targets

def _exports(tmp_path, n=40, start='2009-09-20'):
    d = tmp_path / 'ex'
    d.mkdir()
    dates = pd.date_range(start, periods=n).strftime('%Y-%m-%d')
    per = np.arange(1, n + 1)
    obs_h = np.where(per % 10 == 0, 700.0 + per * 0.01, np.nan)
    p0 = pd.DataFrame({'obsnme': ['h_p0_%04d' % k for k in per],
                       'point': 'P0', 'layer': 1, 'period': per,
                       'date': dates, 'sim_head_m': 699.0 + per * 0.01,
                       'sim_hcorr_m': np.nan, 'obs_head_m': obs_h})
    # a soft-data well (the synthetic quarterly cycle): every 20 days here
    c1 = p0.assign(obsnme=['h_c1_%04d' % k for k in per], point='C1',
                   sim_head_m=750.0,
                   obs_head_m=np.where(per % 20 == 0, 751.0, np.nan))
    pd.concat([p0, c1]).to_csv(d / 'obs_heads.csv', index=False)
    rows = []
    for pt in ('P0', 'SM', 'O1'):
        for sl in (1, 2):
            rows.append(pd.DataFrame({
                'obsnme': 'x', 'point': pt, 'soil_layer': sl, 'period': per,
                'date': dates, 'sim_theta': 0.15,
                'obs_theta': np.where(per % 5 == 0, 0.17, np.nan)}))
    pd.concat(rows).to_csv(d / 'obs_sm.csv', index=False)
    pd.DataFrame({'obsnme': 'x', 'site': 'outlet', 'period': per,
                  'date': dates, 'sim_q_m3d': 0.0,
                  'sim_q_mmd': np.linspace(0.0, 2.0, n),
                  'obs_q_mmd': np.where(per > 5, 0.5, np.nan),
                  'obs_q_m3d': np.nan}).to_csv(d / 'obs_sfr.csv', index=False)
    rows = []
    for site in ('P0', 'EC', 'catchment'):
        rows.append(pd.DataFrame({'obsnme': 'x', 'site': site, 'period': per,
                                  'date': dates, 'sim_et_total_mmd': 1.0 + 0
                                  * per, 'obs_et_mmd': np.nan}))
    pd.concat(rows).to_csv(d / 'obs_et.csv', index=False)
    return str(d)


def test_targets_are_built_weighted_and_collected(tmp_path):
    ex = _exports(tmp_path)
    modis = pd.DataFrame({'start': ['2009-09-22', '2009-10-02', '2010-12-01'],
                          'end': ['2009-09-29', '2009-10-09', '2010-12-08'],
                          'et_mm': [9.0, 10.0, 3.0], 'year': [2009] * 3,
                          'doy': [265, 275, 335]})
    ec = pd.DataFrame({'date': [dt.date(2009, 10, 5)], 'et_mmd': [1.2]})
    tg = po.build_targets(ex, dt.date(2009, 10, 1), modis=modis, ec=ec)
    g = tg.groupby('obgnme')['obsnme'].size().to_dict()
    assert g['hds_p0'] == 4 and g['sm_p0'] == 16 and g['sm_sm'] == 16
    assert g['hds_c'] == 2
    assert 'O1' not in ''.join(tg['key'])              # no probe there
    assert g['et_rs'] == 2                             # 3rd not in the run
    assert g['et_ec'] == 1 and tg.loc[tg['obgnme'] == 'et_ec',
                                      'weight'].iloc[0] == 0.0
    # the warm-up is kept, unweighted
    w = tg.set_index('obsnme')['weight']
    assert w['h_p0_0010'] == 0.0 and w['h_p0_0020'] == 1.0
    assert w['etm_2009265'] == 0.0 and w['etm_2009275'] == 1.0
    sim = po.collect(ex, tg)
    v = sim.set_index('obsnme')['value']
    assert v['h_p0_0020'] == pytest.approx(699.2)
    assert v['etm_2009275'] == pytest.approx(8.0)       # 8 days x 1 mm/d
    q20 = np.log10(np.linspace(0.0, 2.0, 40)[19] + po.Q_EPS)
    assert v['q_0020'] == pytest.approx(q20)
    tg['weight'] = po.balance_weights(tg, sim, total=100.0)
    r = (tg['obsval'] - sim['value']) * tg['weight']
    phi = (r ** 2).groupby(tg['obgnme']).sum()
    assert phi['hds_p0'] == pytest.approx(35.0)       # the heads' share
    assert phi['et_rs'] == pytest.approx(10.0)
    assert phi['et_ec'] == 0.0
    # soft data (user, 2026-10-12): a C-well head weighs 25 % of a P0 head
    w = tg.set_index('obsnme')['weight']
    assert w['h_c1_0020'] == pytest.approx(0.25 * w['h_p0_0020'])
    assert w['h_c1_0040'] == pytest.approx(0.25 * w['h_p0_0020'])
    sd = tg.set_index('obsnme')['standard_deviation']
    assert sd['h_c1_0020'] == pytest.approx(4 * sd['h_p0_0020'])


def test_a_value_the_run_did_not_make_is_a_sentinel_not_a_gap(tmp_path):
    ex = _exports(tmp_path)
    tg = po.build_targets(ex, dt.date(2009, 10, 1))
    tg.loc[0, 'key'] = 'P0|999'
    sim = po.collect(ex, tg)
    assert sim['value'].iloc[0] == 1.0e10 and np.isfinite(sim['value']).all()


# --------------------------------------------------------------- scalars

def _cfg():
    cfg = mcfg.RunConfig.from_dict({})
    cfg.soil.horizon = [mcfg.SoilHorizon(zone=1, smax=0.15, sfc=0.125,
                                         sr=0.075, si=0.1, ks=1.5)]
    return cfg


def test_the_scalars_reach_the_configuration_and_the_layer_files(tmp_path):
    cfg = _cfg()
    t = pst.scalar_table(cfg)
    assert t.set_index('name')['value']['sfc_f_h1'] == pytest.approx(2 / 3)
    t.loc[t['name'] == 'sfc_f_h1', 'value'] = 0.5
    t.loc[t['name'] == 'ks_h1', 'value'] = 4.0
    t.loc[t['name'] == 'ss_l1', 'value'] = 0.3
    t.loc[t['name'] == 'crr_beta', 'value'] = 0.9
    pfw.apply_scalars(cfg, t, mult_dir=str(tmp_path))
    h = cfg.soil.horizon[0]
    assert h.sfc == pytest.approx(0.075 + 0.5 * 0.075) and h.ks == 4.0
    assert cfg.crr.beta == 0.9
    assert float(open(tmp_path / 'ss_l1.txt').read()) == 0.3
    assert os.path.exists(tmp_path / 'drn_l1.txt')


def test_the_bounds_are_the_drafts():
    assert pst.scalar_bounds('sfc_f_h1', 0.5) == ('none', 0.4, 0.6)
    assert pst.scalar_bounds('ks_h2', 2.0) == ('log', pytest.approx(0.2),
                                               20.0)
    assert pst.scalar_bounds('drn_l2', 1.0) == ('log', 0.1, 100.0)
    assert pst.scalar_bounds('ss_l1', 1.0) == ('log', 0.1, 10.0)
    assert pst.scalar_bounds('uzf_eps', 3.5)[1:] == (3.5, 8.0)
    # capped (2026-10-12): at 1 m/d the model crawls and fits no better
    assert pst.scalar_bounds('sfr_rhk', 0.1) == ('log', 0.001, 0.5)


def test_the_prior_starts_inside_the_ranges():
    """eps and beta start inside their ranges: the first prior piled 64 %
    and 57 % of its draws on their bounds."""
    t = pst.scalar_table(_cfg()).set_index('name')['value']
    assert t['uzf_eps'] == 4.5 and t['crr_beta'] == 0.95


def test_the_setup_switches_the_extra_points_on(tmp_path):
    f = tmp_path / 'inputObs.txt'
    f.write_text('# Name X Y lay hi h0 RC STO\n'
                 'P0 739525 4555875 1 2.00 739.6 50 0.025   # Ptb2\n'
                 '##W1 741177 4554666 1 2.00 780.3 50 0.025   #Pvs0\n'
                 '#C4 740764 4555342 1 2.00 767.2 50 0.025\n'
                 '##H3 740730 4553470 1 2.00 807.1 50 0.025\n',
                 encoding='utf-8')
    assert sorted(pst._obs_table_with(str(f))) == ['C4', 'W1']
    txt = f.read_text(encoding='utf-8')
    assert '\nW1 741177' in txt and '\nC4 740764' in txt
    assert '##H3' in txt and txt.startswith('# Name')


def test_the_workspace_may_not_be_inside_ws_root_or_the_repository(tmp_path):
    from types import SimpleNamespace
    mp = SimpleNamespace(WS_ROOT=str(tmp_path / 'ws'), REPO=str(tmp_path /
                                                                'repo'))
    with pytest.raises(SystemExit):
        pst._check_ws(str(tmp_path / 'ws' / 'pest'), mp)
    with pytest.raises(SystemExit):
        pst._check_ws('', mp)
    pst._check_ws(str(tmp_path / 'pest'), mp)
