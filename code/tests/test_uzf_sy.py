# -*- coding: utf-8 -*-
"""The unsaturated zone drains what the aquifer stores: thts - thtr = Sy.

The run of 2026-09-23: Sy 0.01 in the aquifer, thts - thtr = 0.40 in UZF.
A rising water table hands UZF's water above thtr to the aquifer, which
needs only Sy per metre of rise -- so every rise released ten times the
water it took, and on the first transient day the median water table went
from 6.2 m below ground to 0.26 m ABOVE it in 15,269 of 15,750 cells
(5.7e6 m3 in one step). UZF1 never did this: the ini said SPECIFYTHTR 0,
and UZF1 then derived thtr = thts - Sy. UZF6 has no such option; MODFLOW 6
only asks for consistency with STO's Sy, so the build does it.
"""

import importlib.util
import os
import sys
import types

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'app')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import marmites_config as mcfg                                 # noqa: E402
import marmites_props as props                                 # noqa: E402


# ------------------------------------------------------------ the rule
def test_sy_derives_thtr_and_raises_thti_as_for_la_mata():
    thtr, thti, notes = props.uzf_water_contents(0.05, 0.45, 0.15, 0.01)
    assert float(thtr) == pytest.approx(0.44)
    assert float(thti) == pytest.approx(0.44), 'thti below thtr not raised'
    assert any('thts - Sy' in n for n in notes)
    assert any('1 raised to thtr' in n for n in notes)


def test_a_cell_by_cell_sy_gives_a_cell_by_cell_thtr():
    sy = np.array([[0.01, 0.05], [0.10, 0.20]])
    thtr, _thti, _ = props.uzf_water_contents(0.05, 0.45, 0.30, sy)
    assert np.allclose(0.45 - thtr, sy)


def test_sy_at_or_above_thts_is_refused():
    with pytest.raises(props.PropertyError) as exc:
        props.uzf_water_contents(0.05, 0.30, 0.15, np.array([0.1, 0.3]))
    assert 'not positive in 1 cell' in str(exc.value)


def test_inactive_cells_are_neither_touched_nor_judged():
    sy = np.array([0.01, 0.9])                     # 0.9: an inactive cell
    thtr, thti, _ = props.uzf_water_contents(
        0.05, 0.45, 0.15, sy, active=np.array([True, False]))
    assert thtr[1] == pytest.approx(0.05) and thti[1] == pytest.approx(0.15)
    assert thtr[0] == pytest.approx(0.44)


def test_source_refuses_a_zone_that_drains_more_than_sy_quoting_mf6():
    with pytest.raises(props.PropertyError) as exc:
        props.uzf_water_contents(0.05, 0.45, 0.15, 0.01, thtr_from='source')
    msg = str(exc.value)
    assert 'Storage Package' in msg and 'runs away' in msg
    assert 'uzf.thtr_from = "sy"' in msg


def test_source_accepts_a_consistent_thtr_as_given():
    thtr, thti, notes = props.uzf_water_contents(0.44, 0.45, 0.445, 0.01,
                                                 thtr_from='source')
    assert float(thtr) == pytest.approx(0.44)
    assert float(thti) == pytest.approx(0.445) and not notes


def test_a_thti_above_thts_is_lowered_and_counted():
    _thtr, thti, notes = props.uzf_water_contents(0.44, 0.45, 0.9, 0.01,
                                                  thtr_from='source')
    assert float(thti) == pytest.approx(0.45)
    assert any('1 lowered to thts' in n for n in notes)


def test_an_unknown_mode_is_refused():
    with pytest.raises(props.PropertyError):
        props.uzf_water_contents(0.05, 0.45, 0.15, 0.01, thtr_from='value')


# ------------------------------------------------------ the configuration
def test_the_default_is_the_uzf1_derivation():
    assert mcfg.RunConfig.from_dict({}).uzf.thtr_from == 'sy'


def test_lamata_derives_it():
    cfg = mcfg.load_run_config(os.path.join(CODE, 'configs', 'lamata.toml'))
    assert cfg.uzf.thtr_from == 'sy'


def test_an_unknown_thtr_from_is_a_configuration_error():
    with pytest.raises(mcfg.ConfigError) as exc:
        mcfg.RunConfig.from_dict({'uzf': {'thtr_from': 'raster'}})
    assert 'thtr_from' in str(exc.value)


def test_with_sy_a_thti_below_the_given_thtr_is_not_an_error():
    """thtr is not read then, and thti is clipped at the build."""
    cfg = mcfg.RunConfig.from_dict({'uzf': {'thtr_from': 'sy',
                                            'thtr': {'value': 0.3},
                                            'thti': {'value': 0.15}}})
    assert cfg.uzf.thtr_from == 'sy'


def test_with_source_the_old_range_rules_still_hold():
    with pytest.raises(mcfg.ConfigError):
        mcfg.RunConfig.from_dict({'uzf': {'thtr_from': 'source',
                                          'thtr': {'value': 0.3},
                                          'thti': {'value': 0.15}}})


# ------------------------------------------------------------- the panel
def test_the_question_is_on_the_uzf_tab_with_its_choices():
    from lib import schema
    rows = [k for row in schema.UZF_ROWS for k in row if k]
    assert 'uzf.thtr_from' in rows
    assert rows.index('uzf.thtr_from') < rows.index('uzf.thtr')
    src = open(os.path.join(CODE, 'app', 'lib', 'schema.py'),
               encoding='utf-8').read()
    assert "'uzf.thtr_from': lambda: ['sy', 'source']" in src
    page = open(os.path.join(CODE, 'app', 'pages',
                             '4_4_-_Unsaturated_zone_and_groundwater.py'),
                encoding='utf-8').read()
    assert "cfg.uzf.thtr_from == 'sy'" in page


# --------------------------------------------------- the validation panel
def _checks():
    name = '_mm_checks_uzf'
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(CODE, 'app', 'lib', 'checks.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _asc(path, value, nrow=2, ncol=3):
    with open(path, 'w') as fh:
        fh.write('ncols %d\nnrows %d\nxllcorner 0\nyllcorner 0\n'
                 'cellsize 10\nNODATA_value -9999\n' % (ncol, nrow))
        for _ in range(nrow):
            fh.write(' '.join(str(value) for _ in range(ncol)) + '\n')


def _cfg(tmp_path, **uzf):
    ds = tmp_path / 'ds'
    ds.mkdir()
    for L in (1, 2):
        _asc(str(ds / ('Sy_l%d.asc' % L)), 0.01)
        _asc(str(ds / ('ib_l%d.asc' % L)), 1)
    cfg = mcfg.RunConfig.from_dict({'layers': {'nlay': 2}})
    cfg.layers.sy = mcfg.VectorSource(raster='Sy_l%d.asc')
    cfg.layers.ibound = mcfg.VectorSource(raster='ib_l%d.asc')
    for k, v in uzf.items():
        setattr(cfg.uzf, k, v)
    return cfg, str(ds)


def test_the_validation_panel_stops_an_inconsistent_source(tmp_path):
    chk = _checks()
    cfg, ds = _cfg(tmp_path, thtr_from='source')
    got = list(chk.check_uzf_sy(cfg, dataset_dir=ds))
    assert [c.level for c in got] == [chk.ERROR]
    assert got[0].key == 'uzf.thtr_from' and got[0].panel == 4
    assert 'Storage Package' in got[0].detail


def test_the_validation_panel_says_what_the_derivation_will_do(tmp_path):
    chk = _checks()
    cfg, ds = _cfg(tmp_path, thtr_from='sy')
    got = list(chk.check_uzf_sy(cfg, dataset_dir=ds))
    assert got and all(c.level == chk.INFO for c in got)
    text = ' '.join(c.title for c in got)
    assert '0.44' in text and 'raised to thtr' in text


def test_the_check_is_one_of_those_the_launch_asks():
    chk = _checks()
    assert chk.check_uzf_sy in chk.CHECKS


# ------------------------------------------------------- the run warning
coupler = pytest.importorskip('marmites_coupler')


def test_the_warning_counts_every_cell_above_ground_not_only_past_1_m():
    """The run's own warning said '1176 cells, 1 of 60 stress periods':
    it only counted heads MORE than 1 m above ground."""
    top = np.zeros(4)
    heads = np.array([[1.5, 0.5, 0.2, -1.0],       # SP 1
                      [0.3, 0.3, -0.5, -1.0],      # SP 2
                      [0.2, -0.2, -0.5, -1.0]])    # SP 3
    lines = coupler.surface_excess_lines(heads, top, 4)
    text = '\n'.join(lines)
    assert 'up to 1.5 m' in text
    assert 'more than 0.1 m above: 3 of 4 cells' in text
    assert 'in 3 of 3 stress periods' in text
    assert 'at worst 3 cells at once (SP 1)' in text
    assert 'more than 1 m above:   1 cells, in 1 stress period' in text


def test_the_warning_is_silent_below_1_m():
    assert coupler.surface_excess_lines(np.full((2, 3), 0.5),
                                        np.zeros(3), 3) == []


def test_the_rejection_note_turns_mm_per_day_back_into_m3():
    """rejinf_hist is mm/d per cell; the note added it up as m3 and said
    46.1 % where MF6's UZF budget said 2.6 %."""
    src = open(os.path.join(CODE, 'marmites_coupler.py'),
               encoding='utf-8').read()
    assert 'rejected = float(np.sum(self.rejinf_hist))' not in src
    assert ('self.rejinf_hist / self.conv_fact\n'
            '                                    * self.area[None, :]') in src
    # and the arithmetic: 2 cells, 10 mm/d on 100 m2 and 1000 mm/d on 1 m2
    rej_mm = np.array([[10.0, 1000.0]])
    area = np.array([100.0, 1.0])
    assert float(np.sum(rej_mm / 1000.0 * area[None, :])) == \
        pytest.approx(1.0 + 1.0)                    # m3, not 1010


def test_the_seepage_estimate_takes_each_cell_on_its_own_area():
    src = open(os.path.join(CODE, 'marmites_coupler.py'),
               encoding='utf-8').read()
    assert 'np.nanmax(self.exf_hist) / self.conv_fact' not in src
    assert 'self.exf_hist / self.conv_fact' in src


def test_the_warning_measures_from_the_land_not_the_soil_base(capsys):
    """2026-10-04, the LAK run: 'the water table rises up to 1.4 m ABOVE the
    land surface' -- measured from the MF6 top, which is the BASE of the
    MMsoil column. The highest head stood 0.08 m above the land, under a
    full pond. A water table up in the soil column is exfiltration, not a
    flood: the guard takes the land surface."""
    import importlib.util
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(
        't_coupler_mock_land', os.path.join(here, 'test_coupler_mock.py'))
    cm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cm)
    cpl, api, ctx = cm._setup(nper=2, heads0=699.0)
    nrow, ncol = ctx.cMF.nrow, ctx.cMF.ncol
    cpl.mf6b.top = np.full((nrow, ncol), 697.5)          # soil base
    cpl.mf6b._land_surface = lambda: np.full((nrow, ncol), 699.5)
    cpl2 = coupler.MF6Coupler(cpl.mm, ctx, cpl.state, cpl.mf6b,
                              conv_fact=1000.0)
    assert np.allclose(cpl2.top_cell, 697.5)
    assert np.allclose(cpl2.land_cell, 699.5)
    cpl2.run(api)
    out = capsys.readouterr().out
    assert 'ABOVE the land surface' not in out           # 1.5 m into the soil
    # and with the land where the soil base is, it is reported
    cpl.mf6b._land_surface = lambda: np.full((nrow, ncol), 697.5)
    cpl3 = coupler.MF6Coupler(cpl.mm, ctx, cpl.mm.init_state(ctx), cpl.mf6b,
                              conv_fact=1000.0)
    cpl3.run(cm.FakeApi('toy', ctx.cMF.nlay, nrow, ncol, ctx.ncell,
                        nuzf=ctx.ncell + 3, heads0=699.0))
    assert 'up to 1.5 m ABOVE the land surface' in capsys.readouterr().out


def test_cells_without_a_drain_are_not_blamed_on_the_drain():
    src = open(coupler.__file__, encoding='utf-8').read()
    assert 'Only cells WITHOUT a seepage drain stand' in src
    assert '[:, drained]' in src
