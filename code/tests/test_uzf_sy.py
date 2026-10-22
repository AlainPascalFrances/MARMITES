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

Since 2026-10-20 (user) both come from the aquifer layer each UZF object
sits in: thtr = its specific retention Sr (layers.sr), thts = Sy + Sr, the
porosity -- consistent by construction. One panel thts for every layer had
given thtr = 0.45 - 0.01 = 0.44 in the saprolite and the fissured layer
alike.
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
def test_thtr_is_sr_thts_is_sy_plus_sr_and_thti_a_place_between():
    thtr, thts, thti, notes = props.uzf_water_contents(0.275, 0.01)
    assert float(thtr) == pytest.approx(0.275)
    assert float(thts) == pytest.approx(0.285)
    assert float(thti) == pytest.approx(0.275), 'a cold start is empty'
    assert any('thtr = Sr' in n and 'Sy + Sr' in n for n in notes)
    assert any('thtr (empty)' in n for n in notes)
    for fill, want in ((1.0, 0.285), (0.4, 0.279)):
        assert float(props.uzf_water_contents(0.275, 0.01, fill)[2]) == \
            pytest.approx(want)


def test_per_layer_and_per_cell_sy_and_sr():
    sr = np.array([0.275, 0.04])[:, None, None]
    sy = np.array([[[0.01, 0.03]], [[0.005, 0.01]]])
    thtr, thts, _i, _ = props.uzf_water_contents(sr, sy, 0.2)
    assert np.allclose(thts - thtr, sy), 'the zone does not drain Sy'
    assert np.allclose(thtr[0], 0.275) and np.allclose(thtr[1], 0.04)


def test_a_zero_sr_or_a_porosity_above_one_is_refused():
    with pytest.raises(props.PropertyError) as exc:
        props.uzf_water_contents(np.array([0.1, 0.0]), 0.01, 0.15)
    assert 'layers.sr' in str(exc.value) and 'in 1 cell' in str(exc.value)
    with pytest.raises(props.PropertyError) as exc:
        props.uzf_water_contents(0.95, 0.1, 0.5)
    assert 'Sy + Sr' in str(exc.value) and 'above 1' in str(exc.value)


def test_inactive_cells_are_not_judged():
    thtr, thts, thti, _ = props.uzf_water_contents(
        np.array([0.275, 0.0]), 0.01, active=np.array([True, False]))
    assert thtr[0] == pytest.approx(0.275) and thtr[1] == 0.0


def test_one_place_fits_every_layer_where_one_value_could_not():
    """La Mata: [0.275, 0.285] and [0.04, 0.05] -- no absolute value lies
    in both; the place 0.5 is in each (user, 2026-10-22)."""
    sr = np.array([0.275, 0.04])
    thtr, thts, thti, _ = props.uzf_water_contents(sr, 0.01, 0.5)
    assert np.all((thtr < thti) & (thti < thts))
    assert np.allclose(thti, [0.28, 0.045])


# ------------------------------------------------------ the configuration
def test_sr_is_a_layer_property_27_5_and_4_percent():
    c = mcfg.RunConfig.from_dict({})
    assert c.layers.sr.values(2) == [0.275, 0.04]
    assert mcfg.SR_DEFAULT == (0.275, 0.04)
    for gone in ('thts', 'thtr', 'thtr_from'):
        assert not hasattr(c.uzf, gone), gone


def test_an_old_file_with_the_uzf_contents_still_loads():
    c = mcfg.RunConfig.from_dict({'uzf': {'thtr_from': 'sy',
                                          'thtr': {'value': 0.05},
                                          'thts': {'value': 0.45},
                                          'thti': {'value': 0.15},
                                          'vks_from': 'layer'}})
    assert c.layers.sr.values(2) == [0.275, 0.04]
    assert {'thts', 'thtr', 'thtr_from', 'thti', 'vks_from', 'vks'} <= \
        set(mcfg.RETIRED['uzf'])


def test_a_bad_sr_is_a_configuration_error():
    for bad in ([0.0, 0.04], [0.3, 1.2]):
        with pytest.raises(mcfg.ConfigError) as exc:
            mcfg.RunConfig.from_dict({'layers': {'sr': {'value': bad}}})
        assert 'layers.sr' in str(exc.value)


# ------------------------------------------------------------- the panel
def test_sr_is_asked_on_the_layers_tab_one_box_per_layer():
    from lib import schema
    rows = [k for row in schema.GEOMETRY_ROWS for k in row if k]
    assert rows.index('layers.sr') == rows.index('layers.sy') + 1
    assert 'layers.sr' in schema.LAYER_VALUES
    assert schema.PER_LAYER_DEFAULT['layers.sr'] == (0.275, 0.04)
    assert schema.describe('layers.sr')[0] == 'Specific retention'
    uzf_rows = [k for row in schema.UZF_ROWS for k in row if k]
    for gone in ('uzf.thtr_from', 'uzf.thtr', 'uzf.thts'):
        assert gone not in uzf_rows
    page = open(os.path.join(CODE, 'app', 'pages',
                             '4_4_-_Unsaturated_zone_and_groundwater.py'),
                encoding='utf-8').read()
    assert 'thtr_from' not in page


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


def _cfg(tmp_path, sr=None, sy=0.01):
    ds = tmp_path / 'ds'
    ds.mkdir()
    for L in (1, 2):
        _asc(str(ds / ('Sy_l%d.asc' % L)), sy)
        _asc(str(ds / ('ib_l%d.asc' % L)), 1)
    cfg = mcfg.RunConfig.from_dict({'layers': {'nlay': 2}})
    cfg.layers.sy = mcfg.VectorSource(raster='Sy_l%d.asc')
    cfg.layers.ibound = mcfg.VectorSource(raster='ib_l%d.asc')
    if sr is not None:
        cfg.layers.sr = mcfg.VectorSource(value=sr)
    return cfg, str(ds)


def test_the_validation_panel_stops_a_porosity_above_one(tmp_path):
    chk = _checks()
    cfg, ds = _cfg(tmp_path, sr=[0.275, 0.04], sy=0.8)
    got = list(chk.check_uzf_sy(cfg, dataset_dir=ds))
    assert [c.level for c in got] == [chk.ERROR]
    assert got[0].key == 'layers.sr' and got[0].panel == 4
    assert 'above 1' in got[0].detail


def test_the_validation_panel_says_what_uzf_will_take(tmp_path):
    chk = _checks()
    cfg, ds = _cfg(tmp_path)
    got = list(chk.check_uzf_sy(cfg, dataset_dir=ds))
    assert got and all(c.level == chk.INFO for c in got)
    text = ' '.join(c.title for c in got)
    assert '0.04..0.275' in text and '0.05..0.285' in text
    assert 'thtr (empty)' in text


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
