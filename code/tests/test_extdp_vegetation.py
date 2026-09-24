# -*- coding: utf-8 -*-
"""WP2 2.2: UZF's extinction depth per vegetation zone.

La Mata's cover is a patchwork -- grass, holm oak and Pyrenean oak in
percent of each cell -- and every type carries its rooting depth. UZF has
one extinction depth per column, measured from the model's top, which is
the BASE of the MMsoil column: a type reaches UZF with the part of its root
zone below the soil, and the cell takes the cover-weighted mean.
"""

import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'app')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


props = _load('marmites_props_extdp', os.path.join(CODE, 'ppMF6',
                                                   'marmites_props.py'))
import marmites_config as mcfg                                # noqa: E402

ZR = [0.4, 15.0, 10.0]          # grass, holm oak, Pyrenean oak (lamata.toml)


def _cover(*cells):
    """(nveg, 1, n) from per-cell [grass, ilex, pyr] percentages."""
    return np.asarray(cells, dtype=float).T[:, None, :]


# ------------------------------------------------------------ the blend
def test_the_depth_is_the_cover_weighted_root_depth_below_the_soil():
    ext, _ = props.extdp_by_vegetation(_cover([60, 30, 0]), ZR, bare=2.0,
                                       soil_thick=0.5)
    # 60 % grass reaches 0 below 0.5 m of soil; 30 % oak 14.5 m; 10 % bare 2 m
    assert ext[0, 0] == pytest.approx(0.6 * 0.0 + 0.3 * 14.5 + 0.1 * 2.0)


def test_grass_shallower_than_the_soil_gives_uzf_nothing():
    """The soil column already takes what the grass draws: UZF gets 0."""
    ext, _ = props.extdp_by_vegetation(_cover([100, 0, 0]), ZR, bare=2.0,
                                       soil_thick=0.6)
    assert ext[0, 0] == 0.0


def test_a_thin_soil_lets_the_grass_reach_uzf():
    ext, _ = props.extdp_by_vegetation(_cover([100, 0, 0]), ZR, bare=2.0,
                                       soil_thick=0.1)
    assert ext[0, 0] == pytest.approx(0.3)


def test_bare_ground_keeps_the_depth_as_given():
    ext, _ = props.extdp_by_vegetation(_cover([0, 0, 0]), ZR, bare=2.0,
                                       soil_thick=0.5)
    assert ext[0, 0] == pytest.approx(2.0)


def test_a_sprinkle_of_oak_does_not_dry_the_whole_cell_to_15_m():
    """The reason for the mean, not the maximum."""
    ext, _ = props.extdp_by_vegetation(_cover([99, 1, 0]), ZR, bare=2.0)
    assert ext[0, 0] < 1.0


def test_the_bare_depth_and_the_soil_can_vary_per_cell():
    ext, _ = props.extdp_by_vegetation(
        _cover([0, 100, 0], [0, 100, 0], [0, 0, 0]), ZR,
        bare=np.array([[1.0, 1.0, 3.0]]),
        soil_thick=np.array([[0.5, 1.5, 0.5]]))
    assert np.allclose(ext[0], [14.5, 13.5, 3.0])


def test_nodata_soil_counts_as_no_soil():
    ext, _ = props.extdp_by_vegetation(_cover([0, 100, 0]), ZR, bare=2.0,
                                       soil_thick=np.array([[np.nan]]))
    assert ext[0, 0] == pytest.approx(15.0)


# ------------------------------------------------------- irrigated fields
def test_an_irrigated_cell_takes_its_crops_by_the_days_in_the_ground():
    """Field 1 grows crop 1 (1.2 m) for 30 days, lies fallow for 10 (bare,
    2 m) and grows crop 2 (0.6 m) for 60. The natural cover underneath is
    what MMsoil ignores there, so it must not count."""
    sched = np.array([[1, 0, 2]])
    ext, notes = props.extdp_by_vegetation(
        _cover([0, 100, 0], [0, 100, 0]), ZR, bare=2.0, soil_thick=0.2,
        irr=np.array([[1, 0]]), crop_by_period=sched,
        crop_root=[1.2, 0.6], perlen=[30.0, 10.0, 60.0])
    want = (30 * 1.0 + 10 * 2.0 + 60 * 0.4) / 100.0
    assert ext[0, 0] == pytest.approx(want)
    assert ext[0, 1] == pytest.approx(14.8), 'a dry cell took a crop'
    assert '1 irrigated cell(s)' in ' '.join(notes)


def test_only_the_periods_the_run_covers_count():
    sched = np.array([[1, 1, 2, 2]])
    ext, _ = props.extdp_by_vegetation(
        _cover([0, 0, 0]), ZR, bare=2.0, irr=np.array([[1]]),
        crop_by_period=sched, crop_root=[1.2, 0.6], perlen=[1.0, 1.0])
    assert ext[0, 0] == pytest.approx(1.2)


# ------------------------------------------------------------- refusals
def test_one_root_depth_per_vegetation_type():
    with pytest.raises(props.PropertyError):
        props.extdp_by_vegetation(_cover([50, 50, 0]), [0.4, 15.0], 2.0)


def test_a_root_depth_must_be_positive():
    with pytest.raises(props.PropertyError):
        props.extdp_by_vegetation(_cover([50, 50, 0]), [0.4, 0.0, 10.0], 2.0)


def test_a_crop_the_schedule_names_must_exist():
    with pytest.raises(props.PropertyError, match='crop 3'):
        props.extdp_by_vegetation(
            _cover([0, 0, 0]), ZR, 2.0, irr=np.array([[1]]),
            crop_by_period=np.array([[3]]), crop_root=[1.2, 0.6],
            perlen=[1.0])


def test_the_summary_is_area_weighted():
    _, notes = props.extdp_by_vegetation(
        _cover([0, 100, 0], [0, 0, 0]), ZR, bare=1.0,
        areas=np.array([[1.0, 3.0]]))
    assert 'catchment mean 4.5 m (area-weighted)' in notes[1]


# ------------------------------------------------------ the configuration
def test_the_default_is_the_depth_as_given():
    """A configuration that never asked keeps behaving as it did."""
    assert mcfg.RunConfig.from_dict({}).et.extdp_from == 'source'


def test_lamata_takes_it_per_vegetation_zone():
    cfg = mcfg.load_run_config(os.path.join(CODE, 'configs', 'lamata.toml'))
    assert cfg.et.extdp_from == 'vegetation'
    assert [v.root_depth for v in cfg.surface.vegetation] == ZR


def test_an_unknown_mode_is_refused():
    with pytest.raises(mcfg.ConfigError, match='et.extdp_from'):
        mcfg.RunConfig.from_dict({'et': {'extdp_from': 'species'}})


def test_vegetation_mode_needs_every_root_depth():
    cfg = mcfg.load_run_config(os.path.join(CODE, 'configs', 'lamata.toml'))
    cfg.surface.vegetation[0].root_depth = 0.0
    msgs = [e for e in cfg.problems() if 'et.extdp_from' in e]
    assert msgs and cfg.surface.vegetation[0].name in msgs[0]


def test_the_panel_asks_it():
    from lib import schema
    flat = [r for row in schema.UZF_ET_ROWS for r in row]
    assert 'et.extdp_from' in flat
    assert schema.CHOICES['et.extdp_from']() == ['source', 'vegetation']
    assert 'BELOW the soil' in schema.describe('et.extdp_from')[2]


# ------------------------------------------------------------- the driver
def test_the_driver_subtracts_the_soil_it_built_the_top_from():
    """setup_lamata puts the model top at elevation - soil thickness; the
    extinction depth must be measured from the same top."""
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'),
               encoding='utf-8').read()
    assert 'cMF.top = cMF.elev - np.ma.masked_values(gridSOILthick' in src
    body = src[src.index('def extdp_grid'):src.index('def active_area')]
    assert 'soil_thick=soil' in body and 'ctx.gridSOILthick' in body
    assert 'b.uzf_extdp = extdp_grid(cfg, cMF, ctx, DS)' in src
