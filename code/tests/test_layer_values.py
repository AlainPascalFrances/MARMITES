# -*- coding: utf-8 -*-
"""One value PER LAYER for the aquifer layer properties (user,
2026-10-19): the anisotropy Kh/Kv about 2 in the saprolite and 10 in the
fissured layer (docs/MM-MF6/WP7_granite_anisotropy_review.md) -- given on
panel 4, built per layer, and calibrated per layer by PEST++ (bounds 1-3
and 2-30)."""

import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in (TRUNK, os.path.join(TRUNK, 'tools'), os.path.join(TRUNK, 'ppMF6')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg                                  # noqa: E402


def _cfg(nlay=2, **layers):
    cfg = mcfg.RunConfig.from_dict({'layers': dict({'nlay': nlay},
                                                   **layers)})
    return cfg


def test_a_layer_value_may_be_one_number_per_layer(tmp_path):
    cfg = _cfg(k33={'value': [2.0, 10.0]})
    assert cfg.layers.k33.values(2) == [2.0, 10.0]
    assert _cfg(k33={'value': 2.0}).layers.k33.values(2) == [2.0, 2.0]
    f = str(tmp_path / 'c.toml')
    cfg.write_toml(f)
    assert mcfg.load_run_config(f).layers.k33.value == [2.0, 10.0]
    # the panels' path: an override with a bracketed list
    cfg = _cfg(k33={'value': 2.0})
    cfg.apply_overrides(['layers.k33.value=[2.0, 10.0]'], echo=False)
    assert cfg.layers.k33.value == [2.0, 10.0]
    cfg.apply_overrides(['layers.k33.value=3.0'], echo=False)
    assert cfg.layers.k33.value == 3.0


@pytest.mark.parametrize('layers, frag', [
    ({'k33': {'value': [2.0, 10.0, 10.0]}}, 'gives 3 value(s), one per '
                                             'layer, for 2 layer(s)'),
    ({'k': {'value': [0.05, -1.0]}}, 'layers.k must not be negative'),
    ({'k33': {'value': [2.0, 0.0]}}, 'the ratio k/k33, so it must be > 0'),
])
def test_per_layer_values_are_checked(layers, frag):
    with pytest.raises(mcfg.ConfigError, match=frag.replace('(', r'\(')
                       .replace(')', r'\)').replace('/', '.')):
        _cfg(**layers)


def test_one_value_per_layer_is_for_the_layer_properties_only():
    cfg = _cfg()
    cfg.soil.thickness.value = [0.5, 1.0]
    assert any('one number, not one per layer' in e for e in cfg.problems())


def test_the_build_takes_one_value_per_layer(tmp_path):
    import marmites_props as props
    src = mcfg.VectorSource(value=[2.0, 10.0])
    assert props.resolve_source(src, 2, str(tmp_path), 'layers.k33') == \
        [2.0, 10.0]
    with pytest.raises(props.PropertyError, match='2 value'):
        props.resolve_source(src, 3, str(tmp_path), 'layers.k33')
    assert props.value_text(src) == '2 / 10'
    assert props.value_text(mcfg.VectorSource(value=0.05), ' m') == '0.05 m'


def test_the_anisotropy_is_calibrated_per_layer(tmp_path):
    """PEST: Kh/Kv starts at 2 (1-3) in layer 1 and 10 (2-30) in layer 2,
    whatever the configuration holds, and the forward run writes each
    layer's ratio into layers.k33."""
    import pest_lamata_setup as pls
    import pest_lamata_forward as pfw
    (tmp_path / 'MF_ws').mkdir()
    for L in (1, 2):
        (tmp_path / 'MF_ws' / ('hk_l%d.asc' % L)).write_text(
            'ncols 1\nnrows 1\nxllcorner 0\nyllcorner 0\ncellsize 50\n'
            'NODATA_value -9999\n0.05\n')
    cfg = _cfg(k={'raster': 'MF_ws/hk_l%d.asc'}, sy={'value': 0.01},
               k33={'value': 2.0})
    rows = pls.default_pars(cfg, str(tmp_path))
    k33 = {r['layer']: (r['start'], r['lower'], r['upper'], r['log'])
           for r in rows if r['family'] == 'k33'}
    assert k33 == {1: (2.0, 1.0, 3.0, True), 2: (10.0, 2.0, 30.0, True)}
    specs = {s[0]: s[1:] for s in pls.scalar_specs(cfg, rows)}
    assert specs['k33_l2'] == (10.0, 'log', 2.0, 30.0)
    cfg.pest.par = [mcfg.PestPar(**r) for r in rows]
    assert not [e for e in cfg.problems() if 'k33' in e]
    t = pd.DataFrame({'name': ['k33_l2'], 'value': [12.5]})
    pfw.apply_scalars(cfg, t, mult_dir=str(tmp_path))
    assert cfg.layers.k33.value == [2.0, 12.5]
    # a raster, or k33 as a conductivity: no ratio to calibrate
    cfg.layers.k33_as_ratio = False
    assert not [r for r in pls.default_pars(cfg, str(tmp_path))
                if r['family'] == 'k33']
    assert any('Kh/Kv' in e and 'k33_as_ratio' in e for e in cfg.problems())


def test_the_panel_offers_two_and_ten_for_an_unanswered_k33():
    sys.path.insert(0, os.path.join(TRUNK, 'app'))
    from lib import schema
    assert 'layers.k33' in schema.LAYER_VALUES
    assert tuple(schema.PER_LAYER_DEFAULT['layers.k33']) == (2.0, 10.0)
    assert 'ONE VALUE PER LAYER' in schema.describe('layers.k33')[2]
