# -*- coding: utf-8 -*-
"""WP1d -- the SFR network from the MAPPED LINES.

Until now SFR found its cells by looking for ``inputSTREAMw.asc > 0``, and
that raster was never a channel map: its values came from ``Soil_type.shp``,
where ``PONDw`` is 1.5 m on the two alluvium polygons and 0 elsewhere. So the
"stream network" was the alluvium footprint with one width for the whole
catchment.

These tests pin the replacement: the network is what the modeller mapped, a
cell is a stream cell when a line actually crosses it, and the width producer
is resolved AFTER routing -- which is the only point at which contributing
area is known.
"""

import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
REPO = os.path.abspath(os.path.join(CODE, '..'))
DS = os.path.join(REPO, 'example', 'LaMata')
for _p in (CODE, os.path.join(CODE, 'ppMF6')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


mch = _load('marmites_channel_t', os.path.join(CODE, 'ppMF6',
                                               'marmites_channel.py'))
mv = _load('marmites_vector_c', os.path.join(CODE, 'marmites_vector.py'))
cfgmod = _load('marmites_config_c', os.path.join(CODE, 'marmites_config.py'))


class _Net(object):
    """The little of SFRNetwork channel_width needs."""

    def __init__(self, cells, recv, acc):
        self.cells = cells
        self.recv = recv
        self.acc = acc


# ------------------------------------------------------------- reading

def test_read_stream_lines_orders_vertices_by_seq(tmp_path):
    p = tmp_path / 's.csv'
    p.write_text('# provenance\nseg_id,seq,x,y\n'
                 '0,1,10,0\n0,0,0,0\n0,2,20,0\n1,0,0,5\n1,1,10,5\n')
    lines, params = mch.read_stream_lines(str(p))
    assert sorted(lines) == [0, 1]
    assert lines[0] == [(0.0, 0.0), (10.0, 0.0), (20.0, 0.0)]
    assert params == {}


def test_read_stream_lines_picks_up_the_parameter_table(tmp_path):
    (tmp_path / 's.csv').write_text('seg_id,seq,x,y\n0,0,0,0\n0,1,10,0\n')
    (tmp_path / 'p.csv').write_text('# head\nseg_id,width_m,manning\n0,2.5,0.04\n')
    _lines, params = mch.read_stream_lines(str(tmp_path / 's.csv'),
                                           str(tmp_path / 'p.csv'))
    assert params[0]['width_m'] == '2.5'


def test_a_missing_table_says_how_to_make_it(tmp_path):
    with pytest.raises(mch.ChannelError) as e:
        mch.read_stream_lines(str(tmp_path / 'nope.csv'))
    assert 'gis_to_dataset' in str(e.value)


# -------------------------------------------------------------- burning

@pytest.fixture
def grid4():
    """2 x 2 cells of 10 m, origin (0, 0)."""
    return mv.TargetGrid.structured([10.0, 10.0], [10.0, 10.0], 0.0, 0.0)


def test_burn_channel_marks_only_the_cells_a_line_crosses(grid4):
    lines = {0: [(0.0, 15.0), (20.0, 15.0)]}       # across the northern row
    present, seg, length = mch.burn_channel(lines, grid4, (2, 2))
    assert present.tolist() == [[1.0, 1.0], [0.0, 0.0]]
    assert seg.tolist() == [[0, 0], [-1, -1]]
    assert length[0, 0] == pytest.approx(10.0)
    assert length.sum() == pytest.approx(20.0)


def test_burn_channel_gives_a_cell_to_its_dominant_segment(grid4):
    lines = {7: [(0.0, 15.0), (12.0, 15.0)],       # 10 m in cell 0, 2 m in 1
             9: [(12.0, 15.0), (20.0, 15.0)]}      # 8 m in cell 1
    _p, seg, _l = mch.burn_channel(lines, grid4, (2, 2))
    assert seg[0, 0] == 7
    assert seg[0, 1] == 9                          # 8 m beats 2 m


# ------------------------------------------------------- arbolate sum

def test_arbolate_sum_accumulates_downstream():
    """Two headwaters joining, then one trunk reach."""
    a, b, j, out = (0, 0), (0, 2), (1, 1), (2, 1)
    net = _Net(cells=[a, b, j, out],
               recv={a: j, b: j, j: out, out: None}, acc={})
    ln = np.zeros((3, 3))
    ln[a], ln[b], ln[j], ln[out] = 10.0, 20.0, 5.0, 7.0
    arb = mch.arbolate_sum(net, ln)
    assert arb[a] == pytest.approx(10.0)
    assert arb[b] == pytest.approx(20.0)
    assert arb[j] == pytest.approx(35.0)           # 10 + 20 + its own 5
    assert arb[out] == pytest.approx(42.0)         # + its own 7


# ------------------------------------------------------------- width

def _simple_net():
    cells = [(0, 0), (0, 1), (0, 2)]
    return _Net(cells, recv={cells[0]: cells[1], cells[1]: cells[2],
                             cells[2]: None},
                acc={cells[0]: 1, cells[1]: 2, cells[2]: 3})


def test_width_value_producer_is_uniform():
    src = cfgmod.ParamSource(value=2.0)
    w = mch.channel_width(_simple_net(), src, 2500.0, verbose=False)
    assert set(w.values()) == {2.0}


def test_width_column_producer_reads_the_segment_table():
    net = _simple_net()
    seg = np.array([[3, 3, 9]])
    params = {3: {'width_m': '1.2'}, 9: {'width_m': '4.0'}}
    src = cfgmod.ParamSource(column='width_m')
    w = mch.channel_width(net, src, 2500.0, seg_of_cell=seg, params=params,
                          verbose=False)
    assert w[(0, 0)] == pytest.approx(1.2)
    assert w[(0, 2)] == pytest.approx(4.0)


def test_width_column_producer_names_what_the_table_has():
    net = _simple_net()
    src = cfgmod.ParamSource(column='bankfull')
    with pytest.raises(mch.ChannelError) as e:
        mch.channel_width(net, src, 2500.0, seg_of_cell=np.array([[3, 3, 3]]),
                          params={3: {'width_m': '1.2'}}, verbose=False)
    assert 'width_m' in str(e.value)


def test_width_arbolate_producer_grows_downstream():
    """CdL's actual method: headwater value at the top, outlet value at the
    trunk, scaled by the normalised arbolate sum."""
    net = _simple_net()
    ln = np.array([[10.0, 10.0, 10.0]])
    src = cfgmod.ParamSource(drainage={'w_min': 1.5, 'w_max': 3.0,
                                       'power': 2.0})
    w = mch.channel_width(net, src, 2500.0, cell_length=ln, verbose=False)
    assert w[(0, 0)] < w[(0, 1)] < w[(0, 2)]
    assert w[(0, 2)] == pytest.approx(3.0)          # the trunk gets w_max
    assert w[(0, 0)] >= 1.5


def test_width_power_law_takes_area_in_km2():
    """w = 0.5 * A**0.35 with A in m2 gives 7.7 m for ONE 2500 m2 cell, where
    the mapped La Mata channel is 1.5 m. The unit is km2."""
    net = _simple_net()
    src = cfgmod.ParamSource(drainage={'a': 0.5, 'b': 0.35})
    w = mch.channel_width(net, src, 2500.0, verbose=False)
    assert max(w.values()) < 1.0


def test_width_with_no_producer_is_refused():
    with pytest.raises(mch.ChannelError) as e:
        mch.channel_width(_simple_net(), cfgmod.ParamSource(), 2500.0,
                          verbose=False)
    assert 'no producer' in str(e.value)


def test_width_is_floored():
    src = cfgmod.ParamSource(value=0.01)
    w = mch.channel_width(_simple_net(), src, 2500.0, minimum=0.5,
                          verbose=False)
    assert set(w.values()) == {0.5}


# -------------------------------------------------------------- depth

def test_depth_value_producer():
    d = mch.channel_depth(_simple_net(), cfgmod.ParamSource(value=1.25),
                          verbose=False)
    assert set(d.values()) == {1.25}


def test_depth_rejects_a_producer_it_cannot_resolve():
    with pytest.raises(mch.ChannelError) as e:
        mch.channel_depth(_simple_net(), cfgmod.ParamSource(column='d_m'),
                          verbose=False)
    assert 'single value' in str(e.value)


# ------------------------------------------------------- the real data

def test_the_committed_network_burns_onto_the_legacy_grid():
    """97 mapped segments, on the 65 x 60 @ 50 m grid."""
    stream = os.path.join(DS, 'inputSTREAM.csv')
    if not os.path.exists(stream):
        pytest.skip('inputSTREAM.csv not present')
    lines, params = mch.read_stream_lines(
        stream, os.path.join(DS, 'inputSTREAM_param.csv'))
    assert len(lines) == 97 and len(params) == 97
    grid = mv.TargetGrid.structured([50.0] * 60, [50.0] * 65,
                                    739300.0, 4553050.0)
    present, seg, length = mch.burn_channel(lines, grid, (65, 60))
    assert 300 < int(present.sum()) < 400
    assert 14000 < float(length.sum()) < 15000
    assert set(np.unique(seg)) - {-1} <= set(lines)


# ------------------------------------------------ the schema accepts both

@pytest.mark.parametrize('drainage', [
    {'w_min': 1.5, 'w_max': 3.0, 'power': 2.0},
    {'w_min': 1.0, 'w_max': 8.0},
    {'a': 0.5, 'b': 0.35},
])
def test_both_drainage_forms_validate(drainage):
    c = cfgmod.RunConfig.from_dict({'sfr': {'width': {'drainage': drainage}}})
    assert c.sfr.width.producer() == 'drainage'


@pytest.mark.parametrize('drainage, frag', [
    ({'w_min': 3.0, 'w_max': 1.0}, 'w_max must be'),
    ({'w_min': 1.0}, 'needs either'),
    ({'a': 0.5}, 'needs either'),
    ({'a': 0.5, 'b': 0.35, 'c': 1.0}, 'unknown key'),
])
def test_a_malformed_drainage_block_is_refused(drainage, frag):
    with pytest.raises(cfgmod.ConfigError) as e:
        cfgmod.RunConfig.from_dict({'sfr': {'width': {'drainage': drainage}}})
    assert frag in str(e.value)
