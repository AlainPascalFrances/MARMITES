# -*- coding: utf-8 -*-
"""Decision 1: collapse numerical layers onto hydrogeological units.

The La Mata ini declares Mnlay=2 with Mlay=[1,1,1,2,2,2]: six numerical layers
representing two hydrogeological units. Aggregation must preserve the physics
that matters:

  * transmissivity  -> horizontal K is a thickness-weighted ARITHMETIC mean
  * vertical resistance -> vertical K is a thickness-weighted HARMONIC mean
  * the water-table storage -> sy comes from the unit's UPPERMOST layer
  * the model footprint -> a unit is active if ANY of its layers is
"""
import importlib.util
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

import matplotlib  # noqa: E402
matplotlib.use('agg')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L = _load('marmites_layers', os.path.join(TRUNK, 'ppMF6', 'marmites_layers.py'))


def _toy(nrow=3, ncol=4, Mlay=(1, 1, 2, 2)):
    """4 layers -> 2 units, with contrasting K so the means are testable."""
    nlay = len(Mlay)
    o = np.ones((nlay, nrow, ncol))
    return SimpleNamespace(
        nlay=nlay, nrow=nrow, ncol=ncol,
        Mlay=list(Mlay), Mnlay=len(set(Mlay)),
        ibound=np.ones((nlay, nrow, ncol), dtype=int),
        botm=np.stack([o[0] * b for b in (90.0, 80.0, 60.0, 30.0)]),
        thick=np.stack([o[0] * t for t in (10.0, 10.0, 20.0, 30.0)]),
        strt=np.stack([o[0] * s for s in (95.0, 94.0, 93.0, 92.0)]),
        hk_actual=np.stack([o[0] * k for k in (1.0, 3.0, 5.0, 15.0)]),
        vka_actual=np.stack([o[0] * k for k in (0.5, 1.5, 2.5, 7.5)]),
        ss_actual=np.stack([o[0] * s for s in (1e-5, 3e-5, 1e-4, 3e-4)]),
        sy_actual=np.stack([o[0] * s for s in (0.10, 0.20, 0.30, 0.40)]),
        layvka=[0] * nlay, laytyp=[1] * nlay, layavg=[0] * nlay,
        laywet=[0] * nlay, laycbd=[0] * nlay,
        h_plt=[1] * nlay, h_lbl=['L%d' % k for k in range(nlay)],
        iuzfbnd=np.ones((nrow, ncol), dtype=int),
        layer_row_column_elevation_cond={0: [[3, 1, 1, 55.0, 0.02]]},
    )


def test_layer_groups_from_mlay():
    assert L.layer_groups([1, 1, 1, 2, 2, 2], 6) == [[0, 1, 2], [3, 4, 5]]
    assert L.layer_groups([1, 2], 2) == [[0], [1]]


def test_layer_groups_rejects_non_contiguous():
    with pytest.raises(ValueError, match='contiguous'):
        L.layer_groups([1, 2, 1], 3)


def test_layer_groups_rejects_wrong_length():
    with pytest.raises(ValueError, match='Mlay'):
        L.layer_groups([1, 1], 3)


def test_geometry_preserved():
    c = _toy()
    L.aggregate_layers(c, verbose=False)
    assert c.nlay == 2
    # unit bottoms = bottom of each unit's lowest layer
    assert np.allclose(c.botm[0], 80.0) and np.allclose(c.botm[1], 30.0)
    # thickness is summed
    assert np.allclose(c.thick[0], 20.0) and np.allclose(c.thick[1], 50.0)


def test_horizontal_k_is_thickness_weighted_arithmetic():
    """Preserves transmissivity: (1*10 + 3*10)/20 = 2."""
    c = _toy()
    L.aggregate_layers(c, verbose=False)
    assert np.allclose(c.hk_actual[0], 2.0)
    assert np.allclose(c.hk_actual[1], (5.0 * 20 + 15.0 * 30) / 50.0)


def test_vertical_k_is_thickness_weighted_harmonic():
    """Preserves vertical resistance: 20 / (10/0.5 + 10/1.5)."""
    c = _toy()
    L.aggregate_layers(c, verbose=False)
    expect0 = 20.0 / (10.0 / 0.5 + 10.0 / 1.5)
    expect1 = 50.0 / (20.0 / 2.5 + 30.0 / 7.5)
    assert np.allclose(c.vka_actual[0], expect0)
    assert np.allclose(c.vka_actual[1], expect1)
    # harmonic <= arithmetic, always
    assert c.vka_actual[0].max() <= (0.5 * 10 + 1.5 * 10) / 20 + 1e-12


def test_layvka_reset_so_vka_is_no_longer_a_ratio():
    """La Mata has layvka=1 (VKA is the hk/vk RATIO). After aggregation the
    stored value is a true vertical K, so layvka must be 0 or downstream code
    would invert it a second time."""
    c = _toy()
    c.layvka = [1] * c.nlay
    c.vka_actual = np.stack([np.ones((c.nrow, c.ncol)) * r for r in (2.0, 2.0, 4.0, 4.0)])
    L.aggregate_layers(c, verbose=False)
    assert c.layvka == [0, 0]
    # unit 1: k33 per layer = hk/ratio = 1/2 and 3/2 -> harmonic over 10/10 m
    expect = 20.0 / (10.0 / 0.5 + 10.0 / 1.5)
    assert np.allclose(c.vka_actual[0], expect)


def test_sy_taken_from_uppermost_layer():
    """sy acts at the water table, so averaging it down the column is wrong."""
    c = _toy()
    L.aggregate_layers(c, verbose=False)
    assert np.allclose(c.sy_actual[0], 0.10)
    assert np.allclose(c.sy_actual[1], 0.30)


def test_ss_is_thickness_weighted():
    c = _toy()
    L.aggregate_layers(c, verbose=False)
    assert np.allclose(c.ss_actual[0], (1e-5 * 10 + 3e-5 * 10) / 20)


def test_unit_active_if_any_layer_active():
    c = _toy()
    c.ibound[0, 0, 0] = 0          # only the top layer of unit 1 is inactive
    c.ibound[1, 0, 0] = 1
    c.ibound[2, 1, 1] = 0          # both layers of unit 2 inactive at (1,1)
    c.ibound[3, 1, 1] = 0
    L.aggregate_layers(c, verbose=False)
    assert c.ibound[0, 0, 0] != 0, 'unit stays active if any sub-layer is'
    assert c.ibound[1, 1, 1] == 0, 'unit inactive only if all sub-layers are'


def test_drn_cell_list_remapped_to_units():
    """DRN records are keyed by layer index and must follow the collapse."""
    c = _toy()
    L.aggregate_layers(c, verbose=False)
    rec = c.layer_row_column_elevation_cond[0][0]
    assert rec[0] == 1, 'layer 3 belongs to unit 2 -> index 1'
    assert rec[1:] == [1, 1, 55.0, 0.02]


def test_no_op_when_already_aggregated():
    c = _toy(Mlay=(1, 2, 3, 4))
    rep = L.aggregate_layers(c, verbose=False)
    assert rep['aggregated'] is False and c.nlay == 4


# --------------------------------------------------------------------- #
# real La Mata
# --------------------------------------------------------------------- #

def test_lamata_six_to_two():
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s3L.ini')):
        pytest.skip('La Mata dataset not present')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s3L.ini',
                   xllcorner=739300.0, yllcorner=4553050.0)
    assert c.nlay == 6 and c.Mnlay == 2 and c.Mlay == [1, 1, 1, 2, 2, 2]
    before = int((np.abs(np.asarray(c.ibound)) > 0).sum())
    rep = L.aggregate_layers(c, verbose=False)
    assert rep['aggregated'] and c.nlay == 2
    after = int((np.abs(np.asarray(c.ibound)) > 0).sum())
    assert after < before, 'aggregation must reduce the active-cell count'
    # the surface footprint is unchanged: every column active before is active now
    outcrop = np.zeros((c.nrow, c.ncol), dtype=int)
    for k in range(c.nlay):
        outcrop += ((outcrop == 0) & (np.abs(np.asarray(c.ibound))[k] != 0)) * (k + 1)
    assert int((outcrop > 0).sum()) == 1954
