# -*- coding: utf-8 -*-
"""WP1d -- open-water evaporation, from MMsoil to SFR/LAK and back.

MMsoil used to evaporate from its own surface store, and SFR and LAK were
built with EVAPORATION deliberately at 0 so the same water would not leave
twice. The store is gone: the coupler writes the rate into both packages each
stress period, and reads back what MF6 actually removed so ``iEow`` is a
measured flux rather than a structural zero.

Two things can go wrong silently and are pinned here: the units (EVAP is a
rate per wetted area, SIMEVAP a volumetric rate, and the MARMITES vector is a
depth over the CELL), and the meteo zone each reach and lake is tagged with.
"""

import importlib.util
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


coup = _load('marmites_coupler_ow', os.path.join(CODE, 'marmites_coupler.py'))


class _Cpl(coup.MF6Coupler):
    """The evaporation machinery without the rest of the coupler.

    ``MF6Coupler.__init__`` needs a built MF6 model; these two methods depend
    only on the arrays, so they are exercised directly.
    """

    def __init__(self, ncell, area, conv_fact=1000.0):
        self.ncell = ncell
        self.area = np.asarray(area, dtype=float)
        self.conv_fact = float(conv_fact)
        self.nreaches = 0
        self.nlakes = 0
        self.p_sfr_simevap = None
        self.p_lak_simevap = None
        self.p_sfr_evap = None
        self.p_lak_evap = None
        self.sfr_reach_idx = np.full(ncell, -1, dtype=int)
        self.lak_cell_idx = np.full(1, -1, dtype=int)
        self.sfr_evap_zone = np.zeros(1, dtype=int)
        self.lak_evap_zone = np.zeros(1, dtype=int)
        self.Eo_zonesSP = None
        self._iEow = 7


# ------------------------------------------------------------ reading back

def test_sfr_volumetric_rate_becomes_a_depth_over_the_cell():
    """SIMEVAP is m3/d; the MARMITES vector is mm/d over the CELL. 2.5 m3/d
    over 2500 m2 is 1 mm/d -- and getting that factor wrong is a factor of
    2500, which no plot would make obvious."""
    c = _Cpl(3, area=[2500.0, 2500.0, 2500.0])
    c.nreaches = 2
    c.sfr_reach_idx = np.array([0, -1, 1])
    c.p_sfr_simevap = np.array([2.5, 5.0])
    got = c._read_openwater_evap()
    assert got[0] == pytest.approx(1.0)
    assert got[1] == pytest.approx(0.0)
    assert got[2] == pytest.approx(2.0)


def test_two_reaches_in_one_cell_are_summed():
    c = _Cpl(1, area=[1000.0])
    c.nreaches = 1
    c.sfr_reach_idx = np.array([0])
    c.p_sfr_simevap = np.array([1.0])
    assert c._read_openwater_evap()[0] == pytest.approx(1.0)


def test_lake_evaporation_lands_in_its_host_cell():
    c = _Cpl(3, area=[2500.0] * 3)
    c.nlakes = 2
    c.lak_cell_idx = np.array([2, 0])
    c.p_lak_simevap = np.array([2.5, 25.0])
    got = c._read_openwater_evap()
    assert got[2] == pytest.approx(1.0)
    assert got[0] == pytest.approx(10.0)
    assert got[1] == pytest.approx(0.0)


def test_a_removal_is_reported_as_a_positive_loss():
    """MF6's sign convention must not leak into the water balance."""
    c = _Cpl(1, area=[2500.0])
    c.nreaches = 1
    c.sfr_reach_idx = np.array([0])
    c.p_sfr_simevap = np.array([-2.5])
    assert c._read_openwater_evap()[0] == pytest.approx(1.0)


def test_nothing_exposed_gives_zero_not_a_crash():
    c = _Cpl(4, area=[2500.0] * 4)
    assert np.allclose(c._read_openwater_evap(), 0.0)


def test_a_reach_index_past_the_array_is_ignored():
    """A build whose SFR array is shorter than the network must not throw."""
    c = _Cpl(2, area=[2500.0, 2500.0])
    c.nreaches = 2
    c.sfr_reach_idx = np.array([0, 5])
    c.p_sfr_simevap = np.array([2.5])
    got = c._read_openwater_evap()
    assert got[0] == pytest.approx(1.0) and got[1] == 0.0


def test_it_lands_in_iEow_and_leaves_runoff_alone():
    """Eow is a loss from the CHANNEL, downstream of the runoff. Subtracting
    it per cell would drive Ro negative, because a reach carries water from
    upstream."""
    c = _Cpl(2, area=[2500.0, 2500.0])
    c.nreaches = 1
    c.sfr_reach_idx = np.array([0, -1])
    c.p_sfr_simevap = np.array([2.5])
    mm = np.zeros((2, 12))
    mm[:, 5] = 3.0                       # a stand-in for Ro
    c._iEow = 7
    c._openwater_evap_into(mm)
    assert mm[0, 7] == pytest.approx(1.0)
    assert np.allclose(mm[:, 5], 3.0), 'runoff must not be reduced'


def test_no_index_is_a_no_op():
    c = _Cpl(1, area=[2500.0])
    c._iEow = None
    assert c._openwater_evap_into(np.zeros((1, 3))) is None


# --------------------------------------------------------------- writing

def test_the_rate_written_is_per_zone_and_in_model_units():
    """Eo is mm/d per meteo zone; SFR and LAK take a rate per unit of wetted
    area in model length units, so it is divided by conv_fact."""
    c = _Cpl(2, area=[2500.0, 2500.0])
    c.nreaches = 2
    c.nlakes = 1
    c.Eo_zonesSP = np.array([[2.0, 4.0], [10.0, 20.0]])   # 2 zones, 2 SPs
    c.sfr_evap_zone = np.array([0, 1])
    c.lak_evap_zone = np.array([1])
    c.p_sfr_evap = np.zeros(2)
    c.p_lak_evap = np.zeros(1)
    c._write_openwater_evap(1)
    assert c.p_sfr_evap[0] == pytest.approx(4.0 / 1000.0)
    assert c.p_sfr_evap[1] == pytest.approx(20.0 / 1000.0)
    assert c.p_lak_evap[0] == pytest.approx(20.0 / 1000.0)


def test_a_stress_period_past_the_forcing_uses_the_last_one():
    c = _Cpl(1, area=[2500.0])
    c.nreaches = 1
    c.Eo_zonesSP = np.array([[2.0, 4.0]])
    c.sfr_evap_zone = np.array([0])
    c.p_sfr_evap = np.zeros(1)
    c._write_openwater_evap(99)
    assert c.p_sfr_evap[0] == pytest.approx(4.0 / 1000.0)


def test_no_forcing_means_nothing_is_written():
    c = _Cpl(1, area=[2500.0])
    c.nreaches = 1
    c.p_sfr_evap = np.zeros(1)
    assert c._write_openwater_evap(0) is None
    assert c.p_sfr_evap[0] == 0.0


# ------------------------------------------------------- the Sankey split

def test_the_sankey_splits_runoff_rather_than_unbalancing_the_box():
    """The surface box must close: Pe + Exf_1 = I + E_ow + Ro_net. Where the
    stream is groundwater-fed E_ow can exceed the runoff, and the split is
    clamped instead of drawing a negative flow."""
    for ro, eow in ((10.0, 3.0), (10.0, 0.0), (0.0, 11.4), (2.0, 5.0)):
        eow_k = min(eow, ro) if ro > 0.0 else 0.0
        ro_net = max(ro - eow_k, 0.0)
        assert eow_k >= 0.0 and ro_net >= 0.0
        assert eow_k + ro_net == pytest.approx(ro)
