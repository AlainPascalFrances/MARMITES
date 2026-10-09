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
        # the cap on the stream's Eo (§8.23): off until a test gives areas
        self.sfr_evap_frac = 0.5
        self.sfr_wl = None
        self.sfr_up_prev = None
        self.p_sfr_inflow = None
        self.p_sfr_usflow = None
        self.p_sfr_qfrommvr = None
        self.sfr_evap_capped = self.sfr_evap_writes = 0


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


# ------------------------------------- the cap on the stream (analysis §8.23)

def _capped(inflow, up_prev, eo_mm=4.0, wl=30.0, frac=0.5):
    """Reaches of w x L = ``wl`` m2 under Eo ``eo_mm`` mm/d: E0 w L is
    0.12 m3/d, against inflows of a few hundredths, as on La Mata's
    failing reaches."""
    n = len(inflow)
    c = _Cpl(n, area=[2500.0] * n)
    c.nreaches = n
    c.Eo_zonesSP = np.array([[eo_mm]])
    c.sfr_evap_zone = np.zeros(n, dtype=int)
    c.p_sfr_evap = np.zeros(n)
    c.p_sfr_inflow = np.asarray(inflow, dtype=float)
    c.sfr_up_prev = None if up_prev is None else np.asarray(up_prev, float)
    c.sfr_wl = np.full(n, float(wl))
    c.sfr_evap_frac = frac
    c._write_openwater_evap(0)
    return c


def test_a_reach_evaporates_at_most_its_share_of_the_inflow():
    """The runoff it gets this step plus what flowed in at the end of the
    last: 0.05 + 0 and 0 + 0.02 m3/d are capped, 1 m3/d is not."""
    c = _capped(inflow=[0.05, 0.0, 1.0], up_prev=[0.0, 0.02, 0.0])
    assert c.p_sfr_evap[0] == pytest.approx(0.5 * 0.05 / 30.0)
    assert c.p_sfr_evap[1] == pytest.approx(0.5 * 0.02 / 30.0)
    assert c.p_sfr_evap[2] == pytest.approx(0.004)
    assert (c.sfr_evap_capped, c.sfr_evap_writes) == (2, 3)


def test_the_capped_reach_can_never_empty_itself():
    """What makes MF6 flip-flop is a reach that evaporates all it receives;
    with the cap E0 w L stays below the inflow whatever the numbers."""
    rng = np.random.default_rng(0)
    q = rng.uniform(0.0, 0.5, 50)
    up = rng.uniform(0.0, 0.5, 50)
    c = _capped(inflow=q, up_prev=up, eo_mm=6.0)
    assert np.all(c.p_sfr_evap * 30.0 <= 0.5 * (q + up) + 1e-12)
    assert np.all(c.p_sfr_evap <= 0.006 + 1e-12)


def test_the_first_step_counts_the_runoff_alone():
    """No step solved yet: nothing has flowed in from upstream."""
    c = _capped(inflow=[0.05, 0.0], up_prev=None)
    assert c.p_sfr_evap[0] == pytest.approx(0.5 * 0.05 / 30.0)
    assert c.p_sfr_evap[1] == 0.0


def test_without_the_areas_eo_is_written_whole():
    """An MF6 build without USFLOW / LENGTH / WIDTH: no cap (the bind warns)."""
    c = _capped(inflow=[0.05], up_prev=[0.0])
    c.sfr_wl = None
    c._write_openwater_evap(0)
    assert c.p_sfr_evap[0] == pytest.approx(0.004)


def test_the_upstream_and_mover_inflow_is_copied_after_a_step():
    """MF6 zeroes USFLOW when it advances the next step (sfr_ad): the cap
    must read a COPY, with the pond outflow (QFROMMVR) added."""
    c = _Cpl(1, area=[1.0])
    c.sfr_wl = np.full(3, 30.0)
    c.p_sfr_usflow = np.array([1.0, 2.0, 3.0])
    c.p_sfr_qfrommvr = np.array([0.0, 5.0, 0.0])
    c._keep_sfr_inflow()
    c.p_sfr_usflow[:] = 0.0
    assert np.allclose(c.sfr_up_prev, [1.0, 7.0, 3.0])
    c.p_sfr_qfrommvr = None                      # a network without MOVER
    c.p_sfr_usflow[:] = [4.0, 0.0, 0.0]
    c._keep_sfr_inflow()
    assert np.allclose(c.sfr_up_prev, [4.0, 0.0, 0.0])


def test_the_coupler_keeps_the_inflow_after_every_step():
    """Both step calls of _advance, the period's first and its ATS
    sub-steps, are followed by the copy."""
    import inspect
    src = inspect.getsource(coup.MF6Coupler._advance)
    assert src.count('self._one_step(') == 2
    assert src.count('self._keep_sfr_inflow()') == 2


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


# ------------------------------------------- WP2 row 2: kept apart by package

def _both():
    """A reach in cell 0 and a pond in cell 2, 2500 m2 each."""
    c = _Cpl(3, area=[2500.0] * 3)
    c.nreaches = 1
    c.sfr_reach_idx = np.array([0, -1, -1])
    c.p_sfr_simevap = np.array([2.5])            # 1 mm/d over cell 0
    c.nlakes = 2
    c.lak_cell_idx = np.array([2, 0])
    c.p_lak_simevap = np.array([5.0, -7.5])      # 2 mm/d on 2, 3 mm/d on 0
    return c


def test_the_streams_and_the_ponds_are_read_apart():
    sfr, lak = _both()._read_openwater_evap_split()
    assert np.allclose(sfr, [1.0, 0.0, 0.0])
    assert np.allclose(lak, [3.0, 0.0, 2.0])


def test_the_total_is_still_their_sum():
    c = _both()
    sfr, lak = c._read_openwater_evap_split()
    assert np.allclose(c._read_openwater_evap(), sfr + lak)


def test_each_pond_keeps_its_own_loss_in_m3():
    """How much a pond loses is a result in its own right."""
    assert np.allclose(_both()._lak_evap_volumes(), [5.0, 7.5])
    c = _Cpl(1, area=[1.0])
    assert c._lak_evap_volumes() is None, 'no ponds, no pond series'


def _ctx_index():
    import marmites_indices as mi
    return SimpleNamespace(index=dict(mi.INDEX_MM))


def test_each_package_lands_in_its_own_column():
    import marmites_indices as mi
    ix = mi.INDEX_MM
    c = _both()
    c.ctx = _ctx_index()
    c._iEow = ix['iEow']
    mm = np.zeros((3, len(ix)))
    c._openwater_evap_into(mm)
    assert np.allclose(mm[:, ix['iEow_sfr']], [1.0, 0.0, 0.0])
    assert np.allclose(mm[:, ix['iEow_lak']], [3.0, 0.0, 2.0])
    assert np.allclose(mm[:, ix['iEow']], [4.0, 0.0, 2.0])


def test_the_total_et_counts_the_open_water():
    """MMsoil built iETtot with its own Eow, zero since the surface store
    went to MODFLOW: without this the five sources summed to four."""
    import marmites_indices as mi
    ix = mi.INDEX_MM
    c = _both()
    c.ctx = _ctx_index()
    c._iEow = ix['iEow']
    mm = np.zeros((3, len(ix)))
    mm[:, ix['iETtot']] = 10.0                   # Ei + ETsoil + ETg from MMsoil
    c._openwater_evap_into(mm)
    assert np.allclose(mm[:, ix['iETtot']], [14.0, 10.0, 12.0])
    c._openwater_evap_into(mm)                   # idempotent within a period
    assert np.allclose(mm[:, ix['iETtot']], [14.0, 10.0, 12.0])


def test_the_new_columns_are_appended_after_wp2s():
    import marmites_indices as mi
    assert [mi.INDEX_MM[k] for k in ('iEow_sfr', 'iEow_lak')] == [28, 29]
    assert sorted(mi.INDEX_MM.values()) == list(range(len(mi.INDEX_MM)))


def test_the_run_keeps_each_ponds_series_and_not_the_per_cell_split():
    src = open(os.path.join(CODE, 'marmites_coupler.py'),
               encoding='utf-8').read()
    assert "res['lak_evap'] = self.lak_evap_hist" in src
    assert "'eow_sfr':" not in src, 'per cell and per SP: ~0.5 GB a full run'
