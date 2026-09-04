# -*- coding: utf-8 -*-
"""Decision 2/3: groundwater seepage as a smoothed land-surface drain.

UZF's SIMULATE_GWSEEP is deprecated in MF6 and switches discharge on and off
discontinuously. The alternative, used by the CdL reference model, is a DRN at
the land surface with AUXDEPTHNAME, so MF6 ramps the discharge in cubically
over DDRN as the head approaches the surface.

Unlike the reference model, MARMITES does NOT move these flows to SFR with
MVR: the water has to come back into the soil column (Eq. 1 / sec. 2.3), so
the coupler reads the drain SIMVALS and feeds them in as exfiltration.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..', 'trunk'))
DS = os.path.abspath(os.path.join(HERE, '..', 'DataSet_LaMata'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

flopy = pytest.importorskip('flopy')
import matplotlib  # noqa: E402
matplotlib.use('agg')

from test_mf6_build import cmf, mf6mod, _load  # noqa: E402,F401  (pytest fixture)

MC = _load('marmites_coupler', os.path.join(TRUNK, 'marmites_coupler.py'))


def _built(cmf, tmp, seep='drn', cond=None, sfr_cells=()):
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp), daily=True)
    b.seep = seep
    if cond is not None:
        b.drn_seep_cond = cond
    b.sfr_cells = set(sfr_cells)
    b.build()
    return b


def _pkg(b, name):
    try:
        return b.gwf.get_package(name)
    except Exception:
        return None


# ------------------------- package construction ----------------------- #

def test_uzf_mode_builds_no_seep_drain(cmf, tmp_path):
    b = _built(cmf, tmp_path, seep='uzf')
    assert _pkg(b, 'drn_seep') is None
    assert b.ndrnseep == 0 and b.drnseep_id == {}


def test_drn_mode_disables_uzf_gwseep(cmf, tmp_path):
    """Both mechanisms at once would count the discharge twice."""
    b = _built(cmf, tmp_path, seep='drn')
    uzf = _pkg(b, 'uzf')
    assert bool(uzf.simulate_gwseep.get_data()) is False


def test_uzf_mode_keeps_uzf_gwseep(cmf, tmp_path):
    b = _built(cmf, tmp_path, seep='uzf')
    assert bool(_pkg(b, 'uzf').simulate_gwseep.get_data()) is True


def test_seep_drain_on_every_surface_cell(cmf, tmp_path):
    b = _built(cmf, tmp_path, seep='drn')
    assert b.ndrnseep == b.ncell
    assert len(b.drnseep_id) == b.ncell
    assert sorted(b.drnseep_id.values()) == list(range(b.ncell))


def test_drain_sits_at_the_land_surface_with_the_given_conductance(cmf, tmp_path):
    b = _built(cmf, tmp_path, seep='drn', cond=10000.0)
    spd = _pkg(b, 'drn_seep').stress_period_data.get_data(0)
    for rec, (i, j, k) in zip(spd, b.surf_cells):
        assert rec['elev'] == pytest.approx(float(b.top[i, j]))
        assert rec['cond'] == pytest.approx(10000.0)
        assert rec['ddrn'] == pytest.approx(b.drn_seep_ddrn)


def test_default_conductance_is_free_draining(cmf, tmp_path):
    """A seepage face must discharge without the head having to climb above
    ground. La Mata peaks at ~2060 m3/d of seepage per cell, so a conductance
    of order 10 m2/d needs a 200 m head excess to pass it -- the first coupled
    run overshot by 205 m. Anything below ~2000 m2/d is not free-draining."""
    b = _built(cmf, tmp_path, seep='drn')
    assert b.drn_seep_cond >= 2000.0


def test_ddrn_auxiliary_enables_cubic_smoothing(cmf, tmp_path):
    """Without AUXDEPTHNAME the drain switches on discontinuously."""
    b = _built(cmf, tmp_path, seep='drn')
    drn = _pkg(b, 'drn_seep')
    # flopy returns the option as a ('auxiliary', 'ddrn') record
    assert 'ddrn' in str(drn.auxiliary.get_data()).lower()
    assert str(drn.auxdepthname.get_data()).lower() == 'ddrn'
    assert b.drn_seep_ddrn > 0


def test_sfr_cells_get_no_seep_drain(cmf, tmp_path):
    """Decision 4: where SFR exists it takes the seepage, not a drain."""
    ref = _built(cmf, tmp_path, seep='drn')
    skip = {(ref.surf_cells[0][0], ref.surf_cells[0][1]),
            (ref.surf_cells[1][0], ref.surf_cells[1][1])}
    b = _built(cmf, tmp_path, seep='drn', sfr_cells=skip)
    assert b.ndrnseep == ref.ndrnseep - len(skip)
    for c in skip:
        assert c not in b.drnseep_id
    # the remaining indices stay a dense 0..n-1 range
    assert sorted(b.drnseep_id.values()) == list(range(b.ndrnseep))


def test_legacy_outlet_drn_still_present_alongside(cmf, tmp_path):
    b = _built(cmf, tmp_path, seep='drn')
    assert _pkg(b, 'drn') is not None, 'the outlet DRN must not be overwritten'
    assert _pkg(b, 'drn_seep') is not None


# --------------------------- coupler read-back ------------------------ #

class _Ptr(np.ndarray):
    pass


def _coupler_with_drnseep(sim_vals, mapping):
    """Minimal MF6Coupler shell exercising only _read_heads_exf."""
    c = MC.MF6Coupler.__new__(MC.MF6Coupler)
    n = len(mapping)
    c.ncell = n
    c.conv_fact = 1000.0
    c.area = np.full(n, 100.0)
    c.x_index = np.arange(n)
    c.p_x = np.zeros(n)
    c.p_gwd = None
    c.p_drnseep = np.asarray(sim_vals, dtype=float)
    c.drnseep_idx = np.asarray(mapping, dtype=int)
    return c


def test_drnseep_sign_flip_and_units():
    """DRN SIMVALS are negative (water out of the aquifer); exfiltration is
    positive, in mm/d."""
    c = _coupler_with_drnseep([-200.0, -50.0], [0, 1])
    _, exf = c._read_heads_exf()
    # 200 m3/d over 100 m2 = 2 m/d = 2000 mm/d
    assert exf == pytest.approx([2000.0, 500.0])


def test_drnseep_positive_simvals_are_not_infiltration():
    """A drain can only remove water; a positive value must not be read as
    water entering the soil."""
    c = _coupler_with_drnseep([+300.0, -100.0], [0, 1])
    _, exf = c._read_heads_exf()
    assert exf[0] == 0.0 and exf[1] == pytest.approx(1000.0)


def test_cells_without_a_drain_get_zero():
    c = _coupler_with_drnseep([-100.0], [-1, 0, -1])
    _, exf = c._read_heads_exf()
    assert exf[0] == 0.0 and exf[2] == 0.0
    assert exf[1] == pytest.approx(1000.0)


def test_boundary_order_is_not_assumed_to_be_cell_order():
    """The DRN arrays are ordered by boundary. Reading them positionally would
    silently attribute each cell's seepage to a different cell."""
    c = _coupler_with_drnseep([-100.0, -900.0], [1, 0])   # reversed mapping
    _, exf = c._read_heads_exf()
    assert exf == pytest.approx([9000.0, 1000.0])


def test_drnseep_takes_precedence_over_gwd():
    c = _coupler_with_drnseep([-100.0, -100.0], [0, 1])
    c.p_gwd = np.array([1e6, 1e6])       # stale/absent UZF discharge
    _, exf = c._read_heads_exf()
    assert exf == pytest.approx([1000.0, 1000.0])
