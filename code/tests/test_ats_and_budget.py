# -*- coding: utf-8 -*-
"""Adaptive time stepping, and proving MF6 actually solved the problem.

Background (Phase-5 La Mata run): the first transient stress period -- a step
change from the steady state onto a free-draining seepage boundary, in one
1-day step -- did not converge. MF6 moved on anyway. That single period
discharged 9.6e6 m3 through the seepage drain, eight times the recharge of the
entire 1949-day simulation, and left a cumulative mass-balance discrepancy of
-132%. The head field still looked ordinary and the iteration counts looked
modest, so nothing downstream noticed.

Two defences, both tested here:
  * ATS, so MF6 may subdivide a period it cannot solve in one step;
  * check_solution(), which reads the listing files and FAILS the run on a
    non-converged period or a budget that does not close.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

import matplotlib  # noqa: E402
matplotlib.use('agg')

from test_coupler_mock import FakeApi, _setup, coup  # noqa: E402


# --------------------------------------------------------------------- #
# ATS: the coupler must follow MF6 across sub-steps
# --------------------------------------------------------------------- #

class ClockApi(FakeApi):
    """FakeApi with a clock, splitting each period into ``nsub`` sub-steps."""

    def __init__(self, *a, nsub=1, **kw):
        super().__init__(*a, **kw)
        self.nsub = int(nsub)
        self.t = 0.0
        self.steps = 0

    def get_current_time(self):
        return self.t

    def finalize_time_step(self):
        super().finalize_time_step()
        self.steps += 1
        self.t += 1.0 / self.nsub          # each period is 1.0 long


def _unused_clock_for(api, nsub):  # kept out of the way
    return ClockApi(api.name, api.X.size // (api.BOUND.shape[0] or 1) or 1, 1, 1,
                    api.Q.size, api.FINF.size, nsub=nsub)


def _run(nsub, nper=3, mode='lagged'):
    cpl, api, ctx = _setup(nper=nper, mode=mode)
    clock = ClockApi('toy', 1, 1, api.X.size, api.Q.size, api.FINF.size, nsub=nsub)
    clock.X = api.X.copy()
    cpl.run(clock)
    return cpl, clock


def test_single_substep_behaves_as_before():
    cpl, api = _run(nsub=1, nper=3)
    assert api.steps == 4               # 1 steady + 3 transient
    assert cpl.substeps == 0


def test_coupler_follows_mf6_when_ats_splits_a_period():
    """The whole point: without this loop MF6 would fall behind MARMITES and
    the two models would silently be simulating different days."""
    cpl, api = _run(nsub=4, nper=3)
    assert api.steps == 4 * 4, 'each period must be advanced to its end'
    assert cpl.substeps == 3 * 4         # 3 extra sub-steps per period


def test_mf6_clock_reaches_the_end_of_every_period():
    cpl, api = _run(nsub=3, nper=5)
    assert api.get_current_time() == pytest.approx(6.0)   # 1 steady + 5 days


def test_iterative_mode_also_completes_the_period():
    cpl, api = _run(nsub=3, nper=3, mode='iterative')
    assert api.get_current_time() == pytest.approx(4.0)


def test_soil_model_still_steps_once_per_day_when_ats_splits():
    """Sub-steps are a groundwater-solver detail; re-running the soil model per
    sub-step would advance the soil state several times within one day. The
    fluxes are period rates, so they must be written once and left alone."""
    cpl, api = _run(nsub=5, nper=4)
    assert cpl.perc_hist.shape[0] == 4
    # FINF is snapshotted at every finalize_time_step; within a period every
    # sub-step must see the same value
    snaps = np.array(api.finf_at_advance)
    for k in range(1, 5):                       # the 4 transient periods
        blk = snaps[k * 5:(k + 1) * 5]
        assert np.allclose(blk, blk[0]), 'fluxes were rewritten mid-period'


def test_runaway_substepping_is_capped():
    cpl, api, ctx = _setup(nper=2)
    cpl.max_substeps = 10
    clock = ClockApi('toy', 1, 1, api.X.size, api.Q.size, api.FINF.size,
                     nsub=10_000)                   # never reaches the end
    clock.X = api.X.copy()
    with pytest.raises(coup.CouplingError, match='did not reach the end'):
        cpl.run(clock)


def test_missing_clock_degrades_to_one_step_per_period():
    """An MF6 build without get_current_time must still run (no ATS then)."""
    cpl, api, ctx = _setup(nper=3)
    cpl.run(api)                                    # FakeApi has no clock
    assert cpl.substeps == 0


def test_max_outer_defaults_to_the_solver_limit():
    """If the manual outer loop caps below MF6's OUTER_MAXIMUM, MF6 never
    registers a step as failed and ATS cannot retry it -- that was the cause of
    the first transient day advancing unconverged. The default must match the
    solver's own limit, not a lower hard-coded 200."""
    cpl, api, ctx = _setup(nper=2)
    cpl.mf6b.outer_maximum = 500
    rebuilt = coup.MF6Coupler(cpl.mm, cpl.ctx, cpl.state, cpl.mf6b,
                              conv_fact=1000.0, mode='lagged')
    assert rebuilt.max_outer == 500
    # an explicit override still wins
    ov = coup.MF6Coupler(cpl.mm, cpl.ctx, cpl.state, cpl.mf6b,
                         conv_fact=1000.0, mode='lagged', max_outer=42)
    assert ov.max_outer == 42


def test_construction_accepts_a_numpy_perlen():
    """La Mata's perlen is an ndarray; `mf6b.perlen or cMF.perlen` used to raise
    'truth value of an array is ambiguous' before it was made array-safe."""
    cpl, api, ctx = _setup(nper=3)
    cpl.mf6b.perlen = np.array([1.0, 2.0, 3.0])
    rebuilt = coup.MF6Coupler(cpl.mm, cpl.ctx, cpl.state, cpl.mf6b,
                              conv_fact=1000.0, mode='lagged')
    assert rebuilt.perlen == [1.0, 2.0, 3.0]


# --------------------------------------------------------------------- #
# the solution check
# --------------------------------------------------------------------- #

def _cpl_with_ws(tmp_path, mfsim='', lst=''):
    cpl, api, ctx = _setup(nper=2)
    cpl.sim_ws = str(tmp_path)
    (tmp_path / 'mfsim.lst').write_text(mfsim)
    (tmp_path / 'toy.lst').write_text(lst)
    return cpl


CLEAN = """
     PERCENT DISCREPANCY =           0.01     PERCENT DISCREPANCY =          -0.00
"""
BAD = """
     PERCENT DISCREPANCY =        -131.69     PERCENT DISCREPANCY =          -0.00
"""


def test_clean_run_passes(tmp_path, capsys):
    cpl = _cpl_with_ws(tmp_path, mfsim='normal termination\n', lst=CLEAN)
    rep = cpl.check_solution()
    assert rep['ok'] and rep['discrepancy'] == pytest.approx(0.01)
    assert 'solution check: converged' in capsys.readouterr().out


def test_budget_discrepancy_fails_the_run(tmp_path):
    """The exact La Mata failure: -131.69% cumulative."""
    cpl = _cpl_with_ws(tmp_path, mfsim='normal termination\n', lst=BAD)
    with pytest.raises(coup.CouplingError, match='not conserving water'):
        cpl.check_solution()


def test_nonconvergence_in_the_listing_fails_the_run(tmp_path):
    cpl = _cpl_with_ws(
        tmp_path,
        mfsim='Solution 1 did not converge for stress period 2 and time step 1\n',
        lst=CLEAN)
    with pytest.raises(coup.CouplingError, match='did not converge'):
        cpl.check_solution()


def test_nonconvergence_seen_through_the_api_fails_the_run(tmp_path):
    cpl = _cpl_with_ws(tmp_path, mfsim='normal termination\n', lst=CLEAN)
    cpl.n_nonconverged = 2
    with pytest.raises(coup.CouplingError, match='did not converge'):
        cpl.check_solution()


def test_failure_can_be_downgraded_to_a_warning(tmp_path, capsys):
    cpl = _cpl_with_ws(tmp_path, mfsim='normal termination\n', lst=BAD)
    rep = cpl.check_solution(raise_on_fail=False)
    assert rep['ok'] is False
    assert 'not conserving water' in capsys.readouterr().out


def test_discrepancy_limit_is_configurable(tmp_path):
    cpl = _cpl_with_ws(tmp_path, mfsim='normal termination\n', lst=BAD)
    assert cpl.check_solution(max_discrepancy=200.0)['ok'] is True


# --------------------------------------------------------------------- #
# ATS in the written model
# --------------------------------------------------------------------- #

def test_ats_block_written_for_transient_periods_only(tmp_path):
    pytest.importorskip('flopy')
    from test_mf6_build import cmf, mf6mod  # noqa: F401
    DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s3L.ini')):
        pytest.skip('La Mata dataset not present')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s3L.ini',
                   xllcorner=739300.0, yllcorner=4553050.0)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]
    b = mf6mod.clsMF6(c, top=np.asarray(c.elev, float), botm=np.asarray(c.botm, float),
                      sim_ws=str(tmp_path), daily=True)
    b.verbose = False
    b.build()
    b.write()
    ats = (tmp_path / 'lamatamm.tdis.ats').read_text()
    periods = [int(l.split()[0]) for l in ats.splitlines()
               if l.strip() and l.strip()[0].isdigit()]
    assert periods == [2, 3, 4], 'the steady-state period must not get ATS'


def test_ats_can_be_disabled(tmp_path):
    pytest.importorskip('flopy')
    from test_mf6_build import cmf, mf6mod  # noqa: F401
    DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s3L.ini')):
        pytest.skip('La Mata dataset not present')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s3L.ini',
                   xllcorner=739300.0, yllcorner=4553050.0)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]
    b = mf6mod.clsMF6(c, top=np.asarray(c.elev, float), botm=np.asarray(c.botm, float),
                      sim_ws=str(tmp_path), daily=True)
    b.verbose = False
    b.ats = False
    b.build()
    b.write()
    assert not (tmp_path / 'lamatamm.tdis.ats').exists()
