# -*- coding: utf-8 -*-
"""The non-convergence guard (user, 2026-10-20): MF6 stalling is a sign the
solution is not valid, and the run stops instead of crawling for hours.

It counts NEWTON (outer) iterations per stress period, the tries ATS gave
up included, and stops on a stall that LASTS: [run] stall_periods periods
in a row past stall_iterations, or one past max_period_iterations. La
Mata's base run of 2026-10-19 needed at most 31 in a day; realisation 67
got through two hard days (937 and 777) before it stalled nine days in a
row (870-2700 each: a water table rising into a saturated UZF wetting
front, MF6's uz_rise without a Newton derivative). The first guard -- on
linear iterations, one day -- killed six realisations at their first hard
day and one that had converged after an expensive linear solve."""

import inspect
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in (TRUNK, os.path.join(TRUNK, 'tools')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg                                  # noqa: E402
import marmites_coupler as mc                                   # noqa: E402


class _Stepper(mc.MF6Coupler):
    """_advance over scripted (Newton iterations, step length) sub-steps,
    MF6's clock moving by the step length."""

    def __init__(self, stall=300, periods=3, cap=2000):
        self._steps = []
        self._t = 0.0
        self.max_substeps = 5000
        self.stall_iterations = stall
        self.stall_periods = periods
        self.max_period_iterations = cap
        self._stalled = []
        self.substeps = self.n_nonconverged = 0
        self.ctx = None

    def period(self, steps, n):
        """One stress period (0-based ``n``) of one day."""
        self._steps = list(steps)
        return self._advance(None, t_end=self._t + 1.0, period=n)

    def _now(self, api):
        return self._t

    def _one_step(self, api, write_cb=None):
        k, dt = self._steps.pop(0)
        self._t += dt
        return k, k, True

    def _keep_sfr_inflow(self):
        pass

    def _rates_add(self, acc, dt):
        pass

    def _rates_close(self, acc, nsub):
        pass


def test_a_stall_that_lasts_stops_the_run_within_its_last_period():
    cpl = _Stepper()
    cpl.period([(60, 1.0)], 245)
    cpl.period([(871, 1.0)], 246)                    # stalled, 1 of 3
    cpl.period([(100, 0.25)] * 2 + [(1016, 0.5)], 247)  # 1216: 2 of 3
    with pytest.raises(mc.NonConvergence) as e:
        cpl.period([(90, 0.05)] * 20, 248)           # 1800 in the day
    msg = str(e.value)
    assert 'NON-CONVERGENCE: stress period 249 stalls, 3 period(s) in a row'\
        in msg
    assert 'past 300 Newton iterations (871, 1216, 360;' in msg
    assert 'run.stall_iterations, run.stall_periods' in msg
    # stopped within the period, as soon as it passed the limit
    assert len(cpl._steps) == 20 - 4                 # 4 x 90 = 360


def test_a_hard_day_passes_and_a_quiet_one_resets_the_count():
    cpl = _Stepper()
    # realisation 67: 937 on one day, 777 eleven days later, then quiet
    for n, k in enumerate([60, 937, 70, 777, 50, 400, 420, 30, 500]):
        cpl.period([(k, 1.0)], n)
    assert cpl._stalled == [500]
    # the base run's worst day, and the guard off
    assert _Stepper().period([(31, 1.0)], 922) == 31
    off = _Stepper(stall=0, cap=0)
    for n in range(5):
        off.period([(90, 0.025)] * 40, n)            # 3600 a day
    assert off.substeps == 5 * 39


def test_one_period_past_the_cap_stops_at_once():
    cpl = _Stepper()
    with pytest.raises(mc.NonConvergence) as e:
        cpl.period([(100, 0.04)] * 25, 1811)
    msg = str(e.value)
    assert 'stress period 1812 needed more than 2000 Newton' in msg
    assert '(2100 in 21 ATS sub-step(s)' in msg
    assert 'run.max_period_iterations' in msg


def test_one_period_in_a_row_is_a_plain_limit():
    cpl = _Stepper(periods=1, cap=0)
    with pytest.raises(mc.NonConvergence) as e:
        cpl.period([(301, 1.0)], 0)
    assert '1 period(s) in a row past 300 Newton iterations (301;' \
        in str(e.value)


def test_the_count_takes_the_tries_ats_gave_up():
    # every try resets IOUTTOT_TIMESTEP; each failed one ran MXITER
    # iterations and cut the step by dtfailadj (5): 1 -> 0.2 -> 0.04 d
    c = mc.MF6Coupler.__new__(mc.MF6Coupler)
    c.p_iouttot = np.array([7], dtype=np.int32)
    c.p_mxiter = np.array([100], dtype=np.int32)
    c.ats_dtfailadj, c.max_outer = 5.0, 500
    assert c._step_newton(1.0, 0.04) == 7 + 2 * 100
    assert c._step_newton(1.0, 1.0) == 7
    assert c._step_newton(0.5, 0.1) == 107
    assert c._step_newton(None, None) == 7
    c.p_mxiter = None                                   # the solver's own
    assert c._step_newton(1.0, 0.2) == 507
    c.p_iouttot = None                                  # not exposed
    assert c._step_newton(1.0, 0.04) == 0


def test_the_guard_is_a_run_field_with_its_checks():
    c = mcfg.RunConfig.from_dict({})
    g = c.run
    # 300 and a cap of 2000 stopped slow days that converged and single
    # hard days (user, 2026-10-20): 600, no cap
    assert (g.stall_iterations, g.stall_periods,
            g.max_period_iterations) == (600, 3, 0)
    assert not hasattr(c.pest, 'max_solver_iterations')
    assert not any('stall' in e or 'max_period' in e for e in c.problems())
    for k, v, m in (('stall_iterations', -1, 'run.stall_iterations must'),
                    ('stall_periods', 0, 'run.stall_periods must be >= 1'),
                    ('max_period_iterations', -2,
                     'run.max_period_iterations must be >= 0'),
                    ('max_period_iterations', 600,
                     'run.max_period_iterations (600) must be above')):
        c = mcfg.RunConfig.from_dict({})
        setattr(c.run, k, v)
        assert any(m in e for e in c.problems()), (k, v)
    c = mcfg.RunConfig.from_dict({})
    c.run.stall_iterations = c.run.max_period_iterations = 0   # off
    assert not any('stall' in e or 'max_period' in e for e in c.problems())


def test_a_file_with_the_first_guard_still_loads():
    c = mcfg.RunConfig.from_dict({'run': {'max_solver_iterations': 0},
                                  'pest': {'max_solver_iterations': 1500}})
    assert c.run.stall_iterations == 600
    assert 'max_solver_iterations' in mcfg.RETIRED['run']
    assert 'max_solver_iterations' in mcfg.RETIRED['pest']


def test_the_driver_passes_it_on_and_the_set_up_keeps_it():
    sys.path.insert(0, os.path.join(TRUNK, 'app'))
    from lib import runs as runlib
    runner = runlib.import_runner()
    src = inspect.getsource(runner)
    for k in ('stall_iterations', 'stall_periods', 'max_period_iterations'):
        assert 'cpl.%s = int(a.%s)' % (k, k) in src
        assert '%s=int(cfg.run.%s)' % (k, k) in src
    import pest_lamata_setup as pls
    assert 'max_solver_iterations' not in inspect.getsource(pls)
