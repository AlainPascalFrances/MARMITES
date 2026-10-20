# -*- coding: utf-8 -*-
"""The non-convergence guard (user, 2026-10-20): a stress period needing
more linear-solver iterations than [run] max_solver_iterations stops the
run -- such a solution is not valid -- instead of crawling for hours. La
Mata's base run of 2026-10-19 took at most 133 in a day; realisation 67
took 3000-4900 on each day it stalled (a water table rising into a
saturated UZF wetting front, MF6's uz_rise without a Newton derivative)."""

import inspect
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in (TRUNK, os.path.join(TRUNK, 'tools')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg                                  # noqa: E402
import marmites_coupler as mc                                   # noqa: E402


class _Stepper(mc.MF6Coupler):
    """_advance over a scripted list of (iterations, step length): one
    sub-step each, MF6's clock moving by the step length."""

    def __init__(self, steps, limit):
        self._steps = list(steps)
        self._t = 0.0
        self.max_substeps = 5000
        self.max_solver_iterations = limit
        self.substeps = self.n_nonconverged = 0
        self.ctx = None

    def _now(self, api):
        return self._t

    def _one_step(self, api, write_cb=None):
        k, dt = self._steps.pop(0)
        self._t += dt
        return k, True

    def _keep_sfr_inflow(self):
        pass

    def _rates_add(self, acc, dt):
        pass

    def _rates_close(self, acc, nsub):
        pass


def test_a_stress_period_past_the_limit_stops_the_run():
    # 40 sub-steps of 0.025 d at 90 iterations: 3600 in the day
    cpl = _Stepper([(90, 0.025)] * 40, limit=1500)
    with pytest.raises(mc.NonConvergence) as e:
        cpl._advance(None, t_end=1.0, period=1811)
    msg = str(e.value)
    assert 'NON-CONVERGENCE: stress period 1812' in msg
    assert 'more than 1500' in msg and 'run.max_solver_iterations' in msg
    # stopped within the period, as soon as the sum passed the limit
    assert len(cpl._steps) == 40 - 17                    # 17 x 90 = 1530


def test_a_healthy_day_and_no_limit_pass():
    cpl = _Stepper([(133, 1.0)], limit=1500)             # the base's worst
    assert cpl._advance(None, t_end=1.0, period=922) == 133
    cpl = _Stepper([(90, 0.025)] * 40, limit=0)          # 0: no limit
    cpl._advance(None, t_end=1.0, period=1811)
    assert cpl.substeps == 39


def test_the_limit_is_a_run_field_and_the_calibrations_own():
    c = mcfg.RunConfig.from_dict({})
    assert c.run.max_solver_iterations == 0               # no limit
    assert c.pest.max_solver_iterations == 1500
    c.run.max_solver_iterations = -1
    assert any('run.max_solver_iterations must be >= 0' in e
               for e in c.problems())
    c = mcfg.RunConfig.from_dict({})
    c.pest.max_solver_iterations = -5
    assert any('pest.max_solver_iterations must be >= 0' in e
               for e in c.problems())


def test_the_driver_and_the_set_up_pass_it_on():
    sys.path.insert(0, os.path.join(TRUNK, 'app'))
    from lib import runs as runlib
    runner = runlib.import_runner()
    src = inspect.getsource(runner)
    assert 'cpl.max_solver_iterations = int(a.max_solver_iterations)' in src
    assert 'max_solver_iterations=int(cfg.run.max_solver_iterations)' in src
    import pest_lamata_setup as pls
    assert ('cfg.run.max_solver_iterations = int(cfg.pest.'
            'max_solver_iterations)') in inspect.getsource(pls.step_org)
