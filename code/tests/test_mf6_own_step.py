# -*- coding: utf-8 -*-
"""The coupler advances MF6 with MF6's OWN time step, so a step the solver
cannot converge is retried with a smaller one (ATS) instead of accepted.

Driven through prepare_solve / solve / finalize_solve, MF6 reports a failed
step ("did not converge") and ACCEPTS it: the ATS retry lives only in its
own do_time_step (Mf6DoTimestep, mf6core.f90). On 2026-09-27 La Mata's
stream cells oscillated by 14-69 m per iteration on some wet days, and those
days went into the record. do_time_step re-applies the period input on every
try -- uzf_ad copies SINF_PVAR to FINF/SINF and PET_PVAR to PET/PETMAX -- so
the coupler writes those too. Checked against the real libmf6 on a toy
(E:/tmp_claude/helpers/ownstep_toy.py): the value survives, and a hopeless
step is retried by ATS and surfaces from finalize_time_step as an error.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


M = _load('test_coupler_mock_os', os.path.join(HERE, 'test_coupler_mock.py'))


class OwnStepApi(M.FakeApi):
    """The fake with MF6's own time step: the period-input arrays exposed,
    do_time_step running ad (PVAR -> FINF/SINF, PET) on every TRY, and a
    step that can be made to fail its first tries (ATS retries) or all of
    them (finalize_time_step then reports failure)."""

    def __init__(self, *a, fail_tries=None, hopeless_steps=(), **kw):
        super().__init__(*a, **kw)
        self.fail_tries = dict(fail_tries or {})   # step index -> failed tries
        self.hopeless = set(hopeless_steps)
        self.ITERTOT = np.zeros(1, dtype=np.int32)
        self.step_no = -1
        self.tries = []
        self.finf_each_try = []
        self.do_calls = 0
        self._failed = False
        self.prepare_solve_calls_outside = 0

    def get_input_var_names(self):
        n = self.name
        return super().get_input_var_names() + [
            f'{n}/UZF/SINF_PVAR', f'{n}/UZF/PET_PVAR', 'SLN_1/ITERTOT_TIMESTEP']

    def get_value_ptr(self, addr):
        leaf = addr.rsplit('/', 1)[1]
        if leaf == 'SINF_PVAR':
            return self.SINF_PVAR
        if leaf == 'PET_PVAR':
            return self.PET_PVAR
        if leaf == 'ITERTOT_TIMESTEP':
            return self.ITERTOT
        return super().get_value_ptr(addr)

    def prepare_time_step(self, dt):
        self.step_no += 1
        self._failed = False

    def do_time_step(self):
        self.do_calls += 1
        k = self.step_no
        fails = 10 ** 6 if k in self.hopeless else self.fail_tries.get(k, 0)
        tries = 0
        while True:
            tries += 1
            self.prepare_solve(1)             # ad: PVAR -> FINF/SINF, PET
            self.finf_each_try.append((k, tries, self.FINF.copy()))
            for _ in range(self.nouter):
                self.solve(1)
            self.ITERTOT[0] = 3 * tries
            if tries > fails:
                break
            if tries >= 4:                    # ATS reached dtmin: give up
                self._failed = True
                break
        self.tries.append(tries)

    def finalize_time_step(self):
        super().finalize_time_step()
        if self._failed:
            raise RuntimeError('simulation failed to converge')


def _setup(api_kw=None, nper=4):
    cpl, _api, ctx = M._setup(nper=nper, mode='lagged')
    api = OwnStepApi('toy', ctx.cMF.nlay, ctx.cMF.nrow, ctx.cMF.ncol,
                     ctx.ncell, nuzf=ctx.ncell + 3, **(api_kw or {}))
    return cpl, api, ctx


def test_the_coupler_uses_mf6s_own_time_step_when_it_can():
    cpl, api, ctx = _setup()
    res = cpl.run(api)
    assert cpl.mf6_step
    assert api.do_calls == 1 + ctx.cMF.nper          # steady + transient
    assert res['iters_kind'] == 'linear'
    assert res['nonconverged'] == 0


def test_the_written_percolation_survives_every_retry():
    """ad re-applies the period input on each try: FINF must be MMsoil's
    percolation on the retries too, not the build-time rate."""
    cpl, api, ctx = _setup({'fail_tries': {2: 2}})   # SP 1: two failed tries
    res = cpl.run(api)
    tries = [t for (k, t, _f) in api.finf_each_try if k == 2]
    assert tries == [1, 2, 3], tries
    for k, t, finf in api.finf_each_try:
        if k == 0:
            continue                                  # the steady period
        assert np.allclose(finf[:ctx.ncell], res['perc'][k - 1]), (k, t)
        assert not np.any(finf[:ctx.ncell] == -1.0), 'the input file rate'
    assert res['nonconverged'] == 0


def test_a_step_mf6_cannot_solve_after_retrying_is_counted_not_fatal():
    cpl, api, ctx = _setup({'hopeless_steps': {2}})
    res = cpl.run(api)
    assert res['nonconverged'] == 1
    assert api.do_calls == 1 + ctx.cMF.nper          # the run went on
    rep = cpl.check_solution(raise_on_fail=False)
    assert not rep['ok']
    assert rep['nonconverged_api'] == 1


def test_the_uzf_demand_goes_to_the_period_input_too():
    cpl, api, ctx = _setup()
    cpl.run(api)
    assert api.PET_PVAR[:ctx.ncell].min() >= 0.0
    assert np.allclose(api.PET_PVAR[:ctx.ncell], cpl._petmax_written)


def test_without_the_period_input_the_old_stepping_stays():
    """An MF6 build that does not expose SINF_PVAR keeps the prepare_solve
    path -- and says a failed step would not be retried."""
    cpl, api, ctx = M._setup(nper=3, mode='lagged')
    cpl.run(api)
    assert not cpl.mf6_step
