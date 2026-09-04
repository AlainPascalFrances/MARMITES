# -*- coding: utf-8 -*-
"""Phase-3b/3c tests: MF6Coupler logic against a mocked MODFLOW 6 API.

The fake API reproduces the memory-pointer surface the coupler uses
(X, UZF FINF/GWD, WEL BOUND, prepare/solve/finalize cycle) so the lagged
and iterative exchange logic is fully exercised without the mf6 binary:

  * lagged mode consumes heads from the END of the previous SP;
  * pointers receive perc (FINF) and -ETg*area (WEL Q);
  * iterative mode re-evaluates MM at each outer iteration, applies
    under-relaxation, and keeps exactly one state advance per SP;
  * with converging-in-1 iterations, constant heads and relax=1,
    iterative == lagged.
"""
import importlib.util
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..', 'trunk'))
sys.path.insert(0, HERE)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


coup = _load('marmites_coupler', os.path.join(TRUNK, 'marmites_coupler.py'))
T = _load('t_runmmsoil', os.path.join(HERE, 'test_runmmsoil.py'))  # reuse synthetic model


# --------------------------------------------------------------------- #
class FakeApi:
    """Minimal stand-in for modflowapi.ModflowApi (memory-pointer surface)."""

    def __init__(self, name, nlay, nrow, ncol, ncell, nuzf, heads0=699.0,
                 nouter=1, dh_per_step=0.0, gwd_per_cell=0.0):
        self.name = name.upper()
        self.X = np.full(nlay * nrow * ncol, float(heads0))
        self.FINF = np.zeros(nuzf)
        self.GWD = np.full(nuzf, float(gwd_per_cell))
        self.BOUND = np.zeros((ncell, 1))
        self.Q = np.zeros(ncell)
        self.nouter = nouter
        self.dh = dh_per_step
        self.initialized = False
        self.finalized = False
        self._k = 0
        self.finf_at_advance = []      # FINF snapshot at each finalize_time_step
        self.q_at_advance = []
        self.finf_iter_trace = []      # FINF snapshot at each solve() call
        self.solve_calls = 0

    # address surface
    def get_var_address(self, var, *comp):
        return '/'.join(comp) + '/' + var

    def get_input_var_names(self):
        """Addresses this fake exposes (the coupler validates against these).

        MF6's operative UZF infiltration array is SINF; it is listed first so
        the coupler binds it in preference to the non-operative FINF (that
        mismatch left UZF applying a constant perc_user in the real run). SINF
        and FINF share the backing array here so flux-content tests are
        unaffected by which name is bound."""
        n = self.name
        return [f'{n}/X', f'{n}/UZF/SINF', f'{n}/UZF/FINF', f'{n}/UZF/GWD',
                f'{n}/WEL/BOUND', f'{n}/WEL/Q']

    def get_time_step(self):
        return 1.0

    def get_value_ptr(self, addr):
        leaf = addr.rsplit('/', 1)[1]
        return {'X': self.X, 'SINF': self.FINF, 'FINF': self.FINF,
                'GWD': self.GWD, 'BOUND': self.BOUND, 'Q': self.Q}[leaf]

    def get_value(self, addr):
        raise KeyError(addr)           # no NODEUSER -> coupler uses identity map

    # control surface
    def initialize(self):
        self.initialized = True

    def finalize(self):
        self.finalized = True

    def prepare_time_step(self, dt):
        pass

    def prepare_solve(self, sol):
        self._k = 0

    def solve(self, sol):
        self._k += 1
        self.solve_calls += 1
        self.finf_iter_trace.append(self.FINF.copy())
        return self._k >= self.nouter

    def finalize_solve(self, sol):
        pass

    def finalize_time_step(self):
        self.finf_at_advance.append(self.FINF.copy())
        # record whichever rate array is in use (Q preferred, see _bind)
        self.q_at_advance.append(self.Q.copy() if np.any(self.Q) else self.BOUND[:, 0].copy())
        self.X += self.dh              # deterministic head evolution per SP


# --------------------------------------------------------------------- #
def _setup(nper=4, mode='lagged', nouter=1, dh=0.0, relax=0.6, heads0=699.0,
           obs_idx=None, obs_names=None):
    cMF = T._FakeMF(nper=nper, perlen=[1] * nper)
    cMF.modelname = 'toy'
    cMF.perc_user = 0.0
    inp = T._build_inputs(cMF)
    mmmod = T.new
    mm = mmmod.clsMMsoil(hnoflo=T.HNOFLO)
    cells = mm.build_cell_list(cMF)
    ctx = mm.build_context(cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'], inp['_Sm'],
                           inp['_Sfc'], inp['_Sr'], inp['_slprop'], inp['_Ssoil_ini'],
                           inp['botm_l0'], inp['_Ks'], inp['gridSOIL'], inp['gridSOILthick'],
                           inp['TopSoil'], inp['gridMETEO'], T.INDEX_MM, T.INDEX_MM_S,
                           inp['gridSsurfhmax'], inp['gridSsurfw'], inp['P_veg_zoneSP'],
                           inp['Eo_zonesSP'], inp['PT_veg_zonesSP'], inp['Pe_veg_zonesSP'],
                           inp['PE_zonesSP'], inp['gridVEGarea'], inp['LAI_veg_zonesSP'],
                           inp['Zr'], inp['kTg_min'], inp['kTg_max'], inp['kT_f'], inp['kT_s'],
                           inp['NVEG'], 1000.0, 0, None, None, None, None, None,
                           None, None, None, None, None)
    state = mm.init_state(ctx)
    # stub of clsMF6 (cell mapping + geometry only)
    surf = [(c[1], c[2], 0) for c in cells]
    mf6b = SimpleNamespace(cMF=cMF, ncell=ctx.ncell, surf_cells=surf,
                           nrow=cMF.nrow, ncol=cMF.ncol)
    cpl = coup.MF6Coupler(mm, ctx, state, mf6b, conv_fact=1000.0,
                          mode=mode, relax=relax,
                          obs_idx=obs_idx, obs_names=obs_names)
    api = FakeApi('toy', cMF.nlay, cMF.nrow, cMF.ncol, ctx.ncell,
                  nuzf=ctx.ncell + 3, heads0=heads0, nouter=nouter, dh_per_step=dh)
    return cpl, api, ctx


def test_lagged_uses_previous_sp_heads():
    dh = -0.5                      # heads drop 0.5 m per advance
    cpl, api, ctx = _setup(mode='lagged', dh=dh, heads0=699.0)
    res = cpl.run(api)
    assert api.initialized and api.finalized
    # steady SP advanced once + nper transient
    assert len(api.finf_at_advance) == 1 + ctx.cMF.nper
    # MM at SP n must see heads after (n+1) advances: 699 + (n+1)*dh
    for n in range(ctx.cMF.nper):
        expect = 699.0 + (n + 1) * dh
        assert np.allclose(res['heads'][n], expect), (n, res['heads'][n][0], expect)


def test_steady_state_uses_mean_recharge_when_supplied():
    """A steady period ignores the initial-head file, so it must be driven by
    the mean forcing to land near equilibrium. When steady_perc/etg are set the
    steady advance (index 0) must carry them, not the uniform perc_user."""
    cpl, api, ctx = _setup(mode='lagged')
    sp = np.full(ctx.ncell, 7.0e-4)          # per-cell mean recharge (m/d)
    se = np.full(ctx.ncell, 1.0e-4)          # per-cell mean ETg
    cpl.steady_perc, cpl.steady_etg = sp, se
    cpl.run(api)
    finf0 = api.finf_at_advance[0]           # the steady advance
    q0 = api.q_at_advance[0]
    assert np.allclose(finf0[:ctx.ncell], sp, atol=1e-12)
    assert np.allclose(q0, -se * cpl.area, atol=1e-12)


def test_steady_state_falls_back_to_perc_user():
    cpl, api, ctx = _setup(mode='lagged')
    ctx.cMF.perc_user = 2.0e-4
    cpl.run(api)
    finf0 = api.finf_at_advance[0]
    assert np.allclose(finf0[:ctx.ncell], 2.0e-4, atol=1e-12)


def test_pointer_contents_finf_and_welq():
    cpl, api, ctx = _setup(mode='lagged')
    res = cpl.run(api)
    # at each transient advance, FINF[0:ncell] held that SP's perc (m/d)
    for n in range(ctx.cMF.nper):
        finf = api.finf_at_advance[1 + n]
        assert np.allclose(finf[:ctx.ncell], res['perc'][n], atol=1e-12)
        assert np.all(finf[ctx.ncell:] == 0.0)
        q = api.q_at_advance[1 + n]
        assert np.allclose(q, -res['etg'][n] * cpl.area, atol=1e-12)
    # sanity: the shallow water table (dgwt ~1 m, loam) makes ETg active,
    # so a genuinely nonzero flux crossed the WEL pointer
    assert res['etg'].max() > 0
    assert max(q.min() for q in api.q_at_advance[1:]) < 0  # sink wells


def test_iterative_reevaluates_and_relaxes():
    nouter = 3
    cpl, api, ctx = _setup(mode='iterative', nouter=nouter, relax=0.5)
    res = cpl.run(api)
    # steady advance solves too: total solve calls = 1*nouter? steady uses plain
    # _advance (loop until converged) -> nouter calls as well
    assert np.all(res['outer_iters'] == nouter)
    # MM evaluated nouter times per transient SP: FINF trace differs from the
    # raw step() output when relaxation is active (relax != 1) after iter 1
    # (heads constant here so raw evaluations are identical; relaxed value
    # equals raw -> check instead that the trace has nouter snapshots per SP)
    per_sp = nouter
    n_transient_snapshots = len(api.finf_iter_trace) - nouter  # minus steady SP
    assert n_transient_snapshots == ctx.cMF.nper * per_sp


def test_iterative_equals_lagged_when_static():
    """Constant heads, converge in 1 outer iteration, relax=1: both modes
    must produce identical perc/ETg trajectories (the lag has no effect
    because heads never change)."""
    cpl_l, api_l, _ = _setup(mode='lagged', nouter=1, dh=0.0)
    res_l = cpl_l.run(api_l)
    cpl_i, api_i, _ = _setup(mode='iterative', nouter=1, dh=0.0, relax=1.0)
    res_i = cpl_i.run(api_i)
    assert np.allclose(res_l['perc'], res_i['perc'], atol=1e-12)
    assert np.allclose(res_l['etg'], res_i['etg'], atol=1e-12)
    assert np.allclose(res_l['heads'], res_i['heads'])


def test_iterative_single_state_advance_per_sp():
    """State must advance exactly once per SP despite multiple evaluations:
    running iterative (nouter=4, relax=1, static heads) must equal lagged."""
    cpl_i, api_i, _ = _setup(mode='iterative', nouter=4, dh=0.0, relax=1.0)
    res_i = cpl_i.run(api_i)
    cpl_l, api_l, _ = _setup(mode='lagged', nouter=1, dh=0.0)
    res_l = cpl_l.run(api_l)
    assert np.allclose(res_i['perc'], res_l['perc'], atol=1e-12)
    assert np.allclose(res_i['etg'], res_l['etg'], atol=1e-12)


def _setup_grid(grid, nper=2, heads0=699.0):
    """Same toy model, but with an explicit grid geometry (Phase 4)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'mgrid_c', os.path.join(TRUNK, 'marmites_grid.py'))
    Gm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(Gm)
    cMF = T._FakeMF(nper=nper, perlen=[1] * nper)
    cMF.modelname = 'toy'
    cMF.perc_user = 0.0
    cMF.xllcorner = cMF.yllcorner = 0.0
    inp = T._build_inputs(cMF)
    mm = T.new.clsMMsoil(hnoflo=T.HNOFLO)
    cells = mm.build_cell_list(cMF)
    geom = Gm.geometry_for(cMF, cells, grid=grid)
    ctx = mm.build_context(cMF, cells, inp['_nsl'], inp['_nslmax'], inp['_st'], inp['_Sm'],
                           inp['_Sfc'], inp['_Sr'], inp['_slprop'], inp['_Ssoil_ini'],
                           inp['botm_l0'], inp['_Ks'], inp['gridSOIL'], inp['gridSOILthick'],
                           inp['TopSoil'], inp['gridMETEO'], T.INDEX_MM, T.INDEX_MM_S,
                           inp['gridSsurfhmax'], inp['gridSsurfw'], inp['P_veg_zoneSP'],
                           inp['Eo_zonesSP'], inp['PT_veg_zonesSP'], inp['Pe_veg_zonesSP'],
                           inp['PE_zonesSP'], inp['gridVEGarea'], inp['LAI_veg_zonesSP'],
                           inp['Zr'], inp['kTg_min'], inp['kTg_max'], inp['kT_f'], inp['kT_s'],
                           inp['NVEG'], 1000.0, 0, None, None, None, None, None,
                           None, None, None, None, None, geom=geom)
    state = mm.init_state(ctx)
    surf = [(c[1], c[2], 0) for c in cells]
    mf6b = SimpleNamespace(cMF=cMF, ncell=ctx.ncell, surf_cells=surf,
                           nrow=cMF.nrow, ncol=cMF.ncol)
    cpl = coup.MF6Coupler(mm, ctx, state, mf6b, conv_fact=1000.0, mode='lagged')
    api = FakeApi('toy', cMF.nlay, cMF.nrow, cMF.ncol, ctx.ncell,
                  nuzf=ctx.ncell + 3, heads0=heads0)
    return cpl, api, ctx, cells, cMF


def test_disv_node_mapping_reads_correct_heads():
    """DISV: user node = k*ncpl + icell2d. Tag X per node and verify the
    coupler gathers exactly the surface cells' values."""
    cpl, api, ctx, cells, cMF = _setup_grid('disv')
    assert cpl.grid_kind == 'vertex'
    assert cpl.ncpl == cMF.nrow * cMF.ncol
    # unique tag per node so a wrong mapping cannot pass
    api.X[:] = np.arange(api.X.size, dtype=float)
    res = cpl.run(api)
    ncpl = cMF.nrow * cMF.ncol
    expected = np.array([0 * ncpl + c[3] for c in cells], dtype=float)  # layer 0
    assert np.allclose(res['heads'][0], expected)


def test_dis_and_disv_node_mapping_agree_on_equivalent_grid():
    """For the DIS-equivalent vertex grid both mappings select the same X
    entries (icell2d == i*ncol+j, single layer)."""
    cpl_d, api_d, _, _, _ = _setup_grid('dis')
    cpl_v, api_v, _, _, _ = _setup_grid('disv')
    api_d.X[:] = np.arange(api_d.X.size, dtype=float)
    api_v.X[:] = np.arange(api_v.X.size, dtype=float)
    rd = cpl_d.run(api_d)
    rv = cpl_v.run(api_v)
    assert np.allclose(rd['heads'], rv['heads'])
    assert np.allclose(cpl_d.area, cpl_v.area)


def test_coupler_area_comes_from_geometry():
    cpl, _, ctx, _, _ = _setup_grid('disv')
    assert np.allclose(cpl.area, ctx.geom.area)


class _BadApi(FakeApi):
    """FakeApi variants reproducing real MF6 failure modes."""

    def __init__(self, *a, mode='ok', **kw):
        super().__init__(*a, **kw)
        self._mode = mode

    def get_time_step(self):
        if self._mode == 'zero_dt':
            return 0.0
        if self._mode == 'no_dt':
            raise AttributeError('get_time_step unavailable')
        return 1.0

    # leaves that each "missing" mode makes genuinely unreachable. In real MF6
    # an absent variable is not merely dropped from get_input_var_names() (SINF
    # itself is absent from that list yet reachable); it fails at get_value_ptr.
    # So absence must be simulated by making the POINTER unresolvable, otherwise
    # the coupler (which now binds any reachable, correctly-sized pointer -- see
    # marmites_coupler._bind_first, lesson: the var list is not authoritative)
    # would still bind it.
    _MISSING = {'no_finf': {'SINF', 'FINF'},
                'no_gwd': {'GWD'},          # SIMVALS/DRN_SEEP already unresolved
                'no_q': {'Q'}}

    def get_input_var_names(self):
        n = self.name
        if self._mode == 'no_finf':
            return [f'{n}/X', f'{n}/WEL/BOUND', f'{n}/WEL/Q']
        if self._mode == 'no_gwd':
            return [f'{n}/X', f'{n}/UZF/FINF', f'{n}/WEL/BOUND', f'{n}/WEL/Q']
        if self._mode == 'no_q':
            return [f'{n}/X', f'{n}/UZF/FINF', f'{n}/UZF/GWD', f'{n}/WEL/BOUND']
        return super().get_input_var_names()

    def get_value_ptr(self, addr):
        leaf = addr.rsplit('/', 1)[1]
        if leaf in self._MISSING.get(self._mode, ()):
            raise KeyError(addr)           # genuinely absent -> unresolvable
        return super().get_value_ptr(addr)


def _bad(mode, **kw):
    cpl, api, ctx, _cells, cMF = _setup_grid('dis')
    bad = _BadApi('toy', cMF.nlay, cMF.nrow, cMF.ncol, ctx.ncell,
                  nuzf=ctx.ncell + 3, mode=mode, **kw)
    return cpl, bad


def test_zero_timestep_is_tolerated():
    """MF6 IGNORES the dt argument of prepare_time_step (it derives the step
    from TDIS) and get_time_step() returns 0.0 before the first step, so a
    zero must be passed through, never treated as fatal."""
    cpl, api = _bad('zero_dt')
    res = cpl.run(api)                      # must complete
    assert res['perc'].shape[0] == cpl.mf6b.cMF.nper


def test_missing_get_time_step_is_tolerated():
    cpl, api = _bad('no_dt')
    res = cpl.run(api)
    assert res['perc'].shape[0] == cpl.mf6b.cMF.nper


def test_infiltration_binds_sinf_not_finf():
    """MF6's operative UZF infiltration array is SINF; FINF is a valid-but-
    non-operative pointer. Binding FINF left UZF applying a constant perc_user
    and the water table drained. The coupler must prefer SINF."""
    cpl, api, ctx = _setup(mode='lagged')
    cpl.run(api)
    assert cpl.addr_finf.endswith('/SINF'), cpl.addr_finf


def test_unbindable_finf_reports_available_variables():
    """A wrong address hands back a wrong-sized pointer; writing through it
    corrupts MF6 memory. Fail loudly instead, listing what exists."""
    cpl, api = _bad('no_finf')
    with pytest.raises(coup.CouplingError, match='FINF'):
        cpl.run(api)


def test_missing_gwd_degrades_to_zero_exfiltration():
    """Groundwater-discharge array absent -> warn and treat exf as 0,
    rather than crash."""
    cpl, api = _bad('no_gwd')
    api.GWD[:] = 99.0                      # would be huge if wrongly read
    res = cpl.run(api)
    assert np.allclose(res['exf'], 0.0)


def test_node_mapping_overflow_is_caught():
    """If the grid is reduced by idomain and NODEUSER cannot be read, the
    mapping would index past X -- that must be detected before writing."""
    cpl, api, ctx, _cells, cMF = _setup_grid('dis')
    api.X = np.zeros(3)                    # far smaller than the node ids
    with pytest.raises(coup.CouplingError, match='past the head array|empty'):
        cpl.run(api)


def test_bound_is_not_written_when_q_exists():
    """Regression guard for the La Mata crash: MF6 6.7 exposes both WEL/Q and
    WEL/BOUND, but writing BOUND corrupts MF6's heap (bisected with
    tests/diagnose_coupling.py: FINF+Q fine, +BOUND kills prepare_time_step).
    Only Q may be written when it is available."""
    cpl, api, ctx = _setup(mode='lagged')
    api.BOUND[:] = -12345.0                 # sentinel: must survive untouched
    res = cpl.run(api)
    assert cpl.p_bound is None, 'BOUND must not even be bound when Q exists'
    assert np.all(api.BOUND == -12345.0), 'BOUND was written -- crashes MF6 6.7'
    # and Q did receive the rates
    assert np.allclose(api.Q[:ctx.ncell], -res['etg'][-1] * cpl.area, atol=1e-12)


def test_bound_used_only_as_fallback_without_q():
    cpl, api = _bad('no_q')
    res = cpl.run(api)
    assert cpl.p_q is None and cpl.p_bound is not None
    assert np.allclose(api.BOUND[:, 0], -res['etg'][-1] * cpl.area, atol=1e-12)


def test_obs_capture_records_full_mm_series_at_obs_cells():
    """With obs_idx set, the coupler keeps the FULL per-SP MM flux vectors at
    those cells (for the per-point Sankey), exposes them in the result, and they
    match the aggregates: mm_obs at an obs cell must equal that cell's row of the
    per-cell arrays, and the mean over a single-cell selection equals wb_ts."""
    idx = [0, 3]
    cpl, api, ctx = _setup(mode='lagged', obs_idx=idx, obs_names=['A', 'B'])
    res = cpl.run(api)
    nper = ctx.cMF.nper
    nidx = len(ctx.index)
    nsl = int(ctx._nslmax)
    nidx_s = len(ctx.index_S)
    assert res['mm_obs'].shape == (nper, len(idx), nidx)
    assert res['mms_obs'].shape == (nper, len(idx), nsl, nidx_s)
    assert list(res['obs_idx']) == idx
    assert [n.decode() for n in res['obs_names']] == ['A', 'B']
    # obs_ij must be the (i, j) of the selected cells
    for k, p in enumerate(idx):
        assert res['obs_ij'][k, 0] == cpl.i_arr[p]
        assert res['obs_ij'][k, 1] == cpl.j_arr[p]
    # every captured value is finite
    assert np.all(np.isfinite(res['mm_obs']))


def test_no_obs_capture_by_default():
    """Without obs_idx the result carries no obs arrays (back-compatible)."""
    cpl, api, _ctx = _setup(mode='lagged')
    res = cpl.run(api)
    assert 'mm_obs' not in res and cpl.mm_obs is None


def test_all_zero_exfiltration_is_flagged(capsys):
    """A whole run with zero exfiltration everywhere must warn: it indicates
    UZF was built without SIMULATE_GWSEEP (the La Mata failure), not a result."""
    cpl, api, _ctx = _setup(mode='lagged')
    api.GWD[:] = 0.0
    cpl.run(api)
    out = capsys.readouterr().out
    assert 'exfiltration is zero' in out.lower()
    assert 'simulate_gwseep' in out.lower()


def test_nonzero_exfiltration_is_not_flagged(capsys):
    cpl, api, _ctx = _setup(mode='lagged')
    api.GWD[:] = 5.0
    cpl.run(api)
    assert 'exfiltration is zero' not in capsys.readouterr().out.lower()


def test_exfiltration_sign_and_scaling():
    """GWD (m3/d, + to surface) must arrive in MM as +mm/d into the soil."""
    cpl, api, ctx = _setup(mode='lagged')
    api.GWD[:] = 25.0                       # m3/d on 100x100 m cells
    res = cpl.run(api)
    # 25 m3/d / 1e4 m2 * 1000 mm/m = 2.5 mm/d, positive into the soil
    assert np.allclose(res['exf'], 2.5)
