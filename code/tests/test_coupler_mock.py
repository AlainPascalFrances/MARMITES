# -*- coding: utf-8 -*-
"""Phase-3b/3c tests: MF6Coupler logic against a mocked MODFLOW 6 API.

The fake API reproduces the memory-pointer surface the coupler uses
(X, UZF FINF/GWD/PETMAX, the two EVT packages' curves, the
prepare/solve/finalize cycle) so the exchange logic is fully exercised
without the mf6 binary:

  * the coupling is lagged: MM consumes heads from the END of the previous
    SP (the iterative mode went on 2026-10-07);
  * pointers receive perc (FINF) and the day's EVT curves; the fake EVT
    takes their full RATE, the worst case for total ET against PET;
  * WEL is never bound or written: the ETg wells went with the WEL route
    (2026-10-07), and WEL is left for real pumping.
"""
import importlib.util
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
# the coupler draws the EVT curves with marmites_evt; the progress tests
# import marmites_coupler by name (they passed only after another module
# had put TRUNK on the path)
for _p in (TRUNK, os.path.join(TRUNK, 'ppMF6')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


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
                 nouter=1, dh_per_step=0.0, gwd_per_cell=0.0, nseg=8):
        self.name = name.upper()
        self.X = np.full(nlay * nrow * ncol, float(heads0))
        # MF6's UZF cell group keeps TWO infiltration arrays, both set from
        # the period input SINF_PVAR by uzf_ad (setdatafinf), which
        # xmi_prepare_solve runs: FINF is what the kinematic wave ROUTES
        # (UzfCellGroup.f90, surflux = finf + mover), SINF what the budget's
        # INFILTRATION line REPORTS (gwf-uzf.f90, appliedinf). They are
        # separate arrays here too: when this fake shared one backing array
        # for both, writing only SINF passed every test -- and on the real
        # model UZF routed the steady input rate all run (2026-09-23).
        self.SINF_PVAR = np.full(nuzf, -1.0)   # the input file's rate
        self.FINF = np.zeros(nuzf)
        self.SINF = np.zeros(nuzf)
        self.GWD = np.full(nuzf, float(gwd_per_cell))
        # WP2: UZF's PET, as MF6 keeps it -- uzf_ad (in prepare_solve) sets
        # PET and PETMAX from the period input PET_PVAR, and every solve
        # iteration resets PET from PETMAX: PETMAX is the operative array.
        # UZET is the actual ET per object [m3/d], a constant per object here.
        self.PET_PVAR = np.full(nuzf, -1.0)
        self.PET = np.zeros(nuzf)
        self.PETMAX = np.zeros(nuzf)
        self.UZET = np.full(nuzf, -2.0)       # MF6 reports it negative (out)
        # when set (land-object areas, m2): UZF takes ALL of its demand, the
        # worst case for total ET against PET
        self.uzet_area = None
        self.uzet_extra = 0.0                 # m/d taken BEYOND the demand
        self.petmax_at_advance = []
        self.pet_used = []                    # PET as the solve saw it
        # GROUNDWATER ET: the two EVT packages' arrays, by package. MF6 keeps
        # PXDP/PETM (nseg - 1) per record. The fake takes the curve's full
        # RATE [m/d] on the cell's area (evt_area, set by _setup) -- the
        # most EVT can take, so the worst case for total ET against PET.
        self.EVT = {pk: {'SURFACE': np.zeros(ncell), 'RATE': np.zeros(ncell),
                         'DEPTH': np.ones(ncell),
                         'PXDP': np.zeros((ncell, nseg - 1)),
                         'PETM': np.zeros((ncell, nseg - 1)),
                         'SIMVALS': np.zeros(ncell)}
                    for pk in ('EVT_EG', 'EVT_TG')}
        self.evt_area = None
        self.evt_at_advance = []      # {pkg: RATE} at each finalize_time_step
        # a WEL the coupler must never touch (real pumping, some day)
        self.BOUND = np.zeros((ncell, 1))
        self.Q = np.zeros(ncell)
        self.nouter = nouter
        self.dh = dh_per_step
        self.initialized = False
        self.finalized = False
        self._k = 0
        self.finf_at_advance = []      # FINF snapshot at each finalize_time_step
        self.sinf_at_advance = []      # SINF snapshot at each finalize_time_step
        self.finf_iter_trace = []      # FINF snapshot at each solve() call
        self.solve_calls = 0

    # address surface
    def get_var_address(self, var, *comp):
        return '/'.join(comp) + '/' + var

    def get_input_var_names(self):
        """Addresses this fake exposes (the coupler validates against these).

        Both UZF infiltration arrays are listed: FINF (routed) and SINF
        (reported). A coupler must write both."""
        n = self.name
        evt = [f'{n}/{pk}/{v}' for pk in self.EVT for v in self.EVT[pk]]
        return [f'{n}/X', f'{n}/UZF/SINF', f'{n}/UZF/FINF', f'{n}/UZF/GWD',
                f'{n}/UZF/PETMAX', f'{n}/UZF/PET', f'{n}/UZF/UZET',
                f'{n}/WEL/BOUND', f'{n}/WEL/Q'] + evt

    def get_time_step(self):
        return 1.0

    def get_value_ptr(self, addr):
        parts = addr.split('/')
        leaf = parts[-1]
        if len(parts) >= 3 and parts[-2] in self.EVT:
            return self.EVT[parts[-2]][leaf]
        return {'X': self.X, 'SINF': self.SINF, 'FINF': self.FINF,
                'GWD': self.GWD, 'BOUND': self.BOUND, 'Q': self.Q,
                'PETMAX': self.PETMAX, 'PET': self.PET,
                'UZET': self.UZET}[leaf]

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
        # uzf_ad: both arrays back to the period input (setdatafinf)
        self.FINF[:] = self.SINF_PVAR
        self.SINF[:] = self.SINF_PVAR
        self.PET[:] = self.PET_PVAR           # setdataet
        self.PETMAX[:] = self.PET_PVAR
        self._k = 0

    def solve(self, sol):
        self.PET[:] = self.PETMAX             # uzf_solve, every iteration
        self.pet_used.append(self.PET.copy())
        self._k += 1
        self.solve_calls += 1
        self.finf_iter_trace.append(self.FINF.copy())
        return self._k >= self.nouter

    def finalize_solve(self, sol):
        pass

    def finalize_time_step(self):
        self.finf_at_advance.append(self.FINF.copy())    # what was ROUTED
        self.sinf_at_advance.append(self.SINF.copy())    # what was REPORTED
        self.petmax_at_advance.append(self.PETMAX.copy())
        if self.uzet_area is not None:
            a = np.asarray(self.uzet_area, dtype=float)
            self.UZET[:] = 0.0
            self.UZET[:a.size] = -(self.PETMAX[:a.size] + self.uzet_extra) * a
        # EVT takes its curve's full rate, out of the aquifer (negative)
        area = 1.0 if self.evt_area is None else np.asarray(self.evt_area,
                                                            dtype=float)
        for arrs in self.EVT.values():
            arrs['SIMVALS'][:] = -arrs['RATE'] * area
        self.evt_at_advance.append({pk: a['RATE'].copy()
                                    for pk, a in self.EVT.items()})
        self.X += self.dh              # deterministic head evolution per SP


# --------------------------------------------------------------------- #
def _setup(nper=4, nouter=1, dh=0.0, heads0=699.0, obs_idx=None,
           obs_names=None):
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
                           inp['P_veg_zoneSP'],
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
                          obs_idx=obs_idx, obs_names=obs_names)
    api = FakeApi('toy', cMF.nlay, cMF.nrow, cMF.ncol, ctx.ncell,
                  nuzf=ctx.ncell + 3, heads0=heads0, nouter=nouter, dh_per_step=dh)
    api.evt_area = cpl.area
    return cpl, api, ctx


def test_lagged_uses_previous_sp_heads():
    dh = -0.5                      # heads drop 0.5 m per advance
    cpl, api, ctx = _setup(dh=dh, heads0=699.0)
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
    steady advance (index 0) must carry them, not the uniform perc_user: the
    mean ETg as a flat EVT curve (it was a fixed WEL rate until 2026-10-07)."""
    cpl, api, ctx = _setup()
    sp = np.full(ctx.ncell, 7.0e-4)          # per-cell mean recharge (m/d)
    se = np.full(ctx.ncell, 1.0e-4)          # per-cell mean ETg
    cpl.steady_perc, cpl.steady_etg = sp, se
    cpl.run(api)
    finf0 = api.finf_at_advance[0]           # the steady advance
    evt0 = api.evt_at_advance[0]
    assert np.allclose(finf0[:ctx.ncell], sp, atol=1e-12)
    assert np.allclose(evt0['EVT_EG'], se, atol=1e-12)
    assert np.allclose(evt0['EVT_TG'], 0.0)


def test_the_steady_evt_curve_is_flat_down_to_the_cell_bottom():
    """marmites_evt.flat_curve: the full rate from the surface down to the
    last segment above the bottom, zero at the bottom -- the mean rate the
    ETg well drew, reduced as the cell dries (AUTO_FLOW_REDUCE did that)."""
    import marmites_evt as me
    s, r, d, x, y = me.flat_curve(1e-4, 800.0, 790.0, 8)
    assert (s, r, d) == (800.0, 1e-4, 10.0)
    assert np.all(np.diff(x) > 0.0) and 0.0 < x[0] and x[-1] < 1.0
    assert np.all(y == 1.0)
    assert me.flat_curve(0.0, 800.0, 790.0, 8)[1] == 0.0
    assert me.flat_curve(1e-4, 800.0, 801.0, 8)[1] == 0.0   # no thickness


def test_steady_state_falls_back_to_perc_user():
    cpl, api, ctx = _setup()
    ctx.cMF.perc_user = 2.0e-4
    cpl.run(api)
    finf0 = api.finf_at_advance[0]
    assert np.allclose(finf0[:ctx.ncell], 2.0e-4, atol=1e-12)
    assert np.allclose(api.evt_at_advance[0]['EVT_EG'], 0.0), \
        'no mean ETg given: the steady EVT takes nothing'


def test_pointer_contents_finf_and_evt():
    cpl, api, ctx = _setup()
    res = cpl.run(api)
    # at each transient advance, FINF[0:ncell] held that SP's perc (m/d)
    for n in range(ctx.cMF.nper):
        finf = api.finf_at_advance[1 + n]
        assert np.allclose(finf[:ctx.ncell], res['perc'][n], atol=1e-12)
        assert np.all(finf[ctx.ncell:] == 0.0)
        # what EVT took (the fake takes the curves' full rate) is the ETg
        evt = api.evt_at_advance[1 + n]
        assert np.allclose(res['etg'][n], evt['EVT_EG'] + evt['EVT_TG'],
                           atol=1e-15)
    # sanity: the shallow water table (dgwt ~1 m, loam) makes ETg active,
    # so a genuinely nonzero curve crossed the EVT pointers
    assert res['etg'].max() > 0


def test_wel_is_never_bound_or_written():
    """The ETg wells went with the WEL route (2026-10-07). A WEL in the
    model is real pumping, and writing its BOUND corrupted MF6 6.7's heap
    (bisected with tests/diagnose_coupling.py): the coupler leaves it be."""
    cpl, api, ctx = _setup()
    api.BOUND[:] = -12345.0
    api.Q[:] = -54321.0
    cpl.run(api)
    assert not hasattr(cpl, 'p_q') and not hasattr(cpl, 'p_bound')
    assert np.all(api.BOUND == -12345.0) and np.all(api.Q == -54321.0)


def test_a_model_without_evt_is_refused():
    """EVT is the coupled run's only groundwater-ET path: a model built
    without the two packages must stop at the binding, naming them."""
    cpl, api = _bad('no_evt')
    with pytest.raises(coup.CouplingError, match='EVT'):
        cpl.run(api)


def test_the_coupling_modes_are_gone():
    """Lagged is the coupling (2026-10-07): neither a mode nor an
    under-relaxation is accepted any more."""
    import inspect
    sig = inspect.signature(coup.MF6Coupler.__init__)
    assert 'mode' not in sig.parameters and 'relax' not in sig.parameters
    assert not hasattr(coup.MF6Coupler, '_iterative_sp')


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
                           inp['P_veg_zoneSP'],
                           inp['Eo_zonesSP'], inp['PT_veg_zonesSP'], inp['Pe_veg_zonesSP'],
                           inp['PE_zonesSP'], inp['gridVEGarea'], inp['LAI_veg_zonesSP'],
                           inp['Zr'], inp['kTg_min'], inp['kTg_max'], inp['kT_f'], inp['kT_s'],
                           inp['NVEG'], 1000.0, 0, None, None, None, None, None,
                           None, None, None, None, None, geom=geom)
    state = mm.init_state(ctx)
    surf = [(c[1], c[2], 0) for c in cells]
    mf6b = SimpleNamespace(cMF=cMF, ncell=ctx.ncell, surf_cells=surf,
                           nrow=cMF.nrow, ncol=cMF.ncol)
    cpl = coup.MF6Coupler(mm, ctx, state, mf6b, conv_fact=1000.0)
    api = FakeApi('toy', cMF.nlay, cMF.nrow, cMF.ncol, ctx.ncell,
                  nuzf=ctx.ncell + 3, heads0=heads0)
    api.evt_area = cpl.area
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
                'no_gwd': {'GWD'}}          # SIMVALS/DRN_SEEP already unresolved

    def get_input_var_names(self):
        n = self.name
        if self._mode == 'no_finf':
            return [f'{n}/X']
        if self._mode == 'no_gwd':
            return [f'{n}/X', f'{n}/UZF/FINF']
        if self._mode == 'no_evt':
            return [v for v in super().get_input_var_names() if '/EVT_' not in v]
        return super().get_input_var_names()

    def get_value_ptr(self, addr):
        leaf = addr.rsplit('/', 1)[1]
        if leaf in self._MISSING.get(self._mode, ()):
            raise KeyError(addr)           # genuinely absent -> unresolvable
        if self._mode == 'no_evt' and '/EVT_' in addr:
            raise KeyError(addr)           # a model built without EVT
        return super().get_value_ptr(addr)


def _bad(mode, **kw):
    cpl, api, ctx, _cells, cMF = _setup_grid('dis')
    bad = _BadApi('toy', cMF.nlay, cMF.nrow, cMF.ncol, ctx.ncell,
                  nuzf=ctx.ncell + 3, mode=mode, **kw)
    bad.evt_area = cpl.area
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


def test_infiltration_is_written_to_the_routed_and_the_reported_array():
    """UZF ROUTES FINF and REPORTS SINF (MF6 6.7 source). The coupler wrote
    SINF only, after the note here said SINF was 'operative': the budget's
    INFILTRATION line followed MMsoil while the kinematic wave routed the
    steady input rate on every day of the run (2026-09-23: UZF outflows a
    constant 960.9 m3/d, the input file's 961.0; UZF budget 92 % out)."""
    cpl, api, ctx = _setup()
    res = cpl.run(api)
    assert cpl.addr_finf.endswith('/FINF'), cpl.addr_finf
    assert cpl.addr_sinf.endswith('/SINF'), cpl.addr_sinf
    for n in range(ctx.cMF.nper):
        routed = api.finf_at_advance[1 + n][:ctx.ncell]
        reported = api.sinf_at_advance[1 + n][:ctx.ncell]
        assert np.allclose(routed, res['perc'][n]), 'the wave got the input rate'
        assert np.array_equal(routed, reported), 'budget and routing disagree'
        assert not np.any(routed == -1.0), 'SINF_PVAR leaked into the routing'


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


def test_obs_capture_records_full_mm_series_at_obs_cells():
    """With obs_idx set, the coupler keeps the FULL per-SP MM flux vectors at
    those cells (for the per-point Sankey), exposes them in the result, and they
    match the aggregates: mm_obs at an obs cell must equal that cell's row of the
    per-cell arrays, and the mean over a single-cell selection equals wb_ts."""
    idx = [0, 3]
    cpl, api, ctx = _setup(obs_idx=idx, obs_names=['A', 'B'])
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
    cpl, api, _ctx = _setup()
    res = cpl.run(api)
    assert 'mm_obs' not in res and cpl.mm_obs is None


def test_all_zero_exfiltration_is_flagged(capsys):
    """A whole run with zero exfiltration everywhere must warn: it indicates
    UZF was built without SIMULATE_GWSEEP (the La Mata failure), not a result."""
    cpl, api, _ctx = _setup()
    api.GWD[:] = 0.0
    cpl.run(api)
    out = capsys.readouterr().out
    assert 'exfiltration is zero' in out.lower()
    assert 'simulate_gwseep' in out.lower()


def test_nonzero_exfiltration_is_not_flagged(capsys):
    cpl, api, _ctx = _setup()
    api.GWD[:] = 5.0
    cpl.run(api)
    assert 'exfiltration is zero' not in capsys.readouterr().out.lower()


def test_exfiltration_sign_and_scaling():
    """GWD (m3/d, + to surface) must arrive in MM as +mm/d into the soil."""
    cpl, api, ctx = _setup()
    api.GWD[:] = 25.0                       # m3/d on 100x100 m cells
    res = cpl.run(api)
    # 25 m3/d / 1e4 m2 * 1000 mm/m = 2.5 mm/d, positive into the soil
    assert np.allclose(res['exf'], 2.5)


# ------------------------------------------- saying where a long run has got
# La Mata's full record is 1949 daily stress periods, ~4.7 h, and the log
# carried no stress-period marker at all: "is it at period 20 or 1900?" could
# only be answered by measuring the size of the .hds file.

def test_a_duration_is_readable():
    from marmites_coupler import _hms
    assert _hms(45) == '45s'
    assert _hms(3 * 60 + 7) == '3m 07s'
    assert _hms(4 * 3600 + 42 * 60) == '4h 42m'
    assert _hms(-1) == '0s'


def test_progress_is_reported_but_not_every_period(capsys):
    """Every 5%, so the line cannot become the noise it exists to cut."""
    import types

    from marmites_coupler import MF6Coupler

    c = types.SimpleNamespace(_t0=0.0, _progress=None)
    c._progress = types.MethodType(MF6Coupler._progress.__func__
                                   if hasattr(MF6Coupler._progress, '__func__')
                                   else MF6Coupler._progress, c)
    import time as _t
    c._t0 = _t.time()
    for n in range(100):
        c._progress(n, 100)
    out = capsys.readouterr().out
    lines = [ln for ln in out.splitlines() if 'stress period' in ln]
    assert 5 <= len(lines) <= 21, 'reported %d time(s) in 100' % len(lines)
    assert '100/100' in out or '96/100' in out, out[-200:]


def test_progress_says_nothing_for_an_empty_run(capsys):
    import types

    from marmites_coupler import MF6Coupler

    c = types.SimpleNamespace(_t0=0.0)
    c._progress = types.MethodType(MF6Coupler._progress, c)
    c._progress(0, 0)
    assert 'stress period' not in capsys.readouterr().out


# ------------------------------------- a condition that holds is not an event

def test_a_recurring_condition_is_counted_once():
    MMsoil = _load('_mmsoil_tally',
                   os.path.join(os.path.dirname(HERE), 'MARMITESsoil',
                                'MARMITESsoil_v3.py'))

    MMsoil.report_tallies(out=False)              # start clean
    for k in range(2000):
        MMsoil.tally('Tg: soil moisture below wilting point', 0.001 * (k % 5))
    MMsoil.tally('something else')
    lines = MMsoil.report_tallies(out=False)
    assert len(lines) == 2, lines
    joined = ' '.join(lines)
    assert '2000 time(s)' in joined, joined
    assert 'out of range by 0.0000..0.0040' in joined, joined
    assert MMsoil.report_tallies(out=False) == [], 'the tally was not cleared'
