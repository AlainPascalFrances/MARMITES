# -*- coding: utf-8 -*-
"""MODFLOW 6 backend for MARMITES (Phase 3a).

Builds a MODFLOW 6 simulation with flopy.mf6 from the parsed MARMITES/MF
configuration (clsMF attributes), per the NWT->MF6 package mapping of
MARMITES_code_review.md section 4.1 and the decisions of section 6:

  * DIS structured grid (DISV arrives in Phase 4);
  * Newton formulation (`newtonoptions under_relaxation`), IMS complexity
    from the NWT options (COMPLEX for La Mata);
  * STO with an initial steady-state stress period (the dum_sssp1
    convention), daily transient SPs by default;
  * two EVT packages (evt_eg, evt_tg), one record per active surface cell:
    groundwater ET at the solved head, the curves written each day by the
    API coupler (marmites_evt). No WEL: the ETg wells went with the WEL
    route on 2026-10-07, and WEL is left for real pumping;
  * DRN / GHB from the legacy cell lists, grid-agnostic cellids;
  * UZF6 (decision: UZF always): one vertical column of UZF objects per
    active map cell, landflag on the outcrop cell, ivertcon chaining
    downward; finf driven by the API coupler. Budget terms: UZF-GWRCH
    (recharge), UZF-GWD (groundwater discharge = MARMITES exfiltration).

The class only *builds and writes* the simulation; running it requires the
mf6 binary (see tests/run_lamata_mf6.py), and the coupled run goes through
marmites_coupler.MF6Coupler using the MODFLOW 6 API (libmf6).
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import os

import numpy as np

try:
    from flopy.mf6 import (MFSimulation, ModflowGwf, ModflowGwfdis, ModflowGwfdisv,
                           ModflowGwfdrn, ModflowGwfghb, ModflowGwfic, ModflowGwfnpf,
                           ModflowGwflak, ModflowGwfmvr, ModflowGwfoc, ModflowGwfsfr,
                           ModflowGwfevt, ModflowGwfsto, ModflowGwfuzf,
                           ModflowIms, ModflowTdis, ModflowUtllaktab)
    HAS_FLOPY = True
    FLOPY_IMPORT_ERROR = None
except Exception as _exc:  # pragma: no cover - keep the true cause for diagnostics
    HAS_FLOPY = False
    FLOPY_IMPORT_ERROR = _exc


class MF6BuildError(Exception):
    """Invalid input while building the MODFLOW 6 simulation."""


# every SFR reach carries this boundname, so an observation named by it is
# MF6's sum over the whole network
SFR_BOUNDNAME = 'network'


class clsMF6:
    """Build a MODFLOW 6 simulation for the MARMITES coupled model.

    Parameters
    ----------
    cMF : ppMODFLOW_flopy_v3.clsMF
        Parsed legacy configuration (after ppMFtime and the driver's
        top/botm adjustment: top = elev - soil thickness).
    top, botm : 2-D (nrow, ncol) and 3-D (nlay, nrow, ncol) arrays [m]
        Aquifer top and layer bottoms (soil thickness already subtracted,
        as the driver does before running MMsoil).
    sim_ws : str
        Simulation workspace (created if absent).
    daily : bool
        True (default, section 4.3 decision): one transient SP per day
        (perlen=1), preceded by one steady SP. False: use the rainfall-
        aggregated cMF.perlen from ppMFtime.
    exe_name : str
        Path/name of the mf6 executable (only needed to run).
    """

    def __init__(self, cMF, top, botm, sim_ws, daily=True, exe_name='mf6',
                 grid='dis', vertices=None, cell2d=None, cell_nodes=None,
                 gwseep=True, strt_from_dem=None):
        if not HAS_FLOPY:  # pragma: no cover
            raise MF6BuildError('flopy.mf6 could not be imported. Original error:\n'
                                '  %r\n'
                                'Check `python -c "import flopy; print(flopy.__version__)"` '
                                'in this environment (needs flopy >= 3.4).' % (FLOPY_IMPORT_ERROR,))
        if grid not in ('dis', 'disv'):
            raise MF6BuildError("grid must be 'dis' or 'disv' (DISU is out of scope)")
        self.cMF = cMF
        self.sim_ws = sim_ws
        self.daily = bool(daily)
        self.exe_name = exe_name
        self.grid = grid
        # Groundwater discharge to land surface = MARMITES exfiltration (Exf_g).
        # UZF6 computes it ONLY with SIMULATE_GWSEEP; without the option the
        # GWD array stays identically zero and exfiltration is impossible.
        # "Groundwater discharge is nonzero when groundwater head is greater
        # than land surface" (mf6io); here the UZF land surface is the aquifer
        # top = elevation - soil thickness, i.e. the bottom of the MARMITES
        # soil column -- exactly the Exf_g condition of the paper (sec. 2.3).
        self.gwseep = bool(gwseep)
        # Seepage mechanism (decision 2). 'uzf' = UZF SIMULATE_GWSEEP;
        # 'drn' = a smoothed drain at the soil base, as in the CdL reference
        # model. SIMULATE_GWSEEP is deprecated (MF6 6.5.0), which recommends
        # the drain; it is smoothed too (UzfCellGroup.f90 gwseep: the same
        # -s^3 + 2s^2 ramp over SURFDEP), but with its conductance fixed to
        # area x vks / SURFDEP. CdL saw single-cell limit cycles with it.
        # With 'drn' the discharge ramps in over DDRN via cubic smoothing.
        # Unlike the reference, the flows are NOT moved to SFR: MARMITES needs
        # them returned to the soil column (Eq. 1 / sec. 2.3).
        # DEFAULT drn: a clsMF6 built without the run passing a
        # choice must not get the deprecated mechanism by accident.
        self.seep = 'drn'
        # Conductance of a seepage-face drain is a NUMERICAL device, not a
        # physical streambed/aquitard property: its job is to pin the head at
        # the land surface and carry off whatever excess arrives. It must
        # therefore be effectively free-draining. La Mata peaks at ~2060 m3/d
        # of seepage per cell, so C = 10 m2/d would need the head to stand
        # 206 m above ground to pass it -- which is exactly what the first
        # coupled run did (205 m overshoot, 1903 of 1954 cells above ground).
        # 10000 m2/d holds the excess within ~0.2 m.
        self.drn_seep_cond = 10000.0     # m2/d per cell
        # m, cubic smoothing depth -- [seep] ddrn on panel 4's DRN tab; the
        # conductance ramps from 0 at the drain (the aquifer top) to full
        # this much above it (MF6 gwf-drn.f90, get_drain_elevations)
        self.drn_seep_ddrn = 3.5
        # [seep] base: the drain base below the soil base [m], 0 = AT it;
        # [seep] cond_from: 'value' (drn_seep_cond, one number per cell) or
        # 'uzf' -- UZF's own seepage conductance, cell area x vks / SURFDEP.
        # base = SURFDEP/2, ddrn = SURFDEP and 'uzf' reproduce SIMULATE_GWSEEP
        # exactly (UzfCellGroup.f90 gwseep vs gwf-drn.f90: the same cubic).
        self.drn_seep_base = 0.0
        self.drn_seep_cond_from = 'value'
        self.drn_seep_lifted = 0
        self.drn_seep_cond_range = None
        self.drnseep_id = {}             # cell (i, j) -> DRN-SEEP boundary index
        self.ndrnseep = 0
        self.sfr_cells = set()           # (i, j) with an SFR reach: no seep drain
        # SFR. Enabled by setting sfr_pondw (the MARMITES channel-width map).
        # When on, the outlet DRN cells are replaced by the SFR outlet reaches
        # (decision 4), so the catchment discharges through the stream instead
        # of through a drain.
        self.sfr = False
        self.sfr_pondw = None            # channel width per cell [m]
        self.sfr_pondhmax = None         # channel depth per cell [m]
        self.sfr_rhk = 0.1               # streambed K [m/d]
        self.sfr_rbth = 0.5              # streambed thickness [m]
        self.sfr_man = 0.035             # Manning's n
        self.sfr_min_slope = 1e-4        # [sfr] min_slope
        self.sfr_monotonic = True        # [sfr] monotonic_bed
        self.sfr_outlet_cellids = set()  # the cellid of each outlet reach
        self.sfr_net = None              # SFRNetwork once built
        self.sfr_reach_of = {}           # (i, j) -> reach number
        # WP1d: the network comes from the MAPPED LINES. sfr_pondw is then the
        # channel PRESENCE map (1 where a line crosses), and the width and
        # incision are resolved after routing from these ParamSources.
        self.sfr_width_source = None     # marmites_config.ParamSource
        self.sfr_depth_source = None
        self.sfr_seg_of_cell = None      # dominant segment id per cell
        self.sfr_seg_params = None       # inputSTREAM_param.csv rows
        self.sfr_cell_length = None      # mapped channel length per cell [m]
        self.cell_area = None            # m2, for the drainage law
        # UZF residual water content: 'sy' derives thtr = thts - Sy (UZF1
        # without SPECIFYTHTR), 'source' takes cMF.thtr and requires it
        # consistent with Sy. Set from uzf.thtr_from by the driver.
        self.uzf_thtr_from = 'sy'
        # LAK. Enabled by setting lak_shapefile (pond polygons). Each pond
        # becomes one EMBEDDEDV lake connected through its host cell, its true
        # area carried by a stage-volume table, the stream cut out of its
        # footprint (see marmites_lak and _build_ponds).
        self.lak_shapefile = None
        self.lak_depth = None            # per-cell pond depth map [m]
        self.lak_bedleak = 1e-3          # 1/d
        self.lak_surfdep = 0.05          # m
        # groundwater ET: two EVT packages at the solved head (marmites_evt)
        self.evt_nseg = 8
        self.evt_ramp = 0.1              # m
        # LAK's own Newton loop (MAXIMUM_ITERATIONS, MAXIMUM_STAGE_CHANGE);
        # None leaves MF6's defaults, 100 and 1e-5 m. CdL's perched ponds
        # needed 200 and 1e-4 -- [lak] maxiter / stagechg on the panel.
        self.lak_maxiter = None
        self.lak_stagechg = None
        self.ponds = []
        self.lak_of_cell = {}            # footprint cell -> lake index
        self.lak_mvr = True              # route the stream through on-channel ponds
        # Adaptive time stepping. On by default: a stress period MF6 cannot
        # solve in one step is not a result, and without ATS it silently
        # becomes one (see the -132% budget discrepancy in the Phase-5 report).
        self.ats = True
        # d, the shortest retry step -- [run] ats_dtmin on the Run panel.
        # Was 1e-4 d (~9 s), which let a hopeless day be tried ~190 times.
        self.ats_dtmin = 0.01
        # multiplier on the UZF unsaturated vks, to offset the forced EPSILON
        # clamp (2.0 -> 3.5). 1.0 = no change. Calibrate so UZF-GWRCH matches
        # the target recharge.
        self.uzf_vks_scale = 1.0
        self.perioddata = None
        # THE IMS SOLVER, from [solver] on the Run panel (the driver sets
        # it; None = marmites_config.Solver's defaults, one source for both).
        # The legacy ini's NWT HEADTOL / MAXITEROUT / OPTIONS are NOT read:
        # 0.05 m left 35..65 m3/d unaccounted in every step (2026-09-23).
        self.solver = None
        self.outer_maximum = int(self._solver()['outer_maximum'])
        # UNSATURATED-ZONE ET (WP2). ALWAYS SIMULATED -- there is no switch,
        # see the perioddata block below. `uzf_et_form` is 'etwc' (a
        # water-content threshold) or 'etae' (Brooks-Corey capillary
        # pressure), which IS a choice and comes from [et].
        self.uzf_et_form = 'etwc'
        # per-cell (nlay, nrow, ncol) or None -> taken from thtr
        self.uzf_extdp = None
        self.uzf_extwc = None
        # None -> use the ini arrays; (a, b) -> head = a*elevation + b
        self.strt_from_dem = strt_from_dem
        # explicit (nlay,nrow,ncol) initial heads, e.g. a spin-up cycle's final
        # head field; overrides both of the above when set
        self.strt_array = None
        self._strt_stats = None
        self.nlay, self.nrow, self.ncol = cMF.nlay, cMF.nrow, cMF.ncol
        self.top = np.asarray(top, dtype=float)
        self.botm = np.asarray(botm, dtype=float)
        ib = np.asarray(cMF.ibound)
        if (ib < 0).any():
            # constant heads existed as ibound<0 in the legacy world; they
            # would need a CHD package here. La Mata has none.
            raise MF6BuildError('ibound<0 (constant heads) present: add a CHD '
                                'package before building (not implemented).')
        self.idomain = (np.abs(ib) > 0).astype(int)
        # active surface (outcrop) cells in MARMITES cell order
        self.outcropL = np.asarray(cMF.outcropL, dtype=int)
        self.surf_cells = [(i, j, self.outcropL[i, j] - 1)
                           for i in range(self.nrow) for j in range(self.ncol)
                           if self.outcropL[i, j] > 0]
        self.ncell = len(self.surf_cells)
        # refined (quadtree) grids: explicit icell2d per surface cell
        self._cell_nodes = None if cell_nodes is None else np.asarray(cell_nodes, dtype=int)
        # DISV geometry: one cell2d per structured cell unless supplied
        self.ncpl = self.nrow * self.ncol
        self.vertices, self.cell2d = vertices, cell2d
        if self.grid == 'disv' and (vertices is None or cell2d is None):
            import sys as _sys
            _here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if _here not in _sys.path:
                _sys.path.insert(0, _here)
            from marmites_grid import disv_from_structured
            self.vertices, self.cell2d, self.ncpl = disv_from_structured(
                cMF.delr, cMF.delc, getattr(cMF, 'xllcorner', 0.0),
                getattr(cMF, 'yllcorner', 0.0))
        elif self.grid == 'disv':
            self.ncpl = len(cell2d)
        # time discretization
        if self.daily:
            ndays = int(np.sum(cMF.perlen))
            self.perlen = np.ones(ndays, dtype=float)
        else:
            self.perlen = np.asarray(cMF.perlen, dtype=float)
        # A STEADY FIRST PERIOD only when the run has no heads of its own to
        # start from: MF6 ignores the initial heads in a steady period, so a
        # periodic spin-up -- each cycle starting where the last one ended --
        # runs without it (clsMF6.steady_first = False).
        self.steady_first = True
        self.nper = len(self.perlen) + 1          # +1 steady SP at the front
        # per (layer, row, col): the UZF initial water content carried from a
        # previous cycle (NaN = the panel's thti), set by the spin-up
        self.uzf_thti_carry = None
        # per lake: the stage it ENDED at in the previous spin-up cycle or
        # in the saved state a run starts from -- its storage, which the
        # water table alone cannot give back (None = start from the WT)
        self.lak_strt_carry = None
        self.sim = None
        self.gwf = None
        # bookkeeping filled by build()
        self.uzf_packagedata = None
        # the period row per land cell, so a test can read the
        # PET demand and the extinction depth without parsing
        # the written file
        self.uzf_perioddata = None
        self.nuzfcells = 0
        self.surfdep_check = True
        self.eps_clamped = None      # (original, clamped) when EPSILON adjusted

    # ------------------------------------------------------------------ #

    def thti_from_wc(self, wc):
        """The UZF initial water content, (nlay, nrow, ncol), from MF6's
        WCNEW per UZF object -- the mean content over each object's
        unsaturated part at the end of a run. An object with none (0: under
        the water table) is left NaN, i.e. to the panel's thti."""
        out = np.full((self.nlay, self.nrow, self.ncol), np.nan)
        wc = np.asarray(wc, dtype=float).ravel()
        for no, kij in enumerate(getattr(self, 'uzf_obj_kij', None) or []):
            if kij is not None and no < wc.size and wc[no] > 0.0:
                out[kij] = wc[no]
        return out

    @staticmethod
    def _lay3d(arr_like, nlay, nrow, ncol):
        """Normalize scalar / per-layer / full arrays to (nlay, nrow, ncol)."""
        a = np.asarray(arr_like, dtype=float)
        out = np.ones((nlay, nrow, ncol), dtype=float)
        if a.ndim == 0:
            out *= float(a)
        elif a.ndim == 1 and a.shape[0] == nlay:
            for k in range(nlay):
                out[k] *= a[k]
        elif a.ndim == 3:
            out[:] = a
        elif a.ndim == 2:
            out[:] = a[np.newaxis]
        else:
            raise MF6BuildError('cannot broadcast array of shape %s' % (a.shape,))
        return out

    def _solver(self):
        """The IMS settings: ``self.solver`` or the configuration defaults."""
        if self.solver is not None:
            return dict(self.solver)
        import dataclasses
        import marmites_config
        return dataclasses.asdict(marmites_config.Solver())

    def _uzf3d(self, name):
        """A UZF packagedata property as (nlay, nrow, ncol).

        One helper for both paths. The parameter file gave one number per
        property (``['0.05']``), the panel can give a raster per layer, and
        ``checkarray`` turns either into what ``_lay3d`` broadcasts -- so a
        uniform answer still produces exactly the array a scalar did.
        """
        raw = getattr(self.cMF, name)
        if isinstance(raw, (list, tuple)):
            raw = self.cMF.cPROCESS.checkarray(list(raw))
        return self._lay3d(np.asarray(raw, dtype=float),
                           self.nlay, self.nrow, self.ncol)

    # ---- grid-agnostic helpers (DIS vs DISV) -------------------------- #

    def _cellid(self, k, i, j):
        """MODFLOW cellid: (lay, row, col) for DIS, (lay, icell2d) for DISV.

        For a DIS-equivalent vertex grid icell2d == i*ncol + j, i.e. exactly
        the 'node' the MARMITES cell list carries. For a genuinely refined
        (quadtree) grid there is no (i, j): the caller supplies `cell_nodes`
        mapping each surface cell to its icell2d, and `i` is the cell id.
        """
        if self.grid == 'dis':
            return (int(k), int(i), int(j))
        if self._cell_nodes is not None:
            return (int(k), int(self._cell_nodes[int(i)]))
        return (int(k), int(i) * self.ncol + int(j))

    def _griddata(self, arr):
        """Reshape a (nlay, nrow, ncol) array to the grid's layout."""
        a = np.asarray(arr)
        if self.grid == 'dis':
            return a
        return a.reshape(a.shape[0], self.ncpl) if a.ndim == 3 else a.reshape(self.ncpl)

    # UZF6 hard constraints (MODFLOW 6 source: uzf module input checks)
    EPS_MIN, EPS_MAX = 3.5, 14.0

    def _validate_uzf_params(self, thtr, thts, thti, eps):
        """Enforce the UZF6 input rules, converting NWT-era values.

        Takes scalars OR arrays: the Brooks-Corey parameters may be given
        per cell, and a rule that only held for the catchment mean would
        let a single bad cell reach MODFLOW. Raises MF6BuildError for
        physically un-fixable values; clamps EPSILON to the MF6 range with
        an explicit warning, because a value outside it simply cannot be
        represented in UZF6.
        """
        import numpy as _np
        a_thtr, a_thts = _np.asarray(thtr, float), _np.asarray(thts, float)
        a_thti, a_eps = _np.asarray(thti, float), _np.asarray(eps, float)
        if not _np.all(a_thtr > 0.0):
            raise MF6BuildError(
                'UZF6 requires THTR > 0 (lowest %g). UZF1 tolerated 0 when '
                'specifythtr=0; set a residual water content on the UZF '
                'panel.' % a_thtr.min())
        if not _np.all(a_thts > a_thtr):
            raise MF6BuildError('UZF6 requires THTS > THTR everywhere '
                                '(worst pair %g, %g).'
                                % (a_thts.min(), a_thtr.max()))
        if not _np.all((a_thtr <= a_thti) & (a_thti <= a_thts)):
            raise MF6BuildError('UZF6 requires THTR <= THTI <= THTS '
                                'everywhere (THTI spans %g to %g).'
                                % (a_thti.min(), a_thti.max()))
        thtr, thts, thti, eps = a_thtr, a_thts, a_thti, a_eps
        if self.surfdep_check and not (float(self.cMF.surfdep) > 0.0):
            raise MF6BuildError('UZF6 requires SURFDEP > 0.')
        if np.any(eps < self.EPS_MIN) or np.any(eps > self.EPS_MAX):
            worst = float(eps.min() if np.any(eps < self.EPS_MIN)
                          else eps.max())
            nbad = int(np.count_nonzero((eps < self.EPS_MIN)
                                        | (eps > self.EPS_MAX)))
            # the panel refuses one number out of range; a map can still
            # carry such cells (UZF1 accepted any EPSILON)
            eps, new = worst, float(np.clip(worst, self.EPS_MIN, self.EPS_MAX))
            print('\nWARNING! UZF6 requires %.1f <= EPSILON <= %.1f but uzf.eps '
                  'gives %g (%d value(s) out of range).\n         EPSILON clamped '
                  'to %.1f there. The Brooks-Corey exponent controls '
                  'unsaturated\n         relative permeability, so drainage '
                  'through the unsaturated zone changes.\n'
                  '         Give uzf.eps values in range on the UZF panel to '
                  'control this explicitly.'
                  % (self.EPS_MIN, self.EPS_MAX, eps, nbad, new))
            self.eps_clamped = (eps, new)
            eps = np.clip(np.asarray(a_eps, float), self.EPS_MIN, self.EPS_MAX)
        return thtr, thts, thti, eps

    def save_heads_asc(self, heads, prefix):
        """Write an (nlay, nrow, ncol) head field as one ESRI-ASCII grid per
        layer: ``<prefix>_l1.asc`` .. Inactive cells get the nodata value.

        This lets a spin-up's equilibrated head field be reused directly as the
        initial condition of later runs (``strt_heads=<prefix>``), so the
        expensive spin-up is done once, not every time.
        """
        cMF = self.cMF
        h = np.asarray(heads, dtype=float).reshape(self.nlay, self.nrow, self.ncol)
        nodata = -9999.0
        idom = np.asarray(self.idomain)
        hdr = ('ncols %d\nnrows %d\nxllcorner %s\nyllcorner %s\ncellsize %s\n'
               'nodata_value %g\n' % (self.ncol, self.nrow, cMF.xllcorner,
                                      cMF.yllcorner, float(np.mean(cMF.delr)), nodata))
        paths = []
        for k in range(self.nlay):
            arr = np.where(idom[k] > 0, h[k], nodata)
            arr = np.where(np.isfinite(arr), arr, nodata)
            fn = '%s_l%d.asc' % (prefix, k + 1)
            with open(fn, 'w') as f:
                f.write(hdr)
                np.savetxt(f, arr, fmt='%.6g')
            paths.append(fn)
        return paths

    def load_heads_asc(self, prefix):
        """Read a per-layer head field written by :meth:`save_heads_asc` into an
        (nlay, nrow, ncol) array (nodata -> nan)."""
        arrs = []
        for k in range(self.nlay):
            fn = '%s_l%d.asc' % (prefix, k + 1)
            if not os.path.exists(fn):
                raise MF6BuildError('initial-heads file not found: %s' % fn)
            a = np.loadtxt(fn, skiprows=6)
            arrs.append(np.where(a <= -9990.0, np.nan, a))
        return np.stack(arrs)

    def initial_heads(self):
        """Initial heads (IC/STRT), (nlay, nrow, ncol).

        Two sources:

        * ``strt_from_dem = None`` (default): ``cMF.strt``, which
          marmites_props.initial_heads fills from layers.strt (La Mata:
          hi_topL1.asc).

        * ``strt_from_dem = (a, b)``: a linear regression on the DEM,
          ``head = a * elevation + b``, the usual first approximation for a
          sedimentary basin where the water table mimics a subdued replica of
          the topography. Starting far above the water table forces MODFLOW to
          drain a large excess through the first stress periods, which shows up
          as slow convergence (many outer iterations) and spurious infiltration
          being rejected; a regression such as (0.9995, -2.0) starts the model
          much closer to its dynamic equilibrium.

        The result is clipped to stay above each layer's bottom, otherwise the
        Newton formulation starts from a dry cell.
        """
        cMF = self.cMF
        # An explicit array (a previous spin-up cycle's final heads) takes
        # precedence: it is already an equilibrated state, so it is used as-is,
        # only clipped above each layer bottom to keep Newton off a dry cell.
        if getattr(self, 'strt_array', None) is not None:
            strt = np.asarray(self.strt_array, dtype=float).reshape(
                self.nlay, self.nrow, self.ncol).copy()
            botm = np.asarray(self.botm, dtype=float)
            strt = np.maximum(strt, botm + 0.01)
            strt = np.where(np.isfinite(strt), strt, np.asarray(cMF.strt, dtype=float))
            return strt
        if self.strt_from_dem is None:
            return np.asarray(cMF.strt, dtype=float)
        a, b = (float(x) for x in self.strt_from_dem)
        dem = np.asarray(np.ma.filled(np.asarray(cMF.elev), np.nan), dtype=float)
        if dem.ndim == 3:                      # already per layer
            dem = dem[0]
        h0 = a * dem + b
        strt = np.repeat(h0[np.newaxis, :, :], self.nlay, axis=0)
        # keep the start above each layer bottom (a small margin over botm)
        botm = np.asarray(self.botm, dtype=float)
        margin = 0.01
        strt = np.maximum(strt, botm + margin)
        # inactive cells: value is irrelevant but must be finite
        strt = np.where(np.isfinite(strt), strt, np.asarray(cMF.strt, dtype=float))
        self._strt_stats = (float(np.nanmin(h0)), float(np.nanmax(h0)))
        print('initial heads from DEM regression: head = %.4f * elevation %+.3f '
              '(range %.1f..%.1f m)' % (a, b, self._strt_stats[0], self._strt_stats[1]))
        return strt

    def _prop3d(self, name):
        """Fetch a layer property from cMF (handles list-of-arrays/scalars)."""
        val = getattr(self.cMF, name)
        if isinstance(val, list):
            layers = [np.asarray(v, dtype=float) * np.ones((self.nrow, self.ncol))
                      for v in val]
            return np.stack(layers)
        return self._lay3d(val, self.nlay, self.nrow, self.ncol)

    # ------------------------------------------------------------------ #

    def topology(self):
        """Shared-face adjacency of this model's grid, or None on DIS (WP1c.5).

        Built once and cached: SFR routing, LAK host cells and the CRR
        neighbour graph must all agree about what "adjacent" means, and
        rebuilding it per package is how they would drift apart.
        """
        if getattr(self, '_topology', 'unset') != 'unset':
            return self._topology
        self._topology = None
        if self.grid == 'disv' and self.vertices is not None:
            from marmites_topology import MeshTopology
            self._topology = MeshTopology(
                {'vertices': self.vertices, 'cell2d': self.cell2d,
                 'ncpl': self.ncpl})
        return self._topology

    def _land_surface(self):
        """The LAND SURFACE on this grid -- not the model top, which is the
        base of the MMsoil column. A channel is incised below the ground, and
        a pond dug into it."""
        z = getattr(self.cMF, 'elev', None)
        if z is None:
            return np.asarray(self.top, dtype=float)
        return np.asarray(np.ma.getdata(z), dtype=float).reshape(
            np.shape(self.top))

    def _on_catchment_edge(self):
        """``f(cell)``: does the cell have a face on the catchment's external
        boundary -- the mesh edge, or an inactive neighbour?"""
        act = np.asarray(self.idomain, dtype=int).max(axis=0) > 0
        topo = self.topology()
        if topo is not None:
            def f(c):
                i = int(c[0])
                return bool(topo.boundary[i]) or any(
                    not act[n, 0] for n in topo.neighbours[i])
            return f
        nr, nc = act.shape

        def g(c):
            i, j = int(c[0]), int(c[1])
            for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                a, b = i + di, j + dj
                if not (0 <= a < nr and 0 <= b < nc) or not act[a, b]:
                    return True
            return False
        return g

    def _segment_values(self, key, default):
        """Per stream cell, the converter's per-segment ``key`` (manning,
        rhk, rbth in inputSTREAM_param.csv) -- or ``default`` everywhere
        when there is no table or no such column.

        A panel VALUE (``sfr_param_fixed``, set by the driver for a 'value'
        producer) is the value everywhere. The converter only copies it into
        the table, and Launch re-converts when a SHAPEFILE changes, not a
        panel number -- so a streambed thickness set to 0.2 m kept reading
        the table's 0.5 m, silently (2026-10-05). Per-segment producers
        (a column, a raster) still come from the table.
        """
        fixed = (getattr(self, 'sfr_param_fixed', None) or {}).get(key)
        if fixed is not None:
            was = set()
            for r in (self.sfr_seg_params or {}).values():
                try:
                    was.add(float(r[key]))
                except (KeyError, TypeError, ValueError):
                    pass
            was = sorted(was)
            if was and was != [float(fixed)] and getattr(self, 'verbose',
                                                        True):
                print('   SFR %s: the panel\'s %g replaces the converted '
                      'table\'s %s (converted before the panel changed)'
                      % (key, float(fixed), '/'.join('%g' % v for v in was)))
            return float(fixed)
        seg, par = self.sfr_seg_of_cell, self.sfr_seg_params
        if seg is None or not par:
            return default
        out = np.full(np.shape(self.sfr_pondw), float(default))
        seg = np.asarray(seg)
        for c in zip(*np.where(np.asarray(self.sfr_pondw) > 0)):
            row = par.get(int(seg[c])) or {}
            try:
                out[c] = float(row[key])
            except (KeyError, TypeError, ValueError):
                pass
        return out

    def _build_sfr_network(self):
        """Route the mapped channel into an SFR reach network (WP3).

        ONE OUTLET, where the network leaves the catchment -- the lowest
        stream cell on its boundary -- and the reach geometry from the MAP:
        each reach as long as the channel mapped inside its cell, its slope
        over the distance to the next reach, and Manning, streambed K and
        thickness per segment as the converter resolved them.

        THE STREAM'S TOTAL DEPTH below the land surface is the soil depth of
        its cell + the channel depth + the streambed thickness (user rule,
        2026-09-27). The channel is cut THROUGH the MMsoil column into the
        aquifer: the streambed top sits one channel depth below the aquifer
        top (land surface - soil thickness), its bottom one streambed
        thickness lower. With the bed measured from the land surface instead
        (WP3), La Mata's 1.5 m of soil equalled the 1.0 m channel + 0.5 m
        streambed, so 2024 of 2663 streambed bottoms sat exactly on the
        aquifer top -- where the seepage drains hold the water table -- and
        MF6's stream-aquifer exchange, which switches at the streambed
        bottom, took 40 to 500 outer iterations a period instead of 5.
        The ROUTING still follows the land surface.
        """
        self.sfr_outlet_cells = set()
        self.sfr_outlet_cellids = set()
        if self.sfr_pondw is None:
            self.sfr = False
            return None
        from marmites_sfr import catchment_outlet, stream_network, build_sfr
        cMF = self.cMF
        land = self._land_surface()
        stream = [(int(i), int(j))
                  for i, j in zip(*np.where(np.asarray(self.sfr_pondw) > 0))]
        outlet = catchment_outlet(stream, land, self._on_catchment_edge())
        # the old rule stays the fallback: outlet DRN cells on the channel
        drn_cells = []
        if getattr(cMF, 'drn_yn', 0) == 1:
            drn_cells = [(int(i), int(j))
                         for (_l, i, j, _e, _c) in cMF.layer_row_column_elevation_cond[0]]
        topo = self.topology()
        coords = ((lambda c: topo.xy[int(c[0])]) if topo is not None else
                  (lambda c: (float(c[1]), -float(c[0]))))
        net = stream_network(self.sfr_pondw, land,
                             outlets=[outlet] if outlet is not None else None,
                             drn_cells=drn_cells, topology=topo,
                             coords=coords)
        # WP1d: width and incision are resolved AFTER routing, because the
        # drainage producer (w = a*A**b) needs contributing area and that is
        # only known once the reaches are ordered. Falling back to sfr_pondw
        # keeps the raster path working while a case still uses it.
        if self.sfr_width_source is not None:
            from marmites_channel import as_array, channel_depth, channel_width
            verbose = getattr(self, 'verbose', True)
            wid = channel_width(net, self.sfr_width_source, self.cell_area,
                                seg_of_cell=self.sfr_seg_of_cell,
                                params=self.sfr_seg_params,
                                cell_length=self.sfr_cell_length,
                                verbose=verbose)
            self.sfr_pondw = as_array(wid, np.shape(self.sfr_pondw))
            if self.sfr_depth_source is not None:
                dep = channel_depth(net, self.sfr_depth_source, verbose=verbose)
                self.sfr_pondhmax = as_array(dep, np.shape(self.sfr_pondw))
        spacing = ((lambda c, r: topo.distance(c[0], r[0]))
                   if topo is not None else None)
        # the bed datum: the AQUIFER top, i.e. land surface - soil depth
        # (the driver builds top = elev - soil thickness per cell)
        aquifer_top = np.asarray(self.top, dtype=float)
        build_sfr(net, aquifer_top, pondhmax=self.sfr_pondhmax,
                  pondw=self.sfr_pondw,
                  delr=cMF.delr, delc=cMF.delc, botm=self.botm,
                  idomain=self.idomain,
                  rbth=self._segment_values('rbth', self.sfr_rbth),
                  rhk=self._segment_values('rhk', self.sfr_rhk),
                  man=self._segment_values('manning', self.sfr_man),
                  minslope=float(self.sfr_min_slope),
                  monotonic=bool(self.sfr_monotonic),
                  reach_length=self.sfr_cell_length, spacing=spacing,
                  cellid=self._cellid, verbose=getattr(self, 'verbose', True))
        self.sfr = True
        self._bind_sfr_net(net)
        return net

    def _bind_sfr_net(self, net):
        """The sets every later step reads off the network -- rebound when
        the ponds cut their reaches out of it."""
        self.sfr_net = net
        self.sfr_cells = set(net.cells)
        self.sfr_outlet_cells = set(net.outlets)
        self.sfr_outlet_cellids = {tuple(net.packagedata[net.rno[c]][1])
                                   for c in net.outlets}
        self.sfr_reach_of = dict(net.rno)

    def _add_sfr_package(self, gwf, name):
        """ModflowGwfsfr from the routed network.

        EVAPORATION is 0 at BUILD time and written by the coupler each
        stress period from the Eo forcing (WP1d). It used to be left at zero
        permanently, because MARMITES evaporated from its own surface store;
        that store is gone, so the evaporation follows the water to the reach.
        RUNOFF is 0 here as well -- the coupler injects the MARMITES runoff of
        each stress period through the API instead.
        """
        net = self.sfr_net
        nreaches = net.nreaches
        # (rno, status, inflow, rainfall, evaporation, runoff, upstream fraction)
        spd0 = [[r, 'INFLOW', 0.0] for r in range(nreaches)]
        # WP3.4: the outlet as continuous observations, <name>.obs.sfr.csv --
        # what leaves the catchment through the stream (< 0, out of SFR), its
        # stage, and its exchange with the aquifer (> 0 when the reach LOSES
        # water to it). MF6 numbers reaches from 1 and flopy writes
        # observation ids as given.
        obs = []
        outs = [net.rno[c] + 1 for c in net.outlets]
        for k, rn in enumerate(outs):
            tag = '' if len(outs) == 1 else '_%d' % (k + 1)
            obs += [('outflow%s' % tag, 'ext-outflow', rn),
                    ('stage%s' % tag, 'stage', rn),
                    ('leakage%s' % tag, 'sfr', rn)]
        # ...and the network's own budget, which MF6 sums over every reach
        # sharing one boundname (WP3.6): the MM runoff fed in, the open-water
        # evaporation, the net exchange with the aquifer, what the ponds take
        # and give back, and what leaves. Exact at every time step and one
        # small file, where the SFR budget file costs a read per period.
        mover = bool(self.lak_mvr and self.ponds)
        bn = SFR_BOUNDNAME
        obs += [('net_inflow', 'ext-inflow', bn),
                ('net_evaporation', 'evaporation', bn),
                ('net_leakage', 'sfr', bn),
                ('net_outflow', 'ext-outflow', bn)]
        if mover:
            obs += [('net_from_mvr', 'from-mvr', bn),
                    ('net_to_mvr', 'to-mvr', bn)]
        packagedata = [list(r) + [bn] for r in net.packagedata]
        self.sfr_obs_csv = f'{name}.obs.sfr.csv'
        ModflowGwfsfr(gwf, nreaches=nreaches, boundnames=True,
                      packagedata=packagedata,
                      connectiondata=net.connectiondata,
                      perioddata={0: spd0},
                      unit_conversion=86400.0,   # Manning, SI, time unit = day
                      # MVR hands reach flow to the on-channel lakes and takes
                      # their spill back, and a package MVR names must declare
                      # MOVER itself -- otherwise MF6 stops on 'MODEL AND
                      # PACKAGE "…/SFR" DOES NOT HAVE MOVER SPECIFIED'. LAK
                      # already declared it; SFR did not (WP1d).
                      mover=mover,
                      pname='sfr', save_flows=True,
                      budget_filerecord=f'{name}.sfr.cbc',
                      stage_filerecord=f'{name}.sfr.stage',
                      observations={self.sfr_obs_csv: obs})

    def _build_ponds(self):
        """Load the pond polygons; give each a footprint, a host cell and a
        lake geometry; and cut the stream out of the footprints (WP4.1).

        The CdL design (cdl_gwf_model_fable_v2 §5b/§6): a pond owns the
        cells whose centre lies inside it, its one EMBEDDEDV connection is
        the cell holding its centroid, its rim is the mean LAND SURFACE over the
        footprint, and the stream runs THROUGH it -- the reaches inside the
        footprint are excised, the reach entering it hands its flow to the
        lake and the lake spills into the reach leaving it (MVR).

        On the grid the model is ACTUALLY on. On a projected mesh
        cMF.delr/delc/nrow/ncol describe the (ncpl, 1) proxy grid of 1 m
        squares; placing the ponds on that put every host cell somewhere
        else, and left 10 of La Mata's 11 on-channel ponds off the stream.
        """
        self.ponds = []
        if not self.lak_shapefile:
            return []
        import marmites_vector as mv
        from marmites_lak import read_pond_polygons, pond_footprints, POND_DEPTH
        grid = mv.TargetGrid.from_cMF(self.cMF)
        active = (np.asarray(self.idomain) > 0).any(axis=0)
        ponds = pond_footprints(read_pond_polygons(self.lak_shapefile), grid,
                                active=active,
                                verbose=getattr(self, 'verbose', True))
        # WP1d: lak_depth is a single value ([lak] depth) now that the pond
        # depth no longer comes from a raster; a per-cell map is still
        # accepted, so a case with a real bathymetry can supply one.
        depth = None if self.lak_depth is None else np.asarray(self.lak_depth, float)
        # THE STREAM'S RULE (user, 2026-09-27): the pond's total depth
        # below the land surface is the soil depth of its footprint + the
        # pond depth. It is dug through the MMsoil column: its bed one pond
        # depth below the AQUIFER top (land surface - soil thickness), its
        # rim -- where it spills -- the land surface (user, 2026-09-26).
        # With the bed measured from the land surface instead, La Mata's
        # 1.5 m of soil equalled the 1.5 m pond, and every pond bed sat on
        # the aquifer top (-0.44..+0.26 m), where the seepage drains hold
        # the water table -- the switch that made the streams crawl.
        land = self._land_surface()
        for p in ponds:
            i, j = p.cell
            d = POND_DEPTH
            if depth is not None:
                v = float(depth) if depth.ndim == 0 else float(depth[i, j])
                if v > 0:
                    d = v
            p.rim = float(np.mean([land[c] for c in p.cells]))
            p.bottom = float(np.mean([self.top[c] for c in p.cells])) - d
            # the lake connects at the topmost active layer, as the outcropping
            # unit is what a pond actually sits on
            k = 0
            while k < self.nlay - 1 and self.idomain[k, i, j] <= 0:
                k += 1
            p.klay = k
            p.inlet_reaches, p.outlet_reaches = [], []
            p.on_channel = any(c in self.sfr_cells for c in p.cells)
        self.ponds = ponds
        self.lak_of_cell = {c: L for L, p in enumerate(ponds) for c in p.cells}
        net = self.sfr_net
        if net is None or not self.sfr:
            return ponds
        from marmites_sfr import excise_reaches
        owner = {c: self.lak_of_cell[c] for c in net.cells
                 if c in self.lak_of_cell}
        if not owner:
            return ponds
        # The stream's passage through a pond runs from where it first
        # enters the footprint to where it LAST leaves it. A mapped line that
        # wiggles out and back in leaves a reach outside that the lake would
        # spill into and that feeds the lake again -- an MVR loop (one pond
        # on La Mata's mesh). Such detours are cut out with the pond.
        grown = True
        while grown:
            grown = False
            for c in list(owner):
                path, r = [], net.recv[c]
                while r is not None and r not in owner:
                    path.append(r)
                    r = net.recv[r]
                if path and r is not None and owner[r] == owner[c]:
                    for q in path:
                        owner[q] = owner[c]
                    grown = True
        ndetour = len(owner) - sum(1 for c in net.cells if c in self.lak_of_cell)
        into, out_of = excise_reaches(net, list(owner))
        self._bind_sfr_net(net)
        for c, r in into:
            ponds[owner[r]].inlet_reaches.append(net.rno[c])
        for r, c in out_of:
            q = ponds[owner[r]]
            if net.rno[c] not in q.outlet_reaches:
                q.outlet_reaches.append(net.rno[c])
        if getattr(self, 'verbose', True):
            print('LAK: the stream runs through %d pond(s): %d reach(es) inside '
                  'their footprints excised%s -> %d reaches; %d reach(es) hand '
                  'flow to a pond, %d take its spill'
                  % (sum(1 for p in ponds if p.on_channel), len(owner),
                     (' (%d on detours out of and back into a pond)' % ndetour
                      if ndetour else ''),
                     net.nreaches, len(into),
                     sum(len(p.outlet_reaches) for p in ponds)))
        return ponds

    def surface_fractions(self, cells_ij, area):
        """The per-cell SURFACE DESCRIPTOR of cookbook §4a (WP4.6).

        Returns ``(f_lake, f_stream)`` aligned with ``cells_ij`` -- the MM
        cell list as ordered -- and f_soil = 1 - f_lake - f_stream is where
        the MM soil column runs. Built once, from the network and the ponds
        as they will be written; the cell list itself is not touched.

        f_stream is the channel's share of its cell, the reach's surface
        ``rwid * rlen`` over the cell area (La Mata's mesh: median 0.31).
        f_lake spreads the pond's own area over its footprint, so the rain a
        pond receives is its area's and not its cells': ~1 on the mesh,
        where the footprint IS the pond, and the pond's fraction of its one
        host cell on a coarse grid, where a pond is sub-grid. Both are
        capped so that no cell is more than all open water.
        """
        n = len(cells_ij)
        area = np.asarray(area, dtype=float)
        pos = {(int(c[0]), int(c[1])): k for k, c in enumerate(cells_ij)}
        f_lake = np.zeros(n, dtype=float)
        for p in self.ponds:
            ks = [pos[c] for c in p.cells if c in pos]
            if ks:
                f_lake[ks] = min(1.0, float(p.area) / float(area[ks].sum()))
        f_stream = np.zeros(n, dtype=float)
        net = self.sfr_net if self.sfr else None
        if net is not None and net.reach_wid is not None:
            for c, r in net.rno.items():
                k = pos.get((int(c[0]), int(c[1])))
                if k is not None and area[k] > 0:
                    f_stream[k] += (float(net.reach_wid[r])
                                    * float(net.reach_len[r]) / area[k])
        f_stream = np.minimum(f_stream, 1.0 - f_lake)
        return f_lake, f_stream

    def _add_lak_package(self, gwf, name, strt_heads=None):
        """ModflowGwflak: one EMBEDDEDV lake per pond.

        No evaporation is specified at BUILD time; the coupler writes it
        each stress period from the Eo forcing (WP1d), now that MARMITES no
        longer has a surface store of its own to evaporate from.
        """
        from marmites_lak import lake_table
        pkg, conn, tables, outlets = [], [], [], []
        wt = None if strt_heads is None else np.asarray(strt_heads, dtype=float)
        carry = getattr(self, 'lak_strt_carry', None)
        if carry is not None:
            carry = np.asarray(carry, dtype=float).ravel()
            if carry.size == 0:
                # a state saved by a run without LAK: no stage to carry
                print('   LAK: the saved state has no lake stages (it was '
                      'saved without LAK); each lake starts from the water '
                      'table')
                carry = None
            elif carry.size != len(self.ponds):
                print('   LAK: %d carried stage(s) for %d lake(s) -- not this '
                      'set of ponds; each starts from the water table'
                      % (carry.size, len(self.ponds)))
                carry = None
            elif getattr(self, 'verbose', True):
                print('   LAK: %d lake(s) start at the stage they ended at '
                      '(%.2f..%.2f m above their beds)'
                      % (carry.size,
                         min(c - p.bottom for c, p in zip(carry, self.ponds)),
                         max(c - p.bottom for c, p in zip(carry, self.ponds))))
        for L, p in enumerate(self.ponds):
            i, j = p.cell
            if carry is not None and np.isfinite(carry[L]):
                # a spin-up cycle, or a run from a saved state: where the
                # lake ENDED. A perched pond holds water the water table
                # under it does not show, so starting it from the WT again
                # threw its storage away every cycle.
                p.strt = float(max(carry[L], p.bottom))
            # start each lake at its local equilibrium stage: a perched pond
            # starts nearly empty, one in contact with the water table nearly
            # full, so neither gets a shock on the first time step
            elif wt is not None:
                h = float(wt[p.klay, i, j]) if wt.ndim == 3 else float(wt[i, j])
                p.strt = float(np.clip(h, p.bottom + 0.1, p.rim))
            else:
                p.strt = p.bottom + 0.1
            pkg.append([L, p.strt, 1, 'pond%s' % p.fid])
            conn.append([L, 0, self._cellid(p.klay, i, j), 'EMBEDDEDV',
                         self.lak_bedleak, 0.0, 0.0,
                         # connlen/connwidth must be > 0; belev/telev unused
                         0.5 * max(p.depth, 0.1), float(np.sqrt(max(p.area, 1.0)))])
            rows = lake_table(p.bottom, p.rim, p.area)
            fn = f'{name}.lak{L + 1}.tab'
            ModflowUtllaktab(gwf, nrow=len(rows), ncol=4, table=rows,
                             filename=fn, pname=f'laktab_{L + 1}')
            tables.append([L, fn])
            # on-channel ponds spill back into the stream through a Manning
            # outlet at the rim; the flow is handed over by MVR
            if self.lak_mvr and p.outlet_reaches:
                outlets.append([len(outlets), L, -1, 'MANNING', p.rim,
                                float(np.sqrt(max(p.area, 1.0))), 0.035, 1e-3])
        n = len(pkg)
        newton = {}
        if getattr(self, 'lak_maxiter', None) is not None:
            newton['maximum_iterations'] = int(self.lak_maxiter)
        if getattr(self, 'lak_stagechg', None) is not None:
            newton['maximum_stage_change'] = float(self.lak_stagechg)
        ModflowGwflak(gwf, pname='lak', boundnames=True,
                      print_stage=True, save_flows=True,
                      budget_filerecord=f'{name}.lak.cbc',
                      stage_filerecord=f'{name}.lak.stage',
                      mover=bool(outlets),
                      length_conversion=1.0, time_conversion=86400.0,
                      surfdep=self.lak_surfdep, **newton,
                      nlakes=n, noutlets=len(outlets), ntables=n,
                      packagedata=pkg, connectiondata=conn,
                      # ntables=n without `tables` writes the COUNT and not the
                      # block, and MF6 stops with 'Required block "TABLES" not
                      # found. Found block "OUTLETS" instead.' It only surfaced
                      # once the pond source became readable at all (WP1d).
                      tables=tables,
                      outlets=outlets or None,
                      perioddata={0: [[L, 'RAINFALL', 0.0] for L in range(n)]})
        self.lak_outlets = outlets
        if getattr(self, 'verbose', True):
            non = sum(1 for p in self.ponds if p.on_channel)
            print('LAK: %d EMBEDDEDV lake(s), %d on-channel, %d outlet(s); '
                  'bedleak %.3g 1/d (evaporation from the Eo forcing)'
                  % (n, non, len(outlets), self.lak_bedleak))
            print('     surfdep %g m; LAK Newton: %s iterations, stage change '
                  '%s m' % (self.lak_surfdep,
                            newton.get('maximum_iterations', '100 (MF6)'),
                            newton.get('maximum_stage_change', '1e-5 (MF6)')))

    def _add_mvr_package(self, gwf, name):
        """Route the stream through the on-channel ponds.

        Every reach that drained into a pond's footprint hands all its flow
        to the lake; the lake's outlet spills into the reach(es) leaving the
        footprint, shared equally (CdL). The reaches inside the footprint
        are gone, so without this the stream would stop at the pond.
        """
        recs = []
        for L, p in enumerate(self.ponds):
            for r in p.inlet_reaches:
                recs.append(['sfr', r, 'lak', L, 'FACTOR', 1.0])
        for k, out in enumerate(getattr(self, 'lak_outlets', [])):
            p = self.ponds[out[1]]
            for r in p.outlet_reaches:
                recs.append(['lak', k, 'sfr', r, 'FACTOR',
                             1.0 / len(p.outlet_reaches)])
        if not recs:
            return
        ModflowGwfmvr(gwf, maxmvr=len(recs), maxpackages=2,
                      packages=[['sfr'], ['lak']],
                      perioddata={0: recs}, pname='mvr',
                      print_flows=False, budget_filerecord=f'{name}.mvr.cbc')
        if getattr(self, 'verbose', True):
            nin = sum(1 for r in recs if r[0] == 'sfr')
            print('MVR: %d stream->pond hand-off(s), %d pond->stream spill(s)'
                  % (nin, len(recs) - nin))

    def _uzf_vks_grid(self, k33):
        """UZF's vks per (layer, row, col): the layer's k33 (iuzfopt 2) or
        the vks raster (iuzfopt 1), times uzf_vks_scale -- the one rule the
        UZF package and the UZF-equivalent seepage drain both use."""
        cMF = self.cMF
        scale = float(getattr(self, 'uzf_vks_scale', 1.0))
        if int(getattr(cMF, 'iuzfopt', 2)) != 1:
            return scale * np.asarray(k33, dtype=float)
        return scale * self._lay3d(np.asarray(cMF.vks_actual, dtype=float),
                                   self.nlay, self.nrow, self.ncol)

    def _cell_area_grid(self):
        """Cell areas (nrow, ncol) [m2], as MODFLOW 6 computes them: delr x
        delc, or a DISV polygon's shoelace relative to its first vertex."""
        if self.grid == 'dis':
            return np.outer(np.asarray(self.cMF.delc, dtype=float),
                            np.asarray(self.cMF.delr, dtype=float))
        from marmites_grid import polygon_area
        vxy = {int(v[0]): (float(v[1]), float(v[2])) for v in self.vertices}
        out = np.zeros((self.nrow, self.ncol))
        for r in self.cell2d:
            ic = int(r[0])
            out[ic, 0] = polygon_area([vxy[int(iv)]
                                       for iv in r[4:4 + int(r[3])]])
        return out

    def build(self):
        cMF = self.cMF
        os.makedirs(self.sim_ws, exist_ok=True)
        name = cMF.modelname.lower()

        sim = MFSimulation(sim_name=name, sim_ws=self.sim_ws, exe_name=self.exe_name,
                           version='mf6')
        # FULL DOUBLE PRECISION in every file MF6 reads. flopy's default is
        # an 8-digit mantissa, so a vertex at y = 4553208.97 m was written to
        # the centimetre, and MF6's cell areas came out up to 1.9 % off the
        # mesh MMsoil and the coupler use -- 0.0127 m2 is La Mata's smallest
        # cell, and 664 cells were off by more than 0.1 %. Every m3 <-> mm
        # conversion between the two sides was off by as much there; UZF's
        # ET read back on the coupler's area exceeded PET by up to 0.011 mm/d
        # in 207,595 cell-periods (2026-09-24).
        sim.simulation_data.float_precision = 16
        sim.simulation_data.float_characters = 24
        # TDIS: steady SP first (unless the run starts from known heads),
        # then transient
        s0 = 1 if self.steady_first else 0
        self.nper = len(self.perlen) + s0
        perioddata = ([(1.0, 1, 1.0)] if s0 else []) +             [(float(p), 1, 1.0) for p in self.perlen]
        # ATS: let MF6 subdivide any stress period it cannot solve in one step.
        # Without it the first transient day -- a step change from the steady
        # state onto a free-draining seepage boundary -- failed to converge and
        # discharged 9.6e6 m3 in a single step, 8x the recharge of the entire
        # run, leaving a cumulative budget discrepancy of -132%.
        ats = None
        if self.ats:
            # dtmin no longer than the period: MF6 wants dt0 >= dtmin, and
            # a floor at the period length already means "no retry"
            recs = [(i, float(p), min(float(self.ats_dtmin), float(p)),
                     float(p), 2.0, 5.0)
                    for i, (p, _, _) in enumerate(perioddata)
                    if i >= s0]                    # not the steady-state period
            ats = {'maxats': len(recs), 'perioddata': recs}
        ModflowTdis(sim, time_units='DAYS', nper=self.nper, perioddata=perioddata,
                    ats_perioddata=ats)
        self.perioddata = perioddata

        # IMS from [solver] (Run panel). Newton + DBD + BICGSTAB are fixed.
        sv = self._solver()
        self.outer_maximum = int(sv['outer_maximum'])
        ims = ModflowIms(sim, print_option='SUMMARY',
                         complexity=str(sv['complexity']).upper(),
                         outer_dvclose=float(sv['outer_dvclose']),
                         outer_maximum=self.outer_maximum,
                         inner_dvclose=float(sv['inner_dvclose']),
                         # flopy takes inner_rclose as a record: a bare
                         # float (a list is refused) writes it with no
                         # rclose_option = MF6's per-cell infinity norm
                         rcloserecord=float(sv['inner_rclose']),
                         under_relaxation='DBD', linear_acceleration='BICGSTAB')

        gwf = ModflowGwf(sim, modelname=name, newtonoptions='UNDER_RELAXATION',
                         save_flows=True)
        sim.register_ims_package(ims, [name])

        if self.grid == 'dis':
            ModflowGwfdis(gwf, nlay=self.nlay, nrow=self.nrow, ncol=self.ncol,
                          delr=cMF.delr, delc=cMF.delc,
                          top=self.top, botm=self.botm, idomain=self.idomain,
                          length_units='METERS',
                          xorigin=float(cMF.xllcorner), yorigin=float(cMF.yllcorner))
        else:
            ModflowGwfdisv(gwf, nlay=self.nlay, ncpl=self.ncpl,
                           nvert=len(self.vertices),
                           vertices=self.vertices, cell2d=self.cell2d,
                           top=self._griddata(self.top),
                           botm=self._griddata(self.botm),
                           idomain=self._griddata(self.idomain),
                           length_units='METERS')

        strt = self.initial_heads()
        ModflowGwfic(gwf, strt=self._griddata(strt))

        # NPF: icelltype from laytyp; k from hk; k33 from vka
        # legacy LAYVKA != 0 -> VKA is the ratio hk/vk  =>  k33 = hk / vka
        hk = self._prop3d('hk_actual')
        vka = self._prop3d('vka_actual')
        layvka = np.asarray(cMF.layvka, dtype=int)
        k33 = np.empty_like(hk)
        for L in range(self.nlay):
            if layvka[L] != 0:
                with np.errstate(divide='ignore', invalid='ignore'):
                    k33[L] = np.where(vka[L] > 0, hk[L] / vka[L], hk[L])
            else:
                k33[L] = vka[L]
        # [solver] cell_averaging (Run panel): MF6's harmonic default, or
        # an ALTERNATIVE_CELL_AVERAGING -- amt-hmk keeps the conductance up
        # as a cell dewaters, as the CdL model runs
        _avg = str(sv.get('cell_averaging', 'harmonic') or 'harmonic').lower()
        ModflowGwfnpf(gwf, icelltype=list(np.asarray(cMF.laytyp, dtype=int)),
                      k=self._griddata(hk), k33=self._griddata(k33),
                      save_specific_discharge=False,
                      alternative_cell_averaging=(None if _avg == 'harmonic'
                                                  else _avg))

        # STO: steady first SP, transient afterwards -- or transient from
        # the start when the run starts from known heads
        ss = self._prop3d('ss_actual')
        sy = self._prop3d('sy_actual')
        _sto = ({'steady_state': {0: True}, 'transient': {1: True}} if s0
                else {'transient': {0: True}})
        ModflowGwfsto(gwf, iconvert=list(np.asarray(cMF.laytyp, dtype=int)),
                      ss=self._griddata(ss), sy=self._griddata(sy), **_sto)

        # GROUNDWATER ET (2026-10-05; the only route since 2026-10-07): taken
        # at the head MF6 solves for, by two EVT packages -- Eg and Tg apart,
        # so the budget keeps them apart -- one record per land column in
        # the cell order, the coupler writing each day's curve (marmites_evt)
        # and, for a steady first period, a flat curve at the mean ETg.
        # Written inert (rate 0), so a standalone run of the files takes no
        # groundwater ET. The ETg wells this replaced (one WEL per surface
        # cell, AUTO_FLOW_REDUCE) are gone: WEL is left for real pumping.
        self.evt_packages = []
        nseg = int(self.evt_nseg)
        x = [float(v) for v in np.linspace(0.0, 1.0, nseg + 1)[1:-1]]
        for pname in ('evt_eg', 'evt_tg'):
            rows = [[self._cellid(k, i, j), float(self.top[i, j]), 0.0,
                     1.0] + x + [0.0] * (nseg - 1)
                    for (i, j, k) in self.surf_cells]
            ModflowGwfevt(gwf, pname=pname, nseg=nseg,
                          maxbound=len(rows), save_flows=True,
                          stress_period_data={0: rows})
            self.evt_packages.append(pname)

        # SFR and the ponds are resolved first: the outlet reaches replace the
        # outlet DRN cells, and the stream cells are excluded from the seepage
        # drains.
        self._build_sfr_network()
        self._build_ponds()

        # DRN from the legacy list [(l, i, j, elev, cond), ...]
        if getattr(cMF, 'drn_yn', 0) == 1:
            drn_spd = [[self._cellid(l, i, j), float(e), float(c)]
                       for (l, i, j, e, c) in cMF.layer_row_column_elevation_cond[0]
                       # decision 4: the catchment now discharges through the
                       # SFR outlet reach, so keeping the outlet drain as well
                       # would give the water two exits -- but ONLY that one
                       # record, the reach's own cell and layer (WP3.3): the
                       # other outlet drains, and the one under the reach in
                       # the layer below, are legitimate boundary drainage
                       if not (self.sfr and tuple(self._cellid(l, i, j))
                               in self.sfr_outlet_cellids)]
            if drn_spd:
                ModflowGwfdrn(gwf, stress_period_data={0: drn_spd}, pname='drn',
                              maxbound=len(drn_spd), save_flows=True)

        # DRN-SEEP: groundwater seepage at the soil base as a smoothed drain
        # (decision 2/3), MF6's recommended replacement for UZF
        # SIMULATE_GWSEEP (deprecated 6.5.0). Both ramp in with the same
        # cubic; the drain sets its conductance and ramp on its own. It sits
        # at the MF6 top (the base of the MARMITES soil column) with
        # AUXDEPTHNAME, so MF6 scales its conductance from 0 there to full
        # DDRN above it.
        #
        # These flows are NOT routed away with MVR (unlike the CdL reference
        # model, which sends them to SFR). MARMITES needs them back in the soil
        # column: the coupler reads the drain SIMVALS and feeds them to the
        # bottom soil layer as exfiltration, where the upward cascade of
        # Eq. 1/1b can turn the excess into Dunnian runoff.
        self.drnseep_id = {}
        if self.seep == 'drn':
            base = float(getattr(self, 'drn_seep_base', 0.0) or 0.0)
            by_uzf = str(getattr(self, 'drn_seep_cond_from', 'value')) == 'uzf'
            if by_uzf:
                # UZF's own seepage conductance (gwseep: Q = area x vks), over
                # its SURFDEP -- per cell, from the numbers UZF gets
                vks_grid = self._uzf_vks_grid(k33)
                area_grid = self._cell_area_grid()
                surfdep_u = float(np.ravel(np.asarray(cMF.surfdep,
                                                      dtype=float))[0])
            seep_spd, lifted, conds = [], 0, []
            # ... and in every POND cell, the stream cells' exemption (user,
            # 2026-10-04): a drain at the soil base there (100 m2/d) would
            # carry the groundwater past the clay bed (1e-3 /d) and into the
            # pond a day later through MMsoil, so bedleak -- the pond's
            # calibration lever -- would no longer govern what the pond and
            # the aquifer exchange. The whole footprint, not only the host
            # cell that holds the lake's connection.
            lake_cells = {tuple(int(v) for v in c)
                          for c in (self.lak_of_cell or {})}
            self.drn_seep_lake_skipped = 0
            for n, (i, j, k) in enumerate(self.surf_cells):
                if (i, j) in self.sfr_cells:      # SFR handles seepage there
                    continue
                if (int(i), int(j)) in lake_cells:  # ... and LAK here
                    self.drn_seep_lake_skipped += 1
                    continue
                elev = float(self.top[i, j]) - base
                floor = float(np.asarray(self.botm)[k, i, j]) + 0.01
                if elev < floor:
                    # MF6 refuses a drain below its cell's bottom
                    elev, lifted = floor, lifted + 1
                cond = (float(area_grid[i, j]) * float(vks_grid[k, i, j])
                        / surfdep_u if by_uzf else float(self.drn_seep_cond))
                conds.append(cond)
                self.drnseep_id[(i, j)] = len(seep_spd)
                seep_spd.append([self._cellid(k, i, j), elev, cond,
                                 float(self.drn_seep_ddrn)])
            self.drn_seep_lifted = lifted
            self.drn_seep_cond_range = ((min(conds), max(conds))
                                        if conds else None)
            if seep_spd:
                ModflowGwfdrn(gwf, stress_period_data={0: seep_spd},
                              auxiliary=['ddrn'], auxdepthname='ddrn',
                              pname='drn_seep', maxbound=len(seep_spd),
                              save_flows=True)
            self.ndrnseep = len(seep_spd)

        # GHB from the legacy list [(l, i, j, head, cond), ...]
        if getattr(cMF, 'ghb_yn', 0) == 1:
            ghb_spd = [[self._cellid(l, i, j), float(h), float(c)]
                       for (l, i, j, h, c) in cMF.layer_row_column_head_cond[0]]
            ModflowGwfghb(gwf, stress_period_data={0: ghb_spd}, pname='ghb',
                          maxbound=len(ghb_spd), save_flows=True)

        # UZF6: vertical column of UZF objects per active map cell.
        # Object order: the first ncell objects are the land-surface cells in
        # MARMITES cell order (so GWD/FINF slices [0:ncell] map 1:1 to cells).
        # --- UZF6 water-content / Brooks-Corey parameters -----------------
        # NWT->MF6 semantic differences (both hit La Mata, see Phase-4 report):
        #  * UZF1 allowed THTR = 0 when specifythtr = 0; UZF6 *requires*
        #    THTR > 0, so the residual water content given in the ini is
        #    always used (the ini supplies it even when specifythtr = 0).
        #  * UZF1 accepted any EPSILON; UZF6 enforces 3.5 <= EPSILON <= 14.0.
        #  * they are per CELL now: the panel takes a raster, a polygon
        #    attribute or one number for each, and a uniform answer still
        #    broadcasts to exactly the array a single value produced.
        thtr = self._uzf3d('thtr')
        thts = self._uzf3d('thts')
        thti = self._uzf3d('thti')
        eps = self._uzf3d('eps')
        # CONSISTENT WITH STO's Sy -- the very array written to STO above.
        # UZF1 without SPECIFYTHTR derived thtr = thts - Sy; UZF6 has no such
        # option and leaves it to us (marmites_props.uzf_water_contents).
        import marmites_props as _props
        try:
            thtr, thti, _notes = _props.uzf_water_contents(
                thtr, thts, thti, self._prop3d('sy_actual'),
                thtr_from=str(getattr(self, 'uzf_thtr_from', 'sy')),
                active=np.asarray(self.idomain) > 0)
        except _props.PropertyError as exc:
            raise MF6BuildError(str(exc))
        for _n in _notes:
            print(_n)
        # THE WATER THE UNSATURATED ZONE HELD at the end of the previous
        # spin-up cycle: MF6's WCNEW, the mean content over each object's
        # unsaturated part. Restarting at the panel's thti (raised to thtr:
        # empty) would throw that water away between cycles.
        if self.uzf_thti_carry is not None:
            _c = np.asarray(self.uzf_thti_carry, dtype=float)
            _c = _c.reshape(thti.shape)
            _ok = np.isfinite(_c)
            thti = np.where(_ok, np.clip(_c, thtr, thts), thti)
            print('UZF thti carried from the previous cycle in %d cell-layer(s)'
                  % int(_ok.sum()))
        surfdep = float(np.ravel(np.asarray(cMF.surfdep, dtype=float))[0])
        thtr, thts, thti, eps = self._validate_uzf_params(thtr, thts, thti, eps)
        # vertical K for UZF: iuzfopt==1 -> vks array; iuzfopt==2 -> layer k33
        # (the rule is _uzf_vks_grid, shared with the seepage drain)
        # UZF unsaturated conductivity scale. UZF6 forbids EPSILON < 3.5, but the
        # NWT model used 2.0; a higher Brooks-Corey exponent lowers the
        # unsaturated relative permeability K(theta)=VKS*Se^eps, throttling
        # recharge through a deep unsaturated zone (La Mata drained to ~18 m).
        # Raising VKS restores the recharge rate the NWT model had at eps=2.0.
        # This scales ONLY the UZF vks column, never the aquifer NPF k33.
        vks_scale = float(getattr(self, 'uzf_vks_scale', 1.0))
        if vks_scale != 1.0:
            print('UZF vks scaled x%.3g to offset the EPSILON 2.0->3.5 clamp '
                  '(unsaturated recharge throttle)' % vks_scale)
        # the vks per object, as one grid -- shared with the seepage drain
        vks_grid_u = self._uzf_vks_grid(k33)
        # build columns
        pkdata = []          # (iuzno, cellid, landflag, ivertcon, surfdep, vks, thtr, thts, thti, eps)
        subs = []            # subsurface objects appended after the land cells
        iuzno_land = {}
        # first pass: land cells occupy iuzno 0..ncell-1
        for n, (i, j, k) in enumerate(self.surf_cells):
            iuzno_land[(i, j)] = n
        next_no = self.ncell
        col_children = {}    # land iuzno -> list of (iuzno, cellid, k)
        for n, (i, j, k) in enumerate(self.surf_cells):
            chain = []
            for kk in range(k + 1, self.nlay):
                if self.idomain[kk, i, j] > 0:
                    chain.append((next_no, kk))
                    next_no += 1
            col_children[n] = chain
        self.nuzfcells = next_no
        # the UZF objects of each land column, land object first (0-based):
        # the coupler sums UZF's actual ET over them (WP2)
        self.uzf_columns = [[n] + [no for no, _kk in col_children[n]]
                            for n in range(self.ncell)]
        # (layer, row, col) of every UZF object, in object order: what maps
        # MF6's per-object arrays (WCNEW) back onto the grid
        self.uzf_obj_kij = [None] * self.nuzfcells
        for n, (i, j, k) in enumerate(self.surf_cells):
            self.uzf_obj_kij[n] = (k, i, j)
            for no, kk in col_children[n]:
                self.uzf_obj_kij[no] = (kk, i, j)
        for n, (i, j, k) in enumerate(self.surf_cells):
            chain = col_children[n]
            ivertcon = chain[0][0] if chain else -1
            vks = float(vks_grid_u[k, i, j])
            pkdata.append((n, self._cellid(k, i, j), 1, ivertcon, surfdep,
                           vks, float(thtr[k, i, j]), float(thts[k, i, j]),
                           float(thti[k, i, j]), float(eps[k, i, j])))
            for ci, (no, kk) in enumerate(chain):
                child_ivert = chain[ci + 1][0] if ci + 1 < len(chain) else -1
                vks_c = float(vks_grid_u[kk, i, j])
                subs.append((no, self._cellid(kk, i, j), 0, child_ivert,
                             surfdep, vks_c, float(thtr[kk, i, j]),
                             float(thts[kk, i, j]), float(thti[kk, i, j]),
                             float(eps[kk, i, j])))
        pkdata = pkdata + subs
        # PERIOD DATA. (iuzno, finf, pet, extdp, extwc, ha, hroot, rootact)
        # on the land cells only.
        #
        # THE ET SPLIT (WP2, the cookbook's ruling). Total ET is ETsoil +
        # ETuzf + ETg: the soil column and the groundwater are MARMITES's,
        # the UNSATURATED ZONE between them is UZF's. MODFLOW 6 supports
        # exactly that division -- "et can be simulated in the uzf cell and
        # not the gwf cell by omitting keywords linear_gwet and square_gwet"
        # (mf6io) -- so simulate_et goes on WITHOUT either gwet keyword and
        # groundwater ET stays with MARMITES, which applies it through EVT.
        #
        # `pet` starts at 0 and the COUPLER writes the demand each step: it
        # is a daily quantity MARMITES computes, not a property of the
        # model. What comes back is ETuzf ACTUAL, read from the UZF budget,
        # and the residual of the demand chain goes to ETg -- never the
        # demand itself, which UZF may not have been able to meet.
        _extdp = (self._lay3d(np.asarray(self.uzf_extdp, dtype=float),
                              self.nlay, self.nrow, self.ncol)
                  if self.uzf_extdp is not None else None)
        _extwc = (self._lay3d(np.asarray(self.uzf_extwc, dtype=float),
                              self.nlay, self.nrow, self.ncol)
                  if self.uzf_extwc is not None else None)
        pdata0 = []
        for n, (i, j, k) in enumerate(self.surf_cells):
            dp = float(_extdp[k, i, j]) if _extdp is not None else 0.0
            # extwc must lie between thtr and thts; thtr is the floor
            # UZF6 already enforces, so it is the honest default.
            wc = (float(_extwc[k, i, j]) if _extwc is not None
                  else float(thtr[k, i, j]))
            pdata0.append((n, float(getattr(cMF, 'perc_user', 0.0)), 0.0,
                           dp, wc, 0.0, 0.0, 0.0))
        # THE OBJECTS BELOW carry the column's depth too. UZF measures extdp
        # from the land surface and hands the unmet PET down (setbelowpet),
        # but each object computes its own ET zone from ITS OWN extdp
        # (setdataet) -- the land row's is not passed on. Without these rows
        # a 15 m holm-oak depth would stop at the bottom of layer 1. They
        # come after every land row, so their own extwc is the one that
        # stands (setdataetwc also copies the land row's down).
        for n, (i, j, k) in enumerate(self.surf_cells):
            dp = float(_extdp[k, i, j]) if _extdp is not None else 0.0
            for no, kk in col_children[n]:
                wc = (float(_extwc[kk, i, j]) if _extwc is not None
                      else float(thtr[kk, i, j]))
                pdata0.append((no, 0.0, 0.0, dp, wc, 0.0, 0.0, 0.0))
        # SIMULATE_ET IS ALWAYS ON (WP2). Total ET has three sources and the
        # deep unsaturated zone is one of them, so a switch for it could only
        # ever be left in the position that evaporates nothing from the deep
        # zone -- which is the behaviour WP2 exists to end.
        _etkw = {'simulate_et': True}
        # ... and NEITHER linear_gwet NOR square_gwet: groundwater ET is
        # MARMITES's, and asking MODFLOW for it as well would remove the same
        # water twice.
        if self.uzf_et_form == 'etae':
            _etkw['unsat_etae'] = True
        else:
            _etkw['unsat_etwc'] = True
        ModflowGwfuzf(gwf, nuzfcells=self.nuzfcells, ntrailwaves=int(getattr(cMF, 'ntrail2', 7)),
                      nwavesets=int(getattr(cMF, 'nsets', 40)),
                      packagedata=pkdata, perioddata={0: pdata0},
                      # exactly one seepage mechanism, never both, or the
                      # discharge would be counted twice
                      simulate_gwseep=(self.gwseep and self.seep == 'uzf'),
                      pname='uzf', save_flows=True,
                      budget_filerecord=f'{name}.uzf.cbc', **_etkw)

        if self.sfr:
            self._add_sfr_package(gwf, name)
        if self.ponds:
            self._add_lak_package(gwf, name, strt_heads=strt)
            self._add_mvr_package(gwf, name)

        ModflowGwfoc(gwf, head_filerecord=f'{name}.hds',
                     budget_filerecord=f'{name}.cbc',
                     saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])

        self.sim, self.gwf = sim, gwf
        self.uzf_packagedata = pkdata
        self.uzf_perioddata = pdata0
        return sim

    # ------------------------------------------------------------------ #

    def write(self):
        if self.sim is None:
            self.build()
        self.sim.write_simulation(silent=True)

    def run(self):  # pragma: no cover (needs mf6 binary)
        ok, buff = self.sim.run_simulation(silent=False)
        if not ok:
            raise MF6BuildError('MODFLOW 6 run failed; check %s' % self.sim_ws)
        return ok


if __name__ == '__main__':
    print('Build MF6 simulations through tests/run_lamata_mf6.py or the coupler.')
