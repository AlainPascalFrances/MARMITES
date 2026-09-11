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
  * WEL with AUTO_FLOW_REDUCE in every active surface cell (ETg sink,
    driven by the API coupler);
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
                           ModflowGwfsto, ModflowGwfuzf, ModflowGwfwel,
                           ModflowIms, ModflowTdis, ModflowUtllaktab)
    HAS_FLOPY = True
    FLOPY_IMPORT_ERROR = None
except Exception as _exc:  # pragma: no cover - keep the true cause for diagnostics
    HAS_FLOPY = False
    FLOPY_IMPORT_ERROR = _exc


class MF6BuildError(Exception):
    """Invalid input while building the MODFLOW 6 simulation."""


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
        # 'drn' = a smoothed land-surface drain, as in the CdL reference model,
        # which replaced SIMULATE_GWSEEP because the latter is deprecated and
        # switches discharge on/off discontinuously (single-cell limit cycles).
        # With 'drn' the discharge ramps in over DDRN via cubic smoothing.
        # Unlike the reference, the flows are NOT moved to SFR: MARMITES needs
        # them returned to the soil column (Eq. 1 / sec. 2.3).
        self.seep = 'uzf'
        # Conductance of a seepage-face drain is a NUMERICAL device, not a
        # physical streambed/aquitard property: its job is to pin the head at
        # the land surface and carry off whatever excess arrives. It must
        # therefore be effectively free-draining. La Mata peaks at ~2060 m3/d
        # of seepage per cell, so C = 10 m2/d would need the head to stand
        # 206 m above ground to pass it -- which is exactly what the first
        # coupled run did (205 m overshoot, 1903 of 1954 cells above ground).
        # 10000 m2/d holds the excess within ~0.2 m.
        self.drn_seep_cond = 10000.0     # m2/d per cell
        self.drn_seep_ddrn = 3.5         # m, cubic smoothing depth
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
        # LAK. Enabled by setting lak_shapefile (pond polygons). Every La Mata
        # pond is smaller than a 50 m cell, so each becomes one EMBEDDEDV lake
        # inside its host cell with its true area carried by a stage-volume
        # table (see marmites_lak).
        self.lak_shapefile = None
        self.lak_depth = None            # per-cell pond depth map [m]
        self.lak_bedleak = 1e-3          # 1/d
        self.lak_surfdep = 0.05          # m
        self.ponds = []
        self.lak_mvr = True              # route the stream through on-channel ponds
        # Adaptive time stepping. On by default: a stress period MF6 cannot
        # solve in one step is not a result, and without ATS it silently
        # becomes one (see the -132% budget discrepancy in the Phase-5 report).
        self.ats = True
        self.ats_dtmin = 1e-4            # d (~9 s) floor on the sub-step
        # multiplier on the UZF unsaturated vks, to offset the forced EPSILON
        # clamp (2.0 -> 3.5). 1.0 = no change. Calibrate so UZF-GWRCH matches
        # the target recharge.
        self.uzf_vks_scale = 1.0
        self.perioddata = None
        self.outer_maximum = min(int(getattr(cMF, 'maxiterout', 500)), 500)
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
        self.nper = len(self.perlen) + 1          # +1 steady SP at the front
        self.sim = None
        self.gwf = None
        # bookkeeping filled by build()
        self.uzf_packagedata = None
        self.nuzfcells = 0
        self.surfdep_check = True
        self.eps_clamped = None      # (original, clamped) when EPSILON adjusted

    # ------------------------------------------------------------------ #

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

        Raises MF6BuildError for physically un-fixable values; clamps
        EPSILON to the MF6 range with an explicit warning, because a value
        outside it simply cannot be represented in UZF6.
        """
        if not (thtr > 0.0):
            raise MF6BuildError(
                'UZF6 requires THTR > 0 (got %g). UZF1 tolerated 0 when '
                'specifythtr=0; set a residual water content in the MF ini '
                '(the "specifythtr thtr" line).' % thtr)
        if not (thts > thtr):
            raise MF6BuildError('UZF6 requires THTS (%g) > THTR (%g).' % (thts, thtr))
        if not (thtr <= thti <= thts):
            raise MF6BuildError('UZF6 requires THTR <= THTI <= THTS '
                                '(got %g, %g, %g).' % (thtr, thti, thts))
        if self.surfdep_check and not (float(self.cMF.surfdep) > 0.0):
            raise MF6BuildError('UZF6 requires SURFDEP > 0.')
        if eps < self.EPS_MIN or eps > self.EPS_MAX:
            new = min(max(eps, self.EPS_MIN), self.EPS_MAX)
            print('\nWARNING! UZF6 requires %.1f <= EPSILON <= %.1f but the MF '
                  'ini specifies %g.\n         EPSILON clamped to %.1f. This is a real '
                  'NWT->MF6 difference:\n         the Brooks-Corey exponent controls '
                  'unsaturated relative permeability,\n         so drainage through the '
                  'unsaturated zone will differ from the NWT model.\n'
                  '         Set a value in range in the MF ini to control this explicitly.'
                  % (self.EPS_MIN, self.EPS_MAX, eps, new))
            self.eps_clamped = (eps, new)
            eps = new
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

        * ``strt_from_dem = None`` (default): the arrays named in the MF ini
          (La Mata: hi_topL1.asc), i.e. the legacy behaviour.

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

    def _build_sfr_network(self):
        """Route the MARMITES channel map (PONDw) into an SFR reach network."""
        self.sfr_outlet_cells = set()
        if self.sfr_pondw is None:
            self.sfr = False
            return None
        from marmites_sfr import stream_network, build_sfr
        cMF = self.cMF
        # outlets: the existing outlet DRN cells that lie on the channel
        drn_cells = []
        if getattr(cMF, 'drn_yn', 0) == 1:
            drn_cells = [(int(i), int(j))
                         for (_l, i, j, _e, _c) in cMF.layer_row_column_elevation_cond[0]]
        net = stream_network(self.sfr_pondw, self.top, drn_cells=drn_cells,
                             topology=self.topology())
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
        build_sfr(net, self.top, pondhmax=self.sfr_pondhmax, pondw=self.sfr_pondw,
                  delr=cMF.delr, delc=cMF.delc, botm=self.botm,
                  idomain=self.idomain, rbth=self.sfr_rbth, rhk=self.sfr_rhk,
                  man=self.sfr_man, cellid=self._cellid,
                  verbose=getattr(self, 'verbose', True))
        self.sfr = True
        self.sfr_net = net
        self.sfr_cells = set(net.cells)
        self.sfr_outlet_cells = set(net.outlets)
        self.sfr_reach_of = dict(net.rno)
        return net

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
        ModflowGwfsfr(gwf, nreaches=nreaches,
                      packagedata=net.packagedata,
                      connectiondata=net.connectiondata,
                      perioddata={0: spd0},
                      unit_conversion=86400.0,   # Manning, SI, time unit = day
                      # MVR hands reach flow to the on-channel lakes and takes
                      # their spill back, and a package MVR names must declare
                      # MOVER itself -- otherwise MF6 stops on 'MODEL AND
                      # PACKAGE "…/SFR" DOES NOT HAVE MOVER SPECIFIED'. LAK
                      # already declared it; SFR did not (WP1d).
                      mover=bool(self.lak_mvr and self.ponds),
                      pname='sfr', save_flows=True,
                      budget_filerecord=f'{name}.sfr.cbc',
                      stage_filerecord=f'{name}.sfr.stage')

    def _build_ponds(self):
        """Load the pond polygons and give each a host cell and a lake geometry."""
        self.ponds = []
        if not self.lak_shapefile:
            return []
        from marmites_lak import (read_pond_polygons, assign_pond_cells,
                                  POND_DEPTH)
        cMF = self.cMF
        ponds = read_pond_polygons(self.lak_shapefile)
        assign_pond_cells(ponds, float(cMF.xllcorner), float(cMF.yllcorner),
                          cMF.delr, cMF.delc, self.nrow, self.ncol,
                          idomain=self.idomain,
                          verbose=getattr(self, 'verbose', True))
        depth = None if self.lak_depth is None else np.asarray(self.lak_depth, float)
        for p in ponds:
            i, j = p.cell
            d = POND_DEPTH
            if depth is not None and depth[i, j] > 0:
                d = float(depth[i, j])
            p.rim = float(self.top[i, j])
            p.bottom = p.rim - d
            # the lake connects at the topmost active layer, as the outcropping
            # unit is what a pond actually sits on
            k = 0
            while k < self.nlay - 1 and self.idomain[k, i, j] <= 0:
                k += 1
            p.klay = k
            p.on_channel = (i, j) in self.sfr_cells
            if p.on_channel and self.sfr_net is not None:
                p.inlet_reach = self.sfr_reach_of.get((i, j))
                down = self.sfr_net.recv.get((i, j))
                p.outlet_reach = None if down is None else self.sfr_reach_of.get(down)
        self.ponds = ponds
        return ponds

    def _add_lak_package(self, gwf, name, strt_heads=None):
        """ModflowGwflak: one EMBEDDEDV lake per pond.

        No evaporation is specified at BUILD time; the coupler writes it
        each stress period from the Eo forcing (WP1d), now that MARMITES no
        longer has a surface store of its own to evaporate from.
        """
        from marmites_lak import lake_table
        pkg, conn, tables, outlets = [], [], [], []
        wt = None if strt_heads is None else np.asarray(strt_heads, dtype=float)
        for L, p in enumerate(self.ponds):
            i, j = p.cell
            # start each lake at its local equilibrium stage: a perched pond
            # starts nearly empty, one in contact with the water table nearly
            # full, so neither gets a shock on the first time step
            if wt is not None:
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
            if self.lak_mvr and p.outlet_reach is not None:
                outlets.append([len(outlets), L, -1, 'MANNING', p.rim,
                                float(np.sqrt(max(p.area, 1.0))), 0.035, 1e-3])
        n = len(pkg)
        ModflowGwflak(gwf, pname='lak', boundnames=True,
                      print_stage=True, save_flows=True,
                      budget_filerecord=f'{name}.lak.cbc',
                      stage_filerecord=f'{name}.lak.stage',
                      mover=bool(outlets),
                      length_conversion=1.0, time_conversion=86400.0,
                      surfdep=self.lak_surfdep,
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

    def _add_mvr_package(self, gwf, name):
        """Route the stream through the on-channel ponds.

        A reach hosting a pond hands its flow to the lake, and the lake's
        outlet spills into the reach downstream of the pond cell. Without this
        the stream would simply bypass the pond.
        """
        recs = []
        for L, p in enumerate(self.ponds):
            if p.inlet_reach is not None:
                recs.append(['sfr', p.inlet_reach, 'lak', L, 'FACTOR', 1.0])
        for k, out in enumerate(getattr(self, 'lak_outlets', [])):
            p = self.ponds[out[1]]
            if p.outlet_reach is not None:
                recs.append(['lak', k, 'sfr', p.outlet_reach, 'FACTOR', 1.0])
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

    def build(self):
        cMF = self.cMF
        os.makedirs(self.sim_ws, exist_ok=True)
        name = cMF.modelname.lower()

        sim = MFSimulation(sim_name=name, sim_ws=self.sim_ws, exe_name=self.exe_name,
                           version='mf6')
        # TDIS: steady SP first, then transient
        perioddata = [(1.0, 1, 1.0)] + [(float(p), 1, 1.0) for p in self.perlen]
        # ATS: let MF6 subdivide any stress period it cannot solve in one step.
        # Without it the first transient day -- a step change from the steady
        # state onto a free-draining seepage boundary -- failed to converge and
        # discharged 9.6e6 m3 in a single step, 8x the recharge of the entire
        # run, leaving a cumulative budget discrepancy of -132%.
        ats = None
        if self.ats:
            recs = [(i, float(p), self.ats_dtmin, float(p), 2.0, 5.0)
                    for i, (p, _, _) in enumerate(perioddata)
                    if i >= 1]                     # not the steady-state period
            ats = {'maxats': len(recs), 'perioddata': recs}
        ModflowTdis(sim, time_units='DAYS', nper=self.nper, perioddata=perioddata,
                    ats_perioddata=ats)
        self.perioddata = perioddata

        # IMS from the NWT settings (COMPLEX option observed on La Mata)
        complexity = 'COMPLEX' if str(getattr(cMF, 'options', 'COMPLEX')).upper().startswith('COMPLEX') else 'MODERATE'
        self.outer_maximum = min(int(getattr(cMF, 'maxiterout', 500)), 500)
        ims = ModflowIms(sim, print_option='SUMMARY', complexity=complexity,
                         outer_dvclose=float(getattr(cMF, 'headtol', 0.05)),
                         outer_maximum=self.outer_maximum,
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
        ModflowGwfnpf(gwf, icelltype=list(np.asarray(cMF.laytyp, dtype=int)),
                      k=self._griddata(hk), k33=self._griddata(k33),
                      save_specific_discharge=False)

        # STO: steady first SP, transient afterwards
        ss = self._prop3d('ss_actual')
        sy = self._prop3d('sy_actual')
        ModflowGwfsto(gwf, iconvert=list(np.asarray(cMF.laytyp, dtype=int)),
                      ss=self._griddata(ss), sy=self._griddata(sy),
                      steady_state={0: True}, transient={1: True})

        # WEL: one well per active surface cell (ETg sink), q=0 initially,
        # AUTO_FLOW_REDUCE replaces the NWT 'SPECIFY 0.05 iunitramp' option
        wel_spd = [[self._cellid(k, i, j), 0.0] for (i, j, k) in self.surf_cells]
        ModflowGwfwel(gwf, stress_period_data={0: wel_spd},
                      auto_flow_reduce=0.05, pname='wel',
                      maxbound=self.ncell, save_flows=True)

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
                       # would give the water two exits
                       if not (self.sfr and (int(i), int(j)) in self.sfr_outlet_cells)]
            if drn_spd:
                ModflowGwfdrn(gwf, stress_period_data={0: drn_spd}, pname='drn',
                              maxbound=len(drn_spd), save_flows=True)

        # DRN-SEEP: groundwater seepage at the land surface as a smoothed drain
        # (decision 2/3). Alternative to UZF SIMULATE_GWSEEP, which is deprecated
        # in MF6 and switches discharge on and off discontinuously. The drain
        # sits at the land surface with AUXDEPTHNAME so MF6 applies its cubic
        # smoothing over DDRN, ramping the discharge in as the head rises
        # instead of snapping it on.
        #
        # These flows are NOT routed away with MVR (unlike the CdL reference
        # model, which sends them to SFR). MARMITES needs them back in the soil
        # column: the coupler reads the drain SIMVALS and feeds them to the
        # bottom soil layer as exfiltration, where the upward cascade of
        # Eq. 1/1b can turn the excess into Dunnian runoff.
        self.drnseep_id = {}
        if self.seep == 'drn':
            seep_spd = []
            for n, (i, j, k) in enumerate(self.surf_cells):
                if (i, j) in self.sfr_cells:      # SFR handles seepage there
                    continue
                self.drnseep_id[(i, j)] = len(seep_spd)
                seep_spd.append([self._cellid(k, i, j), float(self.top[i, j]),
                                 float(self.drn_seep_cond),
                                 float(self.drn_seep_ddrn)])
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
        thtr = float(np.ravel(np.asarray(cMF.thtr, dtype=float))[0])
        thts = float(np.ravel(np.asarray(cMF.thts, dtype=float))[0])
        thti = float(np.ravel(np.asarray(cMF.thti, dtype=float))[0])
        eps = float(np.ravel(np.asarray(cMF.eps, dtype=float))[0])
        surfdep = float(cMF.surfdep)
        thtr, thts, thti, eps = self._validate_uzf_params(thtr, thts, thti, eps)
        # vertical K for UZF: iuzfopt==1 -> vks array; iuzfopt==2 -> layer k33
        use_layer_vk = int(getattr(cMF, 'iuzfopt', 2)) != 1
        if not use_layer_vk:
            vks3d = self._lay3d(np.asarray(cMF.vks_actual, dtype=float),
                                self.nlay, self.nrow, self.ncol)
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
        for n, (i, j, k) in enumerate(self.surf_cells):
            chain = col_children[n]
            ivertcon = chain[0][0] if chain else -1
            vks = vks_scale * (float(k33[k, i, j]) if use_layer_vk else float(vks3d[k, i, j]))
            pkdata.append((n, self._cellid(k, i, j), 1, ivertcon, surfdep, vks,
                           thtr, thts, thti, eps))
            for ci, (no, kk) in enumerate(chain):
                child_ivert = chain[ci + 1][0] if ci + 1 < len(chain) else -1
                vks_c = vks_scale * (float(k33[kk, i, j]) if use_layer_vk else float(vks3d[kk, i, j]))
                subs.append((no, self._cellid(kk, i, j), 0, child_ivert, surfdep,
                             vks_c, thtr, thts, thti, eps))
        pkdata = pkdata + subs
        # period data: finf on land cells only; ET disabled (MARMITES does ET)
        pdata0 = [(n, float(getattr(cMF, 'perc_user', 0.0)), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                  for n in range(self.ncell)]
        ModflowGwfuzf(gwf, nuzfcells=self.nuzfcells, ntrailwaves=int(getattr(cMF, 'ntrail2', 7)),
                      nwavesets=int(getattr(cMF, 'nsets', 40)),
                      packagedata=pkdata, perioddata={0: pdata0},
                      simulate_et=False,
                      # exactly one seepage mechanism, never both, or the
                      # discharge would be counted twice
                      simulate_gwseep=(self.gwseep and self.seep == 'uzf'),
                      pname='uzf', save_flows=True,
                      budget_filerecord=f'{name}.uzf.cbc')

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
