# -*- coding: utf-8 -*-
"""MARMITES <-> MODFLOW 6 API coupler (Phase 3b/3c).

Replaces the removed file-based Picard loop with a single sequential march
over stress periods through the MODFLOW 6 API (libmf6, via modflowapi/xmipy).
The coupling is LAGGED: MMsoil.step() once per SP, using the heads and the
groundwater discharge (exfiltration) MF6 left at the END of the previous SP;
then MF6 advances the SP with its own time step, which ATS retries shorter
when it fails.

An ``iterative`` mode re-evaluated MMsoil at every MF6 outer iteration from
the current head iterate, under-relaxed (MetaSWAP-style). It went on
2026-10-07: what it existed for -- groundwater ET following the head within
the day -- EVT does inside MF6's own Newton solve, and because it drove the
solve by hand, a step it could not solve was accepted, never retried.

Exchange (memory pointers, no files):
  MM -> MF6 : perc (m/d)  -> UZF FINF/SINF on the land-surface UZF objects
              PETuzf      -> UZF PETMAX, less what EVT may take
              Eg, Tg      -> two EVT curves per cell (EVT_EG, EVT_TG), taken
                             by MF6 at the head it solves for (marmites_evt)
  MF6 -> MM : heads (m)   -> X (per surface cell, dry when h < botm_l0)
              exfiltration-> UZF GWD (m3/d) / cell_area * conv_fact (mm/d, +up)

The MF6 side must be built by ppMF6.marmites_mf6.clsMF6 so that the first
`ncell` UZF objects are the land-surface cells in MARMITES cell order and
the EVT records are in the same order. (Until 2026-10-07 groundwater ET
could also go out as a fixed WEL rate, one well per cell; that route is
gone and WEL is left for real pumping.)

The API object is injectable for testing (tests/test_coupler_mock.py runs
the coupler against a fake API); a real run needs libmf6:
    from modflowapi import ModflowApi
    MF6Coupler(...).run(ModflowApi('libmf6.so', working_directory=sim_ws))
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import os
import re
import time

import numpy as np


# a progress line at least this often [s], whatever the 5 % steps say
PROGRESS_EVERY_S = 900.0


def _hms(seconds):
    """A duration a modeller can act on: 4h 42m, not 16920.3."""
    s = int(max(seconds, 0))
    if s >= 3600:
        return '%dh %02dm' % (s // 3600, (s % 3600) // 60)
    if s >= 60:
        return '%dm %02ds' % (s // 60, s % 60)
    return '%ds' % s


def surface_excess_lines(heads, top, ncell):
    """The warning for a water table standing above the land surface.

    ``heads`` is (nper, ncell), ``top`` (ncell,). Silent unless some head
    is more than 1 m above ground. It used to count only cells MORE THAN
    1 m above -- and so reported "1176 cells, 1 of 60 stress periods" for
    a run whose water table sat above the land surface in thousands of
    cells for all 60 (the UZF/Sy runaway of 2026-09-23). A seepage face
    holds the head a little above ground, so the count that says how much
    of the catchment is flooded is > 0.1 m; > 1 m is kept alongside.
    """
    above = np.asarray(heads, float) - np.asarray(top, float)[None, :]
    if not above.size or np.nanmax(above) <= 1.0:
        return []
    nper = above.shape[0]
    wet, far = above > 0.1, above > 1.0
    per_sp = wet.sum(axis=1)
    k = int(np.argmax(per_sp))
    return ['',
            'WARNING: the water table rises up to %.1f m ABOVE the land '
            'surface.' % float(np.nanmax(above)),
            '         more than 0.1 m above: %d of %d cells at some time, in '
            '%d of %d stress periods; at worst %d cells at once (SP %d)'
            % (int(wet.any(axis=0).sum()), ncell, int(wet.any(axis=1).sum()),
               nper, int(per_sp[k]), k + 1),
            '         more than 1 m above:   %d cells, in %d stress period(s)'
            % (int(far.any(axis=0).sum()), int(far.any(axis=1).sum()))]


class CouplingError(Exception):
    """Raised on MF6 API exchange failures."""


class MF6Coupler:
    """Sequential per-SP coupling driver over MMsoil.step().

    Parameters
    ----------
    mm, ctx, state : clsMMsoil instance, build_context() namespace, init_state()
        The Phase-2 soil model with its static context and carry-over state.
        ctx.cells must be in the same order as clsMF6.surf_cells.
    mf6b : ppMF6.marmites_mf6.clsMF6
        The built MF6 model (for names, cell mapping, geometry).
    conv_fact : float
        Length conversion, mm per model length unit (1000 for metres).
    max_outer : int
        Safety cap on outer iterations per step when MF6 is stepped through
        prepare_solve/solve (an MF6 build without the UZF period-input
        arrays; see _bind) -- MF6's own do_time_step needs none.
    """

    # the coupled run's results, in the model workspace. The name kept the
    # coupling mode when there were two; it stays, so the Results page and
    # the post-processing find the runs made before and after 2026-10-07
    RESULTS_H5 = '_coupled_lagged.h5'

    def __init__(self, mm, ctx, state, mf6b, conv_fact=1000.0,
                 max_outer=None, obs_idx=None, obs_names=None, crr=None):
        self.mm, self.ctx, self.state = mm, ctx, state
        # MMsoil hands over what Eg and Tg COULD take, for the day's EVT
        # curves, instead of taking them itself with its own drawdown: that
        # is the coupled run's contract (the uncoupled MMsoil leaves it off)
        ctx.gw_evt = True
        # Observation cells at which to keep the FULL per-SP MM flux vectors
        # (for the native per-point Sankey / time series). Tiny next to the
        # whole-grid arrays, so unlike wb_map/wb_ts these are stored in full.
        self.obs_idx = None if obs_idx is None else list(np.asarray(obs_idx, int))
        self.obs_names = list(obs_names) if obs_names is not None else None
        self.mf6b = mf6b
        self.conv_fact = float(conv_fact)
        # The manual outer-iteration loop must let MF6 reach its OWN
        # OUTER_MAXIMUM before finalizing, otherwise MF6 never registers the
        # step as failed and ATS cannot retry it with a smaller dt. Capping
        # below OUTER_MAXIMUM was why the first transient day advanced
        # unconverged (see the Phase-5 -130% budget). Default to the solver's
        # own limit; a caller can still override.
        if max_outer is None:
            max_outer = int(getattr(mf6b, 'outer_maximum', 500))
        self.max_outer = int(max_outer)
        self.name = mf6b.cMF.modelname.upper()
        cMF = mf6b.cMF
        # cell mapping / geometry
        self.ncell = ctx.ncell
        if mf6b.ncell != self.ncell:
            raise CouplingError('cell count mismatch: MM %d vs MF6 %d'
                                % (self.ncell, mf6b.ncell))
        self.i_arr = np.array([c[1] for c in ctx.cells], dtype=int)
        self.j_arr = np.array([c[2] for c in ctx.cells], dtype=int)
        self.k_arr = np.array([mf6b.surf_cells[n][2] for n in range(self.ncell)], dtype=int)
        # Phase 4: cell area from the shared grid geometry, so the soil model
        # and the coupler cannot disagree on it (DIS and DISV alike).
        geom = getattr(ctx, 'geom', None)
        if geom is not None:
            self.area = np.asarray(geom.area, dtype=float)
        else:  # pragma: no cover - legacy contexts without a geometry
            delr = np.asarray(cMF.delr, dtype=float)
            delc = np.asarray(cMF.delc, dtype=float)
            self.area = delc[self.i_arr] * delr[self.j_arr]
        self.grid_kind = getattr(geom, 'kind', 'structured')
        self.ncpl = int(getattr(geom, 'ncpl', mf6b.nrow * mf6b.ncol))
        # results
        self.heads_hist = None
        self.exf_hist = None
        self.perc_hist = None
        self.etg_hist = None
        self.outer_iters = None
        self._warned_dt = False
        # per-cell mean recharge / groundwater ET to drive the steady SP0 (set
        # from a prior run); None -> uniform perc_user fallback
        self.steady_perc = None
        self.steady_etg = None
        _perlen = getattr(mf6b, 'perlen', None)
        if _perlen is None:
            _perlen = getattr(cMF, 'perlen', [])
        self.perlen = list(np.asarray(_perlen).ravel())
        self.sim_ws = getattr(mf6b, 'sim_ws', None)
        self.max_substeps = 5000       # ATS safety cap per stress period
        self.substeps = 0              # extra ATS sub-steps taken overall
        self.n_nonconverged = 0
        # MF6's own time step (do_time_step), which RETRIES a step it cannot
        # solve with a smaller one (ATS); set by _bind when the period-input
        # arrays it re-applies on every try are reachable
        self.mf6_step = False
        self.p_sinf_pvar = self.p_pet_pvar = self.p_itertot = None
        # When the march started, for the progress line's estimate. Set here
        # so _progress can never be called before it exists, and reset in
        # run() so a spin-up cycle times itself rather than the whole session.
        self._t0 = time.time()
        self.p_gwd = None
        self.p_rejinf = None
        # WP2: UZF's PET demand (written) and its actual ET (read back)
        self.p_petmax = self.p_pet = self.p_uzet = None
        self.p_wcnew = None
        self.uzf_wc_final = None
        # the lake stages when the run ends (LAK XNEWPAK): where a spin-up's
        # next cycle and a saved state start each pond
        self.p_lak_stage = None
        self.lak_stage_final = None
        # a periodic spin-up: the previous cycle's last seepage, rejected
        # infiltration and UZF ET, read on the first day in its place
        self.carry_in = None
        self.carry_out = None
        self.etuzf_prev = np.zeros(self.ncell)
        self.etuzf_hist = None
        # the UZF objects of each land column (land object first): a
        # column's actual ET is the sum over them -- MF6 hands the demand
        # the land object leaves to the object below (setbelowpet)
        cols = getattr(mf6b, 'uzf_columns', None) or [[k] for k in
                                                      range(self.ncell)]
        self._col_obj = np.array([o for c in cols for o in c], dtype=int)
        self._col_own = np.array([k for k, c in enumerate(cols) for _o in c],
                                 dtype=int)
        self.pet_unmet = None           # catchment PT+PE - (ETsoil+ETuzf+ETg)
        self.n_overdraw = 0
        self.max_overdraw = 0.0
        # the period just advanced: its rates' sub-step means (None when it
        # was one step) and its step count -- see SUBSTEP_RATES
        self._period_rates = None
        self._period_nsub = 1
        # UZF's ET beyond the demand within MF6's wave tolerance (_split_etuzf)
        self.n_resid, self.resid_max, self.resid_m3 = 0, 0.0, 0.0
        self._petmax_written = None
        self.p_drnseep = None
        # groundwater ET: two EVT packages at the solved head (marmites_evt)
        self.evt_nseg = int(getattr(mf6b, 'evt_nseg', 8))
        self.evt_ramp = float(getattr(mf6b, 'evt_ramp', 0.1))
        self.p_evt = {}
        self.p_evt_eg_sim = self.p_evt_tg_sim = None
        # the MF6 top per cell: the BASE of the MMsoil column (UZF's top)
        top = getattr(mf6b, 'top', None)
        self.top_cell = (None if top is None else
                         np.asarray(top, dtype=float)[self.i_arr, self.j_arr])
        # ... and the LAND surface, for the water-table plausibility check
        # (falls back to the top for a builder without one)
        land = (mf6b._land_surface() if hasattr(mf6b, '_land_surface')
                else None)
        self.land_cell = (self.top_cell if land is None else
                          np.asarray(land, dtype=float)[self.i_arr, self.j_arr])
        # SFR: reach number per MARMITES cell (-1 = no reach in this cell), so
        # the runoff of a channel cell can be injected as reach inflow.
        self.p_sfr_inflow = None
        reach_of = getattr(mf6b, 'sfr_reach_of', None) or {}
        self.sfr_reach_idx = np.full(self.ncell, -1, dtype=int)
        for n in range(self.ncell):
            r = reach_of.get((self.i_arr[n], self.j_arr[n]))
            if r is not None:
                self.sfr_reach_idx[n] = int(r)
        self.nreaches = len(set(reach_of.values()))
        # DRN-SEEP boundary index per MARMITES cell (-1 = no seepage drain,
        # e.g. a cell where SFR takes the seepage instead). The DRN package
        # arrays are ordered by boundary, not by cell, so the two orders must
        # be reconciled explicitly rather than assumed equal.
        seep_id = getattr(mf6b, 'drnseep_id', None) or {}
        self.drnseep_idx = np.full(self.ncell, -1, dtype=int)
        for n in range(self.ncell):
            b = seep_id.get((self.i_arr[n], self.j_arr[n]))
            if b is not None:
                self.drnseep_idx[n] = int(b)

        # ---- open-water evaporation (WP1d) -------------------------------
        # Eo used to be applied by MMsoil to its own surface store, and SFR and
        # LAK were built with evaporation deliberately left at 0 so the same
        # water would not leave twice. The surface store went with the pond
        # module, so the evaporation follows the water: onto the stream reaches
        # and the lakes. Eo is per METEO ZONE, so each reach and each lake is
        # tagged with the zone of the cell it sits in.
        self.p_sfr_evap = None
        self.p_lak_evap = None
        self.p_lak_runoff = None        # WP4.4: MM runoff into the ponds
        self.p_sfr_simevap = None
        self.p_lak_simevap = None
        # The stream's Eo is capped at a share of each reach's inflow (sfr.
        # evap_inflow_fraction, analysis §8.23): a reach receiving less than
        # it could evaporate, over a water table below its bed, flip-flops
        # in MF6 between evaporating all of it and leaking all of it, and
        # sfr_fn's perturbed derivative straddles the two -- the solver then
        # cannot converge at a tight outer_dvclose.
        self.sfr_evap_frac = float(getattr(mf6b, 'sfr_evap_frac', 0.5))
        self.p_sfr_usflow = self.p_sfr_qfrommvr = None
        self.sfr_wl = None              # each reach's width x length [m2]
        self.sfr_up_prev = None         # its upstream + mover inflow [m3/d]
        self.sfr_evap_capped = 0        # reach-steps whose EVAP was lowered
        self.sfr_evap_writes = 0        # reach-steps written
        gm = getattr(ctx, 'gridMETEO', None)
        gm = None if gm is None else np.asarray(gm)
        self.sfr_evap_zone = np.zeros(max(self.nreaches, 1), dtype=int)
        if gm is not None and self.nreaches:
            for (i, j), r in reach_of.items():
                if 0 <= int(r) < self.nreaches:
                    self.sfr_evap_zone[int(r)] = int(gm[i, j]) - 1
        ponds = list(getattr(mf6b, 'ponds', None) or [])
        self.nlakes = len(ponds)
        self.lak_evap_zone = np.zeros(max(self.nlakes, 1), dtype=int)
        # MARMITES cell index of each lake, for reporting its evaporation back
        # into the per-cell water balance.
        cell_of = {(int(self.i_arr[n]), int(self.j_arr[n])): n
                   for n in range(self.ncell)}
        self.lak_cell_idx = np.full(max(self.nlakes, 1), -1, dtype=int)
        if gm is not None and self.nlakes:
            for L, p in enumerate(ponds):
                i, j = p.cell
                self.lak_evap_zone[L] = int(gm[i, j]) - 1
                self.lak_cell_idx[L] = cell_of.get((int(i), int(j)), -1)

        # ---- WP4.6: the surface descriptor (cookbook §4a) ----------------
        # Where the MM soil column runs and where open water does instead,
        # per cell: MMsoil reads ctx.f_open, so it is set before the first
        # stress period. Built once, from the network and the ponds as
        # written, on the cell list as ordered -- the list is not touched.
        self.f_lake = np.zeros(self.ncell, dtype=float)
        self.f_stream = np.zeros(self.ncell, dtype=float)
        if hasattr(mf6b, 'surface_fractions'):
            self.f_lake, self.f_stream = mf6b.surface_fractions(
                list(zip(self.i_arr, self.j_arr)), self.area)
        f_open = self.f_lake + self.f_stream
        ctx.f_open = f_open if np.any(f_open > 0.0) else None
        # each lake's FOOTPRINT, weighted by the pond area it holds: where
        # its evaporation is booked and its runoff collected (WP4.4). A pond
        # built without a footprint (an older build) is its host cell.
        self.lak_cells = []
        for L, p in enumerate(ponds):
            ks = [cell_of[c] for c in (getattr(p, 'cells', None) or [p.cell])
                  if tuple(int(v) for v in c) in cell_of]
            ks = np.asarray(ks, dtype=int)
            w = (self.f_lake[ks] * self.area[ks]) if ks.size else np.zeros(0)
            if ks.size and w.sum() <= 0.0:
                w = self.area[ks].copy()
            self.lak_cells.append((ks, w / w.sum() if w.sum() > 0 else w))
        if ctx.f_open is not None and getattr(mf6b, 'verbose', True):
            wa = self.area / self.area.sum()
            lk, st = self.f_lake > 0, self.f_stream > 0
            print('surface: the MM soil column runs on %.2f %% of the catchment;'
                  ' %d pond cell(s) (f_lake %.2f-%.2f), %d channel cell(s) '
                  '(f_stream median %.2f)'
                  % (100.0 * float(np.sum(wa * (1.0 - f_open))), int(lk.sum()),
                     float(self.f_lake[lk].min()) if lk.any() else 0.0,
                     float(self.f_lake[lk].max()) if lk.any() else 0.0,
                     int(st.sum()),
                     float(np.median(self.f_stream[st])) if st.any() else 0.0))
        # The zone rasters are 1-based; a cell outside every mapped zone would
        # index -1 and silently take the LAST zone's Eo.
        self.sfr_evap_zone = np.clip(self.sfr_evap_zone, 0, None)
        self.lak_evap_zone = np.clip(self.lak_evap_zone, 0, None)
        eo = getattr(ctx, 'Eo_zonesSP', None)
        self.Eo_zonesSP = None if eo is None else np.atleast_2d(np.asarray(eo,
                                                                          float))
        self.evap_hist = None

        # ---- WP5: the CRR cascade -----------------------------------------
        # Built once, on the surface descriptor above (cookbook §4a: one
        # mechanism for WP4 and WP5) -- a pond cell is a LAK receiver, a
        # channel cell an SFR one, every other cell a soil column.
        self.crr = None
        self.crr_ts = None
        if crr is not None:
            self.crr = self._build_crr(crr)
            if getattr(mf6b, 'verbose', True) and crr.get('report', True):
                for _line in self.crr.summary(self.area):
                    print(_line)

    def _build_crr(self, opts):
        """The CascadeNetwork of this run, from ``opts`` = {'beta', 'sinks',
        and optionally 'elev'}: the elevation the slopes come from, on the
        model grid (nrow, ncol); by default the land surface (cMF.elev)."""
        import marmites_crr as mcrr
        ctx, mf6b = self.ctx, self.mf6b
        if any(int(c[0]) != k for k, c in enumerate(ctx.cells)):
            raise CouplingError('CRR addresses cells by their place in the '
                                'cell list, and the list is not numbered in '
                                'order')
        kind = np.full(self.ncell, mcrr.KIND_SOIL, dtype=int)
        kind[self.sfr_reach_idx >= 0] = mcrr.KIND_SFR
        for ks, _w in self.lak_cells:              # LAK > SFR > soil
            kind[ks] = mcrr.KIND_LAK
        elev = opts.get('elev')
        if elev is None:
            elev = mf6b._land_surface()
        elev = np.asarray(np.ma.filled(np.ma.asarray(elev, dtype=float),
                                       np.nan), dtype=float)
        elev = elev.reshape(-1, 1) if elev.ndim == 1 else elev
        elev = elev[self.i_arr, self.j_arr]
        topo = mf6b.topology() if hasattr(mf6b, 'topology') else None
        if topo is not None:                       # DISV, (ncpl, 1)
            icell = self.i_arr
        else:
            cMF = mf6b.cMF
            topo = mcrr.structured_topology(cMF.delr, cMF.delc)
            icell = self.i_arr * int(np.size(cMF.delr)) + self.j_arr
        return mcrr.network_from_topology(
            topo, icell, elev, kind, beta=float(opts.get('beta', 1.0)),
            sinks=str(opts.get('sinks', 'evaporate')))

    # the CRR terms kept per stress period [m3/d], catchment totals
    CRR_TERMS = ('runoff', 'runon', 'reinf', 'to_sfr', 'to_lak', 'evap',
                 'evap_sink', 'routed')

    def _crr_record(self, n, out):
        """This period's cascade, catchment totals in m3/d (``crr_ts[n]``)."""
        import marmites_crr as mcrr
        res = out.get('crr')
        if res is None or self.crr_ts is None:
            return
        kind = self.crr.kind
        ri = self.ctx.index.get('iReinf')
        reinf = (0.0 if ri is None else
                 float(np.sum(np.asarray(out['MM'])[:, ri] / self.conv_fact
                              * self.area)))
        self.crr_ts[n] = (float(res.ro.sum()), float(res.runon.sum()), reinf,
                          float(res.deliver[kind == mcrr.KIND_SFR].sum()),
                          float(res.deliver[kind == mcrr.KIND_LAK].sum()),
                          float(res.ecrr.sum()), float(res.ecrr_sink.sum()),
                          float(res.routed.sum()))

    def _print_crr(self, nper):
        """What the cascade did over the run, in mm/yr over the catchment."""
        if self.crr_ts is None or not self.crr_ts.size:
            return
        dt = np.array([float(self.perlen[n]) if n < len(self.perlen) else 1.0
                       for n in range(nper)])
        days = float(dt.sum())
        if days <= 0.0:
            return
        tot = (self.crr_ts * dt[:, None]).sum(axis=0)         # m3
        mm = dict(zip(self.CRR_TERMS, tot / float(self.area.sum()) * 1000.0
                      * 365.25 / days))
        # a cell's runoff = its own excess + the run-on it passed along, so
        # the excess the cells GENERATED is runoff - run-on + reinfiltration
        gen = mm['runoff'] - mm['runon'] + mm['reinf']
        print('\nCRR cascade (mm/yr over the catchment): runoff generated '
              '%.1f = reinfiltrated %.1f + to the streams %.1f + to the ponds '
              '%.1f + evaporated %.1f (at sinks %.1f)'
              % (gen, mm['reinf'], mm['to_sfr'], mm['to_lak'], mm['evap'],
                 mm['evap_sink']))
        print('     run-on received by the soil columns %.1f%s'
              % (mm['runon'], '; %.1f sent from sinks to the nearest stream '
                 'or pond' % mm['routed'] if mm['routed'] else ''))

    def _mm_step(self, n, tstart_MF, heads, exf, rej):
        """One MMsoil stress period. WP2: Eg/Tg see what remains after the
        PREVIOUS period's actual UZF ET (cookbook 2b); WP5: the cascade
        routes the runoff when there is one; then the day's EVT curves, and
        the cap they put on UZF's demand."""
        if self.crr is None:
            out = self.mm.step(self.ctx, n, tstart_MF, heads, exf, self.state,
                               rejinf_cell=rej, etuzf_cell=self.etuzf_prev)
        else:
            out = self.mm.step(self.ctx, n, tstart_MF, heads, exf, self.state,
                               rejinf_cell=rej, etuzf_cell=self.etuzf_prev,
                               crr=self.crr)
        if out.get('gw') is not None:
            out['evt'] = self._evt_curves(out['gw'])
            # what EVT can take today, at most: reserved from UZF's demand,
            # so ETsoil + ETuzf + ETg <= PE + PT holds whatever MF6 solves
            out['etg_cap'] = out['evt']['eg']['RATE'] + out['evt']['tg']['RATE']
        return out

    def _evt_curves(self, gw):
        """The day's two EVT curves per cell from MMsoil's groundwater-ET
        potential (marmites_evt): {'eg'|'tg': {SURFACE, RATE, DEPTH, PXDP,
        PETM}} in m and m/d, in the cell order the packages were built in."""
        import marmites_evt as me
        n, ns = self.ncell, self.evt_nseg
        out = {t: {'SURFACE': np.zeros(n), 'RATE': np.zeros(n),
                   'DEPTH': np.ones(n), 'PXDP': np.zeros((n, ns - 1)),
                   'PETM': np.zeros((n, ns - 1))} for t in ('eg', 'tg')}
        for k in range(n):
            curves = (('eg', me.eg_curve(gw['eg_pe'][k], gw['land'][k],
                                         gw['h0'][k], gw['shah'][k], ns)),
                      ('tg', me.tg_curve(gw['tg_rate'][k], gw['tg_tip'][k],
                                         gw['h0'][k], self.evt_ramp, ns)))
            for t, (s, r, d, x, y) in curves:
                c = out[t]
                c['SURFACE'][k], c['RATE'][k], c['DEPTH'][k] = s, r, d
                c['PXDP'][k], c['PETM'][k] = x, y
        return out

    def _evt_flat(self, rate):
        """A STEADY first period's EVT curves: ``rate`` [m/d per cell area],
        each cell's mean ETg, flat down to the bottom of its surface cell
        (marmites_evt.flat_curve) -- the fixed mean rate the ETg well drew
        until 2026-10-07, reduced as the cell dries as AUTO_FLOW_REDUCE did.
        Booked on EVT_EG as a whole: a steady period's mean has no Eg/Tg
        split, and no result is read from it. Tg's curve is inert."""
        import marmites_evt as me
        n, ns = self.ncell, self.evt_nseg
        rate = np.asarray(rate, dtype=float).ravel()[:n]
        surf = (self.land_cell if self.land_cell is not None
                else np.zeros(n))
        surf = np.asarray(surf, dtype=float)
        # never a NaN into MF6: a cell without a land elevation takes the
        # model top, and failing that a curve that is inert (rate 0)
        if self.top_cell is not None:
            surf = np.where(np.isfinite(surf), surf, self.top_cell)
        surf = np.nan_to_num(surf, nan=0.0)
        botm = getattr(self.mf6b, 'botm', None)
        if botm is not None:
            b = np.asarray(botm, dtype=float)
            bot = b[self.k_arr, self.i_arr, self.j_arr]
        else:                       # a builder without bottoms: no floor
            bot = surf - 1.0e3
        out = {t: {'SURFACE': surf.copy(), 'RATE': np.zeros(n),
                   'DEPTH': np.ones(n), 'PXDP': np.zeros((n, ns - 1)),
                   'PETM': np.zeros((n, ns - 1))} for t in ('eg', 'tg')}
        x = np.linspace(0.0, 1.0, ns + 1)[1:-1]
        out['tg']['PXDP'][:] = x
        for k in range(n):
            s, r, d, xx, y = me.flat_curve(rate[k], surf[k], bot[k], ns)
            c = out['eg']
            c['SURFACE'][k], c['RATE'][k], c['DEPTH'][k] = s, r, d
            c['PXDP'][k], c['PETM'][k] = xx, y
        return out

    def _write_evt(self, evt):
        """Write the day's EVT curves (``_evt_curves``) into MF6."""
        if not evt or not self.p_evt:
            return
        n = self.ncell
        for t, arrs in evt.items():
            p = self.p_evt.get(t) or {}
            for var in ('SURFACE', 'RATE', 'DEPTH'):
                if p.get(var) is not None:
                    p[var][:n] = arrs[var]
            for var in ('PXDP', 'PETM'):
                if p.get(var) is not None:
                    p[var][:n] = np.asarray(arrs[var]).reshape(
                        p[var][:n].shape)

    def _evt_into(self, mm_cells):
        """The groundwater ET MF6's EVT took this period, into the MM vector
        (iEg, iTg, iETg, iETtot) [mm/d per cell]; returns ETg [m/d] for the
        run's books. Sub-step means, as for UZF."""
        if self.p_evt_eg_sim is None:
            return np.zeros(self.ncell)
        n = self.ncell
        eg = np.abs(np.asarray(self._rate('p_evt_eg_sim'), dtype=float)
                    .ravel()[:n]) / self.area * self.conv_fact
        tg = np.abs(np.asarray(self._rate('p_evt_tg_sim'), dtype=float)
                    .ravel()[:n]) / self.area * self.conv_fact
        ix = self.ctx.index
        for key, v in (('iEg', eg), ('iTg', tg), ('iETg', eg + tg)):
            if key in ix:
                mm_cells[:, ix[key]] = v
        if 'iETtot' in ix:
            mm_cells[:, ix['iETtot']] += eg + tg
        return (eg + tg) / self.conv_fact

    def _runoff_to_deliver(self, out):
        """The runoff SFR and LAK take this period [mm/d per cell]: with CRR
        what the cascade brought to each stream and pond cell, without it
        each such cell's own runoff (the soil cells' was not delivered)."""
        if self.crr is not None:
            return np.asarray(out['ro_deliver'], dtype=np.float64)
        return np.asarray(out['MM'], dtype=np.float64)[:, self._iRo]

    # ---------------- address / pointer helpers ------------------------ #

    def _addr(self, api, var, comp):
        """Resolve a variable address with clear failure diagnostics."""
        try:
            return api.get_var_address(var, *comp.split('/'))
        except Exception as exc:  # pragma: no cover
            raise CouplingError(f'cannot resolve MF6 variable {var} in {comp}: {exc}') from exc

    @staticmethod
    def available_vars(api, contains=None):
        """All MF6 input variable addresses (optionally filtered)."""
        try:
            names = list(api.get_input_var_names())
        except Exception:  # pragma: no cover - very old API
            return []
        if contains:
            up = contains.upper()
            names = [n for n in names if up in n.upper()]
        return names

    def _bind_first(self, api, candidates, what, required=True, min_size=None):
        """Bind the first address that actually exists.

        MF6 memory-variable names differ between versions and packages, and a
        wrong address can hand back a pointer of the wrong size -- writing
        through it corrupts MF6's memory and kills the process with no error
        message. So a bound pointer is validated by SIZE (``min_size``).

        The variable list from ``get_input_var_names()`` is used as a fast
        filter, but it is NOT authoritative: advanced-package variables such as
        UZF ``SINF`` are reachable through ``get_value_ptr`` yet are absent from
        that list. Requiring membership silently rejected SINF and fell back to
        the non-operative FINF, so the daily recharge never reached UZF. Now a
        candidate that is missing from the list is still accepted if its pointer
        resolves and has a plausible size.
        """
        known = set(self._known_vars)
        tried = []
        for addr_parts in candidates:
            var, comp = addr_parts
            try:
                addr = api.get_var_address(var, *comp.split('/'))
            except Exception as exc:
                tried.append(f'{comp}/{var} (address error: {exc})')
                continue
            try:
                ptr = api.get_value_ptr(addr)
            except Exception as exc:
                tried.append(f'{addr} (pointer error: {exc})')
                continue
            try:
                size = int(np.asarray(ptr).size)
            except Exception:
                size = None
            if min_size is not None and size is not None and size < min_size:
                tried.append(f'{addr} (size {size} < {min_size})')
                continue
            if known and addr not in known:
                # reachable pointer that MF6 does not list as an input var
                # (e.g. UZF SINF) -- accept it, but say so
                print('NOTE: bound %s via %s (size %s; not in the input-var '
                      'list but a valid pointer)' % (what, addr, size))
            return ptr, addr
        if not required:
            return None, None
        hint = self.available_vars(api, contains=candidates[0][1].split('/')[-1])
        raise CouplingError(
            'could not bind %s.\nTried:\n  %s\nVariables exposed by that package:\n  %s\n'
            'Send this list and I will fix the address mapping.'
            % (what, '\n  '.join(tried), '\n  '.join(hint[:40]) or '(none found)'))

    def _bind(self, api):
        """Bind the exchange pointers once after initialize()."""
        name = self.name
        self._known_vars = self.available_vars(api)
        self.p_x, _ = self._bind_first(api, [('X', name)], 'heads (X)')
        # UZF INFILTRATION: TWO ARRAYS, BOTH WRITTEN. MF6's UZF cell group
        # (UzfCellGroup.f90, MF6 6.7 source) keeps FINF, which the kinematic
        # wave ROUTES (surflux = finf + mover), and SINF, which the budget's
        # INFILTRATION line REPORTS (gwf-uzf.f90, appliedinf). uzf_ad sets
        # both from the period input SINF_PVAR (setdatafinf), and
        # xmi_prepare_solve runs the *_ad routines -- so they are written
        # after prepare_solve, together, as setdatafinf would.
        #
        # This used to bind SINF alone, on the belief that SINF was the
        # 'operative' array and FINF 'non-operative' -- drawn from the budget
        # line, which does follow SINF. The routing never did: on the run of
        # 2026-09-23 UZF's outflows were a constant 960.9 m3/d in every step
        # (the input file's steady rate, 961.0), the UZF package budget was
        # 92 % out, and MMsoil's percolation never entered the unsaturated
        # zone. Writing FINF before prepare_solve (the earlier attempt) was
        # reverted by uzf_ad, which is what made FINF look inert.
        self.p_finf, self.addr_finf = self._bind_first(
            api, [('FINF', f'{name}/UZF'), ('FINF', f'{name}/UZF-1')],
            'UZF routed infiltration (FINF)', min_size=self.ncell)
        self.p_sinf, self.addr_sinf = self._bind_first(
            api, [('SINF', f'{name}/UZF'), ('SINF', f'{name}/UZF-1')],
            'UZF reported infiltration (SINF)', min_size=self.ncell)
        print('coupler: UZF infiltration bound to %s (routed) and %s '
              '(budget)' % (self.addr_finf, self.addr_sinf))
        # groundwater discharge to land surface (= MARMITES exfiltration).
        # Name varies by MF6 version; optional -- without it exfiltration is 0.
        self.p_gwd, self.addr_gwd = self._bind_first(
            api, [('GWD', f'{name}/UZF'), ('GWD', f'{name}/UZF-1'),
                  ('SIMVALS', f'{name}/UZF'), ('SIMVALS', f'{name}/UZF-1')],
            'UZF groundwater discharge', required=False)
        # DRN-SEEP alternative (seep='drn'): the land-surface drain package
        # carries the seepage instead of UZF. Its SIMVALS are negative (water
        # leaving the aquifer), so the sign is flipped on read.
        self.p_drnseep, self.addr_drnseep = self._bind_first(
            api, [('SIMVALS', f'{name}/DRN_SEEP')],
            'DRN-SEEP discharge', required=False)
        # GROUNDWATER ET: the two EVT packages' curve arrays, written every
        # period before MF6's own time step (nothing resets them: only a
        # PERIOD block does -- probed on libmf6 6.7, evt_toy2.py), and the
        # ET they took (SIMVALS, negative out of the aquifer). Required: it
        # is the coupled run's only groundwater-ET path since 2026-10-07.
        self.p_evt = {}
        for tag in ('eg', 'tg'):
            pk = 'EVT_%s' % tag.upper()
            arrs = {}
            for var in ('SURFACE', 'RATE', 'DEPTH', 'PXDP', 'PETM',
                        'SIMVALS'):
                arrs[var], _a = self._bind_first(
                    api, [(var, f'{name}/{pk}')], '%s %s' % (pk, var),
                    min_size=self.ncell)
            self.p_evt[tag] = arrs
        self.p_evt_eg_sim = self.p_evt['eg']['SIMVALS']
        self.p_evt_tg_sim = self.p_evt['tg']['SIMVALS']
        print('coupler: groundwater ET by EVT (EVT_EG, EVT_TG), %d '
              'segments, taken at the solved head' % self.evt_nseg)
        # SFR specified inflow, written each SP with the MARMITES runoff
        if self.nreaches:
            self.p_sfr_inflow, self.addr_sfr_inflow = self._bind_first(
                api, [('INFLOW', f'{name}/SFR'), ('INFLOW', f'{name}/SFR-1')],
                'SFR specified inflow', required=False)
            if self.p_sfr_inflow is None:
                print('WARNING: SFR INFLOW array not exposed by this MF6 build; '
                      'MARMITES runoff will not reach the stream.\n'
                      '         (SFR variables seen: %s)'
                      % ', '.join(self.available_vars(api, 'SFR')[:12]))
            elif self.p_sfr_inflow.size < self.nreaches:
                print('WARNING: SFR INFLOW has %d entries but the network has %d '
                      'reaches; runoff routing disabled.'
                      % (self.p_sfr_inflow.size, self.nreaches))
                self.p_sfr_inflow = None
            # WP1d: open-water evaporation, which MMsoil no longer applies.
            # EVAP is the INPUT rate per unit of wetted area [m/d]; SIMEVAP is
            # what MF6 actually removed, as a volumetric rate [m3/d], capped
            # by what the reach holds. Probed, not assumed: writing 0.005 m/d
            # gave SIMEVAP 1.048 m3/d on the widest reach, whose wetted area
            # is 70.71 x 3.0 = 212 m2 -> 1.06.
            self.p_sfr_evap, self.addr_sfr_evap = self._bind_first(
                api, [('EVAP', f'{name}/SFR'), ('EVAP', f'{name}/SFR-1')],
                'SFR evaporation', required=False)
            if self.p_sfr_evap is None:
                print('WARNING: SFR EVAP array not exposed by this MF6 build; '
                      'open-water evaporation will NOT be applied to the '
                      'stream.')
            self.p_sfr_simevap, _a = self._bind_first(
                api, [('SIMEVAP', f'{name}/SFR'), ('SIMEVAP', f'{name}/SFR-1')],
                'SFR simulated evaporation', required=False)
            # the cap on that evaporation (§8.23): what flows into each reach
            # besides its runoff -- the upstream flow and the pond outflow
            # (QFROMMVR, at SFR's own path, only with MOVER) -- and its area
            self.p_sfr_usflow, _a = self._bind_first(
                api, [('USFLOW', f'{name}/SFR'), ('USFLOW', f'{name}/SFR-1')],
                'SFR upstream inflow', required=False, min_size=self.nreaches)
            self.p_sfr_qfrommvr, _a = self._bind_first(
                api, [('QFROMMVR', f'{name}/SFR'),
                      ('QFROMMVR', f'{name}/SFR-1')],
                'SFR inflow from the mover', required=False,
                min_size=self.nreaches)
            _len, _a = self._bind_first(
                api, [('LENGTH', f'{name}/SFR'), ('LENGTH', f'{name}/SFR-1')],
                'SFR reach length', required=False, min_size=self.nreaches)
            _wid, _a = self._bind_first(
                api, [('WIDTH', f'{name}/SFR'), ('WIDTH', f'{name}/SFR-1')],
                'SFR reach width', required=False, min_size=self.nreaches)
            if self.p_sfr_evap is not None and self.sfr_evap_frac <= 0.0:
                print('coupler: no open-water evaporation from the streams '
                      '(sfr.evap_inflow_fraction = 0)')
            elif self.p_sfr_evap is not None:
                if self.p_sfr_usflow is None or _len is None or _wid is None:
                    print('WARNING: SFR USFLOW, LENGTH or WIDTH not exposed by '
                          'this MF6 build; the stream evaporation is NOT '
                          'capped by the reach inflow (sfr.evap_inflow_'
                          'fraction), and a nearly dry reach can stop the '
                          'solver converging (analysis 8.23).')
                else:
                    n = self.nreaches
                    self.sfr_wl = (np.array(_len, dtype=float)[:n]
                                   * np.array(_wid, dtype=float)[:n])
                    print('coupler: stream evaporation capped at %g x each '
                          'reach\'s inflow (sfr.evap_inflow_fraction)'
                          % self.sfr_evap_frac)
        if self.nlakes:
            self.p_lak_evap, self.addr_lak_evap = self._bind_first(
                api, [('EVAPORATION', f'{name}/LAK'),
                      ('EVAPORATION', f'{name}/LAK-1')],
                'LAK evaporation', required=False)
            if self.p_lak_evap is None:
                print('WARNING: LAK EVAPORATION array not exposed by this MF6 '
                      'build; open-water evaporation will NOT be applied to '
                      'the ponds.')
            self.p_lak_simevap, _a = self._bind_first(
                api, [('EVAP', f'{name}/LAK'), ('EVAP', f'{name}/LAK-1')],
                'LAK simulated evaporation', required=False)
            # the stage itself -- XNEWPAK; STAGE is only the input of a
            # constant-stage lake (gwf-lak.f90)
            self.p_lak_stage, _a = self._bind_first(
                api, [('XNEWPAK', f'{name}/LAK'), ('XNEWPAK', f'{name}/LAK-1')],
                'LAK stage', required=False)
            if self.p_lak_stage is not None \
                    and self.p_lak_stage.size < self.nlakes:
                self.p_lak_stage = None
            # WP4.4: the MM runoff each pond captures, a volumetric rate
            # [m3/d], written with the other API inputs after prepare_solve
            self.p_lak_runoff, self.addr_lak_runoff = self._bind_first(
                api, [('RUNOFF', f'{name}/LAK'), ('RUNOFF', f'{name}/LAK-1')],
                'LAK runoff', required=False)
            if self.p_lak_runoff is None:
                print('WARNING: LAK RUNOFF array not exposed by this MF6 '
                      'build; the rain on the ponds and the runoff they '
                      'capture will NOT reach them.')
            elif self.p_lak_runoff.size < self.nlakes:
                print('WARNING: LAK RUNOFF has %d entries but there are %d '
                      'lakes; pond runoff disabled.'
                      % (self.p_lak_runoff.size, self.nlakes))
                self.p_lak_runoff = None
        if (self.p_sfr_simevap is None and self.p_lak_simevap is None
                and (self.nreaches or self.nlakes)):
            print('WARNING: neither SFR SIMEVAP nor LAK EVAP is exposed; '
                  'open-water evaporation will be APPLIED but will not appear '
                  'in the water balance (Eow reported as 0).')
        if self.p_gwd is None and self.p_drnseep is None:
            print('WARNING: neither the UZF groundwater-discharge array nor a '
                  'DRN_SEEP package was found; exfiltration into the soil will '
                  'be treated as 0.\n         (UZF variables seen: %s)'
                  % ', '.join(self.available_vars(api, 'UZF')[:12]))
        # Rejected infiltration: percolation MARMITES applied that UZF refused
        # (cell saturated, or rate above VKS). It is currently NOT returned to
        # the soil model, so it is tracked and reported to keep the coupled
        # water balance honest -- see the Phase-4 report.
        self.p_rejinf, self.addr_rejinf = self._bind_first(
            api, [('REJINF', f'{name}/UZF'), ('REJINF', f'{name}/UZF-1')],
            'UZF rejected infiltration', required=False)
        # WP2 step 4: UZF's PET demand. MF6 resets the land object's PET from
        # PETMAX on EVERY solve iteration (gwf-uzf.f90 uzf_solve), and uzf_ad
        # sets PET, GWPET and PETMAX from the period input PET_PVAR inside
        # prepare_solve -- so PETMAX is the operative array, written after
        # prepare_solve like FINF. Writing PET alone would be undone at the
        # next iteration: the FINF/SINF trap again. The objects below get the
        # demand the land object left, from MF6 itself (setbelowpet).
        self.p_petmax, self.addr_petmax = self._bind_first(
            api, [('PETMAX', f'{name}/UZF'), ('PETMAX', f'{name}/UZF-1')],
            'UZF PET demand (PETMAX)', min_size=self.ncell)
        self.p_pet, _ = self._bind_first(
            api, [('PET', f'{name}/UZF'), ('PET', f'{name}/UZF-1')],
            'UZF PET', min_size=self.ncell)
        # WP2 step 5: its ACTUAL ET, m3/d per UZF object (uzf_cq)
        self.p_uzet, self.addr_uzet = self._bind_first(
            api, [('UZET', f'{name}/UZF'), ('UZET', f'{name}/UZF-1')],
            'UZF actual ET (UZET)', required=False)
        print('coupler: UZF PET bound to %s; actual ET read from %s'
              % (self.addr_petmax, self.addr_uzet or '(not exposed)'))
        # MF6's OWN TIME STEP, with its retry. A step the solver cannot
        # converge is retried with a smaller one (ATS, dtfailadj) -- but only
        # inside MF6's own do_time_step (Mf6DoTimestep, mf6core.f90). Driven
        # through prepare_solve / solve / finalize_solve, as the coupler was,
        # a failed step is reported ("did not converge") and ACCEPTED: on
        # 2026-09-27 La Mata's stream cells oscillated by 14-69 m per
        # iteration on some wet days and those days went into the record.
        # do_time_step re-applies the period input on every try -- uzf_ad
        # copies SINF_PVAR to FINF/SINF and PET_PVAR to PET/PETMAX
        # (gwf-uzf.f90) -- so the exchanged values are written THERE too;
        # the EVT curves, SFR INFLOW/EVAP and LAK RUNOFF/EVAPORATION are
        # re-read only from a PERIOD block, which the build writes for
        # period 1 alone.
        self.p_sinf_pvar, _a = self._bind_first(
            api, [('SINF_PVAR', f'{name}/UZF'), ('SINF_PVAR', f'{name}/UZF-1')],
            'UZF period infiltration (SINF_PVAR)', required=False,
            min_size=self.ncell)
        self.p_pet_pvar, _a = self._bind_first(
            api, [('PET_PVAR', f'{name}/UZF'), ('PET_PVAR', f'{name}/UZF-1')],
            'UZF period PET (PET_PVAR)', required=False, min_size=self.ncell)
        self.p_itertot, _a = self._bind_first(
            api, [('ITERTOT_TIMESTEP', 'SLN_1')],
            'solver linear iterations', required=False)
        self.mf6_step = (self.p_sinf_pvar is not None
                         and (self.p_pet_pvar is not None
                              or self.p_petmax is None)
                         and hasattr(api, 'do_time_step'))
        if self.mf6_step:
            print('coupler: MF6 advances with its own time step '
                  '(do_time_step): a step it cannot solve is retried with a '
                  'smaller one (ATS)')
        else:
            print('WARNING: the UZF period-input arrays are not reachable, so '
                  'MF6 is stepped through prepare_solve/solve: a step it '
                  'cannot solve is NOT retried, it is accepted unconverged.')
        # each UZF object's mean water content, read once at the end so a
        # periodic spin-up can start the next cycle from it
        self.p_wcnew, _a = self._bind_first(
            api, [('WCNEW', f'{name}/UZF'), ('WCNEW', f'{name}/UZF-1')],
            'UZF water content (WCNEW)', required=False)
        # (The ETg wells' rates were bound here until 2026-10-07: WEL/Q only.
        # MF6 6.7 also exposes WEL/BOUND, and writing it corrupted MF6's heap
        # and killed the process in the next prepare_time_step with no
        # message -- bisected on La Mata, tests/diagnose_coupling.py. Keep
        # that in mind when WEL comes back for real pumping.)
        # node mapping: X is over reduced nodes; NODEUSER maps reduced -> user
        # DIS : user node = k*nrow*ncol + i*ncol + j
        # DISV: user node = k*ncpl + icell2d   (icell2d == the cell list 'node')
        if self.grid_kind == 'vertex':
            icell2d = np.array([c[3] for c in self.ctx.cells], dtype=int)
            user_nodes = self.k_arr * self.ncpl + icell2d
        else:
            nrc = self.mf6b.nrow * self.mf6b.ncol
            user_nodes = self.k_arr * nrc + self.i_arr * self.mf6b.ncol + self.j_arr
        try:
            nodeuser = np.asarray(api.get_value(self._addr(api, 'NODEUSER', f'{name}/DIS')))
            if nodeuser.size == self.p_x.size and nodeuser.size > 0:
                lut = {int(u) - 1: r for r, u in enumerate(nodeuser)}
                self.x_index = np.array([lut[int(u)] for u in user_nodes], dtype=int)
            else:
                self.x_index = user_nodes
        except Exception:
            # unreduced grid (idomain all >0) or API without NODEUSER
            self.x_index = user_nodes

    # ---------------- MM <-> pointer exchange -------------------------- #

    def _check_sizes(self):
        """Validate pointer sizes once, before anything is written.

        Writing past the end of an MF6 array corrupts its memory manager and
        the process dies during the next solve with no message at all -- the
        failure mode is silent, so it is checked up front.
        """
        if self.p_x.size < 1:
            raise CouplingError('MF6 head array X is empty')
        for _p, _nm in ((self.p_finf, 'FINF'), (self.p_sinf, 'SINF')):
            if _p.size < self.ncell:
                raise CouplingError(
                    'UZF %s has %d entries but MARMITES has %d surface cells: '
                    'the bound address is wrong or the UZF package was built '
                    'differently.' % (_nm, _p.size, self.ncell))
        if self.p_gwd is not None and self.p_gwd.size < self.ncell:
            print('WARNING: UZF discharge array is shorter (%d) than the cell '
                  'count (%d); exfiltration disabled.' % (self.p_gwd.size, self.ncell))
            self.p_gwd = None
        if int(self.x_index.max()) >= self.p_x.size:
            raise CouplingError(
                'node mapping points past the head array (max index %d, X size %d). '
                'The grid is reduced by idomain and NODEUSER could not be read.'
                % (int(self.x_index.max()), self.p_x.size))

    # ------------------------- solution validity ----------------------- #

    @staticmethod
    def listing_given_up(txt):
        """Stress periods of ``mfsim.lst`` in which MF6 gave up on a step.

        A failed try prints "Solution N did not converge for stress period
        P ..." and, when ATS retries it, "Failed solution ... will be
        retried" right after. Only a failure with no retry is final.
        """
        lines = txt.splitlines()
        rx = re.compile(r'^\s*Solution \d+ did not converge for stress '
                        r'period (\d+)')
        out = set()
        for i, line in enumerate(lines):
            m = rx.match(line)
            if m and not any('will be retried' in ln
                             for ln in lines[i + 1:i + 6]):
                out.add(int(m.group(1)))
        return out

    def check_solution(self, max_discrepancy=1.0, raise_on_fail=True):
        """Verify MF6 actually solved the problem. Call after the run.

        Heads that look plausible and iteration counts that look modest are
        NOT evidence of a solution. On La Mata a single non-converged stress
        period discharged 9.6e6 m3 -- eight times the recharge of the entire
        simulation -- through the seepage drain, leaving a cumulative budget
        discrepancy of -132% while the head field still looked ordinary. Both
        facts were sitting in the MF6 listing files the whole time.

        Returns a dict; raises CouplingError on failure unless told not to.
        """
        rep = {'nonconverged_api': int(self.n_nonconverged),
               'nonconverged_lst': 0, 'discrepancy': None, 'ok': True,
               'substeps': int(self.substeps), 'messages': []}
        ws = self.sim_ws
        if ws and os.path.isdir(ws):
            fn = os.path.join(ws, 'mfsim.lst')
            if os.path.exists(fn):
                with open(fn, 'r', errors='replace') as fh:
                    txt = fh.read()
                # periods MF6 GAVE UP on -- not every failed try: with
                # ATS each try that is retried prints "did not converge"
                # too (twice), and counting those reported 462 for a run
                # whose steps failed for good in 10 periods (2026-09-27)
                rep['nonconverged_lst'] = len(self.listing_given_up(txt))
            for cand in os.listdir(ws):
                if not cand.endswith('.lst') or cand == 'mfsim.lst':
                    continue
                with open(os.path.join(ws, cand), 'r', errors='replace') as fh:
                    lst = fh.read()
                vals = re.findall(r'PERCENT DISCREPANCY\s*=\s*(-?[\d.]+)', lst)
                if vals:
                    # the cumulative column is the first of each pair
                    rep['discrepancy'] = max(abs(float(v)) for v in vals[-2:])
                break
        if rep['nonconverged_api'] or rep['nonconverged_lst']:
            rep['ok'] = False
            rep['messages'].append(
                '%d stress period(s) did not converge (API count %d, listing '
                'count %d). A non-converged step is not a result: MF6 moves on '
                'with whatever the solver last held, and the water it invents '
                'or destroys stays in the record.'
                % (max(rep['nonconverged_api'], rep['nonconverged_lst']),
                   rep['nonconverged_api'], rep['nonconverged_lst']))
        d = rep['discrepancy']
        if d is not None and d > max_discrepancy:
            rep['ok'] = False
            rep['messages'].append(
                'cumulative mass-balance discrepancy is %.2f%% (limit %.2f%%). '
                'The model is not conserving water.' % (d, max_discrepancy))

        # Infiltration-fidelity check. The recharge coupling is only real if the
        # daily percolation MARMITES writes actually reaches UZF. Binding the
        # wrong UZF variable (FINF vs SINF) left UZF applying a CONSTANT
        # perc_user while MARMITES varied day to day -- the water table drained
        # and nothing flagged it. So: if the written percolation varies but the
        # UZF budget INFILTRATION is (nearly) constant, the coupling is broken.
        infil = self._uzf_infiltration_cv(ws)
        rep['uzf_infiltration'] = infil
        if infil is not None:
            written_cv, uzf_cv, uzf_vals = infil
            if written_cv > 0.2 and uzf_cv < 0.02:
                rep['ok'] = False
                rep['messages'].append(
                    'UZF INFILTRATION is nearly constant (CV %.3f over %s m3/d) '
                    'while the written percolation varies (CV %.2f). The daily '
                    'recharge is NOT reaching UZF -- check the SINF/FINF binding.'
                    % (uzf_cv, [round(v, 1) for v in uzf_vals[:4]], written_cv))
        # The UZF PACKAGE's own budget. The check above reads the budget's
        # INFILTRATION line, which follows SINF -- the REPORTED array -- and
        # so passed while UZF routed a constant FINF: 92.5 % of the UZF
        # budget unaccounted (2026-09-23), with the GWF budget at 0.04 %.
        # Routed and reported water disagreeing is exactly what this sees.
        uzf_d, gwet = self._uzf_budget_checks(ws)
        rep['uzf_discrepancy'], rep['uzf_gwet'] = uzf_d, gwet
        if uzf_d is not None and uzf_d > max_discrepancy:
            rep['ok'] = False
            rep['messages'].append(
                'the UZF package budget is %.2f%% out (limit %.2f%%): what UZF '
                'routes is not what it reports -- check which UZF arrays the '
                'coupler writes (FINF and SINF, PETMAX and PET).'
                % (uzf_d, max_discrepancy))
        # WP2 2.6, the GWET guard: groundwater ET belongs to MARMITES (the
        # EVT packages), so UZF must simulate none -- no configuration key
        # can ask for linear_gwet or square_gwet, and this catches a
        # regression.
        # TOTAL ET NEVER ABOVE PET. uzf_demand makes it hold by construction;
        # a cell-period above it is a fault in the chain, not a result.
        n_over = int(getattr(self, 'n_overdraw', 0))
        rep['et_above_pet'] = n_over
        if n_over:
            rep['ok'] = False
            rep['messages'].append(
                'ETsoil + ETuzf + ETg exceeded PE + PT in %d cell-period(s), '
                'by up to %.3g mm/d. Total ET can never be above PET.'
                % (n_over, float(getattr(self, 'max_overdraw', 0.0))))
        if gwet:
            rep['ok'] = False
            rep['messages'].append(
                'UZF simulated groundwater ET (UZF-GWET %.4g m3): it is '
                'MARMITES\', taken by the EVT packages, and would be counted '
                'twice. The UZF package must be built without linear_gwet / '
                'square_gwet.'
                % gwet)
        if not rep['ok']:
            msg = 'MF6 solution is not valid:\n  - ' + '\n  - '.join(rep['messages'])
            if raise_on_fail:
                raise CouplingError(msg)
            print('\nWARNING: ' + msg)
        elif d is not None:
            print('solution check: converged, cumulative discrepancy %.4f%%, '
                  '%d extra ATS sub-step(s)' % (d, rep['substeps']))
        return rep

    @staticmethod
    def _uzf_budget_checks(ws):
        """``(uzf_discrepancy, uzf_gwet)`` from the model listing file.

        uzf_discrepancy: the |cumulative PERCENT DISCREPANCY| of the LAST
        'UZF BUDGET FOR ENTIRE MODEL' table (None when there is none).
        uzf_gwet: the largest cumulative UZF-GWET / GWET volume reported,
        0.0 when the term is absent -- as it must be.
        """
        if not ws or not os.path.isdir(ws):
            return None, 0.0
        for cand in os.listdir(ws):
            if not cand.endswith('.lst') or cand == 'mfsim.lst':
                continue
            with open(os.path.join(ws, cand), 'r', errors='replace') as fh:
                lst = fh.read()
            d = None
            blocks = lst.split('UZF BUDGET FOR ENTIRE MODEL')
            if len(blocks) > 1:
                last = blocks[-1].split('VOLUME BUDGET FOR ENTIRE MODEL')[0]
                v = re.findall(r'PERCENT DISCREPANCY\s*=\s*(-?[\d.E+-]+)', last)
                if v:
                    d = abs(float(v[0]))           # the cumulative column
            g = [abs(float(x)) for x in re.findall(
                r'(?:UZF-GWET|\bGWET)\s*=\s*(-?[\d.E+-]+)', lst)]
            return d, (max(g) if g else 0.0)
        return None, 0.0

    def _uzf_infiltration_cv(self, ws):
        """Coefficient of variation of written percolation vs UZF INFILTRATION.

        Returns (written_cv, uzf_cv, [sampled uzf infiltration m3/d]) or None if
        the UZF budget cannot be read. A near-zero uzf_cv against a varying
        written_cv means the daily recharge is not reaching UZF.
        """
        if self.perc_hist is None or not ws:
            return None
        uzf_cbc = None
        for cand in os.listdir(ws):
            if cand.endswith('.uzf.cbc'):
                uzf_cbc = os.path.join(ws, cand)
                break
        if uzf_cbc is None:
            return None
        try:
            import flopy
            cbc = flopy.utils.CellBudgetFile(uzf_cbc)
            kk = [k for k in cbc.get_kstpkper() if k[1] >= 1]
            if len(kk) < 3:
                return None
            idx = np.linspace(0, len(kk) - 1, min(10, len(kk))).round().astype(int)
            uzf_vals, wrote = [], []
            for i in sorted(set(idx)):
                k = kk[i]
                d = cbc.get_data(kstpkper=k, text='INFILTRATION')
                if not d:
                    continue
                a = d[0]
                v = a['q'].sum() if hasattr(a, 'dtype') and a.dtype.names else float(np.sum(a))
                uzf_vals.append(abs(float(v)))
                sp = k[1] - 1
                if sp < self.perc_hist.shape[0]:
                    wrote.append(float(self.perc_hist[sp].sum()))
            if len(uzf_vals) < 3:
                return None
            uzf_vals = np.asarray(uzf_vals)
            wrote = np.asarray(wrote)
            cv = lambda x: float(np.std(x) / np.mean(x)) if np.mean(x) > 0 else 0.0
            return cv(wrote), cv(uzf_vals), list(uzf_vals)
        except Exception:                          # pragma: no cover
            return None

    # The rates MF6 reports per TIME STEP that the coupler reads once a
    # stress period is over. When ATS splits the period they are those of
    # its LAST sub-step: on 2026-09-27 a day of 35 sub-steps was booked with
    # the UZF ET, seepage and rejected infiltration of its final 20 minutes,
    # and UZF's wave tolerance -- a depth per step -- divided by that
    # sub-step's length read as ET above PET in 42 cell-periods.
    SUBSTEP_RATES = ('p_uzet', 'p_drnseep', 'p_gwd', 'p_rejinf',
                     'p_sfr_simevap', 'p_lak_simevap',
                     'p_evt_eg_sim', 'p_evt_tg_sim')

    def _rate(self, name):
        """A reported rate for the period just advanced: its mean over the
        sub-steps, weighted by their length, when MF6 split the period;
        MF6's own array otherwise."""
        m = getattr(self, '_period_rates', None)
        if m is not None and name in m:
            return m[name]
        return getattr(self, name, None)

    def _rates_add(self, acc, dt):
        """Add one time step's reported rates to ``acc``, times its length."""
        dt = float(dt) if dt is not None and np.isfinite(dt) and dt > 0 \
            else 1.0
        for name in self.SUBSTEP_RATES:
            ptr = getattr(self, name, None)
            if ptr is None:
                continue
            v = np.asarray(ptr, dtype=float).ravel() * dt
            if name in acc and acc[name].shape == v.shape:
                acc[name] += v
            else:
                acc[name] = v
        acc['__dt'] = acc.get('__dt', 0.0) + dt

    def _rates_close(self, acc, nsub):
        """End of a period: keep the sub-step means (only when it WAS split
        -- a single step's arrays are MF6's own) and how many steps."""
        tot = acc.pop('__dt', 0.0)
        self._period_rates = ({k: v / tot for k, v in acc.items()}
                              if nsub > 1 and tot > 0.0 else None)
        self._period_nsub = int(nsub)

    @staticmethod
    def _now(api):
        try:
            return float(api.get_current_time())
        except Exception:                               # noqa: BLE001
            return None

    def _read_heads_exf(self):
        """Heads [m] and groundwater exfiltration [mm/d, + into the soil].

        The seepage comes from whichever mechanism the model was built with:
        UZF SIMULATE_GWSEEP (GWD, positive out of the aquifer) or the
        land-surface DRN_SEEP package (SIMVALS, negative out of the aquifer).
        """
        heads = np.asarray(self.p_x, dtype=float)[self.x_index]
        if self.p_drnseep is not None:
            sim = np.asarray(self._rate('p_drnseep'), dtype=float).ravel()
            q = np.zeros(self.ncell)
            has = self.drnseep_idx >= 0
            idx = self.drnseep_idx[has]
            if idx.size and int(idx.max()) < sim.size:
                # drains remove water -> SIMVALS <= 0; exfiltration is positive
                q[has] = np.maximum(-sim[idx], 0.0)
            return heads, q / self.area * self.conv_fact
        if self.p_gwd is None:
            return heads, np.zeros(self.ncell)
        gwd = np.asarray(self._rate('p_gwd'),
                         dtype=float).ravel()[:self.ncell]  # m3/d, + to surface
        exf = gwd / self.area * self.conv_fact                          # mm/d, + into soil
        return heads, exf

    def _read_rejinf(self):
        """Rejected infiltration [mm/d, positive] per cell.

        Percolation MARMITES applied that UZF could not accept. It is fed
        back into the soil model's SURFACE store on the next step, where the
        existing capacity logic spills any excess to runoff -- otherwise this
        water would silently leave the coupled water balance.
        """
        if self.p_rejinf is None:
            return np.zeros(self.ncell)
        rej = np.abs(np.asarray(self._rate('p_rejinf'),
                                dtype=float).ravel()[:self.ncell])
        return rej / self.area * self.conv_fact

    def _write_runoff(self, ro):
        """Deliver the MARMITES runoff of this SP to the streams and ponds.

        ``ro`` is per cell in mm/d. A channel cell's runoff goes to its
        reach (SFR INFLOW), a pond footprint's to its lake (LAK RUNOFF,
        WP4.4) -- both volumetric rates [m3/d]. Since WP4.6 it includes the
        rain on the cell's open fraction, which MMsoil hands over as runoff.
        With the CRR cascade (WP5) ``ro`` is what the cascade brought to each
        stream and pond cell -- its own runoff and its upslope cells' -- and
        runoff generated on a soil cell reaches SFR and LAK that way;
        without it, a soil cell's runoff is not delivered at all.
        """
        q = np.maximum(np.asarray(ro, dtype=float), 0.0) / self.conv_fact \
            * self.area                                              # m3/d
        if self.p_sfr_inflow is not None and self.sfr_reach_idx is not None:
            inflow = np.zeros(self.p_sfr_inflow.size)
            has = self.sfr_reach_idx >= 0
            np.add.at(inflow, self.sfr_reach_idx[has], q[has])
            self.p_sfr_inflow[:] = inflow
        if getattr(self, 'p_lak_runoff', None) is not None and self.nlakes:
            rn = np.zeros(self.p_lak_runoff.size)
            for L, (ks, _w) in enumerate(self.lak_cells):
                if L < rn.size and ks.size:
                    rn[L] = float(q[ks].sum())
            self.p_lak_runoff[:] = rn

    def _delivers_runoff(self):
        """Is there a stream or a pond bound to take MARMITES runoff?"""
        return ((self.p_sfr_inflow is not None
                 or getattr(self, 'p_lak_runoff', None) is not None)
                and getattr(self, '_iRo', None) is not None)

    def _write_openwater_evap(self, n):
        """Apply the open-water evaporation of SP ``n`` to SFR and LAK (WP1d).

        Both packages take a rate per unit of wetted area, in model length
        units per day; ``Eo`` is mm/d, so it is divided by ``conv_fact``. MF6
        only removes what the reach or lake actually holds, which is what the
        MMsoil version did too (``Eow`` was capped by the surface store).

        Returns the rate written, for the water-balance record.
        """
        if self.Eo_zonesSP is None:
            return None
        col = min(int(n), self.Eo_zonesSP.shape[1] - 1)
        eo = self.Eo_zonesSP[:, col] / self.conv_fact             # m/d
        wrote = None
        if self.p_sfr_evap is not None and self.nreaches:
            v = self._cap_sfr_evap(eo[np.minimum(self.sfr_evap_zone,
                                                 len(eo) - 1)])
            self.p_sfr_evap[:min(self.p_sfr_evap.size, self.nreaches)] = \
                v[:min(self.p_sfr_evap.size, self.nreaches)]
            wrote = float(v.mean())
        if self.p_lak_evap is not None and self.nlakes:
            v = eo[np.minimum(self.lak_evap_zone, len(eo) - 1)]
            self.p_lak_evap[:min(self.p_lak_evap.size, self.nlakes)] = \
                v[:min(self.p_lak_evap.size, self.nlakes)]
            wrote = float(v.mean()) if wrote is None else wrote
        return wrote

    def _cap_sfr_evap(self, eo):
        """Each reach's Eo [m/d], capped at ``sfr_evap_frac`` of its inflow.

        EVAP <= f x qin / (w L), qin = this step's runoff into the reach
        (INFLOW, written just before by _write_runoff) + the upstream flow
        and pond outflow MF6 left at the end of the previous step
        (_keep_sfr_inflow). Its own groundwater discharge is not counted:
        a reach the aquifer feeds does not flip-flop, and one it stops
        feeding during the step must not be left evaporating it.

        Why (analysis §8.23): MF6 takes a reach's evaporation from the depth
        its previous solve stored. With less inflow than E0 x w x L and the
        head below the bed, the reach alternates between evaporating all of
        it and leaking all of it, sfr_fn's 1e-4 m perturbation lands in the
        other state, and GWF gets a derivative of qin / 1e-4 m where the true
        one is ~0. Kept below the inflow, the evaporation cannot empty the
        reach and the switch is never reached.
        """
        eo = np.asarray(eo, dtype=float)
        if self.sfr_evap_frac <= 0.0:
            # 0: the streams do not evaporate -- no area needed for that
            return np.zeros_like(eo)
        if self.sfr_wl is None:
            return eo
        n = min(eo.size, self.sfr_wl.size)
        qin = np.zeros(n)
        if self.p_sfr_inflow is not None:
            qin += np.asarray(self.p_sfr_inflow, dtype=float)[:n]
        if self.sfr_up_prev is not None:
            qin += self.sfr_up_prev[:n]
        cap = (self.sfr_evap_frac * np.maximum(qin, 0.0)
               / np.maximum(self.sfr_wl[:n], 1e-12))
        out = eo.copy()
        out[:n] = np.minimum(eo[:n], cap)
        self.sfr_evap_capped += int(np.count_nonzero(out[:n] < eo[:n]))
        self.sfr_evap_writes += n
        return out

    def _keep_sfr_inflow(self):
        """Copy what flowed into each reach at the end of the step just
        solved, besides its runoff: USFLOW + QFROMMVR [m3/d].

        MF6 zeroes USFLOW when it advances the next step (sfr_ad, in
        prepare_solve -- which the older prepare_solve route runs BEFORE the
        coupler writes), so the cap of _cap_sfr_evap reads this copy.
        """
        if self.sfr_wl is None:
            return
        n = self.sfr_wl.size
        up = np.array(self.p_sfr_usflow, dtype=float)[:n]
        if self.p_sfr_qfrommvr is not None:
            up = up + np.asarray(self.p_sfr_qfrommvr, dtype=float)[:n]
        self.sfr_up_prev = up

    def _read_openwater_evap(self):
        """What SFR and LAK actually evaporated this SP, per cell, in mm/d.

        MMsoil used to compute this itself as ``Eow``, from its own surface
        store; WP1d moved both the water and the evaporation to MODFLOW, and
        this brings the number back so ``iEow`` is a measured flux again
        rather than a structural zero.

        SFR ``SIMEVAP`` and LAK ``EVAP`` are VOLUMETRIC rates (m3/d) capped by
        what the reach or lake holds -- a dry reach evaporates nothing, which
        is the same limit the old surface store imposed. Dividing by the CELL
        area makes it a depth over the cell, which is the unit every other
        MARMITES flux uses.
        """
        sfr, lak = self._read_openwater_evap_split()
        return sfr + lak

    def _read_openwater_evap_split(self):
        """The same, KEPT APART by package: (SFR, LAK), each per cell in mm/d.

        WP2 row 2. The stream loses water in transit and a pond loses its
        storage, and a model that only knows their sum cannot say which --
        so each is read from its own package and carried in its own column.
        """
        sfr = np.zeros(self.ncell, dtype=float)
        lak = np.zeros(self.ncell, dtype=float)
        if self.p_sfr_simevap is not None and self.nreaches:
            se = np.asarray(self._rate('p_sfr_simevap'), dtype=float)
            has = self.sfr_reach_idx >= 0
            idx = self.sfr_reach_idx[has]
            ok = idx < se.size
            np.add.at(sfr, np.nonzero(has)[0][ok], np.abs(se[idx[ok]]))
        vol = self._lak_evap_volumes()
        foot = getattr(self, 'lak_cells', None)
        if vol is not None:
            for L in range(vol.size):
                if foot is not None and L < len(foot) and foot[L][0].size:
                    # over the pond's footprint, by the pond area each cell
                    # holds (WP4.6) -- the host cell alone would carry a
                    # 2000 m2 pond's evaporation on 25 m2
                    ks, w = foot[L]
                    np.add.at(lak, ks, float(vol[L]) * w)
                    continue
                n = int(self.lak_cell_idx[L])
                if n >= 0:
                    lak[n] += float(vol[L])
        # MF6 reports a removal as a positive rate here; carry it as a
        # positive loss, like every other MARMITES flux.
        return (sfr / self.area * self.conv_fact,                  # mm/d
                lak / self.area * self.conv_fact)

    def _lak_evap_volumes(self):
        """Each pond's actual evaporation this SP [m3/d], positive; None
        when there are no ponds or LAK does not expose it."""
        if self.p_lak_simevap is None or not self.nlakes:
            return None
        le = np.abs(np.asarray(self._rate('p_lak_simevap'), dtype=float))
        return le[:min(self.nlakes, le.size)].copy()

    def _openwater_evap_into(self, mm_cells):
        """Put this SP's open-water evaporation into ``iEow`` of the vector,
        and each package's share into ``iEow_sfr`` and ``iEow_lak``.

        Deliberately does NOT reduce ``iRo``. Per cell it could not: a reach
        carries water from upstream, so its evaporation can exceed the runoff
        that cell generated, and subtracting would drive Ro negative. Eow is a
        loss from the CHANNEL, downstream of the runoff -- which is what the
        water-budget figures now say.

        The TOTAL ET follows. MMsoil built ``iETtot`` with its own Eow, which
        is zero since the surface store went to MODFLOW, so the measured
        value is added here -- or the five sources would sum to four.
        """
        if self._iEow is None:
            return None
        sfr, lak = self._read_openwater_evap_split()
        eow = sfr + lak
        ix = self.ctx.index if getattr(self, 'ctx', None) is not None else {}
        if 'iETtot' in ix:
            mm_cells[:, ix['iETtot']] += eow - mm_cells[:, self._iEow]
        mm_cells[:, self._iEow] = eow
        if 'iEow_sfr' in ix:
            mm_cells[:, ix['iEow_sfr']] = sfr
        if 'iEow_lak' in ix:
            mm_cells[:, ix['iEow_lak']] = lak
        self._eow_split = (sfr, lak)
        return eow

    # mm/d: what the per-cell PET check tolerates -- floating rounding only
    ET_TOL = 1e-6

    @staticmethod
    def uzf_demand(petuzf, etg_cap):
        """UZF's PET demand [m/d]: what the soil left, LESS what groundwater
        ET may take today.

        TOTAL ET CAN NEVER EXCEED PET. MMsoil hands on PETuzf = PE + PT the
        soil did not use, and sizes groundwater ET seeing UZF's ET of the
        PREVIOUS period -- the only one known before MF6 solves. Writing all
        of PETuzf to UZF let this period's UZF ET and ETg together exceed
        it: 125,030 cell-periods, up to 2.05 mm/d, on 2026-09-24. With UZF
        capped at PETuzf - ``etg_cap``, the most the day's EVT curves can
        take (their RATE: EVT only falls from it as the head does),
        ETsoil + ETuzf + ETg <= PE + PT holds by construction, since MF6
        removes at most PETMAX from a column (setbelowpet hands down only
        the unmet part).

        UZF keeps its place in the chain: ETg was sized on PETuzf LESS last
        period's ETuzf, so the cap leaves UZF at least what it took then.
        """
        d = np.asarray(petuzf, dtype=float) - np.asarray(etg_cap, dtype=float)
        return np.maximum(d, 0.0)

    def _write_fluxes(self, perc, petuzf=None, etg_cap=None):
        # etg_cap: the most EVT can take today, which UZF's demand must
        # leave room for -- ETg itself is not known until MF6 solves
        # both UZF arrays, as MF6's setdatafinf sets them: the land cells
        # carry the percolation, the objects below them nothing
        # ... and the period input uzf_ad copies them from on every try of
        # MF6's own time step (see _bind)
        for _p in (self.p_finf, self.p_sinf,
                   getattr(self, 'p_sinf_pvar', None)):
            if _p is None:
                continue
            _p[:self.ncell] = perc                                    # m/d
            if _p.shape[0] > self.ncell:
                _p[self.ncell:] = 0.0
        # WP2: the deep unsaturated zone's PET demand [m/d] on the land
        # objects, after groundwater ET; MF6 passes what they leave to the
        # objects below
        if petuzf is not None and self.p_petmax is not None:
            dem = self.uzf_demand(petuzf, 0.0 if etg_cap is None else etg_cap)
            for _p in (self.p_petmax, self.p_pet,
                       getattr(self, 'p_pet_pvar', None)):
                if _p is not None:
                    _p[:self.ncell] = dem
            self._petmax_written = np.array(dem, dtype=float)

    def _read_etuzf(self):
        """UZF's ACTUAL ET per land cell [mm/d]: UZET (m3/d per object)
        summed over each column's objects, on the cell's own area."""
        if self.p_uzet is None:
            return np.zeros(self.ncell)
        q = np.abs(np.asarray(self._rate('p_uzet'), dtype=float).ravel())
        if self._col_obj.size and self._col_obj.max() >= q.size:
            return np.zeros(self.ncell)
        tot = np.bincount(self._col_own, weights=q[self._col_obj],
                          minlength=self.ncell)
        return tot / self.area * self.conv_fact

    # MF6's own tolerance on UZF ET [m of water per m of unsaturated zone].
    # UNSAT_ETWC books ET as the change in a column's storage, and after
    # taking it merges wave pairs whose water contents differ by less than
    # DEM6 = 1e-6 (UzfCellGroup.f90), which moves up to 1e-6 x the depth
    # between them. Replaying the one-year run of 2026-09-24 against its own
    # MF6 output: UZF took more than the PETMAX written in 612 cell-periods,
    # all by less than 0.73 x 1e-6 m per metre of unsaturated zone.
    UZF_WAVE_TOL = 1e-6

    def _split_etuzf(self, et, n):
        """UZF's ET as ``(booked, residual)`` [mm/d] per cell.

        TOTAL ET CAN NEVER EXCEED PET, and what the coupler wrote to PETMAX
        is what PET left for UZF. What MF6 removed beyond it, WITHIN its wave
        tolerance (UZF_WAVE_TOL x the unsaturated thickness per step), is
        MF6's numerical storage loss, not evapotranspiration: it is booked
        apart (iETuzf_num), so UZF's balance still matches MF6's budget and
        ET stays within the demand. Beyond the tolerance nothing is split --
        it stays ET, and the PET check fails the run on it.

        The tolerance is a depth PER TIME STEP, so a period ATS split into
        n sub-steps may carry n of them (its ET is the sub-steps' mean).
        """
        et = np.asarray(et, dtype=float)
        dem = getattr(self, '_petmax_written', None)
        top = self.top_cell
        if dem is None or top is None or self.heads_hist is None:
            return et, np.zeros_like(et)
        over = np.maximum(et - np.asarray(dem, dtype=float)
                          * self.conv_fact, 0.0)
        uz = np.maximum(top - np.asarray(self.heads_hist[n], dtype=float), 0.0)
        p_x = getattr(self, 'p_x', None)
        if p_x is not None and getattr(self, 'x_index', None) is not None:
            h_end = np.asarray(p_x, dtype=float)[self.x_index]
            uz = np.maximum(uz, top - h_end)
        perlen = float(self.perlen[n]) if n < len(self.perlen) else 1.0
        nsub = max(int(getattr(self, '_period_nsub', 1) or 1), 1)
        tol = (self.UZF_WAVE_TOL * uz * nsub / max(perlen, 1e-12)
               * self.conv_fact)
        resid = np.where((over > 0.0) & (over <= tol), over, 0.0)
        return et - resid, resid

    def _etuzf_into(self, mm_cells, et, n, resid=None):
        """This period's ACTUAL UZF ET into the per-cell vector, the total
        ET made whole, and the PET balance of the period recorded."""
        ix = self.ctx.index
        if 'iETuzf' in ix:
            mm_cells[:, ix['iETuzf']] = et
        if resid is not None and 'iETuzf_num' in ix:
            mm_cells[:, ix['iETuzf_num']] = resid
        if 'iETtot' in ix:
            mm_cells[:, ix['iETtot']] += et
        if not all(k in ix for k in ('iPT', 'iPE', 'iETsoil', 'iETg')):
            return
        demand = mm_cells[:, ix['iPT']] + mm_cells[:, ix['iPE']]
        used = mm_cells[:, ix['iETsoil']] + et + mm_cells[:, ix['iETg']]
        over = used - demand
        bad = over > self.ET_TOL
        self.n_overdraw += int(bad.sum())
        if bad.any():
            self.max_overdraw = max(self.max_overdraw, float(over[bad].max()))
        self.pet_unmet[n] = float(np.average(demand - used, weights=self.area))

    def _print_pet_balance(self, nper):
        """The run's PET balance, catchment, area-weighted [mm/yr]."""
        ix = self.ctx.index
        need = ('iPT', 'iPE', 'iETsoil', 'iETuzf', 'iETg')
        if self.wb_ts is None or not all(k in ix for k in need):
            return
        ts = np.asarray(self.wb_ts)[:nper]
        y = 365.0 / max(float(np.sum(self.perlen[:nper]) or nper), 1e-9)
        tot = {k: float(ts[:, ix[k]].sum()) * y for k in need}
        demand = tot['iPT'] + tot['iPE']
        used = tot['iETsoil'] + tot['iETuzf'] + tot['iETg']
        print('\nPET balance (catchment, mm/yr): demand PT+PE %.1f -> ETsoil '
              '%.1f, ETuzf %.1f, ETg %.1f; unmet %.1f (%.0f %%)'
              % (demand, tot['iETsoil'], tot['iETuzf'], tot['iETg'],
                 demand - used,
                 100.0 * (demand - used) / demand if demand > 0 else 0.0))
        if self.n_overdraw:
            print('      ERROR: ET above PET in %d cell-period(s), at most '
                  '%.3g mm/d -- this must never happen (check_solution '
                  'fails the run)' % (self.n_overdraw, self.max_overdraw))
        if getattr(self, 'n_resid', 0):
            print('      UZF: in %d cell-period(s) MF6 took up to %.3g mm/d '
                  'beyond the demand, within its own wave tolerance (1e-6 m '
                  'per metre of unsaturated zone): %.3g m3 in all, booked as '
                  'UZF numerical loss (iETuzf_num), not ET'
                  % (self.n_resid, self.resid_max, self.resid_m3))
        # the open water, each package on its own line (WP2 row 2)
        if (self.nreaches or self.nlakes) and 'iEow' in ix:
            ow = {k: float(ts[:, ix[k]].sum()) * y
                  for k in ('iEow', 'iEow_sfr', 'iEow_lak') if k in ix}
            print('      open water (mm/yr over the catchment): Eow %.2f = '
                  'streams (SFR) %.2f + ponds (LAK) %.2f'
                  % (ow.get('iEow', 0.0), ow.get('iEow_sfr', 0.0),
                     ow.get('iEow_lak', 0.0)))

    # ---------------- main drive --------------------------------------- #

    def _progress(self, n, nper, last=False):
        """Say where the run is, occasionally.

        A DETACHED RUN OF HOURS MUST SAY WHERE IT IS. La Mata's full record
        is 1949 daily stress periods and took ~4.7 h, and the log carried no
        stress-period marker at all -- so "is it at period 20 or 1900?" could
        only be answered by measuring the size of the .hds file. Nothing was
        wrong with the run; there was simply no way to tell.

        Reported every 5% and never more often than that, so the line cannot
        become the noise it exists to cut through. This is also where WP2.5b
        puts the PET balance, which is why it takes the whole stress period
        rather than just its number.
        """
        if nper <= 0:
            return
        if last and getattr(self, '_prog_from', None) == n + 1:
            return                 # the loop already reported this period
        now = time.time()
        t_last = getattr(self, '_t_last', None)
        t_last = self._t0 if t_last is None else t_last
        # every 5 %, AND whenever PROGRESS_EVERY_S has gone by without a
        # line: a run whose periods crawl (a solve failing and retried
        # through every outer iteration) said nothing for over an hour
        # before its first 5 % (2026-10-09)
        step = max(1, nper // 20)
        due = ((n > 0 and n % step == 0)
               or now - t_last >= PROGRESS_EVERY_S)
        if not last and not due:
            return
        done, left = n + 1, nper - (n + 1)
        # the time these periods took, not only the time left (user,
        # 2026-10-09)
        nblk = max(done - int(getattr(self, '_blk_from', 0)), 1)
        blk = now - t_last
        rate = (now - self._t0) / max(done, 1)
        eta = ('' if last or not left
               else ', ~%s left' % _hms(rate * left))
        self._t_last, self._blk_from = now, done
        print('   stress period %d/%d (%.0f%%): %d period(s) in %s (%s each),'
              ' %s elapsed%s%s'
              % (done, nper, 100.0 * done / nper, nblk, _hms(blk),
                 _hms(blk / nblk), _hms(now - self._t0), eta,
                 self._pet_window(n) if hasattr(self, '_pet_window')
                 else ''))

    def _pet_window(self, n):
        """WP2.5b: the PET balance since the last progress line, catchment
        means in mm/d -- so a starved or over-drawn chain shows up while the
        run goes, not only in the balance at the end."""
        ix = self.ctx.index
        need = ('iPT', 'iPE', 'iETsoil', 'iETuzf', 'iETg')
        ts = getattr(self, 'wb_ts', None)
        start = int(getattr(self, '_prog_from', 0))
        self._prog_from = n + 1
        if ts is None or not all(k in ix for k in need) or n < start:
            return ''
        m = {k: float(np.mean(np.asarray(ts)[start:n + 1, ix[k]]))
             for k in need}
        dem = m['iPT'] + m['iPE']
        used = m['iETsoil'] + m['iETuzf'] + m['iETg']
        return ('  |  PET %.2f mm/d -> soil %.2f, UZF %.2f, groundwater %.2f;'
                ' unmet %.0f%%'
                % (dem, m['iETsoil'], m['iETuzf'], m['iETg'],
                   100.0 * (dem - used) / dem if dem > 0 else 0.0))

    def run(self, api, on_sp=None):
        """Drive the coupled model with an initialized-able MF6 API object.

        api : modflowapi.ModflowApi or compatible (tests inject a fake).
        on_sp : optional callback(n, out_dict) after each stress period.
        """
        cMF = self.mf6b.cMF
        nper_mm = int(cMF.nper)
        self._t0 = time.time()
        self._t_last = self._t0           # the time of the last progress line
        self._blk_from = 0                # and the periods it had reported
        self._prog_from = 0
        self.heads_hist = np.zeros((nper_mm, self.ncell))
        self.exf_hist = np.zeros((nper_mm, self.ncell))
        self.perc_hist = np.zeros((nper_mm, self.ncell))
        self.etg_hist = np.zeros((nper_mm, self.ncell))
        self.rejinf_hist = np.zeros((nper_mm, self.ncell))
        self.etuzf_hist = np.zeros((nper_mm, self.ncell))
        # the soil column's water at the end of every period, per cell [mm
        # over the cell]: the soil-storage map series, on the head series'
        # days (user, 2026-10-10) -- a time mean of its change is ~0 once a
        # spin-up has converged
        self.ssoil_hist = np.zeros((nper_mm, self.ncell))
        self.pet_unmet = np.zeros(nper_mm)
        self.etuzf_prev = np.zeros(self.ncell)
        self.n_overdraw, self.max_overdraw = 0, 0.0
        self.n_resid, self.resid_max, self.resid_m3 = 0, 0.0, 0.0
        # a fresh MF6 session: no step solved yet, nothing flowed in
        self.sfr_up_prev = None
        self.sfr_evap_capped = self.sfr_evap_writes = 0
        self.outer_iters = np.zeros(nper_mm, dtype=int)
        # Water-budget aggregates (compact: the full per-cell/per-SP MM array
        # would be ~370 MB). wb_ts = catchment mean of every MM flux per SP;
        # wb_map = time mean of every MM flux per cell. Together these support
        # the budget time series, totals and maps without storing everything.
        nidx = len(self.ctx.index)
        nidx_s = len(self.ctx.index_S)
        nsl = int(self.ctx._nslmax)
        self._iRo = int(self.ctx.index['iRo'])
        self._iEow = self.ctx.index.get('iEow')
        self._iEow = None if self._iEow is None else int(self._iEow)
        self.runoff_hist = np.zeros((nper_mm, self.ncell))
        self.crr_ts = (None if self.crr is None else
                       np.zeros((nper_mm, len(self.CRR_TERMS))))
        self.evap_hist = np.zeros((nper_mm, self.ncell))
        # ... and by package (WP2 row 2); each pond's own loss in m3/d
        self.evap_sfr_hist = np.zeros((nper_mm, self.ncell))
        self.evap_lak_hist = np.zeros((nper_mm, self.ncell))
        self.lak_evap_hist = np.zeros((nper_mm, max(int(self.nlakes), 0)))
        self.wb_ts = np.zeros((nper_mm, nidx))
        self.wb_map = np.zeros((self.ncell, nidx))
        self.wb_ts_soil = np.zeros((nper_mm, nsl, nidx_s))
        self.wb_map_soil = np.zeros((self.ncell, nsl, nidx_s))
        # full per-SP MM flux vectors at the observation cells (per-point Sankey
        # / time series). nobs is a handful, so these stay in memory in full.
        if self.obs_idx:
            nobs = len(self.obs_idx)
            self.mm_obs = np.zeros((nper_mm, nobs, nidx))
            self.mms_obs = np.zeros((nper_mm, nobs, nsl, nidx_s))
        else:
            self.mm_obs = None
            self.mms_obs = None

        api.initialize()
        try:
            self._bind(api)
            self._check_sizes()
            # --- steady initial SP ---
            # A steady-state period ignores the initial-head file, so the ONLY
            # way to make it land near the dynamic equilibrium is to drive it
            # with the mean actual forcing. When per-cell mean recharge / ETg
            # are supplied (from a prior run) they are used; otherwise fall back
            # to the uniform perc_user, which produces a too-wet near-surface
            # table the transient then has to drain from.
            #
            # WITHOUT ONE (mf6b.steady_first False) the run starts from the
            # heads it was given -- the previous spin-up cycle's last ones --
            # and every period is transient: a periodic spin-up.
            steady = bool(getattr(self.mf6b, 'steady_first', True))
            if not steady:
                print('no steady period: the run starts from the heads it was '
                      'given')
            elif self.steady_perc is not None:
                s_perc = np.asarray(self.steady_perc, dtype=float).ravel()[:self.ncell]
                s_etg = (np.zeros(self.ncell) if self.steady_etg is None else
                         np.asarray(self.steady_etg, dtype=float).ravel()[:self.ncell])
                print('steady state driven by mean recharge %.4g m/d, mean ETg '
                      '%.4g m/d (per cell)' % (s_perc.mean(), s_etg.mean()))
            else:
                s_perc = np.full(self.ncell, float(getattr(cMF, 'perc_user', 0.0)))
                s_etg = np.zeros(self.ncell)
            # end time of every stress period, so the step loop knows when a
            # period is complete even when ATS subdivides it: t_end[n + 1] is
            # the end of MARMITES period n either way
            perlen = list(self.perlen) or [1.0] * nper_mm
            t_end = np.cumsum([1.0 if steady else 0.0]
                              + [float(p) for p in perlen])
            if steady:
                # write the steady fluxes AFTER prepare_time_step (via
                # _advance's callback), or UZF rp reverts SINF to the
                # build-time perc_user; the mean ETg as a flat EVT curve
                s_evt = self._evt_flat(s_etg)

                def _write_steady():
                    self._write_fluxes(s_perc)
                    self._write_evt(s_evt)

                self._advance(api, t_end[0], write_cb=_write_steady)

            # --- transient march, one MM SP per MF6 SP (lagged) ---
            for n in range(nper_mm):
                tstart_MF = n
                heads, exf = self._read_heads_exf()             # end of SP n-1
                rej = self._read_rejinf()                       # returned to surface
                if n == 0 and not steady and self.carry_in:
                    # a periodic cycle's first day: what the previous
                    # cycle's last day left -- MF6 has solved nothing yet
                    exf = np.asarray(self.carry_in['exf'], dtype=float)
                    rej = np.asarray(self.carry_in['rej'], dtype=float)
                    self.etuzf_prev = np.asarray(self.carry_in['etuzf'],
                                                 dtype=float)
                # WP2: Eg/Tg see what remains after the PREVIOUS period's
                # actual UZF ET -- lagged: MM cannot know this period's
                # before MF6 solves (cookbook 2b)
                out = self._mm_step(n, tstart_MF, heads, exf, rej)
                # ALL API inputs (UZF SINF, the EVT curves, SFR INFLOW) must
                # be written after prepare_solve, or MF6 reverts them to the
                # build-time values -- so runoff-to-SFR rides the same
                # callback rather than being written after the advance.
                _o = out
                _ro = (self._runoff_to_deliver(out)
                       if self._delivers_runoff() else None)

                def _write(o=_o, ro=_ro, sp=n):
                    self._write_fluxes(o['perc'], o.get('petuzf'),
                                       etg_cap=o.get('etg_cap'))
                    self._write_evt(o.get('evt'))
                    if ro is not None:
                        self._write_runoff(ro)
                    self._write_openwater_evap(sp)

                self.outer_iters[n] = self._advance(api, t_end[n + 1], write_cb=_write)
                self.heads_hist[n] = heads
                self.exf_hist[n] = exf
                self.perc_hist[n] = out['perc']
                self.rejinf_hist[n] = rej
                mm_cells = np.asarray(out['MM'], dtype=np.float64)
                if self.p_sfr_inflow is not None:
                    self.runoff_hist[n] = mm_cells[:, self._iRo]
                self._crr_record(n, out)
                # WP1d: the open-water evaporation MF6 just simulated goes back
                # into the per-cell vector, so iEow is a measured flux again.
                # It is read AFTER the advance -- it is what the stress period
                # actually removed -- and the runoff injected into SFR above
                # was the gross value, which is correct: MF6 is what
                # evaporates it.
                eow = self._openwater_evap_into(mm_cells)
                if eow is not None:
                    self.evap_hist[n] = eow
                    self.evap_sfr_hist[n], self.evap_lak_hist[n] = \
                        self._eow_split
                    _vol = self._lak_evap_volumes()
                    if _vol is not None:
                        self.lak_evap_hist[n, :_vol.size] = _vol
                # the groundwater ET MF6's EVT took at the head it solved
                # for -- before the PET check reads iETg
                self.etg_hist[n] = self._evt_into(mm_cells)
                # WP2 step 5: what UZF actually took this period -- into the
                # water balance now, and into MM's groundwater ET next period
                et, resid = self._split_etuzf(self._read_etuzf(), n)
                self.etuzf_hist[n] = et
                self.etuzf_prev = et
                self._etuzf_into(mm_cells, et, n, resid)
                if resid.any():
                    self.n_resid += int(np.count_nonzero(resid))
                    self.resid_max = max(self.resid_max, float(resid.max()))
                    self.resid_m3 += float(np.sum(resid / self.conv_fact
                                                  * self.area)
                                           * float(self.perlen[n]
                                                   if n < len(self.perlen)
                                                   else 1.0))
                mms_cells = np.asarray(out['MM_S'], dtype=np.float64)
                _iss = self.ctx.index_S.get('iSsoil')
                if _iss is not None and mms_cells.ndim == 3:
                    self.ssoil_hist[n] = mms_cells[:, :, int(_iss)].sum(axis=1)
                # CATCHMENT means weight each cell by its AREA. A plain mean
                # over cells is the catchment only on a uniform grid: on the
                # Voronoi mesh half the cells cover 5 % of La Mata, refined
                # along the streams, and the plain mean put runoff at 339
                # mm/yr against 61 area-weighted, exfiltration at 310 against
                # 16 (MF6's seepage drain: 15.9) -- 2026-09-23.
                self.wb_ts[n] = np.average(mm_cells, axis=0,
                                           weights=self.area)
                self.wb_map += mm_cells
                self.wb_ts_soil[n] = np.average(mms_cells, axis=0,
                                                weights=self.area)
                self.wb_map_soil += mms_cells
                if self.mm_obs is not None:
                    self.mm_obs[n] = mm_cells[self.obs_idx]
                    self.mms_obs[n] = mms_cells[self.obs_idx]
                self._progress(n, nper_mm)
                if on_sp is not None:
                    on_sp(n, out)
            self._progress(nper_mm - 1, nper_mm, last=True)
            if self.sfr_evap_writes:
                print('stream evaporation: capped at %g x the inflow on %d '
                      'reach-step(s), %.2f %% of those written '
                      '(sfr.evap_inflow_fraction)'
                      % (self.sfr_evap_frac, self.sfr_evap_capped,
                         100.0 * self.sfr_evap_capped / self.sfr_evap_writes))
            # what the unsaturated zone holds at the end: the next spin-up
            # cycle starts from it (clsMF6.uzf_thti_carry) ...
            self.uzf_wc_final = (None if self.p_wcnew is None else
                                 np.array(self.p_wcnew, dtype=float))
            # ... and each pond's stage (a perched pond's storage is not in
            # the heads)
            self.lak_stage_final = (
                None if getattr(self, 'p_lak_stage', None) is None else
                np.array(self.p_lak_stage, dtype=float)[:self.nlakes])
            # ... and what its first day reads as "the previous period"
            _h, _exf = self._read_heads_exf()
            self.carry_out = {'exf': np.array(_exf, dtype=float),
                              'rej': np.array(self._read_rejinf(), dtype=float),
                              'etuzf': np.array(self.etuzf_prev, dtype=float)}
            self.wb_map /= float(nper_mm)
            self.wb_map_soil /= float(nper_mm)
            # Physical-plausibility guard. Exfiltration identically zero over a
            # whole simulation is almost always a configuration fault, not a
            # result: UZF6 computes groundwater discharge ONLY when the UZF
            # package was built with SIMULATE_GWSEEP. This check exists because
            # that option was once omitted and the resulting zero exfiltration
            # looked superficially plausible for a dry period.
            if self.exf_hist.size and not np.any(self.exf_hist > 0.0):
                print('\nWARNING: exfiltration is zero in every cell and every '
                      'stress period.\n         If the water table ever rises '
                      'into the soil zone this is wrong.\n         Check that '
                      + ('the DRN_SEEP package exists and that its drain '
                         'elevations sit at the\n         land surface '
                         '(clsMF6(..., seep="drn")).'
                         if self.p_drnseep is not None else
                         'the UZF package was built with SIMULATE_GWSEEP '
                         '(clsMF6(..., gwseep=True))\n         and that the UZF '
                         'budget contains a groundwater-discharge record.'))
            # Plausibility guard on the water table. A head standing well above
            # the land surface is not a result, it means the seepage boundary
            # cannot discharge fast enough: the first DRN-SEEP run used a
            # conductance of 10 m2/d and the head rose 205 m above ground to
            # force the flux through. Report it rather than let it look like a
            # wet spin-up.
            #
            # Against the LAND surface: the MF6 top is the base of the MMsoil
            # column, and a water table up in the soil column is MARMITES'
            # exfiltration, not a flood. Measured from the top, the LAK run of
            # 2026-10-04 reported "1.4 m ABOVE the land surface" where the
            # highest head stood 0.08 m above it, under a full pond.
            land = getattr(self, 'land_cell', None)
            if land is None:
                land = self.top_cell
            if land is not None and self.heads_hist.size:
                for _line in surface_excess_lines(
                        self.heads_hist, land, self.ncell):
                    print(_line)
                above = self.heads_hist - land[None, :]
                drained = np.asarray(getattr(self, 'drnseep_idx',
                                             np.zeros(0)), int) >= 0
                if drained.size != self.ncell:
                    drained = np.ones(self.ncell, bool)
                if above.max() > 1.0 and self.p_drnseep is not None \
                        and not (above[:, drained] > 1.0).any():
                    # the drain is not the limit where there is none
                    print('         Only cells WITHOUT a seepage drain stand '
                          'that high -- stream and pond cells, where SFR '
                          'and LAK take\n         the seepage through their '
                          'bed: the bed conductance limits them, not the '
                          'drain.')
                elif above.max() > 1.0:
                    if self.p_drnseep is not None:
                        # per CELL: its own exfiltration on its own area.
                        # The maximum rate times the maximum area came from
                        # two different cells -- on the mesh, 0.01..4121 m2.
                        need = float(np.nanmax(
                            (self.exf_hist / self.conv_fact
                             * self.area[None, :])[:, drained]))
                        print('         The seepage drain cannot discharge fast '
                              'enough. Peak seepage is\n         %.0f m3/d per '
                              'cell; a free-draining face needs a conductance of\n'
                              '         about %.0f m2/d to hold the excess head '
                              'near 0.1 m.' % (need, need / 0.1))
                    else:
                        print('         Check the seepage mechanism and the '
                              'initial heads.')

            # water-balance honesty check on the applied percolation
            # rejinf_hist is mm/d PER CELL (_read_rejinf), so it goes back to
            # m3 on each cell's own area. It used to be summed as it stood --
            # mm/d added over cells, reported as m3: 46.1 % on the run of
            # 2026-09-23, where MF6's own UZF budget said 2.6 %.
            applied = float(np.sum(self.perc_hist * self.area[None, :]))
            rejected = float(np.sum(self.rejinf_hist / self.conv_fact
                                    * self.area[None, :]))
            if applied > 0 and rejected > 0:
                print('\nNOTE: UZF rejected %.4g of %.4g m3 of applied percolation '
                      '(%.1f%%).\n      This water re-enters the MARMITES soil '
                      'column from below\n      on the next stress period, and '
                      'whatever the soil cannot hold becomes runoff.'
                      % (rejected, applied, 100.0 * rejected / applied))
            # WP2 2.5b: PET spent once. Unmet demand is normal (a dry soil
            # cannot evaporate what it has not got); ET above the demand is
            # not -- it would be the one-period lag of UZF's ET.
            self._print_pet_balance(nper_mm)
            self._print_crr(nper_mm)
        finally:
            # MF6 can fault inside finalize() when the run is aborted early
            # (the library expects a completed simulation). Never let that
            # mask the original error.
            try:
                api.finalize()
            except Exception as exc:
                print('NOTE: MF6 finalize() failed (%r) -- ignored; this is a '
                      'consequence of stopping early, not the root cause.' % (exc,))
        res = {'heads': self.heads_hist, 'exf': self.exf_hist,
               'perc': self.perc_hist, 'etg': self.etg_hist, 'rejinf': self.rejinf_hist,
               'etuzf': self.etuzf_hist, 'pet_unmet': self.pet_unmet,
               'ssoil': self.ssoil_hist,
               'runoff': self.runoff_hist, 'outer_iters': self.outer_iters,
               # what outer_iters counts: MF6's own step exposes only its
               # linear solves (ITERTOT_TIMESTEP), not the outer iterations
               'iters_kind': ('linear' if self.mf6_step else 'outer'),
               'nonconverged': int(self.n_nonconverged),
               'wb_ts': self.wb_ts, 'wb_map': self.wb_map,
               'wb_ts_soil': self.wb_ts_soil, 'wb_map_soil': self.wb_map_soil}
        # The SFR / LAK split per cell is in wb_ts and wb_map (iEow_sfr,
        # iEow_lak); per cell AND per period it would add ~0.5 GB to a full
        # La Mata run. Each pond's own loss is small and is kept whole.
        if self.lak_evap_hist.size:
            res['lak_evap'] = self.lak_evap_hist      # m3/d per pond
        if self.crr_ts is not None:
            # WP5: the cascade per period, catchment totals in m3/d
            res['crr_ts'] = self.crr_ts
            res['crr_terms'] = np.asarray(self.CRR_TERMS, dtype='S16')
        if self.mm_obs is not None:
            res['mm_obs'] = self.mm_obs
            res['mms_obs'] = self.mms_obs
            res['obs_idx'] = np.asarray(self.obs_idx, int)
            if self.obs_names is not None:
                res['obs_names'] = np.asarray(self.obs_names, dtype='S16')
            # (i, j) of each obs cell, for the per-point aquifer flux lookup
            res['obs_ij'] = np.array([[self.i_arr[p], self.j_arr[p]]
                                      for p in self.obs_idx], int)
        return res

    def _timestep(self, api):
        """Time-step length passed to prepare_time_step().

        NOTE: MF6's BMI implementation *ignores* this argument -- it derives
        the step from TDIS itself -- and get_time_step() legitimately returns
        0.0 before the first step has been prepared. So a non-positive value
        is reported once and passed through, never treated as fatal.
        """
        try:
            dt = float(api.get_time_step())
        except Exception as exc:
            if not self._warned_dt:
                print('NOTE: MF6 get_time_step() unavailable (%r); passing 0.0 '
                      '(MF6 derives the step from TDIS).' % (exc,))
                self._warned_dt = True
            return 0.0
        if not np.isfinite(dt) or dt <= 0.0:
            if not self._warned_dt:
                print('NOTE: MF6 reported dt=%r before the first step; passing it '
                      'through (MF6 derives the step from TDIS).' % dt)
                self._warned_dt = True
            return 0.0 if not np.isfinite(dt) else dt
        return dt

    def _one_step(self, api, write_cb=None):
        """Prepare, solve and finalize a single MF6 time step.

        ``write_cb`` (if given) writes the exchanged fluxes AFTER
        prepare_solve, which is the only point where the write survives. UZF
        re-derives its infiltration (SINF) from the stored period data during
        prepare_time_step AND again during prepare_solve, overwriting anything
        written before them; a value written after prepare_solve is the one the
        solve actually uses. This was pinned down with tests/diag_sinf.py:
        writing after prepare_time_step reverted to perc_user (constant
        recharge, draining water table), writing after prepare_solve sticks.
        """
        if getattr(self, 'mf6_step', False):
            return self._mf6_step(api, write_cb)
        dt = self._timestep(api)
        api.prepare_time_step(dt)
        api.prepare_solve(1)
        if write_cb is not None:
            write_cb()
        kiter = 0
        while kiter < self.max_outer:
            kiter += 1
            if api.solve(1):
                break
        converged = kiter < self.max_outer
        api.finalize_solve(1)
        api.finalize_time_step()
        return kiter, converged

    def _mf6_step(self, api, write_cb=None):
        """One time step by MF6's OWN routine, retried by ATS if it fails.

        The fluxes are written after prepare_time_step -- the period input
        a new PERIOD block would reload is read there -- and before
        do_time_step, whose every try re-applies them (see _bind). A step
        still unconverged when ATS reaches dtmin makes finalize_time_step
        report failure: it is counted, and the run goes on, as before.
        Returns (linear iterations of the last try, converged).
        """
        dt = self._timestep(api)
        api.prepare_time_step(dt)
        if write_cb is not None:
            write_cb()
        api.do_time_step()
        try:
            api.finalize_time_step()
            converged = True
        except Exception:                                 # noqa: BLE001
            converged = False
        kiter = 0
        if self.p_itertot is not None:
            try:
                kiter = int(np.asarray(self.p_itertot).ravel()[0])
            except (TypeError, ValueError, IndexError):   # pragma: no cover
                kiter = 0
        return kiter, converged

    def _advance(self, api, t_end=None, write_cb=None):
        """Advance MF6 to the end of the current stress period.

        With ATS a stress period is NOT one time step: MF6 subdivides any
        period it cannot solve in a single step. Advancing once per period
        would then leave MF6 behind MARMITES by a growing amount -- the two
        models would silently be simulating different days. So the step loop
        runs until MF6's own clock reaches the end of the period.

        ``write_cb`` re-applies the exchanged fluxes after every
        prepare_time_step (the fluxes are stress-period rates, so re-writing
        them on each sub-step is idempotent and keeps them from being reverted
        by rp).
        """
        acc = {}
        t0 = self._now(api)
        kiter, converged = self._one_step(api, write_cb)
        self._keep_sfr_inflow()
        now = self._now(api)
        self._rates_add(acc, None if t0 is None or now is None else now - t0)
        nsub = 1
        if t_end is not None:
            while nsub < self.max_substeps:
                if now is None:                    # pragma: no cover
                    break
                if now >= t_end - 1e-9:
                    break
                k, ok = self._one_step(api, write_cb)
                self._keep_sfr_inflow()
                t1 = self._now(api)
                self._rates_add(acc, None if t1 is None else t1 - now)
                now = t1
                kiter = max(kiter, k)
                converged = converged and ok
                nsub += 1
            else:                                  # pragma: no cover
                raise CouplingError(
                    'MF6 did not reach the end of the stress period after %d '
                    'sub-steps; ATS is shrinking the step without converging.'
                    % self.max_substeps)
        self._rates_close(acc, nsub)
        self.substeps += nsub - 1
        if not converged:
            self.n_nonconverged += 1
        return kiter


if __name__ == '__main__':
    print('Use tests/run_lamata_mf6.py to build and (with mf6 installed) run La Mata.')
