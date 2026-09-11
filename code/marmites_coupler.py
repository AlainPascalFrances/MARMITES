# -*- coding: utf-8 -*-
"""MARMITES <-> MODFLOW 6 API coupler (Phase 3b/3c).

Replaces the removed file-based Picard loop with a single sequential march
over stress periods through the MODFLOW 6 API (libmf6, via modflowapi/xmipy).
Two modes share one exchange layer (MARMITES_code_review.md section 4.3):

  * ``lagged``    -- MMsoil.step() once per SP, using heads and UZF
                     groundwater discharge (exfiltration) from the END of
                     the previous SP; then MF6 advances the SP.
  * ``iterative`` -- MMsoil.step() re-evaluated at every MF6 outer (Newton)
                     iteration from the CURRENT head iterate, with
                     under-relaxation on the exchanged fluxes; removes the
                     one-SP lag entirely (MetaSWAP-style coupling).

Exchange (memory pointers, no files):
  MM -> MF6 : perc (m/d)  -> UZF FINF on the land-surface UZF objects
              ETg  (m/d)  -> WEL Q = -ETg * cell_area (m3/d), AUTO_FLOW_REDUCE
  MF6 -> MM : heads (m)   -> X (per surface cell, dry when h < botm_l0)
              exfiltration-> UZF GWD (m3/d) / cell_area * conv_fact (mm/d, +up)

The MF6 side must be built by ppMF6.marmites_mf6.clsMF6 so that the first
`ncell` UZF objects are the land-surface cells in MARMITES cell order and
the WEL list is in the same order.

The API object is injectable for testing (tests/test_coupler_mock.py runs
both modes against a fake API); a real run needs libmf6:
    from modflowapi import ModflowApi
    MF6Coupler(...).run(ModflowApi('libmf6.so', working_directory=sim_ws))
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import copy
import os
import re

import numpy as np


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
    mode : 'lagged' | 'iterative'
    relax : float
        Under-relaxation factor on the exchanged fluxes in iterative mode
        (new = relax*eval + (1-relax)*previous); 0.5-0.7 recommended.
    max_outer : int
        Safety cap on outer iterations per SP in iterative mode.
    """

    def __init__(self, mm, ctx, state, mf6b, conv_fact=1000.0,
                 mode='lagged', relax=0.6, max_outer=None,
                 obs_idx=None, obs_names=None):
        if mode not in ('lagged', 'iterative'):
            raise CouplingError("mode must be 'lagged' or 'iterative'")
        self.mm, self.ctx, self.state = mm, ctx, state
        # Observation cells at which to keep the FULL per-SP MM flux vectors
        # (for the native per-point Sankey / time series). Tiny next to the
        # whole-grid arrays, so unlike wb_map/wb_ts these are stored in full.
        self.obs_idx = None if obs_idx is None else list(np.asarray(obs_idx, int))
        self.obs_names = list(obs_names) if obs_names is not None else None
        self.mf6b = mf6b
        self.conv_fact = float(conv_fact)
        self.mode = mode
        self.relax = float(relax)
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
        self.p_q = None
        self.p_bound = None
        self.p_gwd = None
        self.p_rejinf = None
        self.p_drnseep = None
        # land-surface elevation per cell, for the water-table plausibility check
        top = getattr(mf6b, 'top', None)
        self.top_cell = (None if top is None else
                         np.asarray(top, dtype=float)[self.i_arr, self.j_arr])
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
        self.p_sfr_simevap = None
        self.p_lak_simevap = None
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
        # The zone rasters are 1-based; a cell outside every mapped zone would
        # index -1 and silently take the LAST zone's Eo.
        self.sfr_evap_zone = np.clip(self.sfr_evap_zone, 0, None)
        self.lak_evap_zone = np.clip(self.lak_evap_zone, 0, None)
        eo = getattr(ctx, 'Eo_zonesSP', None)
        self.Eo_zonesSP = None if eo is None else np.atleast_2d(np.asarray(eo,
                                                                          float))
        self.evap_hist = None

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
        # UZF infiltration. The operative array in MF6's memory manager is
        # SINF (specified infiltration); FINF is only the flopy input keyword.
        # Binding FINF returns a valid-but-non-operative pointer, so the daily
        # percolation written to it is silently ignored and UZF keeps applying
        # the static perc_user from the build -- the recharge coupling looks
        # wired but delivers a constant rate. SINF is therefore tried first.
        self.p_finf, self.addr_finf = self._bind_first(
            api, [('SINF', f'{name}/UZF'), ('SINF', f'{name}/UZF-1'),
                  ('FINF', f'{name}/UZF'), ('FINF', f'{name}/UZF-1')],
            'UZF infiltration (SINF/FINF)', min_size=self.ncell)
        print('coupler: UZF infiltration bound to %s' % self.addr_finf)
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
        # WEL rates: write Q ONLY.
        #
        # MF6 6.7 exposes both WEL/Q and WEL/BOUND. Q is the rate array the
        # package uses. BOUND is *not* a plain writable array here: bisection
        # on La Mata (tests/diagnose_coupling.py) showed that writing FINF+Q
        # runs cleanly, while additionally writing BOUND corrupts MF6's heap
        # and kills the process inside the next prepare_time_step() with no
        # error message. BOUND is therefore only used when Q is unavailable.
        self.p_q, self.addr_q = self._bind_first(
            api, [('Q', f'{name}/WEL'), ('Q', f'{name}/WEL-1')],
            'WEL rates (Q)', required=False)
        self.p_bound = None
        if self.p_q is None:
            self.p_bound, self.addr_bound = self._bind_first(
                api, [('BOUND', f'{name}/WEL'), ('BOUND', f'{name}/WEL-1')],
                'WEL boundary array (BOUND)', required=True)
            print('NOTE: WEL/Q not exposed; falling back to BOUND.')
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
        if self.p_finf.size < self.ncell:
            raise CouplingError(
                'UZF FINF has %d entries but MARMITES has %d surface cells: the '
                'bound address is wrong or the UZF package was built differently.'
                % (self.p_finf.size, self.ncell))
        if self.p_q is not None and self.p_q.size < self.ncell:
            raise CouplingError('WEL Q holds %d entries but MARMITES has %d wells.'
                                % (self.p_q.size, self.ncell))
        if self.p_bound is not None:
            nb = self.p_bound.shape[0] if self.p_bound.ndim == 1 else \
                max(self.p_bound.shape)
            if nb < self.ncell:
                raise CouplingError(
                    'WEL BOUND holds %d entries but MARMITES has %d wells.'
                    % (nb, self.ncell))
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
                rep['nonconverged_lst'] = txt.count('did not converge')
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
        if not rep['ok']:
            msg = 'MF6 solution is not valid:\n  - ' + '\n  - '.join(rep['messages'])
            if raise_on_fail:
                raise CouplingError(msg)
            print('\nWARNING: ' + msg)
        elif d is not None:
            print('solution check: converged, cumulative discrepancy %.4f%%, '
                  '%d extra ATS sub-step(s)' % (d, rep['substeps']))
        return rep

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

    def _read_heads_exf(self):
        """Heads [m] and groundwater exfiltration [mm/d, + into the soil].

        The seepage comes from whichever mechanism the model was built with:
        UZF SIMULATE_GWSEEP (GWD, positive out of the aquifer) or the
        land-surface DRN_SEEP package (SIMVALS, negative out of the aquifer).
        """
        heads = np.asarray(self.p_x, dtype=float)[self.x_index]
        if self.p_drnseep is not None:
            sim = np.asarray(self.p_drnseep, dtype=float).ravel()
            q = np.zeros(self.ncell)
            has = self.drnseep_idx >= 0
            idx = self.drnseep_idx[has]
            if idx.size and int(idx.max()) < sim.size:
                # drains remove water -> SIMVALS <= 0; exfiltration is positive
                q[has] = np.maximum(-sim[idx], 0.0)
            return heads, q / self.area * self.conv_fact
        if self.p_gwd is None:
            return heads, np.zeros(self.ncell)
        gwd = np.asarray(self.p_gwd, dtype=float).ravel()[:self.ncell]  # m3/d, + to surface
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
        rej = np.abs(np.asarray(self.p_rejinf, dtype=float).ravel()[:self.ncell])
        return rej / self.area * self.conv_fact

    def _write_runoff(self, ro):
        """Deliver the MARMITES runoff of this SP to the stream reaches.

        ``ro`` is per cell in mm/d. Only cells that host a reach contribute;
        runoff generated off-channel is handled by the CRR cascade, not here.
        SFR INFLOW is a volumetric rate (m3/d).
        """
        if self.p_sfr_inflow is None or self.sfr_reach_idx is None:
            return
        q = np.asarray(ro, dtype=float) / self.conv_fact * self.area   # m3/d
        inflow = np.zeros(self.p_sfr_inflow.size)
        has = self.sfr_reach_idx >= 0
        np.add.at(inflow, self.sfr_reach_idx[has], np.maximum(q[has], 0.0))
        self.p_sfr_inflow[:] = inflow

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
            v = eo[np.minimum(self.sfr_evap_zone, len(eo) - 1)]
            self.p_sfr_evap[:min(self.p_sfr_evap.size, self.nreaches)] = \
                v[:min(self.p_sfr_evap.size, self.nreaches)]
            wrote = float(v.mean())
        if self.p_lak_evap is not None and self.nlakes:
            v = eo[np.minimum(self.lak_evap_zone, len(eo) - 1)]
            self.p_lak_evap[:min(self.p_lak_evap.size, self.nlakes)] = \
                v[:min(self.p_lak_evap.size, self.nlakes)]
            wrote = float(v.mean()) if wrote is None else wrote
        return wrote

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
        out = np.zeros(self.ncell, dtype=float)
        if self.p_sfr_simevap is not None and self.nreaches:
            se = np.asarray(self.p_sfr_simevap, dtype=float)
            has = self.sfr_reach_idx >= 0
            idx = self.sfr_reach_idx[has]
            ok = idx < se.size
            np.add.at(out, np.nonzero(has)[0][ok], se[idx[ok]])
        if self.p_lak_simevap is not None and self.nlakes:
            le = np.asarray(self.p_lak_simevap, dtype=float)
            for L in range(min(self.nlakes, le.size)):
                n = int(self.lak_cell_idx[L])
                if n >= 0:
                    out[n] += float(le[L])
        # MF6 reports a removal as a positive rate here; carry it as a
        # positive loss, like every other MARMITES flux.
        return np.abs(out) / self.area * self.conv_fact          # mm/d

    def _openwater_evap_into(self, mm_cells):
        """Put this SP's open-water evaporation into ``iEow`` of the vector.

        Deliberately does NOT reduce ``iRo``. Per cell it could not: a reach
        carries water from upstream, so its evaporation can exceed the runoff
        that cell generated, and subtracting would drive Ro negative. Eow is a
        loss from the CHANNEL, downstream of the runoff -- which is what the
        water-budget figures now say.
        """
        if self._iEow is None:
            return None
        eow = self._read_openwater_evap()
        mm_cells[:, self._iEow] = eow
        return eow

    def _write_fluxes(self, perc, etg):
        self.p_finf[:self.ncell] = perc                               # m/d
        if self.p_finf.shape[0] > self.ncell:
            self.p_finf[self.ncell:] = 0.0
        q = -np.asarray(etg, dtype=float) * self.area                 # m3/d, sink
        if self.p_q is not None:
            self.p_q[:self.ncell] = q
            return                       # never also write BOUND (see _bind)
        b = self.p_bound
        if b is None:
            return
        if b.ndim == 1:
            b[:self.ncell] = q
        elif b.shape[0] >= self.ncell:                                # (maxbound, naux+1)
            b[:self.ncell, 0] = q
        else:                                                         # (naux+1, maxbound)
            b[0, :self.ncell] = q

    @staticmethod
    def _clone_state(state):
        return copy.deepcopy(state)

    @staticmethod
    def _restore_state(state, bak):
        state.Ssoil_ini[:] = bak.Ssoil_ini

    # ---------------- main drive --------------------------------------- #

    def run(self, api, on_sp=None):
        """Drive the coupled model with an initialized-able MF6 API object.

        api : modflowapi.ModflowApi or compatible (tests inject a fake).
        on_sp : optional callback(n, out_dict) after each stress period.
        """
        cMF = self.mf6b.cMF
        nper_mm = int(cMF.nper)
        self.heads_hist = np.zeros((nper_mm, self.ncell))
        self.exf_hist = np.zeros((nper_mm, self.ncell))
        self.perc_hist = np.zeros((nper_mm, self.ncell))
        self.etg_hist = np.zeros((nper_mm, self.ncell))
        self.rejinf_hist = np.zeros((nper_mm, self.ncell))
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
        self.evap_hist = np.zeros((nper_mm, self.ncell))
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
            if self.steady_perc is not None:
                s_perc = np.asarray(self.steady_perc, dtype=float).ravel()[:self.ncell]
                s_etg = (np.zeros(self.ncell) if self.steady_etg is None else
                         np.asarray(self.steady_etg, dtype=float).ravel()[:self.ncell])
                print('steady state driven by mean recharge %.4g m/d, mean ETg '
                      '%.4g m/d (per cell)' % (s_perc.mean(), s_etg.mean()))
            else:
                s_perc = np.full(self.ncell, float(getattr(cMF, 'perc_user', 0.0)))
                s_etg = np.zeros(self.ncell)
            # end time of every stress period, so the step loop knows when a
            # period is complete even when ATS subdivides it
            perlen = list(self.perlen) or [1.0] * nper_mm
            t_end = np.cumsum([1.0] + [float(p) for p in perlen])
            # write the steady fluxes AFTER prepare_time_step (via _advance's
            # callback), or UZF rp reverts SINF to the build-time perc_user
            self._advance(api, t_end[0],
                          write_cb=lambda: self._write_fluxes(s_perc, s_etg))

            # --- transient march, one MM SP per MF6 SP ---
            for n in range(nper_mm):
                tstart_MF = n
                if self.mode == 'lagged':
                    heads, exf = self._read_heads_exf()         # end of SP n-1
                    rej = self._read_rejinf()                   # returned to surface
                    out = self.mm.step(self.ctx, n, tstart_MF, heads, exf, self.state,
                                       rejinf_cell=rej)
                    # ALL API inputs (UZF SINF, WEL Q, SFR INFLOW) must be
                    # written after prepare_solve, or MF6 reverts them to the
                    # build-time values -- so runoff-to-SFR rides the same
                    # callback rather than being written after the advance.
                    _o = out
                    _ro = (np.asarray(out['MM'], dtype=np.float64)[:, self._iRo]
                           if self.p_sfr_inflow is not None else None)

                    def _write(o=_o, ro=_ro, sp=n):
                        self._write_fluxes(o['perc'], o['etg'])
                        if ro is not None:
                            self._write_runoff(ro)
                        self._write_openwater_evap(sp)

                    self.outer_iters[n] = self._advance(api, t_end[n + 1], write_cb=_write)
                else:
                    out, heads, exf, rej = self._iterative_sp(api, n, tstart_MF,
                                                             t_end[n + 1])
                self.heads_hist[n] = heads
                self.exf_hist[n] = exf
                self.perc_hist[n] = out['perc']
                self.etg_hist[n] = out['etg']
                self.rejinf_hist[n] = rej
                mm_cells = np.asarray(out['MM'], dtype=np.float64)
                if self.p_sfr_inflow is not None:
                    self.runoff_hist[n] = mm_cells[:, self._iRo]
                # WP1d: the open-water evaporation MF6 just simulated goes back
                # into the per-cell vector, so iEow is a measured flux again.
                # It is read AFTER the advance -- it is what the stress period
                # actually removed -- and the runoff injected into SFR above
                # was the gross value, which is correct: MF6 is what
                # evaporates it.
                eow = self._openwater_evap_into(mm_cells)
                if eow is not None:
                    self.evap_hist[n] = eow
                mms_cells = np.asarray(out['MM_S'], dtype=np.float64)
                self.wb_ts[n] = mm_cells.mean(axis=0)
                self.wb_map += mm_cells
                self.wb_ts_soil[n] = mms_cells.mean(axis=0)
                self.wb_map_soil += mms_cells
                if self.mm_obs is not None:
                    self.mm_obs[n] = mm_cells[self.obs_idx]
                    self.mms_obs[n] = mms_cells[self.obs_idx]
                if on_sp is not None:
                    on_sp(n, out)
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
            if self.top_cell is not None and self.heads_hist.size:
                above = self.heads_hist - self.top_cell[None, :]
                if above.max() > 1.0:
                    ncell_over = int((above > 1.0).any(axis=0).sum())
                    nsp_over = int((above > 1.0).any(axis=1).sum())
                    print('\nWARNING: the water table rises up to %.1f m ABOVE the '
                          'land surface\n         (%d of %d cells, %d of %d stress '
                          'periods).' % (above.max(), ncell_over, self.ncell,
                                         nsp_over, self.heads_hist.shape[0]))
                    if self.p_drnseep is not None:
                        need = float(np.nanmax(self.exf_hist) / self.conv_fact
                                     * np.max(self.area))
                        print('         The seepage drain cannot discharge fast '
                              'enough. Peak seepage is\n         %.0f m3/d per '
                              'cell; a free-draining face needs a conductance of\n'
                              '         about %.0f m2/d to hold the excess head '
                              'near 0.1 m.' % (need, need / 0.1))
                    else:
                        print('         Check the seepage mechanism and the '
                              'initial heads.')

            # water-balance honesty check on the applied percolation
            applied = float(np.sum(self.perc_hist * self.area[None, :]))
            rejected = float(np.sum(self.rejinf_hist))
            if applied > 0 and rejected > 0:
                print('\nNOTE: UZF rejected %.4g of %.4g m3 of applied percolation '
                      '(%.1f%%).\n      This water re-enters the MARMITES soil '
                      'column from below\n      on the next stress period, and '
                      'whatever the soil cannot hold becomes runoff.'
                      % (rejected, applied, 100.0 * rejected / applied))
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
               'runoff': self.runoff_hist, 'outer_iters': self.outer_iters,
               'wb_ts': self.wb_ts, 'wb_map': self.wb_map,
               'wb_ts_soil': self.wb_ts_soil, 'wb_map_soil': self.wb_map_soil}
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
        kiter, converged = self._one_step(api, write_cb)
        nsub = 1
        if t_end is not None:
            while nsub < self.max_substeps:
                try:
                    now = float(api.get_current_time())
                except Exception:                  # pragma: no cover
                    break
                if now >= t_end - 1e-9:
                    break
                k, ok = self._one_step(api, write_cb)
                kiter = max(kiter, k)
                converged = converged and ok
                nsub += 1
            else:                                  # pragma: no cover
                raise CouplingError(
                    'MF6 did not reach the end of the stress period after %d '
                    'sub-steps; ATS is shrinking the step without converging.'
                    % self.max_substeps)
        self.substeps += nsub - 1
        if not converged:
            self.n_nonconverged += 1
        return kiter

    def _iterative_sp(self, api, n, tstart_MF, t_end=None):
        """One SP with MM embedded in the MF6 outer-iteration loop.

        Each outer iteration: restore the soil state, re-evaluate step()
        from the current head iterate, under-relax the exchanged fluxes.
        The state advance of the LAST evaluation is kept.

        The MM coupling is done on the FIRST time step of the period. If ATS
        subdivides the period, the remaining sub-steps are advanced with the
        fluxes already converged on -- re-running the soil model per sub-step
        would advance the soil state several times within one MARMITES day.
        """
        bak = self._clone_state(self.state)
        dt = self._timestep(api)
        api.prepare_time_step(dt)
        api.prepare_solve(1)
        perc_prev = None
        etg_prev = None
        out = None
        heads = exf = None
        kiter = 0
        while kiter < self.max_outer:
            kiter += 1
            heads, exf = self._read_heads_exf()                # current iterate
            rej = self._read_rejinf()
            self._restore_state(self.state, bak)
            out = self.mm.step(self.ctx, n, tstart_MF, heads, exf, self.state,
                               rejinf_cell=rej)
            perc = np.asarray(out['perc'], dtype=float)
            etg = np.asarray(out['etg'], dtype=float)
            if perc_prev is not None:                          # under-relaxation
                perc = self.relax * perc + (1.0 - self.relax) * perc_prev
                etg = self.relax * etg + (1.0 - self.relax) * etg_prev
            perc_prev, etg_prev = perc, etg
            self._write_fluxes(perc, etg)
            self._write_openwater_evap(n)
            if api.solve(1):
                break
        converged = kiter < self.max_outer
        api.finalize_solve(1)
        api.finalize_time_step()
        if not converged:
            self.n_nonconverged += 1
        # finish the period if ATS split it
        nsub = 1
        if t_end is not None:
            while nsub < self.max_substeps:
                try:
                    if float(api.get_current_time()) >= t_end - 1e-9:
                        break
                except Exception:                  # pragma: no cover
                    break
                # re-apply the converged fluxes each sub-step, or rp reverts
                # SINF to the build-time value on the next prepare_time_step
                def _re_apply(p=perc_prev, e=etg_prev, sp=n):
                    self._write_fluxes(p, e)
                    self._write_openwater_evap(sp)

                k, ok = self._one_step(api, write_cb=_re_apply)
                kiter = max(kiter, k)
                if not ok:
                    self.n_nonconverged += 1
                nsub += 1
        self.substeps += nsub - 1
        self.outer_iters[n] = kiter
        return out, heads, exf, rej


if __name__ == '__main__':
    print('Use tests/run_lamata_mf6.py to build and (with mf6 installed) run La Mata.')
