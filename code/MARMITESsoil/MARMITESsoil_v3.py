# -*- coding: utf-8 -*-
"""
MARMITES is a distributed depth-wise lumped-parameter model for
spatio-temporal assessment of water fluxes in the unsaturated zone.
MARMITES is a French word to design a big cooking pot used by sorcerers
for all kinds of experiments!
The main objective of MARMITES development is to partition rainfall
into the several fluxes of the unsaturated and saturated zones.
It applies the concepts enunciated by Lubczynski (2009):
ET = Es + Ei + ETu + ETg;  ETu = Eu + Tu;  ETg = Eg + Tg
MARMITES computes on a daily basis: interception, surface storage and
runoff, evaporation and transpiration from soil and groundwater
(coupled with MODFLOW 6), soil moisture storage and gross recharge.
Input driving forces are rainfall and potential evapo(transpi)ration
daily time series.

References:
Lubczynski (2009), Francés & Lubczynski (2023, Front. Water 5:1055934)

Phase-1 port (2026): Python 3.12, float64 arithmetic (Decimal removed),
bug fixes; structure kept close to v0.3 for regression comparison.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"
__date__ = "2026"

import os
import sys
from types import SimpleNamespace


# ----------------------------------------------------------------- tallies
# A CONDITION THAT IS NORMAL FOR MONTHS IS NOT NEWS EVERY TIME IT HOLDS.
# The soil sitting at wilting point is expected in a semi-arid summer; it
# printed two lines per cell per stress period, 8205 of them in the first
# ~240 periods of La Mata, which is a third of the whole log. Counted here
# and reported once, with how far out of range it went -- the same treatment
# the mesh and drain warnings already get.
_TALLIES = {}


def tally(what, amount=None):
    """Record one occurrence of a recurring condition."""
    rec = _TALLIES.setdefault(what, {'n': 0, 'lo': None, 'hi': None})
    rec['n'] += 1
    if amount is None:
        return
    a = float(amount)
    rec['lo'] = a if rec['lo'] is None else min(rec['lo'], a)
    rec['hi'] = a if rec['hi'] is None else max(rec['hi'], a)


def report_tallies(out=None, clear=True):
    """Print what was tallied, once. Returns the lines, for a test."""
    lines = []
    for what in sorted(_TALLIES):
        rec = _TALLIES[what]
        if rec['lo'] is None:
            lines.append('%s: %d time(s)' % (what, rec['n']))
        else:
            lines.append('%s: %d time(s), out of range by %.4f..%.4f'
                         % (what, rec['n'], rec['lo'], rec['hi']))
    if lines and out is not False:
        print('\nconditions met during the run (counted, not repeated):')
        for ln in lines:
            print('   %s' % ln)
    if clear:
        _TALLIES.clear()
    return lines

import numpy as np


def _grid_module():
    """Import marmites_grid robustly (code/ is not always on sys.path:
    the test-suite loads modules by file path)."""
    if 'marmites_grid' in sys.modules:
        return sys.modules['marmites_grid']
    try:
        import marmites_grid
        return marmites_grid
    except ModuleNotFoundError:
        import importlib.util
        trunk = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        spec = importlib.util.spec_from_file_location(
            'marmites_grid', os.path.join(trunk, 'marmites_grid.py'))
        mod = importlib.util.module_from_spec(spec)
        sys.modules['marmites_grid'] = mod
        spec.loader.exec_module(mod)
        return mod


class MarmitesSoilError(Exception):
    """Raised on invalid states in the MARMITES soil-zone computation."""


class clsMMsoil:
    """
    SOIL class: compute the water balance in the soil (depth-wise, each
    reservoir corresponds to a soil horizon).

    INPUTS
        PARAMETERS
            Sm      max. soil moisture storage (porosity) [-]
            Sfc     soil moisture storage at field capacity [-]
            Sr      residual soil moisture storage [-]
            Ssoil_ini  initial soil moisture [-]
            Ks      saturated hydraulic conductivity [mm/d]
        STATE VARIABLES
            P (rainfall), PT (pot. transpiration), PE (pot. evaporation) [mm/d]
    OUTPUTS
        Pe (effective rainfall), ETsoil, Ssoil, Rp (percolation),
        Ssurf (ponding), Ro (runoff), Eg, Tg [mm/d]
    """

    def __init__(self, hnoflo):
        self.hnoflo = hnoflo

        # Shah, Nachabe & Ross (2007), Ground Water 45(3):329-338, table 3
        # dll [cm], y0 [-], b [cm^-1], ext_d [cm]
        self.paramEg = {
            'sand':             {'dll': 16.0,  'y0': 0.000, 'b': 0.171, 'ext_d': 50.0},
            'loamy sand':       {'dll': 21.0,  'y0': 0.002, 'b': 0.130, 'ext_d': 70.0},
            'sandy loam':       {'dll': 30.0,  'y0': 0.004, 'b': 0.065, 'ext_d': 130.0},
            'sandy clay loam':  {'dll': 30.0,  'y0': 0.006, 'b': 0.046, 'ext_d': 200.0},
            'sandy clay':       {'dll': 20.0,  'y0': 0.005, 'b': 0.042, 'ext_d': 210.0},
            'loam':             {'dll': 33.0,  'y0': 0.004, 'b': 0.028, 'ext_d': 265.0},
            'silty clay':       {'dll': 37.0,  'y0': 0.007, 'b': 0.046, 'ext_d': 335.0},
            'clay loam':        {'dll': 33.0,  'y0': 0.008, 'b': 0.027, 'ext_d': 405.0},
            'silt loam':        {'dll': 38.0,  'y0': 0.006, 'b': 0.019, 'ext_d': 420.0},
            'silt':             {'dll': 31.0,  'y0': 0.007, 'b': 0.021, 'ext_d': 430.0},
            'silty clay loam':  {'dll': 40.0,  'y0': 0.007, 'b': 0.021, 'ext_d': 450.0},
            'clay':             {'dll': 45.0,  'y0': 0.006, 'b': 0.019, 'ext_d': 620.0},
            'sandy loam field_enrico': {'dll': 100.0, 'y0': 0.000, 'b': 0.013, 'ext_d': 475.0},
            'sandy loam field': {'dll': 115.3, 'y0': 0.023, 'b': 0.013, 'ext_d': 1000.0},
        }

    # ------------------------------------------------------------------ #

    @staticmethod
    def _perc(s_tmp, Sm, Sfc, Ks, perlen):
        """Percolation of gravitational water [mm/d]. All args float."""
        if not np.isfinite(s_tmp):
            raise MarmitesSoilError(f'Non-finite soil storage in percolation: {s_tmp!r}')
        if s_tmp <= Sfc:
            return 0.0
        Sg = (s_tmp - Sfc) / (Sm - Sfc)  # gravitational water fraction
        if Ks * Sg * perlen > (s_tmp - Sfc):
            return (s_tmp - Sfc) / perlen
        return Ks * Sg

    @staticmethod
    def _evp(s_tmp, Sm, Sr, pet, perlen):
        """Actual evaporation/transpiration from a soil reservoir [mm/d]."""
        if not np.isfinite(s_tmp):
            raise MarmitesSoilError(f'Non-finite soil storage in evapotranspiration: {s_tmp!r}')
        if s_tmp <= Sr:
            return 0.0
        # effective saturation, at most 1: a store above porosity (rounding
        # in the upward cascade) must not evaporate MORE than its demand --
        # total ET can never be above PET
        Se = min((s_tmp - Sr) / (Sm - Sr), 1.0)
        if pet * Se * perlen > (s_tmp - Sr):
            return (s_tmp - Sr) / perlen
        return pet * Se

    @staticmethod
    def _ktg(theta, Sr, Sm, kTg_min, kTg_max, kT_f, kT_s):
        """The kTg factor of groundwater transpiration at the soil moisture
        ``theta`` of one layer: kTg_min at wilting point, kTg_max at
        saturation, a logistic curve between."""
        if theta > Sr:
            if theta < Sm:
                Ssoil_norm = (theta - Sr) / (Sm - Sr)
                return kTg_max - (kTg_max - kTg_min) / (
                    1.0 + np.exp((Ssoil_norm - kT_f) / kT_s))
            if theta > Sm:
                tally('Tg: soil moisture above porosity', theta - Sm)
            return kTg_max
        if theta < Sr:
            # COUNTED, NOT PRINTED. A semi-arid catchment sits at wilting
            # point for months on end, so this fired 8205 times in the first
            # ~240 stress periods of La Mata -- a third of every line in the
            # log, burying the run's own output. Reported once, with the
            # count and the range, by report_tallies().
            tally('Tg: soil moisture below wilting point', Sr - theta)
        return kTg_min

    @staticmethod
    def _eg(PE, dgwt, head, sy, p):
        """Groundwater evaporation, eq. 17 of Shah et al. (2007), and the
        water table it leaves: ``(Eg [mm], dgwt [cm], head [cm])``.

        ``dgwt`` and ``head`` are in cm, as Shah's parameters ``p`` are
        (``dll``, ``ext_d``); Eg over ``sy`` is the drawdown it causes.

        THE WATER TABLE DROPS BY WHAT EG TAKES, in both branches. When that
        drawdown would pass the extinction depth, Eg is cut to what brings
        the table exactly there -- and the head drops by that same distance,
        ``ext_d - dgwt``. It used to drop by ``ext_d`` itself, the full
        extinction depth (0.5 to 10 m below the table in Shah's table), so
        the groundwater transpiration after it saw a water table metres too
        deep and roots cut off from water they were standing in.
        """
        y0, b, dll, ext_d = p['y0'], p['b'], p['dll'], p['ext_d']
        if dgwt <= dll:
            Eg = PE
        elif dgwt < ext_d:
            # at most PE: every tabulated y0 is > 0, so the fitted curve
            # starts ABOVE 1 just past dll (up to 1.023) -- and total ET can
            # never be above PET
            Eg = PE * min(y0 + np.exp(-b * (dgwt - dll)), 1.0)
        else:
            Eg = 0.0
        if Eg > 0.0:
            # water-table drawdown feedback (0.1: mm -> cm)
            if (dgwt + 0.1 * Eg / sy) > ext_d:
                drop = ext_d - dgwt
                Eg = 10.0 * drop * sy
                dgwt = ext_d
                head -= drop
            else:
                dgwt += 0.1 * Eg / sy
                head -= 0.1 * Eg / sy
        return Eg, dgwt, head

    # ------------------------------------------------------------------ #

    def flux(self, cMF, perleni, Pe, PT, PE, Zr_elev, VEGarea,
             HEADSini, TopSoilLay, BotSoilLay, Tl, nsl, Sm, Sfc, Sr, Ks,
             Ssoil_ini, EXF_ini, dgwt, st, i, j, n,
             kTg_min, kTg_max, kT_f, kT_s, NVEG, LAIveg, REJINF_ini=0.0,
             ETUZF_prev=0.0, RUNON=0.0, IRR=0.0, GW_EVT=False):
        """Soil water balance of one cell for one stress period.

        All storages in mm, all fluxes in mm/d. PT and LAIveg are 1-D
        arrays of length NVEG (scalars per vegetation type); PT is
        consumed (mutated) by soil and groundwater transpiration.

        REJINF_ini : rejected infiltration returned by the groundwater model
            [mm/d] -- percolation MARMITES delivered that the unsaturated zone
            could not accept (cell saturated, or rate above VKS).

            It enters the BOTTOM soil layer, exactly like groundwater
            exfiltration, and is carried upward by the saturation-excess
            cascade. This follows Frances & Lubczynski (2023) Eq. 1, whose
            only subsurface->surface input is Exf_g1 -- exfiltration arriving
            *from the topmost soil layer* -- and section 2.3, where water the
            subsurface cannot accept "eventually creat[es] saturation-excess
            overland flow (Dunnian flow) if all the soil layers turn
            saturated". So it fills the soil from below, then the surface
            store, and the excess becomes runoff.

        ETUZF_prev : the deep unsaturated zone's ACTUAL evapotranspiration
            [mm/d] -- UZF's, read from the MF6 budget -- of the PREVIOUS
            stress period (WP2, step 5 of the demand chain).

        WP2, THE DEMAND CHAIN (cookbook 2b): PET is spent once, in order --
        interception, open water, the soil (Esoil, Tsoil), the deep
        unsaturated zone (ETuzf, in UZF), groundwater (Eg, Tg). What the soil
        leaves is returned as PETuzf, UZF's demand for this stress period;
        Eg and Tg then see what remains after UZF's ACTUAL uptake -- never
        the demand written to UZF, which a dry deep zone cannot meet, or ETg
        would be starved exactly when the deep-rooted trees draw on the water
        table. LAGGED: MM cannot know this period's UZF uptake before MF6
        solves, so the previous period's actual is used; at daily stress
        periods the lag is one day.

        WP1d: MMsoil no longer has a surface RESERVOIR. Ponding capacity and
        open-water evaporation moved to MODFLOW with the water -- SFR for the
        channels, LAK for the charcas, both evaporating from the Eo forcing --
        so nothing is carried between stress periods and everything the soil
        cannot take becomes runoff in the same step.

            Zero for the uncoupled/file-based path, so legacy behaviour is
            unchanged.

        RUNON : runoff the CRR cascade brings from upslope cells [mm/d, per
            unit of soil area] (WP5). It joins the surface water AFTER the
            cell's own rain and exfiltration and infiltrates TOP-DOWN by the
            same law, Eq. 1b, I = min[Ssurf, D1 (phi1 - theta1)]; what the
            soil cannot take leaves as runoff and runs on downslope. The
            bottom-up return of REJINF_ini above is a different flux and is
            untouched (cookbook D4). The REINFILTRATION -- returned last --
            is the part of the run-on the soil took: the capacity the cell's
            own water left, at most the run-on.

        IRR : the SPRINKLER IRRIGATION in Pe [mm over the step, per unit of
            soil area], when soil.irr_infiltration = 'column'; 0 otherwise.
            Sprinklers run for hours at a rate the soil absorbs and the water
            soaks down the column within the day, but Eq. 1b, applied once a
            day, lets it into the TOP horizon only: 2026-10-04, 25 mm on a
            La Mata field whose top horizon had ~19 mm of room ran off as a
            stream pulse. So the irrigation the top horizon cannot take fills
            the horizons below, top-down, each up to saturation -- the mirror
            of the upward saturation-excess cascade -- before any of it runs
            off. Rain, exfiltration and run-on keep Eq. 1b. Booked as
            infiltration (I) and as percolation between the horizons it
            crosses (Rp), so every layer's balance still closes.

        GW_EVT : the coupled run (ctx.gw_evt, set by MF6Coupler). Groundwater
            ET is then MF6's, taken at the head it solves for by EVT
            (marmites_evt): Eg and Tg are returned as 0 and what each source
            COULD take today -- the same demand, soil moisture, root tips and
            Shah curve, at the start-of-day head, with no drawdown of
            MMsoil's own -- comes back last, as a dict. Off (the uncoupled
            path), MMsoil takes Eg and Tg itself, with its drawdown, and the
            dict is None.
        """

        if EXF_ini < 0.0:
            print('WARNING!\nEXFg < 0.0, value %.6f corrected to 0.0.' % EXF_ini)
            EXF_ini = 0.0
        # water returned by the subsurface enters the soil column at its base,
        # whether it is groundwater exfiltration or refused infiltration
        EXF_ini = float(EXF_ini) + abs(float(REJINF_ini))

        # INITIALIZATION
        SAT = np.zeros(nsl, dtype=bool)
        perlen = float(cMF.perlen[n])
        Pe = float(Pe)
        PE = float(PE)
        EXF_ini = float(EXF_ini)
        Tl = np.asarray(Tl, dtype=np.float64)
        Sm = np.asarray(Sm, dtype=np.float64)
        Sfc = np.asarray(Sfc, dtype=np.float64)
        Sr = np.asarray(Sr, dtype=np.float64)
        Ks = np.asarray(Ks, dtype=np.float64)
        PT = np.asarray(PT, dtype=np.float64).reshape(NVEG).copy()
        LAIveg = np.asarray(LAIveg, dtype=np.float64).reshape(NVEG)
        Zr_elev = np.asarray(Zr_elev, dtype=np.float64).reshape(NVEG)
        if VEGarea is None:
            VEGarea = np.zeros(NVEG, dtype=np.float64)
        else:
            VEGarea = np.asarray(VEGarea, dtype=np.float64).reshape(NVEG)

        # Surface water WITHIN this step. Not a storage: nothing is
        # carried over, it is only how Pe and exfiltration reach the soil.
        Ssurf_tmp = Pe

        # Soil storages [mm]
        Ssoil_tmp = np.asarray(Ssoil_ini, dtype=np.float64)[:nsl].copy()
        Rexf_tmp = np.zeros(nsl, dtype=np.float64)

        # EXFILTRATION from MF (groundwater seepage into the soil column)
        Ssoil_tmp[nsl - 1] += EXF_ini * perleni

        # SOIL EXF: saturation excess cascades upward
        # (bug fix vs v0.3: Rexf_tmp was a Python list, `Rexf_tmp /= perlen`
        # raised TypeError; now a numpy array)
        if EXF_ini > 0.0:
            for l in range(nsl - 1, -1, -1):
                if Ssoil_tmp[l] >= Sm[l] * Tl[l]:
                    Rexf_tmp[l] = Ssoil_tmp[l] - Sm[l] * Tl[l]
                    Ssoil_tmp[l] = Sm[l] * Tl[l]
                if l != 0:
                    Ssoil_tmp[l - 1] += Rexf_tmp[l]
            Rexf_tmp /= perlen
            # excess reaching the surface
            Ssurf_tmp += Rexf_tmp[0] * perlen

        IRR = float(IRR)

        def _taken(surf):
            """What the soil takes of ``surf`` [mm]: Eq. 1b into the top
            horizon, then the irrigation in the excess into those below."""
            room0 = Sm[0] * Tl[0] - Ssoil_tmp[0]
            took = room0 if surf > room0 else surf
            down = min(surf - took, IRR) if IRR > 0.0 else 0.0
            for l in range(1, nsl):
                if down <= 0.0:
                    break
                t = min(down, max(Sm[l] * Tl[l] - Ssoil_tmp[l], 0.0))
                took += t
                down -= t
            return took

        # RUN-ON from upslope (WP5 CRR): on top of the cell's own water, so
        # the reinfiltration is what the soil takes beyond that water
        REinf = 0.0
        if RUNON > 0.0:
            own = Ssurf_tmp
            Ssurf_tmp += RUNON * perlen
            REinf = (_taken(Ssurf_tmp) - _taken(own)) / perlen

        # INFILTRATION I into the first soil layer
        if Ssurf_tmp > (Sm[0] * Tl[0] - Ssoil_tmp[0]):
            I = Sm[0] * Tl[0] - Ssoil_tmp[0]
            Ssurf_tmp -= I
        else:
            I = Ssurf_tmp
            Ssurf_tmp = 0.0
        Ssoil_tmp[0] += I
        I /= perlen

        # SAT flags (saturation overland flow) and head correction
        HEADSini_corr = HEADSini * 1.0
        dgwt_corr = dgwt
        if EXF_ini > 0.0:
            for l in range(nsl - 1, -1, -1):
                if Ssoil_tmp[l] == Sm[l] * Tl[l]:
                    HEADSini_corr += float(Tl[l])
                    dgwt_corr -= Tl[l]
                    SAT[l] = True
                else:
                    break

        # SPRINKLER IRRIGATION DOWN THE COLUMN (see IRR): what of it the top
        # horizon could not take fills the horizons below, top-down. After
        # the SAT flags, which stay "saturated from below".
        Rcas = None
        if IRR > 0.0 and Ssurf_tmp > 0.0:
            down = min(Ssurf_tmp, IRR)
            Rcas = np.zeros(nsl, dtype=np.float64)
            for l in range(1, nsl):
                if down <= 0.0:
                    break
                t = min(down, max(Sm[l] * Tl[l] - Ssoil_tmp[l], 0.0))
                Ssoil_tmp[l] += t
                Rcas[:l] += t          # it crossed every horizon above l
                down -= t
            Ssurf_tmp -= Rcas[0]
            I += Rcas[0] / perlen
            Rcas /= perlen

        # RUNOFF. Whatever the soil could not take leaves the cell in the
        # same step; SFR and LAK receive it and evaporate it (WP1d).
        Ro_tmp = Ssurf_tmp / perlen
        Ssurf_tmp = 0.0
        Eow_tmp = 0.0

        Rp_tmp = np.zeros(nsl, dtype=np.float64)
        Tsoil_tmpZr = np.zeros((nsl, NVEG), dtype=np.float64)
        Tsoil_tmp = np.zeros(nsl, dtype=np.float64)
        Esoil_tmp = np.zeros(nsl, dtype=np.float64)
        Ssoil_pc_tmp = np.zeros(nsl, dtype=np.float64)

        # soil layers: transpiration, evaporation, percolation
        for l in range(nsl):
            # Tsoil
            for v in range(NVEG):
                if PT[v] > 0.0 and LAIveg[v] > 0.0 and VEGarea[v] > 0.0:
                    if BotSoilLay[l] > Zr_elev[v]:
                        Tsoil_tmpZr[l, v] = self._evp(Ssoil_tmp[l], Sm[l] * Tl[l],
                                                      Sr[l] * Tl[l], PT[v], perlen)
                    elif TopSoilLay[l] > Zr_elev[v]:
                        PTc = PT[v] * (TopSoilLay[l] - Zr_elev[v]) / Tl[l]
                        Tsoil_tmpZr[l, v] = self._evp(Ssoil_tmp[l], Sm[l] * Tl[l],
                                                      Sr[l] * Tl[l], PTc, perlen)
                    Tsoil_tmp[l] += Tsoil_tmpZr[l, v] * VEGarea[v] * 0.01
                    PT[v] -= Tsoil_tmpZr[l, v]
                    if PT[v] < 0.0:
                        PT[v] = 0.0
            Ssoil_tmp[l] -= Tsoil_tmp[l] * perlen
            # Esoil (only from a dry surface)
            if Ssurf_tmp == 0.0 and PE > 0.0:
                Esoil_tmp[l] = self._evp(Ssoil_tmp[l], Sm[l] * Tl[l],
                                         Sr[l] * Tl[l], PE, perlen)
                Ssoil_tmp[l] -= Esoil_tmp[l] * perlen
                PE -= Esoil_tmp[l]
                if PE < 0.0:
                    PE = 0.0
            # Rp percolation
            if l < (nsl - 1):
                if not SAT[l + 1]:
                    Rp_tmp[l] = self._perc(Ssoil_tmp[l], Sm[l] * Tl[l],
                                           Sfc[l] * Tl[l], Ks[l], perlen)
                    if Rp_tmp[l] * perlen > (Sm[l + 1] * Tl[l + 1] - Ssoil_tmp[l + 1]):
                        Rp_tmp[l] = (Sm[l + 1] * Tl[l + 1] - Ssoil_tmp[l + 1]) / perlen
                    Ssoil_tmp[l + 1] += Rp_tmp[l] * perlen
            elif EXF_ini == 0.0:
                Rp_tmp[l] = self._perc(Ssoil_tmp[l], Sm[l] * Tl[l],
                                       Sfc[l] * Tl[l], Ks[l], perlen)
            Ssoil_tmp[l] -= Rp_tmp[l] * perlen

        if Rcas is not None:
            # the sprinkler water that went down crossed the horizons above
            # where it stopped: percolation in each layer's balance (already
            # moved, so only booked here; never below the bottom horizon)
            Rp_tmp[:nsl - 1] += Rcas[:nsl - 1]

        Ssoil_pc_tmp[:] = Ssoil_tmp / Tl[:nsl]

        # WP2 step 4: what the soil left of PE and PT is the deep unsaturated
        # zone's demand (UZF PET), per unit cell area [mm/d]
        PT_left = float(sum(PT[v] * VEGarea[v] * 0.01 for v in range(NVEG)
                            if LAIveg[v] > 0.0 and VEGarea[v] > 0.0))
        PETuzf = max(float(PE), 0.0) + PT_left
        # WP2 step 5: groundwater ET sees what remains after UZF's ACTUAL
        # uptake. Both remainders shrink in proportion, so the reduction is
        # exactly ETuzf (at most what was left) and E and T keep their ratio.
        ETuzf_used = min(max(float(ETUZF_prev), 0.0), PETuzf)
        if ETuzf_used > 0.0 and PETuzf > 0.0:
            f = 1.0 - ETuzf_used / PETuzf
            PE = PE * f
            PT = PT * f

        if GW_EVT:
            # the coupled run: what each groundwater source COULD take today
            # at the start-of-day head; MF6's EVT takes it at the head it
            # solves for. Nothing is drawn down here.
            gw = {'eg_pe': float(PE) if (Ssurf_tmp == 0.0 and PE > 0.0)
                  else 0.0,
                  'shah': self.paramEg[st], 'h0': float(HEADSini_corr),
                  'tg_rate': np.zeros(NVEG, dtype=np.float64),
                  'tg_tip': np.asarray(Zr_elev, dtype=np.float64).copy()}
            for v in Zr_elev.argsort():
                if not HEADSini_corr > Zr_elev[v]:
                    continue
                for l in range(nsl):
                    if np.isclose(Ssoil_pc_tmp[l], cMF.hnoflo):
                        continue
                    kTg = self._ktg(Ssoil_pc_tmp[l], Sr[l], Sm[l],
                                    float(np.ravel(kTg_min)[v]),
                                    float(np.ravel(kTg_max)[v]),
                                    float(np.ravel(kT_f)[v]),
                                    float(np.ravel(kT_s)[v]))
                    t = PT[v] * kTg
                    PT[v] -= t
                    gw['tg_rate'][v] += t * VEGarea[v] * 0.01
            return (Eow_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp,
                    Ssoil_tmp, Ssoil_pc_tmp, 0.0, 0.0, HEADSini_corr,
                    dgwt_corr, SAT, Rexf_tmp, I, PETuzf, REinf, gw)

        sy_tmp = float(cMF.cPROCESS.float2array(cMF.sy_actual)[cMF.outcropL[i, j] - 1, i, j])
        dgwt_corr_tmp = float(dgwt_corr) * 0.1        # mm -> cm (Shah et al. parameters)
        HEADSini_corr_tmp = float(HEADSini_corr) * 0.1

        # GW evaporation Eg, eq. 17 of Shah et al. (2007)
        Eg_tmp = 0.0
        if Ssurf_tmp == 0.0 and PE > 0.0:
            Eg_tmp, dgwt_corr_tmp, HEADSini_corr_tmp = self._eg(
                PE, dgwt_corr_tmp, HEADSini_corr_tmp, sy_tmp,
                self.paramEg[st])
        Eg_tmp = float(Eg_tmp)
        dgwt_corr_tmp *= 10.0        # cm -> mm
        HEADSini_corr_tmp *= 10.0

        # Groundwater transpiration Tg (phenomenological kTg function)
        Tg_tmp = 0.0
        # ETg is applied through the WEL package, which is always present:
        # computing it IS what the well is for.
        order = Zr_elev.argsort()
        for Zr_elev_, v, kTg_min_, kTg_max_, kT_f_, kT_s_ in zip(
                Zr_elev[order], np.arange(NVEG)[order],
                np.asarray(kTg_min, dtype=np.float64).reshape(NVEG)[order],
                np.asarray(kTg_max, dtype=np.float64).reshape(NVEG)[order],
                np.asarray(kT_f, dtype=np.float64).reshape(NVEG)[order],
                np.asarray(kT_s, dtype=np.float64).reshape(NVEG)[order]):
            if HEADSini_corr_tmp > Zr_elev_:
                for l in range(nsl):
                    if np.isclose(Ssoil_pc_tmp[l], cMF.hnoflo):
                        continue
                    kTg = self._ktg(Ssoil_pc_tmp[l], Sr[l], Sm[l], kTg_min_,
                                    kTg_max_, kT_f_, kT_s_)
                    Tg_tmp_Zr = PT[v] * kTg
                    PT[v] -= Tg_tmp_Zr
                    Tg_tmp1 = Tg_tmp_Zr * VEGarea[v] * 0.01
                    # limit Tg so the corrected head does not drop below root bottom
                    if Tg_tmp1 > 0.0 and (HEADSini_corr_tmp - Tg_tmp1 / sy_tmp) < Zr_elev_:
                        # never negative: once the head sits at the root tip a
                        # rounding residue below it would turn Tg into a source
                        Tg_tmp1 = max((HEADSini_corr_tmp - float(Zr_elev[v]))
                                      * sy_tmp, 0.0)
                    Tg_tmp += Tg_tmp1
                    dgwt_corr_tmp += Tg_tmp1 / sy_tmp
                    HEADSini_corr_tmp -= Tg_tmp1 / sy_tmp

        return (Eow_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp,
                Ssoil_tmp, Ssoil_pc_tmp, Eg_tmp, Tg_tmp, HEADSini_corr,
                dgwt_corr, SAT, Rexf_tmp, I, PETuzf, REinf, None)

    # ------------------------------------------------------------------ #

    # ------------------------------------------------------------------ #
    # Grid-agnostic cell model (Phase 2)
    #
    # MARMITES operates on a 1-D list of active soil columns ("cells").
    # For a structured grid a cell maps to (i, j) with node = i*ncol + j;
    # for a DISV grid (Phase 4) a cell maps to a single icell2d/node and
    # the (i, j) attributes become the same node. The per-stress-period
    # step() method below is the interface the MODFLOW 6 API driver
    # (Phase 3) calls: it takes heads/exfiltration per cell (ncell,) and
    # returns percolation and ETg per cell, advancing the soil state.
    # ------------------------------------------------------------------ #

    @staticmethod
    def build_cell_list(cMF):
        """Return the list of active cells as (cid, i, j, node) tuples.

        Active = cMF.outcropL > 0 (a soil column outcrops there). node is
        the flat index i*ncol + j, i.e. the future DISV node number.
        """
        cells = []
        cid = 0
        for i in range(cMF.nrow):
            for j in range(cMF.ncol):
                if cMF.outcropL[i, j] > 0:
                    cells.append((cid, i, j, i * cMF.ncol + j))
                    cid += 1
        return cells

    def build_context(self, cMF, cells, _nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop,
                      _Ssoil_ini, botm_l0, _Ks, gridSOIL, gridSOILthick, TopSoil, gridMETEO,
                      index, index_S,
                      P_veg_zoneSP, Eo_zonesSP, PT_veg_zonesSP, Pe_veg_zonesSP, PE_zonesSP,
                      gridVEGarea, LAI_veg_zonesSP, Zr, kTg_min, kTg_max, kT_f, kT_s, NVEG,
                      conv_fact, irr_yn, P_irr_zoneSP, PT_irr_zonesSP, Pe_irr_zoneSP,
                      crop_irr_SP, gridIRR, Zr_c, kTg_min_c, kTg_max_c, kT_f_c, kT_s_c,
                      geom=None):
        """Bundle all static (time-invariant) inputs into a namespace so the
        per-cell kernel and step() have a compact signature.

        geom : marmites_grid.CellGeometry or None
            Per-cell area/width provider (Phase 4). When None a
            StructuredGeometry is derived from cMF.delr/delc, reproducing the
            legacy DIS behaviour exactly; pass a VertexGeometry for DISV.
        """
        if geom is None:
            geom = _grid_module().StructuredGeometry.from_cells(cMF, cells)
        return SimpleNamespace(
            cMF=cMF, cells=cells, ncell=len(cells), geom=geom,
            _nsl=_nsl, _nslmax=_nslmax, _st=_st, _Sm=_Sm, _Sfc=_Sfc, _Sr=_Sr,
            _slprop=_slprop, _Ssoil_ini=_Ssoil_ini, botm_l0=botm_l0, _Ks=_Ks,
            gridSOIL=gridSOIL, gridSOILthick=gridSOILthick, TopSoil=TopSoil, gridMETEO=gridMETEO,
            index=index, index_S=index_S,
            P_veg_zoneSP=P_veg_zoneSP, Eo_zonesSP=Eo_zonesSP, PT_veg_zonesSP=PT_veg_zonesSP,
            Pe_veg_zonesSP=Pe_veg_zonesSP, PE_zonesSP=PE_zonesSP, gridVEGarea=gridVEGarea,
            LAI_veg_zonesSP=LAI_veg_zonesSP, Zr=Zr, kTg_min=kTg_min, kTg_max=kTg_max,
            kT_f=kT_f, kT_s=kT_s, NVEG=NVEG, conv_fact=conv_fact, irr_yn=irr_yn,
            P_irr_zoneSP=P_irr_zoneSP, PT_irr_zonesSP=PT_irr_zonesSP, Pe_irr_zoneSP=Pe_irr_zoneSP,
            crop_irr_SP=crop_irr_SP, gridIRR=gridIRR, Zr_c=Zr_c, kTg_min_c=kTg_min_c,
            kTg_max_c=kTg_max_c, kT_f_c=kT_f_c, kT_s_c=kT_s_c,
        )

    def init_state(self, ctx):
        """Per-cell carry-over state, flattened to (ncell,) / (ncell, nslmax).

        This is the (nper, ncpl) data model: state is indexed by cell id,
        not by (row, col)."""
        return SimpleNamespace(
            Ssoil_ini=np.zeros((ctx.ncell, ctx._nslmax), dtype=np.float64),
            # True once the state holds a previous run's END (a periodic
            # spin-up cycle): the first period then starts from it, not from
            # the panel's initial soil moisture
            carried=False,
        )

    def _cell_step(self, ctx, cell, n, tstart_MF, h_MF_ini_tmp, exf_MF_ini_tmp, state,
                   rejinf_cell=0.0, etuzf_cell=0.0, runon=0.0, gw_sink=None):
        """Soil water balance of one active cell for one stress period.

        Returns (MM_tmp, MM_S_tmp (nsl, nindex_S), nsl, perc_vol, etg_vol,
        petuzf_vol) and writes the next-SP state for this cell in `state`.
        ``etuzf_cell`` is the previous stress period's ACTUAL UZF ET [mm/d].

        WP4.6, THE SURFACE DESCRIPTOR (cookbook §4a). A cell may be part open
        water -- ``ctx.f_open`` = its pond (f_lake) and channel (f_stream)
        share -- and the soil column runs on the rest, f_soil. The column is
        computed exactly as before, per unit of SOIL area, and every flux it
        makes enters the cell's balance times f_soil. Over the open fraction
        there is no interception, no soil or groundwater ET, no percolation:
        the rain, and any groundwater seeping up, go straight to the water
        body as runoff (the coupler hands it to SFR / LAK), and its
        evaporation is MF6's, from the Eo forcing (iEow). Every output is
        per unit of CELL area, as before; the carried soil state stays the
        column's own.

        ``runon`` [mm/d per CELL area] is what the CRR cascade brings from
        upslope (WP5); it lands on the soil column, so per SOIL area it is
        runon / f_soil.
        """
        cid, i, j, _node = cell
        fo = getattr(ctx, 'f_open', None)
        f_open = 0.0 if fo is None else min(max(float(fo[cid]), 0.0), 1.0)
        f_soil = 1.0 - f_open
        # what reaches the column from below is per CELL area: the soil share
        # of it, over the soil area, is the same depth for exfiltration --
        # the open share goes to the water body -- but UZF's rejected
        # infiltration and ET only ever came from under the soil fraction
        col_rej = rejinf_cell / f_soil if f_soil > 0.0 else 0.0
        col_etuzf = etuzf_cell / f_soil if f_soil > 0.0 else 0.0
        col_runon = runon / f_soil if f_soil > 0.0 else 0.0
        cMF = ctx.cMF
        # unpack static context to locals (kernel body kept close to v0.3)
        _nsl, _slprop, _st, _Sm, _Sfc, _Sr, _Ks, _Ssoil_ini = (
            ctx._nsl, ctx._slprop, ctx._st, ctx._Sm, ctx._Sfc, ctx._Sr, ctx._Ks, ctx._Ssoil_ini)
        gridSOIL, gridMETEO, gridSOILthick, gridVEGarea = (
            ctx.gridSOIL, ctx.gridMETEO, ctx.gridSOILthick, ctx.gridVEGarea)
        gridIRR, TopSoil, botm_l0 = (ctx.gridIRR, ctx.TopSoil, ctx.botm_l0)
        NVEG, irr_yn, index, index_S, conv_fact = (
            ctx.NVEG, ctx.irr_yn, ctx.index, ctx.index_S, ctx.conv_fact)
        Zr, kTg_min, kTg_max, kT_f, kT_s = ctx.Zr, ctx.kTg_min, ctx.kTg_max, ctx.kT_f, ctx.kT_s
        Zr_c, kTg_min_c, kTg_max_c, kT_f_c, kT_s_c = (
            ctx.Zr_c, ctx.kTg_min_c, ctx.kTg_max_c, ctx.kT_f_c, ctx.kT_s_c)

        SOILzone_tmp = gridSOIL[i, j] - 1
        METEOzone_tmp = gridMETEO[i, j] - 1
        slprop = _slprop[SOILzone_tmp]
        nsl = _nsl[SOILzone_tmp]
        # thickness of soil layers [mm]
        Tl = np.asarray(gridSOILthick[i, j] * slprop * 1000.0, dtype=np.float64)
        # elevation of top and bottom of soil layers
        TopSoilLay = np.zeros(nsl, dtype=np.float64)
        BotSoilLay = np.zeros(nsl, dtype=np.float64)
        for l in range(nsl):
            if l == 0:
                TopSoilLay[0] = TopSoil[i, j]
                BotSoilLay[0] = TopSoil[i, j] - Tl[0]
            else:
                TopSoilLay[l] = BotSoilLay[l - 1]
                BotSoilLay[l] = TopSoilLay[l] - Tl[l]
        first = (n == 0 and not getattr(state, 'carried', False))
        if first:
            Ssoil_ini_tmp = np.asarray(_Ssoil_ini[SOILzone_tmp][:nsl], dtype=np.float64)
            perleni = 1.0
        elif n == 0:
            # a periodic spin-up cycle: the soil the previous cycle ended with
            Ssoil_ini_tmp = np.asarray(state.Ssoil_ini[cid, :nsl], dtype=np.float64)
            perleni = float(cMF.perlen[-1])
        else:
            Ssoil_ini_tmp = np.asarray(state.Ssoil_ini[cid, :nsl], dtype=np.float64)
            perleni = float(cMF.perlen[n - 1])
        IRRfield = 0
        if irr_yn == 1:
            IRRfield = int(gridIRR[i, j])
        if irr_yn == 1 and IRRfield > 0:
            NVEG_tmp = 1
            IRRfield -= 1
            LAIveg_tmp = np.ones(NVEG_tmp, dtype=np.float64)
            # crop id of this SP (bug fix vs v0.3: was an elementwise
            # `!= None` comparison on a numpy array)
            CROP_tmp = int(ctx.crop_irr_SP[IRRfield, tstart_MF])
            P_tmp = float(ctx.P_irr_zoneSP[METEOzone_tmp, IRRfield, tstart_MF])
            Pe_zonesSP_tmp = np.array([ctx.Pe_irr_zoneSP[METEOzone_tmp, IRRfield, tstart_MF]])
            PT_zonesSP_tmp = np.array([ctx.PT_irr_zonesSP[METEOzone_tmp, IRRfield, tstart_MF]])
            Zr_tmp = np.asarray(Zr_c, dtype=np.float64)
            kTg_min_tmp = np.asarray(kTg_min_c, dtype=np.float64)
            kTg_max_tmp = np.asarray(kTg_max_c, dtype=np.float64)
            kT_f_tmp = np.asarray(kT_f_c, dtype=np.float64)
            kT_s_tmp = np.asarray(kT_s_c, dtype=np.float64)
        else:
            NVEG_tmp = NVEG
            CROP_tmp = None
            P_tmp = float(ctx.P_veg_zoneSP[METEOzone_tmp][tstart_MF])
            PT_zonesSP_tmp = np.zeros(NVEG_tmp, dtype=np.float64)
            Pe_zonesSP_tmp = np.zeros(NVEG_tmp, dtype=np.float64)
            LAIveg_tmp = np.zeros(NVEG_tmp, dtype=np.float64)
            VEGarea_tmp = np.zeros(NVEG_tmp, dtype=np.float64)
            Zr_tmp = np.zeros(NVEG_tmp, dtype=np.float64)
            kTg_min_tmp = np.zeros(NVEG_tmp, dtype=np.float64)
            kTg_max_tmp = np.zeros(NVEG_tmp, dtype=np.float64)
            kT_f_tmp = np.zeros(NVEG_tmp, dtype=np.float64)
            kT_s_tmp = np.zeros(NVEG_tmp, dtype=np.float64)
            for v in range(NVEG_tmp):
                PT_zonesSP_tmp[v] = ctx.PT_veg_zonesSP[METEOzone_tmp, v, tstart_MF]
                Pe_zonesSP_tmp[v] = ctx.Pe_veg_zonesSP[METEOzone_tmp, v, tstart_MF]
                LAIveg_tmp[v] = ctx.LAI_veg_zonesSP[v, tstart_MF]
                VEGarea_tmp[v] = gridVEGarea[v, i, j]
                Zr_tmp[v] = float(Zr[v])
                kTg_min_tmp[v] = float(kTg_min[v])
                kTg_max_tmp[v] = float(kTg_max[v])
                kT_f_tmp[v] = float(kT_f[v])
                kT_s_tmp[v] = float(kT_s[v])
        PE_zonesSP_tmp = float(ctx.PE_zonesSP[METEOzone_tmp, SOILzone_tmp, tstart_MF])
        Eo_zonesSP_tmp = float(ctx.Eo_zonesSP[METEOzone_tmp][tstart_MF])
        st = _st[SOILzone_tmp]
        Sm = np.asarray(_Sm[SOILzone_tmp], dtype=np.float64)
        Sfc = np.asarray(_Sfc[SOILzone_tmp], dtype=np.float64)
        Sr = np.asarray(_Sr[SOILzone_tmp], dtype=np.float64)
        Ks = np.asarray(_Ks[SOILzone_tmp], dtype=np.float64)
        shapeFactor = 1.126847784
        # Phase 4: characteristic cell width from the grid geometry provider
        # (DIS: delr[j], legacy-identical; DISV: sqrt(cell area))
        cell_w = float(ctx.geom.width[cid])

        # vegetation patchwork: PT/Pe totals and rooting depths
        SOILarea = 100.0
        if CROP_tmp is not None:
            Zr_elev = np.array([TopSoilLay[0] - float(Zr_tmp[CROP_tmp - 1]) * 1000.0])
            VEGarea_tmp = np.array([100.0]) if CROP_tmp > 0 else np.array([0.0])
            kTg_min_tmp = np.array([kTg_min_tmp[CROP_tmp - 1]])
            kTg_max_tmp = np.array([kTg_max_tmp[CROP_tmp - 1]])
            kT_f_tmp = np.array([kT_f_tmp[CROP_tmp - 1]])
            kT_s_tmp = np.array([kT_s_tmp[CROP_tmp - 1]])
        else:
            Zr_elev = TopSoilLay[0] - Zr_tmp * 1000.0
        Pe_tot = 0.0
        PT_tot = 0.0
        for v in range(NVEG_tmp):
            if LAIveg_tmp[v] > 1.0E-5:
                Pe_tot += Pe_zonesSP_tmp[v] * VEGarea_tmp[v] * 0.01
                PT_tot += PT_zonesSP_tmp[v] * VEGarea_tmp[v] * 0.01
                SOILarea -= VEGarea_tmp[v]
        Pe_tot += P_tmp * SOILarea * 0.01
        INTER_tot = P_tmp - Pe_tot
        # SPRINKLER IRRIGATION down the whole column (soil.irr_infiltration =
        # 'column', ctx.irr_column): an irrigated cell's input is rain AND
        # irrigation (MMsurf adds them), its zone's plain series the rain, so
        # the irrigation is the difference -- its share of what reaches the
        # soil, interception shared pro rata
        irr_col = 0.0
        if CROP_tmp is not None and getattr(ctx, 'irr_column', False) \
                and P_tmp > 0.0:
            rain = float(ctx.P_veg_zoneSP[METEOzone_tmp][tstart_MF])
            if P_tmp > rain:
                irr_col = Pe_tot * (P_tmp - rain) / P_tmp
        PE_tot = PE_zonesSP_tmp * SOILarea * 0.01
        # A TYPE THE DEMAND DOES NOT COUNT TRANSPIRES NOTHING. Dormant (LAI
        # <= 1e-5: La Mata's grass 101 days a year, lai_dry = 0) its PT is
        # left out of PT_tot and its area evaporates as bare soil in PE_tot
        # -- yet flux() let it draw groundwater through Tg on the same area:
        # ET above PET, and the area used twice.
        PT_flux = np.where(np.asarray(LAIveg_tmp) > 1.0E-5,
                           np.asarray(PT_zonesSP_tmp, dtype=np.float64), 0.0)
        # dry cell handling
        # legacy (NWT): dry cells carry the hdry sentinel value;
        # MF6 (cMF.hdry is None): no sentinel exists -- a cell is dry when
        # the head drops below the layer bottom (Newton keeps it active)
        if cMF.hdry is not None:
            dry = np.abs(h_MF_ini_tmp - cMF.hdry) < 1.0E-5
        else:
            dry = h_MF_ini_tmp < botm_l0[i, j]
        if dry:
            HEADSini_drycell = botm_l0[i, j] * 1000.0
        else:
            HEADSini_drycell = h_MF_ini_tmp * 1000.0
        # dgwt and uzthick
        if exf_MF_ini_tmp <= 0.0:
            dgwt = TopSoilLay[0] - HEADSini_drycell
        else:
            dgwt = float(np.sum(Tl[:nsl]))
            HEADSini_drycell = BotSoilLay[nsl - 1]
        uzthick = BotSoilLay[nsl - 1] - HEADSini_drycell
        # for the first SP, Ssoil_ini is in % and has to be converted to mm
        if first:
            Ssoil_ini_tmp = Ssoil_ini_tmp * Tl[:nsl]

        # MAIN SUB-ROUTINE fluxes
        (Eow_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp, Ssoil_tmp,
         Ssoil_pc_tmp, Eg_tmp, Tg_tmp, HEADSini_MM, dgwt_tmp, SAT_tmp, Rexf_tmp,
         I, PETuzf, REinf, GW) = self.flux(cMF, perleni, Pe_tot, PT_flux,
                        PE_zonesSP_tmp * SOILarea * 0.01, Zr_elev,
                        VEGarea_tmp, HEADSini_drycell, TopSoilLay, BotSoilLay,
                        Tl, nsl, Sm, Sfc, Sr, Ks, Ssoil_ini_tmp,
                        exf_MF_ini_tmp, dgwt, st, i, j, n,
                        kTg_min_tmp, kTg_max_tmp, kT_f_tmp, kT_s_tmp,
                        NVEG_tmp, LAIveg_tmp, REJINF_ini=col_rej,
                        ETUZF_prev=col_etuzf, RUNON=col_runon, IRR=irr_col,
                        GW_EVT=bool(getattr(ctx, 'gw_evt', False)))
        if GW is not None and gw_sink is not None:
            # the coupled run: per unit of CELL area and in m, as MF6's
            # EVT takes them -- the column's rates times f_soil
            gw_sink[cid] = {
                'eg_pe': f_soil * GW['eg_pe'] / 1000.0,
                'shah': GW['shah'], 'h0': GW['h0'] / 1000.0,
                'land': float(TopSoilLay[0]) / 1000.0,
                'tg_rate': f_soil * np.asarray(GW['tg_rate']) / 1000.0,
                'tg_tip': np.asarray(GW['tg_tip']) / 1000.0}
        Ssoil_pc_tot = float(np.sum(Ssoil_pc_tmp)) / nsl
        perc = Rp_tmp[-1]
        ETg = Eg_tmp + Tg_tmp
        # WP1d: no surface reservoir, so nothing is carried between stress
        # periods. iSsurf and idSsurf keep their slots in the flux index and
        # are structurally zero; open-water evaporation is now applied by SFR
        # and LAK, from the same Eo forcing.
        dSsurf = 0.0

        # water mass balance (MB) in the soil zone. The water entering from
        # below is groundwater exfiltration AND the rejected infiltration UZF
        # returned (flux adds both at the base); counting the first alone
        # left the soil balance -- and the Sankey -- short by the second
        # (-5.6 %, 36 mm/yr, on 2026-09-24).
        rej_in = abs(float(col_rej))
        exf_in = exf_MF_ini_tmp / cMF.perlen[n] + rej_in
        Esoil_MB = float(np.sum(Esoil_tmp))
        Tsoil_MB = float(np.sum(Tsoil_tmp))
        dSsoil = (Ssoil_tmp - Ssoil_ini_tmp) / cMF.perlen[n]
        dSsoil_tot = float(np.sum(dSsoil))
        ETsoil_tot = Esoil_MB + Tsoil_MB

        MB_l = np.zeros(nsl, dtype=np.float64)
        # Eq. 1: dSsurf/dt = Pe + Exf_g1 - I - Eow - Ro, with dSsurf and
        # Eow now identically zero -- the surface term closes as
        # Pe + Exf_g1 = I + Ro
        # (Exf_g1 = Rexf[0], which now carries any rejected infiltration that
        #  saturated the soil column from below)
        # ... + the CRR run-on (WP5) when the cascade brings any
        if col_runon > 0.0:
            MBsurf = (Pe_tot + col_runon + Rexf_tmp[0]) - (
                Eow_tmp + Ro_tmp + I + dSsurf)
        else:
            MBsurf = (Pe_tot + Rexf_tmp[0]) - (Eow_tmp + Ro_tmp + I + dSsurf)
        if nsl > 1:
            # surficial soil layer
            MB_l[0] = (I + Rexf_tmp[1]) - (Rp_tmp[0] + Rexf_tmp[0]
                                           + Esoil_tmp[0] + Tsoil_tmp[0] + dSsoil[0])
            # intermediate soil layers
            for l in range(1, nsl - 1):
                MB_l[l] = (Rp_tmp[l - 1] + Rexf_tmp[l + 1]) - (
                    Rp_tmp[l] + Rexf_tmp[l] + Esoil_tmp[l] + Tsoil_tmp[l] + dSsoil[l])
            # last soil layer
            l = nsl - 1
            MB_l[l] = (Rp_tmp[l - 1] + exf_in) - (
                Rp_tmp[l] + Rexf_tmp[l] + Esoil_tmp[l] + Tsoil_tmp[l] + dSsoil[l])
        else:
            MB_l[0] = (I + exf_in) - (
                Rp_tmp[0] + Rexf_tmp[0] + Esoil_tmp[0] + Tsoil_tmp[0] + dSsoil[0])
        # total mass balance for the soil
        MB = (I + exf_in) - (
            Rexf_tmp[0] + Esoil_MB + Tsoil_MB + dSsoil_tot + Rp_tmp[-1])

        # export arrays, per unit of CELL area: the column's fluxes times
        # f_soil, and the open fraction's rain and seepage as runoff to the
        # water body (WP4.6). The states (soil moisture, depth to water) are
        # the column's own.
        fs = f_soil
        ro_open = f_open * (P_tmp + max(exf_MF_ini_tmp, 0.0) / cMF.perlen[n])
        if runon > 0.0 and f_soil <= 0.0:
            ro_open += runon        # no column to land on: it runs through
        _vals = {'iP': P_tmp, 'iPT': fs * PT_tot, 'iPE': fs * PE_tot,
                 'iPe': fs * Pe_tot + f_open * P_tmp,
                 'iSsurf': fs * Ssurf_tmp, 'iRo': fs * Ro_tmp + ro_open,
                 'iEXFg': exf_MF_ini_tmp,
                 'iEow': fs * Eow_tmp, 'iMB': fs * MB, 'iEi': fs * INTER_tot,
                 'iEo': Eo_zonesSP_tmp, 'iEg': fs * Eg_tmp, 'iTg': fs * Tg_tmp,
                 'idSsurf': fs * dSsurf, 'iETg': fs * ETg,
                 'iETsoil': fs * ETsoil_tot,
                 'iSsoil_pc': Ssoil_pc_tot, 'idSsoil': fs * dSsoil_tot,
                 'iperc': fs * perc, 'ihcorr': HEADSini_MM * 0.001,
                 'idgwt': -dgwt_tmp * 0.001, 'iuzthick': uzthick * 0.001,
                 'iI': fs * I, 'iMBsurf': fs * MBsurf,
                 # WP2: UZF's demand; its ACTUAL and the total are written by
                 # the coupler once MF6 has solved (iETuzf, iETtot)
                 'iPETuzf': fs * PETuzf, 'iETuzf': 0.0,
                 'iETtot': fs * (INTER_tot + Eow_tmp + ETsoil_tot + ETg),
                 'iRejInf': fs * rej_in}
        if runon > 0.0:
            # WP5: what the cascade brought, and how much of it the soil
            # took (both already inside iI / iRo / iMBsurf above)
            _vals['iRunon'] = fs * col_runon
            _vals['iReinf'] = fs * REinf
        MM_tmp = np.zeros(len(index), dtype=np.float64)
        for _k, _v in _vals.items():
            if _k in index:
                MM_tmp[index[_k]] = float(np.asarray(_v).ravel()[0]
                                          if np.ndim(_v) else _v)
        MM_S_tmp = np.zeros([nsl, len(index_S)], dtype=np.float64)
        for l in range(nsl):
            MM_S_tmp[l, :] = [Esoil_tmp[l], Tsoil_tmp[l], Ssoil_pc_tmp[l],
                              Rp_tmp[l], Rexf_tmp[l], dSsoil[l], Ssoil_tmp[l],
                              SAT_tmp[l], MB_l[l]]
        # write next-SP state for this cell (percolation-driven soil storage)
        # -- the COLUMN's own, before it is expressed per cell area
        state.Ssoil_ini[cid, :nsl] = MM_S_tmp[:nsl, index_S.get('iSsoil')]
        if f_open > 0.0:
            for k in ('iEsoil', 'iTsoil', 'iRsoil', 'iExf', 'idSsoil_s',
                      'iSsoil', 'iMB_s'):
                if k in index_S:
                    MM_S_tmp[:, index_S[k]] *= fs

        # volumetric recharge (UZF finf) and groundwater ET (WEL) rates: from
        # under the soil fraction only -- a pond cell percolates nothing, and
        # the lake exchanges with the aquifer through its own bed instead
        perc_vol = MM_S_tmp[nsl - 1, index_S.get('iRsoil')] / conv_fact
        etg_vol = MM_tmp[index.get('iETg')] / conv_fact
        petuzf_vol = fs * PETuzf / conv_fact              # m/d, UZF PET
        return MM_tmp, MM_S_tmp, nsl, perc_vol, etg_vol, petuzf_vol

    def step(self, ctx, n, tstart_MF, heads_cell, exf_cell, state, rejinf_cell=None,
             etuzf_cell=None, crr=None):
        """Advance the soil water balance one stress period over all cells.

        Parameters
        ----------
        ctx : namespace from build_context()
        n : stress-period index
        tstart_MF : first MODFLOW time step of this SP (index into zone series)
        heads_cell : (ncell,) representative head [m] per active cell (from MF)
        exf_cell : (ncell,) exfiltration into the soil [mm/d] per active cell,
            sign convention: positive = seepage up into the soil column
        state : namespace from init_state(), mutated in place

        Returns dict with per-cell arrays: 'MM' (ncell, nindex),
        'MM_S' (ncell, nslmax, nindex_S), 'perc' (ncell,), 'etg' (ncell,).
        This is grid-agnostic; the Phase-3 API driver scatters perc/etg to
        MODFLOW node arrays, the file-based driver scatters to structured h5.

        crr : marmites_crr.CascadeNetwork or None (WP5). With one, runoff is
            routed downslope cell to cell: the cells are visited in its
            descending-elevation order (the cell list does not move), each
            receiving the run-on of the cells above it, and two more entries
            are returned -- 'ro_deliver', the runoff reaching each stream or
            pond cell [mm/d per cell area], what the coupler writes into SFR
            and LAK instead of iRo; and 'crr', the pass's volumes [m3/d]
            (marmites_crr.CascadeNetwork.route). The runoff the cascade
            evaporates is iEcrr, booked at the cell it left and counted in
            iETtot. Without one, nothing changes -- not a bit.
        """
        nindex = len(ctx.index)
        nindex_S = len(ctx.index_S)
        # DOUBLE PRECISION: the coupler checks ETsoil + ETuzf + ETg <= PE + PT
        # per cell to 1e-6 mm/d, and UZF's demand is PETuzf - ETg -- in
        # float32 each term of a 4 mm/d balance is off by ~5e-7 mm/d
        MM_cells = np.zeros((ctx.ncell, nindex), dtype=np.float64)
        MM_S_cells = np.zeros((ctx.ncell, ctx._nslmax, nindex_S), dtype=np.float64)
        perc_cell = np.zeros(ctx.ncell, dtype=np.float64)
        etg_cell = np.zeros(ctx.ncell, dtype=np.float64)
        rej = (np.zeros(ctx.ncell) if rejinf_cell is None
               else np.asarray(rejinf_cell, dtype=np.float64))
        # WP2: the previous stress period's ACTUAL UZF ET per cell [mm/d]
        etu = (np.zeros(ctx.ncell) if etuzf_cell is None
               else np.asarray(etuzf_cell, dtype=np.float64))
        petuzf_cell = np.zeros(ctx.ncell, dtype=np.float64)

        # the coupled run (MF6Coupler sets ctx.gw_evt): each cell's
        # groundwater-ET potential, for the day's EVT curves
        gw_sink = {} if getattr(ctx, 'gw_evt', False) else None

        def solve(cell, runon=0.0):
            cid = cell[0]
            MM_tmp, MM_S_tmp, nsl, perc_vol, etg_vol, petuzf_vol = self._cell_step(
                ctx, cell, n, tstart_MF, float(heads_cell[cid]), float(exf_cell[cid]), state,
                rejinf_cell=float(rej[cid]), etuzf_cell=float(etu[cid]), runon=runon,
                gw_sink=gw_sink)
            MM_cells[cid, :] = MM_tmp
            MM_S_cells[cid, :nsl, :] = MM_S_tmp
            perc_cell[cid] = perc_vol
            etg_cell[cid] = etg_vol
            petuzf_cell[cid] = petuzf_vol
            return MM_tmp

        out = {'MM': MM_cells, 'MM_S': MM_S_cells, 'perc': perc_cell,
               'etg': etg_cell, 'petuzf': petuzf_cell}
        if crr is None:
            for cell in ctx.cells:
                solve(cell)
            if gw_sink is not None:
                out['gw'] = self._gw_arrays(ctx, gw_sink)
            return out

        # WP5 -- the cascade. It addresses cells by their place in the list,
        # which is their id (build_cell_list numbers them in order).
        iro = ctx.index['iRo']
        area = np.asarray(ctx.geom.area, dtype=np.float64)
        conv = float(ctx.conv_fact)
        res = crr.route(area, lambda k, runon: solve(ctx.cells[k], runon)[iro],
                        conv=conv)
        ecrr = res.ecrr / area * conv                       # mm/d, cell area
        if 'iEcrr' in ctx.index:
            MM_cells[:, ctx.index['iEcrr']] = ecrr
        if 'iETtot' in ctx.index:
            MM_cells[:, ctx.index['iETtot']] += ecrr
        out['ro_deliver'] = res.deliver / area * conv
        out['crr'] = res
        if gw_sink is not None:
            out['gw'] = self._gw_arrays(ctx, gw_sink)
        return out

    @staticmethod
    def _gw_arrays(ctx, gw_sink):
        """The coupled run: each cell's groundwater-ET potential as
        arrays in cell order -- Eg's potential, the land surface and the
        start-of-day head [m, m/d per cell area] with the soil's Shah
        parameters, and per vegetation type Tg's rate and root tip (padded
        to the widest cell: an irrigated crop is one type)."""
        n = ctx.ncell
        nv = max((len(g['tg_rate']) for g in gw_sink.values()), default=1)
        a = {'eg_pe': np.zeros(n), 'land': np.zeros(n), 'h0': np.zeros(n),
             'shah': [None] * n, 'tg_rate': np.zeros((n, nv)),
             'tg_tip': np.full((n, nv), np.nan)}
        for cid, g in gw_sink.items():
            a['eg_pe'][cid] = g['eg_pe']
            a['land'][cid] = g['land']
            a['h0'][cid] = g['h0']
            a['shah'][cid] = g['shah']
            k = len(g['tg_rate'])
            a['tg_rate'][cid, :k] = g['tg_rate']
            a['tg_tip'][cid, :k] = g['tg_tip']
        return a

    def runMMsoil(self, _nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _Ssoil_ini, botm_l0, _Ks,
                  gridSOIL, gridSOILthick, TopSoil, gridMETEO,
                  index, index_S,
                  P_veg_zoneSP, Eo_zonesSP, PT_veg_zonesSP, Pe_veg_zonesSP, PE_zonesSP, gridVEGarea,
                  LAI_veg_zonesSP, Zr, kTg_min, kTg_max, kT_f, kT_s, NVEG,
                  cMF, conv_fact, h5_MF, h5_MM, irr_yn,
                  P_irr_zoneSP=None, PT_irr_zonesSP=None, Pe_irr_zoneSP=None,
                  crop_irr_SP=None, gridIRR=None,
                  Zr_c=None, kTg_min_c=None, kTg_max_c=None, kT_f_c=None, kT_s_c=None,
                  verbose=0, report=None, report_fn=None, stdout=None):
        """Whole-run soil water balance (file-based coupling).

        Thin driver over the grid-agnostic step(): reads heads/exfiltration
        from the MODFLOW HDF5, calls step() per stress period, scatters the
        per-cell results back to the structured (nper, nrow, ncol) HDF5
        datasets consumed by the export/plot code. The Phase-3 MODFLOW 6 API
        driver calls step() directly instead, with heads passed in memory.
        """
        cells = self.build_cell_list(cMF)
        ctx = self.build_context(
            cMF, cells, _nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _Ssoil_ini, botm_l0, _Ks,
            gridSOIL, gridSOILthick, TopSoil, gridMETEO, index, index_S,
            P_veg_zoneSP, Eo_zonesSP, PT_veg_zonesSP, Pe_veg_zonesSP, PE_zonesSP, gridVEGarea,
            LAI_veg_zonesSP, Zr, kTg_min, kTg_max, kT_f, kT_s, NVEG, conv_fact, irr_yn,
            P_irr_zoneSP, PT_irr_zonesSP, Pe_irr_zoneSP, crop_irr_SP, gridIRR,
            Zr_c, kTg_min_c, kTg_max_c, kT_f_c, kT_s_c)
        state = self.init_state(ctx)

        try:
            h_MF_ini = h5_MF['heads4MM'][:, :, :]
            h_MF_ini_mem = 'fast'
        except (MemoryError, OSError):
            h_MF_ini = None
            h_MF_ini_mem = 'slow'
            print('\nRAM memory too small compared to the size of the heads array -> slow computing.')
        exf_MF_ini = None
        exf_MF_ini_mem = 'slow'
        if cMF.uzf_yn == 1:
            try:
                # Phase 4: kept volumetric here; converted to mm/d per cell
                # below using the grid geometry (supports refined/DISV grids;
                # the legacy code assumed a single uniform cell area)
                exf_MF_ini = h5_MF['exf4MM'][:, :, :]
                exf_MF_ini_mem = 'fast'
            except (MemoryError, OSError):
                print('\nRAM memory too small compared to the size of the exfiltration array -> slow computing.')

        i_arr = np.array([c[1] for c in cells], dtype=int)
        j_arr = np.array([c[2] for c in cells], dtype=int)
        area_arr = np.asarray(ctx.geom.area, dtype=np.float64)   # per-cell [L^2]
        tstart_MM = 0
        tstart_MF = 0
        for n in range(cMF.nper):
            if n > 0:
                tstart_MM += cMF.perlen[n - 1]
                tstart_MF += cMF.nstp[n - 1]
            tend_MM = tstart_MM + cMF.perlen[n]
            if h_MF_ini_mem == 'slow':
                h_MF_ini = h5_MF['heads4MM'][tstart_MF:tstart_MF + cMF.nstp[n], :, :]
                h_row = 0
            else:
                h_row = tstart_MF
            # gather per-cell heads/exfiltration for this SP (ncell,)
            heads_cell = h_MF_ini[h_row, i_arr, j_arr].astype(np.float64)
            if cMF.uzf_yn == 1:
                if exf_MF_ini_mem == 'slow':
                    exf_sp = h5_MF['exf4MM'][tstart_MF:tstart_MF + cMF.nstp[n], :, :]
                    exf_raw = exf_sp[0, i_arr, j_arr].astype(np.float64)
                else:
                    exf_raw = exf_MF_ini[h_row, i_arr, j_arr].astype(np.float64)
                # volumetric -> mm/d, positive into the soil, per-cell area
                exf_cell = -exf_raw * conv_fact / area_arr
            else:
                exf_cell = np.zeros(ctx.ncell, dtype=np.float64)

            out = self.step(ctx, n, tstart_MF, heads_cell, exf_cell, state)

            # scatter per-cell results to structured HDF5 datasets
            MM = np.full((cMF.perlen[n], cMF.nrow, cMF.ncol, len(index)), cMF.hnoflo, dtype=np.float32)
            MM_S = np.full((cMF.perlen[n], cMF.nrow, cMF.ncol, _nslmax, len(index_S)), cMF.hnoflo, dtype=np.float32)
            MM_perc_MF = np.zeros((cMF.nrow, cMF.ncol), dtype=np.float32)
            MM_wel_MF = np.zeros((cMF.nrow, cMF.ncol), dtype=np.float32)
            # broadcast the SP value over all days of the SP (as in v0.3)
            MM[:, i_arr, j_arr, :] = out['MM'][np.newaxis, :, :]
            MM_S[:, i_arr, j_arr, :, :] = out['MM_S'][np.newaxis, :, :, :]
            MM_perc_MF[i_arr, j_arr] = out['perc']
            MM_wel_MF[i_arr, j_arr] = out['etg']

            h5_MM['MM'][tstart_MM:tend_MM, :, :, :] = MM
            h5_MM['MM_S'][tstart_MM:tend_MM, :, :, :, :] = MM_S
            h5_MM['perc'][n, :, :] = MM_perc_MF
            h5_MM['ETg'][n, :, :] = MM_wel_MF
            if n > 0 and n % 100 == 0:
                print('Processed data up to stress period %d of %d' % (n, cMF.nper))
                if verbose == 0:
                    sys.stdout = stdout
                    report.close()
                    stdout = sys.stdout
                    report = open(report_fn, 'a')
                    sys.stdout = report
        h5_MM.close()


class SATFLOW:
    """
    SATFLOW: SATurated FLOW
    Water-level fluctuations as a function of recharge.

    INPUT PARAMETERS: hi, h0 (base level), RC (recession constant),
    STO (storage capacity); STATE: Rg (daily gross recharge).
    OUTPUT: h (daily water level).
    """

    def runSATFLOW(self, Rg, hi, h0, RC, STO):
        Rg = np.asarray(Rg, dtype=np.float64)
        h1 = np.zeros(len(Rg), dtype=np.float64)
        h1[0] = hi * 1000.0 + Rg[0] / STO - hi * 1000.0 / RC
        for t in range(1, len(Rg)):
            h1[t] = h1[t - 1] + Rg[t] / STO - h1[t - 1] / RC
        return (h1 + h0 * 1000.0) * 0.001


if __name__ == '__main__':
    print('\nWARNING!\nStart MARMITES-MODFLOW models using the script startMARMITES_v3.py\n')

# EOF
