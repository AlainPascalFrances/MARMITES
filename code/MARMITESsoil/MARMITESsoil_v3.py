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
            Ssurf_max  max. surface (ponding) storage [mm]
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
        Se = (s_tmp - Sr) / (Sm - Sr)  # effective saturation
        if pet * Se * perlen > (s_tmp - Sr):
            return (s_tmp - Sr) / perlen
        return pet * Se

    # ------------------------------------------------------------------ #

    def flux(self, cMF, perleni, Pe, PT, PE, Eosurf_max, Zr_elev, VEGarea,
             HEADSini, TopSoilLay, BotSoilLay, Tl, nsl, Sm, Sfc, Sr, Ks,
             Ssurf_max, Ssoil_ini, Ssurf_ini, EXF_ini, dgwt, st, i, j, n,
             kTg_min, kTg_max, kT_f, kT_s, NVEG, LAIveg, REJINF_ini=0.0):
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
            store, and only the excess above Ssurf_max becomes runoff.

            Zero for the uncoupled/file-based path, so legacy behaviour is
            unchanged.
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
        Eosurf_max = float(Eosurf_max)
        Ssurf_max = float(Ssurf_max)
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

        # Surface storage
        Ssurf_tmp = Pe + float(Ssurf_ini)

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

        # SURFACE storage, runoff Ro, open-water evaporation Eow
        if Ssurf_tmp > Ssurf_max:
            Ro_tmp = (Ssurf_tmp - Ssurf_max) / perlen
            Ssurf_tmp = Ssurf_max
        else:
            Ro_tmp = 0.0
        if Ssurf_tmp > Eosurf_max:
            Eow_tmp = Eosurf_max
            Ssurf_tmp -= Eosurf_max * perlen
        else:
            Eow_tmp = Ssurf_tmp / perlen
            Ssurf_tmp = 0.0

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

        Ssoil_pc_tmp[:] = Ssoil_tmp / Tl[:nsl]

        sy_tmp = float(cMF.cPROCESS.float2array(cMF.sy_actual)[cMF.outcropL[i, j] - 1, i, j])
        dgwt_corr_tmp = float(dgwt_corr) * 0.1        # mm -> cm (Shah et al. parameters)
        HEADSini_corr_tmp = float(HEADSini_corr) * 0.1

        # GW evaporation Eg, eq. 17 of Shah et al. (2007)
        Eg_tmp = 0.0
        if cMF.wel_yn == 1 and Ssurf_tmp == 0.0 and PE > 0.0:
            p = self.paramEg[st]
            y0, b, dll, ext_d = p['y0'], p['b'], p['dll'], p['ext_d']
            if dgwt_corr_tmp <= dll:
                Eg_tmp = PE
            elif dgwt_corr_tmp < ext_d:
                Eg_tmp = PE * (y0 + np.exp(-b * (dgwt_corr_tmp - dll)))
            else:
                Eg_tmp = 0.0
            if Eg_tmp > 0.0:
                # water-table drawdown feedback (0.1: mm -> cm)
                if (dgwt_corr_tmp + 0.1 * Eg_tmp / sy_tmp) > ext_d:
                    Eg_tmp = 10.0 * (ext_d - dgwt_corr_tmp) * sy_tmp
                    dgwt_corr_tmp = ext_d
                    HEADSini_corr_tmp -= ext_d
                else:
                    dgwt_corr_tmp += 0.1 * Eg_tmp / sy_tmp
                    HEADSini_corr_tmp -= 0.1 * Eg_tmp / sy_tmp
        Eg_tmp = float(Eg_tmp)
        dgwt_corr_tmp *= 10.0        # cm -> mm
        HEADSini_corr_tmp *= 10.0

        # Groundwater transpiration Tg (phenomenological kTg function)
        Tg_tmp = 0.0
        if cMF.wel_yn == 1:
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
                        if Ssoil_pc_tmp[l] > Sr[l]:
                            if Ssoil_pc_tmp[l] < Sm[l]:
                                Ssoil_norm = (Ssoil_pc_tmp[l] - Sr[l]) / (Sm[l] - Sr[l])
                                kTg = kTg_max_ - (kTg_max_ - kTg_min_) / (
                                    1.0 + np.exp((Ssoil_norm - kT_f_) / kT_s_))
                            else:
                                if Ssoil_pc_tmp[l] > Sm[l]:
                                    print('WARNING!\nComputing of Tg: soil moisture higher than porosity!'
                                          '\nSoil moisture = %.4f, phi = %.4f' % (Ssoil_pc_tmp[l], Sm[l]))
                                kTg = kTg_max_
                        else:
                            if Ssoil_pc_tmp[l] < Sr[l]:
                                print('WARNING!\nComputing of Tg: soil moisture lower than wilting point!'
                                      '\nSoil moisture = %.4f, WP = %.4f' % (Ssoil_pc_tmp[l], Sr[l]))
                            kTg = kTg_min_
                        Tg_tmp_Zr = PT[v] * kTg
                        PT[v] -= Tg_tmp_Zr
                        Tg_tmp1 = Tg_tmp_Zr * VEGarea[v] * 0.01
                        # limit Tg so the corrected head does not drop below root bottom
                        if Tg_tmp1 > 0.0 and (HEADSini_corr_tmp - Tg_tmp1 / sy_tmp) < Zr_elev_:
                            Tg_tmp1 = (HEADSini_corr_tmp - float(Zr_elev[v])) * sy_tmp
                        Tg_tmp += Tg_tmp1
                        dgwt_corr_tmp += Tg_tmp1 / sy_tmp
                        HEADSini_corr_tmp -= Tg_tmp1 / sy_tmp

        return (Eow_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp,
                Ssoil_tmp, Ssoil_pc_tmp, Eg_tmp, Tg_tmp, HEADSini_corr,
                dgwt_corr, SAT, Rexf_tmp, I)

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
                      index, index_S, gridSsurfhmax, gridSsurfw,
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
            index=index, index_S=index_S, gridSsurfhmax=gridSsurfhmax, gridSsurfw=gridSsurfw,
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
            Ssurf_ini=np.zeros(ctx.ncell, dtype=np.float64),
        )

    def _cell_step(self, ctx, cell, n, tstart_MF, h_MF_ini_tmp, exf_MF_ini_tmp, state,
                   rejinf_cell=0.0):
        """Soil water balance of one active cell for one stress period.

        Returns (MM_tmp list, MM_S_tmp (nsl, nindex_S), nsl, perc_vol, etg_vol)
        and writes the next-SP state for this cell in `state`.
        """
        cid, i, j, _node = cell
        cMF = ctx.cMF
        # unpack static context to locals (kernel body kept close to v0.3)
        _nsl, _slprop, _st, _Sm, _Sfc, _Sr, _Ks, _Ssoil_ini = (
            ctx._nsl, ctx._slprop, ctx._st, ctx._Sm, ctx._Sfc, ctx._Sr, ctx._Ks, ctx._Ssoil_ini)
        gridSOIL, gridMETEO, gridSOILthick, gridVEGarea = (
            ctx.gridSOIL, ctx.gridMETEO, ctx.gridSOILthick, ctx.gridVEGarea)
        gridSsurfhmax, gridSsurfw, gridIRR, TopSoil, botm_l0 = (
            ctx.gridSsurfhmax, ctx.gridSsurfw, ctx.gridIRR, ctx.TopSoil, ctx.botm_l0)
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
        if n == 0:
            Ssoil_ini_tmp = np.asarray(_Ssoil_ini[SOILzone_tmp][:nsl], dtype=np.float64)
            Ssurf_ini_tmp = 0.0
            perleni = 1.0
        else:
            Ssoil_ini_tmp = np.asarray(state.Ssoil_ini[cid, :nsl], dtype=np.float64)
            Ssurf_ini_tmp = float(state.Ssurf_ini[cid])
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
        Ssurf_max = (np.power(cell_w, 3) * gridSsurfhmax[i, j] * gridSsurfw[i, j]
                     * shapeFactor / np.power(100.0, 2) / 10.0)
        Eosurf_max = (cell_w * gridSsurfw[i, j] * shapeFactor
                      * Eo_zonesSP_tmp / np.power(100.0, 2))

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
        PE_tot = PE_zonesSP_tmp * SOILarea * 0.01
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
        if n == 0:
            Ssoil_ini_tmp = Ssoil_ini_tmp * Tl[:nsl]

        # MAIN SUB-ROUTINE fluxes
        (Eow_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp, Ssoil_tmp,
         Ssoil_pc_tmp, Eg_tmp, Tg_tmp, HEADSini_MM, dgwt_tmp, SAT_tmp, Rexf_tmp,
         I) = self.flux(cMF, perleni, Pe_tot, PT_zonesSP_tmp,
                        PE_zonesSP_tmp * SOILarea * 0.01, Eosurf_max, Zr_elev,
                        VEGarea_tmp, HEADSini_drycell, TopSoilLay, BotSoilLay,
                        Tl, nsl, Sm, Sfc, Sr, Ks, Ssurf_max, Ssoil_ini_tmp,
                        Ssurf_ini_tmp, exf_MF_ini_tmp, dgwt, st, i, j, n,
                        kTg_min_tmp, kTg_max_tmp, kT_f_tmp, kT_s_tmp,
                        NVEG_tmp, LAIveg_tmp, REJINF_ini=rejinf_cell)
        Ssoil_pc_tot = float(np.sum(Ssoil_pc_tmp)) / nsl
        perc = Rp_tmp[-1]
        ETg = Eg_tmp + Tg_tmp
        dSsurf = (Ssurf_tmp - Ssurf_ini_tmp) / cMF.perlen[n]

        # water mass balance (MB) in the soil zone
        Esoil_MB = float(np.sum(Esoil_tmp))
        Tsoil_MB = float(np.sum(Tsoil_tmp))
        dSsoil = (Ssoil_tmp - Ssoil_ini_tmp) / cMF.perlen[n]
        dSsoil_tot = float(np.sum(dSsoil))
        ETsoil_tot = Esoil_MB + Tsoil_MB

        MB_l = np.zeros(nsl, dtype=np.float64)
        # Eq. 1: dSsurf/dt = Pe + Exf_g1 - I - Eow - Ro
        # (Exf_g1 = Rexf[0], which now carries any rejected infiltration that
        #  saturated the soil column from below)
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
            MB_l[l] = (Rp_tmp[l - 1] + exf_MF_ini_tmp / cMF.perlen[n]) - (
                Rp_tmp[l] + Rexf_tmp[l] + Esoil_tmp[l] + Tsoil_tmp[l] + dSsoil[l])
        else:
            MB_l[0] = (I + exf_MF_ini_tmp / cMF.perlen[n]) - (
                Rp_tmp[0] + Rexf_tmp[0] + Esoil_tmp[0] + Tsoil_tmp[0] + dSsoil[0])
        # total mass balance for the soil
        MB = (I + exf_MF_ini_tmp / cMF.perlen[n]) - (
            Rexf_tmp[0] + Esoil_MB + Tsoil_MB + dSsoil_tot + Rp_tmp[-1])

        # export arrays
        MM_tmp = [P_tmp, PT_tot, PE_tot, Pe_tot, Ssurf_tmp, Ro_tmp,
                  exf_MF_ini_tmp, Eow_tmp, MB, INTER_tot, Eo_zonesSP_tmp,
                  Eg_tmp, Tg_tmp, dSsurf, ETg, ETsoil_tot, Ssoil_pc_tot,
                  dSsoil_tot, perc, HEADSini_MM * 0.001, -dgwt_tmp * 0.001,
                  uzthick * 0.001, I, MBsurf]
        MM_S_tmp = np.zeros([nsl, len(index_S)], dtype=np.float32)
        for l in range(nsl):
            MM_S_tmp[l, :] = [Esoil_tmp[l], Tsoil_tmp[l], Ssoil_pc_tmp[l],
                              Rp_tmp[l], Rexf_tmp[l], dSsoil[l], Ssoil_tmp[l],
                              SAT_tmp[l], MB_l[l]]

        # volumetric recharge (UZF finf) and groundwater ET (WEL) rates
        perc_vol = MM_S_tmp[nsl - 1, index_S.get('iRsoil')] / conv_fact
        etg_vol = MM_tmp[index.get('iETg')] / conv_fact
        # write next-SP state for this cell (percolation-driven soil storage)
        state.Ssoil_ini[cid, :nsl] = MM_S_tmp[:nsl, index_S.get('iSsoil')]
        state.Ssurf_ini[cid] = MM_tmp[index.get('iSsurf')]
        return MM_tmp, MM_S_tmp, nsl, perc_vol, etg_vol

    def step(self, ctx, n, tstart_MF, heads_cell, exf_cell, state, rejinf_cell=None):
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
        """
        nindex = len(ctx.index)
        nindex_S = len(ctx.index_S)
        MM_cells = np.zeros((ctx.ncell, nindex), dtype=np.float32)
        MM_S_cells = np.zeros((ctx.ncell, ctx._nslmax, nindex_S), dtype=np.float32)
        perc_cell = np.zeros(ctx.ncell, dtype=np.float32)
        etg_cell = np.zeros(ctx.ncell, dtype=np.float32)
        rej = (np.zeros(ctx.ncell) if rejinf_cell is None
               else np.asarray(rejinf_cell, dtype=np.float64))
        for cell in ctx.cells:
            cid = cell[0]
            MM_tmp, MM_S_tmp, nsl, perc_vol, etg_vol = self._cell_step(
                ctx, cell, n, tstart_MF, float(heads_cell[cid]), float(exf_cell[cid]), state,
                rejinf_cell=float(rej[cid]))
            MM_cells[cid, :] = MM_tmp
            MM_S_cells[cid, :nsl, :] = MM_S_tmp
            perc_cell[cid] = perc_vol
            etg_cell[cid] = etg_vol
        return {'MM': MM_cells, 'MM_S': MM_S_cells, 'perc': perc_cell, 'etg': etg_cell}

    def runMMsoil(self, _nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _Ssoil_ini, botm_l0, _Ks,
                  gridSOIL, gridSOILthick, TopSoil, gridMETEO,
                  index, index_S, gridSsurfhmax, gridSsurfw,
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
            gridSOIL, gridSOILthick, TopSoil, gridMETEO, index, index_S, gridSsurfhmax, gridSsurfw,
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
