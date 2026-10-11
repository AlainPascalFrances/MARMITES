# -*- coding: utf-8 -*-
"""MMsoil's stress period over all cells, compiled with numba (WP7 speed).

THE SAME ARITHMETIC AS clsMMsoil._cell_step + flux, statement for statement,
for the coupled run (MF6, groundwater ET by EVT) -- and the CRR cascade
around it, so a day is ONE compiled pass over the cells. The profile of
2026-10-11 put the soil column at 3.1 s of a 4 s day on La Mata (4524 cells,
under the profiler), most of it Python's own overhead: small arrays, dicts
and attribute look-ups per cell per day. A PEST++ calibration runs the model
hundreds of times; that is what this is for.

clsMMsoil.step stays the REFERENCE: it is what runs when numba is missing,
when MM_SOIL_PYTHON is set, or on a path this kernel does not cover (the
uncoupled NWT run: hdry sentinels, MMsoil's own Eg/Tg and drawdown). The
tests hold the two to each other cell by cell.

The per-cell inputs that do not change during a run are PACKED once into
flat arrays (``pack``), keyed on the context, the cover and the cascade.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import os
from types import SimpleNamespace

import numpy as np

try:                                                    # pragma: no cover
    import numba as _numba
    _jit = _numba.njit(cache=True)
except ImportError:                                     # pragma: no cover
    _numba = None
    _jit = None

__all__ = ['available', 'pack', 'step', 'MM_NAMES']

# the MM outputs the kernel writes, in this order (positions resolved
# against ctx.index at pack time; a name the index lacks is skipped)
MM_NAMES = ('iP', 'iPT', 'iPE', 'iPe', 'iSsurf', 'iRo', 'iEXFg', 'iEow',
            'iMB', 'iEi', 'iEo', 'iEg', 'iTg', 'idSsurf', 'iETg', 'iETsoil',
            'iSsoil_pc', 'idSsoil', 'iperc', 'ihcorr', 'idgwt', 'iuzthick',
            'iI', 'iMBsurf', 'iPETuzf', 'iETuzf', 'iETtot', 'iRejInf',
            'iRunon', 'iReinf')
_IRUNON = MM_NAMES.index('iRunon')
_IREINF = MM_NAMES.index('iReinf')
# MM_S columns, POSITIONAL in the reference too:
# [Esoil, Tsoil, Ssoil_pc, Rp, Rexf, dSsoil, Ssoil, SAT, MB_l]
_S_SCALED = (0, 1, 3, 4, 5, 6, 8)   # iEsoil iTsoil iRsoil iExf idSsoil_s iSsoil iMB_s

# the context attributes pack() reads
_PACKED = ('cells', 'gridSOIL', 'gridMETEO', 'gridSOILthick', 'TopSoil',
           'botm_l0', '_slprop', '_Sm', '_Sfc', '_Sr', '_Ks', '_Ssoil_ini',
           '_st', '_nsl', 'gridVEGarea', 'gridIRR', 'P_veg_zoneSP',
           'Eo_zonesSP', 'PT_veg_zonesSP', 'Pe_veg_zonesSP', 'PE_zonesSP',
           'LAI_veg_zonesSP', 'P_irr_zoneSP', 'Pe_irr_zoneSP',
           'PT_irr_zonesSP', 'crop_irr_SP', 'Zr', 'kTg_min', 'kTg_max',
           'kT_f', 'kT_s', 'Zr_c', 'kTg_min_c', 'kTg_max_c', 'kT_f_c',
           'kT_s_c', 'cMF', 'index')

_TALLY_NAMES = ('Tg: soil moisture above porosity',
                'Tg: soil moisture below wilting point')


def available():
    """numba is there and not switched off for a check run."""
    return _jit is not None and not os.environ.get('MM_SOIL_PYTHON')


# ------------------------------------------------------------------ #
# the kernel
# ------------------------------------------------------------------ #

def _np_sum(a, n):
    """numpy's float sum of a[:n], in numpy's order (pairwise above 8)."""
    if n < 8:
        s = 0.0
        for i in range(n):
            s += a[i]
        return s
    r0, r1, r2, r3 = a[0], a[1], a[2], a[3]
    r4, r5, r6, r7 = a[4], a[5], a[6], a[7]
    m = n - n % 8
    for i in range(8, m, 8):
        r0 += a[i]
        r1 += a[i + 1]
        r2 += a[i + 2]
        r3 += a[i + 3]
        r4 += a[i + 4]
        r5 += a[i + 5]
        r6 += a[i + 6]
        r7 += a[i + 7]
    s = ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7))
    for i in range(m, n):
        s += a[i]
    return s


def _evp(s, Sm, Sr, pet, perlen, err):
    if not np.isfinite(s):
        err[0], err[1] = 2.0, s
        return 0.0
    if s <= Sr:
        return 0.0
    Se = min((s - Sr) / (Sm - Sr), 1.0)
    if pet * Se * perlen > (s - Sr):
        return (s - Sr) / perlen
    return pet * Se


def _perc(s, Sm, Sfc, Ks, perlen, err):
    if not np.isfinite(s):
        err[0], err[1] = 1.0, s
        return 0.0
    if s <= Sfc:
        return 0.0
    Sg = (s - Sfc) / (Sm - Sfc)
    if Ks * Sg * perlen > (s - Sfc):
        return (s - Sfc) / perlen
    return Ks * Sg


def _tally(tal, k, amount):
    tal[k, 0] += 1.0
    if tal[k, 0] == 1.0:
        tal[k, 1] = amount
        tal[k, 2] = amount
    else:
        tal[k, 1] = min(tal[k, 1], amount)
        tal[k, 2] = max(tal[k, 2], amount)


def _ktg(theta, Sr, Sm, kmin, kmax, kf, ks, tal):
    if theta > Sr:
        if theta < Sm:
            norm = (theta - Sr) / (Sm - Sr)
            return kmax - (kmax - kmin) / (1.0 + np.exp((norm - kf) / ks))
        if theta > Sm:
            _tally(tal, 0, theta - Sm)
        return kmax
    if theta < Sr:
        _tally(tal, 1, Sr - theta)
    return kmin


def _isclose(a, b):
    return abs(a - b) <= 1e-08 + 1e-05 * abs(b)


def _taken(surf, IRR, nsl, SmTl, Ssoil):
    room0 = SmTl[0] - Ssoil[0]
    took = room0 if surf > room0 else surf
    down = min(surf - took, IRR) if IRR > 0.0 else 0.0
    for l in range(1, nsl):
        if down <= 0.0:
            break
        t = min(down, max(SmTl[l] - Ssoil[l], 0.0))
        took += t
        down -= t
    return took


def _cell(c, n, t, first, perleni, perlen, perlen_n, h, exf, rej, etu,
          runon,
          # packed statics
          sz, mz, irf, nsl_c, Tl_c, Top_c, Bot_c, Sm_c, Sfc_c, Sr_c, Ks_c,
          Ssi_c, shah_c, VEG_c, hdry_mm_c, botm_c, botm32, f_open,
          sumTl_c,
          Zr, kmin, kmax, kf, kss, Zr_c, kmin_c, kmax_c, kf_c, ks_c,
          P_veg, PT_veg, Pe_veg, LAI, PE_z, Eo_z, P_irr, Pe_irr, PT_irr,
          crop, irr_column, hnoflo, conv, state, mmpos,
          # outputs
          MM, MMS, perc_o, etg_o, petuzf_o, eg_pe_o, land_o, h0_o, tgr_o,
          tgt_o, negexf, tal, err):
    fo = min(max(f_open[c], 0.0), 1.0)
    fs = 1.0 - fo
    col_rej = rej / fs if fs > 0.0 else 0.0
    col_etuzf = etu / fs if fs > 0.0 else 0.0
    col_runon = runon / fs if fs > 0.0 else 0.0
    z = sz[c]
    m = mz[c]
    nsl = nsl_c[c]
    Tl = Tl_c[c, :nsl]
    Top = Top_c[c, :nsl]
    Bot = Bot_c[c, :nsl]
    if first:
        Ssi = Ssi_c[c, :nsl].copy()
    else:
        Ssi = state[c, :nsl].copy()
    Sm = Sm_c[c, :nsl]
    Sfc = Sfc_c[c, :nsl]
    Sr = Sr_c[c, :nsl]
    Ks = Ks_c[c, :nsl]
    nveg = Zr.size
    f = irf[c]
    crop_id = -1
    if f >= 0:
        NV = 1
        LAIv = np.ones(1)
        crop_id = crop[f, t]
        P = P_irr[m, f, t]
        Pe_z = np.empty(1)
        Pe_z[0] = Pe_irr[m, f, t]
        PT_z = np.empty(1)
        PT_z[0] = PT_irr[m, f, t]
        nc = Zr_c.size
        k = crop_id - 1 if crop_id - 1 >= 0 else nc + (crop_id - 1)
        Zr_elev = np.empty(1)
        Zr_elev[0] = Top[0] - Zr_c[k] * 1000.0
        VEG = np.empty(1)
        VEG[0] = 100.0 if crop_id > 0 else 0.0
        kmn = np.empty(1)
        kmx = np.empty(1)
        kff = np.empty(1)
        kst = np.empty(1)
        kmn[0], kmx[0], kff[0], kst[0] = kmin_c[k], kmax_c[k], kf_c[k], \
            ks_c[k]
    else:
        NV = nveg
        P = P_veg[m, t]
        PT_z = np.empty(NV)
        Pe_z = np.empty(NV)
        LAIv = np.empty(NV)
        VEG = np.empty(NV)
        for v in range(NV):
            PT_z[v] = PT_veg[m, v, t]
            Pe_z[v] = Pe_veg[m, v, t]
            LAIv[v] = LAI[v, t]
            VEG[v] = VEG_c[c, v]
        Zr_elev = Top[0] - Zr * 1000.0
        kmn, kmx, kff, kst = kmin, kmax, kf, kss
    PEz = PE_z[m, z, t]
    Eo = Eo_z[m, t]
    SOILarea = 100.0
    Pe_tot = 0.0
    PT_tot = 0.0
    for v in range(NV):
        if LAIv[v] > 1.0E-5:
            Pe_tot += Pe_z[v] * VEG[v] * 0.01
            PT_tot += PT_z[v] * VEG[v] * 0.01
            SOILarea -= VEG[v]
    Pe_tot += P * SOILarea * 0.01
    INTER = P - Pe_tot
    irr_col = 0.0
    if f >= 0 and irr_column and P > 0.0:
        rain = P_veg[m, t]
        if P > rain:
            irr_col = Pe_tot * (P - rain) / P
    PE_tot = PEz * SOILarea * 0.01
    PT = np.empty(NV)
    for v in range(NV):
        PT[v] = PT_z[v] if LAIv[v] > 1.0E-5 else 0.0
    # dry cell (MF6: below the layer bottom). numpy compares a Python
    # float with a float32 grid value IN FLOAT32 (NEP 50), so the reference
    # does when botm_l0 is float32
    if botm32:
        dry = np.float32(h) < np.float32(botm_c[c])
    else:
        dry = h < botm_c[c]
    if dry:
        HEADS = hdry_mm_c[c]
    else:
        HEADS = h * 1000.0
    # the head of a dry cell is botm_l0 * 1000.0: a float32 when the grid
    # is, and numpy keeps it float32 through * 1.0, += Tl and * 0.001 below
    h32 = botm32 and dry
    if exf <= 0.0:
        dgwt = Top[0] - HEADS
    else:
        dgwt = sumTl_c[c]
        HEADS = Bot[nsl - 1]
        h32 = False
    uzthick = Bot[nsl - 1] - HEADS
    if first:
        Ssi = Ssi * Tl

    # ---------------- flux (GW_EVT) ----------------
    PE = PE_tot
    EXF_ini = exf
    if EXF_ini < 0.0:
        negexf[c] = EXF_ini
        EXF_ini = 0.0
    EXF = EXF_ini + abs(col_rej)
    SAT = np.zeros(nsl, dtype=np.bool_)
    Pe = Pe_tot
    Ssurf = Pe
    Ssoil = Ssi.copy()
    Rexf = np.zeros(nsl)
    SmTl = Sm * Tl
    SrTl = Sr * Tl
    SfcTl = Sfc * Tl
    Ssoil[nsl - 1] += EXF * perleni
    if EXF > 0.0:
        for l in range(nsl - 1, -1, -1):
            if Ssoil[l] >= SmTl[l]:
                Rexf[l] = Ssoil[l] - SmTl[l]
                Ssoil[l] = SmTl[l]
            if l != 0:
                Ssoil[l - 1] += Rexf[l]
        for l in range(nsl):
            Rexf[l] /= perlen
        Ssurf += Rexf[0] * perlen
    IRR = irr_col
    REinf = 0.0
    if col_runon > 0.0:
        own = Ssurf
        Ssurf += col_runon * perlen
        REinf = (_taken(Ssurf, IRR, nsl, SmTl, Ssoil)
                 - _taken(own, IRR, nsl, SmTl, Ssoil)) / perlen
    if Ssurf > (SmTl[0] - Ssoil[0]):
        I = SmTl[0] - Ssoil[0]
        Ssurf -= I
    else:
        I = Ssurf
        Ssurf = 0.0
    Ssoil[0] += I
    I /= perlen
    HEADS_corr = HEADS * 1.0
    dgwt_corr = dgwt
    if EXF > 0.0:
        for l in range(nsl - 1, -1, -1):
            if Ssoil[l] == SmTl[l]:
                if h32:
                    HEADS_corr = np.float64(np.float32(HEADS_corr)
                                            + np.float32(Tl[l]))
                else:
                    HEADS_corr += Tl[l]
                dgwt_corr -= Tl[l]
                SAT[l] = True
            else:
                break
    has_rcas = False
    Rcas = np.zeros(nsl)
    if IRR > 0.0 and Ssurf > 0.0:
        has_rcas = True
        down = min(Ssurf, IRR)
        for l in range(1, nsl):
            if down <= 0.0:
                break
            tt = min(down, max(SmTl[l] - Ssoil[l], 0.0))
            Ssoil[l] += tt
            for q in range(l):
                Rcas[q] += tt
            down -= tt
        Ssurf -= Rcas[0]
        I += Rcas[0] / perlen
        for l in range(nsl):
            Rcas[l] /= perlen
    Ro = Ssurf / perlen
    Ssurf = 0.0
    Eow = 0.0
    Rp = np.zeros(nsl)
    Tsoil = np.zeros(nsl)
    Esoil = np.zeros(nsl)
    for l in range(nsl):
        for v in range(NV):
            if PT[v] > 0.0 and LAIv[v] > 0.0 and VEG[v] > 0.0:
                TZr = 0.0
                if Bot[l] > Zr_elev[v]:
                    TZr = _evp(Ssoil[l], SmTl[l], SrTl[l], PT[v], perlen, err)
                elif Top[l] > Zr_elev[v]:
                    PTc = PT[v] * (Top[l] - Zr_elev[v]) / Tl[l]
                    TZr = _evp(Ssoil[l], SmTl[l], SrTl[l], PTc, perlen, err)
                Tsoil[l] += TZr * VEG[v] * 0.01
                PT[v] -= TZr
                if PT[v] < 0.0:
                    PT[v] = 0.0
        Ssoil[l] -= Tsoil[l] * perlen
        if Ssurf == 0.0 and PE > 0.0:
            Esoil[l] = _evp(Ssoil[l], SmTl[l], SrTl[l], PE, perlen, err)
            Ssoil[l] -= Esoil[l] * perlen
            PE -= Esoil[l]
            if PE < 0.0:
                PE = 0.0
        if l < nsl - 1:
            if not SAT[l + 1]:
                Rp[l] = _perc(Ssoil[l], SmTl[l], SfcTl[l], Ks[l], perlen, err)
                if Rp[l] * perlen > (SmTl[l + 1] - Ssoil[l + 1]):
                    Rp[l] = (SmTl[l + 1] - Ssoil[l + 1]) / perlen
                Ssoil[l + 1] += Rp[l] * perlen
        elif EXF == 0.0:
            Rp[l] = _perc(Ssoil[l], SmTl[l], SfcTl[l], Ks[l], perlen, err)
        Ssoil[l] -= Rp[l] * perlen
    if has_rcas:
        for l in range(nsl - 1):
            Rp[l] += Rcas[l]
    Spc = np.empty(nsl)
    for l in range(nsl):
        Spc[l] = Ssoil[l] / Tl[l]
    PT_left = 0.0
    for v in range(NV):
        if LAIv[v] > 0.0 and VEG[v] > 0.0:
            PT_left += PT[v] * VEG[v] * 0.01
    PETuzf = max(PE, 0.0) + PT_left
    ETuzf_used = min(max(col_etuzf, 0.0), PETuzf)
    if ETuzf_used > 0.0 and PETuzf > 0.0:
        ff = 1.0 - ETuzf_used / PETuzf
        PE = PE * ff
        for v in range(NV):
            PT[v] = PT[v] * ff
    eg_pe = PE if (Ssurf == 0.0 and PE > 0.0) else 0.0
    tg_rate = np.zeros(NV)
    order = np.argsort(Zr_elev, kind='mergesort')
    for oi in range(NV):
        v = order[oi]
        if not HEADS_corr > Zr_elev[v]:
            continue
        for l in range(nsl):
            if _isclose(Spc[l], hnoflo):
                continue
            kTg = _ktg(Spc[l], Sr[l], Sm[l], kmn[v], kmx[v], kff[v], kst[v],
                       tal)
            tq = PT[v] * kTg
            PT[v] -= tq
            tg_rate[v] += tq * VEG[v] * 0.01
    # ---------------- back in _cell_step ----------------
    eg_pe_o[c] = fs * eg_pe / 1000.0
    h0_o[c] = HEADS_corr / 1000.0
    land_o[c] = Top[0] / 1000.0
    for v in range(NV):
        tgr_o[c, v] = fs * tg_rate[v] / 1000.0
        tgt_o[c, v] = Zr_elev[v] / 1000.0
    Spc_tot = _np_sum(Spc, nsl) / nsl
    perc = Rp[nsl - 1]
    rej_in = abs(col_rej)
    exf_in = exf / perlen_n + rej_in
    Esoil_MB = _np_sum(Esoil, nsl)
    Tsoil_MB = _np_sum(Tsoil, nsl)
    dS = np.empty(nsl)
    for l in range(nsl):
        dS[l] = (Ssoil[l] - Ssi[l]) / perlen_n
    dS_tot = _np_sum(dS, nsl)
    ETsoil_tot = Esoil_MB + Tsoil_MB
    MB_l = np.zeros(nsl)
    if col_runon > 0.0:
        MBsurf = (Pe_tot + col_runon + Rexf[0]) - (Eow + Ro + I + 0.0)
    else:
        MBsurf = (Pe_tot + Rexf[0]) - (Eow + Ro + I + 0.0)
    if nsl > 1:
        MB_l[0] = (I + Rexf[1]) - (Rp[0] + Rexf[0] + Esoil[0] + Tsoil[0]
                                   + dS[0])
        for l in range(1, nsl - 1):
            MB_l[l] = (Rp[l - 1] + Rexf[l + 1]) - (
                Rp[l] + Rexf[l] + Esoil[l] + Tsoil[l] + dS[l])
        l = nsl - 1
        MB_l[l] = (Rp[l - 1] + exf_in) - (
            Rp[l] + Rexf[l] + Esoil[l] + Tsoil[l] + dS[l])
    else:
        MB_l[0] = (I + exf_in) - (
            Rp[0] + Rexf[0] + Esoil[0] + Tsoil[0] + dS[0])
    MB = (I + exf_in) - (Rexf[0] + Esoil_MB + Tsoil_MB + dS_tot
                         + Rp[nsl - 1])
    ro_open = fo * (P + max(exf, 0.0) / perlen_n)
    if runon > 0.0 and fs <= 0.0:
        ro_open += runon
    vals = np.zeros(30)
    vals[0] = P
    vals[1] = fs * PT_tot
    vals[2] = fs * PE_tot
    vals[3] = fs * Pe_tot + fo * P
    vals[4] = fs * Ssurf
    vals[5] = fs * Ro + ro_open
    vals[6] = exf
    vals[7] = fs * Eow
    vals[8] = fs * MB
    vals[9] = fs * INTER
    vals[10] = Eo
    vals[11] = fs * 0.0
    vals[12] = fs * 0.0
    vals[13] = fs * 0.0
    vals[14] = fs * 0.0
    vals[15] = fs * ETsoil_tot
    vals[16] = Spc_tot
    vals[17] = fs * dS_tot
    vals[18] = fs * perc
    if h32:
        vals[19] = np.float64(np.float32(HEADS_corr) * np.float32(0.001))
    else:
        vals[19] = HEADS_corr * 0.001
    vals[20] = -dgwt_corr * 0.001
    vals[21] = uzthick * 0.001
    vals[22] = fs * I
    vals[23] = fs * MBsurf
    vals[24] = fs * PETuzf
    vals[25] = 0.0
    vals[26] = fs * (INTER + Eow + ETsoil_tot + 0.0)
    vals[27] = fs * rej_in
    for q in range(28):
        if mmpos[q] >= 0:
            MM[c, mmpos[q]] = vals[q]
    if runon > 0.0:
        if mmpos[28] >= 0:
            MM[c, mmpos[28]] = fs * col_runon
        if mmpos[29] >= 0:
            MM[c, mmpos[29]] = fs * REinf
    for l in range(nsl):
        MMS[c, l, 0] = Esoil[l]
        MMS[c, l, 1] = Tsoil[l]
        MMS[c, l, 2] = Spc[l]
        MMS[c, l, 3] = Rp[l]
        MMS[c, l, 4] = Rexf[l]
        MMS[c, l, 5] = dS[l]
        MMS[c, l, 6] = Ssoil[l]
        MMS[c, l, 7] = 1.0 if SAT[l] else 0.0
        MMS[c, l, 8] = MB_l[l]
        state[c, l] = Ssoil[l]
    if fo > 0.0:
        for l in range(nsl):
            MMS[c, l, 0] *= fs
            MMS[c, l, 1] *= fs
            MMS[c, l, 3] *= fs
            MMS[c, l, 4] *= fs
            MMS[c, l, 5] *= fs
            MMS[c, l, 6] *= fs
            MMS[c, l, 8] *= fs
    perc_o[c] = MMS[c, nsl - 1, 3] / conv
    etg_o[c] = (MM[c, mmpos[14]] if mmpos[14] >= 0 else 0.0) / conv
    petuzf_o[c] = fs * PETuzf / conv
    return MM[c, mmpos[5]] if mmpos[5] >= 0 else vals[5]


def _day(n, t, first, perleni, perlen_n, heads, exf, rej, etu,
         sz, mz, irf, nsl_c, Tl_c, Top_c, Bot_c, Sm_c, Sfc_c, Sr_c, Ks_c,
         Ssi_c, shah_c, VEG_c, hdry_mm_c, botm_c, botm32, f_open,
         sumTl_c,
         Zr, kmin, kmax, kf, kss, Zr_c, kmin_c, kmax_c, kf_c, ks_c,
         P_veg, PT_veg, Pe_veg, LAI, PE_z, Eo_z, P_irr, Pe_irr, PT_irr,
         crop, irr_column, hnoflo, conv, state, mmpos,
         use_crr, order, kind, sink, route_to, beta, area,
         s_ptr, s_idx, s_a, w_ptr, w_idx, w_a,
         MM, MMS, perc_o, etg_o, petuzf_o, eg_pe_o, land_o, h0_o, tgr_o,
         tgt_o, negexf, tal, err,
         ro, runon_v, deliver, ecrr, ecrr_sink, routed):
    """One stress period: every cell in the cascade's order (or in list
    order without it), the run-on handed on as the reference route does."""
    ncell = heads.size
    perlen = perlen_n
    nmax = 1
    if use_crr:
        for k in range(ncell):
            nmax = max(nmax, s_ptr[k + 1] - s_ptr[k], w_ptr[k + 1] - w_ptr[k])
    vs = np.empty(nmax)
    vw = np.empty(nmax)
    for oi in range(ncell):
        k = order[oi] if use_crr else oi
        ron = runon_v[k] / area[k] * conv if use_crr else 0.0
        r = _cell(k, n, t, first, perleni, perlen, perlen_n, heads[k],
                  exf[k], rej[k], etu[k], ron,
                  sz, mz, irf, nsl_c, Tl_c, Top_c, Bot_c, Sm_c, Sfc_c, Sr_c,
                  Ks_c, Ssi_c, shah_c, VEG_c, hdry_mm_c, botm_c, botm32,
                  f_open, sumTl_c, Zr, kmin, kmax, kf, kss, Zr_c, kmin_c, kmax_c,
                  kf_c, ks_c, P_veg, PT_veg, Pe_veg, LAI, PE_z, Eo_z, P_irr,
                  Pe_irr, PT_irr, crop, irr_column, hnoflo, conv, state,
                  mmpos, MM, MMS, perc_o, etg_o, petuzf_o, eg_pe_o, land_o,
                  h0_o, tgr_o, tgt_o, negexf, tal, err)
        if err[0] != 0.0:
            err[2] = k
            return
        if not use_crr:
            continue
        q = max(r, 0.0) * area[k] / conv
        ro[k] = q
        if q == 0.0:
            continue
        if kind[k] != 0:
            deliver[k] += q
            continue
        if not sink[k]:
            ns = s_ptr[k + 1] - s_ptr[k]
            nw = w_ptr[k + 1] - w_ptr[k]
            for i in range(ns):
                vs[i] = q * s_a[s_ptr[k] + i]
                runon_v[s_idx[s_ptr[k] + i]] += vs[i]
            for i in range(nw):
                vw[i] = q * w_a[w_ptr[k] + i]
                deliver[w_idx[w_ptr[k] + i]] += vw[i]
            ecrr[k] += q - (_np_sum(vs, ns) + _np_sum(vw, nw))
        elif route_to[k] >= 0:
            v = q * beta
            deliver[route_to[k]] += v
            routed[k] += v
            ecrr[k] += q - v
        else:
            ecrr[k] += q
            ecrr_sink[k] += q


if _jit is not None:
    _np_sum = _jit(_np_sum)
    _evp = _jit(_evp)
    _perc = _jit(_perc)
    _tally = _jit(_tally)
    _ktg = _jit(_ktg)
    _isclose = _jit(_isclose)
    _taken = _jit(_taken)
    _cell = _jit(_cell)
    _day = _jit(_day)


# ------------------------------------------------------------------ #
# packing the static inputs, once a run
# ------------------------------------------------------------------ #

def _f64(a, shape=None):
    a = np.ascontiguousarray(np.asarray(a, dtype=np.float64))
    return a.reshape(shape) if shape is not None else a


def pack(mm, ctx, crr=None):
    """The context as flat arrays for the kernel, or None when this path is
    not the kernel's (then the reference runs). Cached on the context."""
    cMF = ctx.cMF
    if getattr(cMF, 'hdry', None) is not None or not getattr(ctx, 'gw_evt',
                                                              False):
        return None
    if list(ctx.index_S) != ['iEsoil', 'iTsoil', 'iSsoil_pc_s', 'iRsoil',
                             'iExf', 'idSsoil_s', 'iSsoil', 'iSAT', 'iMB_s']:
        return None
    fo = getattr(ctx, 'f_open', None)
    # every object the pack reads, by identity: a context whose forcing or
    # grids were REPLACED is packed again (one mutated in place is not --
    # nothing in the model does that)
    key = (id(fo), id(crr), bool(getattr(ctx, 'irr_column', False))) + \
        tuple(id(getattr(ctx, a, None)) for a in _PACKED)
    K = getattr(ctx, '_nb_pack', None)
    if K is not None and K.key == key:
        return K
    cells = ctx.cells
    n = len(cells)
    I = np.array([c[1] for c in cells], dtype=int)
    J = np.array([c[2] for c in cells], dtype=int)
    if any(int(c[0]) != k for k, c in enumerate(cells)):
        return None                       # the cell id must be its position
    gS = np.asarray(ctx.gridSOIL)
    gM = np.asarray(ctx.gridMETEO)
    sz = (gS[I, J].astype(np.int64) - 1)
    mz = (gM[I, J].astype(np.int64) - 1)
    nslmax = int(ctx._nslmax)
    nsl_c = np.array([int(ctx._nsl[z]) for z in sz], dtype=np.int64)
    Tl = np.zeros((n, nslmax))
    Top = np.zeros((n, nslmax))
    Bot = np.zeros((n, nslmax))
    Sm = np.zeros((n, nslmax))
    Sfc = np.zeros((n, nslmax))
    Sr = np.zeros((n, nslmax))
    Ks = np.zeros((n, nslmax))
    Ssi = np.zeros((n, nslmax))
    sumTl = np.zeros(n)
    hdry_mm = np.zeros(n)
    botm = np.zeros(n)
    shah = np.zeros((n, 4))
    for c in range(n):
        i, j, z, nsl = I[c], J[c], sz[c], nsl_c[c]
        # the reference's own expressions, on its own objects, so the
        # dtypes (float32 grids) round exactly as they do there
        tl = np.asarray(ctx.gridSOILthick[i, j] * ctx._slprop[z] * 1000.0,
                        dtype=np.float64)
        Tl[c, :nsl] = tl[:nsl]
        top = np.zeros(nsl, dtype=np.float64)
        bot = np.zeros(nsl, dtype=np.float64)
        for l in range(nsl):
            if l == 0:
                top[0] = ctx.TopSoil[i, j]
                bot[0] = ctx.TopSoil[i, j] - tl[0]
            else:
                top[l] = bot[l - 1]
                bot[l] = top[l] - tl[l]
        Top[c, :nsl], Bot[c, :nsl] = top, bot
        Sm[c, :nsl] = np.asarray(ctx._Sm[z], dtype=np.float64)[:nsl]
        Sfc[c, :nsl] = np.asarray(ctx._Sfc[z], dtype=np.float64)[:nsl]
        Sr[c, :nsl] = np.asarray(ctx._Sr[z], dtype=np.float64)[:nsl]
        Ks[c, :nsl] = np.asarray(ctx._Ks[z], dtype=np.float64)[:nsl]
        Ssi[c, :nsl] = np.asarray(ctx._Ssoil_ini[z][:nsl], dtype=np.float64)
        sumTl[c] = float(np.sum(tl[:nsl]))
        hdry_mm[c] = ctx.botm_l0[i, j] * 1000.0
        botm[c] = ctx.botm_l0[i, j]
        p = mm.paramEg[ctx._st[z]]
        shah[c] = (p['y0'], p['b'], p['dll'], p['ext_d'])
    nveg = int(ctx.NVEG)
    gV = ctx.gridVEGarea
    VEG = np.zeros((n, nveg))
    for v in range(nveg):
        VEG[:, v] = np.asarray(gV[v], dtype=np.float64)[I, J]
    irr_yn = int(ctx.irr_yn) if ctx.irr_yn is not None else 0
    irf = np.full(n, -1, dtype=np.int64)
    if irr_yn == 1:
        g = np.asarray(ctx.gridIRR)[I, J].astype(np.int64)
        irf = np.where(g > 0, g - 1, -1).astype(np.int64)
    nm = int(np.max(mz)) + 1 if n else 1
    P_veg = _f64([np.asarray(ctx.P_veg_zoneSP[k], dtype=np.float64)
                  for k in range(nm)])
    Eo = _f64([np.asarray(ctx.Eo_zonesSP[k], dtype=np.float64)
               for k in range(nm)])
    if irr_yn == 1 and (irf >= 0).any():
        P_irr = _f64(ctx.P_irr_zoneSP)
        Pe_irr = _f64(ctx.Pe_irr_zoneSP)
        PT_irr = _f64(ctx.PT_irr_zonesSP)
        crop = np.ascontiguousarray(np.asarray(ctx.crop_irr_SP)
                                    .astype(np.int64))
        Zr_c = _f64(np.ravel(ctx.Zr_c))
        kmin_c = _f64(np.ravel(ctx.kTg_min_c))
        kmax_c = _f64(np.ravel(ctx.kTg_max_c))
        kf_c = _f64(np.ravel(ctx.kT_f_c))
        ks_c = _f64(np.ravel(ctx.kT_s_c))
    else:
        P_irr = Pe_irr = PT_irr = np.zeros((1, 1, 1))
        crop = np.zeros((1, 1), dtype=np.int64)
        Zr_c = kmin_c = kmax_c = kf_c = ks_c = np.zeros(1)
    f_open = (np.zeros(n) if fo is None else _f64(np.ravel(fo)[:n]))
    pos = np.array([ctx.index.get(k, -1) for k in MM_NAMES], dtype=np.int64)
    K = SimpleNamespace(
        key=key, n=n, nveg=nveg, nslmax=nslmax, I=I, J=J,
        sz=sz, mz=mz, irf=irf, nsl=nsl_c, Tl=Tl, Top=Top, Bot=Bot, Sm=Sm,
        Sfc=Sfc, Sr=Sr, Ks=Ks, Ssi=Ssi, shah=shah, VEG=VEG,
        hdry_mm=hdry_mm, botm=botm, f_open=f_open, sumTl=sumTl,
        botm32=bool(np.asarray(ctx.botm_l0).dtype == np.float32),
        Zr=_f64([float(x) for x in np.ravel(ctx.Zr)[:nveg]]),
        kmin=_f64([float(x) for x in np.ravel(ctx.kTg_min)[:nveg]]),
        kmax=_f64([float(x) for x in np.ravel(ctx.kTg_max)[:nveg]]),
        kf=_f64([float(x) for x in np.ravel(ctx.kT_f)[:nveg]]),
        kss=_f64([float(x) for x in np.ravel(ctx.kT_s)[:nveg]]),
        Zr_c=Zr_c, kmin_c=kmin_c, kmax_c=kmax_c, kf_c=kf_c, ks_c=ks_c,
        P_veg=P_veg, PT_veg=_f64(ctx.PT_veg_zonesSP),
        Pe_veg=_f64(ctx.Pe_veg_zonesSP), LAI=_f64(ctx.LAI_veg_zonesSP),
        PE_z=_f64(ctx.PE_zonesSP), Eo_z=Eo, P_irr=P_irr, Pe_irr=Pe_irr,
        PT_irr=PT_irr, crop=crop,
        irr_column=bool(getattr(ctx, 'irr_column', False)),
        hnoflo=float(cMF.hnoflo), conv=float(ctx.conv_fact),
        perlen=_f64(np.ravel(cMF.perlen)), mmpos=pos,
        paramEg=[mm.paramEg[ctx._st[z]] for z in sz])
    # the cascade, as compressed rows
    K.use_crr = crr is not None
    if crr is not None:
        sp, si, sa, wp, wi, wa = [0], [], [], [0], [], []
        for k in range(n):
            rs, a_s, rw, a_w = crr._split[k]
            si.extend(int(x) for x in rs)
            sa.extend(float(x) for x in a_s)
            wi.extend(int(x) for x in rw)
            wa.extend(float(x) for x in a_w)
            sp.append(len(si))
            wp.append(len(wi))
        K.order = np.ascontiguousarray(crr.order, dtype=np.int64)
        K.kind = np.ascontiguousarray(crr.kind, dtype=np.int64)
        K.sink = np.ascontiguousarray(crr.sink, dtype=np.bool_)
        K.route_to = np.ascontiguousarray(crr.route_to, dtype=np.int64)
        K.beta = float(crr.beta)
        K.s_ptr, K.w_ptr = np.array(sp, np.int64), np.array(wp, np.int64)
        K.s_idx = np.array(si, np.int64) if si else np.zeros(1, np.int64)
        K.s_a = np.array(sa, np.float64) if sa else np.zeros(1)
        K.w_idx = np.array(wi, np.int64) if wi else np.zeros(1, np.int64)
        K.w_a = np.array(wa, np.float64) if wa else np.zeros(1)
    else:
        K.order = np.arange(n, dtype=np.int64)
        K.kind = np.zeros(n, np.int64)
        K.sink = np.zeros(n, np.bool_)
        K.route_to = np.full(n, -1, np.int64)
        K.beta = 1.0
        K.s_ptr = K.w_ptr = np.zeros(n + 1, np.int64)
        K.s_idx = K.w_idx = np.zeros(1, np.int64)
        K.s_a = K.w_a = np.zeros(1)
    ctx._nb_pack = K
    return K


# ------------------------------------------------------------------ #
# the stress period
# ------------------------------------------------------------------ #

def step(mm, ctx, n, tstart_MF, heads_cell, exf_cell, state, rej, etu,
         crr=None):
    """clsMMsoil.step, compiled. Returns its dict, or None when the
    reference has to run (see :func:`pack`)."""
    if not available():
        return None
    K = pack(mm, ctx, crr)
    if K is None:
        return None
    import MARMITESsoil_v3 as ref
    nc, nveg = K.n, K.nveg
    first = (n == 0 and not getattr(state, 'carried', False))
    if first:
        perleni = 1.0
    elif n == 0:
        perleni = float(ctx.cMF.perlen[-1])
    else:
        perleni = float(ctx.cMF.perlen[n - 1])
    perlen_n = float(ctx.cMF.perlen[n])
    nindex, nindex_S = len(ctx.index), len(ctx.index_S)
    MM = np.zeros((nc, nindex))
    MMS = np.zeros((nc, K.nslmax, nindex_S))
    perc = np.zeros(nc)
    etg = np.zeros(nc)
    petuzf = np.zeros(nc)
    eg_pe, land, h0 = np.zeros(nc), np.zeros(nc), np.zeros(nc)
    nv = 1 if (K.irf >= 0).all() and nc else nveg
    nv = max(nv, 1)
    tgr = np.zeros((nc, max(nveg, 1)))
    tgt = np.full((nc, max(nveg, 1)), np.nan)
    negexf = np.zeros(nc)
    tal = np.zeros((2, 3))
    err = np.zeros(3)
    ro, runon_v, deliver = np.zeros(nc), np.zeros(nc), np.zeros(nc)
    ecrr, ecrr_sink, routed = np.zeros(nc), np.zeros(nc), np.zeros(nc)
    area = (np.asarray(ctx.geom.area, dtype=np.float64) if K.use_crr
            else np.ones(nc))
    st = np.ascontiguousarray(state.Ssoil_ini, dtype=np.float64)
    _day(int(n), int(tstart_MF), bool(first), perleni, perlen_n,
         _f64(heads_cell), _f64(exf_cell), _f64(rej), _f64(etu),
         K.sz, K.mz, K.irf, K.nsl, K.Tl, K.Top, K.Bot, K.Sm, K.Sfc, K.Sr,
         K.Ks, K.Ssi, K.shah, K.VEG, K.hdry_mm, K.botm, K.botm32,
         K.f_open, K.sumTl,
         K.Zr, K.kmin, K.kmax, K.kf, K.kss, K.Zr_c, K.kmin_c, K.kmax_c,
         K.kf_c, K.ks_c, K.P_veg, K.PT_veg, K.Pe_veg, K.LAI, K.PE_z, K.Eo_z,
         K.P_irr, K.Pe_irr, K.PT_irr, K.crop, K.irr_column, K.hnoflo,
         K.conv, st, K.mmpos, K.use_crr, K.order, K.kind, K.sink,
         K.route_to, K.beta, area, K.s_ptr, K.s_idx, K.s_a, K.w_ptr,
         K.w_idx, K.w_a, MM, MMS, perc, etg, petuzf, eg_pe, land, h0, tgr,
         tgt, negexf, tal, err, ro, runon_v, deliver, ecrr, ecrr_sink,
         routed)
    if err[0] != 0.0:
        what = 'percolation' if err[0] == 1.0 else 'evapotranspiration'
        raise ref.MarmitesSoilError('Non-finite soil storage in %s: %r'
                                    % (what, float(err[1])))
    state.Ssoil_ini[...] = st
    # what the reference prints and tallies, in its order
    for k in (K.order if K.use_crr else range(nc)):
        if negexf[k] < 0.0:
            print('WARNING!\nEXFg < 0.0, value %.6f corrected to 0.0.'
                  % negexf[k])
    for q, what in enumerate(_TALLY_NAMES):
        if tal[q, 0] > 0:
            rec = ref._TALLIES.setdefault(what, {'n': 0, 'lo': None,
                                                 'hi': None})
            rec['n'] += int(tal[q, 0])
            lo, hi = float(tal[q, 1]), float(tal[q, 2])
            rec['lo'] = lo if rec['lo'] is None else min(rec['lo'], lo)
            rec['hi'] = hi if rec['hi'] is None else max(rec['hi'], hi)
    out = {'MM': MM, 'MM_S': MMS, 'perc': perc, 'etg': etg,
           'petuzf': petuzf}
    # the coupled run's groundwater-ET potential, as _gw_arrays gives it
    tgr_w, tgt_w = tgr[:, :nv], tgt[:, :nv]
    out['gw'] = {'eg_pe': eg_pe, 'land': land, 'h0': h0,
                 'shah': K.shah, 'tg_rate': tgr_w, 'tg_tip': tgt_w}
    if K.use_crr:
        res = SimpleNamespace(ro=ro, runon=runon_v, deliver=deliver,
                              ecrr=ecrr, ecrr_sink=ecrr_sink, routed=routed)
        crr.check(res)
        conv = float(ctx.conv_fact)
        ec = res.ecrr / area * conv
        if 'iEcrr' in ctx.index:
            MM[:, ctx.index['iEcrr']] = ec
        if 'iETtot' in ctx.index:
            MM[:, ctx.index['iETtot']] += ec
        out['ro_deliver'] = res.deliver / area * conv
        out['crr'] = res
    return out
