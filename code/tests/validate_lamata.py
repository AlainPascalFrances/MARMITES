# -*- coding: utf-8 -*-
"""La Mata end-to-end validation of the refactored MMsoil (Phases 1-2).

Runs the refactored soil model against the REAL MODFLOW heads/exfiltration
(_h5_MF.h5) and compares its output to the REFERENCE _h5_MM.h5 produced by
the original Python-2 / Decimal code. Because the port replaced Decimal with
float64, results are expected to match within a small tolerance, not bit-for-bit.

Replicates the driver's setup between clsMF() and the runMMsoil() call, then
truncates to the first --nsp stress periods for runtime (state carry-over
makes the first K SPs directly comparable to the reference).

Usage:  python tests/validate_lamata.py --nsp 25
"""
import argparse
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy'):
    sys.path.insert(0, os.path.join(TRUNK, p))

import matplotlib  # noqa: E402
matplotlib.use('agg')
import h5py  # noqa: E402
import MARMITESutilities as MMutils  # noqa: E402
import ppMODFLOW_flopy_v3 as ppMF  # noqa: E402
import MARMITESsoil_v3 as MMsoil  # noqa: E402
from marmites_indices import INDEX_MM, INDEX_MM_SOIL  # noqa: E402


def read_list(cUTIL, ws, fn):
    return cUTIL.readFile(ws, fn)


def main(nsp):
    cUTIL = MMutils.clsUTILITIES(verbose=1)
    MF_ws = os.path.join(DS, 'MF_ws')
    cMF = ppMF.clsMF(cUTIL, MM_ws=DS, MM_ws_out=DS, MF_ws=MF_ws,
                     MF_ini_fn='__inputMF_flopy_v3_2s3L.ini',
                     xllcorner=739300.0, yllcorner=4553050.0)

    # conv_fact from length unit (lenuni==2 -> metres -> 1000 mm)
    conv_fact = {1: 304.8, 2: 1000.0, 3: 10.0}[cMF.lenuni]

    # ---- parse the MMsurf output file (names of daily series + veg params) ----
    # readFile returns tokens un-stripped (the driver strips at each call site)
    inp = [x.strip() for x in read_list(cUTIL, DS, '__inputMMsurf4MMsoil.txt')]
    l = 0
    NMETEO = int(inp[l]); l += 1
    NVEG = int(inp[l]); l += 1
    NSOIL = int(inp[l]); l += 1
    inputDate_fn = inp[l]; l += 1
    P_veg_fn = inp[l]; l += 1
    Pe_veg_fn = inp[l]; l += 1
    PT_fn = inp[l]; l += 1
    LAI_fn = inp[l]; l += 1
    PE_fn = inp[l]; l += 1
    Eo_fn = inp[l]; l += 1
    Zr = []; kTg_min = []; kTg_max = []; kT_f = []; kT_s = []
    _ = inp[l].split(); l += 1  # VegName
    Zr = [float(x) for x in inp[l].split()]; l += 1
    kTg_min = [float(x) for x in inp[l].split()]; l += 1
    kTg_max = [float(x) for x in inp[l].split()]; l += 1
    kT_f = [float(x) for x in inp[l].split()]; l += 1
    kT_s = [1.0 / float(x) for x in inp[l].split()]; l += 1
    # irrigation block
    NCROP = int(inp[l]); l += 1
    NFIELD = int(inp[l]); l += 1
    P_irr_fn = inp[l]; l += 1
    Pe_irr_fn = inp[l]; l += 1
    PT_irr_fn = inp[l]; l += 1
    Zr_c = np.array([float(x) for x in inp[l].split()]); l += 1
    kTg_min_c = np.array([float(x) for x in inp[l].split()]); l += 1
    kTg_max_c = np.array([float(x) for x in inp[l].split()]); l += 1
    kT_f_c = np.array([float(x) for x in inp[l].split()]); l += 1
    kT_s_c = np.array([1.0 / float(x) for x in inp[l].split()]); l += 1
    crop_irr_fn = inp[l]; l += 1

    # ---- MF time discretization (must reproduce the reference nper) ----
    cMF.ppMFtime(inputDate_fn, P_veg_fn, Pe_veg_fn, PT_fn, LAI_fn, PE_fn, Eo_fn,
                 NMETEO, NVEG, NSOIL, P_irr_fn, Pe_irr_fn, PT_irr_fn, crop_irr_fn, NFIELD)
    print('ppMFtime -> nper=%d, sum(perlen)=%d' % (cMF.nper, int(np.sum(cMF.perlen))))

    # ---- outcrop layer + masks (driver lines ~531-546) ----
    cMF.outcropL = np.zeros((cMF.nrow, cMF.ncol), dtype=int)
    mask_Lsup = np.zeros((cMF.nrow, cMF.ncol), dtype=int)
    for L in range(cMF.nlay):
        mask_Lsup += np.asarray(np.abs(cMF.ibound))[L, :, :]
        iboundBOL = (np.asarray(np.abs(cMF.ibound))[L, :, :] != 0)
        cMF.outcropL += (((cMF.outcropL == 0) & (iboundBOL == 1))) * (L + 1)

    # ---- read grids ----
    gridMETEO = cMF.cPROCESS.inputEsriAscii(grid_fn='inputMETEOzones.asc', datatype=int)
    gridSOIL = cMF.cPROCESS.inputEsriAscii(grid_fn='inputSOILzones.asc', datatype=int)
    gridSOILthick = cMF.cPROCESS.inputEsriAscii(grid_fn='inputSOILthick.asc', datatype=float)
    gridIRR = cMF.cPROCESS.inputEsriAscii(grid_fn='inputIRRzones.asc', datatype=int)

    # ---- input stress-period series ----
    (gridVEGarea, P_veg_zoneSP, Eo_zonesSP, PT_veg_zonesSP, Pe_veg_zonesSP, LAI_veg_zonesSP,
     PE_zonesSP, P_irr_zoneSP, Pe_irr_zoneSP, PT_irr_zonesSP, crop_irr_SP) = cMF.cPROCESS.inputSP(
        NMETEO=NMETEO, NVEG=NVEG, NSOIL=NSOIL, nper=cMF.nper,
        inputZON_SP_P_veg_fn=cMF.inputZON_SP_P_veg_fn, inputZON_SP_Pe_veg_fn=cMF.inputZON_SP_Pe_veg_fn,
        inputZON_SP_LAI_veg_fn=cMF.inputZON_SP_LAI_veg_fn, inputZON_SP_PT_fn=cMF.inputZON_SP_PT_fn,
        inputZON_SP_PE_fn=cMF.inputZON_SP_PE_fn, inputZON_SP_Eo_fn=cMF.inputZON_SP_Eo_fn,
        NFIELD=NFIELD, inputZON_SP_P_irr_fn=cMF.inputZON_SP_P_irr_fn,
        inputZON_SP_Pe_irr_fn=cMF.inputZON_SP_Pe_irr_fn, inputZON_SP_PT_irr_fn=cMF.inputZON_SP_PT_irr_fn,
        input_SP_crop_irr_fn=cMF.input_SP_crop_irr_fn)

    # ---- soil parameters ----
    SOILparam_fn = os.path.join('MF_ws', 'inputSOILparam.txt')
    _nsl, _nam, _st, _slprop, _Sm, _Sfc, _Sr, _S_ini, _Ks = cMF.cPROCESS.inputSoilParam(
        SOILparam_fn=SOILparam_fn, NSOIL=NSOIL)
    _nslmax = max(_nsl)
    for z in range(NSOIL):
        _slprop[z] = np.asarray(_slprop[z])

    # ---- soil-layer top/botm (driver lines ~634-646) ----
    cMF.elev = np.ma.masked_values(np.asarray(cMF.elev), cMF.hnoflo, atol=0.09)
    cMF.top = np.ma.masked_values(cMF.elev, cMF.hnoflo, atol=0.09) - np.ma.masked_values(gridSOILthick, cMF.hnoflo, atol=0.09)
    cMF.botm = np.asarray(cMF.botm)
    cMF.LandSurface = np.zeros(cMF.top.shape, dtype=np.float32)
    for L in range(cMF.nlay):
        cMF.botm[L, :, :] = np.ma.masked_values(cMF.botm[L, :, :], cMF.hnoflo, atol=0.09) - np.ma.masked_values(gridSOILthick, cMF.hnoflo, atol=0.09)
        cMF.LandSurface += cMF.top * (np.asarray(cMF.iuzfbnd) == L + 1)
    cMF.botm = np.ma.masked_values(cMF.botm, cMF.hnoflo, atol=0.09)
    botm_l0 = np.asarray(cMF.botm)[0, :, :]
    for L in range(cMF.nlay):
        cMF.iuzfbnd[cMF.ibound[L] <= 0] = 0

    # ---- truncate to first nsp stress periods for runtime ----
    nsp = min(nsp, cMF.nper)
    cMF.nper = nsp
    cMF.perlen = np.asarray(cMF.perlen)[:nsp]
    cMF.nstp = np.asarray(cMF.nstp)[:nsp]
    ndays = int(np.sum(cMF.perlen))
    print('Validating first %d stress periods (%d days), %d active cells' % (
        nsp, ndays, int((cMF.outcropL > 0).sum())))

    # ---- open real MF h5 (read), build fresh MM h5 (write) ----
    h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
    out_fn = os.path.join(HERE, '_h5_MM_validate.h5')
    if os.path.exists(out_fn):
        os.remove(out_fn)
    h5_MM = h5py.File(out_fn, 'w')
    nd = int(np.sum(cMF.perlen))
    h5_MM.create_dataset('MM', shape=(nd, cMF.nrow, cMF.ncol, len(INDEX_MM)), dtype=np.float32)
    h5_MM.create_dataset('MM_S', shape=(nd, cMF.nrow, cMF.ncol, _nslmax, len(INDEX_MM_SOIL)), dtype=np.float32)
    h5_MM.create_dataset('perc', shape=(cMF.nper, cMF.nrow, cMF.ncol), dtype=np.float32)
    h5_MM.create_dataset('ETg', shape=(cMF.nper, cMF.nrow, cMF.ncol), dtype=np.float32)

    MM_SOIL = MMsoil.clsMMsoil(hnoflo=cMF.hnoflo)
    MM_SOIL.runMMsoil(_nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _S_ini, botm_l0, _Ks,
                      gridSOIL, gridSOILthick, cMF.elev * 1000.0, gridMETEO,
                      INDEX_MM, INDEX_MM_SOIL,
                      P_veg_zoneSP, Eo_zonesSP, PT_veg_zonesSP, Pe_veg_zonesSP, PE_zonesSP, gridVEGarea,
                      LAI_veg_zonesSP, Zr, kTg_min, kTg_max, kT_f, kT_s, NVEG,
                      cMF, conv_fact, h5_MF, h5_MM, irr_yn=1,
                      P_irr_zoneSP=P_irr_zoneSP, PT_irr_zonesSP=PT_irr_zonesSP, Pe_irr_zoneSP=Pe_irr_zoneSP,
                      crop_irr_SP=crop_irr_SP, gridIRR=gridIRR,
                      Zr_c=Zr_c, kTg_min_c=kTg_min_c, kTg_max_c=kTg_max_c, kT_f_c=kT_f_c, kT_s_c=kT_s_c,
                      verbose=1)
    h5_MF.close()

    # ---- compare to reference over the first ndays days ----
    ref = h5py.File(os.path.join(DS, '_h5_MM.h5'), 'r')
    new = h5py.File(out_fn, 'r')
    print('\n=== comparison (first %d days) ===' % ndays)
    active = cMF.outcropL > 0
    results = {}
    for name in ('MM', 'MM_S', 'perc', 'ETg'):
        a = new[name][:]
        if name in ('MM', 'MM_S'):
            b = ref[name][:a.shape[0]]
            m = active[None, :, :, None] if name == 'MM' else active[None, :, :, None, None]
        else:
            b = ref[name][:a.shape[0]]
            m = active[None, :, :]
        a = np.where(np.broadcast_to(m, a.shape), a, 0.0)
        b = np.where(np.broadcast_to(m, b.shape), b, 0.0)
        finite = np.isfinite(a) & np.isfinite(b) & (np.abs(b) < 1e6)
        diff = np.abs(a[finite] - b[finite])
        denom = np.maximum(np.abs(b[finite]), 1e-6)
        rel = diff / denom
        results[name] = (float(diff.max()), float(np.median(diff)), float(rel.max()),
                         float(np.percentile(rel, 99)))
        print('%-5s  max|abs|=%.4g  median|abs|=%.4g  max_rel=%.3g  p99_rel=%.3g'
              % (name, *results[name]))
    ref.close(); new.close()

    # verdict: per-cell fluxes should agree within Decimal-rounding tolerance
    mm_ok = results['MM'][0] < 0.5 and results['MM'][3] < 0.05
    perc_ok = results['perc'][0] < 1e-3
    print('\nVERDICT:', 'PASS' if (mm_ok and perc_ok) else 'REVIEW',
          '(MM p99 rel %.2g, perc max abs %.2g)' % (results['MM'][3], results['perc'][0]))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--nsp', type=int, default=25, help='number of stress periods to validate')
    a = ap.parse_args()
    main(a.nsp)
