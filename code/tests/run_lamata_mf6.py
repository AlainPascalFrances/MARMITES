# -*- coding: utf-8 -*-
"""La Mata MODFLOW 6 coupled run (Phase 3 entry point).

Builds the MF6 simulation from the La Mata dataset and, when the MODFLOW 6
library is available, drives the coupled MARMITES-MF6 model through the API
(lagged or iterative mode). Without libmf6 it stops after writing the
simulation files (still useful: run `mf6` manually in the workspace to
check the groundwater model alone).

Requirements on the executing machine:
    pip install flopy modflowapi
    MODFLOW 6 >= 6.4 binaries:  mf6[.exe] and libmf6[.dll|.so]
    (easiest: `get-modflow :flopy` or download from
     github.com/MODFLOW-USGS/executables)

Usage:
    python tests/run_lamata_mf6.py --build-only
    python tests/run_lamata_mf6.py --libmf6 C:/path/to/libmf6.dll --mode lagged --nsp 60
    python tests/run_lamata_mf6.py --libmf6 ... --mode iterative --relax 0.6
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))

# The repository holds ONLY code, input data and docs. Everything a run
# produces goes to a workspace outside it, laid out as
#     <WS_ROOT>/MF6_ws/                 the MODFLOW 6 model + its output
#     <WS_ROOT>/MMsurf_ws/              MMsurf output
#     <WS_ROOT>/out_<stamp>_<tag>/      MM results (postproc/ + figures/)
# Override with --ws-root or the MARMITES_WS_ROOT environment variable.
WS_ROOT = os.environ.get('MARMITES_WS_ROOT',
                         os.path.join('E:' + os.sep, '00code_ws', 'LaMata_MM-MF6'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

import matplotlib  # noqa: E402
matplotlib.use('agg')
import h5py  # noqa: E402
import MARMITESutilities as MMutils  # noqa: E402
import ppMODFLOW_flopy_v3 as ppMF  # noqa: E402
import MARMITESsoil_v3 as MMsoil  # noqa: E402
from marmites_indices import INDEX_MM, INDEX_MM_SOIL  # noqa: E402
from marmites_mf6 import clsMF6  # noqa: E402
from marmites_coupler import MF6Coupler  # noqa: E402


def _asc(fn):
    """Read an ESRI ASCII grid, nodata -> 0."""
    a = np.loadtxt(fn, skiprows=6)
    return np.where(a <= -9999.0, 0.0, a)


def _write_cell_grid(vals, cells, nrow, ncol, cMF, fn, nodata=-9999.0):
    """Scatter a per-MM-cell field onto the grid and write it as ESRI-ASCII."""
    g = np.full((nrow, ncol), nodata)
    for v, c in zip(vals, cells):
        g[c[1], c[2]] = v
    hdr = ('ncols %d\nnrows %d\nxllcorner %s\nyllcorner %s\ncellsize %s\n'
           'nodata_value %g\n' % (ncol, nrow, cMF.xllcorner, cMF.yllcorner,
                                  float(np.mean(cMF.delr)), nodata))
    with open(fn, 'w') as f:
        f.write(hdr)
        np.savetxt(f, g, fmt='%.6g')


def _read_cell_grid(fn, cells):
    """Read an ESRI-ASCII grid and gather it back to a per-MM-cell vector."""
    g = np.loadtxt(fn, skiprows=6)
    return np.array([g[c[1], c[2]] for c in cells], dtype=float)


def _state_in(a, pref, probe):
    """Resolve a saved-state prefix for READING.

    Run state (equilibrated heads, steady means) is written to the workspace,
    but a baseline set is committed in the repo's DataSet_LaMata/MF_ws. Prefer
    the workspace copy, fall back to the repo one, so `--strt-heads hi_spinup`
    works on a fresh clone and picks up a newer spin-up once one exists.
    ``probe`` is the suffix that identifies the set (e.g. '_l1.asc').
    """
    if os.path.isabs(pref):
        return pref
    cand = os.path.join(a.state_dir, pref)
    if os.path.exists(cand + probe):
        return cand
    return os.path.join(DS, 'MF_ws', pref)


def _state_out(a, pref):
    """Resolve a saved-state prefix for WRITING -- always the workspace, never
    the repository (which holds input data only)."""
    if os.path.isabs(pref):
        return pref
    os.makedirs(a.state_dir, exist_ok=True)
    return os.path.join(a.state_dir, pref)


# La Mata is parameterised at two vertical resolutions. Both are authoritative
# MODFLOW input sets maintained by the modeller -- the 2-layer file is NOT a
# derived aggregation of the 6-layer one; its Sy (0.01 uniform) and its
# confining layer-2 Ss (~1e-7) are hand-set and cannot be reproduced by any
# thickness-weighted mean of the 6-layer values. So for a 2-layer run we read
# the 2-layer file directly rather than aggregating the 6-layer one.
INI_BY_NLAY = {2: '__inputMF_flopy_v3_2s1L.ini', 6: '__inputMF_flopy_v3_2s3L.ini'}


def setup_lamata(daily=True, nsp=None, grid='dis', nlay=None, aggregate=False):
    """Replicate the driver setup; returns (cMF, mm, ctx, state, top, botm).

    ``nlay`` selects the parameter set: 2 reads the 2-layer ini directly, 6 (or
    None) the 6-layer ini. ``aggregate=True`` forces the old behaviour of
    deriving the 2-layer model by collapsing the 6-layer one (kept for
    comparison only).
    """
    if aggregate:
        ini_fn = INI_BY_NLAY[6]
    else:
        ini_fn = INI_BY_NLAY.get(int(nlay) if nlay else 6)
        if ini_fn is None:
            raise SystemExit('--nlay %s has no parameter file; available: %s '
                             '(or use --aggregate to derive it from the 6-layer '
                             'model)' % (nlay, sorted(INI_BY_NLAY)))
    cUTIL = MMutils.clsUTILITIES(verbose=1)
    cMF = ppMF.clsMF(cUTIL, MM_ws=DS, MM_ws_out=DS, MF_ws=os.path.join(DS, 'MF_ws'),
                     MF_ini_fn=ini_fn,
                     xllcorner=739300.0, yllcorner=4553050.0)
    print('parameter set: %s (%d layer(s))' % (ini_fn, cMF.nlay))
    conv_fact = {1: 304.8, 2: 1000.0, 3: 10.0}[cMF.lenuni]

    inp = [x.strip() for x in cUTIL.readFile(DS, '__inputMMsurf4MMsoil.txt')]
    l = 0
    NMETEO = int(inp[l]); l += 1
    NVEG = int(inp[l]); l += 1
    NSOIL = int(inp[l]); l += 1
    inputDate_fn = inp[l]; l += 1
    P_veg_fn, Pe_veg_fn, PT_fn, LAI_fn, PE_fn, Eo_fn = inp[l:l + 6]; l += 6
    _ = inp[l]; l += 1  # VegName
    Zr = [float(x) for x in inp[l].split()]; l += 1
    kTg_min = [float(x) for x in inp[l].split()]; l += 1
    kTg_max = [float(x) for x in inp[l].split()]; l += 1
    kT_f = [float(x) for x in inp[l].split()]; l += 1
    kT_s = [1.0 / float(x) for x in inp[l].split()]; l += 1
    NCROP = int(inp[l]); l += 1
    NFIELD = int(inp[l]); l += 1
    P_irr_fn, Pe_irr_fn, PT_irr_fn = inp[l:l + 3]; l += 3
    Zr_c = np.array([float(x) for x in inp[l].split()]); l += 1
    kTg_min_c = np.array([float(x) for x in inp[l].split()]); l += 1
    kTg_max_c = np.array([float(x) for x in inp[l].split()]); l += 1
    kT_f_c = np.array([float(x) for x in inp[l].split()]); l += 1
    kT_s_c = np.array([1.0 / float(x) for x in inp[l].split()]); l += 1
    crop_irr_fn = inp[l]

    if daily:
        cMF.nper = 1     # perlenmax=1 -> ppMFtime produces daily SPs (decision 4.3)
    cMF.ppMFtime(inputDate_fn, P_veg_fn, Pe_veg_fn, PT_fn, LAI_fn, PE_fn, Eo_fn,
                 NMETEO, NVEG, NSOIL, P_irr_fn, Pe_irr_fn, PT_irr_fn, crop_irr_fn, NFIELD)
    print('time discretization: nper=%d over %d days (daily=%s)'
          % (cMF.nper, int(np.sum(cMF.perlen)), daily))

    # Comparison-only path: derive the 2-layer model by collapsing the 6-layer
    # one. The direct 2-layer parameter file (loaded above when nlay=2) is the
    # authoritative representation; aggregation is kept so the two can be
    # compared. Done BEFORE the outcrop map and the MARMITES cell list, so soil
    # and groundwater models always share the same layering.
    if aggregate and nlay is not None and int(nlay) != cMF.nlay:
        from marmites_layers import aggregate_layers
        if int(nlay) != int(cMF.Mnlay):
            raise SystemExit('--nlay %s requested but the ini declares Mnlay=%d '
                             '(Mlay=%s)' % (nlay, cMF.Mnlay, cMF.Mlay))
        aggregate_layers(cMF)

    # outcrop / masks
    cMF.outcropL = np.zeros((cMF.nrow, cMF.ncol), dtype=int)
    for L in range(cMF.nlay):
        ib = (np.abs(np.asarray(cMF.ibound))[L] != 0)
        cMF.outcropL += ((cMF.outcropL == 0) & ib) * (L + 1)

    gridMETEO = cMF.cPROCESS.inputEsriAscii(grid_fn='inputMETEOzones.asc', datatype=int)
    gridSOIL = cMF.cPROCESS.inputEsriAscii(grid_fn='inputSOILzones.asc', datatype=int)
    gridSOILthick = cMF.cPROCESS.inputEsriAscii(grid_fn='inputSOILthick.asc', datatype=float)
    gridSsurfhmax = cMF.cPROCESS.inputEsriAscii(grid_fn='inputPONDhmax.asc', datatype=float)
    gridSsurfw = cMF.cPROCESS.inputEsriAscii(grid_fn='inputPONDw.asc', datatype=float)
    gridIRR = cMF.cPROCESS.inputEsriAscii(grid_fn='inputIRRzones.asc', datatype=int)

    (gridVEGarea, P_veg_zoneSP, Eo_zonesSP, PT_veg_zonesSP, Pe_veg_zonesSP, LAI_veg_zonesSP,
     PE_zonesSP, P_irr_zoneSP, Pe_irr_zoneSP, PT_irr_zonesSP, crop_irr_SP) = cMF.cPROCESS.inputSP(
        NMETEO=NMETEO, NVEG=NVEG, NSOIL=NSOIL, nper=cMF.nper,
        inputZON_SP_P_veg_fn=cMF.inputZON_SP_P_veg_fn, inputZON_SP_Pe_veg_fn=cMF.inputZON_SP_Pe_veg_fn,
        inputZON_SP_LAI_veg_fn=cMF.inputZON_SP_LAI_veg_fn, inputZON_SP_PT_fn=cMF.inputZON_SP_PT_fn,
        inputZON_SP_PE_fn=cMF.inputZON_SP_PE_fn, inputZON_SP_Eo_fn=cMF.inputZON_SP_Eo_fn,
        NFIELD=NFIELD, inputZON_SP_P_irr_fn=cMF.inputZON_SP_P_irr_fn,
        inputZON_SP_Pe_irr_fn=cMF.inputZON_SP_Pe_irr_fn, inputZON_SP_PT_irr_fn=cMF.inputZON_SP_PT_irr_fn,
        input_SP_crop_irr_fn=cMF.input_SP_crop_irr_fn)

    _nsl, _nam, _st, _slprop, _Sm, _Sfc, _Sr, _S_ini, _Ks = cMF.cPROCESS.inputSoilParam(
        SOILparam_fn=os.path.join('MF_ws', 'inputSOILparam.txt'), NSOIL=NSOIL)
    _nslmax = max(_nsl)
    for z in range(NSOIL):
        _slprop[z] = np.asarray(_slprop[z])

    # driver top/botm adjustment (aquifer sits below the soil column)
    cMF.elev = np.ma.masked_values(np.asarray(cMF.elev), cMF.hnoflo, atol=0.09)
    cMF.top = cMF.elev - np.ma.masked_values(gridSOILthick, cMF.hnoflo, atol=0.09)
    cMF.botm = np.asarray(cMF.botm)
    for L in range(cMF.nlay):
        cMF.botm[L] = np.ma.masked_values(cMF.botm[L], cMF.hnoflo, atol=0.09) - \
            np.ma.masked_values(gridSOILthick, cMF.hnoflo, atol=0.09)
    botm_l0 = np.asarray(cMF.botm)[0]
    for L in range(cMF.nlay):
        cMF.iuzfbnd[cMF.ibound[L] <= 0] = 0

    if nsp:
        cMF.nper = min(nsp, cMF.nper)
        cMF.perlen = np.asarray(cMF.perlen)[:cMF.nper]
        cMF.nstp = np.asarray(cMF.nstp)[:cMF.nper]

    # MF6 semantics: no hdry sentinel; MMsoil switches to h < botm dryness
    cMF.hdry = None

    mm = MMsoil.clsMMsoil(hnoflo=cMF.hnoflo)
    cells = mm.build_cell_list(cMF)
    # Phase 4: cell geometry provider (DIS = legacy delr/delc; DISV = polygons)
    from marmites_grid import geometry_for
    geom = geometry_for(cMF, cells, grid=grid)
    ctx = mm.build_context(cMF, cells, _nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _S_ini,
                           botm_l0, _Ks, gridSOIL, gridSOILthick, cMF.elev * 1000.0, gridMETEO,
                           INDEX_MM, INDEX_MM_SOIL, gridSsurfhmax, gridSsurfw,
                           P_veg_zoneSP, Eo_zonesSP, PT_veg_zonesSP, Pe_veg_zonesSP, PE_zonesSP,
                           gridVEGarea, LAI_veg_zonesSP, Zr, kTg_min, kTg_max, kT_f, kT_s, NVEG,
                           conv_fact, 1, P_irr_zoneSP, PT_irr_zonesSP, Pe_irr_zoneSP,
                           crop_irr_SP, gridIRR, Zr_c, kTg_min_c, kTg_max_c, kT_f_c, kT_s_c,
                           geom=geom)
    state = mm.init_state(ctx)
    top = np.asarray(cMF.top, dtype=float)
    botm = np.asarray(cMF.botm, dtype=float)
    return cMF, mm, ctx, state, top, botm, conv_fact


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--build-only', action='store_true')
    ap.add_argument('--libmf6', default=None, help='path to libmf6.dll/.so (enables the coupled run)')
    ap.add_argument('--mode', choices=['lagged', 'iterative'], default='lagged')
    ap.add_argument('--relax', type=float, default=0.6)
    ap.add_argument('--nsp', type=int, default=None, help='truncate to N stress periods')
    ap.add_argument('--daily', action='store_true', default=True)
    ap.add_argument('--aggregated', dest='daily', action='store_false')
    ap.add_argument('--nlay', type=int, default=None, choices=[2, 6],
                    help='vertical resolution: 2 reads the 2-layer parameter file '
                         '(_2s1L.ini) directly, 6 (default) the 6-layer file')
    ap.add_argument('--aggregate', action='store_true',
                    help='derive the 2-layer model by collapsing the 6-layer one '
                         'instead of reading _2s1L.ini (comparison only)')
    ap.add_argument('--grid', choices=['dis', 'disv'], default='dis',
                    help="'disv' builds the vertex-grid equivalent (Phase 4)")
    ap.add_argument('--strt-dem', nargs=2, type=float, metavar=('A', 'B'), default=None,
                    help='initial heads from the DEM: head = A*elevation + B '
                         '(e.g. --strt-dem 0.9995 -2.0); default uses the ini arrays')
    ap.add_argument('--sfr', action='store_true',
                    help='build the SFR stream network from inputPONDw.asc '
                         '(the SFR outlet then replaces the outlet DRN cells)')
    ap.add_argument('--sfr-rhk', type=float, default=0.1,
                    help='streambed hydraulic conductivity [m/d]')
    ap.add_argument('--lak', nargs='?', const='GIS/lm_ponds.shp', default=None,
                    metavar='SHAPEFILE',
                    help='build EMBEDDEDV lakes from a pond-polygon shapefile '
                         '(default GIS/lm_ponds.shp, relative to the dataset)')
    ap.add_argument('--lak-bedleak', type=float, default=1e-3,
                    help='lakebed leakance [1/d]')
    ap.add_argument('--uzf-vks-scale', type=float, default=1.0, metavar='F',
                    help='multiply the UZF unsaturated vks by F to offset the '
                         'EPSILON 2.0->3.5 clamp (UZF6 minimum), which otherwise '
                         'throttles recharge and drains the water table. '
                         'Calibrate so UZF-GWRCH matches the target (try ~3)')
    ap.add_argument('--seep', choices=['uzf', 'drn'], default='uzf',
                    help='groundwater-seepage mechanism: uzf = UZF SIMULATE_GWSEEP, '
                         'drn = smoothed land-surface drain')
    ap.add_argument('--seep-cond', type=float, default=10000.0,
                    help='DRN-SEEP conductance per cell [m2/d] (--seep drn only). '
                         'This is a free-draining seepage face, not a physical '
                         'aquitard: too small and the water table stands above '
                         'the land surface instead of discharging')
    ap.add_argument('--no-ats', dest='ats', action='store_false', default=True,
                    help='disable adaptive time stepping (not recommended: a '
                         'stress period MF6 cannot solve in one step then '
                         'silently becomes a non-solution)')
    ap.add_argument('--max-discrepancy', type=float, default=1.0,
                    help='fail the run above this cumulative mass-balance '
                         'discrepancy [%%]')
    ap.add_argument('--allow-bad-budget', action='store_true',
                    help='report instead of failing when the solution check '
                         'does not pass')
    ap.add_argument('--strt-heads', default=None, metavar='PREFIX',
                    help='use a saved per-layer head field <PREFIX>_l1.asc.. as '
                         'the initial condition (e.g. an equilibrated spin-up '
                         'result); overrides --strt-dem')
    ap.add_argument('--save-strt', nargs='?', const='hi_spinup', default=None,
                    metavar='PREFIX',
                    help='after the run, save the final head field as '
                         '<PREFIX>_l1.asc.. in the MF workspace (default prefix '
                         'hi_spinup) so it can seed later runs via --strt-heads')
    ap.add_argument('--steady-means', default=None, metavar='PREFIX',
                    help='drive the steady SP0 with saved per-cell mean recharge '
                         'and ETg (<PREFIX>_perc.asc, <PREFIX>_etg.asc) instead '
                         'of uniform perc_user, so SP0 lands near the dynamic '
                         'equilibrium')
    ap.add_argument('--save-means', nargs='?', const='hi_spinup', default=None,
                    metavar='PREFIX',
                    help='after the run, save per-cell mean recharge/ETg as '
                         '<PREFIX>_perc.asc / _etg.asc for reuse via --steady-means')
    ap.add_argument('--spinup', type=int, default=1, metavar='N',
                    help='repeat the forcing up to N times, each cycle starting '
                         'from the previous final heads, until the water table '
                         'stabilises (equilibrium IC). Default 1 = single run')
    ap.add_argument('--spinup-tol', type=float, default=0.05, metavar='M',
                    help='spin-up converged when the mean between-cycle water '
                         'table change is below this [m] (default 0.05)')
    ap.add_argument('--postproc', action='store_true',
                    help='after the run, write obs/head/budget figures + CSVs '
                         'to <ws>/postproc/')
    ap.add_argument('--preproc', action='store_true',
                    help='write the input parameter maps to <out-dir>/_input/')
    ap.add_argument('--gis-ws', default=None,
                    help='GIS workspace for the general map (default: '
                         'MARMITES_GIS_WS, else the built-in path)')
    ap.add_argument('--postproc-only', action='store_true',
                    help='re-draw the figures from a run already on disk: reads '
                         '<ws>/_coupled_<mode>.h5 and the MF6 output, runs no '
                         'MODFLOW (implies --postproc)')
    ap.add_argument('--sankey-min-flux', type=float, default=0.05, metavar='MM',
                    help='native water-balance Sankey: hide flows below this '
                         'magnitude (mm/y) on the core diagram (default 0.05)')
    ap.add_argument('--sankey-obs-years', action='store_true',
                    help='also produce a per-hydrological-year Sankey at each '
                         'observation point (default: whole period only; the '
                         'catchment always gets both)')
    ap.add_argument('--map-days', type=int, default=6, metavar='N',
                    help='head maps on N evenly spaced days (0 = mean only)')
    ap.add_argument('--no-sankey-full', dest='sankey_full', action='store_false',
                    default=True, help='skip the all-flux Sankey (keep only the '
                         'decluttered core diagram)')
    ap.add_argument('--standalone', default=None, metavar='MF6EXE',
                    help='build, then run the model with the mf6 executable '
                         '(no API) -- isolates model faults from coupling faults')
    ap.add_argument('--probe', action='store_true',
                    help='list the MF6 memory variables and exit (diagnostics)')
    ap.add_argument('--ws', default=None,
                    help='MODFLOW 6 workspace (default <ws-root>/MF6_ws)')
    ap.add_argument('--ws-root', default=WS_ROOT, metavar='DIR',
                    help='root for ALL run output, outside the repo '
                         '(default %s, or $MARMITES_WS_ROOT)' % WS_ROOT)
    ap.add_argument('--run-tag', default=None, metavar='TAG',
                    help='label for this run\'s results folder '
                         '<ws-root>/out_<YYYYMMDDHHMM>_<TAG>')
    a = ap.parse_args()
    if a.postproc_only:
        a.postproc = True
    if a.ws is None:
        a.ws = os.path.join(a.ws_root, 'MF6_ws' if a.grid == 'dis' else 'MF6_ws_disv')
    os.makedirs(a.ws, exist_ok=True)
    # results folder for this run, in the legacy out_<stamp>_<tag> style
    tag = a.run_tag or ('%dlay_%s' % (a.nlay or 6, a.mode))
    a.out_dir = os.path.join(a.ws_root,
                             'out_%s_%s' % (time.strftime('%Y%m%d%H%M'), tag))
    # Saved run state (equilibrated heads, steady means) is run OUTPUT, so it is
    # written to the workspace; reading falls back to the baseline committed in
    # the repo's DataSet_LaMata/MF_ws so `--strt-heads hi_spinup` keeps working.
    a.state_dir = a.ws

    cMF, mm, ctx, state, top, botm, conv_fact = setup_lamata(
        daily=a.daily, nsp=a.nsp, grid=a.grid, nlay=a.nlay, aggregate=a.aggregate)
    if a.postproc_only:
        # Re-draw the figures from a run that already happened: everything the
        # post-processing reads is on disk (the coupled HDF5 for the MM side,
        # the .hds/.cbc/.grb for the aquifer side), so MODFLOW need not run
        # again. Iterating on a figure costs seconds instead of the full run.
        h5_fn = os.path.join(a.ws, '_coupled_%s.h5' % a.mode)
        if not os.path.exists(h5_fn):
            raise SystemExit('--postproc-only needs a previous run: %s not found'
                             % h5_fn)
        with h5py.File(h5_fn, 'r') as f:
            res = {k: f[k][:] for k in f.keys()}
        print('re-using %s (%d stress period(s))'
              % (h5_fn, res['wb_ts'].shape[0]))
        _run_postproc(a, cMF, ctx, res)
        return
    b = clsMF6(cMF, top=top, botm=botm, sim_ws=a.ws, daily=True, grid=a.grid,
               strt_from_dem=(tuple(a.strt_dem) if a.strt_dem else None))
    b.seep = a.seep
    b.ats = a.ats
    b.drn_seep_cond = float(a.seep_cond)
    b.uzf_vks_scale = float(a.uzf_vks_scale)
    if a.sfr:
        # the channel map MARMITES already uses for its ponds/surface storage
        b.sfr_pondw = _asc(os.path.join(DS, 'inputPONDw.asc'))
        b.sfr_pondhmax = _asc(os.path.join(DS, 'inputPONDhmax.asc'))
        b.sfr_rhk = float(a.sfr_rhk)
    if a.lak:
        shp = a.lak if os.path.isabs(a.lak) else os.path.join(DS, a.lak)
        b.lak_shapefile = shp
        b.lak_bedleak = float(a.lak_bedleak)
        b.lak_depth = _asc(os.path.join(DS, 'inputPONDhmax.asc'))
    if a.strt_heads:
        # a saved (equilibrated) head field seeds the IC directly, so the
        # spin-up need not be repeated. Resolve relative to the MF workspace.
        pref = _state_in(a, a.strt_heads, '_l1.asc')
        b.strt_array = b.load_heads_asc(pref)
        print('initial heads loaded from %s_l*.asc' % pref)
    b.build()
    b.write()
    print('MF6 (%s) simulation written to %s  (%d SPs incl. steady, %d UZF cells, %d wells)'
          % (a.grid.upper(), a.ws, b.nper, b.nuzfcells, b.ncell))
    if b.seep == 'drn':
        print('   seepage: DRN_SEEP, %d drains, cond %.4g m2/d, ddrn %.4g m '
              '(UZF SIMULATE_GWSEEP off)'
              % (b.ndrnseep, b.drn_seep_cond, b.drn_seep_ddrn))
    else:
        print('   seepage: UZF SIMULATE_GWSEEP')

    if a.standalone:
        # Decisive discriminator: run the SAME model with the mf6 executable,
        # no Python API involved. If this crashes, the fault is in the model
        # we generate; if it runs, the fault is in the coupling exchange.
        import subprocess
        exe = os.path.abspath(a.standalone.strip('"').strip("'"))
        if not os.path.isfile(exe):
            sys.exit('mf6 executable not found: %s' % exe)
        ws = os.path.abspath(a.ws)
        print('\nRunning MF6 standalone in %s\n%s\n' % (ws, '-' * 60))
        proc = subprocess.run([exe], cwd=ws, capture_output=True, text=True)
        print(proc.stdout[-4000:] if proc.stdout else '(no stdout)')
        if proc.stderr:
            print('--- stderr ---\n%s' % proc.stderr[-2000:])
        print('%s\nmf6 exit code: %d  (0 = normal termination)' % ('-' * 60, proc.returncode))
        lst = os.path.join(ws, 'mfsim.lst')
        if os.path.exists(lst):
            with open(lst, errors='replace') as f:
                tail = f.readlines()[-25:]
            print('\n--- mfsim.lst (tail) ---\n%s' % ''.join(tail))
        return


    if a.build_only or not a.libmf6:
        # preproc only needs the built model, so it can run without libmf6;
        # postproc needs MF6 output, so it is skipped here
        if a.preproc:
            from marmites_postprocess import run_preproc
            run_preproc(a.ws, DS, name=cMF.modelname.lower(),
                        gis_ws=a.gis_ws)
        print('\nNo --libmf6 given: stopping after build.\n'
              'To run coupled:  python tests/run_lamata_mf6.py --libmf6 <path> --mode %s' % a.mode)
        return

    # --- robust libmf6 handling ---------------------------------------
    # WinError 87 from CDLL(winmode=0x08) happens when the path is not
    # fully qualified (relative, or polluted with quotes by IDE run
    # configs); dependency DLLs also need the bin dir on the search path.
    lib = os.path.abspath(a.libmf6.strip('"').strip("'"))
    if not os.path.isfile(lib):
        sys.exit('libmf6 not found at: %s' % lib)
    if hasattr(os, 'add_dll_directory'):
        os.add_dll_directory(os.path.dirname(lib))
    from modflowapi import ModflowApi
    api = ModflowApi(lib, working_directory=os.path.abspath(a.ws))

    if a.probe:
        # List what this MF6 build exposes, then stop. Use this when the
        # coupler cannot bind an address (MF6 renames memory variables
        # between versions).
        api.initialize()
        try:
            names = list(api.get_input_var_names())
        finally:
            try:
                api.finalize()
            except Exception:
                pass
        out = os.path.join(os.path.abspath(a.ws), 'mf6_variables.txt')
        with open(out, 'w') as f:
            f.write('\n'.join(names))
        print('MF6 exposes %d variables -> %s' % (len(names), out))
        for key in ('/UZF', '/WEL', '/DIS', '/X'):
            sel = [n for n in names if key in n.upper()]
            print('\n%s (%d):' % (key, len(sel)))
            for n in sel[:25]:
                print('   ', n)
        return
    import flopy
    hds_fn = os.path.join(a.ws, cMF.modelname.lower() + '.hds')

    def _final_heads():
        hf = flopy.utils.HeadFile(hds_fn)
        H = hf.get_data(kstpkper=hf.get_kstpkper()[-1])   # (nlay, nrow, ncol)
        return np.where(np.abs(H) > 1e29, np.nan, np.asarray(H, dtype=float))

    # Spin-up: repeat the whole forcing, each cycle starting from the previous
    # cycle's FINAL heads, until the water table stops moving between cycles.
    # This turns the DEM-regression start (a convenient but non-equilibrium IC,
    # too high in the uplands) into an equilibrated state before the reported
    # run. A single run is just ncyc = 1.
    # mean recharge / ETg to drive the steady SP0 near dynamic equilibrium
    steady_perc = steady_etg = None
    if a.steady_means:
        mp = _state_in(a, a.steady_means, '_perc.asc')
        steady_perc = _read_cell_grid(mp + '_perc.asc', ctx.cells)
        steady_etg = _read_cell_grid(mp + '_etg.asc', ctx.cells)
        print('steady-state means loaded from %s_{perc,etg}.asc' % mp)

    # observation cells at which to keep the full MM flux series (per-point
    # Sankey / time series). Resolved once; the coupler captures them each cycle.
    obs_idx, obs_names = [], []
    if a.postproc:
        try:
            from marmites_postprocess import resolve_obs_cells
            obs_idx, obs_names = resolve_obs_cells(cMF, ctx, DS)
        except Exception as exc:
            print('obs-cell resolution skipped: %r' % exc)

    ncyc = max(1, int(a.spinup))
    prev = None
    # Bound at the end of every cycle below. Named here because the cyc > 0
    # branch reads it: the loop always runs at least once, so it is never
    # actually unbound there, but nothing in the block says so and pyflakes
    # reports it as an undefined name.
    prev_heads = None
    for cyc in range(ncyc):
        if cyc > 0:
            b.strt_array = prev_heads      # equilibrating IC from last cycle
            b.build()
            b.write()
        st = mm.init_state(ctx)            # fresh soil state each cycle
        cpl = MF6Coupler(mm, ctx, st, b, conv_fact=conv_fact,
                         mode=a.mode, relax=a.relax,
                         obs_idx=obs_idx, obs_names=obs_names)
        cpl.steady_perc, cpl.steady_etg = steady_perc, steady_etg
        res = cpl.run(api)
        # A non-converged / non-conserving cycle is not a usable state to
        # iterate from, so the guard runs every cycle.
        cpl.check_solution(max_discrepancy=a.max_discrepancy,
                           raise_on_fail=not a.allow_bad_budget)
        prev_heads = _final_heads()
        # Feed THIS cycle's mean recharge/ETg into the next cycle's steady SP0.
        # Carrying heads alone does not equilibrate (a steady period ignores
        # STRT); driving SP0 with the dynamic mean is what makes the spin-up
        # actually converge. Skipped if the user pinned the means explicitly.
        if not a.steady_means:
            steady_perc = res['perc'].mean(axis=0)
            steady_etg = res['etg'].mean(axis=0)
        if ncyc > 1:
            wt = np.nanmax(prev_heads, axis=0)          # water table per column
            if prev is not None:
                delta = float(np.nanmean(np.abs(wt - prev)))
                print('spin-up cycle %d/%d: mean |dWT| vs previous = %.3f m '
                      '(tol %.3f)' % (cyc + 1, ncyc, delta, a.spinup_tol))
                if delta < a.spinup_tol:
                    print('spin-up converged after %d cycle(s).' % (cyc + 1))
                    break
            else:
                print('spin-up cycle 1/%d done (baseline).' % ncyc)
            prev = wt
    check = None

    out_fn = os.path.join(a.ws, '_coupled_%s.h5' % a.mode)
    with h5py.File(out_fn, 'w') as f:
        for k, v in res.items():
            f.create_dataset(k, data=v)
        f.create_dataset('cell_ij', data=np.array([(c[1], c[2]) for c in ctx.cells]))
        # true grid size: cannot be inferred from active cells alone
        f.create_dataset('grid_shape', data=np.array([cMF.nrow, cMF.ncol]))
    print('\nCoupled run finished. Results: %s' % out_fn)
    print('perc  mean %.4g m/d   ETg mean %.4g m/d   outer iters mean %.1f'
          % (res['perc'].mean(), res['etg'].mean(), res['outer_iters'].mean()))

    # aquifer recharge/discharge balance -- the number to watch when calibrating
    # --uzf-vks-scale: recharge reaching the water table should ~match discharge
    try:
        import flopy
        _lst = flopy.utils.Mf6ListBudget(os.path.join(a.ws, cMF.modelname.lower() + '.lst'))
        _r = _lst.get_dataframes()[0].iloc[1:]
        _A = float(np.sum(b.idomain[0] > 0)) * float(np.mean(cMF.delr)) * float(np.mean(cMF.delc))
        _mmyr = lambda c: float(_r.get(c, 0).mean()) * 365.0 / _A * 1000.0
        rch = _mmyr('UZF-GWRCH_IN')
        dis = _mmyr('DRN2_OUT') + _mmyr('WEL_OUT') + _mmyr('DRN_OUT')
        print('aquifer balance: recharge to WT %.1f mm/yr  vs  discharge %.1f '
              'mm/yr  (deficit %.1f)' % (rch, dis, dis - rch))
        if dis - rch > 5.0:
            print('   -> still draining; raise --uzf-vks-scale (currently %.3g) '
                  'to lift recharge' % a.uzf_vks_scale)
    except Exception as _exc:
        pass

    # Save the final head field for reuse as an IC. Auto-save after a spin-up
    # (so it is never lost), or on explicit --save-strt for a single run.
    save_pref = a.save_strt or ('hi_spinup' if ncyc > 1 else None)
    if save_pref:
        pref = _state_out(a, save_pref)
        paths = b.save_heads_asc(prev_heads, pref)
        print('equilibrated heads saved: %s' % ', '.join(os.path.basename(p) for p in paths))
        print('   reuse with:  --strt-heads %s   (skips the spin-up)' % save_pref)

    # Save per-cell mean recharge / ETg so the steady state of later runs can be
    # driven by the dynamic mean (auto after spin-up, or on explicit --save-means).
    mean_pref = a.save_means or ('hi_spinup' if ncyc > 1 else None)
    if mean_pref:
        mp = _state_out(a, mean_pref)
        _write_cell_grid(res['perc'].mean(axis=0), ctx.cells, cMF.nrow, cMF.ncol,
                         cMF, mp + '_perc.asc')
        _write_cell_grid(res['etg'].mean(axis=0), ctx.cells, cMF.nrow, cMF.ncol,
                         cMF, mp + '_etg.asc')
        print('steady-state means saved: %s_{perc,etg}.asc' % os.path.basename(mp))
        print('   reuse with:  --steady-means %s' % mean_pref)

    _run_postproc(a, cMF, ctx, res)


def _run_postproc(a, cMF, ctx, res):
    """Draw every figure for a completed run.

    Shared by the normal path and by --postproc-only, so re-drawing from an
    existing run takes exactly the same route as drawing at the end of one.
    """
    if a.postproc or a.preproc:
        # All results go to <ws-root>/out_<stamp>_<tag>/, never into the
        # repository and not into the model workspace either.
        from marmites_postprocess import run_preproc, run_postproc, native_suite
        os.makedirs(a.out_dir, exist_ok=True)
        if a.preproc:
            run_preproc(a.ws, DS, name=cMF.modelname.lower(), out_root=a.out_dir,
                        gis_ws=a.gis_ws,
                        cMF=cMF, ctx=ctx, res=res)
        if a.postproc:
            run_postproc(a.ws, DS, name=cMF.modelname.lower(), out_root=a.out_dir)
            # native MARMITESplot figures, driven by the in-memory coupled data
            native_suite(os.path.join(a.out_dir, '_output'), cMF, ctx, res,
                         ds_ws=DS, sim_ws=a.ws, sankey=True,
                         sankey_full=a.sankey_full,
                         sankey_min_flux=a.sankey_min_flux, map_days=a.map_days,
                         sankey_obs_years=a.sankey_obs_years)
            # 01-07 water-budget figures incl. 06_heads/07_coupling and the
            # NWT-vs-MF6 comparison (into <out_dir>/figures_nwt_comparison/)
            try:
                import plot_water_budget as pwb
                pwb.make_figures(a.ws, mode=a.mode, out_dir=a.out_dir)
            except Exception as exc:
                print('   plot_water_budget skipped: %r' % exc)
        print('results written to %s' % a.out_dir)


if __name__ == '__main__':
    main()
