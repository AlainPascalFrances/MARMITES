# -*- coding: utf-8 -*-
"""La Mata MODFLOW 6 coupled run (Phase 3 entry point).

Builds the MF6 simulation from the La Mata dataset and, when the MODFLOW 6
library is available, drives the coupled MARMITES-MF6 model through the API
(lagged: MMsoil once per stress period, then MF6's own time step).
Without libmf6 it stops after writing the simulation files (still useful:
run `mf6` manually in the workspace to check the groundwater model alone).

Requirements on the executing machine:
    pip install flopy modflowapi
    MODFLOW 6 >= 6.4 binaries:  mf6[.exe] and libmf6[.dll|.so]
    (easiest: `get-modflow :flopy` or download from
     github.com/MODFLOW-USGS/executables)

Usage:
    python tests/run_lamata_mf6.py --config configs/lamata.toml
    python tests/run_lamata_mf6.py --config configs/lamata.toml --set run.build_only=true
    python tests/run_lamata_mf6.py --config configs/lamata.toml --set run.nsp=60 --run-tag try
"""
import argparse
import dataclasses
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
# WS_ROOT is mm_paths.WS_ROOT (panel 0, MM_WS_ROOT). It had a literal default
# of its own here, E:\00code_ws\LaMata_MM-MF6, which the run did not read
# but a diagnostic script did (2026-10-07).
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
import marmites_props as props  # noqa: E402
from marmites_coupler import MF6Coupler  # noqa: E402
import marmites_config as mcfg  # noqa: E402
import marmites_surface as msurf  # noqa: E402
import mm_paths  # noqa: E402

WS_ROOT = str(mm_paths.WS_ROOT)


# The BMI/XMI SHARED LIBRARY, not the executable. The coupler steps MODFLOW
# one stress period at a time from inside python and exchanges arrays with it
# between steps; mf6.exe can only run a whole simulation start to finish,
# which is the one thing the coupling cannot use. Same build, two products,
# side by side in the same bin folder -- so pointing at the wrong one is easy
# and has to be said rather than reported as "not found".
LIB_NAMES = ('libmf6.dll', 'libmf6.so', 'libmf6.dylib')


class LibMF6Error(Exception):
    """The MODFLOW 6 library cannot be used as given."""


def resolve_libmf6(given):
    """``paths.libmf6`` -> an absolute library file, or ''. Never raises.

    Blank stops after the build, by design. Everything else is completed
    here rather than at the API call, which happens AFTER the whole model is
    built -- a typo should not cost a build. A FOLDER is completed with the
    library name, because the natural thing to paste is the bin directory.
    """
    given = (given or '').strip().strip('"').strip("'")
    if not given:
        return ''
    if given.lower() == 'auto':
        return mm_paths.LIBMF6 if os.path.isfile(mm_paths.LIBMF6) else ''
    path = os.path.abspath(given)
    if os.path.isdir(path):
        for name in LIB_NAMES:
            cand = os.path.join(path, name)
            if os.path.isfile(cand):
                return cand
        return path          # let check_libmf6 say what is wrong with it
    return path


def check_libmf6(path):
    """Say what is wrong with it, in the words that lead to the fix.

    Returns the path. Raises :class:`LibMF6Error` naming the actual mistake:
    a directory with no library in it, the executable instead of the
    library, or a path that is simply not there.
    """
    if not path:
        return path
    if os.path.isdir(path):
        raise LibMF6Error(
            '%s is a folder and holds no %s. paths.libmf6 is the MODFLOW 6 '
            'LIBRARY, not the folder and not mf6.exe.'
            % (path, ' / '.join(LIB_NAMES)))
    base = os.path.basename(path).lower()
    if not os.path.isfile(path):
        raise LibMF6Error(
            'libmf6 not found at: %s\nSet paths.libmf6 to "auto" to use %s.'
            % (path, mm_paths.LIBMF6))
    if base.startswith('mf6') and base.endswith('.exe'):
        lib = os.path.join(os.path.dirname(path), 'libmf6.dll')
        raise LibMF6Error(
            '%s is the EXECUTABLE. The coupler drives MODFLOW through the '
            'API one stress period at a time, which the executable cannot '
            'do; it needs the shared library beside it%s.'
            % (path, ' -- %s' % lib if os.path.isfile(lib) else ''))
    return path


def _forcing(cfg):
    """The daily forcing: run MMsurf, or check what is already there (WP1d).

    ``run.surface = 1``   MMsurf runs from the configuration and writes the
                          series into the WORKSPACE -- they are run output.
    ``run.surface = 0``   the series must already exist; they are checked for
                          presence AND shape, and the run stops if they are
                          not there rather than proceeding on a stale file.

    With the switch off the series are read from the dataset, which is where
    the committed La Mata forcing lives, so nothing moves until MMsurf is
    actually used.
    """
    cfg = cfg or mcfg.RunConfig.from_dict({})
    if cfg.run.surface:
        # SAID BEFORE IT STARTS, and naming the key. MMsurf is the noisiest
        # thing in the log and the longest part of a build, so "why is this
        # running when I turned it off" has to be answerable from the log:
        # the switch is a panel widget, and what runs is the saved file.
        print('run.surface is ON: MMsurf runs and writes the daily forcing. '
              'Set [run] surface = false to use the series already there.')
        out_ws = msurf.surface_ws(cfg, mm_paths.WS_ROOT, cfg.paths.case)
        return msurf.run(cfg, DS, out_ws, config_hash=cfg.config_hash())
    spec = msurf.forcing_spec(cfg, DS)
    ndays = msurf.check_forcing(spec, must_exist=True,
                                nper=cfg.run.nsp or None)
    print('forcing: %d day(s) from %s (run.surface is off)' % (ndays, DS))
    return spec



def _apply_dem(cMF, cfg, dataset_dir, cache_dir=None):
    """Take the land surface from the DEM panel 1 names. ``(applied, note)``.

    ``[grid] dem`` names the raster on panel 1 and the converter copies it
    into the dataset. It is REQUIRED (2026-10-07): marmites_props.
    land_surface already took the surface from it on the dataset grid, and
    there is no parameter-file elevation raster to fall back on any more.

    Wrapped onto THE GRID THIS RUN USES, after any mesh projection: on the
    dataset grid itself this changes nothing, and on a mesh it puts the
    survey on the mesh cells instead of the 50 m average projected there.

    THE SURFACE MOVES AND THE THICKNESSES DO NOT: top and every botm are
    shifted by the same delta as elev. The DEM refines where the ground is;
    it says nothing about how thick the aquifer below it is, and re-deriving
    the layer geometry from it would invent an answer the raster does not
    hold. Cells the DEM does not reach keep the elevation they had.
    """
    import marmites_dem as mdem

    if not getattr(cfg, 'grid', None) or not cfg.grid.dem:
        return False, ('no [grid] dem: the land surface stays as the dataset '
                       'grid had it')
    path = mdem.dem_path(dataset_dir)
    if not os.path.exists(path):
        return False, ('[grid] dem is %r but %s is not in the dataset -- run '
                       'code/tools/gis_to_dataset.py.'
                       % (cfg.grid.dem, os.path.basename(path)))
    gp = getattr(cMF, 'mesh_gridprops', None)
    if gp is None:
        from marmites_grid import disv_from_structured
        verts, cell2d, ncpl = disv_from_structured(
            cMF.delr, cMF.delc, float(getattr(cMF, 'xllcorner', 0.0)),
            float(getattr(cMF, 'yllcorner', 0.0)))
        gp = {'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl,
              'nlay': int(cMF.nlay)}
    wrapped, info = mdem.wrap_to_grid(path, gp, cache_dir=cache_dir,
                                      warn=lambda m: print('WARNING: %s' % m))

    shape = np.asarray(cMF.elev).shape
    new = np.ma.masked_values(np.ma.filled(wrapped, -9999.0).reshape(shape),
                              -9999.0, atol=1e-6)
    old = np.ma.masked_values(np.asarray(cMF.elev), cMF.hnoflo, atol=0.09)
    both = (~np.ma.getmaskarray(new)) & (~np.ma.getmaskarray(old))
    delta = np.zeros(shape, dtype=float)
    delta[both] = np.ma.getdata(new)[both] - np.ma.getdata(old)[both]

    cMF.elev = np.ma.array(np.where(both, np.ma.getdata(new),
                                    np.ma.getdata(old)),
                           mask=np.ma.getmaskarray(old))
    cMF.top = np.ma.array(np.ma.getdata(cMF.top) + delta,
                          mask=np.ma.getmaskarray(cMF.top))
    botm = np.asarray(cMF.botm, dtype=float)
    for L in range(int(cMF.nlay)):
        botm[L] = botm[L] + delta
    cMF.botm = botm

    # ANYTHING ANCHORED TO THE LAYER GEOMETRY MOVES WITH IT. A DRN elevation
    # under `drn.at_layer_base` is the layer bottom plus 10 mm -- it is a
    # position IN the stack, not a height above sea level -- and the stack
    # has just moved. Leaving the drains behind put 2928 of La Mata's 6059
    # below their own cell bottom, by up to 4.72 m where the DEM lowered the
    # surface most, and MODFLOW refused the whole package:
    #
    #   DRN BOUNDARY (2178) ELEVATION (771.773) IS LESS THAN CELL BOTTOM
    #
    # aborting the process from inside the library. GHB is deliberately NOT
    # shifted: its HEAD is a boundary condition on the water table, an
    # absolute elevation that means the same thing wherever the layer sits.
    _moved = 0
    _recs = getattr(cMF, 'layer_row_column_elevation_cond', None)
    for _spd in ((_recs or {}).values() if isinstance(_recs, dict)
                 else (_recs or ())):
        for _r in _spd:
            _i, _j = int(_r[1]), int(_r[2])
            try:
                _d = float(delta[_i, _j])
            except (IndexError, ValueError):
                continue
            if _d:
                _r[3] = float(_r[3]) + _d
                _moved += 1
    if _moved:
        print('   %d drain elevation(s) moved with the land surface, so they '
              'keep their height above their own layer bottom' % _moved)

    d = delta[both]
    return True, ('land surface from %s (%g m), wrapped onto %d %s cell(s)%s: '
                  'moved by mean %+.3f m, rms %.3f m, |max| %.2f m over %d '
                  'cell(s); %d cell(s) not covered kept their elevation'
                  % (cfg.grid.dem, info['cellsize'], info['ncpl'],
                     cfg.grid_kind, ' (cached)' if info['cached'] else '',
                     d.mean() if d.size else 0.0,
                     float(np.sqrt((d * d).mean())) if d.size else 0.0,
                     float(np.abs(d).max()) if d.size else 0.0, int(both.sum()),
                     int(both.size - both.sum())))


def _crr_options(a, cMF, dataset_dir, cache_dir=None):
    """The CRR cascade's settings for the coupler (WP5), or None when off.

    The slopes come from ``[crr] dem``. When that is the raster the land
    surface was already taken from (the default: both are the dataset's
    inputDEM.asc, the sink-filled survey), the cascade uses that surface
    as it stands -- same file, same wrapping. Any other raster is wrapped
    onto the run's grid the same way, at its own resolution.
    """
    if not getattr(a, 'crr', False):
        return None
    opts = {'beta': float(a.crr_beta), 'sinks': str(a.crr_sinks)}
    name = (getattr(a, 'crr_dem', None) or '').strip()
    path = (name if os.path.isabs(name) else
            os.path.join(str(dataset_dir), name)) if name else None
    same = bool(path is not None and getattr(cMF, 'land_dem', None) and
                os.path.normcase(os.path.abspath(path)) ==
                os.path.normcase(os.path.abspath(cMF.land_dem)))
    if path is None or same:
        print('CRR: slopes from the land surface%s, beta %g, sinks %s -- '
              'from the panel' % (' (%s)' % name if name else '',
                                  opts['beta'], opts['sinks']))
        return opts
    if not os.path.exists(path):
        raise SystemExit('CONFIG ERROR: [crr] dem = %r is not in the dataset '
                         '(%s). Name the sink-filled DEM the cascade should '
                         'follow, or leave it as inputDEM.asc.' % (name, path))
    import marmites_dem as mdem
    gp = getattr(cMF, 'mesh_gridprops', None)
    if gp is None:
        from marmites_grid import disv_from_structured
        verts, cell2d, ncpl = disv_from_structured(
            cMF.delr, cMF.delc, float(getattr(cMF, 'xllcorner', 0.0)),
            float(getattr(cMF, 'yllcorner', 0.0)))
        gp = {'vertices': verts, 'cell2d': cell2d, 'ncpl': ncpl,
              'nlay': int(cMF.nlay)}
    wrapped, info = mdem.wrap_to_grid(path, gp, cache_dir=cache_dir,
                                      warn=lambda m: print('WARNING: %s' % m))
    elev = np.asarray(np.ma.filled(np.ma.asarray(wrapped, dtype=float),
                                   np.nan), dtype=float)
    opts['elev'] = elev.reshape(np.shape(cMF.elev))
    print('CRR: slopes from %s (%g m), wrapped onto %d cell(s)%s; beta %g, '
          'sinks %s -- from the panel'
          % (name, info['cellsize'], info['ncpl'],
             ' (cached)' if info['cached'] else '', opts['beta'],
             opts['sinks']))
    return opts


def ctx_geom_area(cMF):
    """Mean cell area [m2], DIS or DISV. The drainage law w = a*A**b needs it."""
    proj = getattr(cMF, 'mesh_proj', None)
    if proj is not None:
        import marmites_vector as mv
        return mv.TargetGrid.from_gridprops(cMF.mesh_gridprops).area
    return np.outer(np.asarray(cMF.delc, float), np.asarray(cMF.delr, float))


def extdp_grid(cfg, cMF, ctx, dataset_dir, verbose=True):
    """UZF's extinction depth ON THE GRID THE RUN USES (WP2 2.2).

    `et.extdp` as given -- one value per layer, or a raster read and, on a
    mesh, area-averaged onto the cells. (resolve_source hands back raster
    PATHS, which the builder cannot broadcast: a raster never reached UZF.)

    With `et.extdp_from = 'vegetation'` each cell takes the cover-weighted
    rooting depth of its vegetation, and `et.extdp` -- at the surface
    layer -- is the depth of the fraction nothing covers. The cover and the
    depths are MMsoil's own (ctx), so the two sides see one community.
    """
    nlay = int(cMF.nlay)
    src = props.resolve_source(cfg.et.extdp, nlay, dataset_dir, 'et.extdp')
    if isinstance(src[0], str):
        proj = getattr(cMF, 'mesh_proj', None)
        lays = []
        for fn in src:
            arr = _asc(fn)
            if proj is not None:
                arr = np.ma.filled(np.asarray(
                    proj.sample2d(arr, fill=0.0, dtype=float, how='area'),
                    dtype=float), 0.0)
            lays.append(np.asarray(arr, dtype=float).reshape(
                int(cMF.nrow), int(cMF.ncol)))
        base = np.stack(lays)
    else:
        base = np.asarray(src, dtype=float)
    if cfg.et.extdp_from != 'vegetation':
        if verbose:
            print('UZF ET: always on (%s), extinction depth from %s'
                  % (cfg.et.unsat_form, cfg.et.extdp.producer()))
        return base
    # the bare fraction's depth: the source at each column's SURFACE layer
    top_lay = np.clip(np.asarray(cMF.outcropL, dtype=int) - 1, 0, nlay - 1)
    if base.ndim == 1:
        bare = base[top_lay]
    else:
        bare = np.take_along_axis(base, top_lay[None], axis=0)[0]
    act = np.asarray(cMF.outcropL) > 0
    areas = np.asarray(ctx_geom_area(cMF), dtype=float).reshape(act.shape)
    irr = getattr(ctx, 'gridIRR', None) if getattr(ctx, 'irr_yn', 0) else None
    # UZF's top is the land surface minus the soil column (setup_lamata),
    # so each root zone reaches UZF only below the soil
    soil = np.ma.filled(np.ma.masked_values(
        np.asarray(ctx.gridSOILthick, dtype=float), cMF.hnoflo, atol=0.09),
        0.0)
    ext, notes = props.extdp_by_vegetation(
        ctx.gridVEGarea, ctx.Zr, bare, soil_thick=soil, irr=irr,
        crop_by_period=getattr(ctx, 'crop_irr_SP', None),
        crop_root=getattr(ctx, 'Zr_c', None),
        perlen=np.asarray(cMF.perlen, dtype=float)[:int(cMF.nper)],
        areas=np.where(act, areas, 0.0))
    if verbose:
        print('UZF ET: always on (%s), %s' % (cfg.et.unsat_form, notes[0]))
        for n in notes[1:]:
            print(n)
    return ext


def active_area(cMF, b):
    """The catchment's area [m2]: the TRUE areas of the active surface cells.

    The balance line used to take (active cells in layer 1) x mean(delr) x
    mean(delc) -- right on the 50 m grid only. On the Voronoi mesh the cells
    run from ~1 m2 by the drains to ~4000 m2 in the far field, and the line
    printed a recharge of 2,263,794 mm/yr.
    """
    ar = np.asarray(ctx_geom_area(cMF), dtype=float)
    if ar.size != b.nrow * b.ncol:
        raise ValueError('%d cell areas for a %d x %d grid'
                         % (ar.size, b.nrow, b.ncol))
    ar = ar.reshape(b.nrow, b.ncol)
    return float(sum(ar[i, j] for i, j, _lay in b.surf_cells))


# The flows into and out of the AQUIFER, by budget-term prefix. Only the
# terms present in a run's list file count; the rest are simply absent.
RECHARGE_TERMS = ('UZF-GWRCH_IN',)
DISCHARGE_PREFIXES = ('DRN', 'GHB', 'WEL', 'EVT')
# the aquifer's exchange with the streams and the ponds, both ways: seepage
# from them into it (_IN) and groundwater exfiltrating into them (_OUT)
SURFACE_PREFIXES = ('SFR', 'LAK')


def _prefix(c):
    """'DRN_SEEP_OUT' -> 'DRN', 'SFR-1_IN' -> 'SFR', 'DRN2_OUT' -> 'DRN'."""
    return c.split('_')[0].split('-')[0].rstrip('0123456789')


def aquifer_balance(cum, times, area, steady_first=True):
    """Recharge to, and discharge from, the aquifer in mm/yr over the run.

    ``cum`` is the list file's CUMULATIVE volumes [m3], one row per saved
    step; ``times`` their times [d]. The first period is left out, as it
    always was (it is the initial state); the rates are volumes over elapsed
    time, NOT a plain mean of per-step rates -- which weighted a 1-day step
    as much as a 30-day one. Also returns the share of the recharge carried
    by the single largest step, and which step, because one pulse can make
    the whole-run figure meaningless.

    Without a steady first period (a periodic spin-up cycle) there is no
    initial state to leave out: the volumes count from zero at time zero.
    """
    t = np.asarray(times, dtype=float)
    cols = list(cum.columns)
    dis_cols = [c for c in cols if c.endswith('_OUT')
                and _prefix(c) in DISCHARGE_PREFIXES]
    rch_cols = [c for c in cols if c in RECHARGE_TERMS]
    # WP3/WP4: the streams and the ponds exchange with the aquifer directly
    # (MF6 has no unsaturated zone beneath them) -- left out, a losing
    # stream's seepage read as a storage deficit
    sin_cols = [c for c in cols if c.endswith('_IN')
                and _prefix(c) in SURFACE_PREFIXES]
    sout_cols = [c for c in cols if c.endswith('_OUT')
                 and _prefix(c) in SURFACE_PREFIXES]
    bin_cols = [c for c in cols if c.endswith('_IN') and _prefix(c) == 'GHB']
    i0 = 1 if (steady_first and len(t) > 1) else 0
    span = t[-1] - (t[0] if i0 else 0.0)
    if span <= 0 or area <= 0:
        raise ValueError('no elapsed time (%g d) or no area (%g m2)'
                         % (span, area))
    base = cum.iloc[0] if i0 else 0.0

    def vol(c):
        return float(cum[c].iloc[-1] - (base[c] if i0 else 0.0))

    def mmyr(v):
        return v / span / area * 1000.0 * 365.0

    v_rch = sum(vol(c) for c in rch_cols)
    v_dis = sum(vol(c) for c in dis_cols)
    step = np.zeros(len(t))
    for c in rch_cols:
        step += np.diff(np.concatenate([[0.0], cum[c].to_numpy(float)]))
    step = step[i0:]
    k = int(np.argmax(step)) if step.size else 0
    share = float(step[k] / v_rch) if v_rch > 0 else 0.0
    return {'recharge': mmyr(v_rch), 'discharge': mmyr(v_dis),
            'discharge_terms': {c: mmyr(vol(c)) for c in dis_cols},
            'surface_in': mmyr(sum(vol(c) for c in sin_cols)),
            'surface_out': mmyr(sum(vol(c) for c in sout_cols)),
            'boundary_in': mmyr(sum(vol(c) for c in bin_cols)),
            'surface_terms': {c: mmyr(vol(c)) for c in sin_cols + sout_cols},
            'area_km2': area / 1e6, 'days': span,
            'peak_share': share, 'peak_time': float(t[i0 + k]),
            'peak_mm': float(step[k]) / area * 1000.0 if step.size else 0.0}


def balance_lines(bal):
    """The balance as printed at the end of a run."""
    s_in = bal.get('surface_in', 0.0)
    s_out = bal.get('surface_out', 0.0)
    b_in = bal.get('boundary_in', 0.0)
    # the deficit is what storage made up: every way in against every way out
    deficit = (bal['discharge'] + s_out) - (bal['recharge'] + s_in + b_in)
    out = ['aquifer balance over %.0f d on %.3f km2: recharge to WT %.1f '
           'mm/yr  vs  discharge %.1f mm/yr  (deficit %.1f)'
           % (bal['days'], bal['area_km2'], bal['recharge'],
              bal['discharge'], deficit),
           '   discharge by term: ' + ', '.join(
               '%s %.1f' % (c, v) for c, v in sorted(
                   bal['discharge_terms'].items()))]
    if bal.get('surface_terms'):
        out.append('   streams and ponds: seepage into the aquifer %.1f, '
                   'groundwater into them %.1f mm/yr (net %+.1f to the aquifer)'
                   ' -- %s' % (s_in, s_out, s_in - s_out, ', '.join(
                       '%s %.1f' % (c, v) for c, v in sorted(
                           bal['surface_terms'].items()))))
    if b_in:
        out.append('   inflow through the GHB boundary %.1f mm/yr' % b_in)
    if bal['peak_share'] > 0.5:
        out.append('   WARNING: %.0f %% of that recharge arrived in ONE step '
                   '(t = %g d, %.0f mm over the catchment) -- a pulse, not a '
                   'rate; the whole-run figure says little about the rest '
                   'of the run.' % (100.0 * bal['peak_share'],
                                    bal['peak_time'], bal['peak_mm']))
    return out


def _outcrop_bottom(cMF):
    """The bottom of each cell's OUTCROP layer: the layer whose head MMsoil
    reads (the coupler's surface cell) -- below it the cell is dry. Not
    layer 1's: where layer 1 pinches out its bottom IS the soil base, and
    every water table in layer 2 read as a dry cell at the soil base -- full
    groundwater ET 1-6 m above the real table (2026-10-06)."""
    b = np.asarray(cMF.botm, dtype=float)
    k = np.clip(np.asarray(cMF.outcropL, dtype=int) - 1, 0, b.shape[0] - 1)
    return np.take_along_axis(b, k[None], axis=0)[0]


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


def _asc_on_grid(fn, cMF):
    """Read a Tier-A raster and put it on the model's grid (WP1c.5).

    Channel rasters need a resampling rule of their own. Averaging
    inputSTREAMw over ALL the source cells a mesh cell covers would smear a
    2 m channel across a 100 m cell -- every cell downstream of a channel would
    acquire a small positive width and become a stream cell. So the average is
    taken over the STREAM source cells only (`valid = arr > 0`), which gives a
    mesh cell the width of the channel actually crossing it, and leaves a cell
    that no channel crosses at zero.

    This whole raster path is transitional: WP3 builds the reach table from
    `inputSTREAM.csv`, which carries the channel as GRID-INDEPENDENT lines and
    needs no resampling at all.
    """
    arr = _asc(fn)
    proj = getattr(cMF, 'mesh_proj', None)
    if proj is None:
        return arr
    out = proj.sample2d(arr, fill=0.0, dtype=float, how='area',
                        valid=(np.asarray(arr) > 0))
    out = np.ma.filled(np.asarray(out), 0.0)
    # A channel can only exist where the model has a cell. A mesh cell that
    # overhangs the active domain can overlap a channel source cell while
    # being inactive itself, and MODFLOW then rejects the reach with "Cellid
    # is outside of the active model grid". Dropping those is the right
    # answer, but it TRUNCATES the network, so say how many went.
    act = np.asarray(cMF.outcropL) > 0
    dropped = int(np.count_nonzero((out > 0) & ~act))
    if dropped:
        print('   %s: %d channel cell(s) fall outside the active domain on '
              'this mesh and were dropped' % (os.path.basename(fn), dropped))
    return np.where(act, out, 0.0)


def _read_cell_grid(fn, cells):
    """Read an ESRI-ASCII grid and gather it back to a per-MM-cell vector.

    REFUSES A GRID THE CELLS CANNOT ADDRESS, naming both shapes. The cells
    belong to the grid THIS run uses; the file belongs to the grid that
    wrote it. Gather one through the other and numpy raises

        IndexError: index 65 is out of bounds for axis 0 with size 65

    from inside a list comprehension, which says nothing about saved state
    being tied to the grid that produced it. Worse, a file merely LARGER
    than the cells need would not raise at all: it would be silently
    gathered through the wrong cells and drive the steady period with
    another grid's recharge.
    """
    g = np.loadtxt(fn, skiprows=6)
    if g.ndim != 2:
        raise ValueError('%s is not a 2-D grid (shape %r)'
                         % (os.path.basename(fn), g.shape))
    rows = max(int(c[1]) for c in cells) + 1 if len(cells) else 0
    cols = max(int(c[2]) for c in cells) + 1 if len(cells) else 0
    if rows > g.shape[0] or cols > g.shape[1]:
        raise ValueError(
            '%s is %d x %d, but this run\'s cells reach (%d, %d): the file '
            'was written on a different grid.\nSaved state belongs to the '
            'grid that produced it. Regenerate it on this grid, or clear '
            'the [spinup] key that names it.'
            % (os.path.basename(fn), g.shape[0], g.shape[1], rows - 1,
               cols - 1))
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


# NO PARAMETER FILE (user, 2026-10-07). A run used to start by parsing
# example/<case>/MF_ws/__inputMF_flopy_v3_2s1L.ini and let the panels
# override it section by section. A new catchment has no such file, so the
# model description is built from the configuration and the dataset alone
# (clsMF.from_config + marmites_props); the file, where one still lies
# around, is not opened.


def setup_lamata(daily=True, nsp=None, grid='dis', nlay=None,
                 cfg=None, mesh_ws=None):
    """Replicate the driver setup; returns (cMF, mm, ctx, state, top, botm).

    ``cfg`` is REQUIRED: everything the model is comes from it and from the
    dataset (``DS``) -- the legacy MODFLOW parameter file is not read.
    ``nlay`` is accepted for the old call sites and must agree with
    ``cfg.layers.nlay``.

    When ``cfg.grid_kind`` is a genuine mesh producer -- 'quadtree',
    'voronoi' -- the model is built on the dataset's structured raster grid
    and then PROJECTED onto the mesh (``marmites_mesh.project_model``). The
    returned ``cMF`` then carries ``mesh_gridprops`` and ``mesh_proj``, and
    its arrays are ``(ncpl, 1)``; see ``marmites_mesh`` for why that shape
    makes the rest of the code work unchanged.
    """
    if cfg is None:
        raise SystemExit('setup_lamata needs the run configuration: the '
                         'model is built from it, not from a parameter file')
    if nlay is not None and int(nlay) != int(cfg.layers.nlay):
        raise SystemExit('nlay %d disagrees with layers.nlay = %d'
                         % (int(nlay), int(cfg.layers.nlay)))
    cUTIL = MMutils.clsUTILITIES(verbose=1)
    # THE GRID IS THE DATASET'S: the rectangle the converter wrote every
    # raster onto, from the Grid panel. The rasters carry it in their own
    # headers; a dataset with none has not been converted yet.
    _rect, _names, _others = props.dataset_grid(DS)
    if _rect is None:
        raise SystemExit('the dataset %s holds no raster to take the grid '
                         'from -- run the converter (Launch on the Grid '
                         'panel) first' % DS)
    cMF = ppMF.clsMF.from_config(
        cUTIL, MM_ws=DS, MM_ws_out=DS, MF_ws=os.path.join(DS, 'MF_ws'),
        grid=_rect, nlay=int(cfg.layers.nlay),
        hnoflo=float(cfg.layers.hnoflo),
        modelname=cfg.meta.model_name(cfg.paths.case),
        steady_recharge=float(cfg.spinup.steady_recharge))
    print('model: %s, %d layer(s), from the configuration and the dataset '
          '%s (no parameter file)' % (cMF.modelname, cMF.nlay, DS))
    # ONE SENTINEL, NOT TWO -- cMF's and the raster reader's. from_config
    # sets both from the panel already; this is the one function that keeps
    # them equal (setting only cMF's copy once turned La Mata's twelve
    # drains into 7800 at the aquifer floor), and it says nothing when they
    # agree.
    props.apply_hnoflo(cfg, cMF)
    # the raster grid, and any raster that does not sit on it
    props.check_grid(cMF, DS)
    # THE LAND SURFACE first: botm is elevation minus the cumulative
    # thickness. From the DEM panel 1 names, wrapped onto the dataset grid
    # (the parameter file named a 50 m elevation raster).
    props.land_surface(cfg, cMF, DS,
                       cache_dir=(os.path.join(mesh_ws, '_dem50')
                                  if mesh_ws else None))
    # THE LAYER PROPERTIES ARE THE PANEL'S -- ibound, thickness, k, k33, Ss,
    # Sy -- every one required: a blank one stops the run naming it.
    props.apply_layer_properties(cfg, cMF, DS, required=True)
    props.check_land_surface(cMF)
    props.initial_heads(cfg, cMF, DS)
    props.uzf_footprint(cMF)
    # THE CATCHMENT IS THE GEOGRAPHIC REFERENCE. It does not decide which
    # cells are active -- layers.ibound does, per layer -- but the active
    # cells have to sit inside it, and a model whose cells fall outside is
    # in a different coordinate system.
    props.check_catchment(cfg, cMF, mm_paths.GIS)
    # GHB and DRN from [ghb] and [drn] -- AFTER the properties, because a
    # drain taken at the base of its layer reads botm.
    props.apply_boundaries(cfg, cMF, DS)
    props.apply_uzf(cfg, cMF, DS)
    conv_fact = 1000.0

    # --- the forcing (WP1d) ------------------------------------------------
    # Replaces the positional parsing of __inputMMsurf4MMsoil.txt. That file
    # was written by MMsurf and read back here, and it was AUTHORITATIVE: with
    # MMsurf not running, editing Zr or kT* in the ini changed nothing, and
    # the two disagreed. Everything now comes from the configuration, and
    # kT_s arrives as the slope -- so the old `1.0 / x` is gone with the file.
    spec = _forcing(cfg)
    NMETEO, NVEG, NSOIL = spec.nmeteo, spec.nveg, spec.nsoil
    NCROP, NFIELD = spec.ncrop, spec.nfield
    # ABSOLUTE paths: with run.surface on, the series are in the workspace,
    # not the dataset. os.path.join(MM_ws, <absolute>) returns the absolute
    # one, so ppMFtime finds them wherever they are and cMF.MM_ws -- which
    # also locates the rasters and receives the stress-period files -- does
    # not have to move.
    inputDate_fn = spec.path('date')
    P_veg_fn, Pe_veg_fn = spec.path('rf_veg'), spec.path('tf_veg')
    PT_fn, LAI_fn = spec.path('pt_veg'), spec.path('lai_veg')
    PE_fn, Eo_fn = spec.path('pe'), spec.path('eo')
    P_irr_fn, Pe_irr_fn = spec.path('rf_irr'), spec.path('tf_irr')
    PT_irr_fn, crop_irr_fn = spec.path('pt_irr'), spec.path('crop_irr')
    Zr, kTg_min, kTg_max = spec.Zr, spec.kTg_min, spec.kTg_max
    kT_f, kT_s = spec.kT_f, spec.kT_s
    Zr_c = np.array(spec.Zr_c)
    kTg_min_c = np.array(spec.kTg_min_c)
    kTg_max_c = np.array(spec.kTg_max_c)
    kT_f_c = np.array(spec.kT_f_c)
    kT_s_c = np.array(spec.kT_s_c)

    # ppMFtime reads cMF.nper as the LONGEST a stress period may be, not as
    # a count: 1 gives one period per day, and anything larger lets it average
    # dry days together. Both numbers come from the front-end now -- the
    # MODFLOW ini's own value is no longer what decides it.
    if daily:
        cMF.nper = 1     # perlenmax=1 -> ppMFtime produces daily SPs (decision 4.3)
    elif getattr(cfg, 'run', None) is not None:
        cMF.nper = int(cfg.run.perlen_max)
    cMF.ppMFtime(inputDate_fn, P_veg_fn, Pe_veg_fn, PT_fn, LAI_fn, PE_fn, Eo_fn,
                 NMETEO, NVEG, NSOIL, P_irr_fn, Pe_irr_fn, PT_irr_fn, crop_irr_fn, NFIELD)
    print('time discretization: nper=%d over %d days (daily=%s)'
          % (cMF.nper, int(np.sum(cMF.perlen)), daily))

    # outcrop / masks
    cMF.outcropL = np.zeros((cMF.nrow, cMF.ncol), dtype=int)
    for L in range(cMF.nlay):
        ib = (np.abs(np.asarray(cMF.ibound))[L] != 0)
        cMF.outcropL += ((cMF.outcropL == 0) & ib) * (L + 1)

    gridMETEO = cMF.cPROCESS.inputEsriAscii(grid_fn='inputMETEOzones.asc', datatype=int)
    # THE SOIL ZONES AND THICKNESS FROM THE SOIL PANEL -- raster, polygon
    # layer or one value, whichever it names. They came from
    # inputSOILzones.asc and inputSOILthick.asc, filenames written right
    # here, so the panel's soil.zones and soil.thickness changed nothing.
    # On a MESH the polygon inputs are overlaid again on the mesh cells
    # below, and those are the ones the model uses: the 50 m values are an
    # intermediate, so their summary is not printed -- two lines with two
    # sets of numbers read as a contradiction (2026-10-07).
    _on_mesh = cfg.grid_kind in ('quadtree', 'voronoi')
    gridSOIL = props.soil_grid(cfg, cMF, DS, 'zones', kind='int')
    gridSOILthick = props.soil_grid(cfg, cMF, DS, 'thickness')
    print('soil zones: %s; soil thickness: %s -- from the panel'
          % (cfg.soil.zones.producer(), cfg.soil.thickness.producer()))
    gridIRR = cMF.cPROCESS.inputEsriAscii(grid_fn='inputIRRzones.asc', datatype=int)

    # THE VEGETATION COVER FROM THE SOIL PANEL: its vegetation layer, class
    # column and class table, put onto the grid by exact area overlay. It came
    # from inputVEG1area.asc .. inputVEG3area.asc, filenames written into
    # MARMITESprocess, so the panel changed nothing. The two agree on the
    # trees to 0.01 %; on grass they do not -- the polygons give ~89 % where
    # the old raster gave 25 % -- and the polygons are the reference: the
    # summer is carried by the grass wilting in the seasonal forcing, not by
    # the cover map.
    _veg = props.veg_cover(
        cfg, cMF, DS, NVEG,
        cache_dir=(os.path.join(os.path.dirname(mesh_ws), '_overlay')
                   if mesh_ws else None),
        verbose=not _on_mesh)
    (gridVEGarea, P_veg_zoneSP, Eo_zonesSP, PT_veg_zonesSP, Pe_veg_zonesSP, LAI_veg_zonesSP,
     PE_zonesSP, P_irr_zoneSP, Pe_irr_zoneSP, PT_irr_zonesSP, crop_irr_SP) = cMF.cPROCESS.inputSP(
        NMETEO=NMETEO, NVEG=NVEG, NSOIL=NSOIL, nper=cMF.nper,
        inputZON_SP_P_veg_fn=cMF.inputZON_SP_P_veg_fn, inputZON_SP_Pe_veg_fn=cMF.inputZON_SP_Pe_veg_fn,
        inputZON_SP_LAI_veg_fn=cMF.inputZON_SP_LAI_veg_fn, inputZON_SP_PT_fn=cMF.inputZON_SP_PT_fn,
        inputZON_SP_PE_fn=cMF.inputZON_SP_PE_fn, inputZON_SP_Eo_fn=cMF.inputZON_SP_Eo_fn,
        NFIELD=NFIELD, inputZON_SP_P_irr_fn=cMF.inputZON_SP_P_irr_fn,
        inputZON_SP_Pe_irr_fn=cMF.inputZON_SP_Pe_irr_fn, inputZON_SP_PT_irr_fn=cMF.inputZON_SP_PT_irr_fn,
        input_SP_crop_irr_fn=cMF.input_SP_crop_irr_fn, gridVEGarea=_veg)

    # THE SOIL COLUMN FROM THE SOIL PANEL. It came from MF_ws/inputSOILparam.txt,
    # at a path hard-coded right here, so the panel's soil.params was never
    # read and editing it changed nothing -- the audit behind the cookbook's
    # Appendix B found it. Same arrays, same order; only the source moved.
    if cfg is not None:
        _nsl, _nam, _st, _slprop, _Sm, _Sfc, _Sr, _S_ini, _Ks = \
            props.soil_parameters(cfg, nsoil=NSOIL)
        print('soil: %d zone(s), %s horizon(s), from the panel'
              % (len(_nsl), '/'.join(str(n) for n in _nsl)))
    else:
        _nsl, _nam, _st, _slprop, _Sm, _Sfc, _Sr, _S_ini, _Ks = \
            cMF.cPROCESS.inputSoilParam(
                SOILparam_fn=os.path.join('MF_ws', 'inputSOILparam.txt'),
                NSOIL=NSOIL)
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
    botm_l0 = _outcrop_bottom(cMF)
    for L in range(cMF.nlay):
        cMF.iuzfbnd[cMF.ibound[L] <= 0] = 0

    if nsp:
        cMF.nper = min(nsp, cMF.nper)
        cMF.perlen = np.asarray(cMF.perlen)[:cMF.nper]
        cMF.nstp = np.asarray(cMF.nstp)[:cMF.nper]

    # MF6 semantics: no hdry sentinel; MMsoil switches to h < botm dryness
    cMF.hdry = None

    # ---- WP1c.1: project onto an unstructured mesh, if one is configured.
    # Everything above ran on the structured raster, which is the point: the
    # mesh path reuses the whole validated setup and only changes the
    # discretisation the model is expressed on.
    grids = {'gridMETEO': gridMETEO, 'gridSOIL': gridSOIL,
             'gridSOILthick': gridSOILthick, 'gridIRR': gridIRR,
             'gridVEGarea': gridVEGarea}
    mesh_kind = (cfg.grid_kind if cfg is not None else 'structured')
    if mesh_kind in ('quadtree', 'voronoi'):
        import marmites_mesh
        import marmites_meshes
        gp, info = marmites_meshes.build_mesh(
            cfg, cMF, cache_dir=mesh_ws, dataset_dir=DS, model_ws=mesh_ws)
        print('mesh: %s, ncpl=%d, signature %s%s'
              % (info['kind'], info['ncpl'], info['signature'],
                 ' (from cache)' if info['cached'] else ''))
        if 'size_equiv' in info:
            print('      cell area %.4g..%.4g m2, mean %.4g (equivalent '
                  'square side %.1f m)'
                  % (info['area_min'], info['area_max'], info['area_mean'],
                     info['size_equiv']))
        cMF, grids, proj = marmites_mesh.project_model(
            cMF, gp, grids, how=cfg.grid.resample)
        _cov = proj.overlap_report()
        print('      resample=%s, %.1f source cell(s) per mesh cell, '
              'cell coverage %.3f, mesh tiles %.3f%% of the grid rectangle'
              % (cfg.grid.resample, _cov['src_per_cell_mean'],
                 _cov['coverage_mean'], 100.0 * _cov['domain_ratio']))
        cMF.mesh_gridprops, cMF.mesh_proj = gp, proj
        # which mesh: a saved state records it, and is refused on another
        cMF.mesh_signature = info['signature']
        gridMETEO = grids['gridMETEO']; gridSOIL = grids['gridSOIL']
        gridSOILthick = grids['gridSOILthick']
        gridIRR = grids['gridIRR']; gridVEGarea = grids['gridVEGarea']
        # THE POLYGON INPUTS ON THE MESH CELLS THEMSELVES. Above they were
        # overlaid on the 50 m grid and then resampled onto the mesh, so a
        # mesh cell smaller than 50 m -- half the La Mata Voronoi cells are
        # under 46 m2, the riparian corridor the 15,586 vegetation polygons
        # describe in detail -- inherited its 50 m cell's average cover and
        # majority soil zone. A RASTER source stays raster -> mesh: the
        # raster is the data.
        if cfg is not None:
            _cells = props.mesh_cells(cMF)
            _ovl = os.path.join(os.path.dirname(mesh_ws), '_overlay') \
                if mesh_ws else None
            gridVEGarea = props.veg_cover(cfg, cMF, DS, NVEG, cache_dir=_ovl,
                                          cells=_cells,
                                          cells_key=info['signature'])
            if cfg.soil.zones.producer() == 'layer':
                _z = np.asarray(props.soil_grid(cfg, cMF, DS, 'zones',
                                                kind='int', cells=_cells))
                _none = np.abs(_z - cMF.hnoflo) < 1.0      # no polygon reached
                _keep = _none & (np.abs(np.asarray(gridSOIL) - cMF.hnoflo)
                                 >= 1.0)
                gridSOIL = np.where(_keep, gridSOIL, _z).astype(int)
                print('soil zones: majority of the soil polygons on each mesh '
                      'cell%s' % ('' if not _keep.any() else
                                  ' (%d cell(s) no polygon reaches keep the '
                                  'resampled zone)' % int(_keep.sum())))
            if cfg.soil.thickness.producer() == 'layer':
                # top and every botm were set from the RESAMPLED thickness
                # (top = elev - thick): they move by the difference, so the
                # aquifer keeps its own thickness under the new soil column
                _t = np.asarray(props.soil_grid(cfg, cMF, DS, 'thickness',
                                                cells=_cells), dtype=float)
                _ok = (np.abs(_t - cMF.hnoflo) > 0.09) & \
                    (np.abs(np.asarray(gridSOILthick) - cMF.hnoflo) > 0.09)
                _d = np.where(_ok, _t - np.asarray(gridSOILthick, float), 0.0)
                cMF.top = np.asarray(cMF.top) - _d
                cMF.botm = np.asarray(cMF.botm) - _d[None, :, :]
                gridSOILthick = np.where(_ok, _t, gridSOILthick)
                print('soil thickness: area mean of the soil polygons on each '
                      'mesh cell (top and bottoms moved by up to %.2f m)'
                      % float(np.abs(_d).max()))
        # derived from the PROJECTED arrays, never carried over from the raster
        botm_l0 = _outcrop_bottom(cMF)
        print('projected onto the mesh: %d active cell(s) of %d'
              % (int(np.count_nonzero(cMF.outcropL > 0)), info['ncpl']))

    # ---- the land surface, from the raster panel 1 names (WP1d).
    # HERE, after the projection: the DEM is wrapped onto the grid the run
    # actually uses, so the survey is resampled once instead of twice.
    _applied, _note = _apply_dem(cMF, cfg, DS,
                                 cache_dir=os.path.join(mesh_ws or '', '_dem')
                                 if mesh_ws else None)
    print('elevation: %s' % _note)
    # which raster the land surface now is -- the CRR cascade reuses it
    # rather than wrapping the same file a second time (_crr_options)
    import marmites_dem as _mdem
    cMF.land_dem = _mdem.dem_path(DS) if _applied else None
    if _applied:
        botm_l0 = _outcrop_bottom(cMF)
    # THE BOUNDARY LINES, on the grid the run uses and its FINAL bottoms: a
    # drain at the base of its layer is placed after the land surface moved
    # the layers, not before, so no re-anchoring is needed
    if cfg is not None:
        props.apply_line_boundaries(cfg, cMF, DS)
    # how many DRN / GHB cells each layer HAS, now that they are built on
    # the run's grid -- the figures draw a layer's term only where it is > 0
    props.boundary_cell_counts(cMF)

    mm = MMsoil.clsMMsoil(hnoflo=cMF.hnoflo)
    cells = mm.build_cell_list(cMF)
    # Phase 4: cell geometry provider (DIS = legacy delr/delc; DISV = polygons)
    from marmites_grid import geometry_for
    if getattr(cMF, 'mesh_gridprops', None) is not None:
        # icell2d == the row index under the (ncpl, 1) convention
        from marmites_grid import VertexGeometry
        geom = VertexGeometry.from_vertices(
            cMF.mesh_gridprops['vertices'], cMF.mesh_gridprops['cell2d'],
            np.array([c[3] for c in cells], dtype=int), nlay=cMF.nlay)
    else:
        geom = geometry_for(cMF, cells, grid=grid)
    ctx = mm.build_context(cMF, cells, _nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _S_ini,
                           botm_l0, _Ks, gridSOIL, gridSOILthick, cMF.elev * 1000.0, gridMETEO,
                           INDEX_MM, INDEX_MM_SOIL,
                           P_veg_zoneSP, Eo_zonesSP, PT_veg_zonesSP, Pe_veg_zonesSP, PE_zonesSP,
                           gridVEGarea, LAI_veg_zonesSP, Zr, kTg_min, kTg_max, kT_f, kT_s, NVEG,
                           conv_fact, 1, P_irr_zoneSP, PT_irr_zonesSP, Pe_irr_zoneSP,
                           crop_irr_SP, gridIRR, Zr_c, kTg_min_c, kTg_max_c, kT_f_c, kT_s_c,
                           geom=geom)
    state = mm.init_state(ctx)
    top = np.asarray(cMF.top, dtype=float)
    botm = np.asarray(cMF.botm, dtype=float)
    return cMF, mm, ctx, state, top, botm, conv_fact


def _state_sidecar(a, prefix):
    """Path of the scope sidecar written beside a saved state prefix."""
    return mcfg.state_sidecar(a.state_dir, prefix)


def save_run_state(pref, b, cpl, st):
    """The rest of where a run ENDED, beside its heads: ``<pref>_state.npz``.

    Heads alone do not say where a run ended. The water the unsaturated
    zone holds (MF6's WCNEW per UZF object, as the next thti), the soil
    moisture (MMsoil's state) and the previous day's seepage, rejected
    infiltration and UZF ET (what the first day reads before MF6 has solved
    anything) are the rest -- what a periodic spin-up cycle carries, and
    what a run started from these heads needs to start where this one
    stopped. Covered by the same scope sidecar as the heads.
    """
    shape = (int(b.nlay), int(b.nrow), int(b.ncol))
    wc = (b.thti_from_wc(cpl.uzf_wc_final)
          if getattr(cpl, 'uzf_wc_final', None) is not None
          else np.full(shape, np.nan))
    carry = getattr(cpl, 'carry_out', None) or {}
    n = int(np.asarray(st.Ssoil_ini).shape[0])
    fn = pref + '_state.npz'
    lak = getattr(cpl, 'lak_stage_final', None)
    np.savez(fn, uzf_wc=wc, soil=np.asarray(st.Ssoil_ini, dtype=float),
             exf=np.asarray(carry.get('exf', np.zeros(n)), dtype=float),
             rej=np.asarray(carry.get('rej', np.zeros(n)), dtype=float),
             etuzf=np.asarray(carry.get('etuzf', np.zeros(n)), dtype=float),
             # each pond's stage; empty when the run had no LAK
             lak_stage=(np.zeros(0) if lak is None
                        else np.asarray(lak, dtype=float)))
    return fn


LAST_CYCLE = '_lastcycle'


def save_cycle_state(a, b, cMF, ctx, base, heads, cpl, st, res, grid=None):
    """Save a spin-up cycle that PASSED its check as ``<base>_lastcycle``:
    heads, the rest of the state, and its mean recharge / ETg -- everything a
    run started from it needs. Overwritten by every cycle that passes.

    The final state is only written once the spin-up ends, so a later cycle
    that failed its check (2026-10-04: one sub-step of SP74 in cycle 2)
    threw away a cycle that had converged, ~2.75 h of run. Written AS the
    cycle passes rather than when a later one fails, so it also survives a
    killed process or a crashed machine. Returns the prefix name.
    """
    name = base + LAST_CYCLE
    pref = _state_out(a, name)
    b.save_heads_asc(heads, pref)
    save_run_state(pref, b, cpl, st)
    _write_cell_grid(res['perc'].mean(axis=0), ctx.cells, cMF.nrow, cMF.ncol,
                     cMF, pref + '_perc.asc')
    _write_cell_grid(res['etg'].mean(axis=0), ctx.cells, cMF.nrow, cMF.ncol,
                     cMF, pref + '_etg.asc')
    _write_state_scope(a, a.config, name, grid=grid)
    return name


def drop_cycle_state(a, base, nlay):
    """Remove ``<base>_lastcycle`` once the spin-up's own state is saved: it
    is then the same cycle under a second name."""
    pref = _state_out(a, base + LAST_CYCLE)
    for fn in (['%s_l%d.asc' % (pref, k + 1) for k in range(int(nlay))]
               + [pref + s for s in ('_state.npz', '_perc.asc', '_etg.asc')]
               + [_state_sidecar(a, base + LAST_CYCLE)]):
        if os.path.exists(fn):
            os.remove(fn)


def load_run_state(pref, b, ctx):
    """What :func:`save_run_state` wrote for ``pref``, checked against this
    model -- or None, saying why, when there is none or it does not fit
    (a state saved before 2026-09-25 has heads only)."""
    fn = pref + '_state.npz'
    if not os.path.exists(fn):
        print('   no %s beside the heads: the unsaturated zone and the soil '
              'start from the panel\'s initial values'
              % os.path.basename(fn))
        return None
    z = np.load(fn)
    shape = (int(b.nlay), int(b.nrow), int(b.ncol))
    soil = np.asarray(z['soil'], dtype=float)
    want = (int(ctx.ncell), int(ctx._nslmax))
    bad = [what for what, ok in (
        ('unsaturated zone %s, model %s' % (z['uzf_wc'].shape, shape),
         z['uzf_wc'].shape == shape),
        ('soil %s, model %s' % (soil.shape, want), soil.shape == want),
        ('exchange terms %d cells, model %d' % (z['exf'].size, want[0]),
         z['exf'].size == want[0]))
        if not ok]
    if bad:
        print('   %s does not fit this model (%s): the unsaturated zone and '
              'the soil start from the panel\'s initial values'
              % (os.path.basename(fn), '; '.join(bad)))
        return None
    print('   the unsaturated zone, the soil and the previous day\'s exchange '
          'terms from %s' % os.path.basename(fn))
    # the ponds' stages, when the state was saved with LAK on (the builder
    # checks the count against the ponds it builds)
    lak = (np.asarray(z['lak_stage'], dtype=float)
           if 'lak_stage' in z.files else np.zeros(0))
    return {'uzf_wc': np.asarray(z['uzf_wc'], dtype=float), 'soil': soil,
            'carry': {k: np.asarray(z[k], dtype=float)
                      for k in ('exf', 'rej', 'etuzf')},
            'lak_stage': lak if lak.size else None}


def _run_grid(cMF):
    """The grid a saved state has to match: its shape as the state files
    are written -- ``(ncpl, 1)`` on a mesh -- and the mesh signature."""
    return {'shape': (int(cMF.nrow), int(cMF.ncol)),
            'signature': getattr(cMF, 'mesh_signature', None)}


def _write_state_scope(a, cfg, prefix, grid=None):
    """Record WHICH configuration produced a saved state (WP0.6) -- and,
    since 2026-09-27, which GRID: the scope names only the grid kind, so a
    state saved on one voronoi mesh passed as belonging to another."""
    import json
    payload = {'state_hash': cfg.state_hash(), 'scope': cfg.state_scope()}
    if grid:
        payload['grid'] = {'shape': list(grid['shape']),
                           'signature': grid.get('signature')}
    with open(_state_sidecar(a, prefix), 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, default=str)


def _check_state_scope(a, cfg):
    """Report on saved state that cannot be reused here. Never refuses.

    Saved state belongs to the grid and layer set that produced it: handing
    a structured-grid field to a mesh model would give MODFLOW an array of
    the wrong length. This USED TO REFUSE, with "CONFIG ERROR: ... has no
    scope sidecar". That was the wrong answer to "this state does not fit",
    because a run can always start from the land surface -- which is where
    a spin-up starts from anyway. props.resolve_initial_heads decides, and
    the run says which it chose.

    `spinup.steady_means` is DROPPED when it cannot be used here, and the
    run says so. It used only to be reported -- and then read anyway, which
    is how a note about state belonging to another grid was followed three
    lines later by

        IndexError: index 65 is out of bounds for axis 0 with size 65

    from gathering a 65x60 structured .asc through voronoi cell indices. A
    warning that does not prevent the failure it predicts is not a warning.
    Nothing is lost by dropping it: with no means pinned, the spin-up takes
    the steady period's averages from the cycle it just ran, which is where
    they came from in the first place.
    """
    why = mcfg.state_problem(cfg, a.state_dir)
    if why:
        print('note: %s' % why)
    for what, prefix in (('spinup.strt_heads', cfg.spinup.strt_heads),
                         ('spinup.steady_means', cfg.spinup.steady_means)):
        prefix = (prefix or '').strip()
        if prefix and not os.path.exists(_state_sidecar(a, prefix)):
            print('note: %s = %r has no scope sidecar (written before WP0)'
                  % (what, prefix))
    # strt_heads has a fallback of its own -- resolve_initial_heads starts the
    # water table from the land surface and says so -- so only the means are
    # dropped here.
    if why and a.steady_means:
        print('   %s is NOT used: the steady period takes its averages from '
              'the spin-up cycle instead.' % 'spinup.steady_means')
        a.steady_means = None


def _args_from_config(cfg, probe=False):
    """Map a RunConfig onto the legacy attribute names the driver body uses.

    WP0 is a REFACTOR, not a rewrite: the body of main() below is untouched, so
    a configuration that mirrors the old flags produces byte-identical MF6
    input. This function is the whole of the translation, and it is the only
    place the old flag vocabulary survives.

    Blank strings and 0 in the schema mean "not set" for the flags whose
    argparse default was ``None``.
    """
    def _or_none(s):
        s = (s or '').strip()
        return s or None

    lak_source = None
    if cfg.lak.enable:
        # The builder fits an embedded lake to each pond FOOTPRINT, so it
        # needs the polygons -- [lak] geometry, the derived GeoJSON the
        # converter writes -- not [lak] source, which is the centroid/DEM
        # table. Passing source here is what made LAK fail on a missing .dbf.
        lak_source = cfg.lak.geometry or 'inputPONDS.geojson'

    libmf6 = resolve_libmf6((cfg.paths.libmf6 or '').strip())

    return argparse.Namespace(
        # run
        nsp=(cfg.run.nsp or None), daily=cfg.run.daily, ats=cfg.run.ats,
        ats_dtmin=float(cfg.run.ats_dtmin),
        build_only=cfg.run.build_only, standalone=_or_none(cfg.run.standalone),
        probe=bool(probe or cfg.run.probe),
        max_discrepancy=cfg.run.max_discrepancy,
        allow_bad_budget=cfg.run.allow_bad_budget,
        # grid / layers. Two separate notions (WP1c.1): `grid` is the MF6
        # DISCRETISATION ('dis' or 'disv', all clsMF6 understands), `mesh_kind`
        # is the PRODUCER that made it. Conflating them is what made 'quadtree'
        # reach clsMF6 as a grid type it rejects.
        grid=('dis' if cfg.grid_kind == 'structured' else 'disv'),
        mesh_kind=cfg.grid_kind,
        nlay=cfg.layers.nlay,
        # packages
        uzf_vks_scale=cfg.uzf.vks_scale,
        uzf_et_form=cfg.et.unsat_form,
        seep=cfg.seep.kind, seep_cond=cfg.seep.cond,
        seep_ddrn=float(cfg.seep.ddrn),
        seep_base=float(cfg.seep.base), seep_cond_from=str(cfg.seep.cond_from),
        sfr=cfg.sfr.enable,
        sfr_rhk=(cfg.sfr.rhk.value if cfg.sfr.rhk.value is not None else 0.1),
        lak=lak_source, lak_bedleak=cfg.lak.bedleak,
        # sprinkler irrigation down the whole soil column (2026-10-04)
        irr_column=(cfg.soil.irr_infiltration == 'column'),
        # WP5, the runoff cascade
        crr=bool(cfg.crr.enable), crr_beta=float(cfg.crr.beta),
        crr_sinks=str(cfg.crr.sinks), crr_dem=_or_none(cfg.crr.dem),
        # spin-up / initial state
        spinup=cfg.spinup.cycles, spinup_tol=cfg.spinup.tol,
        strt_heads=_or_none(cfg.spinup.strt_heads),
        steady_means=_or_none(cfg.spinup.steady_means),
        save_strt=_or_none(cfg.spinup.save_strt),
        save_means=_or_none(cfg.spinup.save_means),
        strt_dem=(list(cfg.spinup.strt_dem) or None),
        # post-processing
        # BOTH: the panel offers run.plot as the group's switch and
        # postproc.enable as a field, and each promises to stop the figures.
        # run.plot was read by nothing, so only one of the two kept its word.
        postproc=bool(cfg.postproc.enable and cfg.run.plot),
        # the input stage: [postproc] input_maps, the one switch for it
        preproc=cfg.postproc.input_maps,
        postproc_only=cfg.postproc.only, gis_ws=_or_none(cfg.paths.gis_ws),
        sankey_min_flux=cfg.postproc.sankey_min_flux,
        sankey_full=cfg.postproc.sankey_full,
        sankey_obs_years=cfg.postproc.sankey_obs_years,
        map_days=cfg.postproc.map_days,
        # where things go
        ws=_or_none(cfg.paths.ws),
        ws_root=str(mm_paths.WS_ROOT),
        run_tag=(cfg.meta.name or None),
        libmf6=(libmf6 or None),
        # carried for provenance
        config=cfg,
    )


def main():
    ap = argparse.ArgumentParser(
        description='Run the MARMITES / MODFLOW 6 coupled model for one case '
                    'study. Everything that used to be a flag is now a key in '
                    'the configuration file; see code/configs/lamata.toml.',
        epilog='example:  python code/tests/run_lamata_mf6.py '
               '--config code/configs/lamata.toml --set run.nsp=365')
    ap.add_argument('--config', required=True, metavar='FILE',
                    help='run configuration (TOML). The single source of truth '
                         'for every model and run setting.')
    ap.add_argument('--set', action='append', dest='overrides', default=[],
                    metavar='SECTION.KEY=VALUE',
                    help='one-off override of a configuration key, repeatable. '
                         'Every override is echoed at startup, because a '
                         'setting that changes a run without appearing '
                         'anywhere is how a multi-hour run gets wasted.')
    ap.add_argument('--run-tag', default=None, metavar='TAG',
                    help='label this run, overriding meta.name -> '
                         '<ws-root>/out_<YYYYMMDDHHMM>_<TAG>')
    ap.add_argument('--probe', action='store_true',
                    help='list the MF6 memory variables and exit (diagnostic)')
    ns = ap.parse_args()

    cfg = mcfg.load_run_config(ns.config)
    cfg.apply_overrides(ns.overrides)
    if ns.run_tag:
        cfg.meta.name = ns.run_tag
    cfg.require_implemented_grid()
    print('config: %s  (hash %s)' % (os.path.abspath(ns.config), cfg.config_hash()))
    # a renamed or RETIRED key in the file: the panels say so, and so does
    # the run -- a file still asking for the WEL route or the iterative
    # coupling (2026-10-07) runs without them, and its log must tell
    for _line in getattr(cfg, 'migrated', ()):
        # ASCII: the run's console or log may not be UTF-8
        print('NOTE (config): %s' % _line.replace('—', '--'))
    # THE DATASET IS THE CASE'S: <example_root>/<paths.case>, panel 0's
    # folder (mm_paths.dataset_dir). It was this file's own location,
    # <repo>/example/LaMata, whatever paths.case or panel 0 said -- one more
    # thing a second catchment could not change (2026-10-07).
    global DS
    DS = str(mm_paths.dataset_dir(cfg.paths.case))

    a = _args_from_config(cfg, probe=ns.probe)
    if a.postproc_only:
        a.postproc = True
    if a.ws is None:
        # One workspace per MESH, not per discretisation: a quadtree and a
        # DISV-from-DIS model are both 'disv' but share no file, and letting
        # them overwrite each other is a grid-cache bug waiting to happen.
        a.ws = mcfg.state_workspace(cfg, a.ws_root)
    os.makedirs(a.ws, exist_ok=True)
    # results folder for this run, in the legacy out_<stamp>_<tag> style
    tag = a.run_tag or ('%dlay' % (a.nlay or 6))
    a.out_dir = os.path.join(a.ws_root,
                             'out_%s_%s' % (time.strftime('%Y%m%d%H%M'), tag))
    # PROVENANCE (WP0.5): the RESOLVED configuration -- after --set -- is copied
    # into the run folder, so every result says exactly what produced it.
    _prov_dir = os.path.join(a.out_dir, '_input')
    os.makedirs(_prov_dir, exist_ok=True)
    cfg.write_toml(os.path.join(_prov_dir, 'resolved_config.toml'))
    # Saved run state (equilibrated heads, steady means) is run OUTPUT, so it is
    # written to the workspace; reading falls back to the baseline committed in
    # the repo's example/LaMata/MF_ws so `spinup.strt_heads` keeps working.
    a.state_dir = a.ws
    # THE LIBRARY IS CHECKED BEFORE THE BUILD. It is not used until the model
    # has been written, which on La Mata is minutes away, and a run that
    # spends them only to exit on a mistyped path has wasted all of them.
    # Blank is not a mistake -- it means "build and stop" -- so only a path
    # that was GIVEN is judged.
    try:
        check_libmf6(a.libmf6)
    except LibMF6Error as exc:
        sys.exit('paths.libmf6: %s' % exc)
    # run.model OFF: produce the forcing and stop. MMsurf is a run of its
    # own -- that is what the switch on the driving-forces panel promises --
    # and everything below this line is the model. The switch was read by
    # NOTHING before WP1d: turning it off changed the file and not the run.
    if not cfg.run.model:
        _forcing(cfg)
        print('run.model is off: the forcing is done and the model is not '
              'built. Turn it on to run MMsoil + MODFLOW 6.')
        return

    # STATE GUARD (WP0.6): saved state is only valid for the grid and layer set
    # it was produced on, so the sidecar carries THAT scope rather than the whole
    # configuration -- a full-config hash would trip on an unrelated key such as
    # postproc.map_days. A mismatch stops the run naming the offending key,
    # which is the failure the CdL grid-design cache produced once.
    #
    # AFTER the switch, not before: this is about MODFLOW's initial state, so
    # it has no business stopping a forcing-only run -- which is exactly what
    # it did, and the reason MMsurf never started.
    _check_state_scope(a, cfg)

    cMF, mm, ctx, state, top, botm, conv_fact = setup_lamata(
        daily=a.daily, nsp=a.nsp, grid=a.grid, nlay=a.nlay,
        cfg=cfg, mesh_ws=os.path.join(a.ws, '_mesh'))
    # how the irrigated fields' irrigation enters the soil (MMsoil reads it)
    ctx.irr_column = bool(getattr(a, 'irr_column', False))
    _irr = getattr(ctx, 'gridIRR', None)
    _n_irr = (0 if _irr is None else
              sum(1 for c in ctx.cells if int(_irr[c[1], c[2]]) > 0))
    print('irrigation: %d irrigated cell(s); the irrigation enters %s -- from '
          'the panel' % (_n_irr, 'the whole soil column, top-down (sprinklers)'
                         if ctx.irr_column else
                         'the top horizon only, like rain'))
    if a.postproc_only:
        # Re-draw the figures from a run that already happened: everything the
        # post-processing reads is on disk (the coupled HDF5 for the MM side,
        # the .hds/.cbc/.grb for the aquifer side), so MODFLOW need not run
        # again. Iterating on a figure costs seconds instead of the full run.
        h5_fn = os.path.join(a.ws, MF6Coupler.RESULTS_H5)
        if not os.path.exists(h5_fn):
            raise SystemExit('--postproc-only needs a previous run: %s not found'
                             % h5_fn)
        with h5py.File(h5_fn, 'r') as f:
            res = {k: f[k][:] for k in f.keys()}
        print('re-using %s (%d stress period(s))'
              % (h5_fn, res['wb_ts'].shape[0]))
        _run_postproc(a, cMF, ctx, res)
        return
    _gp = getattr(cMF, 'mesh_gridprops', None)
    _grid = _run_grid(cMF)
    b = clsMF6(cMF, top=top, botm=botm, sim_ws=a.ws, daily=True, grid=a.grid,
               vertices=(_gp['vertices'] if _gp else None),
               cell2d=(_gp['cell2d'] if _gp else None),
               strt_from_dem=(tuple(a.strt_dem) if a.strt_dem else None))
    b.seep = a.seep
    b.ats = a.ats
    b.ats_dtmin = float(getattr(a, 'ats_dtmin', b.ats_dtmin))
    if b.ats:
        print('ATS: a step MF6 cannot solve is retried at 1/5 of its length, '
              'down to %g d -- from the panel' % b.ats_dtmin)
    b.drn_seep_cond = float(a.seep_cond)
    b.drn_seep_ddrn = float(getattr(a, 'seep_ddrn', b.drn_seep_ddrn))
    b.drn_seep_base = float(getattr(a, 'seep_base', 0.0))
    b.drn_seep_cond_from = str(getattr(a, 'seep_cond_from', 'value'))
    b.uzf_vks_scale = float(a.uzf_vks_scale)
    if cfg is not None:
        b.uzf_thtr_from = str(cfg.uzf.thtr_from)
        b.solver = dataclasses.asdict(cfg.solver)
        b.outer_maximum = int(cfg.solver.outer_maximum)
        print('solver: %s, outer_dvclose %g m (max %d), inner_dvclose %g m, '
              'inner_rclose %g m3/d, cell averaging %s -- from the panel'
              % (cfg.solver.complexity.upper(), cfg.solver.outer_dvclose,
                 cfg.solver.outer_maximum, cfg.solver.inner_dvclose,
                 cfg.solver.inner_rclose, cfg.solver.cell_averaging))
    # UNSATURATED-ZONE ET. The extinction depth follows the usual rule --
    # a raster, a column of the vegetation layer, or one value -- resolved
    # on the grid the run uses; or, per vegetation zone, the rooting depth
    # of each cell's cover (et.extdp_from, WP2 2.2).
    b.uzf_et_form = str(getattr(a, 'uzf_et_form', 'etwc'))
    # groundwater ET: the builder makes the two EVT packages, MMsoil hands
    # over what Eg and Tg could take (the coupler sets ctx.gw_evt), MF6
    # takes it at the head it solves for. The WEL route went on 2026-10-07.
    if cfg is not None:
        b.evt_nseg = int(cfg.et.evt_nseg)
        b.evt_ramp = float(cfg.et.evt_ramp)
        print('groundwater ET: EVT (Eg, Tg) at the head MF6 solves for, %d '
              'segments, Tg ramp %g m above the root tips -- from the panel'
              % (b.evt_nseg, b.evt_ramp))
    if cfg is not None:
        b.uzf_extdp = extdp_grid(cfg, cMF, ctx, DS)
    if a.sfr:
        # WP1d: the network is the hydrography the modeller MAPPED, burned onto
        # whichever grid panel 1 produced -- not inputSTREAMw.asc, which was
        # the alluvium footprint of Soil_type.shp with one width for the whole
        # catchment. Width and incision are resolved after routing.
        import marmites_channel as mch
        import marmites_vector as mv
        lines, seg_params = mch.read_stream_lines(
            os.path.join(DS, 'inputSTREAM.csv'),
            os.path.join(DS, 'inputSTREAM_param.csv'))
        vgrid = mv.TargetGrid.from_cMF(cMF)
        present, seg_of_cell, ch_len = mch.burn_channel(
            lines, vgrid, (cMF.nrow, cMF.ncol))
        act = np.asarray(cMF.outcropL) > 0
        dropped = int(np.count_nonzero((present > 0) & ~act))
        if dropped:
            print('   %d channel cell(s) fall outside the active domain and '
                  'were dropped' % dropped)
        b.sfr_pondw = np.where(act, present, 0.0)
        b.sfr_pondhmax = np.zeros_like(b.sfr_pondw)
        b.sfr_seg_of_cell = seg_of_cell
        b.sfr_seg_params = seg_params
        b.sfr_cell_length = ch_len
        b.sfr_width_source = cfg.sfr.width if cfg else None
        b.sfr_depth_source = cfg.sfr.depth if cfg else None
        b.cell_area = float(np.mean(np.asarray(ctx_geom_area(cMF))))
        print('   stream network: %d segment(s) -> %d cell(s), %.0f m mapped'
              % (len(lines), int((b.sfr_pondw > 0).sum()), float(ch_len.sum())))
        b.sfr_rhk = float(a.sfr_rhk)
        # WP3: the panel's bed rules, and its values where the converter's
        # per-segment table has none
        if cfg is not None:
            b.sfr_min_slope = float(cfg.sfr.min_slope)
            b.sfr_monotonic = bool(cfg.sfr.monotonic_bed)
            if cfg.sfr.manning.value is not None:
                b.sfr_man = float(cfg.sfr.manning.value)
            if cfg.sfr.rbth.value is not None:
                b.sfr_rbth = float(cfg.sfr.rbth.value)
            # A single panel VALUE is the value on every segment, whatever
            # the converted table still holds: Launch re-converts when a
            # shapefile changes, not when a panel number does
            b.sfr_param_fixed = {
                key: float(src.value)
                for key, src in (('manning', cfg.sfr.manning),
                                 ('rhk', cfg.sfr.rhk), ('rbth', cfg.sfr.rbth))
                if src.producer() == 'value'}
    if a.lak:
        shp = a.lak if os.path.isabs(a.lak) else os.path.join(DS, a.lak)
        b.lak_shapefile = shp
        b.lak_bedleak = float(a.lak_bedleak)
        # WP1d: the pond depth used to be read from inputSTREAMhmax.asc, which
        # is retired -- and was never a pond map anyway (1.0 m on the alluvium
        # polygons). It comes from [lak] depth now; None falls back to
        # marmites_lak.POND_DEPTH.
        _d = (cfg.lak.depth if cfg else None)
        b.lak_depth = None if _d is None else float(_d)
        # the settings that stabilised CdL's perched ponds (cookbook 4.1):
        # they were on the panel and reached nothing -- MF6 ran its own
        # defaults, 100 iterations and a 1e-5 m stage change
        if cfg is not None:
            b.lak_surfdep = float(cfg.lak.surfdep)
            b.lak_maxiter = int(cfg.lak.maxiter)
            b.lak_stagechg = float(cfg.lak.stagechg)
    # WHERE THE RUN STARTS. The saved state when it exists and belongs to
    # this grid and layer set; otherwise a cold start from layers.strt, or
    # from the land surface when that is blank. Never an array nobody chose.
    _kind, _payload, _why = props.resolve_initial_heads(cfg, a.state_dir,
                                                        grid=_grid)
    saved_state = None
    if _kind == 'saved':
        pref = _state_in(a, str(_payload), '_l1.asc')
        b.strt_array = b.load_heads_asc(pref)
        print('initial heads loaded from %s_l*.asc' % pref)
        # A RUN FROM SAVED HEADS STARTS FROM THEM: no steady period, which
        # MF6 solves ignoring the initial heads -- the saved heads were only
        # ever a first guess for it. The rest of the saved state, where
        # there is one, starts the unsaturated zone and the soil.
        b.steady_first = False
        saved_state = load_run_state(pref, b, ctx)
        if saved_state is not None:
            b.uzf_thti_carry = saved_state['uzf_wc']
            b.lak_strt_carry = saved_state.get('lak_stage')
            if a.lak and b.lak_strt_carry is None:
                print('   the state was saved without LAK: each lake starts '
                      'from the water table under it')
        if a.steady_means:
            print('   spinup.steady_means is not used: a run from saved heads '
                  'has no steady period to drive')
    elif _kind == 'strt':
        # layers.strt, on the grid the run uses (projected with the rest),
        # kept above each cell's bottom like any explicit start
        b.strt_array = np.asarray(cMF.strt, dtype=float)
    else:
        b.strt_from_dem = tuple(_payload)
    b.build()
    b.write()
    print('MF6 (%s) simulation written to %s  (%d SPs%s, %d UZF cells, %d EVT '
          'records per package)'
          % (a.grid.upper(), a.ws, b.nper,
             ' incl. steady' if b.steady_first else ', no steady period',
             b.nuzfcells, b.ncell))
    if b.seep == 'drn':
        _rng = b.drn_seep_cond_range or (b.drn_seep_cond, b.drn_seep_cond)
        print('   seepage: DRN_SEEP, %d drains %s, cond %s, ramped in over '
              '%.4g m above it -- from the panel (UZF SIMULATE_GWSEEP off)'
              % (b.ndrnseep,
                 ('%.4g m below the soil base' % b.drn_seep_base)
                 if b.drn_seep_base > 0 else 'at the soil base',
                 ('%.4g m2/d' % _rng[0]) if b.drn_seep_cond_from != 'uzf'
                 else ('UZF\'s area x vks / SURFDEP, %.4g..%.4g m2/d'
                       % _rng),
                 b.drn_seep_ddrn))
        if b.drn_seep_lifted:
            print('   seepage: %d drain(s) lifted to 1 cm above their cell '
                  'bottom (MF6 refuses a drain below it)' % b.drn_seep_lifted)
        if getattr(b, 'drn_seep_lake_skipped', 0):
            print('   seepage: none in %d pond cell(s), as in the stream '
                  'cells -- the pond and the aquifer exchange through its bed '
                  '(lak.bedleak)' % b.drn_seep_lake_skipped)
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


    # WP5: the runoff cascade, settled once for every spin-up cycle
    crr_opts = _crr_options(a, cMF, DS, cache_dir=os.path.join(a.ws, '_mesh',
                                                               '_dem'))

    if a.build_only or not a.libmf6:
        if crr_opts is not None:
            # the cascade is built from the model as written, so a build-only
            # check sees its sinks before anything runs (cookbook 5.3)
            MF6Coupler(mm, ctx, mm.init_state(ctx), b, conv_fact=conv_fact,
                       crr=crr_opts)
        # preproc only needs the built model, so it can run without libmf6;
        # postproc needs MF6 output, so it is skipped here
        if a.preproc:
            from marmites_postprocess import run_preproc
            run_preproc(a.ws, DS, name=cMF.modelname.lower(),
                        gis_ws=a.gis_ws)
        print('\nNo --libmf6 given: stopping after build.\n'
              'To run coupled: set paths.libmf6 (the Run panel) and run.build_only '
              'off.')
        return

    # --- robust libmf6 handling ---------------------------------------
    # WinError 87 from CDLL(winmode=0x08) happens when the path is not
    # fully qualified (relative, or polluted with quotes by IDE run
    # configs); dependency DLLs also need the bin dir on the search path.
    # The path itself was settled and checked BEFORE the build -- see
    # check_libmf6 -- so by here it is a file that exists.
    lib = os.path.abspath(a.libmf6.strip('"').strip("'"))
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
        for key in ('/UZF', '/EVT', '/DIS', '/X'):
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
    if a.steady_means and cfg is not None:
        # the same question as for the heads, now that the grid is known:
        # means saved on another mesh are dropped, not gathered through
        # this one's cell indices
        _why_m = mcfg.state_problem(cfg, a.state_dir, grid=_grid,
                                    keys=('spinup.steady_means',))
        if _why_m:
            print('%s\n   spinup.steady_means is NOT used: the steady period '
                  'takes a uniform recharge instead.' % _why_m)
            a.steady_means = None
    if a.steady_means:
        mp = _state_in(a, a.steady_means, '_perc.asc')
        steady_perc = _read_cell_grid(mp + '_perc.asc', ctx.cells)
        steady_etg = _read_cell_grid(mp + '_etg.asc', ctx.cells)
        print('steady-state means loaded from %s_{perc,etg}.asc' % mp)

    # observation cells at which to keep the full MM flux series (per-point
    # Sankey / time series). Resolved once; the coupler captures them each cycle.
    obs_idx, obs_names = [], []
    # THE OBSERVATION FILES FROM THE STATE VARIABLES PANEL, before anything
    # reads one. The post-processing named them as literals, so the panel's
    # [obs] changed nothing -- see marmites_postprocess.OBS. Set in this
    # process, which is also the one the post-processing runs in.
    if getattr(a, 'config', None) is not None:
        import marmites_postprocess as _pp
        _o = _pp.use_observations(a.config)
        print('observations: points %s, heads %s_*, soil moisture %s_*, '
              'runoff %s_* -- from the panel'
              % (_o['table'], _o['heads'], _o['sm'], _o['ro']))
    # ON EVERY RUN, figures or not: the observation exports (WP6.4) are what
    # a calibration's forward runs are read by, and those draw nothing
    try:
        from marmites_postprocess import resolve_obs_cells
        obs_idx, obs_names = resolve_obs_cells(cMF, ctx, DS,
                                               verbose=bool(a.postproc))
    except Exception as exc:
        print('obs-cell resolution skipped: %r' % exc)

    ncyc = max(1, int(a.spinup))
    prev = None
    # Bound at the end of every cycle below. Named here because the cyc > 0
    # branch reads it: the loop always runs at least once, so it is never
    # actually unbound there, but nothing in the block says so and pyflakes
    # reports it as an undefined name.
    prev_heads = None
    spin_converged, delta = False, float('nan')
    # where a cycle that passed its check is kept while later ones run
    cycle_base = (a.save_strt or 'hi_spinup') if ncyc > 1 else None
    last_good = None                     # (cycle number, saved name)
    for cyc in range(ncyc):
        if cyc > 0:
            # PERIODIC SPIN-UP: the next cycle starts where this one ended --
            # its heads, the water its unsaturated zone holds and its soil --
            # and runs WITHOUT a steady period, which ignores the heads and
            # restarts from a mean-forcing equilibrium the year never
            # reaches (2026-09-25: every cycle ended 0.6 m below its start,
            # so the cycles agreed with each other and not with themselves).
            b.strt_array = prev_heads
            b.steady_first = False
            b.uzf_thti_carry = (None if cpl.uzf_wc_final is None
                                else b.thti_from_wc(cpl.uzf_wc_final))
            # ... and each pond where it ended, not re-derived from the WT
            b.lak_strt_carry = getattr(cpl, 'lak_stage_final', None)
            b.build()
            b.write()
            st.carried = True              # the soil the last cycle ended with
        else:
            st = mm.init_state(ctx)
            if saved_state is not None:
                # a run from a saved state: its soil, not the panel's
                st.Ssoil_ini[:] = saved_state['soil']
                st.carried = True
        _carry = (cpl.carry_out if cyc > 0 else
                  saved_state['carry'] if saved_state is not None else None)
        cpl = MF6Coupler(mm, ctx, st, b, conv_fact=conv_fact,
                         obs_idx=obs_idx, obs_names=obs_names,
                         crr=(None if crr_opts is None else
                              dict(crr_opts, report=(cyc == 0))))
        cpl.carry_in = _carry
        cpl.steady_perc, cpl.steady_etg = steady_perc, steady_etg
        try:
            res = cpl.run(api)
            # What held for months rather than what happened once: the soil
            # at wilting point printed 8205 lines in 240 stress periods.
            MMsoil.report_tallies()
            # A non-converged / non-conserving cycle is not a usable state to
            # iterate from, so the guard runs every cycle.
            cpl.check_solution(max_discrepancy=a.max_discrepancy,
                               raise_on_fail=not a.allow_bad_budget)
        except BaseException:
            # interrupted too: what an earlier cycle reached is on disk
            if last_good is not None:
                print('\nCycle %d did not finish as a usable state. Cycle %d, '
                      'the last that passed its check, is saved as "%s".\n'
                      '   continue from it with:  spinup.strt_heads = "%s" '
                      '(and spinup.steady_means = "%s")'
                      % (cyc + 1, last_good[0], last_good[1], last_good[1],
                         last_good[1]))
            raise
        prev_heads = _final_heads()
        # THIS cycle's mean recharge/ETg: no longer the next cycle's steady
        # period (it has none -- see the periodic spin-up above), but what is
        # saved as the steady means a later run can start from. Skipped if
        # the user pinned the means explicitly.
        if not a.steady_means:
            steady_perc = res['perc'].mean(axis=0)
            steady_etg = res['etg'].mean(axis=0)
        if cycle_base and cyc + 1 < ncyc:
            # not after the last cycle: the final save below writes that one
            last_good = (cyc + 1, save_cycle_state(
                a, b, cMF, ctx, cycle_base, prev_heads, cpl, st, res,
                grid=_grid))
            print('spin-up cycle %d passed its check: saved as "%s" until the '
                  'spin-up ends' % last_good)
        if ncyc > 1:
            wt = np.nanmax(prev_heads, axis=0)          # water table per column
            if prev is not None:
                delta = float(np.nanmean(np.abs(wt - prev)))
                print('spin-up cycle %d/%d: mean |dWT| vs previous = %.3f m '
                      '(tol %.3f)' % (cyc + 1, ncyc, delta, a.spinup_tol))
                if delta < a.spinup_tol:
                    print('spin-up converged after %d cycle(s).' % (cyc + 1))
                    spin_converged = True
                    break
            else:
                print('spin-up cycle 1/%d done (baseline).' % ncyc)
            prev = wt
    # A spin-up that ran out of cycles is NOT an equilibrium, and the heads
    # it leaves were saved as 'equilibrated' all the same (2026-09-23: 6
    # cycles, |dWT| 2.46 -> 0.38 m against a 0.05 m tolerance).
    if ncyc > 1 and not spin_converged:
        print('\nWARNING: the spin-up did NOT converge in %d cycle(s): the last '
              'mean |dWT| was %.3f m against a tolerance of %.3f m.\n'
              '         The heads saved below are the last cycle\'s, not an '
              'equilibrium. Raise spinup.cycles, or start the next run from '
              'them to continue.' % (ncyc, delta, a.spinup_tol))
    check = None

    out_fn = os.path.join(a.ws, MF6Coupler.RESULTS_H5)
    with h5py.File(out_fn, 'w') as f:
        for k, v in res.items():
            f.create_dataset(k, data=v)
        f.create_dataset('cell_ij', data=np.array([(c[1], c[2]) for c in ctx.cells]))
        # true grid size: cannot be inferred from active cells alone
        f.create_dataset('grid_shape', data=np.array([cMF.nrow, cMF.ncol]))
        # each cell's area [m2], in the order of cell_ij: a catchment mean
        # is AREA-weighted, and on a mesh the cells differ by 10^5
        f.create_dataset('cell_area', data=np.asarray(cpl.area, dtype=float))
    print('\nCoupled run finished. Results: %s' % out_fn)
    _w = np.asarray(cpl.area, dtype=float)
    print('perc  mean %.4g m/d   ETg mean %.4g m/d (catchment, area-weighted)'
          '   %s iters mean %.1f'
          % (float(np.average(res['perc'], axis=1, weights=_w).mean()),
             float(np.average(res['etg'], axis=1, weights=_w).mean()),
             res.get('iters_kind', 'outer'), res['outer_iters'].mean()))
    # WP1d: open-water evaporation is MF6's now, read back from SFR SIMEVAP and
    # LAK EVAP and carried in the MM vector as iEow, so it is a measured flux
    # in the water balance rather than a structural zero.
    # KEPT APART by package (WP2 row 2): the streams lose water in transit,
    # a pond loses its storage -- and each pond's own loss is in the h5.
    _ev = getattr(cpl, 'evap_hist', None)
    if _ev is not None and np.size(_ev):
        _cells = int(np.count_nonzero(_ev.sum(axis=0)))
        if _cells:
            _m = {k: float(np.average(getattr(cpl, h), axis=1,
                                      weights=_w).mean())
                  for k, h in (('sfr', 'evap_sfr_hist'),
                               ('lak', 'evap_lak_hist'))}
            print('E_ow  %.4g mm/d catchment mean (area-weighted) = streams '
                  '(SFR SIMEVAP) %.4g + ponds (LAK EVAP) %.4g, from %d '
                  'cell(s) with open water'
                  % (_m['sfr'] + _m['lak'], _m['sfr'], _m['lak'], _cells))
            if 'lak_evap' in res:
                _lv = np.asarray(res['lak_evap'], dtype=float)
                _pl = np.asarray(cMF.perlen, dtype=float)[:_lv.shape[0]]
                print('      each pond, m3 over the run: %s'
                      % ', '.join('%d: %.4g' % (L + 1, v) for L, v in
                                  enumerate((_lv * _pl[:, None]).sum(axis=0))))
        else:
            print('E_ow  zero -- no open water evaporated (dry channels, or '
                  'the simulated-evaporation arrays are not exposed)')

    # aquifer recharge/discharge balance -- the number to watch when calibrating
    # --uzf-vks-scale: recharge reaching the water table should ~match discharge
    try:
        import flopy
        _lst = flopy.utils.Mf6ListBudget(os.path.join(a.ws, cMF.modelname.lower() + '.lst'))
        _bal = aquifer_balance(_lst.get_dataframes(diff=False)[1],
                               _lst.get_times(), active_area(cMF, b),
                               steady_first=bool(getattr(b, 'steady_first',
                                                         True)))
        for _line in balance_lines(_bal):
            print(_line)
        # A deficit is advice only over a whole year: a summer window drains
        # by nature (the 60-day June-July runs), and the lever is the panel's
        # uzf.vks_scale, not the command-line flag it used to name.
        if _bal['discharge'] - _bal['recharge'] > 5.0:
            if _bal['days'] >= 365:
                print('   -> the aquifer drains over the run; if that is not '
                      'expected, uzf.vks_scale (now %.3g) lifts recharge'
                      % a.uzf_vks_scale)
            else:
                print('   (a %.0f-day window: a deficit says nothing about '
                      'the long-term balance)' % _bal['days'])
    except Exception as _exc:                           # noqa: BLE001
        print('aquifer balance: not computed (%s: %s)'
              % (type(_exc).__name__, _exc))

    # Save the final head field for reuse as an IC. Auto-save after a spin-up
    # (so it is never lost), or on explicit --save-strt for a single run.
    save_pref = a.save_strt or ('hi_spinup' if ncyc > 1 else None)
    if save_pref:
        pref = _state_out(a, save_pref)
        paths = b.save_heads_asc(prev_heads, pref)
        # ... and the rest of where the run ended, so a run started from
        # these heads starts where this one stopped
        paths = list(paths) + [save_run_state(pref, b, cpl, st)]
        _write_state_scope(a, a.config, save_pref,      # WP0.6 scope sidecar
                           grid=_grid)
        print('%s heads saved: %s'
              % ('final' if ncyc <= 1 else 'equilibrated' if spin_converged
                 else 'NOT-converged spin-up', ', '.join(
                     os.path.basename(p) for p in paths)))
        print('   reuse with:  spinup.strt_heads = "%s"   (a run starts where '
              'this one ended, with no steady period; spinup.cycles = 1 runs '
              'it once)' % save_pref)

    # Save per-cell mean recharge / ETg so the steady state of later runs can be
    # driven by the dynamic mean (auto after spin-up, or on explicit --save-means).
    mean_pref = a.save_means or ('hi_spinup' if ncyc > 1 else None)
    if mean_pref:
        mp = _state_out(a, mean_pref)
        _write_cell_grid(res['perc'].mean(axis=0), ctx.cells, cMF.nrow, cMF.ncol,
                         cMF, mp + '_perc.asc')
        _write_cell_grid(res['etg'].mean(axis=0), ctx.cells, cMF.nrow, cMF.ncol,
                         cMF, mp + '_etg.asc')
        _write_state_scope(a, a.config, mean_pref,      # WP0.6 scope sidecar
                           grid=_grid)
        print('steady-state means saved: %s_{perc,etg}.asc' % os.path.basename(mp))
        print('   reuse with:  spinup.steady_means = "%s"' % mean_pref)
    # the spin-up's own state is on disk: the cycle kept while it ran is a
    # duplicate now (or an older cycle)
    if cycle_base and last_good is not None:
        drop_cycle_state(a, cycle_base, b.nlay)

    _run_postproc(a, cMF, ctx, res)


def _run_postproc(a, cMF, ctx, res):
    """Draw every figure for a completed run.

    Shared by the normal path and by --postproc-only, so re-drawing from an
    existing run takes exactly the same route as drawing at the end of one.
    """
    # THE OBSERVATION EXPORTS (WP6.4), figures or not: one row per stress
    # period for heads, soil moisture, streamflow and ET -- what WP7's
    # forward runs are read by. Into the results folder, and into the model
    # workspace under a name that does not change from run to run.
    try:
        from marmites_postprocess import export_observations
        export_observations(a.ws, cMF.modelname.lower(), DS, cMF, ctx, res,
                            [os.path.join(a.out_dir, '_output'),
                             os.path.join(a.ws, 'obs_exports')])
    except Exception as exc:                            # noqa: BLE001
        print('   observation exports skipped: %r' % exc)
    if a.postproc or a.preproc:
        # All results go to <ws-root>/out_<stamp>_<tag>/, never into the
        # repository and not into the model workspace either.
        from marmites_postprocess import run_preproc, run_postproc, native_suite
        os.makedirs(a.out_dir, exist_ok=True)
        # THE PLOTS PANEL'S FIGURE SETTINGS. The hydrological year, the tick
        # density and the water-balance unit are set on cMF, where the
        # figures look; the switches go to the calls below. Every one of
        # them was read by nothing -- sankey was a literal True here.
        cfg = getattr(a, 'config', None)
        if cfg is not None:
            props.apply_plot_settings(cfg, cMF)
        if a.preproc:
            run_preproc(a.ws, DS, name=cMF.modelname.lower(), out_root=a.out_dir,
                        gis_ws=a.gis_ws,
                        cMF=cMF, ctx=ctx, res=res,
                        input_maps=(cfg.postproc.input_maps
                                    if cfg is not None else True))
        if a.postproc:
            run_postproc(a.ws, DS, name=cMF.modelname.lower(), out_root=a.out_dir)
            # native MARMITESplot figures, driven by the in-memory coupled data
            native_suite(os.path.join(a.out_dir, '_output'), cMF, ctx, res,
                         ds_ws=DS, sim_ws=a.ws,
                         sankey=(cfg.postproc.sankey if cfg is not None
                                 else True),
                         sankey_full=a.sankey_full,
                         sankey_min_flux=a.sankey_min_flux, map_days=a.map_days,
                         sankey_obs_years=a.sankey_obs_years,
                         obs_series=(cfg.postproc.obs_series
                                     if cfg is not None else True),
                         result_maps=(cfg.postproc.result_maps
                                      if cfg is not None else True))
            # 01-07 water-budget figures incl. 06_heads/07_coupling and the
            # NWT-vs-MF6 comparison (into <out_dir>/figures_nwt_comparison/)
            try:
                import plot_water_budget as pwb
                # The MODFLOW-NWT reference is a 65 x 60 structured run. On a
                # mesh it used to be switched off altogether -- the catchment
                # series and totals too, which need no grid at all. The maps
                # now put the mesh on that grid by an area-weighted overlay
                # (plot_water_budget._mesh_to_grid), so it is compared again.
                pwb.make_figures(a.ws, out_dir=a.out_dir)
            except Exception as exc:
                print('   plot_water_budget skipped: %r' % exc)
        print('results written to %s' % a.out_dir)


if __name__ == '__main__':
    main()
