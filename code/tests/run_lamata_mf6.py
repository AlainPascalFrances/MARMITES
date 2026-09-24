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
import marmites_props as props  # noqa: E402
from marmites_coupler import MF6Coupler  # noqa: E402
import marmites_config as mcfg  # noqa: E402
import marmites_surface as msurf  # noqa: E402
import mm_paths  # noqa: E402


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

    THE FRONT-END FIELD IS THE SWITCH. ``[grid] dem`` names the raster on
    panel 1; when it is set and the converter has copied it into the dataset,
    the surface comes from there. Blank -- a case that has no DEM -- and the
    model keeps the elevation raster the MF parameter file names, exactly as
    before. There is no separate flag to forget to set.

    Wrapped onto THE GRID THIS RUN USES, after any mesh projection, so the
    survey is resampled once instead of twice. On La Mata the difference
    between the two is 0.65 m rms and 4.4 m at worst.

    THE SURFACE MOVES AND THE THICKNESSES DO NOT: top and every botm are
    shifted by the same delta as elev. The DEM refines where the ground is;
    it says nothing about how thick the aquifer below it is, and re-deriving
    the layer geometry from it would invent an answer the raster does not
    hold. Cells the DEM does not reach keep the elevation they had.
    """
    import marmites_dem as mdem

    if not getattr(cfg, 'grid', None) or not cfg.grid.dem:
        return False, 'no [grid] dem: the land surface is the MF ini raster'
    path = mdem.dem_path(dataset_dir)
    if not os.path.exists(path):
        return False, ('[grid] dem is %r but %s is not in the dataset -- run '
                       'code/tools/gis_to_dataset.py. The land surface is the '
                       'MF ini raster.' % (cfg.grid.dem,
                                           os.path.basename(path)))
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


def ctx_geom_area(cMF):
    """Mean cell area [m2], DIS or DISV. The drainage law w = a*A**b needs it."""
    proj = getattr(cMF, 'mesh_proj', None)
    if proj is not None:
        import marmites_vector as mv
        return mv.TargetGrid.from_gridprops(cMF.mesh_gridprops).area
    return np.outer(np.asarray(cMF.delc, float), np.asarray(cMF.delr, float))


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
DISCHARGE_PREFIXES = ('DRN', 'GHB', 'WEL')


def aquifer_balance(cum, times, area):
    """Recharge to, and discharge from, the aquifer in mm/yr over the run.

    ``cum`` is the list file's CUMULATIVE volumes [m3], one row per saved
    step; ``times`` their times [d]. The first period is left out, as it
    always was (it is the initial state); the rates are volumes over elapsed
    time, NOT a plain mean of per-step rates -- which weighted a 1-day step
    as much as a 30-day one. Also returns the share of the recharge carried
    by the single largest step, and which step, because one pulse can make
    the whole-run figure meaningless.
    """
    t = np.asarray(times, dtype=float)
    cols = list(cum.columns)
    dis_cols = [c for c in cols if c.endswith('_OUT')
                and c.split('_')[0].split('-')[0].rstrip('0123456789')
                in DISCHARGE_PREFIXES]
    rch_cols = [c for c in cols if c in RECHARGE_TERMS]
    i0 = 1 if len(t) > 1 else 0
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
            'area_km2': area / 1e6, 'days': span,
            'peak_share': share, 'peak_time': float(t[i0 + k]),
            'peak_mm': float(step[k]) / area * 1000.0 if step.size else 0.0}


def balance_lines(bal):
    """The balance as printed at the end of a run."""
    out = ['aquifer balance over %.0f d on %.3f km2: recharge to WT %.1f '
           'mm/yr  vs  discharge %.1f mm/yr  (deficit %.1f)'
           % (bal['days'], bal['area_km2'], bal['recharge'],
              bal['discharge'], bal['discharge'] - bal['recharge']),
           '   discharge by term: ' + ', '.join(
               '%s %.1f' % (c, v) for c, v in sorted(
                   bal['discharge_terms'].items()))]
    if bal['peak_share'] > 0.5:
        out.append('   WARNING: %.0f %% of that recharge arrived in ONE step '
                   '(t = %g d, %.0f mm over the catchment) -- a pulse, not a '
                   'rate; the whole-run figure says little about the rest '
                   'of the run.' % (100.0 * bal['peak_share'],
                                    bal['peak_time'], bal['peak_mm']))
    return out


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


# ONE parameter file. La Mata used to carry two, at two vertical
# resolutions, and --nlay chose between them -- so "how many layers" was
# really "which of the two files exists". The number of layers is a panel
# field now (MODFLOW aquifer layers -> Aquifer layers), and what the file
# still supplies is being emptied section by section.
MF_INI = '__inputMF_flopy_v3_2s1L.ini'


def setup_lamata(daily=True, nsp=None, grid='dis', nlay=None,
                 cfg=None, mesh_ws=None):
    """Replicate the driver setup; returns (cMF, mm, ctx, state, top, botm).

    ``nlay`` is checked against the parameter file rather than choosing
    between two of them: the layer count is a panel field now.

    ``cfg`` (WP1c.1) enables the unstructured path. When ``cfg.grid_kind`` is
    a genuine mesh producer -- 'quadtree', later 'voronoi' -- the model is
    built on the structured raster exactly as before and then PROJECTED onto
    the mesh (``marmites_mesh.project_model``). The returned ``cMF`` then
    carries ``mesh_gridprops`` and ``mesh_proj``, and its arrays are
    ``(ncpl, 1)``; see ``marmites_mesh`` for why that shape makes the rest of
    the code work unchanged. 'structured' and 'disv' never project: they are
    the regression anchor and must stay byte-identical.
    """
    ini_fn = MF_INI
    cUTIL = MMutils.clsUTILITIES(verbose=1)
    # THE ORIGIN COMES FROM THE DATASET. It used to be two literals here,
    # repeated in the parameter file and in every test -- three copies of
    # a number that belongs to the grid the converter wrote the rasters
    # onto. The rasters carry it in their own headers, so that is where it
    # is read from; a dataset with no raster yet falls back to the
    # parameter file, which is the only thing left that knows.
    _rect, _names, _others = props.dataset_grid(DS)
    if _rect is None:
        _xll, _yll = 0.0, 0.0
        print('no raster in the dataset: the origin comes from %s' % ini_fn)
    else:
        _xll, _yll = float(_rect[0]), float(_rect[1])
    cMF = ppMF.clsMF(cUTIL, MM_ws=DS, MM_ws_out=DS, MF_ws=os.path.join(DS, 'MF_ws'),
                     MF_ini_fn=ini_fn, grid=_rect,
                     xllcorner=_xll, yllcorner=_yll)
    print('parameter set: %s (%d layer(s))' % (ini_fn, cMF.nlay))
    # ... and the SHAPE is checked against those same rasters. It cannot be
    # corrected here -- clsMF has already sized and read every array with
    # the parameter file's nrow and ncol -- so a disagreement stops the run
    # rather than producing a model quietly built on the wrong rectangle.
    props.check_grid(cMF, DS)
    # THE LAYER COUNT IS THE PANEL'S. It cannot be substituted here for the
    # same reason the grid shape cannot: the parameter file's per-layer
    # lists -- ibound and strt -- have already been read with its own nlay.
    # Once those two come from the panel as well, this becomes an override
    # instead of a check.
    if nlay is not None and int(nlay) != int(cMF.nlay):
        raise SystemExit(
            'layers.nlay = %d but %s describes %d layer(s). ibound and strt '
            'still come from that file, so the two have to agree; give the '
            'rasters for %d layers there, or set the panel to %d.'
            % (int(nlay), ini_fn, cMF.nlay, int(nlay), cMF.nlay))
    # THE FRONT-END OWNS hnoflo (WP1d, geometry). It is the value that marks
    # a cell as having nothing to report, and MARMITES masks on it as well --
    # so the two have to be the SAME number, which is why it is asked once on
    # the panel rather than twice in two files. Applied right after the ini is
    # parsed, before anything reads it.
    # THE MODEL'S NAME IS THE PANEL'S. MODFLOW 6 names every file it writes
    # after it, and it used to live in the parameter file ("lamataMM") --
    # the model was named where the modeller never looked.
    if cfg is not None:
        _model = cfg.meta.model_name(cfg.paths.case)
        if _model != str(cMF.modelname).lower():
            print('model name: %s (the parameter file said %s)'
                  % (_model, cMF.modelname))
            cMF.modelname = _model
    # ONE SENTINEL, NOT TWO -- into cMF AND the raster reader. Setting only
    # cMF's copy is what turned La Mata's twelve drains into 7800 at the
    # aquifer floor; see props.apply_hnoflo, which the tests call too.
    if cfg is not None:
        props.apply_hnoflo(cfg, cMF)
    # THE FRONT-END OWNS THE LAYER PROPERTIES TOO. thickness, k, Ss and Sy
    # are asked on the panel, and whatever it answers replaces what the ini
    # parsed -- before the cell list, the soil model or any MF6 package has
    # read them. What the panel does not answer is left as the ini had it.
    props.apply_layer_properties(cfg, cMF, DS)
    # THE BOUNDARY PACKAGES TOO. GHB and DRN are rebuilt from [ghb] and
    # [drn] -- AFTER the properties, because a drain taken at the base of
    # its layer reads botm, and botm is the panel thickness now.
    # THE CATCHMENT IS THE GEOGRAPHIC REFERENCE. It does not decide which
    # cells are active -- layers.ibound does, per layer -- but the active
    # cells have to sit inside it, and a model whose cells fall outside is
    # in a different coordinate system.
    props.check_catchment(cfg, cMF, mm_paths.GIS)
    props.apply_boundaries(cfg, cMF, DS)
    props.apply_uzf(cfg, cMF, DS)
    # LENGTHS ARE METRES, always. The panel asks for a projected CRS in
    # metres and every raster is metric, so lenuni is 2 and the conversion
    # to the millimetres MARMITES works in is fixed. It was read from the
    # parameter file, where nothing could have set it to anything else.
    cMF.lenuni = 2
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
    if cfg is not None:
        gridSOIL = props.soil_grid(cfg, cMF, DS, 'zones', kind='int')
        gridSOILthick = props.soil_grid(cfg, cMF, DS, 'thickness')
        print('soil zones: %s; soil thickness: %s -- from the panel'
              % (cfg.soil.zones.producer(), cfg.soil.thickness.producer()))
    else:
        gridSOIL = cMF.cPROCESS.inputEsriAscii(grid_fn='inputSOILzones.asc',
                                               datatype=int)
        gridSOILthick = cMF.cPROCESS.inputEsriAscii(
            grid_fn='inputSOILthick.asc', datatype=float)
    gridIRR = cMF.cPROCESS.inputEsriAscii(grid_fn='inputIRRzones.asc', datatype=int)

    # THE VEGETATION COVER FROM THE SOIL PANEL: its vegetation layer, class
    # column and class table, put onto the grid by exact area overlay. It came
    # from inputVEG1area.asc .. inputVEG3area.asc, filenames written into
    # MARMITESprocess, so the panel changed nothing. The two agree on the
    # trees to 0.01 %; on grass they do not -- the polygons give ~89 % where
    # the old raster gave 25 % -- and the polygons are the reference: the
    # summer is carried by the grass wilting in the seasonal forcing, not by
    # the cover map.
    _veg = None
    if cfg is not None:
        _veg = props.veg_cover(
            cfg, cMF, DS, NVEG,
            cache_dir=(os.path.join(os.path.dirname(mesh_ws), '_overlay')
                       if mesh_ws else None))
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
    botm_l0 = np.asarray(cMF.botm)[0]
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
        botm_l0 = np.asarray(cMF.botm)[0]
        print('projected onto the mesh: %d active cell(s) of %d'
              % (int(np.count_nonzero(cMF.outcropL > 0)), info['ncpl']))

    # ---- the land surface, from the raster panel 1 names (WP1d).
    # HERE, after the projection: the DEM is wrapped onto the grid the run
    # actually uses, so the survey is resampled once instead of twice.
    _applied, _note = _apply_dem(cMF, cfg, DS,
                                 cache_dir=os.path.join(mesh_ws or '', '_dem')
                                 if mesh_ws else None)
    print('elevation: %s' % _note)
    if _applied:
        botm_l0 = np.asarray(cMF.botm)[0]

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


def _write_state_scope(a, cfg, prefix):
    """Record WHICH configuration produced a saved state (WP0.6)."""
    import json
    payload = {'state_hash': cfg.state_hash(), 'scope': cfg.state_scope()}
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
        mode=cfg.run.mode, relax=cfg.run.relax,
        nsp=(cfg.run.nsp or None), daily=cfg.run.daily, ats=cfg.run.ats,
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
        sfr=cfg.sfr.enable,
        sfr_rhk=(cfg.sfr.rhk.value if cfg.sfr.rhk.value is not None else 0.1),
        lak=lak_source, lak_bedleak=cfg.lak.bedleak,
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
        preproc=cfg.postproc.preproc,
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
    tag = a.run_tag or ('%dlay_%s' % (a.nlay or 6, a.mode))
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
    _gp = getattr(cMF, 'mesh_gridprops', None)
    b = clsMF6(cMF, top=top, botm=botm, sim_ws=a.ws, daily=True, grid=a.grid,
               vertices=(_gp['vertices'] if _gp else None),
               cell2d=(_gp['cell2d'] if _gp else None),
               strt_from_dem=(tuple(a.strt_dem) if a.strt_dem else None))
    b.seep = a.seep
    b.ats = a.ats
    b.drn_seep_cond = float(a.seep_cond)
    b.uzf_vks_scale = float(a.uzf_vks_scale)
    if cfg is not None:
        b.uzf_thtr_from = str(cfg.uzf.thtr_from)
        b.solver = dataclasses.asdict(cfg.solver)
        b.outer_maximum = int(cfg.solver.outer_maximum)
        print('solver: %s, outer_dvclose %g m (max %d), inner_dvclose %g m, '
              'inner_rclose %g m3/d -- from the panel'
              % (cfg.solver.complexity.upper(), cfg.solver.outer_dvclose,
                 cfg.solver.outer_maximum, cfg.solver.inner_dvclose,
                 cfg.solver.inner_rclose))
    # UNSATURATED-ZONE ET. The extinction depth follows the usual rule --
    # a raster, a column of the vegetation layer, or one value -- so it is
    # resolved the way every other spatial input is, per layer and then
    # broadcast over the column.
    b.uzf_et_form = str(getattr(a, 'uzf_et_form', 'etwc'))
    if cfg is not None:
        b.uzf_extdp = props.resolve_source(
            cfg.et.extdp, int(cMF.nlay), DS, 'et.extdp')
        print('UZF ET: always on (%s), extinction depth from %s'
              % (b.uzf_et_form, cfg.et.extdp.producer()))
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
    # WHERE THE RUN STARTS. The saved state when it exists and belongs to
    # this grid and layer set; the land surface otherwise. Never the
    # parameter file's array by accident -- that is how a run ends up
    # starting from a state nobody chose.
    _kind, _payload, _why = props.resolve_initial_heads(cfg, a.state_dir)
    if _kind == 'saved':
        pref = _state_in(a, str(_payload), '_l1.asc')
        b.strt_array = b.load_heads_asc(pref)
        print('initial heads loaded from %s_l*.asc' % pref)
    else:
        b.strt_from_dem = tuple(_payload)
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
    spin_converged, delta = False, float('nan')
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
        # What held for months rather than what happened once: the soil at
        # wilting point printed 8205 lines in 240 stress periods before this.
        MMsoil.report_tallies()
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

    out_fn = os.path.join(a.ws, '_coupled_%s.h5' % a.mode)
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
          '   outer iters mean %.1f'
          % (float(np.average(res['perc'], axis=1, weights=_w).mean()),
             float(np.average(res['etg'], axis=1, weights=_w).mean()),
             res['outer_iters'].mean()))
    # WP1d: open-water evaporation is MF6's now, read back from SFR SIMEVAP and
    # LAK EVAP and carried in the MM vector as iEow, so it is a measured flux
    # in the water balance rather than a structural zero.
    _ev = getattr(cpl, 'evap_hist', None)
    if _ev is not None and np.size(_ev):
        _cells = int(np.count_nonzero(_ev.sum(axis=0)))
        if _cells:
            print('E_ow  %.4g mm/d catchment mean, from %d cell(s) with open '
                  'water (SFR SIMEVAP + LAK EVAP)'
                  % (float(np.mean(np.sum(_ev, axis=1))) / _ev.shape[1], _cells))
        else:
            print('E_ow  zero -- no open water evaporated (dry channels, or '
                  'the simulated-evaporation arrays are not exposed)')

    # aquifer recharge/discharge balance -- the number to watch when calibrating
    # --uzf-vks-scale: recharge reaching the water table should ~match discharge
    try:
        import flopy
        _lst = flopy.utils.Mf6ListBudget(os.path.join(a.ws, cMF.modelname.lower() + '.lst'))
        _bal = aquifer_balance(_lst.get_dataframes(diff=False)[1],
                               _lst.get_times(), active_area(cMF, b))
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
        _write_state_scope(a, a.config, save_pref)      # WP0.6 scope sidecar
        print('%s heads saved: %s'
              % ('final' if ncyc <= 1 else 'equilibrated' if spin_converged
                 else 'NOT-converged spin-up', ', '.join(
                     os.path.basename(p) for p in paths)))
        print('   reuse with:  spinup.strt_heads = "%s"   (skips the spin-up)' % save_pref)

    # Save per-cell mean recharge / ETg so the steady state of later runs can be
    # driven by the dynamic mean (auto after spin-up, or on explicit --save-means).
    mean_pref = a.save_means or ('hi_spinup' if ncyc > 1 else None)
    if mean_pref:
        mp = _state_out(a, mean_pref)
        _write_cell_grid(res['perc'].mean(axis=0), ctx.cells, cMF.nrow, cMF.ncol,
                         cMF, mp + '_perc.asc')
        _write_cell_grid(res['etg'].mean(axis=0), ctx.cells, cMF.nrow, cMF.ncol,
                         cMF, mp + '_etg.asc')
        _write_state_scope(a, a.config, mean_pref)      # WP0.6 scope sidecar
        print('steady-state means saved: %s_{perc,etg}.asc' % os.path.basename(mp))
        print('   reuse with:  spinup.steady_means = "%s"' % mean_pref)

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
                pwb.make_figures(a.ws, mode=a.mode, out_dir=a.out_dir)
            except Exception as exc:
                print('   plot_water_budget skipped: %r' % exc)
        print('results written to %s' % a.out_dir)


if __name__ == '__main__':
    main()
