# -*- coding: utf-8 -*-
"""The layer properties, from the TOML instead of the MODFLOW ini.  WP1d.

``__inputMF_flopy_v3_*.ini`` answers the same questions the panel
*MODFLOW aquifer layers* now asks -- the thickness of each layer and its
k, Ss and Sy -- and two files that answer the same question are one file
too many. This module takes the answers from the configuration and puts
them into the parsed ``clsMF`` object, in the SHAPE the ini would have
produced, so that every consumer downstream is unchanged:

    cMF.thick        (nlay, nrow, ncol)   float array
    cMF.hk_actual    list of nlay items, each a float or a 2-D array
    cMF.ss_actual    idem
    cMF.sy_actual    idem
    cMF.botm         (nlay, nrow, ncol)   recomputed from elev - cumulative thickness

The arithmetic is NOT reimplemented here: ``checkarray`` still turns a name
or a number into an array, and the ``botm`` loop is the one from
``clsMF.__init__``, moved rather than copied in spirit. What changes is
only WHERE the answer comes from.

A source that names a shapefile is refused rather than guessed at: wrapping
a polygon attribute onto the layers is the converter's job, and it does not
do it yet. Saying so is better than silently falling back to the ini, which
is exactly the confusion this module exists to end.
"""

import os

import numpy as np

__author__ = "Alain P. Francés <frances.alain@gmail.com>"

__all__ = ['PROPERTIES', 'apply_layer_properties', 'apply_boundaries',
           'apply_uzf', 'resolve_source', 'PropertyError']


class PropertyError(Exception):
    """A layer property cannot be resolved from the configuration."""


# (config field on [layers], the clsMF attribute the ini filled, a label)
PROPERTIES = (
    # ibound is INTEGER and must be read before the thickness: botm is
    # elev - cumulative thickness*|ibound|, so it decides where a layer
    # contributes at all.
    ('ibound', 'ibound', 'active cells'),
    ('thickness', 'thick', 'layer thickness'),
    ('k', 'hk', 'hydraulic conductivity'),
    ('k33', 'vka', 'vertical conductivity'),
    ('ss', 'ss', 'specific storage'),
    ('sy', 'sy', 'specific yield'),
)

# The per-layer FLAGS. The ini gave one integer per layer and every
# parameter set in the repository repeats the same integer nlay times, so
# the panel asks once and this fans it out.
FLAGS = (
    ('convertible', 'laytyp', lambda on: 1 if on else 0),
    ('k33_as_ratio', 'layvka', lambda on: 1 if on else 0),
)


def resolve_source(src, nlay, dataset_dir, what):
    """One ``VectorSource`` as the per-layer list the ini would have given.

    Returns a list of length ``nlay`` whose items are either floats (one
    number for the whole layer) or ABSOLUTE raster paths. Absolute is
    deliberate: ``checkarray`` joins what it gets with ``MF_ws``, and
    ``os.path.join`` returns an absolute second argument unchanged, so the
    rasters can live anywhere under the dataset rather than only in MF_ws.
    """
    producer = src.producer()
    if producer is None:
        return None                      # not answered: the ini still supplies it
    if producer == 'layer':
        raise PropertyError(
            '%s names the shapefile %r, and wrapping a polygon attribute onto '
            'the layers is not implemented yet -- give a raster or a single '
            'value instead.' % (what, src.layer))
    if producer == 'value':
        return [float(src.value)] * nlay
    names = src.rasters(nlay)
    if len(names) != nlay:
        # One name for every layer is allowed; anything else is a pattern
        # that does not fit the stack, and guessing which layer it meant is
        # how a model ends up with layer 2 holding layer 1's numbers.
        if len(names) == 1:
            names = names * nlay
        else:
            raise PropertyError(
                '%s expands to %d raster(s) for %d layer(s)'
                % (what, len(names), nlay))
    out = []
    for name in names:
        path = name if os.path.isabs(name) else os.path.join(dataset_dir, name)
        _must_exist(path, what)
        out.append(path)
    return out


def _must_exist(path, what):
    """Present AND spelled the way it is on disk.

    CASE-STRICT on purpose: Windows opens Ss_l2.asc when the file is called
    ss_l2.asc and Linux does not, so a check that trusted the filesystem
    would let a model run here and fail there -- and a per-layer pattern is
    exactly what makes that mistake easy to write.
    """
    folder, base = os.path.dirname(path), os.path.basename(path)
    try:
        here = os.listdir(folder)
    except OSError:
        raise PropertyError('%s: %s is not a folder' % (what, folder))
    if base in here:
        return
    near = [h for h in here if h.lower() == base.lower()]
    if near:
        raise PropertyError(
            '%s: %s does not exist -- the file on disk is spelled %r. Windows '
            'would open it and Linux would not, so the name has to match.'
            % (what, path, near[0]))
    raise PropertyError('%s: %s does not exist' % (what, path))


def apply_layer_properties(cfg, cMF, dataset_dir, verbose=True):
    """Put the configured layer properties into ``cMF``. Returns what it did.

    Called straight after ``clsMF`` is built, before the cell list, the soil
    model or the MF6 packages read anything from it. Every property the
    configuration does not answer is left exactly as the ini parsed it.
    """
    if cfg is None:
        return []
    done = []

    # The flags first: one answer per model, fanned out over the layers the
    # way the ini spelled it out. They are unconditional -- unlike the
    # sources there is no "not answered" state for a boolean, so the panel
    # is always the authority.
    for field, attr, to_int in FLAGS:
        if hasattr(cfg.layers, field):
            setattr(cMF, attr,
                    [to_int(getattr(cfg.layers, field))] * int(cMF.nlay))

    for field, attr, label in PROPERTIES:
        src = getattr(cfg.layers, field, None)
        if src is None:
            continue
        values = resolve_source(src, int(cMF.nlay), dataset_dir,
                                'layers.%s' % field)
        if values is None:
            continue
        setattr(cMF, attr, values)
        done.append((field, src.producer(),
                     src.raster or ('%g' % src.value if src.value is not None
                                    else '')))

    if not done:
        return done

    # Re-run the SAME conversion the constructor runs, so a raster read here
    # and a raster read there cannot diverge.
    touched = {attr for field, attr, _ in PROPERTIES
               if field in [d[0] for d in done]}
    for attr in ('hk', 'vka', 'ss', 'sy'):
        if attr in touched:
            setattr(cMF, attr + '_actual',
                    cMF.cPROCESS.checkarray(getattr(cMF, attr)))
    if 'ibound' in touched:
        # INTEGER, and reshaped the way the constructor reshapes it, so a
        # one-layer model gets (1, nrow, ncol) rather than (nrow, ncol).
        ib = cMF.cPROCESS.checkarray(cMF.ibound, dtype=int)
        ib = np.asarray(ib)
        if int(cMF.nlay) < 2:
            ib = ib.reshape((1, cMF.nrow, cMF.ncol))
        cMF.ibound = ib
    if 'thick' in touched:
        cMF.thick = cMF.cPROCESS.float2array(
            cMF.cPROCESS.checkarray(cMF.thick))
        _recompute_botm(cMF)

    if verbose:
        for field, producer, what in done:
            print('layers.%s: %s from the panel%s'
                  % (field, producer, (' (%s)' % what) if what else ''))
    return done


def _recompute_botm(cMF):
    """``botm`` from the land surface and the cumulative thickness.

    The loop from ``clsMF.__init__``: each layer's bottom is the elevation
    minus every thickness down to and including its own, and an inactive
    cell contributes nothing, so a layer that is absent there does not push
    the ones below it down.
    """
    elev = np.ma.masked_values(np.asarray(cMF.elev), cMF.hnoflo, atol=0.09)
    ibound = np.abs(np.asarray(cMF.ibound))
    botm, cum = [], None
    for l in range(int(cMF.nlay)):
        layer = cMF.thick[l, :, :] * ibound[l, :, :]
        cum = layer if cum is None else cum + layer
        botm.append(np.ma.masked_values(elev - cum, cMF.hnoflo, atol=0.09))
    cMF.botm = np.asarray(botm)
    if int(cMF.nlay) < 2 and isinstance(cMF.botm, list):
        cMF.botm = np.ma.masked_values(
            np.asarray(cMF.botm).reshape((1, cMF.nrow, cMF.ncol)),
            cMF.hnoflo, atol=0.09)


# =====================================================================
#  The boundary packages: GHB and DRN
# =====================================================================
# The parameter file kept the FOOTPRINT of each boundary inside its
# rasters -- zero everywhere except on the boundary -- and two of its
# conventions only in a comment: that a negative drain elevation means the
# bottom of the layer, and that a conductance equal to hnoflo means "not
# here". The configuration says the first out loud (drn.at_layer_base) and
# keeps the second, because it is how nodata reaches these arrays.
#
# The loops below are the ones from clsMF.__init__, and a test requires the
# list they build to equal the one the parameter file built, entry for
# entry, on the real La Mata drains.


def apply_boundaries(cfg, cMF, dataset_dir, verbose=True):
    """Build GHB and DRN from the configuration. Returns what it did."""
    if cfg is None:
        return []
    done = []
    for name in ('ghb', 'drn'):
        pkg = getattr(cfg, name, None)
        if pkg is None:
            continue
        setattr(cMF, name + '_yn', 1 if pkg.enable else 0)
        if not pkg.enable:
            continue
        builder = _build_ghb if name == 'ghb' else _build_drn
        n = builder(pkg, cMF, dataset_dir)
        done.append((name, n))
        if verbose:
            print('%s: %d cell(s) on layer(s) %s from the panel'
                  % (name.upper(), n,
                     ', '.join(str(L) for L in pkg.layers)))
    return done


def _layer_arrays(src, layers, cMF, dataset_dir, what):
    """{layer index, 0-based: 2-D array} for the layers the package covers.

    A layer the package does not list gets no array and therefore no
    boundary cells -- which is the panel's way of saying what the legacy
    rasters said with a plane of zeros.
    """
    nlay = int(cMF.nlay)
    values = resolve_source(src, nlay, dataset_dir, what)
    if values is None:
        raise PropertyError('%s: nothing to read' % what)
    out = {}
    for L in layers:
        idx = int(L) - 1                      # MODFLOW counts layers from 1
        if not (0 <= idx < nlay):
            raise PropertyError('%s: layer %s is outside 1..%d'
                                % (what, L, nlay))
        v = values[idx]
        plane = np.zeros((cMF.nrow, cMF.ncol), dtype=float)
        if isinstance(v, str):
            plane = cMF.cPROCESS.convASCIIraster2array(v, plane)
        else:
            plane[:, :] = float(v)
        out[idx] = plane
    return out


def _build_ghb(pkg, cMF, dataset_dir):
    """cMF.ghb_head_array, ghb_cond_array and layer_row_column_head_cond."""
    nlay, nrow, ncol = int(cMF.nlay), cMF.nrow, cMF.ncol
    head = _layer_arrays(pkg.head, pkg.layers, cMF, dataset_dir, 'ghb.head')
    cond = _layer_arrays(pkg.cond, pkg.layers, cMF, dataset_dir, 'ghb.cond')
    cMF.ghb_head_array = np.zeros((nlay, nrow, ncol))
    cMF.ghb_cond_array = np.zeros((nlay, nrow, ncol))
    cMF.layer_row_column_head_cond = {0: []}
    for l in sorted(head):
        cMF.ghb_head_array[l, :, :] = head[l]
        cMF.ghb_cond_array[l, :, :] = cond[l]
        for i in range(nrow):
            for j in range(ncol):
                h = cMF.ghb_head_array[l, i, j]
                if h == 0 or h == cMF.hnoflo:
                    continue
                if cMF.ghb_cond_array[l, i, j] == cMF.hnoflo:
                    continue
                cMF.layer_row_column_head_cond[0].append(
                    [l, i, j, h, cMF.ghb_cond_array[l, i, j]])
    return len(cMF.layer_row_column_head_cond[0])


def _build_drn(pkg, cMF, dataset_dir):
    """cMF.drn_elev_array, drn_cond_array and the elevation/cond list.

    ``at_layer_base`` puts every drain just above the bottom of its own
    layer (botm + 0.01 m). The parameter file spelled that as a NEGATIVE
    elevation in the raster, and a negative elevation is still honoured
    cell by cell, so a map that mixes the two keeps working.
    """
    nlay, nrow, ncol = int(cMF.nlay), cMF.nrow, cMF.ncol
    elev = _layer_arrays(pkg.elevation, pkg.layers, cMF, dataset_dir,
                         'drn.elevation')
    cond = _layer_arrays(pkg.cond, pkg.layers, cMF, dataset_dir, 'drn.cond')
    cMF.drn_elev_array = np.zeros((nlay, nrow, ncol))
    cMF.drn_cond_array = np.zeros((nlay, nrow, ncol))
    cMF.layer_row_column_elevation_cond = {0: []}
    for l in sorted(elev):
        cMF.drn_elev_array[l, :, :] = elev[l]
        cMF.drn_cond_array[l, :, :] = cond[l]
        for i in range(nrow):
            for j in range(ncol):
                e = cMF.drn_elev_array[l, i, j]
                if e == 0:
                    continue
                if e != cMF.hnoflo and (pkg.at_layer_base or e < 0):
                    base = (cMF.botm[l] if isinstance(cMF.botm[l], float)
                            else cMF.botm[l, i, j])
                    e_out = base + 0.01
                else:
                    e_out = e + 0.01
                if cMF.drn_cond_array[l, i, j] != cMF.hnoflo:
                    cMF.layer_row_column_elevation_cond[0].append(
                        [l, i, j, e_out, cMF.drn_cond_array[l, i, j]])
                cMF.drn_elev_array[l, i, j] = e_out
    return len(cMF.layer_row_column_elevation_cond[0])


# =====================================================================
#  UZF
# =====================================================================
# The panel answers what ModflowGwfuzf takes, and this puts it where the
# build looks for it -- which is still the legacy attribute names, because
# clsMF6 reads them through getattr. The UZF1 names that MODFLOW 6 has no
# equivalent for are not set here BECAUSE THEY ARE NOT ASKED: see the Uzf
# docstring in marmites_config for where each one went.

# (config field on [uzf], the clsMF attribute the build reads)
# The plain numbers. ntrailwaves/nwavesets are counts and surfdep is one
# depth for the model, so none of them is a map.
UZF_FIELDS = (
    ('surfdep', 'surfdep'),
)
# The Brooks-Corey properties, one value per UZF OBJECT. They reach the
# build as the per-layer list checkarray understands, exactly like the
# aquifer properties -- so a single value still broadcasts to what a
# scalar produced, and a raster gives each cell its own.
UZF_SOURCES = (
    ('eps', 'eps'),
    ('thtr', 'thtr'),
    ('thts', 'thts'),
    ('thti', 'thti'),
)


def apply_uzf(cfg, cMF, dataset_dir, verbose=True):
    """Put the configured unsaturated zone into ``cMF``. Returns what it did.

    ``vks_from`` becomes the legacy ``iuzfopt``: 'raster' is 1 (read the
    map) and 'layer' is 2 (use the layer's own k33). The build tests
    ``iuzfopt != 1``, so the meaning survives the rename.
    """
    if cfg is None or getattr(cfg, 'uzf', None) is None:
        return []
    u = cfg.uzf
    done = []
    for field, attr in UZF_FIELDS:
        if not hasattr(u, field):
            continue
        setattr(cMF, attr, float(getattr(u, field)))
        done.append(field)
    for field, attr in UZF_SOURCES:
        src = getattr(u, field, None)
        if src is None:
            continue
        values = resolve_source(src, int(cMF.nlay), dataset_dir,
                                'uzf.%s' % field)
        if values is None:
            raise PropertyError('uzf.%s: nothing to read' % field)
        setattr(cMF, attr, values)
        setattr(cMF, attr + '_actual', cMF.cPROCESS.checkarray(values))
        done.append(field)
    cMF.ntrail2 = int(u.ntrailwaves)
    cMF.nsets = int(u.nwavesets)
    cMF.iuzfopt = 1 if u.vks_from == 'raster' else 2
    done.append('vks_from')
    if cMF.iuzfopt == 1:
        values = resolve_source(u.vks, int(cMF.nlay), dataset_dir, 'uzf.vks')
        if values is None:
            raise PropertyError('uzf.vks_from is raster but uzf.vks is empty')
        cMF.vks = values
        cMF.vks_actual = cMF.cPROCESS.checkarray(values)
        done.append('vks')
    if verbose:
        print('UZF: %s from the panel (vks from the %s)'
              % (', '.join(done), u.vks_from))
    return done


# =====================================================================
#  The grid the rasters stand on
# =====================================================================
# nrow, ncol, delr, delc and the origin were declared TWICE: in the
# MODFLOW parameter file, and implicitly by every raster in the dataset.
# They agreed in La Mata by luck -- all 42 rasters happen to carry the
# rectangle the ini names -- and nothing checked.
#
# The RASTERS are the authority. They are what the converter wrote onto
# the grid the Grid panel defined, so their header IS that grid; the ini
# only repeated it. Unlike the layer properties this cannot be applied
# after the fact: clsMF sizes and reads every array with the ini's nrow
# and ncol, so a disagreement has already corrupted them by the time this
# runs. It is therefore a CHECK that refuses, plus the origin, which the
# driver had hard-coded.


class GridMismatch(PropertyError):
    """The parameter file and the dataset's rasters describe different grids."""


def dataset_grid(dataset_dir):
    """``(xll, yll, nrow, ncol, cellsize)`` the dataset's rasters declare."""
    import marmites_meshes as meshes
    rect, names, others = meshes.dataset_rectangle(str(dataset_dir))
    if rect is None:
        return None, [], []
    return rect, names, others


def check_grid(cMF, dataset_dir, strict=True, verbose=True):
    """Does the parameter file describe the grid the rasters are on?

    Returns the rectangle the rasters declare, or None when the dataset
    holds no raster to compare against -- a new catchment, where there is
    nothing to disagree with yet.
    """
    rect, names, others = dataset_grid(dataset_dir)
    if rect is None:
        return None
    xll, yll, nrow, ncol, cell = rect
    delr = float(np.ravel(np.asarray(cMF.delr, dtype=float))[0])
    delc = float(np.ravel(np.asarray(cMF.delc, dtype=float))[0])
    bad = []
    if int(cMF.nrow) != int(nrow):
        bad.append('nrow %d vs %d' % (cMF.nrow, nrow))
    if int(cMF.ncol) != int(ncol):
        bad.append('ncol %d vs %d' % (cMF.ncol, ncol))
    if abs(delr - cell) > 1e-6 or abs(delc - cell) > 1e-6:
        bad.append('cell %g x %g vs %g' % (delr, delc, cell))
    if abs(float(cMF.xllcorner) - xll) > 1e-3:
        bad.append('xll %g vs %g' % (cMF.xllcorner, xll))
    if abs(float(cMF.yllcorner) - yll) > 1e-3:
        bad.append('yll %g vs %g' % (cMF.yllcorner, yll))
    if bad and strict:
        raise GridMismatch(
            'the MODFLOW parameter file and the dataset rasters describe '
            'different grids (%s). The rasters win -- they are what the '
            'converter wrote onto the grid the Grid panel defined -- but '
            'every array has already been read with the parameter file\'s '
            'shape, so this cannot be corrected here. Rebuild the dataset '
            'from the Grid panel, or fix the parameter file.\n'
            '  %d raster(s) agree on %r%s'
            % ('; '.join(bad), len(names), rect,
               ''.join('\n  %d other(s) disagree among themselves: %r'
                       % (len(n), r) for r, n in others)))
    if verbose:
        if others:
            print('WARNING: %d raster(s) do not sit on the model grid: %s'
                  % (sum(len(n) for _r, n in others),
                     ', '.join(n[0] for _r, n in others)))
        print('grid: %d x %d cells of %g m at (%g, %g), from %d dataset '
              'raster(s)%s' % (nrow, ncol, cell, xll, yll, len(names),
                               '' if not bad else ' -- MISMATCH: %s'
                               % '; '.join(bad)))
    return rect


# =====================================================================
#  The catchment polygon as the GEOGRAPHIC REFERENCE
# =====================================================================
# It does not decide which cells are active -- `layers.ibound` does, one
# map per layer, because a layer can pinch out inside the catchment and in
# La Mata one does (layer 1 absent in 84 cells where layer 2 is present,
# with a 20-35 m thickness still written there, so the thickness cannot
# express it either).
#
# What the polygon IS, is the thing every input has to agree with about
# where the model is. So it is used to CHECK: the active cells should sit
# inside it, and a model whose active cells fall largely outside it is
# almost certainly in a different coordinate system.


def catchment_mask(cfg, cMF, gis_dir):
    """The cells the catchment polygon touches: (nrow, ncol) of 0/1.

    A cell counts as inside when the polygon touches it AT ALL -- the
    modeller's rule. None when there is no boundary layer to read.
    """
    import marmites_vector as mv

    name = getattr(getattr(cfg, 'grid', None), 'boundary', '') if cfg else ''
    if not name:
        return None
    path = os.path.join(str(gis_dir), name)
    if not os.path.exists(path):
        raise PropertyError(
            'grid.boundary = %r is not in %s. It is the catchment: the '
            'geographic reference every input is checked against.'
            % (name, gis_dir))
    grid = mv.TargetGrid.from_cMF(cMF)
    mask, _report = mv.overlay_polygons(mv.Layer(path), grid, how='presence',
                                        fill=0, dtype=int)
    return np.asarray(mask, dtype=int).reshape(grid.shape)


def check_catchment(cfg, cMF, gis_dir, verbose=True):
    """Do the active cells sit inside the catchment? Returns a report.

    NOT an error by itself: a cell clipped by the boundary is a real
    modelling choice, and the two maps are allowed to differ at the edge.
    What it catches is the case that is never intentional -- active cells
    far outside the polygon, which means the rasters and the catchment are
    not in the same coordinate system.
    """
    mask = catchment_mask(cfg, cMF, gis_dir)
    if mask is None:
        return None
    active = (np.abs(np.asarray(cMF.ibound, dtype=int)) != 0).any(axis=0)
    inside = int((active & (mask != 0)).sum())
    outside = int((active & (mask == 0)).sum())
    total = int(active.sum())
    report = {'active': total, 'inside': inside, 'outside': outside,
              'catchment': int((mask != 0).sum()),
              'fraction_inside': (float(inside) / total) if total else 0.0}
    if verbose:
        print('catchment: %d of %d active cell(s) inside %s (%.1f%%); '
              'the polygon touches %d cell(s)'
              % (inside, total, cfg.grid.boundary,
                 100.0 * report['fraction_inside'], report['catchment']))
    if total and report['fraction_inside'] < 0.5:
        raise PropertyError(
            'only %d of %d active cells (%.1f%%) fall inside %s. The active '
            'cells and the catchment are not in the same place -- check the '
            'coordinate system of the rasters against the shapefile.'
            % (inside, total, 100.0 * report['fraction_inside'],
               cfg.grid.boundary))
    return report


# =====================================================================
#  The state a run starts from
# =====================================================================
# `spinup.strt_heads` names a saved, equilibrated head field. It belongs
# to the grid and layer set that produced it, so reusing it under a
# different grid would hand MODFLOW an array of the wrong length -- which
# is why a guard existed. But the guard REFUSED, and a refusal is the
# wrong answer to "this state does not fit": the run can always start from
# the land surface instead, which is what a spin-up starts from anyway.
#
# So the rule is: use the saved state when it exists AND belongs here;
# otherwise COMPUTE the initial heads from the DEM and say so. Never fall
# back silently to whatever the parameter file named, which is how a run
# ends up starting from a state nobody chose.

# head = a * elevation + b. Slightly below the land surface: starting AT
# it forces MODFLOW to drain a large excess through the first stress
# periods, which reads as slow convergence and spurious rejected
# infiltration rather than as hydrology.
DEFAULT_STRT_DEM = (1.0, -2.0)


def resolve_initial_heads(cfg, state_dir, verbose=True):
    """How this run should start. ``(kind, payload, why)``.

    kind is 'saved' (payload = the prefix), or 'dem' (payload = (a, b)).
    ``why`` is the human sentence explaining the choice -- empty when the
    saved state was simply usable.
    """
    import marmites_config as mcfg

    prefix = (getattr(cfg.spinup, 'strt_heads', '') or '').strip()
    dem = tuple(getattr(cfg.spinup, 'strt_dem', ()) or ()) or DEFAULT_STRT_DEM

    if not prefix:
        why = ('no saved state is named, so the water table starts at '
               'elevation * %g %+g m' % dem)
        if verbose:
            print('initial heads: %s' % why)
        return 'dem', dem, why

    missing = [os.path.basename(p)
               for p in _state_files(state_dir, prefix, int(cfg.layers.nlay))
               if not os.path.exists(p)]
    if missing:
        why = ('%s is named but %s not in %s, so the water table starts '
               'at elevation * %g %+g m instead'
               % (prefix, ', '.join(missing), state_dir, dem[0], dem[1]))
        if verbose:
            print('initial heads: %s' % why)
        return 'dem', dem, why

    problem = mcfg.state_problem(cfg, state_dir)
    if problem:
        why = ('%s   Starting from elevation * %g %+g m instead of refusing '
               'to run.' % (problem, dem[0], dem[1]))
        if verbose:
            print('initial heads: %s' % why)
        return 'dem', dem, why

    if verbose:
        print('initial heads: the saved state %s (it belongs to this grid '
              'and layer set)' % prefix)
    return 'saved', prefix, ''


def _state_files(state_dir, prefix, nlay):
    """The per-layer files a saved head field is stored in."""
    base = prefix if os.path.isabs(prefix) else os.path.join(str(state_dir),
                                                             prefix)
    return [('%s_l%d.asc' % (base, k + 1)) for k in range(max(int(nlay), 1))]
