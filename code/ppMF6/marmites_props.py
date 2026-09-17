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

__all__ = ['PROPERTIES', 'apply_layer_properties', 'resolve_source',
           'PropertyError']


class PropertyError(Exception):
    """A layer property cannot be resolved from the configuration."""


# (config field on [layers], the clsMF attribute the ini filled, a label)
PROPERTIES = (
    ('thickness', 'thick', 'layer thickness'),
    ('k', 'hk', 'hydraulic conductivity'),
    ('ss', 'ss', 'specific storage'),
    ('sy', 'sy', 'specific yield'),
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
    for attr in ('hk', 'ss', 'sy'):
        if attr in touched:
            setattr(cMF, attr + '_actual',
                    cMF.cPROCESS.checkarray(getattr(cMF, attr)))
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
