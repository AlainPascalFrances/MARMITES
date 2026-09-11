# -*- coding: utf-8 -*-
"""WP1d -- the stream network from the MAPPED LINES, not from a raster.

Until now SFR found its cells by looking for ``inputSTREAMw.asc > 0``. That
raster was never a channel map: its values came from ``Soil_type.shp``, where
``PONDw`` is 1.5 m on the two alluvium polygons and 0 everywhere else. So what
the model called "the stream network" was the alluvium footprint, and the
channel width was one number for the whole catchment.

The network now comes from the hydrography the modeller mapped --
``inputSTREAM.csv``, the grid-independent vertex table WP1 writes from
``hydrography.shp`` -- burned onto whichever grid panel 1 produced. A cell is
a stream cell when a mapped line actually crosses it, and the reach takes the
attributes of the segment that contributes the most length to it.

Width
-----
``[sfr] width`` declares where the value comes from, and all three producers
are honoured here:

    value     one number for the whole network
    column    a per-segment attribute, carried by inputSTREAM_param.csv
    drainage  w = a * A**b, with A the contributing area

The drainage law was declared in the configuration from WP0 and never
resolved -- the converter cannot compute contributing area, because that
depends on the routing the model does on its own grid. It is resolved HERE,
after ``stream_network`` has routed the reaches, from ``net.acc``: the number
of stream cells draining through each reach, times the cell area.
"""

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

__all__ = ['ChannelError', 'read_stream_lines', 'burn_channel',
           'channel_width', 'channel_depth']


class ChannelError(Exception):
    """The mapped stream network cannot be put on the grid."""


def _rows(path):
    """CSV rows, skipping the '#' provenance header WP1 writes."""
    import csv

    with open(path, encoding='utf-8-sig') as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith('#')]
    return list(csv.DictReader(lines))


def read_stream_lines(stream_csv, param_csv=None):
    """The mapped network as ``{seg_id: [(x, y), ...]}`` plus its parameters.

    ``inputSTREAM.csv`` carries ``seg_id, seq, x, y`` in MODEL CRS, so this is
    grid-independent: the same file serves a structured grid and a mesh.
    """
    if not os.path.exists(stream_csv):
        raise ChannelError(
            'the stream network table is missing: %s\nRun '
            'code/tools/gis_to_dataset.py to write it from the hydrography '
            'layer.' % stream_csv)
    segs = {}
    for r in _rows(stream_csv):
        segs.setdefault(int(r['seg_id']), []).append(
            (int(r['seq']), float(r['x']), float(r['y'])))
    if not segs:
        raise ChannelError('%s holds no vertices' % stream_csv)
    lines = {k: [(x, y) for _s, x, y in sorted(v)] for k, v in segs.items()}

    params = {}
    if param_csv and os.path.exists(param_csv):
        for r in _rows(param_csv):
            params[int(r['seg_id'])] = r
    return lines, params


def burn_channel(lines, grid, shape):
    """Put the mapped lines on the grid.

    Returns ``(present, seg_of_cell, length)``, each shaped like the model
    grid -- ``(nrow, ncol)``, i.e. ``(ncpl, 1)`` on a mesh. ``seg_of_cell`` is
    the segment contributing the most length to each cell (-1 where none),
    which is the right rule for a per-segment attribute: a cell crossed by two
    reaches takes the one that dominates it.
    """
    import marmites_vector as mv

    seg_ids = sorted(lines)
    best = np.zeros(grid.ncell, dtype=float)
    length = np.zeros(grid.ncell, dtype=float)
    pick = np.full(grid.ncell, -1, dtype=int)
    for sid in seg_ids:
        pts = lines[sid]
        for k in range(len(pts) - 1):
            p, q = pts[k], pts[k + 1]
            bb = (min(p[0], q[0]), min(p[1], q[1]),
                  max(p[0], q[0]), max(p[1], q[1]))
            for ic in grid.candidates(bb):
                ln, _ = mv._clip_segment_to_convex(p, q, grid.polygons[ic])
                if ln <= 0.0:
                    continue
                length[ic] += ln
                if ln > best[ic]:
                    best[ic], pick[ic] = ln, sid
    present = (length > 0.0).astype(float)
    return (present.reshape(shape), pick.reshape(shape),
            length.reshape(shape))


def arbolate_sum(net, cell_length):
    """Cumulative channel length draining through each reach [m].

    ``net.cells`` is in reach order, headwater first, so one forward pass
    down the receiver chain accumulates it.
    """
    ln = np.asarray(cell_length, dtype=float) if cell_length is not None else None
    arb = {}
    for c in net.cells:
        own = float(ln[c]) if ln is not None else 0.0
        arb[c] = arb.get(c, 0.0) + own
        r = net.recv.get(c)
        if r is not None:
            arb[r] = arb.get(r, 0.0) + arb[c]
    return arb


def channel_width(net, source, cell_area, seg_of_cell=None, params=None,
                  cell_length=None, minimum=0.5, verbose=True):
    """Width per stream cell, resolving whichever producer is configured.

    ``source`` is a ``marmites_config.ParamSource``. ``net`` is the ROUTED
    network, because the drainage producer needs contributing area and that
    is only known once the reaches are ordered.
    """
    shape = np.asarray(cell_area).shape if np.ndim(cell_area) else None
    wid = {}
    producer = source.producer() if source is not None else None

    if producer == 'value':
        for c in net.cells:
            wid[c] = float(source.value)
        how = 'uniform %.2f m' % float(source.value)

    elif producer == 'column':
        if seg_of_cell is None or not params:
            raise ChannelError(
                "sfr.width asks for column %r, but no per-segment parameter "
                "table was loaded (inputSTREAM_param.csv)" % source.column)
        col = source.column
        miss = [s for s in params if col not in params[s]]
        if miss:
            raise ChannelError(
                'sfr.width asks for column %r, which inputSTREAM_param.csv '
                'does not have; it has: %s'
                % (col, ', '.join(sorted(next(iter(params.values()))))))
        for c in net.cells:
            sid = int(np.asarray(seg_of_cell)[c])
            try:
                wid[c] = float(params[sid][col])
            except (KeyError, TypeError, ValueError):
                raise ChannelError(
                    'segment %s has no usable %r in inputSTREAM_param.csv'
                    % (sid, col))
        how = 'column %s' % col

    elif producer == 'drainage':
        d = dict(source.drainage)
        if {'w_min', 'w_max'} <= set(d):
            # The CdL method (cdl_gwf_model_fable_v2.py, "SFR width scaled by
            # drainage"): the width grows from the headwater value to the
            # outlet value with the normalised ARBOLATE SUM -- the cumulative
            # channel length draining through the reach -- not with a power of
            # contributing area. Its purpose is that "the outlet trunk is not
            # a 1 m pipe carrying the whole-catchment discharge".
            w0, w1 = float(d['w_min']), float(d['w_max'])
            p = float(d.get('power', 2.0))
            arb = arbolate_sum(net, cell_length)
            amax = max(arb.values()) if arb else 1.0
            amax = amax if amax > 0 else 1.0
            for c in net.cells:
                wid[c] = w0 + (w1 - w0) * (arb.get(c, 0.0) / amax) ** p
            how = ('arbolate sum, %.2f -> %.2f m (power %.3g, trunk drains '
                   '%.1f km)' % (w0, w1, p, amax / 1000.0))
        else:
            # A Hack-type power law on contributing AREA, in km2 -- the unit
            # matters: with A in m2 a single 2500 m2 cell already gives
            # 0.5 * 2500**0.35 = 7.7 m on La Mata, where the mapped channel is
            # 1.5 m wide.
            a = float(d['a'])
            b = float(d['b'])
            area = (float(cell_area) if shape is None
                    else float(np.mean(np.asarray(cell_area))))
            for c in net.cells:
                A = max(float(net.acc.get(c, 1)), 1.0) * area / 1e6   # km2
                wid[c] = a * A ** b
            how = 'w = %.3g * A**%.3g, A in km2' % (a, b)

    else:
        raise ChannelError(
            'sfr.width has no producer set -- give it value, column or '
            'drainage (see [sfr.width] in the configuration)')

    for c in wid:
        wid[c] = max(wid[c], minimum)
    if verbose:
        v = np.array(list(wid.values()))
        print('   channel width from %s: %.2f-%.2f m (mean %.2f)'
              % (how, v.min(), v.max(), v.mean()))
    return wid


def channel_depth(net, source, minimum=0.05, verbose=True):
    """Incision depth per stream cell, from ``[sfr] depth``.

    Replaces ``inputSTREAMhmax.asc``, which -- like the width raster -- held
    one number (1.0 m) on the alluvium polygons and zero elsewhere.
    """
    producer = source.producer() if source is not None else None
    if producer == 'value':
        d = max(float(source.value), minimum)
        if verbose:
            print('   channel incision: %.2f m (uniform)' % d)
        return {c: d for c in net.cells}
    raise ChannelError(
        'sfr.depth must be a single value for now (got producer %r); a '
        'per-segment depth needs a column in inputSTREAM_param.csv'
        % producer)


def as_array(per_cell, shape, fill=0.0):
    """A ``{(i, j): value}`` map as a grid array, for build_sfr."""
    out = np.full(shape, float(fill), dtype=float)
    for c, v in per_cell.items():
        out[c] = v
    return out
