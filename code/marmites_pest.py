# -*- coding: utf-8 -*-
"""The model side of the PEST++ calibration (WP7).

Two things a forward run hands to the model, and nothing PEST-specific
beyond them -- no pyemu here, so the model keeps running where pyemu is not
installed:

``apply_multipliers(cMF, folder)``
    Multiplier files written by the forward run, applied to the aquifer on
    the grid the model RUNS on -- after the 50 m rasters were projected onto
    the mesh -- so a pilot point multiplies the mesh cell it falls in, not a
    50 m cell that is resampled afterwards:

        k_l<L>.txt    horizontal K of layer L: one value per cell, or one
        sy_l<L>.txt   specific yield of layer L: idem
        ss_l<L>.txt   specific storage of layer L: idem
        drn_l<L>.txt  the outlet drains of layer L: one value

    and one that is not a multiplier (2026-10-18, user):

        k_l1_max_ratio.txt  layer 1's K at most this times layer 2's, cell
                            by cell (0.1: the saprolite an order of
                            magnitude below the fissured layer)

    (The specific retention is a layer property of the configuration since
    2026-10-20, layers.sr: the forward run writes it there, as Kh/Kv.)

    One value per line; per cell in the model's own cell order (the mesh's
    icell2d, the structured grid's row-major order). A missing file leaves
    that property as it is. k33 follows K (layers.k33_as_ratio), UZF's vks
    follows k33 (uzf.vks_from = layer), UZF's thts = Sy + Sr follows Sy:
    the multiplied values reach them all.

    When Sy changes, ``cMF.uzf_carry_ref`` keeps the model's (thtr, thts)
    before it: the saved state's UZF water content was left with those, and
    is carried relative to them (clsMF6.uzf_thti_carry_ref) -- for a state
    saved without its own.

``set_param(cfg, dotted, value)``
    A configuration value by its dotted name, list items included --
    ``soil.horizon[2].ks``, ``drn.cond.value`` -- which ``--set`` cannot
    address. The forward run writes the worker's configuration with it.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import os
import re

import numpy as np

__all__ = ['PestError', 'apply_multipliers', 'set_param', 'get_param',
           'MULT_PROPS']

# file stem -> the cMF attribute it multiplies
MULT_PROPS = {'k': 'hk_actual', 'sy': 'sy_actual', 'ss': 'ss_actual'}


class PestError(Exception):
    """A PEST hand-over the model cannot use as given."""


def _read(path):
    v = np.loadtxt(path, dtype=float, ndmin=1).ravel()
    if not v.size or not np.all(np.isfinite(v)) or np.any(v <= 0.0):
        raise PestError('%s: a multiplier must be finite and > 0 (got '
                        '%d value(s), min %s)' % (path, v.size,
                                                  v.min() if v.size else '-'))
    return v


def _layers3d(val, nlay, nrow, ncol):
    """A layer property as (nlay, nrow, ncol) float -- it may be a list of
    per-layer floats or arrays, or an array."""
    if isinstance(val, (list, tuple)):
        return np.stack([np.asarray(v, dtype=float) * np.ones((nrow, ncol))
                         for v in val])
    a = np.asarray(val, dtype=float)
    if a.ndim == 0:
        return np.full((nlay, nrow, ncol), float(a))
    if a.ndim == 2:
        return np.broadcast_to(a, (nlay, nrow, ncol)).copy()
    return a.reshape(nlay, nrow, ncol).copy()


def apply_multipliers(cMF, folder, out=print):
    """Multiply the aquifer properties (and the outlet drains) by the files
    in ``folder``; returns ``{name: (min, max)}`` of what was applied."""
    if not folder:
        return {}
    if not os.path.isdir(folder):
        raise PestError('the PEST multiplier folder %s is not there' % folder)
    nlay, nrow, ncol = int(cMF.nlay), int(cMF.nrow), int(cMF.ncol)
    ncell = nrow * ncol
    done = {}
    if getattr(cMF, 'sr_actual', None) is not None and any(
            os.path.exists(os.path.join(folder, 'sy_l%d.txt' % (L + 1)))
            for L in range(nlay)):
        # what the saved state's UZF water content is relative to: the
        # model's water contents before the Sy multiplier (thtr = Sr,
        # thts = Sy + Sr)
        sr0 = _layers3d(cMF.sr_actual, nlay, nrow, ncol)
        cMF.uzf_carry_ref = (sr0, _layers3d(cMF.sy_actual, nlay, nrow,
                                            ncol) + sr0)
    for stem, attr in MULT_PROPS.items():
        files = [os.path.join(folder, '%s_l%d.txt' % (stem, L + 1))
                 for L in range(nlay)]
        if not any(os.path.exists(f) for f in files):
            continue
        arr = _layers3d(getattr(cMF, attr), nlay, nrow, ncol)
        for L, f in enumerate(files):
            if not os.path.exists(f):
                continue
            m = _read(f)
            if m.size == 1:
                arr[L] *= m[0]
            elif m.size == ncell:
                arr[L] *= m.reshape(nrow, ncol)
            else:
                raise PestError('%s: %d value(s), the grid has %d cell(s) '
                                'per layer (or give one value)'
                                % (f, m.size, ncell))
            done['%s_l%d' % (stem, L + 1)] = (float(m.min()), float(m.max()))
        setattr(cMF, attr, [arr[L] for L in range(nlay)])
    # the order of the layers' K: layer 1 at most `ratio` x layer 2
    rf = os.path.join(folder, 'k_l1_max_ratio.txt')
    if os.path.exists(rf) and nlay >= 2:
        ratio = float(_read(rf)[0])
        hk = _layers3d(cMF.hk_actual, nlay, nrow, ncol)
        over = hk[0] > ratio * hk[1]
        hk[0] = np.where(over, ratio * hk[1], hk[0])
        cMF.hk_actual = [hk[L] for L in range(nlay)]
        if out:
            out('PEST: layer 1 K at most %g x layer 2 K -- %d cell(s) '
                'brought down to it' % (ratio, int(over.sum())))
    # the outlet drains, one multiplier per layer
    recs = getattr(cMF, 'layer_row_column_elevation_cond', None)
    dfiles = {L: os.path.join(folder, 'drn_l%d.txt' % (L + 1))
              for L in range(nlay)}
    if recs and any(os.path.exists(f) for f in dfiles.values()):
        mult = {L: (_read(f)[0] if os.path.exists(f) else 1.0)
                for L, f in dfiles.items()}
        new = []
        for (l, i, j, e, c) in recs[0]:
            new.append((l, i, j, e, float(c) * float(mult.get(int(l), 1.0))))
        recs[0] = new
        for L, f in dfiles.items():
            if os.path.exists(f):
                done['drn_l%d' % (L + 1)] = (mult[L], mult[L])
    if out and done:
        out('PEST multipliers from %s, change from the reference model: %s'
            % (folder, ', '.join(
                '%s %s' % (k, (pct_change(lo) if lo == hi else
                               '%s .. %s' % (pct_change(lo), pct_change(hi))))
                for k, (lo, hi) in sorted(done.items()))))
    return done


def pct_change(m):
    """A multiplier as the change it makes to the reference model, in
    percent (user, 2026-10-21): x0.87 -> '-13.0 %', x21.7 -> '+2070 %'."""
    p = (float(m) - 1.0) * 100.0
    if abs(p) < 0.05:
        return '0 %'
    return ('%+.0f %%' if abs(p) >= 100.0 else '%+.1f %%') % p


_TOKEN = re.compile(r'^([A-Za-z_][A-Za-z_0-9]*)(?:\[(\d+)\])?$')


def _walk(cfg, dotted):
    parts = dotted.split('.')
    obj = cfg
    for k, part in enumerate(parts):
        m = _TOKEN.match(part)
        if not m:
            raise PestError('%r: cannot read %r' % (dotted, part))
        name, idx = m.group(1), m.group(2)
        last = k == len(parts) - 1
        if not hasattr(obj, name):
            raise PestError('%r: no %r in %s' % (dotted, name,
                                                 type(obj).__name__))
        if last and idx is None:
            return obj, name, None
        val = getattr(obj, name)
        if idx is not None:
            i = int(idx)
            if not isinstance(val, list) or i >= len(val):
                raise PestError('%r: %s has no item %d' % (dotted, name, i))
            if last:
                return val, i, 'item'
            val = val[i]
        obj = val
    raise PestError('%r: empty name' % dotted)


def get_param(cfg, dotted):
    obj, key, kind = _walk(cfg, dotted)
    return obj[key] if kind == 'item' else getattr(obj, key)


def set_param(cfg, dotted, value):
    """Set ``dotted`` (``soil.horizon[2].ks``) to ``value``, typed like the
    value it replaces."""
    obj, key, kind = _walk(cfg, dotted)
    old = obj[key] if kind == 'item' else getattr(obj, key)
    if isinstance(old, bool):
        raise PestError('%r is a switch, not a parameter' % dotted)
    if isinstance(old, int):
        value = int(round(float(value)))
    elif old is None or isinstance(old, float):
        value = float(value)
    else:
        raise PestError('%r is %s, not a number' % (dotted,
                                                    type(old).__name__))
    if kind == 'item':
        obj[key] = value
    else:
        setattr(obj, key, value)
    return old
