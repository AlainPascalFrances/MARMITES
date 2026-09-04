# -*- coding: utf-8 -*-
"""Aggregate MODFLOW numerical layers into hydrogeological units.

The La Mata MF ini already declares the intent:

    Mnlay = 2
    Mlay  = [1, 1, 1, 2, 2, 2]

i.e. the six numerical layers represent **two** hydrogeological units. Running
one numerical layer per geologic unit (as the CdL reference model does) cuts
the cell count and the number of UZF objects by a factor of three and makes the
numerical layering match the conceptual one.

Aggregation rules (applied per unit, over the layers it contains):

    botm      bottom of the unit's lowest layer          geometry preserved
    k         thickness-weighted ARITHMETIC mean         preserves transmissivity
    k33       thickness-weighted HARMONIC mean           preserves vertical resistance
    ss        thickness-weighted arithmetic mean         storage preserved
    sy        value of the unit's UPPERMOST layer        sy acts at the water table
    strt      value of the unit's uppermost layer        head at the top of the unit
    ibound    active if ANY constituent layer is active  keeps the footprint
    thick     sum of the constituent thicknesses         geometry preserved

The aggregation is applied to the parsed ``clsMF`` object **before** the
MARMITES cell list is built, so the soil model and the groundwater model always
agree on the layering and on ``outcropL``.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import numpy as np

__all__ = ['aggregate_layers', 'layer_groups']


def layer_groups(Mlay, nlay):
    """[[layer indices of unit 1], [unit 2], ...] from the 1-based Mlay list."""
    Mlay = np.asarray(Mlay, dtype=int).ravel()
    if Mlay.size != nlay:
        raise ValueError('Mlay has %d entries but the model has %d layers'
                         % (Mlay.size, nlay))
    units = sorted(set(int(u) for u in Mlay))
    groups = [[k for k in range(nlay) if int(Mlay[k]) == u] for u in units]
    for u, g in zip(units, groups):
        if not g:
            raise ValueError('hydrogeological unit %d has no layers' % u)
        if g != list(range(g[0], g[-1] + 1)):
            raise ValueError('unit %d is not a contiguous block of layers: %s'
                             % (u, g))
    return groups


def _as3d(val, nlay, nrow, ncol):
    """Normalize scalar / per-layer / 2-D / 3-D input to (nlay, nrow, ncol)."""
    if isinstance(val, list):
        try:
            return np.stack([np.asarray(v, dtype=float) * np.ones((nrow, ncol))
                             for v in val])
        except Exception:
            val = np.asarray(val, dtype=float)
    a = np.asarray(val, dtype=float)
    out = np.ones((nlay, nrow, ncol), dtype=float)
    if a.ndim == 0:
        out *= float(a)
    elif a.ndim == 1 and a.shape[0] == nlay:
        for k in range(nlay):
            out[k] *= a[k]
    elif a.ndim == 2:
        out[:] = a[np.newaxis]
    elif a.ndim == 3:
        out[:] = a
    else:
        raise ValueError('cannot normalize array of shape %s' % (a.shape,))
    return out


def _harmonic(vals, wts):
    """Thickness-weighted harmonic mean (preserves vertical resistance)."""
    vals = np.asarray(vals, dtype=float)
    wts = np.asarray(wts, dtype=float)
    tot = wts.sum(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        res = np.where(vals > 0, wts / np.where(vals > 0, vals, 1.0), 0.0).sum(axis=0)
        out = np.where(res > 0, tot / np.where(res > 0, res, 1.0), 0.0)
    return out


def aggregate_layers(cMF, verbose=True):
    """Collapse cMF's numerical layers onto its hydrogeological units.

    Mutates ``cMF`` in place (nlay, ibound, botm, thick, strt, hk/vka/ss/sy,
    iuzfbnd, Mlay/Mnlay) and returns a short report dict.
    """
    nlay, nrow, ncol = int(cMF.nlay), int(cMF.nrow), int(cMF.ncol)
    groups = layer_groups(cMF.Mlay, nlay)
    nnew = len(groups)
    if nnew == nlay:
        return {'aggregated': False, 'nlay': nlay, 'groups': groups}

    ib = np.abs(_as3d(cMF.ibound, nlay, nrow, ncol))
    botm = _as3d(cMF.botm, nlay, nrow, ncol)
    thick = _as3d(cMF.thick, nlay, nrow, ncol)
    hk = _as3d(cMF.hk_actual, nlay, nrow, ncol)
    vka = _as3d(cMF.vka_actual, nlay, nrow, ncol)
    ss = _as3d(cMF.ss_actual, nlay, nrow, ncol)
    sy = _as3d(cMF.sy_actual, nlay, nrow, ncol)
    strt = _as3d(cMF.strt, nlay, nrow, ncol)

    # layvka semantics: VKA is the hk/vk RATIO when layvka != 0 -> convert to a
    # true vertical K before averaging, otherwise the harmonic mean is meaningless
    layvka = np.asarray(cMF.layvka, dtype=int).ravel()
    k33 = np.empty_like(hk)
    for k in range(nlay):
        if layvka[k] != 0:
            with np.errstate(divide='ignore', invalid='ignore'):
                k33[k] = np.where(vka[k] > 0, hk[k] / vka[k], hk[k])
        else:
            k33[k] = vka[k]

    n_ib = np.zeros((nnew, nrow, ncol))
    n_botm = np.zeros((nnew, nrow, ncol))
    n_thick = np.zeros((nnew, nrow, ncol))
    n_hk = np.zeros((nnew, nrow, ncol))
    n_k33 = np.zeros((nnew, nrow, ncol))
    n_ss = np.zeros((nnew, nrow, ncol))
    n_sy = np.zeros((nnew, nrow, ncol))
    n_strt = np.zeros((nnew, nrow, ncol))

    for u, g in enumerate(groups):
        w = thick[g]                                  # (nsub, nrow, ncol)
        act = ib[g] > 0
        w_eff = np.where(act, np.maximum(w, 0.0), 0.0)
        tot = w_eff.sum(axis=0)
        safe = np.where(tot > 0, tot, 1.0)
        n_ib[u] = act.any(axis=0).astype(float)       # active if ANY sub-layer is
        n_botm[u] = botm[g[-1]]                       # bottom of the lowest layer
        n_thick[u] = tot
        n_hk[u] = (hk[g] * w_eff).sum(axis=0) / safe            # arithmetic
        n_k33[u] = _harmonic(k33[g], w_eff)                     # harmonic
        n_ss[u] = (ss[g] * w_eff).sum(axis=0) / safe            # arithmetic
        n_sy[u] = sy[g[0]]                                      # uppermost layer
        n_strt[u] = strt[g[0]]
        # cells inactive everywhere in the unit: keep a finite, harmless value
        dead = tot <= 0
        if dead.any():
            n_hk[u][dead] = hk[g[0]][dead]
            n_k33[u][dead] = k33[g[0]][dead]
            n_ss[u][dead] = ss[g[0]][dead]
            n_thick[u][dead] = np.maximum(thick[g].sum(axis=0)[dead], 0.0)

    sign = np.sign(np.asarray(_as3d(cMF.ibound, nlay, nrow, ncol)))
    neg = (sign < 0).any(axis=0)
    ibound_new = n_ib.astype(int)
    if neg.any():                                     # preserve constant heads
        for u in range(nnew):
            ibound_new[u][neg & (ibound_new[u] > 0)] = -1

    cMF.nlay = nnew
    cMF.ibound = ibound_new
    cMF.botm = n_botm
    cMF.thick = n_thick
    cMF.strt = n_strt
    cMF.hk_actual = n_hk
    cMF.ss_actual = n_ss
    cMF.sy_actual = n_sy
    # store a TRUE vertical K and mark layvka=0 so downstream code stops
    # treating vka as a ratio
    cMF.vka_actual = n_k33
    cMF.layvka = [0] * nnew
    cMF.laytyp = [int(np.asarray(cMF.laytyp).ravel()[g[0]]) for g in groups]
    cMF.layavg = [int(np.asarray(cMF.layavg).ravel()[g[0]]) for g in groups]
    cMF.laywet = [0] * nnew
    cMF.laycbd = [0] * nnew
    cMF.Mlay = list(range(1, nnew + 1))
    cMF.Mnlay = nnew
    cMF.h_plt = [1] * nnew
    cMF.h_lbl = ['U%d' % (u + 1) for u in range(nnew)]
    if getattr(cMF, 'iuzfbnd', None) is not None:
        iu = np.asarray(cMF.iuzfbnd)
        if iu.ndim == 3:
            cMF.iuzfbnd = iu[[g[0] for g in groups]]
    # DRN / GHB cell lists are keyed by layer index -> remap onto the units
    lut = {}
    for u, g in enumerate(groups):
        for k in g:
            lut[k] = u
    for attr in ('layer_row_column_elevation_cond', 'layer_row_column_head_cond'):
        d = getattr(cMF, attr, None)
        if isinstance(d, dict) and 0 in d:
            seen, remapped = set(), []
            for rec in d[0]:
                l, i, j = int(rec[0]), int(rec[1]), int(rec[2])
                key = (lut.get(l, 0), i, j)
                if key in seen:            # collapsed duplicates: keep the first
                    continue
                seen.add(key)
                remapped.append([lut.get(l, 0), i, j] + list(rec[3:]))
            d[0] = remapped

    report = {'aggregated': True, 'nlay_before': nlay, 'nlay': nnew,
              'groups': groups,
              'active_cells': int((ibound_new != 0).sum())}
    if verbose:
        print('layer aggregation: %d numerical layers -> %d hydrogeological '
              'unit(s) %s' % (nlay, nnew, [[k + 1 for k in g] for g in groups]))
        for u, g in enumerate(groups):
            a = int((ibound_new[u] != 0).sum())
            print('   unit %d: layers %s, %d active cells, k %.3g..%.3g, '
                  'k33 %.3g..%.3g m/d'
                  % (u + 1, [k + 1 for k in g], a,
                     n_hk[u][ibound_new[u] != 0].min() if a else 0,
                     n_hk[u][ibound_new[u] != 0].max() if a else 0,
                     n_k33[u][ibound_new[u] != 0].min() if a else 0,
                     n_k33[u][ibound_new[u] != 0].max() if a else 0))
    return report
