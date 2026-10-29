# -*- coding: utf-8 -*-
"""Maps of ONE day of a finished run (user, 2026-10-29): the Results panel
asks for a day and a flux and draws, for that day alone, the map the run
draws as a time mean (GWmap_*, MMmap_*) -- on the linear ramp and on the
percentile ramp, like the run's own maps.

The run keeps every day of the aquifer budget in MF6's cell budget and
heads files, and the soil column's in the coupled h5. All three sit in the
run's MF6 workspace, and the NEXT run there overwrites them. So a day can
be drawn only while they are still this run's: their files must date from
the run, from its results folder's stamp to its last figure -- the rule
marmites_runinfo applies to a saved state.

One map is the panel's own: the wetness W (user, 2026-10-29), the soil's
saturation with the open water counted as saturated, kept for comparing
the model with imagery.

Streamlit-free: the page only asks and shows. The figures go into
``<out>/_output/day_maps`` as ``daymap_<GW|MM>_<flux>_<YYYY-MM-DD>.png``
and the percentile twin ``..._pct.png``, kept with the run's results.
"""

import glob
import os
import time
from datetime import date, datetime

import numpy as np

__all__ = ['FLUXES', 'DAY_DIR', 'label', 'run_dates', 'model_name',
           'outputs_of', 'day_values', 'draw', 'drawn']

DAY_DIR = 'day_maps'                 # under <out>/_output
_H5 = '_coupled_lagged.h5'           # MF6Coupler.RESULTS_H5
_SLACK_S = 300.0                     # as marmites_runinfo._SLACK_S
_FOLDERS = ('_input', '_output', 'figures_nwt_comparison')

# What can be drawn: (key, family, label, unit, colour map). GW = one panel
# per layer from MF6's outputs, the same terms as the run's GWmap_* (Rn and
# the others of RECHARGE_MAPS, _AQ_MAPS, dSg, the head); MM = the soil
# column, one panel, from the coupled h5. A signed field (Rn, Re, FLF, dSg)
# goes onto plotLAYER's diverging ramp whatever its colour map.
FLUXES = (
    ('Rn', 'GW', 'net recharge (Rg + Rsw - Exf - ETg)', 'mm/d', 'Blues'),
    ('Re', 'GW', 'effective recharge (Rg + Rsw - Exf)', 'mm/d', 'Blues'),
    ('Rg', 'GW', 'recharge to groundwater', 'mm/d', 'Blues'),
    ('Rsw', 'GW', 'recharge from the streams and ponds (SFR->GW + LAK->GW)',
     'mm/d', 'Blues'),
    ('EXFg', 'GW', 'exfiltration: seepage + GW->SFR + GW->LAK', 'mm/d',
     'Reds'),
    ('ETg', 'GW', 'groundwater ET (Eg + Tg)', 'mm/d', 'Reds'),
    ('Eg', 'GW', 'groundwater evaporation (EVT)', 'mm/d', 'Reds'),
    ('Tg', 'GW', 'groundwater transpiration (EVT)', 'mm/d', 'Reds'),
    ('DRN', 'GW', 'boundary drainage', 'mm/d', 'Reds'),
    ('FLF', 'GW', 'flow across the lower face (+ down, - up)', 'mm/d',
     'Blues'),
    ('dSg', 'GW', 'release from groundwater storage', 'mm/d', 'Blues'),
    ('head', 'GW', 'head', 'm', 'Blues'),
    # the wetness the user kept for the comparison with imagery (2026-10-29):
    # the soil's saturation, comparable across soil zones, the open water
    # (ponds, flowing reaches) counted as saturated
    ('W', 'MM', 'wetness (1 - f_ow) Se + f_ow: soil saturation, open water '
     'as 1', '-', 'Blues'),
    # the h5 keeps the column's storage per cell, not each horizon's: the
    # column's moisture is its storage over its thickness (the run's MMmap
    # theta is the horizons' plain mean, the same where they are as thick)
    ('theta', 'MM', 'soil moisture of the soil column (storage / thickness)',
     '% vol', 'Blues'),
    ('Ssoil', 'MM', 'soil water storage', 'mm', 'Blues'),
    ('Rp', 'MM', 'percolation', 'mm/d', 'Blues'),
    ('ETuzf', 'MM', 'unsaturated zone ET (UZF)', 'mm/d', 'Reds'),
    ('Ro', 'MM', 'runoff', 'mm/d', 'Reds'),
)
_BY_KEY = {f[0]: f for f in FLUXES}


def label(key):
    """The flux as the page lists it."""
    _k, fam, lbl, unit, _c = _BY_KEY[key]
    return '%s: %s, %s [%s]' % ('Aquifer' if fam == 'GW' else 'Soil column',
                                key, lbl, unit)


_TRUNK = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def _pp():
    import marmites_postprocess as pp
    return pp


def _on_path():
    import sys
    for p in (_TRUNK, os.path.join(_TRUNK, 'MARMITESutilities')):
        if p not in sys.path:
            sys.path.insert(0, p)


def _utils():
    """MARMITESutilities, whose remapped colour map plotLAYER's diverging
    ramp takes from the ``cMF`` it is given."""
    _on_path()
    import MARMITESutilities as MMutils
    return MMutils


# ------------------------------------------------------------------ #
# which run, which files, which days
# ------------------------------------------------------------------ #

def run_dates(out_dir):
    """The date of each stress period of the run, from a per-period table
    it wrote into ``_output`` (sp_run_timing.csv, else the outlet's flow),
    or None."""
    import pandas as pd
    for f in ('sp_run_timing.csv', 'outlet_streamflow.csv'):
        p = os.path.join(out_dir, '_output', f)
        if not os.path.exists(p):
            continue
        try:
            d = pd.read_csv(p, usecols=['date'])['date']
            return [x.date() for x in pd.to_datetime(d)]
        except Exception:                            # noqa: BLE001
            continue
    return None


def model_name(sim_ws):
    """The groundwater model's name, from the simulation's mfsim.nam."""
    try:
        with open(os.path.join(sim_ws, 'mfsim.nam'), encoding='utf-8',
                  errors='replace') as fh:
            lines = fh.read().splitlines()
    except OSError:
        return None
    inside = False
    for ln in lines:
        t = ln.split('#')[0].split()
        if not t:
            continue
        if t[0].upper() == 'BEGIN' and t[1:2] and t[1].upper() == 'MODELS':
            inside = True
        elif t[0].upper() == 'END':
            inside = False
        elif inside and t[0].upper().startswith('GWF') and len(t) >= 3:
            return t[2].lower()
    return None


def _window(out_dir):
    """(start, end) epoch seconds of the run that wrote ``out_dir``: its
    stamp (out_<YYYYMMDDHHMM>_..., the run's start) and its last result --
    the files of its folders, not the day maps added later."""
    try:
        start = time.mktime(time.strptime(
            os.path.basename(out_dir).split('_')[1], '%Y%m%d%H%M'))
    except (IndexError, ValueError):
        start = None
    end = None
    for sub in _FOLDERS:
        d = os.path.join(out_dir, sub)
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            p = os.path.join(d, f)
            if os.path.isfile(p):
                end = max(end or 0.0, os.path.getmtime(p))
    return start, end


def _when(t):
    return time.strftime('%d %b %Y %H:%M', time.localtime(t))


def outputs_of(out_dir, sim_ws=None):
    """Where the model outputs of the run that wrote ``out_dir`` are and
    whether they are still its own: ``{'sim_ws', 'name', 'dates', 'h5',
    'soil', 'why'}``, ``why`` '' when a day can be drawn, else the reason;
    ``soil`` the run's soil horizons ``(zone, slprop, smax, sr)``, which
    the wetness W needs.

    ``sim_ws``: the MF6 workspace; by default the one the run's resolved
    configuration names, by the driver's rule (``paths.ws``, else the
    mesh's workspace under WS_ROOT, the results folder's parent).
    """
    res = {'sim_ws': sim_ws, 'name': None, 'dates': None, 'h5': None,
           'soil': None, 'why': ''}
    cfg_p = os.path.join(out_dir, '_input', 'resolved_config.toml')
    cfg = None
    if os.path.exists(cfg_p):
        try:
            import marmites_config as mcfg
            cfg = mcfg.load_run_config(cfg_p)
        except Exception as exc:                     # noqa: BLE001
            if sim_ws is None:
                res['why'] = 'its configuration copy does not load (%s)' % exc
                return res
    if sim_ws is None:
        if cfg is None:
            res['why'] = 'it has no copy of its configuration'
            return res
        sim_ws = (str(cfg.paths.ws) if cfg.paths.ws else
                  mcfg.state_workspace(cfg, os.path.dirname(
                      os.path.abspath(out_dir))))
        res['sim_ws'] = sim_ws
    if cfg is not None:
        res['soil'] = [(int(h.zone), float(h.slprop), float(h.smax),
                        float(h.sr)) for h in cfg.soil.horizon]
        # for the comparison with an image (marmites_wetimage): the model's
        # CRS and its irrigated fields, a shapefile in DATA_ROOT/GIS
        if getattr(cfg.grid, 'crs_epsg', 0):
            res['crs'] = 'EPSG:%d' % int(cfg.grid.crs_epsg)
        if cfg.surface.irrigation and cfg.surface.irr_zones.layer:
            res['irr_layer'] = str(cfg.surface.irr_zones.layer)
    name = model_name(sim_ws)
    res['name'] = name
    if not name:
        res['why'] = 'there is no MF6 simulation in %s' % sim_ws
        return res
    dates = run_dates(out_dir)
    res['dates'] = dates
    if not dates:
        res['why'] = ('it did not record the date of its stress periods '
                      '(sp_run_timing.csv, written since 2026-10-26)')
        return res
    files = [os.path.join(sim_ws, '%s.%s' % (name, x)) for x in ('cbc', 'hds')]
    missing = [os.path.basename(p) for p in files if not os.path.exists(p)]
    if missing:
        res['why'] = '%s not in %s' % (', '.join(missing), sim_ws)
        return res
    h5 = os.path.join(sim_ws, _H5)
    res['h5'] = h5 if os.path.exists(h5) else None
    start, end = _window(out_dir)
    t = max(os.path.getmtime(p) for p in files)
    if start is not None and t < start - 60.0:
        res['why'] = ('the model outputs in %s are older than this run (%s): '
                      'it did not run the model (a figures-only re-draw?)'
                      % (sim_ws, _when(t)))
    elif end is not None and t > end + _SLACK_S:
        res['why'] = ('the model outputs in %s were written after this run '
                      'ended (%s), by a later run: run this configuration '
                      'again to draw its days' % (sim_ws, _when(t)))
    return res


# ------------------------------------------------------------------ #
# the values of one day
# ------------------------------------------------------------------ #

def _grid(sim_ws, name):
    """The cells as plotLAYER draws them: polygons in real coordinates,
    areas, centres, the inactive mask (nlay, ncpl) and the column count
    that turns an (i, j) into a cell number."""
    import flopy
    grb = _pp()._grb_path(sim_ws, name)
    mg = flopy.mf6.utils.MfGrdFile(grb, verbose=False).modelgrid
    nlay = int(mg.nlay)
    if getattr(mg, 'nrow', None):                  # structured
        nrow, ncol = int(mg.nrow), int(mg.ncol)
        polys = [np.asarray(mg.get_cell_vertices(i, j), float)
                 for i in range(nrow) for j in range(ncol)]
    else:                                          # a mesh: (ncpl, 1)
        n = int(np.atleast_1d(mg.ncpl)[0])
        polys = [np.asarray(mg.get_cell_vertices(c), float)
                 for c in range(n)]
        nrow, ncol = n, 1
    ncpl = len(polys)
    area = np.array([0.5 * abs(np.dot(p[:, 0], np.roll(p[:, 1], 1))
                               - np.dot(p[:, 1], np.roll(p[:, 0], 1)))
                     for p in polys])
    cen = np.array([p[:, :2].mean(axis=0) for p in polys])
    idom = mg.idomain
    off = (np.zeros((nlay, ncpl), bool) if idom is None
           else np.asarray(idom).reshape(nlay, ncpl) <= 0)
    return {'polys': polys, 'area': area, 'xc': cen[:, 0], 'yc': cen[:, 1],
            'off': off, 'nlay': nlay, 'ncpl': ncpl, 'nrow': nrow,
            'ncol': ncol, 'grb': grb}


def _nodes(ij, g):
    """The cell numbers of MARMITES' (i, j) pairs: (row, column) on a
    structured grid; on a mesh (ncpl, 1) the cell is i -- or j, the
    coupler's (0, cell)."""
    ij = np.asarray(ij, dtype=int).reshape(-1, 2)
    if g['ncol'] == 1 and ij.size and not ij[:, 1].any():
        return ij[:, 0]
    if g['ncol'] == 1:
        return ij[:, 1]
    return ij[:, 0] * g['ncol'] + ij[:, 1]


def _budget(sim_ws, name, k, nper, g):
    """Every aquifer record of stress period ``k`` (0-based, of the run's
    ``nper``), m3/d as (nlay, ncpl), each period the time-weighted mean of
    its steps -- the map pass of the run, on one period."""
    import flopy
    pp = _pp()
    nlay, ncpl = g['nlay'], g['ncpl']
    cbc = flopy.utils.CellBudgetFile(os.path.join(sim_ws, '%s.cbc' % name))
    try:
        kk = cbc.get_kstpkper()
        grp = pp.period_steps(cbc.get_times(), kk,
                              steady=pp.steady_first(sim_ws, name))[-nper:][k]
        recs = set(r.strip() for r in cbc.get_unique_record_names(decode=True))
        out = {}
        for key, text, pak2 in pp._AQ_RECORDS:
            if text not in recs:
                continue
            acc, a_in, a_out = (np.zeros((nlay, ncpl)) for _ in range(3))
            got = False
            for r, w in grp:
                try:
                    d = cbc.get_data(text=text, paknam2=pak2, kstpkper=kk[r],
                                     full3D=True)
                except ValueError:       # a package this run did not have
                    d = None
                if not d:
                    continue
                a = np.ma.filled(np.asarray(d[0], float), 0.0).reshape(
                    nlay, ncpl)
                acc += w * a
                a_in += w * np.maximum(a, 0.0)
                a_out += w * np.minimum(a, 0.0)
                got = True
            if got:
                out[key] = acc
                if key in pp._AQ_SPLIT:
                    out[key + '_IN'], out[key + '_OUT'] = a_in, a_out
        if nlay > 1 and 'FLOW-JA-FACE' in recs:
            src, pos, nodes = pp._ja_down_index(g['grb'], nlay, g['nrow'],
                                                g['ncol'])
            acc = np.zeros(nodes)
            for r, w in grp:
                fj = np.asarray(cbc.get_data(text='FLOW-JA-FACE',
                                             kstpkper=kk[r])[0], float).ravel()
                a = np.zeros(nodes)
                a[src] = -fj[pos]                    # flopy's sign convention
                acc += w * a
            out['FLF'] = acc.reshape(nlay, ncpl)
        return out
    finally:
        cbc.close()


def _head(sim_ws, name, k, nper, g):
    import flopy
    pp = _pp()
    hds = flopy.utils.HeadFile(os.path.join(sim_ws, '%s.hds' % name))
    try:
        kk = pp.period_end_kk(hds, steady=pp.steady_first(sim_ws, name))
        h = np.asarray(hds.get_data(kstpkper=kk[-nper:][k]), float)
    finally:
        hds.close()
    h = h.reshape(g['nlay'], g['ncpl'])
    return np.where(np.abs(h) > 1e29, np.nan, h)


def _thickness(h):
    """Each soil column's horizon thicknesses [mm], (ncell, nsl), from the
    run's means: a horizon's storage over its moisture (a fraction), both
    means of the same constant ratio."""
    _on_path()
    from marmites_indices import INDEX_MM_SOIL as IXS
    ms = np.asarray(h['wb_map_soil'][...], float)
    s, f = ms[:, :, IXS['iSsoil']], ms[:, :, IXS['iSsoil_pc_s']]
    return np.where(f > 0, s / np.where(f > 0, f, 1.0), 0.0)


def _column_storages(tl, soil):
    """Each column's storage at saturation and at residual moisture [mm]:
    its horizons' smax and sr times their thicknesses, NaN for a column of
    no zone. The h5 does not say which zone a column is; MMsoil splits it
    by the zone's slprop, so the split tells -- unless two zones split
    alike and differ in their properties, which is refused."""
    zones = {}
    for z, slp, smax, sr in soil:
        zones.setdefault(int(z), []).append((float(slp), float(smax), float(sr)))
    nsl = tl.shape[1]
    t = tl.sum(axis=1)
    share = tl / np.where(t > 0, t, 1.0)[:, None]
    sig = {z: np.r_[[h[0] for h in hz], np.zeros(max(nsl - len(hz), 0))][:nsl]
           for z, hz in zones.items()}
    keys = sorted(sig)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            if np.allclose(sig[a], sig[b], atol=1e-3) and zones[a] != zones[b]:
                raise ValueError('soil zones %d and %d split their columns '
                                 'alike, so the wetness cannot tell them '
                                 'apart' % (a, b))
    ssat, sres = np.full(t.size, np.nan), np.full(t.size, np.nan)
    for z in keys:
        m = (t > 0) & np.all(np.abs(share - sig[z]) < 1e-3, axis=1)
        ssat[m] = sum(h[1] * tl[m, l] for l, h in enumerate(zones[z]))
        sres[m] = sum(h[2] * tl[m, l] for l, h in enumerate(zones[z]))
    return ssat, sres


def _block(path, block):
    """The rows of an MF6 input file's BEGIN <block> ... END <block>."""
    rows, inside = [], False
    with open(path, encoding='utf-8', errors='replace') as fh:
        for ln in fh:
            t = ln.split('#')[0].split()
            if len(t) >= 2 and t[0].upper() == 'BEGIN' \
                    and t[1].lower() == block:
                inside = True
            elif len(t) >= 2 and t[0].upper() == 'END' \
                    and t[1].lower() == block:
                break
            elif inside and t:
                rows.append(t)
    return rows


def _stages(path, k, nper, steady):
    """An SFR or LAK stage file's values at the end of period ``k``."""
    import flopy
    hf = flopy.utils.HeadFile(path, text='STAGE', precision='double')
    try:
        kk = _pp().period_end_kk(hf, steady=steady)[-nper:][k]
        return np.asarray(hf.get_data(kstpkper=kk), float).ravel()
    finally:
        hf.close()


def _cellid(t, g):
    """The cell number of an MF6 cellid's tokens (layer first): (layer,
    cell) on a mesh, (layer, row, column) on a structured grid."""
    if g['ncol'] == 1:
        return int(t[1]) - 1
    return (int(t[1]) - 1) * g['ncol'] + int(t[2]) - 1


def pond_cells(sim_ws, name, g):
    """``{lake: [cells]}``: the cells each pond (LAK) sits on."""
    lak = os.path.join(sim_ws, '%s.lak' % name)
    out = {}
    if os.path.exists(lak):
        nci = 2 if g['ncol'] == 1 else 3
        for t in _block(lak, 'connectiondata'):
            out.setdefault(int(t[0]) - 1, []).append(_cellid(t[2:2 + nci], g))
    return out


def _open_water(sim_ws, name, k, nper, g):
    """The share of each cell under open water on period ``k``: a pond's
    cells while its stage is above its bottom (the first stage of its
    table), and a reach's channel, length x width, while it carries
    water -- read from the run's MF6 files."""
    f = np.zeros(g['ncpl'])
    nci = 2 if g['ncol'] == 1 else 3     # (layer, cell) or (layer, row, col)

    def node(t):
        return _cellid(t, g)
    steady = _pp().steady_first(sim_ws, name)
    sfr = os.path.join(sim_ws, '%s.sfr' % name)
    if os.path.exists(sfr) and os.path.exists(sfr + '.stage'):
        stage = _stages(sfr + '.stage', k, nper, steady)
        for t in _block(sfr, 'packagedata'):
            if t[1].upper() == 'NONE':
                continue
            r, c = int(t[0]) - 1, node(t[1:1 + nci])
            rlen, rwid, rtp = (float(t[1 + nci]), float(t[2 + nci]),
                               float(t[4 + nci]))
            if r < stage.size and stage[r] - rtp > 1e-3:
                f[c] += rlen * rwid / g['area'][c]
    lak = os.path.join(sim_ws, '%s.lak' % name)
    if os.path.exists(lak) and os.path.exists(lak + '.stage'):
        stage = _stages(lak + '.stage', k, nper, steady)
        bottom = {}
        for t in _block(lak, 'tables'):
            tab = os.path.join(sim_ws, t[-1])
            rows = _block(tab, 'table') if os.path.exists(tab) else []
            if rows:
                bottom[int(t[0]) - 1] = float(rows[0][0])
        for n, cells in pond_cells(sim_ws, name, g).items():
            if n < stage.size and stage[n] - bottom.get(n, -np.inf) > 1e-3:
                f[cells] = 1.0
    return np.clip(f, 0.0, 1.0)


def saturation(h5, k, g, soil):
    """Se of each soil column on period ``k``, on the grid's cells: (S -
    S_r) / (S_sat - S_r), clipped to 0..1; NaN off the soil columns."""
    import h5py
    if not soil:
        raise ValueError("the wetness needs the run's soil horizons, from "
                         "its configuration copy")
    with h5py.File(h5, 'r') as h:
        node = _nodes(h['cell_ij'][...], g)
        ss = np.asarray(h['ssoil'][k, :], float)
        tl = _thickness(h)
    ssat, sres = _column_storages(tl, soil)
    se = np.full(g['ncpl'], np.nan)
    se[node] = np.clip((ss - sres) / np.where(ssat > sres, ssat - sres,
                                              np.nan), 0.0, 1.0)
    return se


def _wetness(sim_ws, name, k, nper, g, h5, soil):
    """W = (1 - f_ow) Se + f_ow (user, 2026-10-29): the effective
    saturation of the soil column, Se = (S - S_r) / (S_sat - S_r), with the
    cell's open water counted as saturated; 0 at residual moisture, 1
    saturated or under water. Dimensionless, like a satellite soil water
    index; the field-capacity version saturated over most of La Mata in
    June and was dropped."""
    se = saturation(h5, k, g, soil)
    f = _open_water(sim_ws, name, k, nper, g)
    w = np.where(np.isfinite(se), (1.0 - f) * se + f,
                 np.where(f > 0, f, np.nan))
    return w[None, :]


def _soil(h5, key, k, g):
    """One soil-column field of period ``k`` on the grid's cells, (1, ncpl),
    NaN off the soil columns."""
    import h5py
    with h5py.File(h5, 'r') as h:
        node = _nodes(h['cell_ij'][...], g)
        if key == 'theta':
            ss = np.asarray(h['ssoil'][k, :], float)
            tl = _thickness(h)
            t = tl.sum(axis=1)
            v = np.where(t > 0, ss / np.where(t > 0, t, 1.0) * 100.0, np.nan)
        else:
            ds, fac = {'Ssoil': ('ssoil', 1.0), 'Rp': ('perc', 1000.0),
                       'ETuzf': ('etuzf', 1.0), 'Ro': ('runoff', 1.0)}[key]
            if ds not in h:
                raise ValueError('this run did not keep %s per cell' % ds)
            v = np.asarray(h[ds][k, :], float) * fac
    out = np.full(g['ncpl'], np.nan)
    out[node] = v
    return out[None, :]


def day_values(sim_ws, name, k, nper, key, g=None, h5=None, soil=None):
    """The map of ``key`` on stress period ``k`` (0-based of ``nper``):
    (nlay, ncpl) for an aquifer term, (1, ncpl) for the soil column, in
    the unit FLUXES gives, NaN where there is nothing to show. ``soil``:
    the run's horizons ``(zone, slprop, smax, sr)``, for W."""
    g = g or _grid(sim_ws, name)
    fam = _BY_KEY[key][1]
    if fam == 'MM':
        if not h5:
            raise ValueError('the coupled results (%s) are not in %s'
                             % (_H5, sim_ws))
        if key == 'W':
            return _wetness(sim_ws, name, k, nper, g, h5, soil)
        return _soil(h5, key, k, g)
    if key == 'head':
        v = _head(sim_ws, name, k, nper, g)
    else:
        pp = _pp()
        m = _budget(sim_ws, name, k, nper, g)

        def rec(r):
            return m[r] if r in m else np.zeros((g['nlay'], g['ncpl']))
        if key in ('Rn', 'Re', 'Rsw', 'EXFg'):
            rm = pp.recharge_maps(m, g['nlay'])
            if key not in rm:
                raise ValueError('the run has no UZF recharge')
            v = rm[key]
        else:
            v = {'Rg': lambda: rec('UZF-GWRCH'),
                 'ETg': lambda: -(rec('WEL') + rec('EVT_EG') + rec('EVT_TG')),
                 'Eg': lambda: -rec('EVT_EG'),
                 'Tg': lambda: -rec('EVT_TG'),
                 'DRN': lambda: -rec('DRN'),
                 'FLF': lambda: rec('FLF'),
                 'dSg': lambda: rec('STO-SS') + rec('STO-SY')}[key]()
        v = np.asarray(v, float) * 1000.0 / g['area'][None, :]  # m3/d -> mm/d
    off = g['off']
    if key == 'FLF':
        # the bottom layer has no lower face: nothing to show there
        v, off = v[:-1], off[:-1]
    return np.where(off, np.nan, v)


# ------------------------------------------------------------------ #
# the figure
# ------------------------------------------------------------------ #

def _fname(fam, key, day, suffix, page=None):
    return 'daymap_%s_%s_%s%s%s.png' % (fam, key, day.isoformat(),
                                        '' if page is None else
                                        '_page%d' % page, suffix)


def draw(out_dir, day, key, info=None, verbose=False):
    """Draw ``key`` on ``day`` (a date) for the run that wrote ``out_dir``,
    on the linear ramp and, for a field that is not constant, the
    percentile one; returns the PNGs written in ``<out>/_output/day_maps``.
    Raises ValueError with the reason when the day cannot be drawn."""
    import matplotlib
    matplotlib.use('agg')
    from types import SimpleNamespace
    pp = _pp()
    info = info or outputs_of(out_dir)
    if info.get('why'):
        raise ValueError(info['why'])
    if isinstance(day, datetime):
        day = day.date()
    dates = list(info['dates'])
    if day not in dates:
        raise ValueError('%s is not a day of this run (%s to %s)'
                         % (day, dates[0], dates[-1]))
    k, nper = dates.index(day), len(dates)
    sim_ws, name = info['sim_ws'], info['name']
    _key, fam, lbl, unit, cname = _BY_KEY[key]
    g = _grid(sim_ws, name)
    V = day_values(sim_ws, name, k, nper, key, g=g, h5=info.get('h5'),
                   soil=info.get('soil'))
    nl = V.shape[0]
    hnoflo = 9999.999
    off = ~np.isfinite(V)
    vals = V[~off]
    if not vals.size:
        raise ValueError('nothing to draw: %s is empty on %s' % (key, day))
    lo, hi = pp._vrange(vals)
    pts = None
    if info.get('h5'):
        import h5py
        with h5py.File(info['h5'], 'r') as h:
            if 'obs_names' in h and 'obs_ij' in h:
                nm = [n.decode() if isinstance(n, bytes) else str(n)
                      for n in h['obs_names'][...]]
                c = _nodes(h['obs_ij'][...], g)
                pts = [nm, [float(g['yc'][i]) for i in c],
                       [float(g['xc'][i]) for i in c], [0] * len(c)]
    contours = None
    if key == 'head' and hi > lo:
        contours = pp._head_levels(lo, hi)
    MMplot = pp._mmplot()
    geo = SimpleNamespace(cUTIL=_utils().clsUTILITIES)  # the diverging remap
    dest = os.path.join(out_dir, '_output', DAY_DIR)
    os.makedirs(dest, exist_ok=True)
    dn = matplotlib.dates.date2num(datetime(day.year, day.month, day.day))
    written = []
    modes = (('linspace', ''),) + ((('percentile', '_pct'),)
                                   if vals.min() < vals.max() else ())
    for itype, suffix in modes:
        title = '%smap_%s_%s%s' % (fam, key, day.isoformat(), suffix)
        for old in glob.glob(os.path.join(dest, '_sp_plt_%s_day*.png'
                                          % title)):
            os.remove(old)
        MMplot.plotLAYER(
            days=[k], str_per=[k + 1], Date=[dn],
            JD=[day.timetuple().tm_yday], ncol=1, nrow=g['ncpl'], nlay=nl,
            nplot=nl, V=np.where(off, hnoflo, V)[None, :, :, None],
            cmap=matplotlib.colormaps[cname],
            CBlabel='%s [%s]%s' % (lbl, unit, ', percentile classes'
                                   if suffix else ''), msg='',
            plt_title=title, MM_ws=dest, interval_type=itype,
            interval_num=5, Vmax=[hi], Vmin=[lo],
            fmt='%5.1f' if unit in ('% vol', 'm') else '%5.2f',
            points=pts, mask=off[:, :, None], hnoflo=hnoflo, cMF=geo,
            polys=g['polys'], contours=contours is not None,
            contour_levels=contours, contour_fmt='%.0f')
        pages = sorted(glob.glob(os.path.join(
            dest, '_sp_plt_%s_day%05d_*_*.png' % (title, k + 1))))
        for n, p in enumerate(pages):
            to = os.path.join(dest, _fname(fam, key, day, suffix,
                                           None if len(pages) == 1 else n + 1))
            os.replace(p, to)
            written.append(to)
    if verbose:
        print('day map %s %s: %s' % (key, day, ', '.join(
            os.path.basename(p) for p in written)))
    return written


def drawn(out_dir):
    """The day maps drawn for this run, ``[(sub, fname, path)]``: the
    latest first, each map before its percentile twin -- and the
    comparisons with an image (``wetimage_``, marmites_wetimage)."""
    d = os.path.join(out_dir, '_output', DAY_DIR)
    if not os.path.isdir(d):
        return []
    fs = [f for f in os.listdir(d)
          if f.startswith(('daymap_', 'wetimage_')) and f.endswith('.png')]
    groups = {}
    for f in fs:
        base = f[:-len('_pct.png')] + '.png' if f.endswith('_pct.png') else f
        groups.setdefault(base, []).append(f)
    sub = '_output/%s' % DAY_DIR

    def newest(base):
        return max(os.path.getmtime(os.path.join(d, f))
                   for f in groups[base])
    out = []
    for base in sorted(groups, key=newest, reverse=True):
        for f in sorted(groups[base], key=lambda f: f.endswith('_pct.png')):
            out.append((sub, f, os.path.join(d, f)))
    return out
