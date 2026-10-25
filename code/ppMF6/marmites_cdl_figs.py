# -*- coding: utf-8 -*-
"""The CdL model's input figures, for any MARMITES model (user, 2026-10-25),
plus a cross-section pair per pond -- drawn by the input stage
(marmites_postprocess.run_preproc, [postproc] input_maps) into
<run>/_input/:

  IN_101_aquifer_params.png   horizontal K per MF6 layer, active cells only,
                              one shared (log) ramp
  IN_102_layer_surfaces.png   top / bottom / thickness of each MF6 layer
  IN_103_ibound_uzf.png       idomain per layer + the layer of each UZF land
                              object; the ponds outlined and numbered
  IN_104_land_cover_extdp.png the vegetation type covering most of each cell
                              + UZF's extinction depth
  IN_105_uzf_params.png       UZF's vks, thts, thti, thtr, eps (land objects)
  IN_106_forcing.png          monthly rainfall and open-water evaporation
  IN_107_ponds_zoom_p<n>.png  each pond: the grid, its SFR and LAK cells, its
                              footprint and the stream (6 per page)
  IN_108_pond_sections_p<n>.png  each pond, two sections ACROSS the reach --
                              through its LAK cell and through the nearest
                              SFR cell: land surface, soil column, MF6 layers,
                              the pond (bottom, rim, initial stage), the
                              streambed, the initial water table (3 ponds a
                              page)

Everything model-side is read from the WRITTEN simulation (as the model map
is); the land surface, the forcing and the vegetation cover come from the
run's own cMF / ctx, the vegetation names and root depths from the
configuration.
"""

import os

import numpy as np

PONDS_PER_ZOOM_PAGE = 6
PONDS_PER_SECTION_PAGE = 3
# the LAK and SFR cells look the same in the pond zooms and the sections
LAK_STYLE = dict(facecolor='orange', alpha=0.75, edgecolor='crimson', lw=1.2)
SFR_STYLE = dict(facecolor='deepskyblue', alpha=0.55, edgecolor='none')


# --------------------------------------------------------------------- #
# reading
# --------------------------------------------------------------------- #

def load_model(sim_ws, name):
    """(gwf, mg) of the written simulation, or (None, None)."""
    import flopy
    if not os.path.exists(os.path.join(sim_ws, 'mfsim.nam')):
        return None, None
    sim = flopy.mf6.MFSimulation.load(
        sim_ws=sim_ws, verbosity_level=0,
        load_only=['dis', 'disv', 'npf', 'sto', 'ic', 'uzf', 'lak', 'sfr'])
    gwf = (sim.get_model(name) if name and name in sim.model_names
           else sim.get_model())
    return gwf, gwf.modelgrid


def _pkg(gwf, key):
    for p in gwf.packagelist:
        if p.package_type.lower() == key:
            return p
    return None


def _node(mg):
    vertex = mg.grid_type == 'vertex'
    ncol = 1 if vertex else int(mg.ncol)

    def node(cid):
        c = [int(v) for v in cid]
        return c[1] if vertex else c[1] * ncol + c[2]
    return node


def _layered(a, mg):
    """An (nlay, ncell) view of a layered array."""
    return np.asarray(a, dtype=float).reshape(int(mg.nlay), -1)


def _polys(mg):
    n = int(np.asarray(mg.xcellcenters).size)
    return [np.asarray(mg.get_cell_vertices(i), dtype=float)[:, :2]
            for i in range(n)]


def _boundary(ds_ws):
    from marmites_postprocess import _geojson_polygons
    rings = []
    for _props, rr in _geojson_polygons(os.path.join(ds_ws or '',
                                                     'inputBOUNDARY.geojson'))[0]:
        rings += rr
    return rings


def uzf_land(gwf, mg):
    """Per cell: the UZF land object's layer (1-based), vks, thtr, thts,
    thti, eps and extdp (period 1) -- NaN where none."""
    node = _node(mg)
    n = int(np.asarray(mg.xcellcenters).size)
    out = {k: np.full(n, np.nan) for k in ('layer', 'vks', 'thtr', 'thts',
                                           'thti', 'eps', 'extdp')}
    u = _pkg(gwf, 'uzf')
    if u is None:
        return out
    pk = u.packagedata.get_data()
    land = {}
    for r in pk:
        if int(r['landflag']) == 1:
            c = node(r['cellid'])
            land[int(r['ifno'])] = c
            out['layer'][c] = int(r['cellid'][0]) + 1
            for k in ('vks', 'thtr', 'thts', 'thti', 'eps'):
                out[k][c] = float(r[k])
    pd = u.perioddata.get_data(0) if u.perioddata.has_data() else None
    if pd is not None:
        for r in pd:
            c = land.get(int(r['ifno']))
            if c is not None:
                try:
                    out['extdp'][c] = float(r['extdp'])
                except (TypeError, ValueError):
                    pass
    return out


def lakes(gwf, mg, sim_ws):
    """Each lake: name (pond<fid>), its connection cell(s), initial stage,
    rim (outlet invert, else the initial stage) and bottom (the first stage
    of its table)."""
    node = _node(mg)
    lk = _pkg(gwf, 'lak')
    if lk is None:
        return []
    out = {}
    for r in lk.packagedata.get_data():
        out[int(r['ifno'])] = {'name': str(r['boundname']) if 'boundname'
                               in r.dtype.names else 'lake%d' % (r['ifno'] + 1),
                               'strt': float(r['strt']), 'cells': [],
                               'rim': float(r['strt']), 'bottom': np.nan}
    for r in lk.connectiondata.get_data():
        out[int(r['ifno'])]['cells'].append(node(r['cellid']))
    if lk.outlets.has_data():
        for r in lk.outlets.get_data():
            L = int(r['lakein'])
            if L in out:
                out[L]['rim'] = float(r['invert'])
    if lk.tables.has_data():
        for r in lk.tables.get_data():
            L = int(r['ifno'])
            fn = os.path.join(sim_ws, str(r['tab6_filename']))
            if L in out and os.path.exists(fn):
                with open(fn) as fh:
                    lines = fh.read().splitlines()
                for i, ln in enumerate(lines):
                    if ln.strip().upper().startswith('BEGIN TABLE'):
                        out[L]['bottom'] = float(lines[i + 1].split()[0])
                        break
    return [out[k] for k in sorted(out)]


def sfr_reaches(gwf, mg):
    """cell -> (reach number, rtp, rbth, rwid) of the first reach in it."""
    node = _node(mg)
    s = _pkg(gwf, 'sfr')
    out = {}
    if s is None:
        return out
    for r in s.packagedata.get_data():
        c = node(r['cellid'])
        if c not in out:
            out[c] = (int(r['ifno']) + 1, float(r['rtp']), float(r['rbth']),
                      float(r['rwid']))
    return out


# --------------------------------------------------------------------- #
# drawing helpers
# --------------------------------------------------------------------- #

def _map(ax, polys, values, rings, cmap='viridis', norm=None, title='',
         vmin=None, vmax=None):
    """A cell-valued map: NaN cells light grey."""
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Polygon as MplPoly
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(v)
    if (~ok).any():
        ax.add_collection(PatchCollection(
            [MplPoly(p) for p, o in zip(polys, ok) if not o],
            facecolor='0.88', edgecolor='none'))
    pc = PatchCollection([MplPoly(p) for p, o in zip(polys, ok) if o],
                         cmap=cmap, norm=norm, edgecolor='none')
    pc.set_array(v[ok])
    if norm is None and vmin is not None:
        pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    for r in rings:
        ax.plot(r[:, 0], r[:, 1], color='k', lw=0.8)
    allxy = np.vstack(polys)
    ax.set_xlim(allxy[:, 0].min(), allxy[:, 0].max())
    ax.set_ylim(allxy[:, 1].min(), allxy[:, 1].max())
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.ticklabel_format(useOffset=False, style='plain')
    return pc


def _log_ticks(cb, lo, hi):
    """Plain 1-2-5 ticks on a log colour bar (one decade gets one tick
    otherwise)."""
    t = [m * 10.0 ** e for e in range(int(np.floor(np.log10(lo))),
                                      int(np.ceil(np.log10(hi))) + 1)
         for m in (1, 2, 5)]
    t = [v for v in t if lo <= v <= hi]
    cb.set_ticks(t)
    cb.set_ticklabels(['%g' % v for v in t])
    cb.ax.minorticks_off()


def _save(fig, out, fn, written):
    import matplotlib.pyplot as plt
    path = os.path.join(out, fn)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    written.append(path)


# --------------------------------------------------------------------- #
# the figures
# --------------------------------------------------------------------- #

def fig_aquifer_params(out, gwf, mg, polys, rings, written):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    npf = _pkg(gwf, 'npf')
    k = _layered(npf.k.array, mg)
    idm = _layered(mg.idomain, mg) > 0
    act = np.where(idm, k, np.nan)
    pos = act[np.isfinite(act) & (act > 0)]
    norm = LogNorm(vmin=pos.min(), vmax=pos.max())
    nlay = int(mg.nlay)
    fig, axes = plt.subplots(1, nlay, figsize=(6 * nlay, 6.5), squeeze=False,
                             layout='constrained')
    pc = None
    for L in range(nlay):
        a = act[L]
        pc = _map(axes[0, L], polys, np.where(a > 0, a, np.nan), rings,
                  norm=norm, title='MF6 layer %d -- Kh, %d active cells '
                  '(median %.3g m/d)' % (L + 1, int(idm[L].sum()),
                                         np.nanmedian(a)))
    cb = fig.colorbar(pc, ax=axes.ravel().tolist(), shrink=0.85, pad=0.02)
    cb.set_label('Kh (m/d) -- shared log ramp; grey: inactive')
    _log_ticks(cb, norm.vmin, norm.vmax)
    fig.suptitle('Horizontal K per MF6 layer, shown only where the layer is '
                 'ACTIVE', fontsize=11)
    _save(fig, out, 'IN_101_aquifer_params.png', written)


def fig_layer_surfaces(out, mg, polys, rings, written):
    import matplotlib.pyplot as plt
    nlay = int(mg.nlay)
    top = np.asarray(mg.top, dtype=float).ravel()
    bot = _layered(mg.botm, mg)
    idm = _layered(mg.idomain, mg) > 0
    fig, axes = plt.subplots(nlay, 3, figsize=(16, 5.2 * nlay), squeeze=False)
    for L in range(nlay):
        t = top if L == 0 else bot[L - 1]
        for c, (arr, ttl, cm, un) in enumerate([
                (t, 'TOP', 'viridis', 'm a.s.l.'),
                (bot[L], 'BOTTOM', 'viridis', 'm a.s.l.'),
                (t - bot[L], 'THICKNESS', 'YlGnBu', 'm')]):
            pc = _map(axes[L, c], polys, np.where(idm[L], arr, np.nan), rings,
                      cmap=cm, title='MF6 layer %d -- %s (%s)' % (L + 1, ttl,
                                                                  un))
            fig.colorbar(pc, ax=axes[L, c], shrink=0.7)
    fig.suptitle('MF6 layers: top / bottom / thickness (layer 1\'s top is the '
                 'base of the soil column); grey: inactive')
    fig.tight_layout()
    _save(fig, out, 'IN_102_layer_surfaces.png', written)


def _pond_rings(ds_ws):
    from marmites_postprocess import _geojson_polygons
    out = []
    for props, rr in _geojson_polygons(os.path.join(ds_ws or '',
                                                    'inputPONDS.geojson'))[0]:
        fid = props.get('fid', props.get('id', ''))
        for r in rr:
            out.append((fid, r))
    return out


def fig_ibound_uzf(out, gwf, mg, polys, rings, ds_ws, written):
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    nlay = int(mg.nlay)
    idm = _layered(mg.idomain, mg)
    land = uzf_land(gwf, mg)['layer']
    cm = ListedColormap(['0.85', '#3182bd'])
    nm = BoundaryNorm([-0.5, 0.5, 1.5], cm.N)
    ponds = _pond_rings(ds_ws)
    panels = [(np.where(idm[L] > 0, 1.0, 0.0), 'idomain -- MF6 layer %d (%d '
               'active)' % (L + 1, int((idm[L] > 0).sum())), cm, nm)
              for L in range(nlay)]
    lcm = ListedColormap(['#9ecae1', '#08519c', '#ff7f00', '#e31a1c'][:nlay])
    lnm = BoundaryNorm(np.arange(0.5, nlay + 1.0), lcm.N)
    panels.append((land, 'UZF land objects -- the layer each sits in', lcm,
                   lnm))
    fig, axes = plt.subplots(1, len(panels), figsize=(5.4 * len(panels), 6),
                             squeeze=False, layout='constrained')
    for ax, (arr, ttl, c_, n_) in zip(axes[0], panels):
        pc = _map(ax, polys, arr, rings, cmap=c_, norm=n_, title=ttl)
        for fid, r in ponds:
            ax.plot(r[:, 0], r[:, 1], color='navy', lw=1.0)
            ax.annotate(str(fid), (r[:, 0].mean(), r[:, 1].mean()),
                        color='navy', fontsize=7, fontweight='bold',
                        ha='center', va='bottom', xytext=(0, 4),
                        textcoords='offset points')
        if ttl.startswith('UZF'):
            cb = fig.colorbar(pc, ax=ax, shrink=0.6,
                              ticks=np.arange(1, nlay + 1))
            cb.set_label('layer')
    fig.suptitle('idomain per MF6 layer (blue active, grey inactive) + the '
                 'layer of each UZF land object; ponds outlined and numbered',
                 fontsize=11)
    _save(fig, out, 'IN_103_ibound_uzf.png', written)


def fig_land_cover_extdp(out, gwf, mg, polys, rings, ctx, cfg, written):
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch
    n = len(polys)
    cover = np.asarray(getattr(ctx, 'gridVEGarea'), dtype=float)
    cover = cover.reshape(cover.shape[0], -1)[:, :n]          # (NVEG, ncell)
    nveg = cover.shape[0]
    vegs = list(getattr(getattr(cfg, 'surface', None), 'vegetation', []) or [])
    names = [v.name for v in vegs][:nveg] or ['type %d' % (v + 1)
                                               for v in range(nveg)]
    zr = [float(v.root_depth) for v in vegs][:nveg]
    tot = cover.sum(axis=0)
    dom = np.where(tot > 0.0, np.argmax(cover, axis=0).astype(float), np.nan)
    colors = ['#c7e9c0', '#7a0177', '#238b45', '#fec44f', '#8c510a',
              '#41b6c4'][:nveg]
    cm = ListedColormap(colors)
    nm = BoundaryNorm(np.arange(-0.5, nveg), cm.N)
    ext = uzf_land(gwf, mg)['extdp']
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    _map(axes[0], polys, dom, rings, cmap=cm, norm=nm,
         title='Vegetation -- the type covering most of each cell')
    share = [100.0 * np.mean(dom == v) for v in range(nveg)]
    axes[0].legend(handles=[Patch(facecolor=colors[v], label='%s -- root %s m '
                                  '(%.0f %% of the cells)' % (
                                      names[v], ('%g' % zr[v]) if v < len(zr)
                                      else '?', share[v]))
                            for v in range(nveg)], loc='upper right',
                   fontsize=7, title='type -- root depth (cells)')
    pc = _map(axes[1], polys, ext, rings, cmap='YlGn',
              title='UZF ET extinction depth (m), mean %.2f m'
              % np.nanmean(ext))
    fig.colorbar(pc, ax=axes[1], shrink=0.7, label='extdp (m)')
    fig.suptitle('Vegetation -> UZF ET extinction depth (each cell: its types\' '
                 'root depth below the soil column, weighted by their cover; '
                 'the uncovered fraction takes the extinction depth set for '
                 'it)', fontsize=10)
    _save(fig, out, 'IN_104_land_cover_extdp.png', written)


def fig_uzf_params(out, gwf, mg, polys, rings, written):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    u = uzf_land(gwf, mg)
    panels = [('vks', 'vks (m/d, log ramp)', 'viridis'),
              ('thts', 'thts (Sy + Sr)', 'YlGnBu'),
              ('thti', 'thti (initial)', 'YlGnBu'),
              ('thtr', 'thtr (specific retention)', 'YlOrBr'),
              ('eps', 'eps (Brooks-Corey)', 'PuBu')]
    fig, axes = plt.subplots(2, 3, figsize=(16, 11), layout='constrained')
    for ax, (k, ttl, cm) in zip(axes.ravel(), panels):
        v = u[k]
        pos = v[np.isfinite(v) & (v > 0)]
        norm = (LogNorm(vmin=pos.min(), vmax=pos.max())
                if k == 'vks' and pos.size and pos.max() > pos.min() else None)
        pc = _map(ax, polys, v, rings, cmap=cm, norm=norm,
                  title='UZF -- %s' % ttl)
        cb = fig.colorbar(pc, ax=ax, shrink=0.7)
        if norm is not None:
            _log_ticks(cb, norm.vmin, norm.vmax)
    axes.ravel()[-1].axis('off')
    fig.suptitle('UZF parameters of the land objects (as written for MF6)')
    _save(fig, out, 'IN_105_uzf_params.png', written)


def fig_forcing(out, cMF, ctx, written):
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    import pandas as pd
    P = np.asarray(getattr(ctx, 'P_veg_zoneSP'), dtype=float)
    Eo = np.asarray(getattr(ctx, 'Eo_zonesSP'), dtype=float)
    P = P.reshape(P.shape[0], -1).mean(axis=0)
    Eo = Eo.reshape(Eo.shape[0], -1).mean(axis=0)
    d = np.asarray(getattr(cMF, 'inputDate'), dtype=float)
    n = min(len(d), len(P), len(Eo), int(getattr(cMF, 'nper', len(P))))
    idx = pd.to_datetime(mdates.num2date(d[:n])).tz_localize(None)
    daily = pd.DataFrame({'P': P[:n], 'Eo': Eo[:n]}, index=idx)
    m = daily.resample('MS').sum()
    # a month the run covers only in part would read as a dry one
    full = daily['P'].resample('MS').count() == m.index.days_in_month
    m = m[full.values] if full.any() else m
    fig, ax = plt.subplots(figsize=(13, 4.5))
    # Eo in the lower half, the rainfall bars hanging in the upper half
    ax.plot(m.index, m['Eo'], color='darkorange', lw=1.5,
            label='open-water evaporation Eo (mm/month)')
    ax.set_ylabel('Eo (mm/month)', color='darkorange')
    ax.set_ylim(0, 2.0 * m['Eo'].max())
    ax2 = ax.twinx()
    ax2.bar(m.index, m['P'], width=25, color='tab:blue',
            label='rainfall (mm/month)')
    ax2.set_ylim(2.0 * m['P'].max(), 0)
    ax2.set_ylabel('rainfall (mm/month)', color='tab:blue')
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h2 + h1, l2 + l1, loc='upper center', fontsize=8, ncol=2,
              bbox_to_anchor=(0.5, -0.08))
    ax.set_title('Climatic forcing: rainfall (blue, top) and open-water '
                 'evaporation (orange, bottom), totals of the full months '
                 '[%s .. %s]'
                 % (m.index[0].strftime('%Y-%m'),
                    m.index[-1].strftime('%Y-%m')), fontsize=10)
    _save(fig, out, 'IN_106_forcing.png', written)


def _stream_lines(ds_ws):
    fn = os.path.join(ds_ws or '', 'inputSTREAM.csv')
    if not os.path.exists(fn):
        return []
    import marmites_channel as mch
    lines, _par = mch.read_stream_lines(fn)
    return [np.asarray(v, dtype=float)[:, :2] for v in lines.values()]


def _pond_ring_of(name, ponds):
    fid = str(name).lower().replace('pond', '')
    for f, r in ponds:
        if str(f) == fid:
            return f, r
    return None, None


def fig_ponds_zoom(out, mg, polys, ds_ws, lks, sfr, written):
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Patch
    from matplotlib.patches import Polygon as MplPoly
    sfr = set(sfr)
    lines = _stream_lines(ds_ws)
    ponds = _pond_rings(ds_ws)
    xc = np.asarray(mg.xcellcenters, dtype=float).ravel()
    yc = np.asarray(mg.ycellcenters, dtype=float).ravel()
    npag = int(np.ceil(len(lks) / float(PONDS_PER_ZOOM_PAGE)))
    for pg in range(npag):
        chunk = lks[pg * PONDS_PER_ZOOM_PAGE:(pg + 1) * PONDS_PER_ZOOM_PAGE]
        fig, axes = plt.subplots(2, 3, figsize=(16, 11))
        for ax in axes.ravel()[len(chunk):]:
            ax.axis('off')
        for ax, lk in zip(axes.ravel(), chunk):
            fid, ring = _pond_ring_of(lk['name'], ponds)
            cx = xc[lk['cells']]
            cy = yc[lk['cells']]
            xs = np.r_[cx, ring[:, 0]] if ring is not None else cx
            ys = np.r_[cy, ring[:, 1]] if ring is not None else cy
            x0, x1 = xs.min() - 25, xs.max() + 25
            y0, y1 = ys.min() - 25, ys.max() + 25
            win = np.flatnonzero((xc >= x0) & (xc <= x1) & (yc >= y0)
                                 & (yc <= y1))
            ax.add_collection(PatchCollection([MplPoly(polys[i]) for i in win],
                                              facecolor='none',
                                              edgecolor='0.75', lw=0.4))
            sw = [i for i in win if i in sfr]
            if sw:
                ax.add_collection(PatchCollection(
                    [MplPoly(polys[i]) for i in sw], **SFR_STYLE))
            ax.add_collection(PatchCollection(
                [MplPoly(polys[i]) for i in lk['cells']], **LAK_STYLE))
            if ring is not None:
                ax.plot(ring[:, 0], ring[:, 1], color='navy', lw=1.8)
            for ln in lines:
                ax.plot(ln[:, 0], ln[:, 1], color='tab:blue', lw=1.2)
            nsl = len(set(lk['cells']) & sfr)
            ax.set_title('pond %s: %d LAK cell(s), SFR∩LAK = %d (expect 0)'
                         % (fid if fid is not None else lk['name'],
                            len(lk['cells']), nsl), fontsize=10)
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            ax.set_aspect('equal')
            ax.tick_params(labelsize=7)
            ax.ticklabel_format(useOffset=False, style='plain')
        fig.legend(handles=[
            Patch(label='SFR reach cells', **SFR_STYLE),
            Patch(label='LAK cells', **LAK_STYLE),
            Patch(facecolor='none', edgecolor='navy', label='pond footprint')],
            loc='upper center', ncol=3, fontsize=10)
        fig.suptitle('Ponds -- LAK vs SFR cells (page %d/%d)' % (pg + 1, npag),
                     y=0.955)
        _save(fig, out, 'IN_107_ponds_zoom_p%d.png' % (pg + 1), written)



# --------------------------------------------------------------------- #
# the cross sections
# --------------------------------------------------------------------- #

def _nearest_on_lines(p, lines):
    """(point, unit direction) of the stream nearest to ``p``."""
    best = None
    for ln in lines:
        a, b = ln[:-1], ln[1:]
        ab = b - a
        L2 = np.maximum((ab ** 2).sum(axis=1), 1e-12)
        t = np.clip(((p - a) * ab).sum(axis=1) / L2, 0.0, 1.0)
        q = a + ab * t[:, None]
        d = np.hypot(*(q - p).T)
        i = int(np.argmin(d))
        if best is None or d[i] < best[0]:
            best = (d[i], q[i], ab[i] / np.sqrt(L2[i]))
    return (None, None) if best is None else (best[1], best[2])


def _cells_along(xs, ys, polys, xc, yc):
    """The cell each point lies in (-1 outside), by point-in-polygon on the
    cells near the line only."""
    from matplotlib.path import Path
    x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()
    pad = 60.0
    near = np.flatnonzero((xc >= x0 - pad) & (xc <= x1 + pad)
                          & (yc >= y0 - pad) & (yc <= y1 + pad))
    pts = np.c_[xs, ys]
    out = np.full(len(xs), -1, dtype=int)
    for i in near:
        inside = Path(polys[i]).contains_points(pts)
        out[inside & (out < 0)] = i
    return out


def _section(ax, centre, normal, half, mg, polys, land, strt, sfr, lk_cells,
             lk, title):
    xc = np.asarray(mg.xcellcenters, dtype=float).ravel()
    yc = np.asarray(mg.ycellcenters, dtype=float).ravel()
    s = np.arange(-half, half + 0.25, 0.5)
    xs, ys = centre[0] + s * normal[0], centre[1] + s * normal[1]
    cell = _cells_along(xs, ys, polys, xc, yc)
    ok = cell >= 0
    top = np.asarray(mg.top, dtype=float).ravel()
    bot = _layered(mg.botm, mg)
    idm = _layered(mg.idomain, mg) > 0
    nan = np.full(len(s), np.nan)

    def along(a):
        v = nan.copy()
        v[ok] = np.asarray(a, dtype=float)[cell[ok]]
        return v
    gl, t1 = along(land), along(top)
    b = [along(np.where(idm[L], bot[L], np.nan)) for L in range(int(mg.nlay))]
    # the initial water table: the head of each column's top active layer
    wt_cell = np.full(top.size, np.nan)
    for L in range(int(mg.nlay) - 1, -1, -1):
        wt_cell = np.where(idm[L], strt[L], wt_cell)
    wt = along(wt_cell)
    # the pond takes the place of the soil column in its LAK cell(s)
    inl = np.isin(cell, lk_cells) if lk is not None else np.zeros(len(s), bool)
    ax.fill_between(s, t1, gl, where=np.isfinite(gl) & np.isfinite(t1) & ~inl,
                    color='#d9b38c', step='mid', label='soil column (MMsoil)')
    ax.fill_between(s, b[0], t1, color='0.85', step='mid', label='MF6 layer 1')
    if len(b) > 1:
        ax.fill_between(s, b[1], b[0], color='0.65', step='mid',
                        label='MF6 layer 2')
    ax.step(s, gl, where='mid', color='saddlebrown', lw=1.4,
            label='land surface')
    ax.step(s, t1, where='mid', color='k', lw=0.9,
            label='MF6 top (soil base)')
    ax.step(s, wt, where='mid', color='royalblue', lw=1.4, ls='--',
            label='initial water table (IC)')
    lows = [np.nanmin(wt), np.nanmin(t1)]
    highs = [np.nanmax(gl)]
    # the pond, over its LAK cell(s)
    if lk is not None:
        if inl.any():
            # the pit, bed to land surface: empty above the water
            ax.fill_between(s, lk['bottom'], np.fmax(gl, lk['rim']),
                            where=inl, step='mid', facecolor='white',
                            edgecolor='none')
            ax.fill_between(s, lk['bottom'], lk['strt'], where=inl,
                            step='mid', facecolor=LAK_STYLE['facecolor'],
                            alpha=LAK_STYLE['alpha'], edgecolor='none',
                            label='pond water, initial stage')
            ax.fill_between(s, lk['bottom'], lk['rim'], where=inl,
                            step='mid', facecolor='none',
                            edgecolor=LAK_STYLE['edgecolor'],
                            lw=LAK_STYLE['lw'],
                            label='LAK cell: the pond, bed to rim')
            ax.hlines(lk['rim'], s[inl].min(), s[inl].max(), color='navy',
                      ls=':', lw=1.4, label='pond rim / outlet')
            base = float(np.mean(top[lk_cells]))
            ground = float(np.mean(np.asarray(land, dtype=float)[lk_cells]))
            ax.text(0.5, 0.97, 'spills at %.2f m, %.1f m below the land; bed '
                    '%.2f m, %.1f m below the land and %.1f m below the soil '
                    'base' % (lk['rim'], ground - lk['rim'], lk['bottom'],
                              ground - lk['bottom'], base - lk['bottom']),
                    transform=ax.transAxes, ha='center', va='top', fontsize=7,
                    color=LAK_STYLE['edgecolor'])
            lows.append(lk['bottom'])
            highs.append(lk['rim'])
    # the streambed, over the SFR cells
    rtp = nan.copy()
    rbt = nan.copy()
    for i, c in enumerate(cell):
        if c in sfr:
            rtp[i], rbt[i] = sfr[c][1], sfr[c][1] - sfr[c][2]
    if np.isfinite(rtp).any():
        ax.fill_between(s, rbt, rtp, where=np.isfinite(rtp), step='mid',
                        label='SFR cells: streambed (rtp - rbth)',
                        **SFR_STYLE)
        lows.append(np.nanmin(rbt))
    lo = np.nanmin(lows) - 3.0
    hi = np.nanmax(highs) + 1.5          # room for the pond's note
    ax.set_ylim(lo, hi)
    ax.set_xlim(-half, half)
    # the section's own cell, at 0
    ax.axvline(0.0, color='0.35', lw=0.7, ls=(0, (2, 3)))
    below = ['layer %d %.1f m' % (L + 1, np.nanmin(b[L]))
             for L in range(len(b)) if np.nanmax(b[L]) < lo]
    if below:
        ax.text(0.01, 0.02, 'bottom below the frame: %s' % ', '.join(below),
                transform=ax.transAxes, fontsize=7, color='0.3')
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('distance across the reach (m)', fontsize=8)
    ax.set_ylabel('m a.s.l.', fontsize=8)
    ax.tick_params(labelsize=7)


def fig_pond_sections(out, gwf, mg, polys, ds_ws, cMF, lks, written):
    import matplotlib.pyplot as plt
    lines = _stream_lines(ds_ws)
    if not lines:
        return
    sfr = sfr_reaches(gwf, mg)
    land = np.asarray(np.ma.filled(np.ma.asarray(cMF.elev, dtype=float),
                                   np.nan), dtype=float).ravel()
    n = int(np.asarray(mg.xcellcenters).size)
    land = land[:n] if land.size >= n else np.asarray(mg.top).ravel()
    strt = _layered(_pkg(gwf, 'ic').strt.array, mg)
    xc = np.asarray(mg.xcellcenters, dtype=float).ravel()
    yc = np.asarray(mg.ycellcenters, dtype=float).ravel()
    area = np.array([0.5 * abs(np.dot(p[:, 0], np.roll(p[:, 1], -1))
                               - np.dot(p[:, 1], np.roll(p[:, 0], -1)))
                     for p in polys])
    ponds = _pond_rings(ds_ws)
    npag = int(np.ceil(len(lks) / float(PONDS_PER_SECTION_PAGE)))
    for pg in range(npag):
        chunk = lks[pg * PONDS_PER_SECTION_PAGE:
                    (pg + 1) * PONDS_PER_SECTION_PAGE]
        fig, axes = plt.subplots(PONDS_PER_SECTION_PAGE, 2,
                                 figsize=(15, 4.6 * PONDS_PER_SECTION_PAGE),
                                 squeeze=False)
        for ax in axes.ravel()[2 * len(chunk):]:
            ax.axis('off')
        for r, lk in enumerate(chunk):
            fid, ring = _pond_ring_of(lk['name'], ponds)
            c_lak = lk['cells'][0]
            p = np.array([xc[c_lak], yc[c_lak]])
            _q, t = _nearest_on_lines(p, lines)
            nrm = np.array([-t[1], t[0]])
            half = max(60.0, 4.0 * np.sqrt(area[c_lak]))
            _section(axes[r, 0], p, nrm, half, mg, polys, land, strt, sfr,
                     lk['cells'], lk, 'pond %s -- across the reach through its '
                     'LAK cell %d' % (fid, c_lak + 1))
            # the SFR cell nearest the pond, outside its LAK cells
            cand = [c for c in sfr if c not in lk['cells']]
            if cand:
                c_s = min(cand, key=lambda c: np.hypot(xc[c] - p[0],
                                                       yc[c] - p[1]))
                ps = np.array([xc[c_s], yc[c_s]])
                _q2, t2 = _nearest_on_lines(ps, lines)
                n2 = np.array([-t2[1], t2[0]])
                _section(axes[r, 1], ps, n2, half, mg, polys, land, strt, sfr,
                         lk['cells'], lk, 'pond %s -- across the reach through '
                         'SFR cell %d (reach %d)' % (fid, c_s + 1,
                                                     sfr[c_s][0]))
        h, lab = [], []
        for ax in axes.ravel():
            for hh, ll in zip(*ax.get_legend_handles_labels()):
                if ll not in lab:
                    h.append(hh)
                    lab.append(ll)
        fig.suptitle('Ponds -- sections across the reach (initial state; the '
                     'dotted line marks the cell named; page %d/%d)'
                     % (pg + 1, npag))
        fig.tight_layout(rect=(0, 0.05, 1, 1))
        fig.legend(h, lab, loc='lower center', ncol=5, fontsize=8)
        _save(fig, out, 'IN_108_pond_sections_p%d.png' % (pg + 1), written)


# --------------------------------------------------------------------- #
# all of them
# --------------------------------------------------------------------- #

def cdl_figures(out, sim_ws, name, ds_ws, cMF=None, ctx=None, cfg=None,
                verbose=True):
    """The CdL input figures and the pond sections into ``out``; the files
    written. A figure that cannot be drawn is skipped with a note."""
    import matplotlib
    matplotlib.use('agg')
    written = []
    gwf, mg = load_model(sim_ws, name)
    if gwf is None:
        return written
    polys = _polys(mg)
    rings = _boundary(ds_ws)
    lks = lakes(gwf, mg, sim_ws)
    sfr = sfr_reaches(gwf, mg)
    jobs = [('aquifer params', lambda: fig_aquifer_params(
                out, gwf, mg, polys, rings, written)),
            ('layer surfaces', lambda: fig_layer_surfaces(
                out, mg, polys, rings, written)),
            ('ibound / UZF', lambda: fig_ibound_uzf(
                out, gwf, mg, polys, rings, ds_ws, written)),
            ('UZF parameters', lambda: fig_uzf_params(
                out, gwf, mg, polys, rings, written))]
    if ctx is not None:
        jobs.append(('land cover / extdp', lambda: fig_land_cover_extdp(
            out, gwf, mg, polys, rings, ctx, cfg, written)))
        if cMF is not None:
            jobs.append(('forcing', lambda: fig_forcing(out, cMF, ctx,
                                                        written)))
    if lks:
        jobs.append(('pond zooms', lambda: fig_ponds_zoom(
            out, mg, polys, ds_ws, lks, sfr, written)))
        if cMF is not None:
            jobs.append(('pond sections', lambda: fig_pond_sections(
                out, gwf, mg, polys, ds_ws, cMF, lks, written)))
    for what, job in jobs:
        try:
            job()
        except Exception as exc:                 # pragma: no cover
            if verbose:
                print('   %s skipped: %r' % (what, exc))
    if verbose:
        print('CdL input figures: %d file(s)' % len(written))
    return written
