# -*- coding: utf-8 -*-
"""Panel 1 -- THE GRID.  WP1d.

The first thing a modeller settles, and the reason it is first: every other
input is WRAPPED ONTO whatever this produces. The soil zones, the vegetation
cover, the stream network and the observation points are vector layers in the
project CRS, projected onto the grid at run time -- so changing the grid does
not mean re-making any of them.

Two buttons, and the difference between them is the whole design of the page:

  Create grid       an EXPERIMENT. Builds from what is on screen, keeps the
                    attempt, writes nothing the model reads. Press it as
                    often as you like, on as many settings as you like.
  Select this grid  the COMMITMENT. Writes the settings that produced the
                    chosen attempt into the configuration, and that is the
                    grid a run will use.

There is deliberately no separate "Validate & save" here: on this panel a
save IS the selection, and two buttons that both write would only differ in
whether the mesh had been looked at first.
"""

import os
import sys
from pathlib import Path

import numpy as np
import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP, os.path.join(CODE, 'ppMF6')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg              # noqa: E402
import mm_paths                             # noqa: E402
from lib import editor, loaders, panelui, schema   # noqa: E402

st.set_page_config(page_title='1 Grid', page_icon='🗺️', layout='wide')
case = st.session_state.get('case', 'LaMata')
cfg, path = panelui.pick_config()
panelui.dataset_banner(cfg)
panel = panelui.header(1)


def _boundary_bbox(cfg):
    """The catchment polygon's bounding box, or None with the reason."""
    from marmites_vector import Layer, VectorError
    p = os.path.join(str(mm_paths.GIS), cfg.grid.boundary)
    try:
        return Layer(p).bbox, None
    except (VectorError, Exception) as exc:          # noqa: B014
        return None, '%r' % exc


def _ws_root(cfg):
    return (str(Path(cfg.paths.ws).parent) if cfg.paths.ws
            else str(mm_paths.WS_ROOT))


def _attempt_dir(cfg, tag):
    """Where one EXPERIMENT's mesh is cached.

    Keyed on the attempt rather than on the kind alone: two voronoi meshes at
    different cell sizes have to coexist, or they cannot be compared -- which
    is what made the old single Create grid button feel like it worked once.
    """
    return os.path.join(_ws_root(cfg), '_grid_attempts', tag)


def _build_grid(cfg, cache_dir, force=True):
    """Build the mesh for THIS configuration object. (ok, [lines], info)."""
    import marmites_meshes as mm

    lines = []

    def warn(msg):
        lines.append('WARNING: %s' % msg)

    bbox, why = _boundary_bbox(cfg)
    if bbox is None and not cfg.grid.override.enable:
        return False, ['The catchment polygon could not be read, so the grid '
                       'rectangle cannot be derived: %s' % why], None
    try:
        stub = mm.grid_stub(cfg, bbox, nlay=cfg.layers.nlay)
    except Exception as exc:
        return False, ['%r' % exc], None
    snap = mm.rectangle_cell_size(cfg)
    lines.append('rectangle: %d rows x %d cols of %g m, origin %.1f, %.1f'
                 % (stub.nrow, stub.ncol, snap, stub.xllcorner, stub.yllcorner))
    os.makedirs(cache_dir, exist_ok=True)
    try:
        _gp, info = mm.build_mesh(cfg, stub, cache_dir=cache_dir,
                                  dataset_dir=str(mm_paths.dataset_dir(
                                      cfg.paths.case)),
                                  model_ws=cache_dir, warn=warn, force=force)
    except Exception as exc:
        return False, ['The %s producer failed: %r'
                       % (cfg.grid_kind, exc)] + lines, None
    head = ('%s mesh: %d cells, mean %.0f m² (%.1f m equivalent side)'
            % (info['kind'], info['ncpl'], info.get('area_mean', 0.0),
               (info.get('area_mean', 0.0) ** 0.5)))
    lines.append('signature %s' % info['signature'])
    lines.append('cached at %s' % cache_dir)
    return True, [head] + lines, info


def _signature(cfg):
    """The mesh signature THIS configuration asks for, or None.

    The cache is keyed on it, so it is also the only honest way for the
    viewer to answer "is the mesh on screen the one these settings describe?"
    -- comparing cell counts would pass a mesh built with the same count and
    a different corridor.
    """
    import marmites_meshes as mm
    bbox, _why = _boundary_bbox(cfg)
    try:
        return mm.mesh_signature(cfg, mm.grid_stub(cfg, bbox,
                                                   nlay=cfg.layers.nlay))
    except Exception:                                    # noqa: BLE001
        return None


def _find_mesh(cfg, attempts, want):
    """The cached mesh for ``want``, wherever it is. (gp, sig, where, fresh).

    Looks through this session's attempts first and the run cache second,
    and reports whether what it found actually MATCHES the signature -- a
    stale mesh drawn without comment is how panel 1 came to show a grid
    nobody had asked for.
    """
    kind = cfg.grid_kind
    seen = []
    for a in reversed(attempts):
        if a['ok'] and a['cfg'].grid_kind == kind:
            gp, sig = loaders.read_mesh_at(a['cache'], kind)
            if gp is not None:
                seen.append((gp, sig, 'attempt %s' % a['tag']))
    gp, sig = loaders.read_mesh(_ws_root(cfg), kind)
    if gp is not None:
        seen.append((gp, sig, 'the run cache'))
    for gp, sig, where in seen:
        if want and (sig or {}).get('signature') == want:
            return gp, sig, where, True
    return (seen[0] + (False,)) if seen else (None, None, '', False)


def _settings_of(trial):
    """The whole ``[grid]`` block of an attempt, as edits to save.

    Taken from the attempt's own configuration rather than from the widgets:
    selecting is the moment the grid on screen becomes the model's, and by
    then the boxes may well have moved on to the next experiment.
    """
    return {k: v for k, v in schema.fields_of(trial, 'grid')
            if k not in schema.GRID_DERIVED and not schema.is_source(v)}


def _have_plotly():
    try:
        import plotly.graph_objects            # noqa: F401
        return True
    except ImportError:
        return False


def _plotly_mesh(polys, areas, kind, ncpl, colour_by, overlays, epsg, title):
    """The mesh as an INTERACTIVE figure: plotly's own zoom, pan and reset.

    One trace per area class rather than one per cell -- 4000 traces would
    make the page unusable -- with the polygons separated by ``None`` inside
    each trace, which is how plotly draws many outlines in one go.
    """
    import numpy as np
    import plotly.graph_objects as go
    from matplotlib import colormaps

    fig = go.Figure()
    if colour_by == 'cell area':
        classes = loaders.area_classes(areas)
        cmap = colormaps['viridis']
        for k, (lo, hi, sel) in enumerate(classes):
            xs, ys = [], []
            for i in np.nonzero(sel)[0]:
                p = polys[i]
                xs.extend([q[0] for q in p] + [p[0][0], None])
                ys.extend([q[1] for q in p] + [p[0][1], None])
            r, g, b, _a = cmap(k / max(len(classes) - 1, 1))
            fig.add_trace(go.Scatter(
                x=xs, y=ys, fill='toself', mode='lines',
                fillcolor='rgb(%d,%d,%d)' % (r * 255, g * 255, b * 255),
                line=dict(color='rgba(68,68,68,0.35)', width=0.4),
                name=('%.0f m² (every cell)' % lo if hi <= lo
                      else '%.0f–%.0f m²' % (lo, hi)),
                hoverinfo='name'))
    else:
        xs, ys = [], []
        for p in polys:
            xs.extend([q[0] for q in p] + [p[0][0], None])
            ys.extend([q[1] for q in p] + [p[0][1], None])
        fig.add_trace(go.Scatter(
            x=xs, y=ys, fill='toself', mode='lines',
            fillcolor='#dfe7f5', line=dict(color='rgba(68,68,68,0.35)',
                                           width=0.4),
            name='cells', hoverinfo='skip'))

    for name, (xs, ys, how) in overlays.items():
        if how in ('lak', 'sfr'):
            # FILLED CELLS, not an outline: the point of drawing them is to
            # see which cells the package will own, so the cell is the thing
            # that has to be visible.
            lak = how == 'lak'
            fig.add_trace(go.Scatter(
                x=xs, y=ys, mode='lines', name=name, fill='toself',
                fillcolor='rgba(255,140,0,0.75)' if lak
                else 'rgba(31,119,180,0.45)',
                line=dict(color='#8a4b00' if lak else '#12537e', width=0.6),
                hoverinfo='name'))
        elif how == 'poly':
            fig.add_trace(go.Scatter(
                x=xs, y=ys, mode='lines', name=name, fill='toself',
                fillcolor='rgba(23,190,207,0.45)',
                line=dict(color='#0b6d78', width=1.4), hoverinfo='name'))
        elif how == 'line':
            fig.add_trace(go.Scatter(x=xs, y=ys, mode='lines', name=name,
                                     line=dict(color='#1f77b4', width=1.6)))
        elif how == 'dash':
            fig.add_trace(go.Scatter(x=xs, y=ys, mode='lines', name=name,
                                     line=dict(color='#d62728', width=1.6,
                                               dash='dash')))
        else:
            fig.add_trace(go.Scatter(
                x=xs, y=ys, mode='markers', name=name,
                marker=dict(size=9, color='#ffcc00', symbol='triangle-up',
                            line=dict(color='black', width=0.6)),
                text=how if isinstance(how, list) else None,
                hovertemplate='%{text}<br>%{x:.0f}, %{y:.0f}<extra></extra>'
                if isinstance(how, list) else None))

    fig.update_layout(
        title=title, height=760, margin=dict(l=10, r=10, t=50, b=10),
        xaxis_title='x [m], EPSG:%d' % epsg, yaxis_title='y [m]',
        showlegend=True, dragmode='pan',
        legend=dict(orientation='h', y=-0.08))
    # Equal aspect, or a mesh looks stretched and cells look like rectangles.
    fig.update_yaxes(scaleanchor='x', scaleratio=1)
    return fig


def _static_mesh(polys, areas, colour_by, overlays, epsg, title, view):
    """The same map as a PICTURE, for when plotly is not installed."""
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    fig, ax = plt.subplots(figsize=(9, 9))
    pc = PolyCollection(polys, edgecolors='#44444455', linewidths=0.3)
    if colour_by == 'cell area':
        pc.set_array(areas)
        pc.set_cmap('viridis')
        fig.colorbar(pc, ax=ax, shrink=0.7, label='cell area [m²]')
    else:
        pc.set_facecolor('#dfe7f5')
    ax.add_collection(pc)

    for name, (xs, ys, how) in overlays.items():
        if how in ('lak', 'sfr'):
            ax.fill(xs, ys,
                    facecolor='#ff8c00bf' if how == 'lak' else '#1f77b473',
                    edgecolor='#8a4b00' if how == 'lak' else '#12537e',
                    lw=0.6, zorder=2, label=name)
        elif how == 'poly':
            ax.fill(xs, ys, facecolor='#17becf70', edgecolor='#0b6d78',
                    lw=1.2, zorder=3, label=name)
        elif how in ('line', 'dash'):
            ax.plot(xs, ys, '--' if how == 'dash' else '-',
                    color='#d62728' if how == 'dash' else '#1f77b4',
                    lw=1.2, zorder=3, label=name)
        elif how == 'points':
            ax.scatter(xs, ys, s=28, c='#17becf', edgecolor='k',
                       linewidth=0.4, zorder=4, label=name)
        else:
            ax.scatter(xs, ys, s=45, marker='^', c='#ffcc00', edgecolor='k',
                       linewidth=0.5, zorder=5, label=name)
            for x, y, nm in zip(xs, ys, how):
                ax.annotate(nm, (x, y), fontsize=7, xytext=(3, 3),
                            textcoords='offset points', zorder=6)

    ax.set_xlim(view[0], view[2])
    ax.set_ylim(view[1], view[3])
    ax.set_aspect('equal')
    ax.set_xlabel('x [m], EPSG:%d' % epsg)
    ax.set_ylabel('y [m]')
    ax.set_title(title)
    return fig


def _provenance(path):
    """The ``# source / # size`` header the converter writes, as a dict."""
    out = {}
    try:
        with open(path, encoding='utf-8') as fh:
            for line in fh:
                if not line.startswith('#'):
                    break
                if ':' in line:
                    k, v = line[1:].split(':', 1)
                    out[k.strip()] = v.strip()
    except OSError:
        return {}
    return out


# What the MESH PRODUCERS read out of the dataset, and the shapefile each is
# derived from. Voronoi triangulates inputWATERSHED.csv and refines on
# inputSTREAM.csv; the quadtree refines on the same streams. So these two are
# the dataset tables a grid actually depends on -- and the reason the
# converter has to run BEFORE a build rather than after a selection.
MESH_INPUTS = (('inputWATERSHED.csv', None),          # None = cfg.grid.boundary
               ('inputSTREAM.csv', 'hydrography.shp'))


def _dataset_stale(cfg):
    """Is what the producers read still what the cartography says? (bool, why).

    Compares the SOURCE and its size recorded in each table's header against
    the shapefile on disk. A missing table, a table made from a different
    file, or a file that has been re-exported since all mean the same thing:
    a mesh built now would be built on the previous cartography.
    """
    ds = str(mm_paths.dataset_dir(cfg.paths.case))
    why = []
    for table, source in MESH_INPUTS:
        out = os.path.join(ds, table)
        name = source or cfg.grid.boundary
        src = (name if os.path.isabs(name)
               else os.path.join(str(mm_paths.GIS), name))
        if not os.path.exists(out):
            why.append('%s has never been written' % table)
            continue
        prov = _provenance(out)
        was = prov.get('source', '')
        if was and os.path.normcase(was) != os.path.normcase(src):
            why.append('%s was made from %s, not %s'
                       % (table, os.path.basename(was), os.path.basename(src)))
            continue
        if not os.path.exists(src):
            continue                      # the check below would be noise
        size = prov.get('size', '')
        try:
            on_disk = os.path.getsize(src)
        except OSError:
            continue
        if size and ('%d bytes' % on_disk) not in size:
            why.append('%s has changed since %s was written'
                       % (os.path.basename(src), table))
    return bool(why), '; '.join(why)


def _trial_of(cfg, edited):
    """The configuration AS EDITED ON SCREEN, or None if it is not valid yet.

    The checks below have to describe the grid the button would build, not
    the one the file still holds -- otherwise changing the kind leaves the
    warnings talking about the previous one.
    """
    try:
        trial, _todo = editor.apply_changes(cfg, edited)
    except (editor.EditError, mcfg.ConfigError):
        return None
    return trial


def _rectangle_check(cfg, bbox):
    """Whether the derived rectangle stands on the dataset's own rasters."""
    import marmites_meshes as mm

    if cfg is None or not bbox:
        return None
    try:
        return mm.rectangle_check(cfg, bbox,
                                  str(mm_paths.dataset_dir(cfg.paths.case)))
    except Exception as exc:                       # never break the page
        return {'status': 'error', 'detail': repr(exc), 'derived': None,
                'source': None, 'overhang': {}, 'names': [], 'others': []}


def _show_rectangle_check(rc, kind):
    """Draw it. Returns True when the legacy rectangle can be adopted.

    STOPGAP (WP1d). The model is still assembled from rasters frozen on one
    lattice, and `project_model` resamples that assembly onto the grid: a
    grid that reaches past them fails inside the projection, cells first,
    with a message about tops and bottoms rather than about rectangles. So
    the disagreement is stated HERE, where the rectangle is chosen. It
    disappears for good when the MODEL panel re-derives those rasters onto
    whatever rectangle this panel produces.
    """
    if rc is None or rc['status'] == 'error':
        return False
    if rc['status'] == 'none':
        st.caption('No raster in the dataset yet, so there is nothing for '
                   'this rectangle to disagree with.')
        return False

    d, s = rc['derived'], rc['source']
    line = ('grid `%.10g, %.10g` %d × %d @ %g m — rasters `%.10g, %.10g` %d × %d @ %g m'
            % (d[0], d[1], d[2], d[3], d[4], s[0], s[1], s[2], s[3], s[4]))
    if rc['status'] == 'ok':
        st.caption('✅ This rectangle stands on the %d raster(s) the model '
                   'reads — %s.' % (len(rc['names']), line))
    elif rc['status'] == 'shifted':
        st.warning('**The grid and the model rasters are not on the same '
                   'lattice.** %s — %s. Every cell is then resampled from '
                   'fractions of four, which blurs the legacy model instead '
                   'of reproducing it.' % (line, rc['detail']))
    else:
        st.warning(
            '**This grid does not stand on the rasters the model reads.** '
            '%s: %s.\n\nThe grid itself will build. What fails is the MODEL '
            'build, later: it assembles the model from those rasters and '
            'resamples it onto these cells, and the cells with nothing '
            'underneath come out with their top at or below their bottom. '
            'The real cure is for the model panel to re-derive its own '
            'rasters onto this rectangle; that is not built yet.'
            % (line, rc['detail']))
    for rect, names in rc['others']:
        st.caption('⚠️ %d raster(s) sit on a THIRD rectangle — `%.10g, %.10g` %d × '
                   '%d @ %g m: %s. Stale exports, most likely, but the model '
                   'would read them as it reads the rest.'
                   % (len(names), rect[0], rect[1], rect[2], rect[3], rect[4],
                      ', '.join(names)))
    return rc['status'] in ('overhang', 'shifted') \
        and kind in ('structured', 'dis', 'disv')


def _run_converter(case, cfg_path, dry):
    """Run the WP1 converter as a subprocess and return its output.

    A subprocess, not an import: it keeps geopandas out of this process's
    import graph, and it is exactly the command a user would type.
    """
    import subprocess
    cmd = [sys.executable, os.path.join(CODE, 'tools', 'gis_to_dataset.py'),
           '--case', case, '--config', cfg_path]
    if dry:
        cmd.append('--dry-run')
    r = subprocess.run(cmd, capture_output=True, text=True,
                       cwd=str(mm_paths.REPO), timeout=900)
    return (r.stdout or '') + (('\n' + r.stderr) if r.stderr else '')


def _describe(cfg):
    """One line naming the grid a configuration asks for."""
    k = cfg.grid_kind
    if k == 'voronoi':
        v = cfg.grid.voronoi
        return ('voronoi, %g m background%s'
                % (v.cell_far,
                   (', %g m corridor to %g m' % (v.cell_near_stream,
                                                 v.stream_buffer))
                   if v.stream_refine else ', unrefined'))
    if k == 'quadtree':
        q = cfg.grid.quadtree
        return ('quadtree, %g m background, %d level(s)%s'
                % (cfg.grid.cell_size, q.refine_level,
                   '' if q.refine_streams else ', unrefined'))
    return '%s, %g m cells' % (k, cfg.grid.cell_size)


ATTEMPTS = st.session_state.setdefault('grid_attempts', [])

tab_domain, tab_mesh = st.tabs(['Catchment & grid', 'Visualize grid'])

# ===================================================================== 1a
with tab_domain:
    st.markdown('#### Catchment and grid definition')
    st.caption('Catchment boundary: polygon shape file with projected '
               'coordinates (metric). It defines the active domain, the mesh '
               'boundary and the model rectangle.')

    # The form CHECKS every layer as it draws it -- the CRS it shows comes
    # from the catchment's check -- and hands the reports back for the
    # read-outs below.
    edited, chosen, reports = panelui.grid_permanent_form(cfg, columns=3)
    rep = reports['grid.boundary']
    bnd = rep['path']
    for dotted, r in reports.items():
        label = schema.describe(dotted)[0]
        for msg in r['errors']:
            st.error('**%s** — %s' % (label, msg))
        for msg in r['warnings']:
            st.warning('**%s** — %s' % (label, msg))
    ok = [schema.describe(d)[0] for d, r in reports.items() if r['ok']]
    if len(ok) == len(reports) and reports:
        st.success('Valid: %s.' % ', '.join(ok).lower())

    if rep['bbox']:
        x0, y0, x1, y1 = rep['bbox']
        # The cell count follows the KIND chosen just above, and the size that
        # kind actually uses -- read live, because the box it comes from is
        # drawn further down the page.
        snap = float(st.session_state.get(
            'grid.voronoi.cell_far' if chosen == 'voronoi' else 'grid.cell_size',
            (cfg.grid.voronoi.cell_far if chosen == 'voronoi'
             else cfg.grid.cell_size)) or 0.0)
        c1, c2, c3 = st.columns(3)
        c1.metric('Polygon area', '%.3f km²' % (rep['area_m2'] / 1e6))
        c2.metric('Extent', '%.0f × %.0f m' % (x1 - x0, y1 - y0))
        c3.metric('Cells at %g m' % snap,
                  '%d' % (((x1 - x0) / snap + 1) * ((y1 - y0) / snap + 1))
                  if snap > 0 else '—')
        st.caption('`%s` — %d %s feature(s), CRS as declared: **%s**%s'
                   % (os.path.basename(bnd), rep['features'],
                      rep['kind'] or '?',
                      rep['crs_name'] or 'UNDECLARED (no .prj)',
                      ' (EPSG:%d)' % rep['epsg'] if rep['epsg'] else ''))
    edited.update(panelui.grid_kind_form(cfg, chosen, edited, columns=3))

    if chosen in ('structured', 'dis'):
        st.caption('A structured grid needs nothing beyond the cell size: it '
                   'is the rectangle above, divided. It is the regression '
                   'anchor, not the default.')
    elif chosen == 'disv':
        st.caption('The structured grid re-expressed as polygons. Nothing '
                   'about the geometry changes, which is what makes it the '
                   'control for the mesh path.')
    elif chosen == 'voronoi':
        st.caption('Sizes are the side of the EQUIVALENT SQUARE, so 100 aims '
                   'at 10 000 m². The build prints the area it actually '
                   'achieved — read that, not this. You give the size AT the '
                   'stream and how fast it may grow; the corridor follows, '
                   'because each band is as wide as its own cells.')
        v = panelui.live_voronoi(cfg, edited)
        bands = v.graded_bands()
        if bands:
            # Every column a STRING: a mixed int/str column cannot be
            # converted to Arrow, and streamlit's repair pass logs a
            # traceback on every render.
            lo, rows = 0.0, []
            for i, (d, s) in enumerate(bands, start=1):
                rows.append({'band': '%d' % i, 'from [m]': '%g' % lo,
                             'to [m]': '%g' % d, 'cell size [m]': '%g' % s})
                lo = d
            rows.append({'band': 'background', 'from [m]': '%g' % lo,
                         'to [m]': '∞', 'cell size [m]': '%g' % v.cell_far})
            st.dataframe(rows, width='stretch', hide_index=True)
        elif v.stream_refine:
            st.info('No bands: the size at the stream must be smaller than '
                    'the background for a corridor to mean anything.')
    if cfg.grid.override.enable:
        st.warning('**Override is ON**: the grid is taken from the origin and '
                   'shape above, not derived from the polygon. That is how the '
                   'legacy 65 × 60 @ 50 m grid is reproduced when something '
                   'needs comparing against it.')

    # ---- experiment --------------------------------------------------
    st.markdown('#### Create the grid')
    st.caption('Builds the mesh from the settings ABOVE AS THEY STAND — no '
               'save needed — and keeps it as an attempt, so two kinds or two '
               'cell sizes can be compared. Nothing a run reads is written '
               'until you select one (panel **Visualize grid**).')

    # Said BEFORE the button, because it is about the rectangle the button
    # would use and the answer changes nothing about how to build it -- only
    # about whether the model can later be put on it.
    _trial = _trial_of(cfg, edited)
    _rc = _rectangle_check(_trial, rep['bbox'])
    if _show_rectangle_check(_rc, _trial.grid_kind if _trial else chosen):
        _s = _rc['source']
        if st.button('Pin the grid to the rasters’ rectangle',
                     key='adoptrect',
                     help='Turns the override ON and fills it with the '
                          'rasters\' own origin, shape and cell size. Legal '
                          'only on a structured or disv grid — on a mesh the '
                          'override has no meaning, so there the answer is '
                          'the model panel.'):
            panelui.park('grid.override.enable', True)
            panelui.park('grid.override.xllcorner', float(_s[0]))
            panelui.park('grid.override.yllcorner', float(_s[1]))
            panelui.park('grid.override.nrow', int(_s[2]))
            panelui.park('grid.override.ncol', int(_s[3]))
            panelui.park('grid.cell_size', float(_s[4]))
            st.rerun()
        st.caption('It pins the grid to a rectangle that was chosen years '
                   'ago, cutting %s of the catchment. Use it to reproduce '
                   'the legacy model, not to build a new one.'
                   % ('the western edge' if 'west' in (_rc['overhang'] or {})
                      else 'an edge'))

    cbuild, cclear, cmsg = st.columns([1, 1, 3])
    if cbuild.button('Create grid', type='primary', key='mkgrid'):
        try:
            trial, _todo = editor.apply_changes(cfg, edited)
        except (editor.EditError, mcfg.ConfigError) as exc:
            cmsg.error('These settings are not valid, so nothing was built:'
                       '\n\n%s' % exc)
        else:
            tag = '%s_%d' % (trial.grid_kind, len(ATTEMPTS) + 1)
            pre = []
            # The producers read the DATASET, not the shapefiles: voronoi
            # triangulates inputWATERSHED.csv and refines on
            # inputSTREAM.csv. So a boundary changed on this page means
            # nothing to the mesh until the converter has run -- which is why
            # it runs HERE, silently, and not on the selection.
            stale, why = _dataset_stale(trial)
            if stale:
                with st.spinner('Re-reading the cartography…'):
                    pre = ['the dataset was out of date (%s), so the '
                           'converter ran first' % why,
                           _run_converter(case, path, dry=False), '']
            with st.spinner('Building the %s grid…' % trial.grid_kind):
                ok, lines, info = _build_grid(trial, _attempt_dir(cfg, tag))
            lines = lines[:1] + pre + lines[1:] if lines else pre
            ATTEMPTS.append({'tag': tag, 'ok': ok, 'lines': lines,
                             'info': info, 'cfg': trial,
                             'label': _describe(trial),
                             'cache': _attempt_dir(cfg, tag)})
    if cclear.button('Clear attempts', key='clrgrid') and ATTEMPTS:
        ATTEMPTS.clear()
        st.rerun()

    if ATTEMPTS:
        st.markdown('##### Attempts this session')
        rows = []
        for a in ATTEMPTS:
            i = a['info'] or {}
            rows.append({
                'attempt': a['tag'],
                'settings': a['label'],
                'cells': i.get('ncpl', '—'),
                'mean cell [m²]': ('%.0f' % i['area_mean']
                                   if i.get('area_mean') else '—'),
                'equivalent side [m]': ('%.1f' % (i['area_mean'] ** 0.5)
                                        if i.get('area_mean') else '—'),
                'built': '✅' if a['ok'] else '❌',
            })
        st.dataframe(rows, width='stretch', hide_index=True)
        st.success('After the grid is created, go to panel **Visualize grid**.')
    else:
        st.info('No attempt yet. Press **Create grid** — it builds from the '
                'settings above without saving anything.')


# ===================================================================== 1b
with tab_mesh:
    ws_root = _ws_root(cfg)
    want = _signature(cfg)

    # What can be shown: this session's attempts, plus whatever is cached for
    # the configured kind. Choosing here rather than on the first tab is the
    # point -- an attempt is chosen by LOOKING at it.
    choices = [('%s — %s' % (a['tag'], a['label']), a) for a in ATTEMPTS
               if a['ok']]
    saved_gp, saved_sig, saved_where, fresh = _find_mesh(cfg, [], want)
    if saved_gp is not None:
        choices.append(('the %s cached for `%s` (%s)'
                        % ('mesh' if cfg.grid_kind not in ('structured', 'dis')
                           else 'grid', cfg.grid_kind, saved_where), None))
    if not choices:
        grid_fn, _s = loaders.mesh_cache_paths(ws_root, cfg.grid_kind)
        st.info(
            'Nothing to show for `%s` yet.\n\n'
            'Press **Create grid** on the first tab — it builds from the '
            'settings on screen without saving anything. A run caches its own '
            'at\n\n`%s`' % (cfg.grid_kind, grid_fn))
        st.stop()

    labels = [c[0] for c in choices]
    pick = st.selectbox('Attempt to inspect or select', labels,
                        index=len(labels) - 1, key='pick_attempt')
    att = dict(choices)[pick]

    if att is not None:
        gp, sig = loaders.read_mesh_at(att['cache'], att['cfg'].grid_kind)
        kind = att['cfg'].grid_kind
        shown_cfg = att['cfg']
        (st.success if att['ok'] else st.error)(att['lines'][0])
        with st.expander('Build log', expanded=not att['ok']):
            st.code('\n'.join(att['lines']), language='text')
        if gp is None:
            st.error('The mesh for %s is no longer in %s.'
                     % (att['tag'], att['cache']))
            st.stop()
    else:
        gp, sig, kind, shown_cfg = saved_gp, saved_sig, cfg.grid_kind, cfg
        if not fresh:
            st.warning(
                'This is **not** the grid the settings describe. It comes '
                'from %s, built under a different `[grid]` block — its '
                'signature is `%s` and the settings ask for `%s`. Press '
                '**Create grid** on the first tab.'
                % (saved_where, (sig or {}).get('signature', 'unknown'),
                   want or '?'))
        else:
            st.caption('From %s — signature `%s`, which is what the current '
                       'settings ask for.' % (saved_where, want))

    st.markdown('**`%s`** — %s' % (kind, _describe(shown_cfg)))

    (polys, areas, centres) = loaders.mesh_polygons(gp)
    ncpl = int(gp['ncpl'])
    m1, m2, m3, m4 = st.columns(4)
    m1.metric('cells (ncpl)', '%d' % ncpl)
    m2.metric('mean cell', '%.0f m²' % areas.mean())
    m3.metric('equivalent side', '%.1f m' % np.sqrt(areas.mean()))
    m4.metric('smallest / largest', '%.0f / %.0f m²' % (areas.min(), areas.max()))
    st.caption('mesh signature `%s` — the cache is keyed on it, so a mesh built '
               'under a different `[grid]` block is rebuilt rather than reused.'
               % (sig or {}).get('signature', '?'))

    with st.expander('Topology check (WP1c.5)', expanded=False):
        try:
            from marmites_topology import MeshTopology
            topo = MeshTopology(gp)
            rep = topo.check(raise_on_error=False)
            a, b, c = st.columns(3)
            a.metric('faces', '%d' % rep['faces'])
            b.metric('neighbours / cell', '%.2f' % rep['mean_neighbours'])
            c.metric('components', '%d' % rep['components'])
            st.write('interior faces **%d**, boundary cells **%d**, '
                     'duplicate vertices merged **%d**, isolated cells **%d**'
                     % (rep['interior_faces'], rep['boundary_cells'],
                        rep['merged_vertices'], rep['isolated_cells']))
            if rep['components'] > 1:
                st.warning(
                    'The full mesh has %d components. That is expected: flopy '
                    'clips Voronoi cells against the domain corners and leaves '
                    'a couple of non-conforming triangles. What matters is '
                    'that the ACTIVE cells form one component, which the run '
                    'checks -- water in a second component could never reach '
                    'an outlet.' % rep['components'])
            else:
                st.success('One connected component — every cell can route to '
                           'an outlet.')
        except Exception as exc:
            st.info('Topology check unavailable: %r' % exc)

    colour_by = st.radio('Colour cells by', ['cell area', 'uniform'],
                         horizontal=True)
    show = st.multiselect(
        'Overlay', ['stream network', 'ponds', 'pond cells (LAK)',
                    'stream cells (SFR)', 'observation points',
                    'catchment boundary'],
        default=['stream network', 'observation points'])

    # ---- how the map is drawn -----------------------------------------
    # Plotly carries its own zoom, pan, box-zoom and reset, which is what a
    # mesh wants. Without it the figure reaches the browser as a picture, so
    # the view becomes a STATE that buttons move and the axes are set from.
    interactive = _have_plotly() and st.toggle(
        'Interactive map — drag to pan, wheel or box to zoom, and the toolbar '
        'has "Reset axes"', value=True, key='interactive_map')

    _vx = [float(v[1]) for v in gp['vertices']]
    _vy = [float(v[2]) for v in gp['vertices']]
    full = (min(_vx), min(_vy), max(_vx), max(_vy))
    view = full
    if not interactive:
        # The full extent is what "Original extent" goes back to -- and what
        # a newly selected mesh resets to, since a window from the previous
        # grid would be meaningless on this one.
        if st.session_state.get('view_of') != (pick, ncpl):
            st.session_state['view_of'] = (pick, ncpl)
            st.session_state['view'] = None
        view = st.session_state.get('view') or full

        def _set(v):
            st.session_state['view'] = v
            st.rerun()

        def _zoom(f):
            x0, y0, x1, y1 = view
            cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
            w, h = (x1 - x0) / 2, (y1 - y0) / 2
            _set((cx - w * f, cy - h * f, cx + w * f, cy + h * f))

        def _pan(dx, dy):
            x0, y0, x1, y1 = view
            sx, sy = (x1 - x0) * dx * 0.25, (y1 - y0) * dy * 0.25
            _set((x0 + sx, y0 + sy, x1 + sx, y1 + sy))

        z1, z2, z3, z4, z5, z6, z7 = st.columns(7)
        if z1.button('🔍 +', key='zin', help='Zoom in'):
            _zoom(1 / 1.6)
        if z2.button('🔍 −', key='zout', help='Zoom out'):
            _zoom(1.6)
        if z3.button('⟲ Original extent', key='zreset'):
            _set(None)
        if z4.button('←', key='pleft'):
            _pan(-1, 0)
        if z5.button('→', key='pright'):
            _pan(1, 0)
        if z6.button('↑', key='pup'):
            _pan(0, 1)
        if z7.button('↓', key='pdown'):
            _pan(0, -1)
        if st.session_state.get('view'):
            x0, y0, x1, y1 = view
            st.caption('view %.0f–%.0f × %.0f–%.0f m (%.0f × %.0f m)'
                       % (x0, x1, y0, y1, x1 - x0, y1 - y0))

    DS = mm_paths.dataset_dir(case)

    def _csv(fn, cols):
        p = os.path.join(str(DS), fn)
        if not os.path.exists(p):
            return None
        rows, header = [], None
        with open(p, encoding='utf-8') as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith('#'):
                    continue
                if header is None:
                    header = [x.strip() for x in s.split(',')]
                    continue
                rec = dict(zip(header, s.split(',')))
                try:
                    rows.append([float(rec[c]) if c not in ('seg_id', 'seq',
                                                            'ring_id', 'fid')
                                 else int(float(rec[c])) for c in cols])
                except (KeyError, ValueError):
                    continue
        return rows or None

    # ---- the overlays, gathered once ----------------------------------
    # Both renderers draw the same things, so they are collected here rather
    # than twice: (xs, ys, how), with None breaks between parts.
    overlays = {}
    if 'stream network' in show:
        segs = {}
        for sid, seq, x, y in (_csv('inputSTREAM.csv',
                                    ['seg_id', 'seq', 'x', 'y']) or []):
            segs.setdefault(sid, []).append((seq, x, y))
        xs, ys = [], []
        for sid in sorted(segs):
            pts = [(x, y) for _s, x, y in sorted(segs[sid])]
            xs.extend([p[0] for p in pts] + [None])
            ys.extend([p[1] for p in pts] + [None])
        if xs:
            overlays['stream network'] = (xs, ys, 'line')
    if 'catchment boundary' in show:
        ring = _csv('inputWATERSHED.csv', ['ring_id', 'seq', 'x', 'y']) or []
        if ring:
            pts = [(x, y) for _r, _s, x, y in sorted(ring, key=lambda t: t[1])]
            overlays['catchment boundary'] = ([p[0] for p in pts],
                                              [p[1] for p in pts], 'dash')
    if 'ponds' in show:
        # The FOOTPRINT, not the centroid: a pond is a polygon the mesh has
        # to cover, and a dot says nothing about whether it does.
        import marmites_meshes as _mmesh
        rings = _mmesh.pond_rings(str(DS))
        if rings:
            xs, ys = [], []
            for ring in rings:
                xs.extend([p[0] for p in ring] + [ring[0][0], None])
                ys.extend([p[1] for p in ring] + [ring[0][1], None])
            overlays['ponds'] = (xs, ys, 'poly')
        else:
            pd_ = _csv('inputPONDS.csv', ['fid', 'x', 'y']) or []
            if pd_:
                overlays['ponds'] = ([p[1] for p in pd_], [p[2] for p in pd_],
                                     'points')
    # ---- the cells the packages will own (WP4 read ahead of time) -----
    # Drawn so that a refinement can be SEEN to have worked. A pond that owns
    # no cell of its own, or a lake that does not touch the network it is
    # supposed to drain into, is a thing to find out HERE and not in the LAK
    # file.
    want_lak = 'pond cells (LAK)' in show
    want_sfr = 'stream cells (SFR)' in show
    if want_lak or want_sfr:
        import marmites_meshes as _mm

        def _cells(idx):
            xs, ys = [], []
            for i in idx:
                p = polys[i]
                xs.extend([q[0] for q in p] + [p[0][0], None])
                ys.extend([q[1] for q in p] + [p[0][1], None])
            return xs, ys

        rings = _mm.pond_rings(str(DS))
        segs = []
        try:
            segs = _mm.stream_lines(os.path.join(str(DS), 'inputSTREAM.csv'))
        except Exception:
            segs = []
        per_pond = _mm.pond_cells(gp, rings) if rings else []
        lake = sorted({i for cells in per_pond for i in cells})
        line = _mm.stream_cells(gp, segs) if segs else []
        # A cell cannot be both a lake and a reach. THE LAKE TAKES IT: the
        # stream is cut where it enters the pond and rejoined below it, and
        # what carries the water across is a mover, not a reach inside the
        # water. So the cells drawn as SFR are the line's cells MINUS the
        # lake's -- the picture says what the packages will do.
        sfr = [i for i in line if i not in set(lake)]
        if want_lak and lake:
            xs, ys = _cells(lake)
            overlays['pond cells (LAK)'] = (xs, ys, 'lak')
        if want_sfr and sfr:
            xs, ys = _cells(sfr)
            overlays['stream cells (SFR)'] = (xs, ys, 'sfr')
        if per_pond:
            loose = [k + 1 for k, cells in enumerate(per_pond)
                     if not _mm.cells_touch(gp, cells, sfr)]
            sizes = [len(c) for c in per_pond]
            taken = len(line) - len(sfr)
            st.caption(
                '%d pond(s) own %d cell(s) between them (%d to %d each), and '
                '%s. The mapped stream runs through %d cell(s); **%d of them '
                'are inside a pond and belong to the LAKE**, so the reaches '
                'stop at the rim and a mover carries the flow across — that '
                'is panel 3, not this one. %d stay as reaches.'
                % (len(per_pond), sum(sizes), min(sizes), max(sizes),
                   'every one touches a reach' if not loose
                   else 'pond(s) %s touch NO reach'
                   % ', '.join(str(k) for k in loose),
                   len(line), taken, len(sfr)))

    obs = []
    if 'observation points' in show:
        try:
            from marmites_postprocess import obs_points
            obs = obs_points(str(DS))
            if obs:
                overlays['observation points'] = (
                    [p['x'] for p in obs], [p['y'] for p in obs],
                    [p['name'] for p in obs])
        except Exception as exc:
            st.caption('observation points not drawn: %r' % exc)

    title = '%s — %s mesh, %d cells' % (case, kind, ncpl)
    if interactive:
        st.plotly_chart(
            _plotly_mesh(polys, areas, kind, ncpl, colour_by, overlays,
                         cfg.grid.crs_epsg or 23029, title),
            width='stretch',
            config={'scrollZoom': True, 'displaylogo': False})
        st.caption('Cells are the model\'s own polygons. Drag to pan, scroll '
                   'or drag a box to zoom, and double-click or **Reset axes** '
                   'in the toolbar to come back to the full extent. An '
                   'observation point sits in the cell whose polygon contains '
                   'it -- on a coarse mesh two nearby points share one, and '
                   'the run says so.')
    else:
        import matplotlib.pyplot as plt
        fig = _static_mesh(polys, areas, colour_by, overlays,
                           cfg.grid.crs_epsg or 23029,
                           title + ('' if st.session_state.get('view') is None
                                    else ' (zoomed)'), view)
        st.pyplot(fig, width='content')
        plt.close(fig)
        st.caption('Cells are the model\'s own polygons. An observation point '
                   'sits in the cell whose polygon contains it -- on a coarse '
                   'mesh two nearby points share one, and the run says so.')

    # ---- commit, AFTER the map ---------------------------------------
    # Below the figure because that is the order of the decision: look at the
    # mesh, then commit to it.
    if att is not None:
        st.markdown('---')
        if st.button('Select this grid for the model', type='primary',
                     key='selgrid', width='stretch'):
            try:
                applied, digest = editor.save(cfg, path,
                                              _settings_of(att['cfg']))
            except (editor.EditError, mcfg.ConfigError) as exc:
                st.error('NOT selected — the configuration would be '
                         'invalid:\n\n%s' % exc)
            else:
                st.session_state['grid_selected'] = att['tag']
                st.success('Selected **%s** — %d change(s), hash %s'
                           % (att['tag'], len(applied), digest))
                # The settings are written; put the mesh they produced where
                # the driver looks, so the run reuses it instead of spending
                # the build again. The signature goes with it, so a later
                # change to [grid] still invalidates it.
                dst, moved = loaders.promote_mesh(att['cache'], ws_root, kind)
                if moved:
                    st.caption('Mesh promoted to `%s` — a run will reuse it '
                               'rather than rebuild.' % dst)
                # Saved spin-up state belongs to the grid it was produced on.
                # Saying so HERE is the difference between a clear message now
                # and a CONFIG ERROR at the start of the next run.
                stale = [n for n in ('strt_heads', 'steady_means')
                         if getattr(cfg.spinup, n)]
                if stale and kind not in ('structured', 'dis', 'disv'):
                    st.warning(
                        'This grid is a mesh, and `spinup.%s` still names '
                        'state produced on the structured grid. A run will '
                        'stop rather than feed MODFLOW an array of the wrong '
                        'length: clear those keys, or re-run the spin-up on '
                        'this mesh.' % '` / `spinup.'.join(stale))
        if kind != cfg.grid_kind:
            st.info('The configuration currently says **%s**; this attempt is '
                    '**%s**.' % (cfg.grid_kind, kind))
