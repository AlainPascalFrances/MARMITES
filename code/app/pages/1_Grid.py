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
from lib import editor, loaders, panelui    # noqa: E402

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
    if cfg.grid_kind in ('structured', 'dis'):
        info = {'kind': 'structured', 'ncpl': int(stub.nrow * stub.ncol),
                'area_mean': snap * snap, 'signature': 'rectangle',
                'nrow': int(stub.nrow), 'ncol': int(stub.ncol)}
        return True, ['Structured grid: %d x %d = %d cells of %g m — nothing '
                      'to build, it IS the rectangle.'
                      % (stub.nrow, stub.ncol, info['ncpl'], snap)] + lines, info
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

tab_domain, tab_mesh = st.tabs(['Catchment & grid', 'The selected mesh'])

# ===================================================================== 1a
with tab_domain:
    st.markdown('#### The catchment and the grid built inside it')
    st.caption('A PROJECTED, metric polygon. It defines the active domain, '
               'the mesh boundary and the model rectangle — so it is checked '
               'here, before anything is built on it.')

    edited, chosen = panelui.grid_permanent_form(cfg, columns=3)

    # ---- is the file a usable catchment polygon? ---------------------
    from marmites_vector import check_polygon_layer          # noqa: E402
    picked = edited.get('grid.boundary', cfg.grid.boundary)
    bnd = (picked if os.path.isabs(picked)
           else os.path.join(str(mm_paths.GIS), picked))
    rep = check_polygon_layer(bnd, expect_epsg=edited.get('grid.crs_epsg',
                                                          cfg.grid.crs_epsg))
    for msg in rep['errors']:
        st.error(msg)
    for msg in rep['warnings']:
        st.warning(msg)

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
    if rep['ok']:
        st.success('Valid catchment polygon.')

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
                   'achieved — read that, not this. The transition bands are '
                   'derived from the corridor and the maximum size ratio.')
    elif chosen == 'quadtree':
        st.caption('GRIDGEN halves a cell per refinement level, so level 2 on '
                   'a %g m background gives %g m along the streams.'
                   % (cfg.grid.cell_size,
                      cfg.grid.cell_size / (2 ** cfg.grid.quadtree.refine_level)))
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
               'until you select one.')

    cbuild, cclear, cmsg = st.columns([1, 1, 3])
    if cbuild.button('Create grid', type='primary', key='mkgrid'):
        try:
            trial, _todo = editor.apply_changes(cfg, edited)
        except (editor.EditError, mcfg.ConfigError) as exc:
            cmsg.error('These settings are not valid, so nothing was built:'
                       '\n\n%s' % exc)
        else:
            tag = '%s_%d' % (trial.grid_kind, len(ATTEMPTS) + 1)
            with st.spinner('Building the %s grid…' % trial.grid_kind):
                ok, lines, info = _build_grid(trial, _attempt_dir(cfg, tag))
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

        names = [a['tag'] for a in ATTEMPTS]
        pick = st.selectbox('Attempt to inspect or select', names,
                            index=len(names) - 1, key='pick_attempt')
        att = ATTEMPTS[names.index(pick)]
        (st.success if att['ok'] else st.error)(att['lines'][0])
        with st.expander('Build log', expanded=not att['ok']):
            st.code('\n'.join(att['lines']), language='text')

        # ---- commit --------------------------------------------------
        st.markdown('#### Select this grid for the model')
        st.caption('Writes the settings that produced **%s** into `%s`. From '
                   'then on, that is the grid a run builds and every other '
                   'panel wraps its layers onto.' % (pick, os.path.basename(path)))
        csel, cnote = st.columns([1, 3])
        if csel.button('Select this grid for the model', type='primary',
                       key='selgrid', disabled=not att['ok']):
            try:
                applied, digest = editor.save(cfg, path, edited)
            except (editor.EditError, mcfg.ConfigError) as exc:
                cnote.error('NOT selected — the configuration would be '
                            'invalid:\n\n%s' % exc)
            else:
                st.session_state['grid_selected'] = pick
                cnote.success('Selected **%s** — %d change(s), hash %s'
                              % (pick, len(applied), digest))
                # The saved spin-up state belongs to the grid it was produced
                # on. Saying so here is the difference between a clear message
                # now and a CONFIG ERROR at the start of the next run.
                stale = [n for n in ('strt_heads', 'steady_means')
                         if getattr(cfg.spinup, n)]
                if stale and att['cfg'].grid_kind not in ('structured', 'disv'):
                    cnote.warning(
                        'This grid is a mesh, and `spinup.%s` still names '
                        'state produced on the structured grid. A run will '
                        'stop rather than feed MODFLOW an array of the wrong '
                        'length: clear those keys, or re-run the spin-up on '
                        'this mesh.' % '` / `spinup.'.join(stale))
                st.rerun()
        if att['ok'] and not att['cfg'].grid_kind == cfg.grid_kind:
            cnote.info('The configuration currently says **%s**; this attempt '
                       'is **%s**.' % (cfg.grid_kind, att['cfg'].grid_kind))
    else:
        st.info('No attempt yet. Press **Create grid** — it builds from the '
                'settings above without saving anything.')

    with st.expander('Update the dataset from the cartography'):
        st.caption('The shapefiles stay in the GIS folder and are read ONLY by '
                   'the converter, which writes grid-independent tables into '
                   'the dataset. Those are what a run reads, and what the '
                   'other panels wrap onto this grid.')
        c1, c2 = st.columns(2)
        if c1.button('Preview (dry run)'):
            st.session_state['conv'] = _run_converter(case, path, dry=True)
        if c2.button('Update dataset', type='primary'):
            st.session_state['conv'] = _run_converter(case, path, dry=False)
        if st.session_state.get('conv'):
            st.code(st.session_state['conv'], language='text')

# ===================================================================== 1b
with tab_mesh:
    st.caption('True cell polygons, straight from the mesh a run would use. '
               'The maps in a run are drawn on a display raster instead; this '
               'is the grid itself.')
    kind = cfg.grid_kind
    st.markdown('**Selected: `%s`** — %s' % (kind, _describe(cfg)))

    ws_root = _ws_root(cfg)
    gp, sig = loaders.read_mesh(ws_root, kind)
    if gp is None:
        # Fall back to the attempt built on this page, so a mesh can be looked
        # at before a run has ever been launched on it.
        for a in reversed(ATTEMPTS):
            if a['ok'] and a['cfg'].grid_kind == kind:
                gp, sig = loaders.read_mesh(a['cache'], kind)
                if gp is None:
                    gp, sig = loaders.read_mesh(os.path.dirname(a['cache']),
                                                kind)
                break
    if gp is None:
        grid_fn, _s = loaders.mesh_cache_paths(ws_root, kind)
        st.info(
            'No mesh cached for `%s` yet.\n\n'
            'Build one with **Create grid** on the first tab, or start a run '
            'from the **Run** panel; a run caches it at\n\n`%s`\n\n'
            'A structured grid has no mesh to show — it is the rectangle.'
            % (kind, grid_fn))
        st.stop()

    @st.cache_data(show_spinner='Reading the mesh...')
    def _mesh(ws_root, kind, _sig):
        g, _s = loaders.read_mesh(ws_root, kind)
        return loaders.mesh_polygons(g), int(g['ncpl'])

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
        'Overlay', ['stream network', 'ponds', 'observation points',
                    'catchment boundary'],
        default=['stream network', 'observation points'])

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

    import matplotlib                                     # noqa: E402
    matplotlib.use('agg')
    import matplotlib.pyplot as plt                       # noqa: E402
    from matplotlib.collections import PolyCollection     # noqa: E402

    fig, ax = plt.subplots(figsize=(9, 9))
    pc = PolyCollection(polys, edgecolors='#44444455', linewidths=0.3)
    if colour_by == 'cell area':
        pc.set_array(areas)
        pc.set_cmap('viridis')
        fig.colorbar(pc, ax=ax, shrink=0.7, label='cell area [m²]')
    else:
        pc.set_facecolor('#dfe7f5')
    ax.add_collection(pc)

    if 'stream network' in show:
        segs = {}
        for sid, seq, x, y in (_csv('inputSTREAM.csv',
                                    ['seg_id', 'seq', 'x', 'y']) or []):
            segs.setdefault(sid, []).append((seq, x, y))
        for sid in sorted(segs):
            pts = [(x, y) for _s, x, y in sorted(segs[sid])]
            ax.plot([p[0] for p in pts], [p[1] for p in pts], '-',
                    color='#1f77b4', lw=1.2, zorder=3)
    if 'ponds' in show:
        pd_ = _csv('inputPONDS.csv', ['fid', 'x', 'y']) or []
        if pd_:
            ax.scatter([p[1] for p in pd_], [p[2] for p in pd_], s=28,
                       c='#17becf', edgecolor='k', linewidth=0.4, zorder=4)
    if 'catchment boundary' in show:
        ring = _csv('inputWATERSHED.csv', ['ring_id', 'seq', 'x', 'y']) or []
        if ring:
            pts = [(x, y) for _r, _s, x, y in sorted(ring, key=lambda t: t[1])]
            ax.plot([p[0] for p in pts], [p[1] for p in pts], '--',
                    color='#d62728', lw=1.2, zorder=3)
    if 'observation points' in show:
        try:
            from marmites_postprocess import obs_points
            pts = obs_points(str(DS))
            ax.scatter([p['x'] for p in pts], [p['y'] for p in pts], s=45,
                       marker='^', c='#ffcc00', edgecolor='k', linewidth=0.5,
                       zorder=5)
            for p in pts:
                ax.annotate(p['name'], (p['x'], p['y']), fontsize=7,
                            xytext=(3, 3), textcoords='offset points', zorder=6)
        except Exception as exc:
            st.caption('observation points not drawn: %r' % exc)

    ax.autoscale_view()
    ax.set_aspect('equal')
    ax.set_xlabel('x [m], EPSG:%d' % (cfg.grid.crs_epsg or 23029))
    ax.set_ylabel('y [m]')
    ax.set_title('%s — %s mesh, %d cells' % (case, kind, ncpl))
    st.pyplot(fig, width='content')
    plt.close(fig)

    st.caption('Cells are the model\'s own polygons. An observation point sits '
               'in the cell whose polygon contains it — on a coarse mesh two '
               'nearby points can share one cell, and the run says so.')
