# -*- coding: utf-8 -*-
"""Panel 1 -- THE GRID.  WP1d.

The first thing a modeller settles, and the reason it is first: every other
input is WRAPPED ONTO whatever this produces. The soil zones, the vegetation
cover, the stream network and the observation points are vector layers in the
project CRS, projected onto the grid at run time -- so changing the grid does
not mean re-making any of them.

Three parts: the catchment and the grid settings, a check that the polygon is
where it should be, and the mesh a run would actually use.
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
from lib import loaders, panelui            # noqa: E402

st.set_page_config(page_title='1 Grid', page_icon='🗺️', layout='wide')
case = st.session_state.get('case', 'LaMata')
cfg, path = panelui.pick_config()
panelui.dataset_banner(cfg)
panel = panelui.header(1)

tab_domain, tab_mesh = st.tabs(['Catchment & grid', 'The mesh a run would use'])

# ===================================================================== 1a
with tab_domain:
    st.markdown('#### The catchment')
    st.caption('A PROJECTED, metric polygon in the GIS folder. It defines the '
               'active domain, the mesh boundary and the model rectangle.')

    gis = Path(mm_paths.GIS)
    shp = sorted(f.name for f in gis.glob('*.shp')) if gis.is_dir() else []
    if shp and cfg.grid.boundary not in shp:
        st.warning('`%s` is not in %s. Present: %s'
                   % (cfg.grid.boundary, gis, ', '.join(shp[:12])))

    edited = panelui.section_form(cfg, 'grid', columns=3)

    # ---- is the polygon where it should be? --------------------------
    bnd = gis / cfg.grid.boundary
    if bnd.exists():
        try:
            sys.path.insert(0, CODE)
            from marmites_vector import Layer, _signed_area
            lay = Layer(str(bnd))
            area = 0.0
            for i in range(len(lay)):
                rings = lay.rings(i)
                if rings:
                    area += abs(_signed_area(max(
                        rings, key=lambda r: abs(_signed_area(r)))))
            x0, y0, x1, y1 = lay.bbox
            c1, c2, c3 = st.columns(3)
            c1.metric('Polygon area', '%.3f km²' % (area / 1e6))
            c2.metric('Extent', '%.0f × %.0f m' % (x1 - x0, y1 - y0))
            c3.metric('Cells at %g m' % cfg.grid.cell_size,
                      '%d' % (((x1 - x0) / cfg.grid.cell_size + 1)
                              * ((y1 - y0) / cfg.grid.cell_size + 1)))
            st.caption('`%s` — %d feature(s), CRS as declared: %s'
                       % (bnd.name, len(lay),
                          (lay.crs_wkt.split('"')[1] if '"' in lay.crs_wkt
                           else 'UNDECLARED (no .prj)')))
            if x1 - x0 < 100 or y1 - y0 < 100:
                st.error('That extent is too small to be metres. The '
                         'catchment polygon must be PROJECTED, not '
                         'latitude/longitude.')
        except Exception as exc:
            st.info('Could not read the polygon: %r' % exc)
    else:
        st.error('The catchment polygon is missing: `%s`' % bnd)

    st.markdown('#### The grid built inside it')
    if cfg.grid.override.enable:
        st.warning('**Override is ON**: the grid is taken from the origin and '
                   'shape below, not derived from the polygon. That is how the '
                   'legacy 65 × 60 @ 50 m grid is reproduced when something '
                   'needs comparing against it.')
    if cfg.grid_kind in ('voronoi', 'quadtree'):
        st.info('A mesh is built the first time a run uses it and then cached. '
                'The sizes below are the side of the EQUIVALENT SQUARE — the '
                'build prints the area it actually achieved, which is what to '
                'read.')

    panelui.save_button(cfg, path, edited)

    with st.expander('What the converter would export for this grid'):
        st.markdown("""
The shapefiles stay in the GIS folder and are read **only** by the converter,
which writes grid-independent tables into the dataset. Those are what a run
reads, and what the panels below wrap onto this grid:

```bash
python code/tools/gis_to_dataset.py --case %s --config %s --dry-run
```
""" % (case, os.path.relpath(path, str(mm_paths.REPO)).replace('\\', '/')))

# ===================================================================== 1b
with tab_mesh:
    st.caption('True cell polygons, straight from the mesh the driver cached. '
               'The maps in a run are drawn on a display raster instead; this '
               'is the grid itself.')
    kind = st.selectbox('Grid kind', list(mcfg.GRID_KINDS),
                        index=list(mcfg.GRID_KINDS).index(cfg.grid_kind),
                        key='mesh_kind')
    ws_root = (str(Path(cfg.paths.ws).parent) if cfg.paths.ws
               else str(mm_paths.WS_ROOT))
    gp, sig = loaders.read_mesh(ws_root, kind)
    if gp is None:
        grid_fn, _s = loaders.mesh_cache_paths(ws_root, kind)
        st.info(
            'No mesh cached for `%s` yet.\n\n'
            'It is built the first time a run uses that grid, and cached at\n\n'
            '`%s`\n\n'
            'Start one from the **Run** panel with `grid.kind=%s`. A '
            'structured grid has no mesh to show.' % (kind, grid_fn, kind))
        st.stop()

    @st.cache_data(show_spinner='Reading the mesh...')
    def _mesh(ws_root, kind, _sig):
        g, _s = loaders.read_mesh(ws_root, kind)
        return loaders.mesh_polygons(g), int(g['ncpl'])

    (polys, areas, centres), ncpl = _mesh(ws_root, kind,
                                          (sig or {}).get('signature', ''))
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
