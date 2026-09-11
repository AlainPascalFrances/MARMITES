# -*- coding: utf-8 -*-
"""Grid — look at the mesh a run would actually use.  WP1c.7.

The figure suite draws mesh results through a display RASTER, so a Voronoi
cell arrives there as a block of pixels. This page draws the TRUE polygons,
which is what you want when the question is "is the grid right?" rather than
"what did the model compute?".

It reads the mesh from the cache the driver writes (`<ws>/MF6_ws_<kind>/_mesh/`)
rather than building it: building needs the .ini, the time discretisation and
Triangle, which is a ~20 s round trip no page should make on every rerun.
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

import marmites_config as mcfg        # noqa: E402
import mm_paths                       # noqa: E402
from lib import loaders               # noqa: E402

st.set_page_config(page_title='Grid', page_icon='🕸️', layout='wide')
CONFIG_DIR = os.path.join(CODE, 'configs')
case = st.session_state.get('case', 'LaMata')

st.title('Grid')
st.caption('True cell polygons, straight from the mesh the driver cached. '
           'The maps in a run are drawn on a display raster instead; this is '
           'the grid itself.')

files = sorted(f for f in os.listdir(CONFIG_DIR) if f.endswith('.toml')) \
    if os.path.isdir(CONFIG_DIR) else []
if not files:
    st.error('No configuration in %s' % CONFIG_DIR)
    st.stop()

c1, c2 = st.columns([2, 3])
chosen = c1.selectbox('Configuration', files,
                      index=files.index('lamata.toml') if 'lamata.toml' in files else 0)
cfg = mcfg.load_run_config(os.path.join(CONFIG_DIR, chosen))
kind = c2.selectbox('Grid kind', list(mcfg.GRID_KINDS),
                    index=list(mcfg.GRID_KINDS).index(cfg.grid_kind))

# `paths.ws` names ONE workspace; the mesh cache lives per grid kind under
# the workspace ROOT, which is mm_paths.WS_ROOT unless paths.ws overrides it.
ws_root = (str(Path(cfg.paths.ws).parent) if cfg.paths.ws
           else str(mm_paths.WS_ROOT))
gp, sig = loaders.read_mesh(ws_root, kind)

if gp is None:
    grid_fn, _s = loaders.mesh_cache_paths(ws_root, kind)
    st.info(
        'No mesh cached for `%s` yet.\n\n'
        'It is built the first time a run uses that grid, and cached at\n\n'
        '`%s`\n\n'
        'Start one from the **Run** page with `grid.kind=%s`, or use '
        '`--set grid.kind=%s`. A structured grid has no mesh to show.'
        % (kind, grid_fn, kind, kind))
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

# ---- topology check --------------------------------------------------- #
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
                'The full mesh has %d components. That is expected here: '
                'flopy clips Voronoi cells against the domain corners and '
                'leaves a couple of non-conforming triangles. What matters '
                'is that the ACTIVE cells form one component, which the run '
                'checks — water in a second component could never reach an '
                'outlet.' % rep['components'])
        else:
            st.success('One connected component — every cell can route to an '
                       'outlet.')
    except Exception as exc:
        st.info('Topology check unavailable: %r' % exc)

# ---- the map ---------------------------------------------------------- #
colour_by = st.radio('Colour cells by', ['cell area', 'uniform'],
                     horizontal=True)
show = st.multiselect(
    'Overlay', ['stream network', 'ponds', 'observation points',
                'catchment boundary'],
    default=['stream network', 'observation points'])

DS = mm_paths.dataset_dir(case)


def _csv(fn, cols):
    """Read one of the WP1 grid-independent tables; None if absent."""
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
                   c='#17becf', edgecolor='k', linewidth=0.4, zorder=4,
                   label='ponds')
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
ax.set_xlabel('x [m], EPSG:23029')
ax.set_ylabel('y [m]')
ax.set_title('%s — %s mesh, %d cells' % (case, kind, ncpl))
st.pyplot(fig, width='content')
plt.close(fig)

st.caption('Cells are the model\'s own polygons. An observation point sits in '
           'the cell whose polygon contains it (WP1c.6) — on a coarse mesh two '
           'nearby points can share one cell, and the run says so when they do.')
