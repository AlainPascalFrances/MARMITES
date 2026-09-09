# -*- coding: utf-8 -*-
"""Inputs — the Tier-A files and the Tier-B cartography.  WP1b, 1b.2/1b.3/1b.9."""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

import mm_paths                       # noqa: E402
from lib import loaders               # noqa: E402

st.set_page_config(page_title='Inputs', page_icon='📂', layout='wide')
case = st.session_state.get('case', 'LaMata')
DS = mm_paths.dataset_dir(case)

st.title('Inputs — %s' % case)
st.caption('Tier A is what MM and MF read directly and lives in the repository. '
           'Tier B is the source cartography and lives outside it, in '
           '`$MM_DATA_ROOT/GIS`. This page is one of only two places that '
           'imports geopandas; the model path has none.')


# ---- caching: key on (path, size, mtime), never on the path alone -----------
@st.cache_data(show_spinner=False)
def _asc(path, _key):
    return loaders.read_asc(path)


@st.cache_data(show_spinner=False)
def _table(path, _key):
    return loaders.read_table(path)


@st.cache_data(show_spinner=False)
def _vector(path, _key):
    gdf, note = loaders.read_vector(path)
    return gdf.to_json(), note


def _norm(arr):
    """Scale a raster to 0-1 for display, ignoring nodata."""
    import numpy as np
    a = np.asarray(arr, dtype=float)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return np.zeros_like(a)
    lo, hi = float(finite.min()), float(finite.max())
    out = (a - lo) / (hi - lo) if hi > lo else np.zeros_like(a)
    return np.nan_to_num(out, nan=0.0)


def _run_converter(case, dry):
    """Run the WP1 converter as a subprocess and return its output.

    A subprocess, not an import: it keeps geopandas out of this process's
    import graph even here, and it is exactly the command a user would type.
    """
    import subprocess
    cmd = [sys.executable, os.path.join(CODE, 'tools', 'gis_to_dataset.py'),
           '--case', case]
    if dry:
        cmd.append('--dry-run')
    r = subprocess.run(cmd, capture_output=True, text=True,
                       cwd=str(mm_paths.REPO), timeout=900)
    return (r.stdout or '') + (('\n' + r.stderr) if r.stderr else '')


tab_a, tab_b = st.tabs(['Tier A — the model inputs', 'Tier B — cartography'])

with tab_a:
    if not DS.is_dir():
        st.error('Case folder not found: %s' % DS)
    else:
        inv = loaders.inventory(DS)
        missing = [f for _g, items in inv for f, _s, ok in items if not ok]
        if missing:
            st.warning('%d listed input(s) missing: %s'
                       % (len(missing), ', '.join(missing[:8])))
        for group, items in inv:
            with st.expander('%s  (%d)' % (group, len(items)),
                             expanded=group.startswith('Generated')):
                for rel, size, ok in items:
                    c1, c2, c3 = st.columns([6, 2, 2])
                    c1.write(('✅ ' if ok else '❌ ') + '`%s`' % rel)
                    c2.write('%.1f kB' % (size / 1024.0) if ok else '—')
                    if ok and c3.button('view', key='v_' + rel):
                        st.session_state['view'] = rel

        rel = st.session_state.get('view')
        if rel:
            p = DS / rel
            st.markdown('---')
            st.subheader(rel)
            if p.suffix.lower() == '.asc' or p.suffix.upper() == '.ASC':
                arr, hdr = _asc(str(p), loaders.digest(str(p)))
                st.write(hdr)
                st.image(_norm(arr), clamp=True, use_container_width=False,
                         caption='%s (%d x %d)' % (rel, arr.shape[0], arr.shape[1]))
            elif p.suffix.lower() == '.csv':
                prov, header, rows = _table(str(p), loaders.digest(str(p)))
                if prov:
                    st.code('\n'.join(prov), language='text')
                st.dataframe([dict(zip(header, r)) for r in rows[:500]],
                             use_container_width=True)
            else:
                st.code(p.read_text(encoding='utf-8', errors='replace')[:20000])

with tab_b:
    gis = str(mm_paths.GIS)
    st.write('`%s`' % gis)
    layers = loaders.vector_layers(gis)
    chosen = [n for n, label, _c, ok in layers
              if ok and st.checkbox('%s  (`%s`)' % (label, n), value=ok,
                                    key='lay_' + n)]
    if st.button('Show map'):
        try:
            import folium
            from streamlit_folium import st_folium
        except Exception:
            st.info('Install `folium` and `streamlit-folium` (see '
                    '`environment-ui.yml`) to draw the map. The CRS report '
                    'below works without them.')
            chosen = chosen
        else:
            m = folium.Map(tiles='OpenStreetMap')
            bounds = None
            for name in chosen:
                gj, note = _vector(os.path.join(gis, name),
                                   loaders.digest(os.path.join(gis, name)))
                colour = dict((n, c) for n, _l, c, _o in layers)[name]
                folium.GeoJson(gj, name=name,
                               style_function=lambda _f, c=colour: {
                                   'color': c, 'weight': 2, 'fillOpacity': 0.2},
                               tooltip=folium.GeoJsonTooltip(
                                   fields=[], aliases=[])).add_to(m)
            folium.LayerControl().add_to(m)
            st_folium(m, height=560, use_container_width=True)

    st.markdown('#### CRS')
    st.caption('Reported, never guessed: the La Mata DEM rasters carry an '
               'equivalent but *unnamed* projection with no EPSG code, and a '
               'layer with no `.prj` is flagged rather than silently misplaced.')
    for name, label, _c, ok in layers:
        if not ok:
            continue
        try:
            _gj, note = _vector(os.path.join(gis, name),
                                loaders.digest(os.path.join(gis, name)))
            st.write('`%-24s` %s' % (name, note))
        except Exception as exc:
            st.write('`%-24s` unreadable: %s' % (name, exc))

    st.markdown('---')
    st.markdown('#### Update the dataset from the cartography  *(WP1.9)*')
    st.caption('Runs `code/tools/gis_to_dataset.py`. Preview first: it '
               'overwrites the generated Tier-A tables.')
    c1, c2 = st.columns(2)
    if c1.button('Preview (dry run)'):
        st.session_state['conv'] = _run_converter(case, dry=True)
    if c2.button('Update dataset', type='primary'):
        st.session_state['conv'] = _run_converter(case, dry=False)
    if st.session_state.get('conv'):
        st.code(st.session_state['conv'], language='text')
