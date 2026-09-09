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


def _bounds_of(geojson_str):
    """[[south, west], [north, east]] of a GeoJSON string, for fit_bounds."""
    import json as _json

    def walk(coords, acc):
        if coords and isinstance(coords[0], (int, float)):
            x, y = coords[0], coords[1]
            acc[0] = min(acc[0], y); acc[1] = min(acc[1], x)
            acc[2] = max(acc[2], y); acc[3] = max(acc[3], x)
        else:
            for c in coords:
                walk(c, acc)
        return acc

    acc = [90.0, 180.0, -90.0, -180.0]
    for feat in _json.loads(geojson_str).get('features', []):
        geom = feat.get('geometry') or {}
        if geom.get('coordinates'):
            walk(geom['coordinates'], acc)
    if acc[0] > acc[2]:
        return None
    return [[acc[0], acc[1]], [acc[2], acc[3]]]


def _merge_bounds(a, b):
    if not a:
        return b
    if not b:
        return a
    return [[min(a[0][0], b[0][0]), min(a[0][1], b[0][1])],
            [max(a[1][0], b[1][0]), max(a[1][1], b[1][1])]]


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
                st.image(_norm(arr), clamp=True, width='content',
                         caption='%s (%d x %d)' % (rel, arr.shape[0], arr.shape[1]))
            elif p.suffix.lower() == '.csv':
                prov, header, rows = _table(str(p), loaders.digest(str(p)))
                if prov:
                    st.code('\n'.join(prov), language='text')
                st.dataframe([dict(zip(header, r)) for r in rows[:500]],
                             width='stretch')
            else:
                st.code(p.read_text(encoding='utf-8', errors='replace')[:20000])

with tab_b:
    gis = str(mm_paths.GIS)
    st.write('`%s`' % gis)
    layers = loaders.vector_layers(gis)
    chosen = [n for n, label, _c, ok in layers
              if ok and st.checkbox('%s  (`%s`)' % (label, n), value=ok,
                                    key='lay_' + n)]
    # A TOGGLE, not a button. st.button is True only on the rerun its own click
    # causes, so with a button the map vanished the moment any checkbox was
    # touched -- the same "Streamlit reruns the whole script" trap the Run page
    # is built around, walked into here in the view layer.
    if st.toggle('Show map', value=False, key='show_map'):
        try:
            import folium
            from streamlit_folium import st_folium
        except Exception as exc:
            st.info('Map needs `folium` and `streamlit-folium` (%s). The CRS '
                    'report below works without them.' % exc)
        else:
            if not chosen:
                st.info('Tick at least one layer.')
            else:
                m = folium.Map(tiles='OpenStreetMap')
                colours = {n: c for n, _l, c, _o in layers}
                bounds = None
                for name in chosen:
                    p = os.path.join(gis, name)
                    gj, _note = _vector(p, loaders.digest(p))
                    folium.GeoJson(
                        gj, name=name,
                        style_function=(lambda _f, c=colours[name]: {
                            'color': c, 'weight': 2,
                            'fillColor': c, 'fillOpacity': 0.25}),
                        marker=folium.CircleMarker(
                            radius=4, color=colours[name], fill=True),
                    ).add_to(m)
                    b = _bounds_of(gj)
                    bounds = b if bounds is None else _merge_bounds(bounds, b)
                folium.LayerControl(collapsed=False).add_to(m)
                if bounds:
                    # Without this the map opens on the whole world: folium
                    # centres on (0, 0) when given no location.
                    m.fit_bounds(bounds)
                st_folium(m, height=560, width=None, returned_objects=[])

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
