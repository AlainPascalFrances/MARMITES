# -*- coding: utf-8 -*-
"""Results — the figures a run wrote.  WP1b, step 1b.8.

SHOW, DO NOT REDRAW. The native suite already writes every figure into
`out_<stamp>_<tag>/_input` and `_output`; this page presents them with a run
picker and a filter. One source of truth per figure, the same rule that
produced the suite in the first place. Interactivity is added only where it
buys something, and that starts in stage 2 (WP6).

ORGANISED BY TABS (user, 2026-10-09): input maps, output maps, time series,
calibration, the water budgets and the ponds' -- lib/results.py says
which figure goes where, with its title. Two figures to a row, aligned on
their tops; each title centred BELOW its figure.

MAPS OF ONE DAY (user, 2026-10-29) are the one thing drawn here: a day and
a flux, from the run's own model outputs, while they are still in its
workspace (ppMF6/marmites_daymap) -- and the wetness W of a day against the
NDVI of a Landsat or Sentinel-2 scene of it, by linear regression
(ppMF6/marmites_wetimage; the scenes fetched by ppMF6/marmites_rs). They are
kept in the run's `_output/day_maps`.
"""

import html
import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP, os.path.join(CODE, 'ppMF6')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_daymap as DM          # noqa: E402
import marmites_rs as RS              # noqa: E402
import marmites_smrs as SM            # noqa: E402
import marmites_wetimage as WI        # noqa: E402
import mm_paths                       # noqa: E402
from lib import results as R          # noqa: E402

st.set_page_config(page_title='Results', page_icon='📊', layout='wide')
WS = str(mm_paths.WS_ROOT)

st.title('Results')
st.caption('The figures the run wrote, from `%s`.' % WS)


@st.cache_data(show_spinner=False)
def _out_dirs(ws, stamp):
    # ``stamp``, the folder's mtime, refreshes the list when a run adds its
    # folder. No leading underscore: st.cache_data leaves such arguments out
    # of its key, and the list then stayed as the app first read it (a run
    # of 2026-10-27 was missing until a restart)
    if not os.path.isdir(ws):
        return []
    return sorted((d for d in os.listdir(ws) if d.startswith('out_')),
                  reverse=True)


def _key_for(ws):
    try:
        return os.stat(ws).st_mtime_ns
    except OSError:
        return 0


dirs = _out_dirs(WS, _key_for(WS))
if not dirs:
    st.info('No `out_*` results folder in the workspace yet.')
    st.stop()

sel = st.selectbox('Run', dirs)
root = os.path.join(WS, sel)

# Provenance: WP0.5 copies the resolved configuration into every run folder.
prov = os.path.join(root, '_input', 'resolved_config.toml')
if os.path.exists(prov):
    with st.expander('Configuration that produced this run', expanded=False):
        st.code(open(prov, encoding='utf-8').read(), language='toml')
else:
    st.caption('No `resolved_config.toml` — this run predates WP0.5.')

pngs = []
for sub in R.FOLDERS:
    d = os.path.join(root, sub)
    if os.path.isdir(d):
        for f in sorted(os.listdir(d)):
            if f.lower().endswith('.png'):
                pngs.append((sub, f, os.path.join(d, f)))
summary = os.path.join(root, 'figures_nwt_comparison', '00_summary.txt')
csvs = []
for sub in ('_output',):
    d = os.path.join(root, sub)
    if os.path.isdir(d):
        csvs += [(f, os.path.join(d, f)) for f in sorted(os.listdir(d))
                 if f.lower().endswith('.csv')]

if not pngs:
    st.warning('No figures in this run folder. Re-draw them without '
               're-running the model: turn **Run MMsoil + MODFLOW 6** off '
               'with the figures on, and launch -- or\n\n`python '
               'code/tests/run_lamata_mf6.py --config code/configs/lamata.toml '
               '--set run.model=false --set run.plot=true --run-tag %s`'
               % sel.split('_', 2)[-1])
    st.stop()

needle = st.text_input('Filter by title or file name', '').strip().lower()
groups = R.arrange(pngs)
if needle:
    groups = {k: [it for it in v if needle in it[3].lower()
                  or needle in it[1].lower()] for k, v in groups.items()}
st.caption('%d of %d figure(s)' % (sum(len(v) for v in groups.values()),
                                   len(pngs)))


def _caption(sub, fname, text):
    """The title, centred under the figure, and where the file is."""
    st.markdown(
        '<div style="text-align:center;margin:-0.4rem 0 1.4rem 0">'
        '<b>%s</b><br><span style="font-size:0.75rem;color:grey">%s/%s'
        '</span></div>' % (html.escape(text), html.escape(sub),
                           html.escape(fname)),
        unsafe_allow_html=True)


def _grid(items):
    """Two figures to a row, aligned on their tops. A map of two layers
    takes the whole row, so its maps show at the size of a one-layer map's
    (user, 2026-10-09); a one-layer figure alone on its row is centred at
    the width of a column."""
    for row in R.rows(items):
        if len(row) == 2:
            cols = st.columns(2, gap='medium', vertical_alignment='top')
        elif R.panels(row[0][2]) > 1:
            cols = [st.container()]
        else:
            cols = [st.columns([1, 2, 1])[1]]
        for col, (sub, fname, path, text) in zip(cols, row):
            with col:
                st.image(path, width='stretch')
                _caption(sub, fname, text)


def _wet_ndvi(root, info):
    """The wetness W of a day against the NDVI of a satellite scene of it
    (user, 2026-10-29): a scene already fetched, or one found and fetched
    here -- Landsat 5/7/8/9 or Sentinel-2, from the Planetary Computer, into
    the RS folder. The regression is per model cell, per vegetation class;
    the scene's true colour is only the visual support."""
    st.markdown('**The wetness W against the NDVI of a satellite scene**')
    dates = info['dates']
    days = set(dates)
    folders = {s: RS.folder(mm_paths.RS, s) for s in RS.SENSORS}
    prods = [p for p in RS.products(*folders.values())
             if RS.as_date(p['date']) in days]
    if prods:
        c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment='bottom')
        lab = ['%s  %s  (%g m, %.0f %% of the window clear)'
               % (p['date'], p['platform'], p['res'],
                  100 * p.get('clear_share_window', 0)) for p in prods]
        i = c1.selectbox('Scene fetched (its NDVI against W of its day)',
                         range(len(prods)), format_func=lab.__getitem__,
                         key='nd_prod_%s' % sel)
        # a cell's class: its main vegetation type when it covers this much
        # of the cell, else bare soil (user, 2026-10-29)
        mc = c2.number_input('Main vegetation, at least % of the cell', 0.0,
                             100.0, WI.MIN_COVER, step=5.0, key='nd_mincover')
        if c3.button('Compare', key='nd_go', type='primary'):
            with st.spinner('NDVI of %s against W...' % lab[i].split('  (')[0]):
                try:
                    WI.compare(root, prods[i], info=info, gis=mm_paths.GIS,
                               min_cover=mc)
                except Exception as exc:              # noqa: BLE001
                    st.error('Not compared: %s' % exc)
    else:
        st.caption('No scene fetched for the days of this run yet: find one '
                   'below.')
    with st.expander('Find and fetch a scene (Planetary Computer, no account)'):
        c1, c2, c3, c4 = st.columns([2, 2, 1, 1], vertical_alignment='bottom')
        sensor = c1.selectbox('Product', list(RS.SENSORS),
                              format_func=lambda s: RS.SENSORS[s]['label'],
                              key='nd_sensor')
        span = c2.date_input('Between', value=(dates[0], dates[-1]),
                             min_value=dates[0], max_value=dates[-1],
                             key='nd_span_%s' % sel)
        cloud = c3.number_input('Scene cloud, at most %', 0, 100, 20,
                                key='nd_cloud')
        if c4.button('Search', key='nd_search') and len(span) == 2:
            with st.spinner('Searching the catalogue...'):
                try:
                    g = DM._grid(info['sim_ws'], info['name'])
                    b = RS.catchment_bounds(g)
                    found = RS.search(sensor, RS.lonlat_bbox(b, info['crs']),
                                      span[0], span[1], cloud)
                    st.session_state['nd_found'] = (sel, sensor, found, b)
                except Exception as exc:              # noqa: BLE001
                    st.error('The catalogue did not answer: %s' % exc)
        hit = st.session_state.get('nd_found')
        if hit and hit[0] == sel and hit[1] == sensor:
            scenes = [s_ for s_ in hit[2] if RS.as_date(s_['date']) in days]
            if not scenes:
                st.info('No %s scene over the catchment in that window.'
                        % RS.SENSORS[sensor]['label'])
            else:
                have = {(p['platform'], p['date']) for p in prods}
                lab = ['%s  %s  scene cloud %s %%%s'
                       % (s_['date'], s_['platform'],
                          '?' if s_['cloud'] is None else '%.0f' % s_['cloud'],
                          '  (fetched)' if (s_['platform'], s_['date']) in have
                          else '') for s_ in scenes]
                c1, c2 = st.columns([4, 1], vertical_alignment='bottom')
                j = c1.selectbox('%d scene(s)' % len(scenes), range(len(scenes)),
                                 format_func=lab.__getitem__, key='nd_scene')
                if c2.button('Fetch', key='nd_fetch'):
                    with st.spinner('Reading %s over the catchment...'
                                    % scenes[j]['id']):
                        try:
                            RS.fetch(scenes[j], info['crs'], hit[3],
                                     folders[sensor])
                            st.rerun()
                        except Exception as exc:      # noqa: BLE001
                            st.error('Not fetched: %s' % exc)
        st.caption('Into `%s`: the NDVI and the true colour on the model\'s '
                   'grid; clouds, their shadows, snow and missing data '
                   'masked.' % folders[sensor])


def _sm_products(root):
    """Satellite soil moisture against the model's saturation Se, catchment
    mean (ppMF6/marmites_smrs; user, 2026-10-29): the products downloaded
    into the SOIL MOISTURE folder of panel 0 (mm_paths.SOIL_MOISTURE, by
    default RS/soil_moisture), compared on request with the run's own
    outputs while they are still in its workspace. The figures join this
    tab (the run's _output/sm_products). A product's folder may hold the
    Copernicus Browser's Statistical Info exports instead (2026-10-30): the
    catchment to draw there is offered as GeoJSON."""
    have = SM.available(mm_paths.SOIL_MOISTURE)
    with st.expander("Satellite soil moisture against the model's saturation "
                     "Se", expanded=False):
        try:
            info = DM.outputs_of(root)
        except Exception as exc:                      # noqa: BLE001
            info = {'why': 'its model outputs could not be read (%s)' % exc}
        if not info.get('why') and info.get('crs'):
            try:
                gj = _catchment_geojson(
                    info['sim_ws'], info['name'], info['crs'],
                    max(_key_for(os.path.join(info['sim_ws'], info['name'] + x))
                        for x in ('.dis', '.disv')))
                st.download_button(
                    'The catchment as GeoJSON (WGS84), to upload as the area '
                    'in the Copernicus Browser', gj, key='sm_aoi',
                    file_name='%s_catchment_wgs84.geojson' % info['name'],
                    mime='application/geo+json')
            except Exception:                         # noqa: BLE001
                pass
        if not have:
            st.caption('No product in `%s` yet (Copernicus SWI 0.1 deg or '
                       '1 km, Surface Soil Moisture 1 km): the SOIL MOISTURE '
                       'folder of the Home panel, which says what each is.'
                       % mm_paths.SOIL_MOISTURE)
            return
        if info.get('why'):
            st.info('No comparison for this run: %s.' % info['why'])
            return
        c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment='bottom')
        p = c1.selectbox('Soil-moisture product', list(have), key='sm_prod',
                         format_func=lambda k: '%s (%d file(s))'
                         % (SM.PRODUCTS[k]['label'], have[k]))
        T = c2.number_input('T, days (SWI T-value / filter)', 1, 100,
                            SM.T_DEFAULT, key='sm_T')
        if c3.button('Compare', key='sm_go', type='primary'):
            with st.spinner('Se against %s...' % SM.PRODUCTS[p]['short']):
                try:
                    SM.compare(root, p, SM.folder(mm_paths.SOIL_MOISTURE, p),
                               info=info, T=int(T))
                    st.rerun()
                except Exception as exc:              # noqa: BLE001
                    st.error('Not compared: %s' % exc)
        st.caption('An SWI is compared as it is (both % of saturation); a '
                   'surface product is filtered to T days, then both are '
                   'scaled to their own 2nd-98th percentiles. Scores over '
                   'the whole hydrological years. A folder holds either the '
                   "product's files or the Browser's Statistical Info "
                   'exports, one CSV per layer (SWI001..SWI100).')


@st.cache_data(show_spinner=False)
def _catchment_geojson(sim_ws, name, crs, stamp):
    # ``stamp``: the grid file's mtime, a new grid there redraws it
    return SM.catchment_geojson({'sim_ws': sim_ws, 'name': name, 'crs': crs})


def _day_maps(root):
    """A day and a flux, drawn from the run's own model outputs while they
    are still in its workspace -- the next run there overwrites them --
    then every map of one day drawn for this run, the latest first."""
    try:
        info = DM.outputs_of(root)
    except Exception as exc:                          # noqa: BLE001
        info = {'why': 'its model outputs could not be read (%s)' % exc}
    if info.get('why'):
        st.info('No new map of one day for this run: %s.' % info['why'])
    else:
        dates = info['dates']
        c1, c2, c3 = st.columns([1, 3, 1], vertical_alignment='bottom')
        # one key per run: another run has other days
        day = c1.date_input('Day', value=dates[-1], min_value=dates[0],
                            max_value=dates[-1], key='dm_day_%s' % sel)
        key = c2.selectbox('Flux', [f[0] for f in DM.FLUXES],
                           format_func=DM.label, key='dm_flux')
        if c3.button('Draw', key='dm_draw', type='primary'):
            with st.spinner('Reading %s for %s (the cell budget takes about '
                            '20 s the first time)...' % (key, day)):
                try:
                    DM.draw(root, day, key, info=info)
                except Exception as exc:              # noqa: BLE001
                    st.error('Not drawn: %s' % exc)
        st.caption('From the model outputs in `%s`, on the linear and the '
                   'percentile colour ramp; kept in this run\'s '
                   '`_output/%s`.' % (info['sim_ws'], DM.DAY_DIR))
        _wet_ndvi(root, info)
    items = DM.drawn(root)
    if items:
        _grid([(sub, f, p, R.title(sub, f)) for sub, f, p in items])
    else:
        st.caption('No map of one day drawn for this run yet.')


# the day maps' tab always, and without a count: a label that changes when a
# map is drawn would send the page back to the first tab
keys = [k for k, _l in R.TABS if groups.get(k) or k == 'daymap'
        or (k == 'budget' and os.path.exists(summary) and not needle)]
tabs = st.tabs([dict(R.TABS)[k] if k == 'daymap' else
                '%s (%d)' % (dict(R.TABS)[k], len(groups.get(k, [])))
                for k in keys])
for key, tab in zip(keys, tabs):
    with tab:
        if key == 'daymap':
            _day_maps(root)
            continue
        if key == 'calib':
            _sm_products(root)
        if key == 'budget' and os.path.exists(summary) and not needle:
            st.code(open(summary, encoding='utf-8', errors='replace').read(),
                    language='text')
            _caption('figures_nwt_comparison', '00_summary.txt',
                     R.title('figures_nwt_comparison', '00_summary.txt'))
        _grid(groups.get(key, []))

if csvs:
    with st.expander('CSV output (%d)' % len(csvs)):
        for fname, fpath in csvs:
            st.write('`%s`' % fname)
            st.download_button('download %s' % fname,
                               data=open(fpath, 'rb').read(),
                               file_name=fname, key='dl_' + fname)
