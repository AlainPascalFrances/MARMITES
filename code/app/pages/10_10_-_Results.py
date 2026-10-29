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
workspace (ppMF6/marmites_daymap) -- and the wetness W of a day against a
georeferenced image of it (ppMF6/marmites_wetimage). They are kept in the
run's `_output/day_maps`.
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


def _wet_image(root, info):
    """The wetness W of a day against a georeferenced RGB image of it: the
    GeoTIFFs in the workspace and the remote-sensing folder, the day the
    file name carries by default."""
    st.markdown('**The wetness W against an image of the same day**')
    found = WI.images(WS, mm_paths.RS)
    if not found:
        st.caption('No GeoTIFF in `%s` or `%s`: put a georeferenced RGB image '
                   'of a day of the run there.' % (WS, mm_paths.RS))
        return
    dates = info['dates']
    c1, c2, c3 = st.columns([3, 1, 1], vertical_alignment='bottom')
    img = c1.selectbox('Image (a georeferenced RGB GeoTIFF)', found,
                       key='wi_img')
    d0 = WI.image_date(img)
    # one key per run and image: the default is the image's own day
    day = c2.date_input('Its day', value=d0 if d0 in dates else dates[-1],
                        min_value=dates[0], max_value=dates[-1],
                        key='wi_day_%s_%s' % (sel, os.path.basename(img)))
    if c3.button('Compare', key='wi_go', type='primary'):
        with st.spinner('Comparing the wetness of %s with %s...'
                        % (day, os.path.basename(img))):
            try:
                WI.compare(root, img, day, info=info, gis=mm_paths.GIS)
            except Exception as exc:                  # noqa: BLE001
                st.error('Not compared: %s' % exc)
    st.caption('Green = excess green above its Otsu threshold; 50 m blocks, '
               'the model\'s wettest blocks over the same area as the green '
               'ones; the irrigated fields (the run\'s irrigation layer) '
               'scored apart.')


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
        _wet_image(root, info)
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
