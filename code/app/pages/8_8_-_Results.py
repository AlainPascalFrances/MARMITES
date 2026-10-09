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
"""

import html
import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

import mm_paths                       # noqa: E402
from lib import results as R          # noqa: E402

st.set_page_config(page_title='Results', page_icon='📊', layout='wide')
WS = str(mm_paths.WS_ROOT)

st.title('Results')
st.caption('The figures the run wrote, from `%s`.' % WS)


@st.cache_data(show_spinner=False)
def _out_dirs(ws, _key):
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
    st.warning('No figures in this run folder. Re-draw them without re-running '
               'the model:\n\n`python code/tests/run_lamata_mf6.py --config '
               'code/configs/lamata.toml --set postproc.only=true --run-tag %s`'
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
    """Two figures to a row, aligned on their tops; an odd last one centred
    at the same width."""
    for k in range(0, len(items), 2):
        pair = items[k:k + 2]
        cols = (st.columns(2, gap='medium', vertical_alignment='top')
                if len(pair) == 2 else [st.columns([1, 2, 1])[1]])
        for col, (sub, fname, path, text) in zip(cols, pair):
            with col:
                st.image(path, width='stretch')
                _caption(sub, fname, text)


keys = [k for k, _l in R.TABS if groups.get(k)
        or (k == 'budget' and os.path.exists(summary) and not needle)]
tabs = st.tabs(['%s (%d)' % (dict(R.TABS)[k], len(groups.get(k, [])))
                for k in keys])
for key, tab in zip(keys, tabs):
    with tab:
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
