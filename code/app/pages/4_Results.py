# -*- coding: utf-8 -*-
"""Results — the figures a run wrote.  WP1b, step 1b.8.

SHOW, DO NOT REDRAW. The native suite already writes every figure into
`out_<stamp>_<tag>/_input` and `_output`; this page presents them with a run
picker and a filter. One source of truth per figure, the same rule that
produced the suite in the first place. Interactivity is added only where it
buys something, and that starts in stage 2 (WP6).
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

import mm_paths                       # noqa: E402

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
for sub in ('_input', '_output', 'figures_nwt_comparison'):
    d = os.path.join(root, sub)
    if os.path.isdir(d):
        for f in sorted(os.listdir(d)):
            if f.lower().endswith('.png'):
                pngs.append((sub, f, os.path.join(d, f)))
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

folders = sorted({s for s, _f, _p in pngs})
pick = st.multiselect('Folder', folders, default=folders)
needle = st.text_input('Filter by name', '')
shown = [(s, f, p) for s, f, p in pngs
         if s in pick and needle.lower() in f.lower()]
st.caption('%d of %d figure(s)' % (len(shown), len(pngs)))

ncol = st.slider('Columns', 1, 4, 2)
cols = st.columns(ncol)
for k, (sub, fname, fpath) in enumerate(shown):
    with cols[k % ncol]:
        st.image(fpath, caption='%s/%s' % (sub, fname), use_container_width=True)

if csvs:
    with st.expander('CSV output (%d)' % len(csvs)):
        for fname, fpath in csvs:
            st.write('`%s`' % fname)
            st.download_button('download %s' % fname,
                               data=open(fpath, 'rb').read(),
                               file_name=fname, key='dl_' + fname)
