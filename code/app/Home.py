# -*- coding: utf-8 -*-
"""Panel 0 -- MARMITES / MODFLOW 6 front-end.  WP1d.

    streamlit run code/app/Home.py

A VIEWER, AN EDITOR AND A LAUNCHER. It reads the inputs, edits the
configuration, starts a run and shows the output. It is never part of the
model path: ``code/app/`` may import the model, the model may never import
``code/app/`` or streamlit, and ``code/tests/test_repo_hygiene.py`` asserts
it -- so the model keeps running headless from Spyder and from a PEST worker
where Streamlit is not installed.

The panels are in the order the work is done, and the grid comes first
because every other input is wrapped onto it.
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

import mm_paths                                  # noqa: E402
from lib import panelui, schema, runs as runlib  # noqa: E402

st.set_page_config(page_title='MARMITES / MF6', page_icon='💧', layout='wide')


def case_selector():
    cases = []
    ex = mm_paths.EXAMPLE_ROOT
    if ex.is_dir():
        cases = sorted(d.name for d in ex.iterdir()
                       if d.is_dir() and not d.name.startswith('.'))
    cases = cases or ['LaMata']
    default = st.session_state.get('case', 'LaMata')
    idx = cases.index(default) if default in cases else 0
    return st.sidebar.selectbox('Case study', cases, index=idx, key='case')


def where_things_are():
    """Panel 0 owns the machine-specific folders -- set here, once.

    They used to be read-only in the sidebar under a note telling the modeller
    to edit ``code/mm_paths.py``. A path is not source code: it says where one
    person's data sits, so it is set here and saved to a machine-local file
    that is not tracked. An ``MM_*`` environment variable still wins, and the
    form says so per row rather than leaving it to be discovered.
    """
    st.markdown('### Where things are, on this machine')
    st.caption('Set once. Written to `code/configs/paths.local.toml`, which '
               'is machine-local and not tracked — nothing here goes into the '
               'repository, and no file has to be edited by hand.')

    current = {k: str(getattr(mm_paths, {'gis': 'GIS',
                                         'example_root': 'EXAMPLE_ROOT'}.get(
                                  k, k.upper()))) for k in mm_paths.SETTABLE}
    pinned = [k for k in mm_paths.SETTABLE
              if 'environment variable' in mm_paths.source_of(k)]
    if pinned:
        st.info('Set by the environment, so they cannot be changed here: %s. '
                'That is deliberate — a batch run or a PEST worker overrides '
                'these for one run.'
                % ', '.join('`%s`' % mm_paths.SETTABLE[k][0] for k in pinned))

    edits, cols = {}, st.columns(2)
    for i, (key, (env, _default, doc)) in enumerate(mm_paths.SETTABLE.items()):
        with cols[i % 2]:
            label = key.replace('_', ' ').upper()
            got = panelui.folder_picker(
                label, current[key], key='path_%s' % key,
                help_='`%s`  \n%s  \nFrom %s.'
                      % (env, doc, mm_paths.source_of(key)),
                want='file' if key in ('nwt_ref', 'python_exe') else 'dir')
            exists = os.path.exists(got) if got else False
            st.caption(('🟢 exists' if exists else '🔴 does not exist')
                       + ('  ·  read-only: %s is set' % env
                          if key in pinned else ''))
            edits[key] = got

    c1, c2 = st.columns([1, 3])
    if c1.button('Save these paths', type='primary', key='save_paths'):
        saved = {k: v for k, v in edits.items() if k not in pinned}
        try:
            mm_paths.save_settings(saved)
        except OSError as exc:
            c2.error('Could not write %s: %r' % (mm_paths.SETTINGS, exc))
        else:
            c2.success('Saved to `%s`. Applied now — no restart.'
                       % mm_paths.settings_label())
            st.rerun()
    missing = [k for k, v in edits.items() if v and not os.path.exists(v)]
    if missing:
        c2.warning('%d path(s) do not exist yet: %s. A run needs the dataset '
                   'and the workspace; the rest are needed only by what uses '
                   'them.' % (len(missing), ', '.join(missing)))


def runs_dir(cfg=None):
    return (cfg.ui.runs_dir if cfg and cfg.ui.runs_dir else
            os.path.join(str(mm_paths.WS_ROOT), 'runs'))


def main():
    case = case_selector()
    cfg, path = panelui.pick_config()
    panelui.dataset_banner(cfg)

    st.title('💧  MARMITES / MODFLOW 6')

    _mc1, _mc2 = st.columns([2, 3])
    with _mc1:
        _edited = panelui.rows_form(cfg, (('meta.model', None),), 'meta',
                                    columns=1)
    _effective = ((_edited.get('meta.model', cfg.meta.model) or '').strip()
                  or case).lower()
    with _mc2:
        st.markdown('')
        st.caption('A soil water balance coupled to MODFLOW 6 through the '
                   'API. Case **%s**, configuration `%s`.'
                   % (case, os.path.basename(path)))
        st.caption('MODFLOW writes `%s.hds`, `%s.cbc`, `%s.lst`. Saving '
                   'renames this file to match.'
                   % (_effective, _effective, _effective))
    panelui.save_button(cfg, path, _edited, key='save_model_name')

    st.markdown('### What this configuration will run')
    cols = st.columns(3)
    for c, (sw, what) in zip(cols, [
            ('run.surface', 'MMsurf — the daily forcing'),
            ('run.model', 'MMsoil + MODFLOW 6'),
            ('run.plot', 'Figures')]):
        section, key = sw.split('.')
        on = bool(getattr(getattr(cfg, section), key))
        c.metric(what, 'ON' if on else 'off')
    if not cfg.run.surface:
        st.caption('MMsurf is off, so the daily forcing must already exist. It '
                   'is checked for presence AND shape before the run starts — '
                   'a stale or truncated file stops the run rather than being '
                   'used.')

    st.markdown('### The panels, in the order to fill them')
    # BUILT from the schema, not typed out: the table went stale the moment
    # the Model panel was split three ways, and a table of contents that
    # disagrees with the sidebar is worse than none.
    rows = ['| Panel | What it settles | Switch |', '|---|---|---|']
    for num, title, _icon, switch, _sections, blurb in schema.PANELS:
        if num == 0:
            continue
        rows.append('| **%d %s** | %s | %s |'
                    % (num, title, blurb.split('.')[0].strip().rstrip('.'),
                       '`%s`' % switch if switch else 'always'))
    rows.append('| **7 Run** | launch it, and follow the log | — |')
    rows.append('| **8 Results** | the figures a run wrote | — |')
    st.markdown(chr(10).join(rows))
    st.markdown("""

**The grid comes first**, because every other input is *wrapped onto it*: the
soil zones, the vegetation cover, the stream network and the observation
points are vector layers, projected onto whichever grid the Grid panel produced.
Change the grid and they follow — they do not have to be re-made.

The switches are not decoration. They are the same `[run]` keys the driver
reads, so turning one off means that half of the model does not execute.
""")

    where_things_are()

    ds = mm_paths.dataset_dir(case)
    c1, c2, c3 = st.columns(3)
    import marmites_meshes as _mm                 # noqa: E402
    c1.metric('Grid', '%s @ %g m' % (cfg.grid_kind,
                                     _mm.rectangle_cell_size(cfg)))
    c2.metric('Vegetation types', len(cfg.surface.vegetation))
    c3.metric('Runs on record', len(runlib.list_runs(runs_dir(cfg))))

    if not ds.is_dir():
        st.error('The dataset folder does not exist: %s' % ds)

    recent = runlib.list_runs(runs_dir(cfg), limit=5)
    if recent:
        st.markdown('### Recent runs')
        for r in recent:
            state = r.get('state')
            icon = ('⏳' if state == 'running'
                    else '✅' if r.get('outcome') == 'completed'
                    else '⛔' if r.get('outcome') == 'stopped' else '❌')
            st.write('%s `%s` — %s, started %s'
                     % (icon, r['run_id'], state, r.get('started')))

    st.markdown('---')
    st.caption('A run is launched as a DETACHED process: it survives closing '
               'this tab, and its log keeps being written. Nothing here runs a '
               'model inline — a coupled run is minutes and a calibration is '
               'hours.')


main()
