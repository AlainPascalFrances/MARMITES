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

    current = {k: getattr(mm_paths, {'gis': 'GIS', 'rs': 'RS',
                                     'example_root': 'EXAMPLE_ROOT'}.get(
                              k, k.upper())) for k in mm_paths.SETTABLE}
    # an optional root left unset is None, which must show as a blank box
    current = {k: '' if v is None else str(v) for k, v in current.items()}
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
            if key == 'pest_ws':
                pest_ws_warning(got)
            if key == 'rs' and exists:
                rs_products(got)
            if key == 'soil_moisture' and exists:
                sm_products(got)
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


def rs_products(path):
    """What the remote-sensing folder holds: one line per product folder
    and the files under it, so a download that stopped short shows here."""
    try:
        subs = sorted(d for d in os.listdir(path)
                      if os.path.isdir(os.path.join(path, d)))
    except OSError:
        return
    if not subs:
        st.caption('No product folder yet: download each product into a '
                   'folder of its own, e.g. `MODIS4LaMata`.')
        return
    rows = []
    for d in subs:
        n = 0
        for _root, _dirs, files in os.walk(os.path.join(path, d)):
            n += sum(1 for f in files if f.lower().endswith(('.hdf', '.tif',
                                                             '.nc', '.h5')))
        rows.append('`%s` (%d file%s)' % (d, n, '' if n == 1 else 's'))
    st.caption('Products: ' + ', '.join(rows))


def sm_products(path):
    """What the soil-moisture folder holds: per product folder its files, or
    the layers of its Copernicus Browser exports; a folder that is not one of
    the products (ppMF6/marmites_smrs) is named as such."""
    import re
    ppmf6 = os.path.join(CODE, 'ppMF6')
    if ppmf6 not in sys.path:
        sys.path.insert(0, ppmf6)
    try:
        import marmites_smrs as SMRS
        known = list(SMRS.PRODUCTS)
    except Exception:                                    # noqa: BLE001
        known = []
    try:
        subs = sorted(d for d in os.listdir(path)
                      if os.path.isdir(os.path.join(path, d)))
    except OSError:
        return
    rows = []
    for d in subs:
        files = [f for _r, _d, fs in os.walk(os.path.join(path, d)) for f in fs
                 if f.lower().endswith(('.nc', '.nc4', '.tif', '.tiff',
                                        '.csv'))]
        csv = [f for f in files if f.lower().endswith('.csv')]
        layers = sorted({m.group(1).upper() for m in
                         (re.search(r'-(SWI_?\d+|SSM)-', f, re.I) for f in csv)
                         if m})
        what = ('Browser exports: %s' % (', '.join(layers) or
                                         '%d CSV' % len(csv)) if csv else
                '%d file%s' % (len(files), '' if len(files) == 1 else 's'))
        if known and d not in known:
            what += ', not a product folder'
        rows.append('`%s` (%s)' % (d, what))
    if not rows:
        st.caption('No product folder yet: one per product, named %s.'
                   % ', '.join('`%s`' % k for k in known))
        return
    st.caption('Products: ' + '; '.join(rows))


def pest_ws_warning(path):
    """The PEST workspace needs room: say how much, and what the drive has.

    A calibration keeps a model copy per worker and every realisation's
    outputs -- CdL's PEST folder is ~40 GB -- so the warning comes with the
    free space of the drive chosen, and a folder inside WS_ROOT or the
    repository is refused in words (a PEST run must never write into the
    workspace the modeller's own runs use)."""
    import shutil
    note = ('**Needs tens of GB.** The calibration copies the model for every '
            'worker and keeps each realisation\'s outputs: CdL\'s PEST folder '
            '(`X:\\3p1p1\\DR3PDA1`) is ~40 GB.')
    if not path:
        st.warning(note + ' Not set: no calibration runs until it is.')
        return
    probe = path
    while probe and not os.path.exists(probe):
        parent = os.path.dirname(probe)
        probe = parent if parent != probe else ''
    free = ''
    if probe:
        try:
            free = ' Free on that drive now: **%.0f GB**.' % (
                shutil.disk_usage(probe).free / 1e9)
        except OSError:
            free = ''
    st.warning(note + free)
    real = os.path.normcase(os.path.abspath(path))
    for name, root in (('WS_ROOT', mm_paths.WS_ROOT),
                       ('the repository', mm_paths.REPO)):
        r = os.path.normcase(os.path.abspath(str(root)))
        if real == r or real.startswith(r + os.sep):
            st.error('It is inside %s (`%s`): choose a folder of its own -- '
                     'a calibration must not write where the model runs or '
                     'the code live.' % (name, root))


def runs_dir(cfg=None):
    return (cfg.ui.runs_dir if cfg and cfg.ui.runs_dir else
            os.path.join(str(mm_paths.WS_ROOT), 'runs'))


def main():
    case = case_selector()
    cfg, path = panelui.pick_config()
    panelui.dataset_banner(cfg)

    st.title('💧  MARMITES / MODFLOW 6')
    st.caption('A soil water balance coupled to MODFLOW 6 through the API.')

    _mc1, _mc2 = st.columns([2, 3])
    with _mc1:
        _edited = panelui.rows_form(cfg, (('meta.model', None),), 'meta',
                                    columns=1)
    _effective = ((_edited.get('meta.model', cfg.meta.model) or '').strip()
                  or case).lower()
    with _mc2:
        st.markdown('')
        st.caption('Case **%s**, configuration `%s`.'
                   % (case, os.path.basename(path)))
        st.caption('MODFLOW writes `%s.hds`, `%s.cbc`, `%s.lst`. Saving '
                   'renames this file to match.'
                   % (_effective, _effective, _effective))
    panelui.remember(_edited)

    # THE FILE, which is what a run reads -- not the toggles, which are
    # widgets until a panel saves them. Where the two disagree the metric
    # says so rather than showing one of them and leaving the other to be
    # discovered by a run.
    st.markdown('### What this configuration will run')
    _pending = dict((sw, (live, saved))
                    for sw, live, saved in panelui.unsaved_switches(cfg))
    cols = st.columns(3)
    for c, (sw, what) in zip(cols, [
            ('run.surface', 'MMsurf — the daily forcing'),
            ('run.model', 'MMsoil + MODFLOW 6'),
            ('run.plot', 'Figures')]):
        section, key = sw.split('.')
        on = bool(getattr(getattr(cfg, section), key))
        c.metric(what, 'ON' if on else 'off',
                 delta=('panel says %s — NOT SAVED'
                        % ('ON' if _pending[sw][0] else 'off')
                        if sw in _pending else None),
                 delta_color='inverse' if sw in _pending else 'normal')
    if _pending:
        st.error('A switch was changed on a panel and not saved. **The run '
                 'reads the file**, so it would ignore the change. Save the '
                 'panel it was changed on.')
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
    rows.append('| **7 Run** | validate the configuration (compulsory), '
                'then launch it and follow the log | — |')
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

    # The one save, as on every panel. Last, so it counts what this page has
    # remembered on this run.
    panelui.sidebar_save(cfg, path)


main()
