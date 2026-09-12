# -*- coding: utf-8 -*-
"""WP1d -- the widget layer shared by the four panels.

The ONLY module under ``code/app/lib`` that imports streamlit, and it does so
for one reason: every panel needs the same "field -> widget" mapping, and
duplicating it five times is how two panels end up disagreeing about what a
value means. Everything with substance still lives in ``schema.py`` and
``editor.py``, which stay importable without streamlit and are tested there.
"""

import dataclasses
import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for _p in (CODE, APP):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import marmites_config as mcfg          # noqa: E402
import mm_paths                         # noqa: E402
from lib import editor, schema          # noqa: E402

CONFIG_DIR = os.path.join(CODE, 'configs')

__all__ = ['pick_config', 'header', 'master_switch', 'section_form',
           'grid_form', 'table_form', 'save_button', 'CONFIG_DIR']


def pick_config():
    """The configuration every panel edits, remembered across pages."""
    files = sorted(f for f in os.listdir(CONFIG_DIR)
                   if f.endswith('.toml')) if os.path.isdir(CONFIG_DIR) else []
    if not files:
        st.error('No configuration in %s' % CONFIG_DIR)
        st.stop()
    default = st.session_state.get('config_file', 'lamata.toml')
    idx = files.index(default) if default in files else 0
    chosen = st.sidebar.selectbox('Configuration', files, index=idx,
                                  key='config_file')
    path = os.path.join(CONFIG_DIR, chosen)
    try:
        cfg = mcfg.load_run_config(path)
    except mcfg.ConfigError as exc:
        st.error('This configuration is not valid, so no panel can show it:\n\n'
                 '%s' % exc)
        st.stop()
    st.sidebar.caption('hash `%s`' % cfg.config_hash())
    return cfg, path


def header(number):
    """Title, blurb and the panel strip, so the order is always visible."""
    panel = next(p for p in schema.PANELS if p[0] == number)
    _n, title, icon, _sw, _sec, blurb = panel
    st.title('%s  %d — %s' % (icon, number, title))
    st.caption(blurb)
    strip = '  '.join(
        ('**%d %s**' % (p[0], p[1])) if p[0] == number else ('%d %s' % (p[0], p[1]))
        for p in schema.PANELS)
    st.caption(strip)
    return panel


def master_switch(cfg, switch):
    """The group's on/off, as the RUN key the driver actually reads."""
    if switch is None:
        return None
    section, key = switch.split('.')
    label, _u, help_ = schema.describe(switch)
    cur = bool(getattr(getattr(cfg, section), key))
    val = st.toggle('**%s**' % label, value=cur, help=help_, key='sw_%s' % switch)
    if val != cur:
        st.caption('changed — press *Validate & save* below to keep it')
    if not val:
        st.info('This group is OFF. The panel still edits its settings; the '
                'run simply will not execute it.')
    return {switch: val}


def _widget(dotted, value, prefix=''):
    """One field, with its label, units and explanation."""
    label, units, help_ = schema.describe(dotted)
    shown = '%s [%s]' % (label, units) if units else label
    key = '%s%s' % (prefix, dotted)
    help_full = '`%s`%s' % (dotted, ('  \n' + help_) if help_ else '')
    if schema.is_source(value):
        return _source_widget(dotted, value, shown, help_full, key)
    options = schema.choices_for(dotted)
    if options:
        cur = str(value)
        if cur not in options:
            options = options + [cur]
        return st.selectbox(shown, options, index=options.index(cur),
                            help=help_full, key=key)
    if isinstance(value, bool):
        return st.checkbox(shown, value=value, help=help_full, key=key)
    if isinstance(value, int) and not isinstance(value, bool):
        return st.number_input(shown, value=int(value), step=1,
                               help=help_full, key=key)
    if isinstance(value, float):
        return st.number_input(shown, value=float(value), format='%g',
                               help=help_full, key=key)
    if isinstance(value, (list, dict)):
        st.caption('%s — `%s`  \n*edit in the TOML tab*' % (shown, value))
        return None
    return st.text_input(shown, value=str(value), help=help_full, key=key)


def _source_widget(dotted, src, shown, help_full, key):
    """A ParamSource / VectorSource: the producer, then its one value.

    Shown as a choice because exactly one producer may be set -- offering
    five boxes invites setting two, which the schema then refuses.
    """
    producers = [n for n in ('raster', 'layer', 'column', 'value', 'drainage')
                 if hasattr(src, n)]
    cur = src.producer()
    st.markdown('**%s**' % shown)
    st.caption(help_full)
    idx = producers.index(cur) if cur in producers else 0
    which = st.selectbox('source', producers, index=idx,
                         key=key + '.__producer', label_visibility='collapsed')
    val = getattr(src, which)
    if which == 'drainage':
        st.caption('`%s` — edit in the TOML tab' % (val or {}))
        return {}
    if isinstance(val, float) or (val is None and which == 'value'):
        out = st.number_input('value', value=float(val if val is not None else 0.0),
                              format='%g', key=key + '.__v',
                              label_visibility='collapsed')
    else:
        out = st.text_input('value', value=str(val or ''), key=key + '.__v',
                            label_visibility='collapsed')
    # Only the chosen producer is written; the others are cleared, so the
    # precedence raster > layer > value cannot be decided by a leftover.
    edits = {}
    for n in producers:
        if n == 'drainage':
            continue
        blank = 0.0 if n == 'value' else ''
        edits['%s.%s' % (dotted, n)] = out if n == which else blank
    if which == 'value':
        edits.pop('%s.value' % dotted, None)
        edits['%s.value' % dotted] = out
    return edits


def section_form(cfg, section, columns=2, skip_subpanels=True, only=None):
    """Every plain field of a section, as widgets. Returns the edits.

    Conditional blocks (``[grid.voronoi]`` and the like) are left out by
    default and drawn by :func:`subpanel_form` under the field that controls
    them -- settings that only apply for one grid kind should not sit beside
    the ones that always apply, looking equally live.
    """
    edited = {}
    rows = schema.fields_of(cfg, section)
    if only is not None:
        rows = [r for r in rows if r[0].startswith(only)]
    elif skip_subpanels:
        rows = [r for r in rows if schema.subpanel_for(r[0]) is None]
    if not rows:
        return edited
    cols = st.columns(columns)
    for k, (dotted, value) in enumerate(rows):
        with cols[k % columns]:
            got = _widget(dotted, value)
            if isinstance(got, dict):
                edited.update(got)
            elif got is not None:
                edited[dotted] = got
    return edited


def subpanel_form(cfg, prefix, chosen, columns=3):
    """A conditional block, shown only when its controlling field selects it.

    ``chosen`` is the value the panel currently has for the controlling field
    -- the WIDGET's value, not the saved one, so the sub-panel follows the
    selector immediately rather than after a save.
    """
    rule = schema.subpanel_for(prefix)
    if rule is None:
        raise KeyError('%s is not a conditional block' % prefix)
    _control, applies = rule
    if chosen not in applies:
        return {}
    section = prefix.split('.')[0]
    label = prefix.split('.')[-1]
    st.markdown('##### `[%s]` — applies to **%s** only' % (prefix, chosen))
    return section_form(cfg, section, columns=columns, only=prefix + '.')


def native_dialog(want='dir', initial='', title='Select'):
    """The operating system's own folder/file dialog. (path or '', reason).

    Run in a SUBPROCESS, not in this one. Streamlit executes the script on a
    worker thread, and driving Tk from a non-main thread is how a server
    wedges with no window to close; a child process also cannot take the app
    down with it if Tk is missing or there is no display. It is the machine
    the SERVER runs on, which is the machine whose paths are being set --
    the same one the in-page browser walks.
    """
    import subprocess
    import sys
    code = (
        'import sys, tkinter as tk\n'
        'from tkinter import filedialog\n'
        'r = tk.Tk(); r.withdraw(); r.attributes("-topmost", True)\n'
        'p = filedialog.%s(initialdir=sys.argv[1], title=sys.argv[2])\n'
        'r.destroy()\n'
        'sys.stdout.write(p or "")\n'
        % ('askdirectory' if want == 'dir' else 'askopenfilename'))
    try:
        r = subprocess.run([sys.executable, '-c', code,
                            str(initial or os.getcwd()), str(title)],
                           capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        return '', 'the dialog was left open for five minutes'
    except OSError as exc:
        return '', '%r' % exc
    if r.returncode != 0:
        return '', (r.stderr or '').strip().splitlines()[-1:] and \
            (r.stderr or '').strip().splitlines()[-1] or 'the dialog failed'
    return (r.stdout or '').strip(), ''


def folder_picker(label, value, key, help_=None, want='dir', native=True):
    """A path: typed, browsed to in the page, or chosen in the OS dialog.

    ``…`` opens the operating system's own dialog, which is what a modeller
    expects. The in-page browser stays beside it because the dialog opens on
    the machine the SERVER runs on -- the same machine either way, but only
    one of the two works when there is no display to open a window on.

    ``want='file'`` picks a file instead of a folder.
    """
    # A widget's value can only be set BEFORE it is instantiated, so a choice
    # made by a button lower down is parked here and applied on the rerun.
    # Writing st.session_state[key] after the text_input exists raises
    # StreamlitWidgetAlreadyInstantiatedError -- which is exactly what the
    # first version of this did.
    pending = key + '.__pending'
    if pending in st.session_state:
        st.session_state[key] = st.session_state.pop(pending)

    c1, c2 = st.columns([6, 1])
    with c1:
        typed = st.text_input(label, value=str(value or ''), key=key,
                              help=help_)
    with c2:
        st.markdown('<div style="height:1.85rem"></div>',
                    unsafe_allow_html=True)
        if native and st.button('…', key=key + '.__native',
                                help='Open the system dialog'):
            start = typed if os.path.isdir(typed) else os.path.dirname(typed)
            got, why = native_dialog(want, start, 'Select %s' % label)
            if got:
                st.session_state[pending] = got
                st.rerun()
            elif why:
                st.session_state[key + '.__why'] = why
    if st.session_state.get(key + '.__why'):
        st.caption('The system dialog could not be used (%s). The browser '
                   'below does the same job.'
                   % st.session_state.pop(key + '.__why'))

    here = st.session_state.get(key + '.__at') or (
        typed if os.path.isdir(typed) else os.path.dirname(typed) or os.getcwd())
    with st.expander('Browse in the page…', expanded=False):
        st.caption('`%s`' % here)
        try:
            entries = sorted(os.listdir(here))
        except OSError as exc:
            st.error('%r' % exc)
            entries = []
        dirs = [d for d in entries if os.path.isdir(os.path.join(here, d))]
        b1, b2 = st.columns([1, 3])
        if b1.button('⬆ up', key=key + '.__up'):
            st.session_state[key + '.__at'] = \
                os.path.dirname(here.rstrip('\\/')) or here
            st.rerun()
        go = b2.selectbox('Subfolder', ['—'] + dirs, key=key + '.__sub')
        if go and go != '—':
            st.session_state[key + '.__at'] = os.path.join(here, go)
            st.session_state.pop(key + '.__sub', None)
            st.rerun()
        if want == 'file':
            files = [f for f in entries
                     if os.path.isfile(os.path.join(here, f))]
            pick = st.selectbox('File', ['—'] + files, key=key + '.__file')
            if pick and pick != '—' and st.button('Use this file',
                                                  key=key + '.__usef'):
                st.session_state[pending] = os.path.join(here, pick)
                st.rerun()
        elif st.button('Use this folder', key=key + '.__use'):
            st.session_state[pending] = here
            st.rerun()
    return typed


def boundary_picker(dotted, value, folder_key='gis_folder'):
    """Choose the catchment polygon from the shapefiles on this machine.

    A dropdown of what is actually THERE rather than a name to type: the
    commonest way this field goes wrong is a file that has been renamed or
    exported somewhere else, and a text box reports that only when the run
    fails. The folder defaults to ``DATA_ROOT/GIS`` -- where the shapefiles
    live and are read by the converter alone -- and can be pointed anywhere.

    Stores a BARE FILENAME when the file is in the GIS folder, so the
    configuration stays portable between machines, and an absolute path
    otherwise. Both are accepted everywhere the boundary is opened.
    """
    from marmites_vector import find_shapefiles

    gis = str(mm_paths.GIS)
    folder = st.text_input('Folder to look in', value=st.session_state.get(
        folder_key, gis), key=folder_key,
        help='Defaults to DATA_ROOT/GIS. The shapefiles stay here and are '
             'read only by the converter -- they never enter the repository.')
    found = find_shapefiles(folder)
    if not os.path.isdir(folder):
        st.error('No such folder: `%s`' % folder)
    elif not found:
        st.warning('No .shp in `%s` (or one level below it).' % folder)

    # What the configuration currently names, resolved the way the model
    # resolves it, so the saved value is always one of the options.
    current = os.path.join(gis, value) if value and not os.path.isabs(value) \
        else (value or '')
    options = list(found)
    if current and current not in options:
        options.insert(0, current)

    label, units, help_ = schema.describe(dotted)
    shown = '%s [%s]' % (label, units) if units else label
    if not options:
        return st.text_input(shown, value=str(value), key=dotted,
                             help='`%s`  \n%s' % (dotted, help_))
    pick = st.selectbox(
        shown, options, index=options.index(current) if current in options
        else 0, key=dotted, help='`%s`  \n%s' % (dotted, help_),
        format_func=lambda p: (os.path.relpath(p, folder)
                               if os.path.isdir(folder)
                               and p.startswith(os.path.abspath(folder))
                               else p))
    # Back to what the file should hold.
    try:
        rel = os.path.relpath(pick, gis)
    except ValueError:                         # different drive
        return pick
    return rel if not rel.startswith('..') else pick


def derived_value(cfg, dotted, edited, values):
    """A derived field's value for the settings currently ON SCREEN.

    The configuration recomputes these in ``validate()``, i.e. on save. The
    panel has to do it a step earlier, or a box that says "derived" would go
    on showing the previous corridor's bands until the modeller saved -- and
    the first thing they would do is doubt the number rather than the box.
    """
    import copy

    if dotted == 'grid.voronoi.trans_levels':
        v = copy.deepcopy(cfg.grid.voronoi)
        for name in ('cell_far', 'cell_near_stream', 'stream_buffer',
                     'grade_ratio', 'stream_refine'):
            key = 'grid.voronoi.%s' % name
            if key in edited:
                setattr(v, name, edited[key])
            elif key in st.session_state:
                setattr(v, name, st.session_state[key])
        try:
            bands = v.bands()
        except Exception:                                # noqa: BLE001
            return str(values.get(dotted, ''))
        if not bands:
            return '(none — the mesh is uniform)'
        return ', '.join('%g' % b for b in bands)
    return str(values.get(dotted, ''))


def grid_form(cfg, columns=3):
    """Panel 1's ``[grid]`` block: what is always asked, then the kind's own.

    The permanent block is the three things that do not depend on the
    producer -- the catchment, its CRS, and which producer -- plus the
    wrapping rule, which applies to every kind. Everything else belongs to
    ONE kind and is drawn underneath it, so a setting that does nothing for
    the chosen kind never sits beside one that does.

    Returns ``(edits, chosen_kind)``; the edits are keyed by dotted path
    exactly like :func:`section_form`.
    """
    edited, chosen = grid_permanent_form(cfg, columns=columns)
    values = dict(schema.fields_of(cfg, 'grid'))
    edited.update(grid_kind_form(cfg, chosen, edited, values, columns=columns))
    return edited, chosen


def grid_permanent_form(cfg, columns=3):
    """The half of ``[grid]`` that applies whatever the producer is.

    Split from the kind's own settings so the panel can put the CATCHMENT
    READ-OUTS -- area, extent, cell count -- between the two: they describe
    what has just been chosen here, and they are what says whether the file
    is the right one before any producer setting matters.
    """
    values = dict(schema.fields_of(cfg, 'grid'))
    edited = {}
    st.markdown('##### Always asked')
    cols = st.columns(columns)
    for k, dotted in enumerate(schema.GRID_PERMANENT):
        with cols[k % columns]:
            if dotted == 'grid.boundary':
                got = boundary_picker(dotted, values[dotted])
            else:
                got = _widget(dotted, values[dotted])
            if got is not None:
                edited[dotted] = got
    raw = edited.get('grid.kind', cfg.grid_kind)
    return edited, mcfg._GRID_ALIAS.get(raw, raw)


def grid_kind_form(cfg, chosen, edited=None, values=None, columns=3):
    """The settings belonging to ONE producer, and only that one."""
    values = dict(schema.fields_of(cfg, 'grid')) if values is None else values
    edited = dict(edited or {})
    out = {}
    st.markdown('##### `%s` — the settings this producer uses' % chosen)

    # A field is BLOCKED when the switch it depends on is off. Drawn greyed
    # and empty rather than hidden, so what the switch would give is visible
    # (panel 1, D4); the value itself is cleared by GridVoronoi.refresh().
    #
    # The switch's own widget is drawn INSIDE the loop below, i.e. after this,
    # so its live value comes from session_state -- where streamlit keeps the
    # widget's current value under its key -- and only falls back to the saved
    # configuration on the very first render.
    off = set()
    for switch, dependents in schema.GRID_GATED.items():
        state = st.session_state.get(switch, values.get(switch, True))
        if not edited.get(switch, state):
            off.update(dependents)

    # Row by row, so a switch sits on its own line above what it controls.
    for row in schema.GRID_SUBPANEL.get(chosen, ()):
        cols = st.columns(max(columns, len(row)))
        for k, dotted in enumerate(row):
            with cols[k]:
                label, units, help_ = schema.describe(dotted)
                shown = '%s [%s]' % (label, units) if units else label
                if dotted in schema.GRID_DERIVED:
                    # Recomputed from what is ON SCREEN, not from the saved
                    # configuration: a derived box that only caught up after a
                    # save would show the PREVIOUS corridor's bands, and the
                    # first thing doubted would be the number, not the box.
                    #
                    # It has to go through session_state. A KEYED widget takes
                    # `value` as its value on the FIRST render only and reads
                    # session_state on every one after, so passing a freshly
                    # computed `value` to a keyed box changes nothing at all.
                    ro = 'ro_%s' % dotted
                    st.session_state[ro] = derived_value(
                        cfg, dotted, dict(edited, **out), values)
                    st.text_input(shown, disabled=True, key=ro,
                                  help='`%s`  \n%s' % (dotted, help_))
                    continue
                if dotted in off:
                    st.text_input(shown, value='', disabled=True,
                                  key='off_%s' % dotted,
                                  help='`%s`  \nCleared while its switch is '
                                       'off.' % dotted)
                    continue
                got = _widget(dotted, values[dotted])
                if got is not None:
                    out[dotted] = got
    return out


def table_form(cfg, dotted, singular, path):
    """One array of tables, as an editable grid with its column help."""
    rows = editor.table_rows(cfg, dotted)
    element_help = []
    for name in (rows[0] if rows else editor.row_defaults(cfg, dotted)):
        label, units, help_ = schema.describe('%s.%s' % (dotted.split('.')[-1],
                                                         name))
        element_help.append('**%s** (`%s`%s)%s'
                            % (label, name, ', ' + units if units else '',
                               ' — ' + help_ if help_ else ''))
    st.markdown('#### %s — %d entr%s'
                % (singular, len(rows), 'y' if len(rows) == 1 else 'ies'))
    with st.expander('What each column means'):
        for line in element_help:
            st.markdown('- ' + line)

    edited = st.data_editor(rows, num_rows='dynamic', width='stretch',
                            key='tbl_%s' % dotted)
    c1, c2 = st.columns([1, 4])
    if c1.button('Save %s' % singular.lower(), key='save_%s' % dotted):
        try:
            new = editor.set_table(cfg, dotted, list(edited))
        except (editor.EditError, mcfg.ConfigError) as exc:
            c2.error('NOT saved:\n\n%s' % exc)
        else:
            new.write_toml(path)
            c2.success('Saved %d row(s) — hash %s'
                       % (len(edited), new.config_hash()))
            st.rerun()


def save_button(cfg, path, edited, label='Validate & save'):
    """The one way a panel writes. Validates first, and says what changed."""
    st.markdown('---')
    c1, c2 = st.columns([1, 3])
    if not c1.button(label, type='primary', key='save_%s' % id(edited)):
        return
    try:
        applied, digest = editor.save(cfg, path, edited)
    except (editor.EditError, mcfg.ConfigError) as exc:
        c2.error('NOT saved — the configuration would be invalid:\n\n%s' % exc)
        return
    if not applied:
        c2.info('Nothing changed.')
        return
    c2.success('Saved %d change(s) — hash %s' % (len(applied), digest))
    with c2.expander('What changed'):
        for a in applied:
            st.code(a, language='ini')
    st.rerun()


def dataset_banner(cfg):
    """Where this panel is actually reading from, on this machine."""
    ds = mm_paths.dataset_dir(cfg.paths.case)
    rows = [('dataset', ds), ('GIS (converter only)', mm_paths.GIS),
            ('workspace', mm_paths.WS_ROOT)]
    for label, p in rows:
        ok = os.path.exists(str(p))
        st.sidebar.markdown('%s **%s**  \n`%s`' % ('🟢' if ok else '🔴', label, p))
    return ds
