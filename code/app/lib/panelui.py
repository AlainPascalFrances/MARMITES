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
           'grid_form', 'grid_permanent_form', 'grid_kind_form',
           'gis_folder_box', 'layer_picker', 'resolve_layer', 'folder_picker',
           'table_form', 'save_button', 'CONFIG_DIR']


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
    """A path: typed, or chosen in the operating system's own dialog.

    ``…`` opens the OS dialog, which is what a modeller expects. There is no
    in-page folder browser beside it any more: with the dialog working it was
    two ways to do one thing, and the text box is the fallback when the
    dialog cannot open a window.

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
        st.caption('The system dialog could not be used (%s) — type the path '
                   'instead.' % st.session_state.pop(key + '.__why'))
    return typed


def gis_folder_box(folder_key='gis_folder'):
    """Where to look for shapefiles. Shared by the layer pickers below."""
    from marmites_vector import find_shapefiles

    gis = str(mm_paths.GIS)
    folder = st.text_input(
        'Folder to look in', value=st.session_state.get(folder_key, gis),
        key=folder_key,
        help='Defaults to DATA_ROOT/GIS. The shapefiles stay here and are '
             'read only by the converter -- they never enter the repository.')
    found = find_shapefiles(folder)
    if not os.path.isdir(folder):
        st.error('No such folder: `%s`' % folder)
    elif not found:
        st.warning('No .shp in `%s` (or one level below it).' % folder)
    return folder, found


NONE = '\u2014 none \u2014'


def layer_picker(dotted, value, folder, found, optional=False):
    """Choose one shapefile from what is actually in the folder.

    A dropdown of what is THERE rather than a name to type, because the
    commonest way these fields go wrong is a file renamed or exported
    somewhere else, and a text box reports that only when the run fails.

    Stores a BARE FILENAME while the file is in the GIS folder, so the
    configuration stays portable between machines, and an absolute path
    otherwise -- both are accepted everywhere a layer is opened. ``optional``
    adds a "none" entry, which is what a catchment with no mapped network or
    no ponds selects, and which is the DEFAULT for a new case.
    """
    gis = str(mm_paths.GIS)
    label, units, help_ = schema.describe(dotted)
    shown = '%s [%s]' % (label, units) if units else label

    current = os.path.join(gis, value) if value and not os.path.isabs(value) \
        else (value or '')
    options = list(found)
    if current and current not in options:
        options.insert(0, current)
    if optional:
        options = [NONE] + options
    if not options:
        return st.text_input(shown, value=str(value), key=dotted,
                             help='`%s`  \n%s' % (dotted, help_))

    pick = st.selectbox(
        shown, options,
        index=options.index(current) if current in options else 0,
        key=dotted, help='`%s`  \n%s' % (dotted, help_),
        format_func=lambda p: (p if p == NONE else
                               (os.path.relpath(p, folder)
                                if os.path.isdir(folder)
                                and p.startswith(os.path.abspath(folder))
                                else p)))
    if pick == NONE:
        return ''
    try:
        rel = os.path.relpath(pick, gis)
    except ValueError:                         # a different drive
        return pick
    return rel if not rel.startswith('..') else pick


def resolve_layer(name):
    """A configured layer name as an absolute path, or '' when unset."""
    if not name:
        return ''
    return (name if os.path.isabs(name)
            else os.path.join(str(mm_paths.GIS), name))


def derived_value(cfg, dotted, edited, values):
    """A derived field's value for the settings currently ON SCREEN.

    The configuration recomputes these in ``validate()``, i.e. on save. The
    panel has to do it a step earlier, or a box that says "derived" would go
    on showing the previous corridor's bands until the modeller saved -- and
    the first thing they would do is doubt the number rather than the box.
    """
    import copy

    if dotted.startswith('grid.voronoi.'):
        v = live_voronoi(cfg, edited)
        try:
            bands = v.graded_bands()
        except Exception:                                # noqa: BLE001
            return str(values.get(dotted, ''))
        if dotted == 'grid.voronoi.stream_buffer':
            return ('%g' % bands[-1][0]) if bands else '0'
        if dotted == 'grid.voronoi.trans_levels':
            if not bands:
                return '(none — the mesh is uniform)'
            return ', '.join('%g' % d for d, _s in bands)
    return str(values.get(dotted, ''))


def live_voronoi(cfg, edited=None):
    """``cfg.grid.voronoi`` with the values currently ON SCREEN applied.

    A copy, so nothing is written by drawing a page. The derived corridor and
    bands are computed from this: the configuration only recomputes them on
    save, and a box that caught up only then would show the previous
    corridor's numbers.
    """
    import copy

    v = copy.deepcopy(cfg.grid.voronoi)
    edited = edited or {}
    for name in ('cell_far', 'cell_near_stream', 'grade_ratio',
                 'stream_refine'):
        key = 'grid.voronoi.%s' % name
        if key in edited:
            setattr(v, name, edited[key])
        elif key in st.session_state:
            setattr(v, name, st.session_state[key])
    return v


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
    edited, chosen, _reps = grid_permanent_form(cfg, columns=columns)
    values = dict(schema.fields_of(cfg, 'grid'))
    edited.update(grid_kind_form(cfg, chosen, edited, values, columns=columns))
    return edited, chosen


def grid_permanent_form(cfg, columns=3):
    """The half of ``[grid]`` that applies whatever the producer is.

    Returns ``(edits, chosen_kind, report)``. The catchment is CHECKED here,
    once, and the report handed back: the CRS shown below the file comes from
    that check, and so do the read-outs the page draws underneath -- reading
    the shapefile twice per render to tell the same thing twice would be one
    read too many.
    """
    from marmites_vector import check_polygon_layer

    values = dict(schema.fields_of(cfg, 'grid'))
    edited, reports = {}, {}

    # --- the three layers a GRID can depend on -------------------------
    # All asked HERE, at the same level, because the refinement options
    # cannot be answered before it is known whether there IS a network or a
    # pond to refine around -- and because a new catchment starts with an
    # empty list and has to be able to say so.
    folder, found = gis_folder_box()
    cols = st.columns(3)
    with cols[0]:
        picked = layer_picker('grid.boundary', values['grid.boundary'],
                              folder, found)
    with cols[1]:
        streams = layer_picker('grid.streams', values['grid.streams'],
                               folder, found, optional=True)
    with cols[2]:
        ponds = layer_picker('grid.ponds', values['grid.ponds'],
                             folder, found, optional=True)
    edited['grid.boundary'] = picked
    edited['grid.streams'] = streams
    edited['grid.ponds'] = ponds

    path = resolve_layer(picked)
    rep = check_polygon_layer(path, expect_epsg=cfg.grid.crs_epsg)
    reports['grid.boundary'] = rep
    # The other two are checked AGAINST the catchment: same CRS, and an
    # extent that actually overlaps it. Nothing reprojects, so a layer from
    # another catchment or in another CRS is a layer the mesh refines nowhere.
    for dotted, name, want in (('grid.streams', streams, 'line'),
                               ('grid.ponds', ponds, 'polygon')):
        if not name:
            continue
        reports[dotted] = check_polygon_layer(
            resolve_layer(name), expect_epsg=(rep['epsg']
                                              or cfg.grid.crs_epsg),
            against=rep['bbox'], want=want)

    # Read-only: the CRS is a PROPERTY OF THE FILE, not a choice. It is taken
    # from the .prj -- by its authority code, or by matching the definition
    # when ArcGIS wrote none -- and only falls back to the configured value
    # when the file says nothing a code can be made of.
    epsg = int(rep['epsg'] or cfg.grid.crs_epsg or 0)
    label, units, _h = schema.describe('grid.crs_epsg')
    c1, _c2 = st.columns(2)
    with c1:
        st.text_input('%s [%s]' % (label, units), disabled=True,
                      key='ro_grid.crs_epsg',
                      value=('EPSG:%d' % epsg) if epsg else 'UNDECLARED',
                      help='`grid.crs_epsg`  \nRead from the layer, not '
                           'chosen: every other layer is assumed to be in '
                           'this CRS and nothing is reprojected.')
        st.caption(
            'from the `.prj`%s'
            % (' (matched — it carries no EPSG code)' if rep.get('epsg_matched')
               else '' if rep['epsg'] else
               ': it declares none, so the configured value stands'))
    if epsg and epsg != int(cfg.grid.crs_epsg):
        edited['grid.crs_epsg'] = epsg

    # --- the producer, and how values are carried onto it --------------
    cols = st.columns(columns)
    for k, dotted in enumerate(('grid.kind', 'grid.resample')):
        with cols[k]:
            got = _widget(dotted, values[dotted])
            if got is not None:
                edited[dotted] = got
    raw = edited.get('grid.kind', cfg.grid_kind)
    return edited, mcfg._GRID_ALIAS.get(raw, raw), reports


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
    # A switch whose LAYER has not been given cannot be answered at all: it
    # is drawn off and disabled, with the reason, rather than offered and
    # then refused by validate(). Read from the edits, so choosing the layer
    # enables it immediately instead of after a save.
    unavailable = {}
    for switch, needs in schema.GRID_NEEDS_LAYER.items():
        have = edited.get(needs, values.get(needs, ''))
        if not have:
            unavailable[switch] = needs

    off = set()
    for switch, dependents in schema.GRID_GATED.items():
        state = st.session_state.get(switch, values.get(switch, True))
        if switch in unavailable or not edited.get(switch, state):
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
                if dotted in unavailable:
                    st.checkbox(shown, value=False, disabled=True,
                                key='na_%s' % dotted,
                                help='`%s`  \nNeeds `%s`, which is not set.'
                                     % (dotted, unavailable[dotted]))
                    st.caption('needs a %s layer'
                               % unavailable[dotted].split('.')[-1])
                    out[dotted] = False
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
