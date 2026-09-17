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

__all__ = ['pick_config', 'header', 'master_switch', 'switch_and_save',
           'section_form',
           'grid_form', 'grid_permanent_form', 'grid_kind_form',
           'gis_folder_box', 'layer_picker', 'resolve_layer', 'folder_picker',
           'table_form', 'save_button', 'park', 'live',
           'surface_folder_box', 'file_picker', 'path_box',
           'rows_form', 'read_only', 'column_box', 'how_box',
           'record_lines',
           'CONFIG_DIR']


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
    # Remembered here because it is the ONE place every panel passes
    # through, and the source widgets below need it without being handed a
    # configuration they otherwise have no use for.
    st.session_state['__dataset_dir'] = str(
        mm_paths.dataset_dir(cfg.paths.case))
    st.sidebar.caption('hash `%s`' % cfg.config_hash())
    for line in getattr(cfg, 'migrated', ()):
        # Said once, where the file is chosen: the name in the file is not
        # the name in the panel until something saves it back.
        st.sidebar.info('This file uses an old name — %s. Saving from any '
                        'panel writes the new one.' % line)
    return cfg, path


def header(number):
    """Title, blurb and the panel strip, so the order is always visible."""
    panel = next(p for p in schema.PANELS if p[0] == number)
    _n, title, icon, _sw, _sec, blurb = panel
    st.title('%s  %d — %s' % (icon, number, title))
    st.caption(blurb)
    return panel


def master_switch(cfg, switch):
    """The group's on/off, as the RUN key the driver actually reads.

    ONE switch can appear on more than one panel -- ``run.model`` is the soil
    water balance and MODFLOW together, which are no longer separable -- and
    then the two have to agree. A widget's own state does not survive moving
    to another page, so the value is kept under a plain session key as well,
    and that is what the next panel starts from. Saving on either writes the
    same field, so they cannot drift.
    """
    if switch is None:
        return None
    section, key = switch.split('.')
    label, _u, help_ = schema.describe(switch)
    saved = bool(getattr(getattr(cfg, section), key))
    live_key = 'live_%s' % switch
    cur = bool(st.session_state.get(live_key, saved))
    val = st.toggle('**%s**' % label, value=cur, help=help_,
                    key='sw_%s' % switch)
    st.session_state[live_key] = bool(val)
    if val != saved:
        st.caption('changed — press *Validate & save* below to keep it')
    if not val:
        st.info('This group is OFF. The panel still edits its settings; the '
                'run simply will not execute it.')
    sharing = [p for p in schema.PANELS if p[3] == switch]
    if len(sharing) > 1:
        st.caption('One switch for panel%s %s — they run, or do not run, '
                   'together.'
                   % ('s' if len(sharing) > 1 else '',
                      ' and '.join('**%d %s**' % (p[0], p[1])
                                   for p in sharing)))
    return {switch: val}


def park(key, value):
    """Set a widget's value from a button drawn BELOW it, on the next run.

    A keyed widget takes ``value=`` on its FIRST render only and reads
    session_state on every one after, and session_state[key] cannot be
    assigned once the widget exists. So a button further down the page parks
    what it wants here and the widget picks it up when the page reruns.
    """
    st.session_state[key + '.__pending'] = value


def live(key, default=None):
    """What a widget will hold on THIS run -- a parked value included.

    Anything deciding whether ANOTHER widget is drawn has to ask this rather
    than session_state: the parked value is only moved across when the widget
    itself is instantiated, which for a switch happens after the fields it
    gates have already been placed. Reading session_state alone left a switch
    turned on by a button with its own dependants still greyed for one run --
    and the values parked into them unconsumed.
    """
    pending = key + '.__pending'
    if pending in st.session_state:
        return st.session_state[pending]
    return st.session_state.get(key, default)


def _widget(dotted, value, prefix=''):
    """One field, with its label, units and explanation."""
    label, units, help_ = schema.describe(dotted)
    shown = '%s [%s]' % (label, units) if units else label
    key = '%s%s' % (prefix, dotted)
    # Anything parked by a button lower down the page (see park) is applied
    # HERE, before the widget is instantiated, which is the only moment it
    # can be. A widget given its value that way must NOT also be given a
    # default: streamlit warns that two things are deciding it -- and it is
    # right, the session_state one silently wins.
    parked = key + '.__pending' in st.session_state
    if parked:
        st.session_state[key] = st.session_state.pop(key + '.__pending')

    def _default(**kw):
        return {} if parked else kw

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
        return st.checkbox(shown, help=help_full, key=key,
                           **_default(value=value))
    if isinstance(value, int) and not isinstance(value, bool):
        return st.number_input(shown, step=1, help=help_full, key=key,
                               **_default(value=int(value)))
    if isinstance(value, float):
        return st.number_input(shown, format='%g', help=help_full, key=key,
                               **_default(value=float(value)))
    if isinstance(value, (list, dict)):
        st.caption('%s — `%s`  \n*edit in the TOML tab*' % (shown, value))
        return None
    return st.text_input(shown, help=help_full, key=key,
                         **_default(value=str(value)))


# A producer that names a FILE, and the folder that kind of file lives in.
# By RULE, not by listing the fields one at a time: a source is a source, and
# the next one added should get the dialog without anyone remembering to add
# it here. `column` names an attribute and `value` is a number, so neither is
# a file and neither gets one.

def _producer_folder(which):
    """Where to start the dialog for this producer, or None if not a file."""
    if which == 'layer':
        return str(mm_paths.GIS)          # cartography, read by the converter
    if which == 'raster':
        # The DATASET copy, grid-independent, written by the converter -- not
        # the GIS original.
        return st.session_state.get('__dataset_dir', '')
    return None


def layer_box(dotted, current, nlay, key=None):
    """Which MODFLOW layers a boundary applies to, 1-based.

    Counted from 1 because that is how MODFLOW, the parameter file and every
    conversation about the model count them; the run subtracts one when it
    indexes an array. Layers the model no longer has are dropped from the
    default rather than silently kept, so shrinking nlay cannot leave a
    boundary pointing at a layer that is not there.
    """
    label, units, help_ = schema.describe(dotted)
    shown = '%s [%s]' % (label, units) if units else label
    options = list(range(1, max(int(nlay), 0) + 1))
    default = [int(v) for v in (current or []) if int(v) in options]
    return st.multiselect(shown, options, default=default,
                          key=key or dotted,
                          help='`%s`  \n%s' % (dotted, help_),
                          format_func=lambda n: 'layer %d' % n)


def boundary_note(cfg, name, value_field, dataset_dir):
    """What this boundary package resolves to, and whether the run reads it.

    The same two things the layers tab says, because they are the same two
    questions: which files a `%d` pattern opens on the layers the boundary
    applies to, and whether answering here changes what a run does. A tab
    that showed neither would be asking the modeller to take it on trust.
    """
    pkg = getattr(cfg, name)
    if not pkg.enable:
        st.caption('Off — `ModflowGwf%s` is not built and nothing above is '
                   'read.' % name)
        return
    rows = []
    for field in (value_field, 'cond'):
        src = getattr(pkg, field)
        if src.producer() != 'raster' or '%d' not in (src.raster or ''):
            continue
        marks = []
        for L in pkg.layers:
            fn = src.raster % int(L)
            there = _spelled_here(os.path.join(str(dataset_dir or ''), fn))
            marks.append('%s `%s`' % ('🟢' if there else '🔴',
                                      fn))
        rows.append('- **%s** — %s' % (field, ', '.join(marks)))
    if rows:
        st.markdown('On layer(s) %s:\n\n%s'
                    % (', '.join(str(L) for L in pkg.layers),
                       '\n'.join(rows)))
    st.success('The run builds `[%s]` from these answers, not from '
               '`MF_ws/__inputMF_flopy_v3_*.ini`. A cell carries a boundary '
               'where the source above produces a value; everywhere else it '
               'falls back to `fill` and there is none.' % name)


def _spelled_here(path):
    """Is the file there, spelled exactly as asked? (case-strict)"""
    try:
        return os.path.basename(path) in os.listdir(os.path.dirname(path))
    except OSError:
        return False


def label_of(dotted):
    """The field's label, for a box whose own label is collapsed."""
    label, units, _h = schema.describe(dotted)
    return '%s [%s]' % (label, units) if units else label


def how_box(dotted, current, key, layer_name):
    """How the layer's values are carried onto a cell.

    The modes that apply depend on the GEOMETRY -- `length` means nothing for
    a polygon and `area_mean` nothing for a line -- so the list is built from
    the layer's own header rather than offering everything and letting the
    run refuse it. `auto` is the default and picks by kind.
    """
    from marmites_vector import layer_kind, overlay_modes_for

    try:
        kind = layer_kind(resolve_layer(str(layer_name or '')))
    except Exception:                       # never break the panel for this
        kind = ''
    modes = list(overlay_modes_for(kind))
    cur = str(current or 'auto')
    if cur not in modes:
        modes.append(cur)
    st.caption('Overlay rule%s' % (' (%s layer)' % kind if kind else ''))
    return st.selectbox(
        'Overlay rule', modes, index=modes.index(cur),
        key=key + '.__how', label_visibility='collapsed',
        help='`%s.how`  \nHow a cell takes its value from the features that '
             'touch it. auto picks by geometry: majority for a polygon class '
             'map, longest for a line. majority is the class covering most '
             'of the cell, area_fraction the percentage covered, area_mean '
             'the area-weighted mean of the column, presence 1 where '
             'anything touches. For lines, longest is the feature with most '
             'length in the cell and length the metres of it.' % dotted)


def _layer_modifiers(dotted, src, key, layer_name):
    """The column and the overlay rule, side by side under their layer.

    Both are properties OF the layer named above them, not alternatives to
    it, and both are read off that layer's own header.
    """
    c1, c2 = st.columns(2)
    with c1:
        col = column_box(dotted, src.column, key + '.__col',
                         layer_name)
    with c2:
        how = how_box(dotted, getattr(src, 'how', 'auto'), key, layer_name) \
            if hasattr(src, 'how') else None
    return col, how


def column_box(dotted, current, widget_key, layer_name):
    """The attribute column, chosen from what the shapefile actually carries.

    A LIST, not a box to type into: the names are in the file, so asking
    someone to remember GRID_CODE against GRIDCODE is asking them to make a
    mistake the panel could have prevented. Read from the .dbf header alone,
    so offering them costs a few hundred bytes rather than the whole layer.

    A layer that is not there yet, or has no attributes, falls back to a text
    box -- a name typed for a file still to be exported should not be thrown
    away, and a picker with nothing in it is worse than a box.
    """
    from marmites_vector import layer_fields

    hint = ('`%s.column`  \nThe attribute of that layer carrying the value. '
            'None uses the layer\'s own geometry -- a zone per polygon, in '
            'file order.' % dotted)
    cur = str(current or '')
    try:
        fields = layer_fields(resolve_layer(str(layer_name or '')))
    except Exception:                       # never break the panel for this
        fields = []
    st.caption('Attribute column')
    if not fields:
        return st.text_input('Attribute column', value=cur,
                             key=widget_key, help=hint,
                             label_visibility='collapsed')
    # A column the file no longer has is SHOWN rather than dropped: it is
    # what the configuration says, and silently replacing it with the first
    # attribute in the file would be a change nobody asked for.
    options = [NONE] + fields + ([cur] if cur and cur not in fields else [])
    pick = st.selectbox('Attribute column', options,
                        index=options.index(cur) if cur in options else 0,
                        key=widget_key, help=hint,
                        label_visibility='collapsed')
    if cur and cur not in fields:
        st.caption('⚠️ `%s` is not an attribute of that layer.' % cur)
    return '' if pick == NONE else pick


def _source_widget(dotted, src, shown, help_full, key):
    """A ParamSource / VectorSource: the producer, then its one value.

    Shown as a choice because exactly one producer may be set -- offering
    five boxes invites setting two, which the schema then refuses.

    `column` is a producer for a ParamSource, where it means a column of the
    SFR source layer. It is NOT one for a VectorSource: there it is an
    attribute OF the layer named beside it, and `producer()` never returns
    it. Offering it as a choice there set the column and cleared the layer,
    leaving a source that produces nothing -- so for those it is drawn UNDER
    the layer box instead, which is where it belongs.
    """
    layered = hasattr(src, 'layer')
    producers = [n for n in ('raster', 'layer', 'column', 'value', 'drainage')
                 if hasattr(src, n) and not (layered and n == 'column')]
    cur = src.producer()
    # The field's name and its explanation ride on the title as a tooltip,
    # the way the layer pickers do. As body text they were three lines of
    # grey under every one of these boxes -- and on a panel that has a dozen
    # of them, the text is most of the panel.
    st.markdown('**%s**' % shown, help=help_full)
    idx = producers.index(cur) if cur in producers else 0
    which = st.selectbox('source', producers, index=idx,
                         key=key + '.__producer', label_visibility='collapsed')
    val = getattr(src, which)
    if which == 'drainage':
        st.caption('`%s` — edit in the TOML tab' % (val or {}))
        return {}
    col = how = None
    folder = _producer_folder(which)
    if folder is not None:
        # It names a FILE -- a shapefile in the cartography folder, a raster
        # in the dataset -- so it is chosen the way every other file on these
        # panels is chosen: typed if you know the name, picked if you do not.
        out = path_box(key + '.__v', label_of(dotted), val, folder,
                       help_=help_full, collapsed=True)
    elif which == 'value' and dotted in schema.INTEGER_VALUE:
        # A zone is a NUMBER of a zone: 1.0 zones is not a thing, and a box
        # that offers decimals invites one. Read defensively -- a file edited
        # by hand, or by an older version of this panel, may hold a string.
        try:
            start = int(float(val))
        except (TypeError, ValueError):
            start = 0
        out = st.number_input('value', value=start, step=1,
                              key=key + '.__v', label_visibility='collapsed')
    elif isinstance(val, float) or (val is None and which == 'value'):
        out = st.number_input('value', value=float(val if val is not None else 0.0),
                              format='%g', key=key + '.__v',
                              label_visibility='collapsed')
    else:
        out = st.text_input('value', value=str(val or ''), key=key + '.__v',
                            label_visibility='collapsed')
    if layered and which == 'layer':
        col, how = _layer_modifiers(dotted, src, key, out)

    # Only the chosen producer is written; the others are CLEARED, so the
    # precedence raster > layer > value cannot be decided by a leftover.
    #
    # `value` is the exception, and clearing it was a bug: it was written as
    # 0.0, which is a legitimate zone NUMBER rather than an absence, and the
    # field's default is None -- so the editor had no type to coerce to and
    # stored the string "0.0". It is written only when it is the producer.
    # Leaving an old number behind is harmless: raster and layer both beat it.
    edits = {}
    for n in producers:
        if n in ('drainage', 'value'):
            continue
        edits['%s.%s' % (dotted, n)] = out if n == which else ''
    if which == 'value':
        # The FIELD is a float even where the box is a count, so the file
        # keeps one type and a saved 1 does not become an int on one machine
        # and a float on the next.
        edits['%s.value' % dotted] = float(out)
    if layered:
        # Written only with the layer it qualifies; cleared with it, or a
        # column left behind would describe a layer that is no longer named.
        edits['%s.column' % dotted] = col if which == 'layer' else ''
        if how is not None:
            edits['%s.how' % dotted] = how
    return edits


def section_form(cfg, section, columns=2, skip_subpanels=True, only=None):
    """Every plain field of a section, as widgets. Returns the edits.

    Conditional blocks (``[grid.voronoi]`` and the like) are left out by
    default and drawn by :func:`subpanel_form` under the field that controls
    them -- settings that only apply for one grid kind should not sit beside
    the ones that always apply, looking equally live.
    """
    edited = {}
    rows = [r for r in schema.fields_of(cfg, section)
            if r[0] not in schema.HIDDEN]
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


def surface_folder_box(default, key='surface_folder'):
    """Where the meteorological record and its companions live.

    The same shape as :func:`gis_folder_box`: a folder stated ONCE, and the
    pickers below it are relative to it. It used to be a caption stating a
    path the modeller could not change -- which is fine until the records are
    somewhere else, and then there is nothing to do about it.
    """
    folder = st.text_input(
        'Folder with SURFACE information',
        value=st.session_state.get(key, str(default)), key=key,
        help='Where MMsurf reads the meteorological record, the irrigation '
             'series and the crop schedules. Defaults to MMsurf_ws beside '
             'the case.')
    if not os.path.isdir(folder):
        st.error('No such folder: `%s`' % folder)
    return folder


def _as_pattern(name):
    """``__inputFIELD1_crop_schedule.txt`` -> ``__inputFIELD%d_...``.

    The field number is the FIRST run of digits: a schedule is picked by
    pointing at one field's file, and what has to be stored is the pattern
    the other fields are read with.
    """
    import re

    if '%d' in name:
        return name, ''
    new, n = re.subn(r'\d+', '%d', name, count=1)
    if not n:
        return name, ('%s has no field number in it, so it cannot be a '
                      'pattern -- the name needs a %%d where the field '
                      'index goes.' % name)
    return new, ''


def path_box(key, label, value, folder, help_=None, pattern=False,
             collapsed=False):
    """A file named relative to ``folder``, or chosen in the OS dialog.

    Stores the BARE NAME while the file is in that folder, so the
    configuration stays portable between machines, and an absolute path
    otherwise -- the same rule the layer pickers follow.
    """
    pending = key + '.__pending'
    if pending in st.session_state:
        st.session_state[key] = st.session_state.pop(pending)

    c1, c2 = st.columns([6, 1])
    with c1:
        typed = st.text_input(
            label, value=str(value or ''), key=key, help=help_,
            label_visibility='collapsed' if collapsed else 'visible')
    with c2:
        if not collapsed:
            st.markdown('<div style="height:1.85rem"></div>',
                        unsafe_allow_html=True)
        if st.button('…', key=key + '.__pick',
                     help='Choose the file on this machine'):
            start = folder if os.path.isdir(folder or '') else os.getcwd()
            got, why = native_dialog('file', start, 'Select %s' % label)
            if got:
                # RELATIVE to the folder wherever it sits below it, not only
                # when it sits directly in it: the soil parameters live in
                # MF_ws/ under the dataset, and storing an absolute path for
                # them would tie the configuration to this machine.
                name = got
                if os.path.isdir(folder or ''):
                    try:
                        rel = os.path.relpath(got, os.path.abspath(folder))
                    except ValueError:           # a different drive
                        rel = ''
                    if rel and not rel.startswith('..'):
                        name = rel.replace(os.sep, '/')
                note = ''
                if pattern:
                    name, note = _as_pattern(name)
                st.session_state[pending] = name
                if note:
                    st.session_state[key + '.__why'] = note
                st.rerun()
            elif why:
                st.session_state[key + '.__why'] = why
    if st.session_state.get(key + '.__why'):
        st.caption(st.session_state.pop(key + '.__why'))
    return typed


def file_picker(dotted, value, folder, pattern=False):
    """One of panel 2's records: the box, its label, and the dialog."""
    label, units, help_ = schema.describe(dotted)
    shown = '%s [%s]' % (label, units) if units else label
    return path_box(dotted, shown, value, folder,
                    help_='`%s`  \n%s' % (dotted, help_), pattern=pattern)


def read_only(dotted, value):
    """One field, shown as it stands and not editable.

    SHOWN, not hidden and not cleared: a switch that is off is a statement
    about the RUN, not about the answers, and a modeller turning irrigation
    back on should find the series still named. It returns nothing, so a save
    leaves the field exactly as the file has it.
    """
    label, units, help_ = schema.describe(dotted)
    shown = '%s [%s]' % (label, units) if units else label
    hint = '`%s`  \n%s' % (dotted, help_)
    if schema.is_source(value):
        which = value.producer() or 'value'
        st.markdown('**%s**' % shown, help=hint)
        st.text_input(shown, value='%s: %s' % (which,
                                               getattr(value, which, '')),
                      disabled=True, key='ro_%s' % dotted,
                      label_visibility='collapsed')
    elif isinstance(value, bool):
        st.checkbox(shown, value=value, disabled=True, help=hint,
                    key='ro_%s' % dotted)
    elif isinstance(value, (int, float)):
        st.number_input(shown, value=value, disabled=True, help=hint,
                        key='ro_%s' % dotted)
    else:
        st.text_input(shown, value=str(value or ''), disabled=True,
                      help=hint, key='ro_%s' % dotted)


def rows_form(cfg, rows, section, columns=2, folder=None, files=(),
              patterns=(), gated=None):
    """A section laid out ROW BY ROW, as written, rather than in field order.

    The order a dataclass happens to declare its fields in is not the order a
    modeller fills them in, and a round-robin across columns puts unrelated
    answers side by side. Fields named in ``files`` get the system dialog.

    ``gated`` maps a switch to the fields it controls: with the switch off
    they are drawn READ-ONLY rather than live, so a panel cannot be left
    saying something the run will not do.
    """
    values = dict(schema.fields_of(cfg, section))
    edited = {}
    off = set()
    for switch, dependents in (gated or {}).items():
        # Read live, so the fields grey on the same run the box is unticked.
        # Every switch here is drawn BEFORE what it controls, which is the
        # layout's job to keep true.
        if not live(switch, values.get(switch, True)):
            off.update(dependents)
    for row in rows:
        cols = st.columns(max(columns, len(row)))
        for k, dotted in enumerate(row):
            if dotted is None or dotted not in values:
                continue          # an empty cell, so the column keeps place
            with cols[k]:
                if dotted in off:
                    read_only(dotted, values[dotted])
                    continue
                if dotted in schema.LAYER_LIST:
                    edited[dotted] = layer_box(dotted, values[dotted],
                                               cfg.layers.nlay)
                    continue
                of = schema.COLUMN_OF.get(dotted)
                if of is not None:
                    # A plain string field that names an ATTRIBUTE of the
                    # layer another field names -- so it is chosen from that
                    # layer, exactly as a source's column is.
                    edited[dotted] = column_box(
                        dotted, values[dotted], dotted,
                        live(of, values.get(of, '')))
                    continue
                if dotted in files:
                    got = file_picker(dotted, values[dotted], folder,
                                      pattern=dotted in patterns)
                else:
                    got = _widget(dotted, values[dotted])
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
    """Where to look for cartography. Shared by the pickers below.

    Returns ``(folder, shapefiles, rasters)``. An ESRI grid is a DIRECTORY of
    .adf files, so the raster list is not just a filtered file listing.
    """
    from marmites_vector import find_rasters, find_shapefiles

    gis = str(mm_paths.GIS)
    folder = st.text_input(
        'Folder with GIS information (defined in panel Home)',
        value=st.session_state.get(folder_key, gis), key=folder_key,
        help='Defaults to DATA_ROOT/GIS, which panel Home sets. The '
             'shapefiles and rasters stay here and are read only by the '
             'converter -- they never enter the repository.')
    found = find_shapefiles(folder)
    if not os.path.isdir(folder):
        st.error('No such folder: `%s`' % folder)
    elif not found:
        st.warning('No .shp in `%s` (or one level below it).' % folder)
    return folder, found, find_rasters(folder)


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
    hint = '`%s`  \n%s' % (dotted, help_)
    # The TITLE on its own line and the geometry under it, so the four
    # pickers' titles sit at the same height across the columns instead of
    # stepping down as the longer ones wrap.
    #
    # ONE block, not a markdown and a caption: two blocks put a paragraph
    # gap between the title and its bracket, and the widget's own label is
    # collapsed below -- which takes its help icon with it. Written here, the
    # icon sits on the title where it is being looked for.
    title, _sep, qualifier = label.partition(' (')
    st.markdown('**%s**  \n<span style="font-size:0.8em;opacity:0.6">%s'
                '</span>' % (title, ('(%s' % qualifier) if qualifier
                             else '&nbsp;'),
                help=hint, unsafe_allow_html=True)

    current = os.path.join(gis, value) if value and not os.path.isabs(value) \
        else (value or '')
    options = list(found)
    if current and current not in options:
        options.insert(0, current)
    if optional:
        options = [NONE] + options
    if not options:
        return st.text_input(shown, value=str(value), key=dotted,
                             label_visibility='collapsed', help=hint)

    pick = st.selectbox(
        shown, options,
        index=options.index(current) if current in options else 0,
        key=dotted, help=hint, label_visibility='collapsed',
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
    folder, found, rasters = gis_folder_box()
    cols = st.columns(4)
    with cols[0]:
        picked = layer_picker('grid.boundary', values['grid.boundary'],
                              folder, found)
    with cols[1]:
        streams = layer_picker('grid.streams', values['grid.streams'],
                               folder, found, optional=True)
    with cols[2]:
        ponds = layer_picker('grid.ponds', values['grid.ponds'],
                             folder, found, optional=True)
    with cols[3]:
        dem = layer_picker('grid.dem', values['grid.dem'], folder, rasters,
                           optional=True)
    edited['grid.boundary'] = picked
    edited['grid.streams'] = streams
    edited['grid.ponds'] = ponds
    edited['grid.dem'] = dem

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
    if dem:
        from marmites_vector import check_raster
        reports['grid.dem'] = check_raster(
            resolve_layer(dem), expect_epsg=(rep['epsg'] or cfg.grid.crs_epsg),
            against=rep['bbox'])

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
        state = live(switch, values.get(switch, True))
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
                note = schema.GRID_NOTES.get(dotted)
                if note is not None:
                    st.caption(note(
                        lambda d, f=None: live(d, values.get(d, f))))
    return out


def record_lines(cfg, dotted, folder):
    """``[(ok, name, what, size)]`` for the records this table reads.

    Stated under the TABLE rather than in a list of its own: a block of green
    and red dots says whether files are there, and nothing about what they
    are for. Under the table that reads them it is one answer.
    """
    fields = schema.SURFACE_TABLE_FILES.get(dotted, ())
    out = []
    for field in fields:
        label = schema.describe(field)[0]
        name = live(field, getattr(cfg.surface, field.split('.')[-1], ''))
        if field == 'surface.irr_ts' or field == 'surface.crop_schedule':
            if not live('surface.irrigation', cfg.surface.irrigation):
                continue
        names = [name]
        if '%d' in str(name):
            n = int(live('surface.nfield', cfg.surface.nfield) or 0)
            names = [name % (f + 1) for f in range(n)]
        for one in names:
            p = one if os.path.isabs(str(one)) else os.path.join(folder, one)
            ok = os.path.exists(p)
            out.append((ok, one, label,
                        '%.1f KB' % (os.path.getsize(p) / 1024.0) if ok
                        else 'not found'))
    return out


def table_form(cfg, dotted, singular, path, records=()):
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
    for ok, name, what, size in records:
        st.caption('%s `%s` — %s, %s' % ('🟢' if ok else '🔴', name, what.lower(),
                                         size))
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


def switch_and_save(cfg, panel):
    """The group's switch and the panel's save, side by side at the TOP.

    Returns ``(edited, slot)``. The save is DRAWN here and FILLED last: it
    writes everything the panel collected, and a panel collects it tab by
    tab as the tabs are drawn, so a button that captured the edits at this
    point would capture an empty dict. The container reserves the place and
    ``save_button(slot=...)`` puts the button in it once every sub-panel has
    had its say.

    One switch, one save, both where the modeller looks first -- rather than
    a button at the foot of a page whose tabs may each be a screenful.
    """
    edited = {}
    c1, c2 = st.columns([2, 3])
    with c1:
        edited.update(master_switch(cfg, panel[3]) or {})
    return edited, c2.container()


def save_button(cfg, path, edited, label='Validate & save', slot=None,
                key='save_panel'):
    """The one way a panel writes. Validates first, and says what changed.

    ``slot`` is a container reserved earlier -- see :func:`switch_and_save`.
    Without one the button is drawn where it stands, under a rule.

    ``key`` is CONSTANT, and has to be. It used to be ``'save_%s' % id(edited)``
    -- the id of a dict built fresh on every run -- so the button was a
    different widget each rerun: the click arrived for a key that no longer
    existed and the new button read False. It appeared to work whenever
    CPython happened to hand the new dict the address the old one had just
    freed, which is most of the time and not all of it. The symptom was a
    Validate & save that silently did nothing.
    """
    box = slot if slot is not None else st.container()
    with box:
        if slot is None:
            st.markdown('---')
        c1, c2 = st.columns([1, 3])
        if not c1.button(label, type='primary', key=key):
            return
        try:
            applied, digest = editor.save(cfg, path, edited)
        except (editor.EditError, mcfg.ConfigError) as exc:
            c2.error('NOT saved — the configuration would be invalid:\n\n%s'
                     % exc)
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
