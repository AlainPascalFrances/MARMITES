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
