# -*- coding: utf-8 -*-
"""Panel 2 -- SURFACE (MMsurf).  WP1d.

MMsurf turns the hourly meteorological record into the daily forcing MMsoil
consumes: rainfall, throughfall, potential transpiration per vegetation type,
potential evaporation per surface soil, open-water evaporation and LAI.

This panel replaces ``MMsurf_ws/__inputMMsurf.ini`` entirely, and with it
``__inputMMsurf4MMsoil.txt`` -- the file MMsurf wrote and the driver read back,
which was AUTHORITATIVE: with MMsurf not running, editing Zr or kT* in the ini
changed nothing, and the two disagreed. There is one source now, and it is
what these widgets edit.

The parameter tables are edited as GRIDS, one row per zone or type, because
the ini's positional lines are the single easiest thing in this model to get
silently wrong.
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_surface as msurf            # noqa: E402
import mm_paths                             # noqa: E402
from lib import panelui, schema             # noqa: E402

st.set_page_config(page_title='2 Surface', page_icon='🌦️', layout='wide')
case = st.session_state.get('case', 'LaMata')
cfg, path = panelui.pick_config()
ds = panelui.dataset_banner(cfg)
panel = panelui.header(2)

edited = {}
edited.update(panelui.master_switch(cfg, panel[3]) or {})

tab_par, tab_tables, tab_forcing = st.tabs(
    ['Records & options', 'Parameter tables', 'The forcing it produces'])

# ------------------------------------------------------- records & options
with tab_par:
    st.markdown('#### The meteorological record and its companions')
    in_ws = panelui.surface_folder_box(os.path.join(str(ds), 'MMsurf_ws'))

    # Laid out row by row rather than in field order: the three RECORDS
    # belong together down the left, and the two spatial layers -- a
    # different kind of answer entirely -- belong together on the right.
    edited.update(panelui.rows_form(
        cfg, schema.SURFACE_ROWS, 'surface', columns=2, folder=in_ws,
        files=schema.SURFACE_FILES, patterns=schema.SURFACE_PATTERNS))

    # From the BOXES: a green light against the file that was named before
    # the last edit is worse than no light at all.
    def live(dotted, fallback):
        return panelui.live(dotted, fallback)

    rows = [(live('surface.meteo_ts', cfg.surface.meteo_ts),
             'meteorological record')]
    if live('surface.irrigation', cfg.surface.irrigation):
        rows.append((live('surface.irr_ts', cfg.surface.irr_ts),
                     'irrigation series'))
        pattern = live('surface.crop_schedule', cfg.surface.crop_schedule)
        for f in range(int(live('surface.nfield', cfg.surface.nfield) or 0)):
            try:
                name = pattern % (f + 1)
            except TypeError:
                name = pattern
            rows.append((name, 'crop schedule, field %d' % (f + 1)))
    st.markdown('#### Are they there?')
    for fn, what in rows:
        p = os.path.join(in_ws, fn)
        ok = os.path.exists(p)
        size = ('%.1f KB' % (os.path.getsize(p) / 1024.0)) if ok else 'missing'
        st.markdown('%s `%s` — %s, %s' % ('🟢' if ok else '🔴', fn, what, size))

    st.info('The meteorological record is ONE file: `Date`, `Time`, then SIX '
            'columns per station in the order **P, Ta, RHa, Pa, wind, '
            'radiation**, hourly. The header row is compulsory and is never '
            'parsed — column ORDER is the contract. The irrigation series is '
            'one column per field, and its dates are never read either, so it '
            'must be row-aligned with the meteorological record.')

    panelui.save_button(cfg, path, edited)

# ------------------------------------------------------- parameter tables
with tab_tables:
    st.caption('One row per zone or type. MMsurf prepends a reference '
               '`grassFAO56` with fixed FAO-56 parameters, so three '
               'vegetation types here means four internally — index 0 is '
               'never one of yours.')
    for dotted, num, singular, count in schema.TABLES:
        if num != 2:
            continue
        st.markdown('---')
        panelui.table_form(cfg, dotted, singular, path)

# --------------------------------------------------------- the forcing out
with tab_forcing:
    spec = msurf.forcing_spec(cfg, str(ds))
    on = cfg.run.surface
    st.markdown('#### %s' % ('MMsurf will WRITE these into the workspace'
                             if on else
                             'MMsurf is off, so these must already exist'))
    if on:
        st.caption('They are run OUTPUT and go to `%s` — never into the '
                   'repository.'
                   % msurf.surface_ws(cfg, mm_paths.WS_ROOT, cfg.paths.case))
    try:
        ndays = msurf.check_forcing(spec, must_exist=not on)
        if ndays:
            st.success('%d day(s) of forcing, and every file has the shape '
                       'its block count implies.' % ndays)
    except msurf.MMsurfError as exc:
        (st.warning if on else st.error)(str(exc))
        ndays = None

    st.markdown('The block count is what makes a truncated file findable: '
                'nothing inside these files says how many blocks they hold.')
    rows = []
    for role in spec.roles():
        p = spec.path(role)
        ok = os.path.exists(p)
        rows.append({'file': spec.files[role],
                     'blocks': spec.blocks(role),
                     'expected values': (ndays * spec.blocks(role)
                                         if ndays else '?'),
                     'present': '🟢' if ok else '🔴',
                     'size KB': round(os.path.getsize(p) / 1024.0, 1) if ok
                                else 0.0})
    st.dataframe(rows, width='stretch', hide_index=True)
    st.caption('The committed `inputZONRFe_veg_d.txt` held 4869 values against '
               'a 1949-day record — not a whole number of blocks — and that '
               'went unnoticed for years. It is checked now, before the run.')
