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

st.set_page_config(page_title='2 Surface and driving forces', page_icon='🌦️', layout='wide')
case = st.session_state.get('case', 'LaMata')
cfg, path = panelui.pick_config()
ds = panelui.dataset_banner(cfg)
panel = panelui.header(2)

edited, save_slot = panelui.switch_and_save(cfg, panel)

tab_par, tab_time, tab_tables, tab_forcing = st.tabs(
    ['Records & options', 'Time discretisation', 'Parameter tables',
     'The forcing it produces'])

# ------------------------------------------------------- records & options
with tab_par:
    in_ws = panelui.surface_folder_box(os.path.join(str(ds), 'MMsurf_ws'))

    # One column per SUBJECT: the meteorology down the left, the irrigation
    # down the right, each in the order it is filled in.
    edited.update(panelui.rows_form(
        cfg, schema.SURFACE_ROWS, 'surface', columns=2, folder=in_ws,
        files=schema.SURFACE_FILES, patterns=schema.SURFACE_PATTERNS,
        gated=schema.SURFACE_GATED))


# --------------------------------------------------- time discretisation
# HERE, and not in the MODFLOW parameter file where the aggregation limit
# used to sit under the name `nper`. The days come from the RECORD -- MMsurf
# reads the hourly series and writes one row per day -- and what the run does
# with them is a question about that record, not about MODFLOW.
with tab_time:
    st.caption('The days themselves come from the meteorological record: '
               'MMsurf writes one row per day of it. What the run does with '
               'them is decided here.')
    edited.update(panelui.rows_form(cfg, schema.TIME_ROWS, 'run', columns=2))

    ndays = msurf.record_days(str(ds))
    if ndays:
        daily = panelui.live('run.daily', cfg.run.daily)
        nsp = int(panelui.live('run.nsp', cfg.run.nsp) or 0)
        if daily:
            st.success('%d day(s) in the record, so %d stress period%s.'
                       % (ndays, nsp or ndays, '' if (nsp or ndays) == 1
                          else 's'))
        else:
            most = int(panelui.live('run.perlen_max', cfg.run.perlen_max) or 1)
            st.info('%d day(s) in the record. Aggregated, the count is not '
                    'known until the rainfall is read: a day with rain gets a '
                    'period of its own and the dry days after it are averaged '
                    'up to %d, so the run reports it as it starts. At most '
                    '%d period(s), at least %d.'
                    % (ndays, most, ndays, -(-ndays // most)))
        if nsp:
            st.warning('**run.nsp = %d**: the run stops after %d stress '
                       'period(s) instead of covering the record. Useful to '
                       'try a change; misleading in a result.' % (nsp, nsp))
    else:
        st.caption('The day count comes from `inputDATE.txt`, which MMsurf '
                   'writes; it is not there yet.')


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
        panelui.table_form(cfg, dotted, singular, path,
                           records=panelui.record_lines(cfg, dotted, in_ws))

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

# Drawn at the TOP, beside the switch, and FILLED here: it writes everything
# the panel collected, and a panel collects it tab by tab as the tabs are
# drawn. Calling it earlier would capture an empty dict.
panelui.save_button(cfg, path, edited, slot=save_slot)
