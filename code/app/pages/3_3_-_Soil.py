# -*- coding: utf-8 -*-
"""Panel 3 -- SOIL: what MMsoil reads.  WP1d.

The soil column -- zones, thickness, the per-zone parameter file --
and the vegetation cover each cell carries. MODFLOW's own inputs
moved to panel 4 when this page was split: they are one model but two
different sets of questions, answered at different times.

The spatial inputs are VECTOR LAYERS, wrapped onto whichever grid
panel 1 produced. Each declares where its value comes from, with an
explicit precedence -- raster beats layer beats a single value --
rather than being decided by which file happens to exist.
"""

import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP, os.path.join(CODE, 'ppMF6')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg              # noqa: E402
import mm_paths                             # noqa: E402
from lib import panelui, schema             # noqa: E402

st.set_page_config(page_title='3 Soil', page_icon='🌱', layout='wide')
case = st.session_state.get('case', 'LaMata')
cfg, path = panelui.pick_config()
ds = panelui.dataset_banner(cfg)
panel = panelui.header(3)

edited = panelui.panel_switch(cfg, panel)

# The CARTOGRAPHY tab that stood here (the converter) moved to the Grid
# panel: it serves panels 1 to 5, and Launch now converts by itself when a
# shapefile has changed (lib.dataset_state).
tab_soil, tab_veg = st.tabs(['Soil column', 'Vegetation characteristics'])

# ------------------------------------------------------------------ soil
with tab_soil:
    edited.update(panelui.rows_form(cfg, schema.SOIL_ROWS, 'soil',
                                    columns=2, folder=str(ds),
                                    files=schema.SOIL_DATASET_FILES))

    st.info('**Precedence is explicit, not decided by which file exists:** a '
            'raster beats a polygon attribute, and a polygon attribute beats '
            'a single value. Nothing set is an error, not a silent zero.')

    # THE SOIL COLUMN, AS TABLES. It was the positional inputSOILparam.txt,
    # and the run read it from a path hard-coded in the driver -- so the
    # field that named it here changed nothing. The column is now edited in
    # place, and it is what MMsoil receives.
    st.markdown('#### The soil column')
    st.caption('One row per soil zone, in zone-CODE order: row N is the zone '
               'whose code is N in the soil-zone layer above. Then one row '
               'per horizon, top to bottom, each naming its zone. The '
               'horizons of a zone share its thickness through `slprop`, '
               'which must sum to 1.')
    for dotted in ('soil.zone', 'soil.horizon'):
        singular = next(s for d, _n, s, _c in schema.TABLES if d == dotted)
        panelui.table_form(cfg, dotted, singular, path)

    # IMPORT, ONCE, from an old parameter file -- the way a case that still
    # has one (CdL) moves over. It writes through the SAME save as the
    # sidebar, validated, so a file the rules refuse is refused here too.
    with st.expander('Import the column from an old inputSOILparam.txt'):
        _default = os.path.join('MF_ws', 'inputSOILparam.txt')
        imp = st.text_input('File, in the dataset folder or absolute',
                            value=_default, key='soil_import_path')
        full = imp if os.path.isabs(imp) else os.path.join(str(ds), imp)
        st.caption('`%s` — %s' % (full, 'found' if os.path.exists(full)
                                  else 'not there'))
        if st.button('Import and replace the tables', key='soil_import',
                     disabled=not os.path.exists(full)):
            try:
                zones, horizons = mcfg.soil_tables_from_param_file(full)
            except mcfg.ConfigError as exc:
                st.error(str(exc))
            else:
                panelui.remember_table('soil.zone', zones)
                panelui.remember_table('soil.horizon', horizons)
                _applied, why, ok = panelui.save_now(cfg, path)
                if not ok:
                    st.error('NOT imported -- the column would be invalid:'
                             '\n\n%s' % why)
                else:
                    # The data editors hold their own copy; drop it so they
                    # redraw from the file that was just written.
                    for k in ('tbl_soil.zone', 'tbl_soil.horizon'):
                        st.session_state.pop(k, None)
                    st.session_state['__saved_note'] = (
                        'Imported %d zone(s) and %d horizon(s) from %s'
                        % (len(zones), len(horizons),
                           os.path.basename(full)))
                    st.rerun()

# ------------------------------------------------------------ vegetation
# [soil] carries the vegetation COVER as well, because they share a file --
# not because they are one question. Shown beside the soil settings the
# second was read as more of the first.
with tab_veg:
    edited.update(panelui.rows_form(cfg, schema.SOIL_VEG_ROWS, 'soil',
                                    columns=2, folder=str(mm_paths.GIS),
                                    files=schema.SOIL_FILES))

    st.markdown('#### Vegetation classes')
    # Concatenated rather than %-formatted: the sentence itself contains a
    # per-cent sign, and one of those in a format string is a bug waiting for
    # the day somebody adds a second placeholder.
    st.caption('What the layer\'s class column holds, and which vegetation '
               'type of ' + schema.panel_name(2) + ' it means. The share of '
               'each cell covered is an exact AREA OVERLAY — a cell 37 % '
               'covered gets 37, not the class that happened to sit under '
               'its centre.')
    # ONLY the vegetation table. This looped over every table of the panel,
    # which was fine while there was one -- the soil column's two tables
    # belong on the Soil column tab, and would have appeared here too.
    for dotted, num, singular, _c in schema.TABLES:
        if dotted == 'soil.veg_class':
            panelui.table_form(cfg, dotted, singular, path)

    names = [v.name for v in cfg.surface.vegetation]
    if names:
        st.caption('Vegetation types available from %s: '
                   % schema.panel_name(2) +
                   ', '.join('**%d** %s' % (k + 1, n)
                             for k, n in enumerate(names)))


# REMEMBERED, not written: the one save is in the sidebar (see panelui).
panelui.remember(edited)
panelui.sidebar_save(cfg, path)
