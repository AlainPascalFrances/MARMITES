# -*- coding: utf-8 -*-
"""Configuration — edit and validate a run configuration.  WP1b 1b.5 + WP1.7.

The widgets are BUILT FROM THE SCHEMA, not hand-written, so a key added to
`RunConfig` appears here by itself. Saving refuses on a validation failure --
the same fail-fast rule as WP0.6 -- and writes to `code/configs/`, never to a
run's `resolved_config.toml`.

The second tab is WP1.7: every MMsurf parameter, ALREADY FILLED from the case
study, with its units and its meaning, instead of the positional wall of
numbers the ini actually is.
"""

import dataclasses
import os
import sys

import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg        # noqa: E402
import mmsurf_config as msc           # noqa: E402
import mm_paths                       # noqa: E402

st.set_page_config(page_title='Configuration', page_icon='⚙️', layout='wide')
case = st.session_state.get('case', 'LaMata')
CONFIG_DIR = os.path.join(CODE, 'configs')

st.title('Configuration')

files = sorted(f for f in os.listdir(CONFIG_DIR) if f.endswith('.toml')) \
    if os.path.isdir(CONFIG_DIR) else []
if not files:
    st.error('No configuration in %s' % CONFIG_DIR)
    st.stop()

chosen = st.selectbox('Configuration file', files,
                      index=files.index('lamata.toml') if 'lamata.toml' in files else 0)
path = os.path.join(CONFIG_DIR, chosen)
try:
    cfg = mcfg.load_run_config(path)
except mcfg.ConfigError as exc:
    st.error(str(exc))
    st.stop()

st.caption('`%s` — hash `%s`' % (path, cfg.config_hash()))

tab_run, tab_mmsurf, tab_raw = st.tabs(
    ['Run configuration', 'MMsurf parameters (read-only)', 'TOML'])

# --------------------------------------------------------------- run config
with tab_run:
    st.caption('Widgets are generated from the schema, so a new key appears '
               'here without touching this page.')
    edited = {}
    sections = [n for n in mcfg._SECTIONS]
    cols = st.columns(2)
    for k, name in enumerate(sections):
        sec = getattr(cfg, name)
        with cols[k % 2].expander('[%s]' % name, expanded=name in
                                  ('meta', 'run', 'layers', 'seep')):
            for f in dataclasses.fields(sec):
                v = getattr(sec, f.name)
                key = '%s.%s' % (name, f.name)
                if dataclasses.is_dataclass(v):
                    st.caption('`%s` — nested, edit in the TOML tab' % key)
                    continue
                if isinstance(v, bool):
                    edited[key] = st.checkbox(key, value=v, key=key)
                elif isinstance(v, int):
                    edited[key] = st.number_input(key, value=int(v), step=1, key=key)
                elif isinstance(v, float):
                    edited[key] = st.number_input(key, value=float(v),
                                                  format='%g', key=key)
                elif isinstance(v, (list, dict)):
                    st.caption('`%s` = `%s` — edit in the TOML tab' % (key, v))
                else:
                    edited[key] = st.text_input(key, value=str(v), key=key)

    st.markdown('---')
    c1, c2 = st.columns([1, 3])
    if c1.button('Validate & save', type='primary'):
        overrides = []
        for key, val in edited.items():
            sec, leaf = key.split('.', 1)
            cur = getattr(getattr(cfg, sec), leaf)
            if str(cur) != str(val):
                overrides.append('%s=%s' % (key, val))
        if not overrides:
            c2.info('Nothing changed.')
        else:
            try:
                cfg.apply_overrides(overrides, echo=False)
            except mcfg.ConfigError as exc:
                c2.error('NOT saved — the configuration would be invalid:\n\n%s'
                         % exc)
            else:
                cfg.write_toml(path)
                c2.success('Saved %d change(s) to %s (hash %s)'
                           % (len(overrides), chosen, cfg.config_hash()))

# ------------------------------------------------------------ MMsurf (1.7)
with tab_mmsurf:
    ini = mm_paths.dataset_dir(case) / 'MMsurf_ws' / '__inputMMsurf.ini'
    st.caption('`%s`' % ini)
    if not ini.exists():
        st.warning('MMsurf ini not found for this case.')
    else:
        try:
            ms = msc.load_mmsurf_ini(str(ini))
        except msc.MMsurfError as exc:
            st.error(str(exc))
        else:
            st.write(ms.counts)
            params = msc.enumerate_parameters(ms)
            blocks = sorted({p.block for p in params},
                            key=lambda b: ['meteo', 'veg', 'crop', 'soil',
                                           'openwater'].index(b))
            for block in blocks:
                rows = [p for p in params if p.block == block]
                members = sorted({p.member for p in rows})
                with st.expander('%s — %s  (%d parameter(s))'
                                 % (block, ', '.join(members), len(rows)),
                                 expanded=(block == 'veg')):
                    st.dataframe(
                        [{'member': p.member, 'parameter': p.name,
                          'value': p.value, 'units': p.units,
                          'description': p.description} for p in rows],
                        width='stretch', hide_index=True)
            st.info('Read-only for now. These migrate into the schema as '
                    '`[mmsurf.*]`; `Zr` is the natural source for the WP2 '
                    '`[et] extdp_source = "veg_zone"` extinction depth.')

# ------------------------------------------------------------------- raw
with tab_raw:
    st.code(open(path, encoding='utf-8').read(), language='toml')
