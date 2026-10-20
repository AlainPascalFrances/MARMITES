# -*- coding: utf-8 -*-
"""Panel 8 -- PEST++ SET-UP: everything that defines a calibration.  WP7
front-end, phase 2 (docs/MM-MF6/WP7_PEST_frontend_plan.md).

WP7 front-end, phase 3: the BASE is a finished run of the registry -- its
configuration copy is the model calibrated, its end state where every
realisation starts (marmites_runinfo) -- or, for a state no recorded run
vouches for, this configuration and a saved state chosen by hand.

Every choice here is a field of the configuration's ``[pest]`` section,
remembered like the other panels' and written by the sidebar's one save.
The PEST tools read them (code/tools/pest_lamata_setup.py at set-up,
pest_lamata_run.py at launch), so a set-up made from this page is the
calibration the page shows. Four tabs, in the order a calibration is
thought through:

  Base and workspace         where every realisation starts, what is fitted
  Parameters                 the table: one row per family and layer
  Observations and weights   the groups, their share of phi, their noise
  Ensemble and PEST++        realisations, iterations, agents, the start

The PEST workspace is set on Home and only shown here (user, 2026-10-17).
Tab ⑤ (phase 4) sets the calibration up and, once a set-up has
succeeded, launches it -- detached jobs, lib/pest_jobs.py; the progress is
panel 9's.
"""

import copy
import glob
import os
import sys

import pandas as pd
import streamlit as st

APP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.abspath(os.path.join(APP, '..'))
for p in (CODE, APP, os.path.join(CODE, 'tools')):
    if p not in sys.path:
        sys.path.insert(0, p)

import marmites_config as mcfg          # noqa: E402
import marmites_runinfo as ri           # noqa: E402
import mm_paths                         # noqa: E402
from lib import editor, panelui, schema  # noqa: E402
from lib import runs as runlib          # noqa: E402

st.set_page_config(page_title='8 PEST++ set-up', page_icon='🎯',
                   layout='wide')
cfg, path = panelui.pick_config()
ds = str(panelui.dataset_banner(cfg))

st.title('🎯  8 — PEST++ set-up')
st.caption('What a calibration of this configuration is: its start, its '
           'parameters, its observations, its ensemble. Saved with the '
           'configuration (`[pest]`) by *Validate & save* in the sidebar; '
           'the PEST tools read it at set-up and at launch.')
WS = str(mm_paths.PEST_WS) if mm_paths.PEST_WS else ''
st.text_input('PEST workspace (set on Home)', value=WS or 'not set',
              disabled=True, key='pest_ws_shown')

try:
    import pest_lamata_setup as pls
    import pest_lamata_run as plr
except Exception as exc:                                 # noqa: BLE001
    st.error('The PEST tools cannot be loaded: %s' % exc)
    st.stop()

edited = {}

# what each observation group is (pest_lamata_obs)
GROUP_WHAT = {
    'hds_p0': "P0's head, near-daily, measured",
    'hds_c': 'C1..C5, W1 heads, quarterly: soft data',
    'sm_p0': "P0's soil moisture, 2 horizons",
    'sm_sm': "SM's soil moisture, 2 horizons",
    'q_log': 'outlet flow, log10(Q + offset), daily',
    'q_dry': "the gauge's dry days: the model's flow against zero",
    'q_year': 'outlet volume per hydrological year',
    'et_rs': 'MODIS catchment ET per 8-day composite',
    'et_pt': 'measured ET at the points (EC tower...): validation',
}


def _what(cfg, r):
    """A parameter row in words, with its units."""
    kind, units, how = mcfg.PEST_FAMILIES.get(r['family'], ('', '', ''))
    if kind == 'layer':
        txt = 'layer %d' % r['layer']
    elif kind == 'horizon' and 1 <= r['layer'] <= len(cfg.soil.horizon):
        h = cfg.soil.horizon[r['layer'] - 1]
        txt = 'horizon %d (zone %d)' % (r['layer'], h.zone)
        if r['family'] == 'sfc':
            txt += ', within sr %g .. smax %g' % (h.sr, h.smax)
    else:
        txt = ''
    return '%s%s — %s' % (txt + ', ' if txt else '', how, units)


tab_base, tab_par, tab_obs, tab_ens, tab_run = st.tabs([
    '① Base and workspace', '② Parameters', '③ Observations and weights',
    '④ Ensemble and PEST++', '⑤ Set up and launch'])

# ------------------------------------------------------------------- base
@st.cache_data(ttl=60, show_spinner=False)
def _base_of(run_dir, ws_root, stamp):
    """marmites_runinfo.base_of, read again when the run's files change
    (``stamp``)."""
    return ri.base_of(run_dir, ws_root)


def _stamp(run_dir):
    return tuple(os.path.getmtime(os.path.join(run_dir, f))
                 if os.path.exists(os.path.join(run_dir, f)) else 0.0
                 for f in ('run.log', ri.INFO, 'status.json'))


def _differs(a, b, skip=('pest', 'ui', 'meta')):
    """[(dotted key, this configuration's value, the run's)] where they
    differ -- the calibration's own settings left out."""
    out = []

    def walk(x, y, pre):
        if isinstance(x, dict) and isinstance(y, dict):
            for k in sorted(set(x) | set(y)):
                walk(x.get(k), y.get(k), '%s.%s' % (pre, k))
        elif x != y:
            out.append((pre, x, y))
    da, db = a.to_dict(), b.to_dict()
    for s in da:
        if s not in skip:
            walk(da[s], db.get(s), s)
    return out


with tab_base:
    st.markdown('#### Where every realisation starts')
    RUNS = ri.runs_dir(cfg, mm_paths.WS_ROOT)
    done = [r['run_id'] for r in runlib.list_runs(RUNS, limit=60)
            if r.get('state') == 'finished'
            and r.get('outcome') == 'completed']
    cur_run = (cfg.pest.base_run or '').strip()
    ids = [''] + done + ([cur_run] if cur_run and cur_run not in done
                         else [])
    label, _u, help_ = schema.describe('pest.base_run')
    kw = {} if 'pest.base_run' in st.session_state else {
        'index': ids.index(cur_run) if cur_run in ids else 0}
    base_run = st.selectbox(
        label, ids, key='pest.base_run',
        format_func=lambda r: r or ('— none: this configuration and a saved '
                                    'state chosen by hand —'),
        help='`pest.base_run`  \n%s  \nThe completed runs of `%s`.'
             % (help_, RUNS), **kw) or ''
    edited['pest.base_run'] = base_run
    base_cfg = cfg                  # the model the calibration would run
    if base_run:
        rd = os.path.join(RUNS, base_run)
        b = _base_of(rd, str(mm_paths.WS_ROOT), _stamp(rd))
        if b['why']:
            st.error('Run `%s` cannot be the base: %s.' % (base_run, b['why']))
        else:
            import time
            st.success(
                '**The model:** the configuration run `%s` ran (`%s`%s).  \n'
                '**The start:** its end state `%s`, saved %s in `%s`%s.'
                % (base_run, b['config'],
                   ', recovered from its results folder'
                   if b['source'] == 'log' else '',
                   b['state'], time.strftime('%d %b %Y %H:%M',
                                             time.localtime(b['state_time'])),
                   b['state_dir'],
                   '' if b.get('converged', True) else
                   ' -- a spin-up saved BEFORE it converged'))
            try:
                base_cfg = mcfg.load_run_config(b['config'])
                base_cfg.pest = copy.deepcopy(cfg.pest)
            except mcfg.ConfigError as exc:
                st.error('Its configuration copy does not load: %s' % exc)
                base_cfg = cfg
            diff = _differs(cfg, base_cfg) if base_cfg is not cfg else []
            if diff:
                st.warning('Its configuration differs from this one in %d '
                           'setting(s): the calibration runs the RUN\'s, '
                           'with this configuration\'s calibration settings.'
                           % len(diff))
                with st.expander('The differences'):
                    st.dataframe(pd.DataFrame(
                        [(k, str(x), str(y)) for k, x, y in diff],
                        columns=['setting', 'this configuration', 'the run']),
                        hide_index=True, width='stretch')
            st.caption('The run covers %s.'
                       % ('the whole record' if not base_cfg.run.nsp else
                          '%d stress periods' % base_cfg.run.nsp))
    else:
        state_dir = mcfg.state_workspace(cfg, mm_paths.WS_ROOT)
        found = {d['name']: d for d in mcfg.saved_states(
            state_dir, int(cfg.layers.nlay), 'spinup.strt_heads')}
        cur = (cfg.pest.base_state or '').strip()
        names = [''] + list(found) + ([cur] if cur and cur not in found
                                      else [])

        def _fmt(name):
            if not name:
                return '— choose a saved state —'
            d = found.get(name)
            if d is None:
                return '%s  ·  not found in the workspace' % name
            import time
            bits = [name, time.strftime('%d %b %Y %H:%M',
                                        time.localtime(d['mtime']))]
            if d['shape']:
                bits.append('%s cells' % mcfg.cells_text(d['shape']))
            if not d['full']:
                bits.append('heads only')
            return '  ·  '.join(bits)

        label, _u, help_ = schema.describe('pest.base_state')
        kw = {} if 'pest.base_state' in st.session_state else {
            'index': names.index(cur) if cur in names else 0}
        edited['pest.base_state'] = st.selectbox(
            label, names, format_func=_fmt, key='pest.base_state',
            help='`pest.base_state`  \n%s  \nFrom `%s`.' % (help_, state_dir),
            **kw) or ''
        if not edited['pest.base_state']:
            st.warning('No saved state chosen: the set-up needs one.')
        else:
            st.warning('The model calibrated is THIS configuration as it is '
                       'now, and nothing vouches that it produced the state '
                       'chosen. A run solution, above, carries both -- '
                       'every run the Run panel launches records them.')
    ds = str(mm_paths.dataset_dir(base_cfg.paths.case))

    st.markdown('#### What is fitted')
    edited.update(panelui.rows_form(cfg, schema.PEST_BASE_ROWS, 'pest'))
    probe = copy.deepcopy(base_cfg)
    probe.pest.calib_start = str(edited.get('pest.calib_start',
                                            cfg.pest.calib_start) or '')
    try:
        st.caption('The first day fitted: **%s**; the observations before '
                   'it are reported, never fitted.'
                   % pls.calib_start(probe, ds))
    except (ValueError, OSError, SystemExit) as exc:
        st.error('No first day can be worked out: %s' % exc)

# ------------------------------------------------------------- parameters
with tab_par:
    st.markdown('#### Pilot points (K, Sy, Ss)')
    edited.update(panelui.rows_form(cfg, schema.PEST_PILOT_ROWS, 'pest',
                                    columns=3))
    st.markdown('#### The parameter table')
    try:
        defaults = pls.default_pars(base_cfg, ds)
    except (SystemExit, OSError, ValueError) as exc:
        defaults = []
        st.error('No default table for this configuration: %s' % exc)
    saved = editor.table_rows(cfg, 'pest.par')
    # what the file says (or, empty, the defaults) -- unless an edit not
    # yet saved is remembered: it outlives a visit to another page
    shown = saved or defaults
    rows = panelui.pending_tables().get('pest.par', shown)
    if not cfg.pest.par:
        st.info('The configuration holds no table yet: these are the '
                'DEFAULTS for it (the calibration of 2026-10-18 -- the '
                'user\'s hydraulic table for K, Sy and Sr; each horizon\'s '
                'field capacity within its wilting point and porosity). '
                'Edit them and save to write the table.')
    if st.button('Reset the table to the defaults', key='pest_par_reset',
                 disabled=not cfg.pest.par or not defaults):
        panelui.remember_table('pest.par', defaults)
        st.session_state.pop('tbl_pest.par', None)
        st.rerun()
    df = pd.DataFrame(rows, columns=['family', 'layer', 'enable', 'start',
                                     'lower', 'upper', 'log'])
    df.insert(2, 'what', [_what(base_cfg, r) for r in rows])
    cc = {'family': st.column_config.TextColumn('Family', disabled=True),
          'layer': st.column_config.NumberColumn('Layer / horizon',
                                                 disabled=True),
          'what': st.column_config.TextColumn('What, units', disabled=True,
                                              width='large'),
          'enable': st.column_config.CheckboxColumn('On'),
          'start': st.column_config.NumberColumn('Start', format='%.6g'),
          'lower': st.column_config.NumberColumn('Lower', format='%.6g'),
          'upper': st.column_config.NumberColumn('Upper', format='%.6g'),
          'log': st.column_config.CheckboxColumn('Log')}
    with st.expander('What each column means'):
        for name in ('family', 'layer', 'enable', 'start', 'lower', 'upper',
                     'log'):
            lab, _u, hlp = schema.describe('par.%s' % name)
            st.markdown('- **%s** (`%s`) — %s' % (lab, name, hlp))
    got = st.data_editor(df, column_config=cc, num_rows='fixed',
                         hide_index=True, width='stretch',
                         key='tbl_pest.par')
    out = [{'family': str(r['family']), 'layer': int(r['layer']),
            'enable': bool(r['enable']), 'start': float(r['start']),
            'lower': float(r['lower']), 'upper': float(r['upper']),
            'log': bool(r['log'])} for r in got.to_dict('records')]
    # remembered only when it says something the file does not: the
    # defaults shown for an empty table are not an edit
    if out != shown:
        panelui.remember_table('pest.par', out)
        st.caption('Changed — *Validate & save* in the sidebar writes it '
                   'with the rest of the configuration.')
    else:
        panelui.forget_table('pest.par')
    probe = copy.deepcopy(base_cfg)
    probe.pest.par = [mcfg.PestPar(**r) for r in out]
    for k in ('pest.pp_space', 'pest.vario_range', 'pest.k_order_decades'):
        if k in edited:
            setattr(probe.pest, k.split('.')[1], float(edited[k]))
    bad = [e for e in probe.problems() if e.startswith('pest.par')]
    for e in bad:
        st.error(e)
    n_on = sum(1 for r in out if r['enable'])
    n_pp = sum(1 for r in out if r['enable']
               and r['family'] in mcfg.PEST_PILOT)
    st.caption('%d row(s) on: %d on pilot points (a point every %g m on '
               'each layer\'s active cells), %d one-value parameter(s). The '
               'prior\'s spread on the mesh is drawn at set-up.'
               % (n_on, n_pp,
                  float(edited.get('pest.pp_space', cfg.pest.pp_space)),
                  n_on - n_pp))
    off = [r['family'] for r in out if r['enable']
           and pls.PACKAGE_OF.get(r['family'])
           and not getattr(base_cfg, pls.PACKAGE_OF[r['family']]).enable]
    if off:
        st.warning('On, but their package is off, so not calibrated: %s'
                   % ', '.join(off))

# ----------------------------------------------------------- observations
with tab_obs:
    st.markdown('#### The groups and their share of phi')
    gsaved = editor.table_rows(cfg, 'pest.group')
    gdf = pd.DataFrame(panelui.pending_tables().get('pest.group', gsaved),
                       columns=['name', 'share', 'noise',
                                        'noise_rel'])
    gdf.insert(1, 'what', [GROUP_WHAT.get(n, '') for n in gdf['name']])
    gcc = {'name': st.column_config.TextColumn('Group', disabled=True),
           'what': st.column_config.TextColumn('What', disabled=True,
                                               width='large'),
           'share': st.column_config.NumberColumn('Share of phi',
                                                  format='%.4g'),
           'noise': st.column_config.NumberColumn('Noise', format='%.4g'),
           'noise_rel': st.column_config.NumberColumn('Relative noise',
                                                      format='%.4g')}
    with st.expander('What each column means'):
        for name in ('name', 'share', 'noise', 'noise_rel'):
            lab, _u, hlp = schema.describe('group.%s' % name)
            st.markdown('- **%s** (`%s`) — %s' % (lab, name, hlp))
    gget = st.data_editor(gdf, column_config=gcc, num_rows='fixed',
                          hide_index=True, width='stretch',
                          key='tbl_pest.group')
    gout = [{'name': str(r['name']), 'share': float(r['share']),
             'noise': float(r['noise']), 'noise_rel': float(r['noise_rel'])}
            for r in gget.to_dict('records')]
    if gout != gsaved:
        panelui.remember_table('pest.group', gout)
        st.caption('Changed — *Validate & save* in the sidebar writes it '
                   'with the rest of the configuration.')
    else:
        panelui.forget_table('pest.group')
    total = sum(g['share'] for g in gout if g['name'] not in mcfg.PEST_SOFT)
    (st.caption if abs(total - 1.0) <= 1e-6 else st.error)(
        'The shares of the measured groups sum to %g%s.'
        % (total, '' if abs(total - 1.0) <= 1e-6 else ' -- they must make 1'))
    edited.update(panelui.rows_form(cfg, schema.PEST_OBS_ROWS, 'pest'))

    st.markdown('#### Evapotranspiration')
    modis = os.path.join(str(mm_paths.RS), 'MODIS4%s' % cfg.paths.case)
    n_hdf = len(glob.glob(os.path.join(modis, '**', '*.hdf'),
                          recursive=True)) if os.path.isdir(modis) else 0
    st.caption('%s **MODIS** (group et_rs), from the RS folder set on Home: '
               '`%s` — %d composite file(s).'
               % ('🟢' if n_hdf else '🔴', modis, n_hdf))
    pts = sorted(os.path.basename(f) for f in glob.glob(os.path.join(
        ds, '%s_*.txt' % cfg.obs.aet_prefix))) if cfg.obs.aet_prefix else []
    st.caption('**Measured ET at the points** (group et_pt: validation, '
               'never fitted), the files `%s_<point>.txt` of the dataset: %s.'
               % (cfg.obs.aet_prefix or '?', ', '.join(pts) or 'none'))

    st.markdown('#### The targets of the last set-up')
    tg_f = os.path.join(WS, 'org', 'pest_targets.csv') if WS else ''
    if tg_f and os.path.exists(tg_f):
        tg = pd.read_csv(tg_f)
        summ = tg.groupby('obgnme').agg(
            targets=('obsnme', 'size'),
            fitted=('weight', lambda w: int((w > 0).sum())),
            weight=('weight', 'max')).reset_index().rename(
            columns={'obgnme': 'group'})
        st.dataframe(summ, hide_index=True, width='stretch')
        st.caption('From `%s`: what the set-up built from its base run, '
                   'with the settings in force then.' % tg_f)
    else:
        st.caption('None yet in this PEST workspace.')

# --------------------------------------------------------------- ensemble
with tab_ens:
    st.markdown('#### Ensemble and PEST++ options')
    edited.update(panelui.rows_form(cfg, schema.PEST_ENSEMBLE_ROWS, 'pest'))
    agents = int(edited.get('pest.agents', cfg.pest.agents))
    if agents > plr.MAX_WORKERS:
        st.error('At most %d agents on this server; %d asked.'
                 % (plr.MAX_WORKERS, agents))
    start = edited.get('pest.start', cfg.pest.start)
    if start == 'warm' and WS:
        jcb = sorted(glob.glob(os.path.join(WS, '*', 'lamata.*.par.jcb')))
        if jcb:
            st.caption('Parameter ensembles in the PEST workspace:  \n'
                       + '  \n'.join('`%s`' % f for f in jcb[-8:]))
    st.caption('Set up and launch on the next tab.')

# ------------------------------------------------------- set up and launch
# Phase 4 (user, 2026-10-20): the set-up and, once one has succeeded, the
# launch -- DETACHED jobs (lib/pest_jobs.py) on the SAVED configuration,
# one at a time; the job's state and log below, refreshed while it runs.
JOB_REFRESH_S = 10


def _job_view(ws):
    """The latest PEST job: what, how it stands, its log; a stop."""
    job = pest_jobs.latest(ws)
    if job is None:
        st.caption('No job started from this page yet.')
        return
    label = {'setup': 'Set-up', 'launch': 'Calibration'}.get(job['kind'],
                                                           job['kind'])
    state = job.get('state')
    head = '**%s** `%s` -- started %s' % (label, job['job_id'],
                                         job.get('started'))
    if state == 'running':
        st.info(head + ', **running**')
    elif job.get('outcome') == 'ok':
        st.success(head + ', finished %s (%s min): **succeeded**'
                   % (job.get('finished'), job.get('minutes')))
    else:
        st.error(head + ', ended %s: **%s**%s'
                 % (job.get('finished') or '', job.get('outcome'),
                    '' if job.get('returncode') is None
                    else ' (exit %s)' % job.get('returncode')))
    st.code(pest_jobs.log_tail(ws, job['job_id'], 40) or '(no output yet)',
            language='text')
    if job['kind'] == 'launch' and state == 'running':
        panelui.panel_link('pages/9_9_-_PEST_progress.py',
                           'Follow it on PEST++ progress', icon='📈')
    if state == 'running':
        sure = st.checkbox('Stop this job, and everything it started',
                           key='pest_stop_ok')
        if st.button('Stop', key='pest_stop_go', disabled=not sure):
            pest_jobs.stop(ws, job['job_id'])
            st.rerun()


with tab_run:
    st.markdown('#### Set up, then launch')
    st.caption('Both run DETACHED in the PEST workspace -- closing the '
               'browser stops nothing -- one PEST job at a time, on the '
               'SAVED configuration.')
    if not WS:
        st.info('No PEST workspace is set yet: choose one on **Home**.')
    else:
        from lib import pest_jobs, pest_progress
        panelui.remember(edited)
        todo = panelui.unsaved_changes(cfg)
        busy = pest_jobs.active(WS)
        pp_alive = pest_progress.pestpp_running()
        py = (str(mm_paths.PYTHON_EXE)
              if os.path.exists(str(mm_paths.PYTHON_EXE)) else sys.executable)
        try:
            _m, _s, _d, base_what = pls.resolve_base(cfg, mm_paths.WS_ROOT)
            base_ok = True
        except SystemExit as exc:
            base_ok, base_what = False, str(exc)
        pest_errs = [e for e in cfg.problems() if e.startswith('pest')]
        common = []
        if todo:
            common.append('%d change(s) not saved: *Validate & save* first'
                          % len(todo))
        if busy:
            common.append('job %s is running' % busy[0]['job_id'])
        elif pp_alive:
            common.append('pestpp-ies is running on this machine')
        c1, c2 = st.columns(2)
        with c1:
            st.markdown('##### ① Set up')
            st.caption('From the base -- %s: the model copy, its base run '
                       '(about an hour), the targets and their weights, the '
                       'template and a prior of %d realisations. Replaces '
                       '`org` and `template` in the workspace.'
                       % (base_what if base_ok else 'not resolved',
                          cfg.pest.nreal))
            no_setup = list(common)
            if not base_ok:
                no_setup.append(base_what)
            no_setup += pest_errs
            if st.button('Set up the calibration', type='primary',
                         key='pest_setup_go', disabled=bool(no_setup)):
                try:
                    pest_jobs.launch('setup', path, WS, python_exe=py)
                except (pest_jobs.JobBusy, FileNotFoundError) as exc:
                    st.error('Not started: %s' % exc)
                else:
                    st.rerun()
            for why in no_setup:
                st.caption('⛔ %s' % why)
        with c2:
            st.markdown('##### ② Launch')
            ready, why_ready = pest_jobs.setup_ready(WS)
            st.caption('%d agent(s), %d iteration(s) (noptmax), from the %s: '
                       'hours to days on this server. A previous `master` '
                       'is moved aside, never deleted.'
                       % (cfg.pest.agents, cfg.pest.noptmax, cfg.pest.start))
            no_launch = list(common)
            if not ready:
                no_launch.append(why_ready)
            if cfg.pest.agents > plr.MAX_WORKERS:
                no_launch.append('at most %d agents' % plr.MAX_WORKERS)
            sure = st.checkbox('I understand: %d agents for hours to days'
                               % cfg.pest.agents, key='pest_launch_ok',
                               disabled=bool(no_launch))
            if st.button('Launch the calibration', type='primary',
                         key='pest_launch_go',
                         disabled=bool(no_launch) or not sure):
                try:
                    pest_jobs.launch('launch', path, WS, python_exe=py)
                except (pest_jobs.JobBusy, FileNotFoundError) as exc:
                    st.error('Not started: %s' % exc)
                else:
                    st.rerun()
            for why in no_launch:
                st.caption('⛔ %s' % why)
            if ready and not no_launch:
                st.caption('✅ %s' % why_ready)
        st.markdown('##### The current job')
        _live = bool(busy)
        st.fragment(run_every=JOB_REFRESH_S if _live else None)(
            _job_view)(WS)

# REMEMBERED, not written: the one save is in the sidebar (see panelui).
panelui.remember(edited)
panelui.sidebar_save(cfg, path)
