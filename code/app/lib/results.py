# -*- coding: utf-8 -*-
"""The Results panel's order: which tab each figure belongs on, its title,
and its place in the tab (user, 2026-10-09).

Six tabs -- input maps, output maps, time series, calibration (the state
variables against what was measured), the water budgets (the catchment
and every observation point), the ponds'
water budget -- plus "Other" for a figure no rule knows, so nothing a run
wrote is ever hidden. Streamlit-free, so it is tested without a page.
"""

import os
import re

__all__ = ['TABS', 'classify', 'title', 'arrange', 'panels', 'rows', 'FOLDERS']

TABS = (
    ('input', 'Input maps'),
    ('output', 'Output maps'),
    ('series', 'Time series'),
    ('calib', 'Calibration (state variables)'),
    ('budget', 'Water budgets'),
    ('ponds', 'Ponds water budget'),
    ('other', 'Other'),
)
# the folders of a run where the figures are, in the order they are read
FOLDERS = ('_input', '_output', 'figures_nwt_comparison')

# --- what each figure is ----------------------------------------------- #
# input fields (marmites_postprocess._native_input_maps, IN_<nnn>_<stem>)
_IN = {
    'general_map': 'General map of the site (GIS layers)',
    'model_map': 'The model: grid, streams, ponds and boundary cells',
    'elev': 'Land surface elevation',
    'top': 'Aquifer top (the soil base)',
    'botm': 'Bottom of each layer',
    'thick': 'Aquifer thickness',
    'strt': 'Initial heads',
    'gridSOILthick': 'Soil thickness',
    'hk': 'Horizontal hydraulic conductivity',
    'Ss': 'Specific storage',
    'Sy': 'Specific yield',
    'T': 'Transmissivity',
    'vka': 'Vertical hydraulic conductivity',
    'drn_cond': 'Drain conductance',
    'drn_elev': 'Drain elevation',
    'ghb_cond': 'General-head boundary conductance',
    'ghb_head': 'General-head boundary head',
    'eps': 'Unsaturated zone: Brooks-Corey exponent (epsilon)',
    'thts': 'Unsaturated zone: saturated water content',
    'thti': 'Unsaturated zone: initial water content',
    'thtr': 'Unsaturated zone: residual water content',
    'vks': 'Unsaturated zone: vertical hydraulic conductivity',
    'ibound': 'Active cells',
    'SOILzones': 'Soil zones',
    'METEOzones': 'Meteorological zones',
    'IRRzones': 'Irrigation zones',
}
# soil-column maps (marmites_postprocess._MAP_FLUXES, MMmap_<stem>)
_MM = {
    'P': 'rainfall', 'Pe': 'effective rainfall', 'I': 'infiltration',
    'Rp': 'percolation', 'theta': 'soil moisture',
    'dgwt': 'depth to the water table', 'Ro': 'runoff',
    'Ei': 'interception', 'Eow': 'open-water evaporation',
    'ETsoil': 'soil evapotranspiration', 'Eg': 'groundwater evaporation',
    'Tg': 'groundwater transpiration', 'ETg': 'groundwater evapotranspiration',
    'EXFg': 'exfiltration', 'uzthick': 'unsaturated zone thickness',
    'dSsoil': 'change in soil storage', 'dSsurf': 'change in surface storage',
    'Runon': 'run-on from upslope (cascade)',
    'Reinf': 'reinfiltrated run-on (cascade)',
    'Ecrr': 'runoff evaporated by the cascade',
}
# aquifer maps (marmites_postprocess._AQ_MAPS and the derived ones)
_GW = {
    'head': 'mean head', 'Rg': 'recharge to groundwater',
    'EXFg': 'seepage to the surface', 'DRN': 'boundary drainage',
    'ETg': 'groundwater evapotranspiration',
    'Eg': 'groundwater evaporation (EVT)',
    'Tg': 'groundwater transpiration (EVT)',
    'FLF': 'flow across the lower face of the layer (+ down, - up)',
    'FUF': 'flow across the upper face of the layer (+ down, - up)',
    'SFRg': 'exchange with the streams (+ into the aquifer)',
    'LAKg': 'exchange with the ponds (+ into the aquifer)',
    'Re': 'effective recharge (Rg - Exf)',
    'Rn': 'net recharge (Rg - Exf - ETg)',
    'dSg': 'release from groundwater storage',
}
# stream maps (marmites_postprocess._native_result_maps, SWmap_<stem>)
_SW = {
    'Q': 'mean streamflow (log10 of 1 + Q, Q in m3/d)',
}
# the comparison with the legacy NWT run (tests/plot_water_budget.py)
_NWT = {
    '01_wb_timeseries': 'Water budget over time, against the NWT run',
    '02_wb_cumulative': 'Cumulative water budget, against the NWT run',
    '03_wb_totals': 'Water budget totals, against the NWT run',
    '04_wb_hydroyear': 'Water budget per hydrological year, against the NWT '
                       'run',
    '06_heads': 'Heads, against the NWT run',
    '07_coupling': 'The coupling: what MARMITES and MODFLOW exchanged',
}
_FIXED = {
    'budget_compartment': 'Groundwater budget by compartment',
    'budget_uzf': 'Unsaturated zone (UZF) budget',
    'budget_sfr': 'Stream (SFR) budget',
    'budget_sfr_ts': 'Stream (SFR) budget over time',
    'budget_lak': 'Ponds (LAK) budget',
    'budget_lak_ts': 'Ponds (LAK) budget over time',
    'lake_stage': 'Pond stages against their beds and rims',
    'lake_budget_years_total': 'All ponds: water budget per hydrological '
                               'year [m3]',
    'lake_budget_years_by_pond': 'Each pond: water budget per hydrological '
                                 'year [m3]',
    'lake_budget_years_total_mm': 'All ponds: water budget per hydrological '
                                  'year [mm over the pond cells]',
    'lake_budget_years_by_pond_mm': 'Each pond: water budget per '
                                    'hydrological year [mm over its cells]',
    'obs_heads': 'Heads at the observation points: simulated and observed',
    'outlet_streamflow': 'Streamflow at the catchment outlet: simulated and '
                         'observed',
}
_CRIT = {'RMSE': 'root mean square error', 'RSR': 'RMSE-to-SD ratio',
         'NSE': 'Nash-Sutcliffe efficiency', 'r': 'correlation'}


def _stem(fname):
    """The name without the extension and the native suite's page suffix
    (``_1_1``)."""
    s = os.path.splitext(fname)[0]
    return re.sub(r'_\d+_\d+$', '', s)


def classify(sub, fname):
    """The tab key of one figure (``sub`` is its folder in the run)."""
    s = _stem(fname)
    if sub == '_input':
        return 'input'
    if (s.startswith(('lake_stage', 'lake_budget')) or
            s.startswith('budget_lak')):
        return 'ponds'
    if ('_sp_plt_GWmap_' in s or '_sp_plt_MMmap_' in s
            or '_sp_plt_SWmap_' in s or s.startswith('05_map_')):
        return 'output'
    if ('calibcrit' in s or s in ('obs_heads', 'outlet_streamflow')
            or s.startswith('sm_depth_')):
        return 'calib'
    if ('WBsankey' in s or s in ('budget_compartment', 'budget_uzf',
                                 'budget_sfr', '03_wb_totals',
                                 '04_wb_hydroyear', '00_summary')):
        return 'budget'
    if (re.match(r'_0[A-Za-z0-9]+_ts', s) or s in (
            'budget_sfr_ts', '01_wb_timeseries', '02_wb_cumulative',
            '06_heads', '07_coupling')):
        return 'series'
    return 'other'


def title(sub, fname):
    """A figure's title for the panel: what it shows, in words."""
    s = _stem(fname)
    m = re.match(r'(?:_sp_plt_)?IN_\d+_(.+)$', s)
    if m:
        key = m.group(1)
        v = re.match(r'VEG(\d+)area$', key)
        if v:
            return 'Vegetation type %s: cover' % v.group(1)
        return _IN.get(key, key)
    m = re.match(r'_sp_plt_MMmap_(.+)$', s)
    if m:
        return 'Soil column: %s' % _MM.get(m.group(1), m.group(1))
    m = re.match(r'_sp_plt_GWmap_head_series_day0*(\d+)$', s)
    if m:
        return 'Aquifer: head on day %s' % m.group(1)
    m = re.match(r'_sp_plt_GWmap_(.+)$', s)
    if m:
        return 'Aquifer: %s' % _GW.get(m.group(1), m.group(1))
    m = re.match(r'_sp_plt_SWmap_(.+)$', s)
    if m:
        return 'Streams: %s' % _SW.get(m.group(1), m.group(1))
    m = re.match(r'05_map_(.+)$', s)
    if m:
        return 'Map of %s, against the NWT run' % m.group(1)
    if s in _NWT:
        return _NWT[s]
    if s in _FIXED:
        return _FIXED[s]
    m = re.match(r'__plt_calibcrit(\w+?)_\d+$', s)
    if m:
        return 'Calibration criterion: %s' % _CRIT.get(m.group(1), m.group(1))
    m = re.match(r'sm_depth_(\w+)$', s)
    if m:
        return '%s: soil moisture at depth, simulated and observed' % m.group(1)
    m = re.match(r'_0(\w+?)_ts(_part2|GW_part3MF)?$', s)
    if m:
        part = {None: 'fluxes and state variables',
                '_part2': 'fluxes and state variables (2)',
                'GW_part3MF': 'groundwater (MODFLOW)'}[m.group(2)]
        return '%s: time series of %s' % (m.group(1), part)
    m = re.match(r'_(catchment|catchment_full|obs_(\w+))_WBsankey_0whole(_pc)?$',
                 s)
    if m:
        who = ('the catchment' if m.group(1) == 'catchment' else
               'the catchment, every flux' if m.group(1) == 'catchment_full'
               else 'point %s' % m.group(2))
        return 'Water balance of %s%s' % (who, ', in % of rainfall'
                                         if m.group(3) else ' [mm/y]')
    if s == '00_summary':
        return 'Water budget against the NWT run (table)'
    return s


# the order inside a tab: the catchment before the points, a family together
_ORDER = {
    'input': (r'IN_000_general', r'IN_000_model', r'IN_'),
    'output': (r'MMmap_', r'GWmap_head$', r'GWmap_head_series', r'GWmap_',
               r'SWmap_', r'05_map_'),
    'series': (r'01_wb', r'02_wb', r'06_heads', r'07_coupling',
               r'budget_sfr_ts', r'_0'),
    'calib': (r'outlet_streamflow', r'obs_heads', r'calibcrit', r'sm_depth_'),
    'budget': (r'00_summary', r'_catchment_WBsankey', r'_catchment_full',
               r'03_wb', r'04_wb', r'budget_compartment', r'budget_uzf',
               r'budget_sfr', r'_obs_'),
    'ponds': (r'lake_budget_years_total$', r'lake_budget_years_by_pond$',
              r'lake_budget_years_total_mm', r'lake_budget_years_by_pond_mm',
              r'budget_lak_ts', r'budget_lak$', r'lake_stage'),
}


def _rank(tab, sub, fname):
    # a point's three time-series sheets in their own order: 1, 2, groundwater
    s = _stem(fname).replace('_tsGW_part3MF', '_ts_part3')
    for k, pat in enumerate(_ORDER.get(tab, ())):
        if re.search(pat, s):
            return (k, s)
    return (len(_ORDER.get(tab, ())), s)


def panels(path):
    """How many map panels a figure holds side by side: what plotLAYER
    writes in the PNG's 'MM-panels' text (2026-10-09), else 1. A figure of
    two layers takes the whole row on the panel, so its maps show at the
    size of a one-layer map's."""
    try:
        from PIL import Image
        with Image.open(path) as im:
            return max(1, int(im.info.get('MM-panels', 1)))
    except Exception:                    # no PIL, not a PNG, no such key
        return 1


def rows(items, path_of=lambda it: it[2]):
    """The page's rows: a figure of several panels alone on its row, the
    others two to a row in their order. A one-panel figure left alone
    before a wide one gets a row of its own."""
    out, pend = [], []
    for it in items:
        if panels(path_of(it)) > 1:
            if pend:
                out.append(pend)
                pend = []
            out.append([it])
            continue
        pend.append(it)
        if len(pend) == 2:
            out.append(pend)
            pend = []
    if pend:
        out.append(pend)
    return out


def arrange(files):
    """``{tab: [(sub, fname, path, title)]}`` for ``(sub, fname, path)``
    items, each tab in its own order; tabs with nothing are left out."""
    out = {}
    for sub, fname, path in files:
        out.setdefault(classify(sub, fname), []).append(
            (sub, fname, path, title(sub, fname)))
    for tab, items in out.items():
        items.sort(key=lambda it: _rank(tab, it[0], it[1]))
    return out
