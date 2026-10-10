# -*- coding: utf-8 -*-
"""The Results panel by tabs (user, 2026-10-09): input maps, output maps,
time series, calibration (state variables), the water budgets and the
ponds' water budget. Every figure a run writes has a tab and a title."""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
spec = importlib.util.spec_from_file_location(
    '_app_results', os.path.join(CODE, 'app', 'lib', 'results.py'))
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)

# what run 20261008114055 wrote (one of each family), with its tab
RUN = {
    ('_input', 'IN_000_general_map.png'): 'input',
    ('_input', 'IN_000_model_map.png'): 'input',
    ('_input', '_sp_plt_IN_007_hk_1_1.png'): 'input',
    ('_input', '_sp_plt_IN_022_VEG1area_1_1.png'): 'input',
    ('_output', '_sp_plt_MMmap_ETsoil_1_1.png'): 'output',
    ('_output', '_sp_plt_GWmap_head_1_1.png'): 'output',
    ('_output', '_sp_plt_GWmap_head_series_day00073_1_1.png'): 'output',
    ('_output', '_sp_plt_GWmap_dSg_1_1.png'): 'output',
    ('figures_nwt_comparison', '05_map_ETg.png'): 'output',
    ('_output', '_0P0_ts.png'): 'series',
    ('_output', '_0P0_ts_part2.png'): 'series',
    ('_output', '_0P0_tsGW_part3MF.png'): 'series',
    ('_output', 'budget_sfr_ts.png'): 'series',
    ('figures_nwt_comparison', '01_wb_timeseries.png'): 'series',
    ('figures_nwt_comparison', '07_coupling.png'): 'series',
    ('_output', '__plt_calibcritNSE_0.png'): 'calib',
    ('_output', 'obs_heads.png'): 'calib',
    ('_output', 'outlet_streamflow.png'): 'calib',
    ('_output', 'sm_depth_C1.png'): 'calib',
    ('_output', '_catchment_WBsankey_0whole.png'): 'budget',
    ('_output', '_obs_C1_WBsankey_0whole_pc.png'): 'budget',
    ('_output', 'budget_compartment.png'): 'budget',
    ('_output', 'budget_uzf.png'): 'budget',
    ('figures_nwt_comparison', '03_wb_totals.png'): 'budget',
    ('_output', 'lake_budget_years_by_pond_mm.png'): 'ponds',
    ('_output', 'lake_budget_years_total.png'): 'ponds',
    ('_output', 'budget_lak_ts.png'): 'ponds',
    ('_output', 'lake_stage.png'): 'ponds',
}


def test_every_figure_has_its_tab():
    got = {k: R.classify(*k) for k in RUN}
    assert got == RUN
    assert [k for k, _l in R.TABS][:6] == ['input', 'output', 'series',
                                           'calib', 'budget', 'ponds']


def test_an_unknown_figure_is_shown_under_other_not_hidden():
    assert R.classify('_output', 'something_new.png') == 'other'
    assert R.title('_output', 'something_new.png') == 'something_new'


def test_every_known_figure_has_a_title_in_words():
    for sub, fname in RUN:
        t = R.title(sub, fname)
        assert t and t != os.path.splitext(fname)[0], fname
    assert R.title('_input', '_sp_plt_IN_007_hk_1_1.png') == \
        'Horizontal hydraulic conductivity'
    assert R.title('_output', '_sp_plt_GWmap_head_series_day00073_1_1.png') \
        == 'Aquifer: head on day 73'
    assert R.title('_output', '_obs_C1_WBsankey_0whole_pc.png') == \
        'Water balance of point C1, in % of rainfall'
    assert R.title('_output', 'lake_budget_years_by_pond_mm.png').endswith(
        '[mm over its cells]')


def test_the_stream_and_exchange_maps_have_their_tab_and_title():
    """User, 2026-10-10: the streamflow on a map, the exchange with the
    streams and the ponds, the flow across the upper face -- and the
    recharge labels with Exf and ETg as the amounts their maps show."""
    q = ('_output', '_sp_plt_SWmap_Q_1_1.png')
    assert R.classify(*q) == 'output'
    assert R.title(*q) == 'Streams: mean streamflow (log10 of 1 + Q, Q in m3/d)'
    for stem, words in (('FUF', 'upper face'), ('FLF', 'lower face'),
                        ('SFRg', 'streams'), ('LAKg', 'ponds')):
        t = R.title('_output', '_sp_plt_GWmap_%s_1_1.png' % stem)
        assert t.startswith('Aquifer: ') and words in t, t
    assert R.title('_output', '_sp_plt_GWmap_Re_1_1.png') == \
        'Aquifer: effective recharge (Rg - Exf)'
    assert R.title('_output', '_sp_plt_GWmap_Rn_1_1.png') == \
        'Aquifer: net recharge (Rg - Exf - ETg)'
    order = R.arrange([('_output', '_sp_plt_SWmap_Q_1_1.png', 'q.png'),
                       ('_output', '_sp_plt_GWmap_Rg_1_1.png', 'rg.png'),
                       ('_output', '_sp_plt_MMmap_Ro_1_1.png', 'ro.png')])
    assert [it[1] for it in order['output']][-1].startswith(
        '_sp_plt_SWmap_'), order['output']


def test_the_soil_storage_series_follows_the_head_series():
    """User, 2026-10-10: the soil storage on the head series' days, right
    after them in the tab, titled by its day."""
    f = '_sp_plt_MMmap_Ssoil_series_day00200_1_1.png'
    assert R.classify('_output', f) == 'output'
    assert R.title('_output', f) == 'Soil column: soil water storage on day 200'
    names = ['_sp_plt_MMmap_Ro_1_1.png', '_sp_plt_GWmap_head_1_1.png',
             '_sp_plt_GWmap_head_series_day00001_1_1.png', f,
             '_sp_plt_GWmap_Rg_1_1.png']
    got = [it[1] for it in R.arrange([('_output', n, n) for n in names])
           ['output']]
    assert got == names, got


def test_a_two_layer_map_takes_the_whole_row(tmp_path):
    """User, 2026-10-09: a map of two layers (side by side) occupies both
    columns, a one-layer map one -- so their maps show at the same size."""
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    paths = {}
    for name, n in (('a', 1), ('b', 1), ('c', 1), ('wide', 2), ('d', 1)):
        fig = plt.figure(figsize=(1, 1))
        fn = str(tmp_path / ('%s.png' % name))
        fig.savefig(fn, metadata={'MM-panels': str(n)})
        plt.close(fig)
        paths[name] = fn
    plain = str(tmp_path / 'plain.png')
    fig = plt.figure(figsize=(1, 1))
    fig.savefig(plain)
    plt.close(fig)
    assert R.panels(paths['wide']) == 2 and R.panels(paths['a']) == 1
    assert R.panels(plain) == 1, 'a figure without the key is one panel'
    items = [('_output', n + '.png', paths[n], n)
             for n in ('a', 'b', 'c', 'wide', 'd')]
    rows = [[it[3] for it in r] for r in R.rows(items)]
    assert rows == [['a', 'b'], ['c'], ['wide'], ['d']]


def test_inside_a_tab_the_catchment_comes_before_the_points():
    files = [(s, f, f) for s, f in RUN]
    g = R.arrange(files)
    series = [f for _s, f, _p, _t in g['series']]
    assert series.index('01_wb_timeseries.png') < series.index('_0P0_ts.png')
    assert series[-3:] == ['_0P0_ts.png', '_0P0_ts_part2.png',
                           '_0P0_tsGW_part3MF.png'], series
    ponds = [f for _s, f, _p, _t in g['ponds']]
    assert ponds[0] == 'lake_budget_years_total.png'
    inputs = [f for _s, f, _p, _t in g['input']]
    assert inputs[:2] == ['IN_000_general_map.png', 'IN_000_model_map.png']
