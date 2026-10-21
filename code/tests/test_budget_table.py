# -*- coding: utf-8 -*-
"""The water-budget table at the end of a run (the Run panel's log): MF6
against the NWT reference, the difference in mm/year AND in percent of the
reference (user, 2026-10-21)."""

import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def _pwb():
    spec = importlib.util.spec_from_file_location(
        '_pwb_table', os.path.join(HERE, 'plot_water_budget.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_difference_in_percent_of_the_reference():
    pwb = _pwb()
    assert pwb.pct_diff('Ro', 96.4, 139.2) == '-30.7 %'
    assert pwb.pct_diff('I', 344.2, 323.8) == '+6.3 %'
    # a storage change, a reference under 1 mm/year: no percentage
    assert pwb.pct_diff('dSsoil', 0.1, -1.3) == '-'
    assert pwb.pct_diff('Eow', 1.6, 0.4) == '-'
    assert pwb.pct_diff('Qout', 153.6, None) == '-'


def test_the_table_has_both_columns(tmp_path, capsys):
    pwb = _pwb()
    new = {'wb_ts': np.zeros((5, 3)), 'outer_iters': np.array([3, 5]),
           'rejinf': np.zeros((5, 2)), 'iters_kind': 'linear'}
    pwb.write_summary(['P', 'Ro', 'dSsoil'], [450.4, 96.4, 0.1],
                      [450.6, 139.2, -1.3], str(tmp_path), new,
                      sw={'gw_to_sfr': 68.1, 'sfr_to_gw': 6.4})
    out = capsys.readouterr().out
    head = [ln for ln in out.splitlines() if ln.startswith('flux')][0]
    assert head.split()[-2:] == ['diff', '%'] and 'NWT (Picard)' in head
    row = [ln for ln in out.splitlines() if ln.startswith('Ro ')][0]
    assert row.split()[-3:] == ['-42.8', '-30.7', '%']
    assert [ln for ln in out.splitlines()
            if ln.startswith('dSsoil')][0].split()[-1] == '-'
    sfr = [ln for ln in out.splitlines() if ln.startswith('Ro+netSFR')][0]
    assert '+13.6 %' in sfr                    # 158.1 against 139.2
