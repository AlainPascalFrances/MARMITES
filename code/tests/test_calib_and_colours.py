# -*- coding: utf-8 -*-
"""Two post-processing failures of the run of 2026-09-23.

* 'h: error' at C1, C2, C3: quarterly piezometers, ONE observation in a
  60-day run, so RSR (rmse / std of the observations) was left None and the
  '%.2f' print raised; the RSR calibration figure then failed on the same
  None ("'>=' not supported between float and NoneType").
* 'GWmap map head_series skipped: 257 color bins ... ncolors = 256': the
  colour levels could outnumber the colours over a wide head range.
"""

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (os.path.join(CODE, 'MARMITESutilities'),
           os.path.join(CODE, 'MARMITESutilities', 'MARMITESplot'), CODE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

HN = 9999.999


@pytest.fixture(scope='module')
def proc():
    mod = pytest.importorskip('MARMITESprocess_v3')
    return mod


def test_one_observation_gives_rmse_and_marks_the_rest_undefined(proc):
    p = proc.clsPROCESS.__new__(proc.clsPROCESS)
    rmse, rsr, nse, r = p.compCalibCrit(np.array([752.6]),
                                        np.array([752.7]), HN)
    assert rmse == pytest.approx(0.1)
    assert rsr == HN and nse == HN and r == HN


def test_undefined_criteria_print_as_na(proc):
    assert proc._crit(HN, HN) == 'n/a'
    assert proc._crit(None, HN) == 'n/a'
    assert proc._crit(2.4412, HN) == '2.44'
    assert proc._crit(1.23, HN, '%.1f') == '1.2'


def test_a_single_observation_point_is_reported_not_errored(proc, capsys):
    p = proc.clsPROCESS.__new__(proc.clsPROCESS)
    p.compCalibCritObs(np.zeros((3, 1)), np.array([752.6, 752.6, 752.6]),
                       np.array([[HN, HN, HN]]),
                       np.array([HN, 752.7, HN]), HN, 'C3', 1)
    out = capsys.readouterr().out
    assert 'h: 0.10 m / n/a / n/a / n/a' in out, out
    assert 'error' not in out


def test_the_calibration_axis_ignores_undefined_criteria():
    mmplot = pytest.importorskip('MARMITESplot_v3')
    assert mmplot._defined_max([0.5, HN, None, 2.0], HN) == 2.0
    assert mmplot._defined_max([HN, None], HN) == 0.0


def test_colour_levels_never_outnumber_the_colours():
    mmplot = pytest.importorskip('MARMITESplot_v3')
    import matplotlib
    cmap = matplotlib.colormaps['Blues']
    for lo, hi in ((733.1, 811.2), (0.0, 1.0), (-3.0, 7.0), (1e-5, 3e-5)):
        levels = mmplot._colour_levels(cmap, lo, hi)
        assert len(levels) - 1 <= cmap.N, (lo, hi, len(levels))
        matplotlib.colors.BoundaryNorm(levels, cmap.N, extend='both')
