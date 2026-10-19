# -*- coding: utf-8 -*-
"""The water table MMsoil's groundwater ET starts from (user, 2026-10-19):
the surface cell's head, or -- that cell dry -- the first layer below that
is not, never above that layer's top (a piezometric level under a dry
layer is no water table); and MMsoil's dry floor, the bottom of the
column's deepest active layer, so a column is dry only when every layer
is. ETg from a deeper layer matters where roots are longer than the
surface layer is thick (not La Mata: layer 1 >= 20 m, roots <= 15 m)."""

import os
import sys
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in (TRUNK, os.path.join(TRUNK, 'app')):
    if p not in sys.path:
        sys.path.insert(0, p)

from marmites_coupler import MF6Coupler                         # noqa: E402

# three layers on one row of columns: top 100, bottoms 80 / 50 / 0
BOT = (80.0, 50.0, 0.0)
NCOL = 7


def _coupler(k_surface, lut=None):
    c = MF6Coupler.__new__(MF6Coupler)
    botm = np.array([[[b] * NCOL] for b in BOT])            # (3, 1, NCOL)
    c.mf6b = SimpleNamespace(botm=botm, nrow=1, ncol=NCOL)
    c.grid_kind = 'structured'
    c.ncpl = NCOL
    c.k_arr = np.asarray(k_surface, dtype=int)
    c.i_arr = np.zeros(NCOL, dtype=int)
    c.j_arr = np.arange(NCOL)
    c._below, c._bot_top = c._layers_below(lut)
    return c


def _x(layers):
    """The head array of the unreduced grid (user node = k * NCOL + j)."""
    return np.asarray([h for row in layers for h in row], dtype=float)


def test_the_water_table_is_the_first_layer_that_is_not_dry():
    #            wet   L2    L2 piezo  L3     all dry  L2 surf  L2 surf dry
    k_surf = [0, 0, 0, 0, 0, 1, 1]
    h1 = [90.0, 75.0, 75.0, 75.0, 75.0, 0.0, 0.0]
    h2 = [70.0, 60.0, 85.0, 45.0, 45.0, 60.0, 45.0]
    h3 = [60.0, 55.0, 55.0, 30.0, -5.0, 40.0, 30.0]
    c = _coupler(k_surf)
    c.p_x = _x([h1, h2, h3])
    surface = np.array([h1[j] if k_surf[j] == 0 else h2[j]
                        for j in range(NCOL)])
    got = c.water_table(surface)
    assert got.tolist() == [90.0,     # layer 1 wet: its own head
                            60.0,     # layer 1 dry: layer 2's water table
                            80.0,     # layer 2 above its top: its top
                            30.0,     # layers 1 and 2 dry: layer 3's
                            -5.0,     # all dry: the deepest head (floor)
                            60.0,     # layer 2 is the surface and wet
                            30.0]     # layer 2 the surface, dry: layer 3
    # untouched: what MF6 reported
    assert surface[1] == 75.0


def test_an_inactive_layer_is_stepped_over():
    """Column 1's layer 2 is no active cell (absent from the reduced
    grid): under its dry layer 1 the water table is layer 3's."""
    user = [k * NCOL + j for k in range(3) for j in range(NCOL)]
    gone = 1 * NCOL + 1
    lut = {u: r for r, u in enumerate(u for u in user if u != gone)}
    c = _coupler([0] * NCOL, lut=lut)
    h = {0: [90.0] * NCOL, 1: [70.0] * NCOL, 2: [60.0] * NCOL}
    h[0][1], h[2][1] = 75.0, 33.0
    c.p_x = np.asarray([h[k][j] for k in range(3) for j in range(NCOL)
                        if k * NCOL + j != gone], dtype=float)
    got = c.water_table(np.asarray(h[0], dtype=float))
    assert got[1] == 33.0 and got[0] == 90.0


def test_a_single_layer_model_is_left_alone():
    c = MF6Coupler.__new__(MF6Coupler)
    c._below, c._bot_top = [], None
    h = np.array([1.0, 2.0])
    assert c.water_table(h) is h


def test_the_dry_floor_is_the_deepest_active_bottom():
    """run_lamata_mf6._dry_floor: the deepest active layer's bottom, also
    where layer 1 pinches out (the outcrop is layer 2) or layer 3 is
    absent."""
    from lib import runs as runlib
    runner = runlib.import_runner()
    ib = np.ones((3, 1, 3), dtype=int)
    ib[0, 0, 1] = 0                 # layer 1 absent: outcrop layer 2
    ib[2, 0, 2] = 0                 # layer 3 absent
    cMF = SimpleNamespace(botm=np.array([[[80.0] * 3], [[50.0] * 3],
                                         [[0.0] * 3]]),
                          ibound=ib, outcropL=np.array([[1, 2, 1]]))
    assert runner._dry_floor(cMF).tolist() == [[0.0, 0.0, 50.0]]
