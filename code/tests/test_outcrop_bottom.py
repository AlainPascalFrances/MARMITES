# -*- coding: utf-8 -*-
"""MMsoil's dry-cell test reads the bottom of the OUTCROP layer.

MMsoil takes a surface cell for dry when its head is below ``botm_l0`` and
then puts the water table at that bottom. The driver passed layer 1's bottom
everywhere: where layer 1 pinches out (thickness 0, idomain 0) that bottom is
the soil base itself, so every water table in layer 2 read as dry AT THE
SOIL BASE -- Eg and Tg at their full potential 1-6 m above the real table
(La Mata, 181 cells near the streams, 2026-10-06).
"""
import importlib.util
import os
import sys
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, HERE, os.path.join(CODE, 'ppMF6'),
           os.path.join(CODE, 'MARMITESutilities'),
           os.path.join(CODE, 'MARMITESsoil'), os.path.join(CODE, 'ppMF_FloPy')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _driver():
    spec = importlib.util.spec_from_file_location(
        '_mm_runner_outcrop', os.path.join(HERE, 'run_lamata_mf6.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_bottom_of_each_cells_outcrop_layer():
    rl = _driver()
    # 1 x 3 cells, 3 layers: layer 1 present / pinched out / all inactive
    botm = np.array([[[100.0, 100.0, 100.0]],
                     [[80.0, 70.0, 60.0]],
                     [[50.0, 40.0, 30.0]]])
    cMF = SimpleNamespace(botm=botm, outcropL=np.array([[1, 2, 0]]))
    b = rl._outcrop_bottom(cMF)
    assert b.shape == (1, 3)
    # layer 1's bottom, layer 2's bottom, and layer 1's for an inactive cell
    assert b.tolist() == [[100.0, 70.0, 100.0]]


def test_a_layer_2_water_table_is_not_dry():
    """The head of a layer-2 surface cell, below layer 1's (pinched-out)
    bottom but above its own: not dry, MMsoil keeps the real head."""
    rl = _driver()
    botm = np.array([[[100.0]], [[70.0]]])
    b = rl._outcrop_bottom(SimpleNamespace(botm=botm,
                                           outcropL=np.array([[2]])))
    h = 96.5                                   # 3.5 m below the soil base
    assert not h < b[0, 0]
    assert h < botm[0, 0, 0]                   # what the old test said: dry


def test_the_driver_passes_it_everywhere():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert 'botm_l0 = np.asarray(cMF.botm)[0]' not in src
    assert src.count('botm_l0 = _outcrop_bottom(cMF)') == 3
