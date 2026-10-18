# -*- coding: utf-8 -*-
"""The UZF water content carried into a run, relative to the water contents
it was left with (2026-10-18).

A state saved with thtr/thts 0.44/0.45 (the reference model: thts 0.45,
Sy 0.01) restarting a PEST realisation whose thtr/thts are 0.275/0.30 was
clipped to thts: the unsaturated zone started full and drained into a
recharge pulse of metres of head. What is carried now is the FILL of the
mobile range, (wc - thtr0) / (thts0 - thtr0), mapped onto the new range."""

import importlib.util
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
for p in (TRUNK, os.path.join(TRUNK, 'ppMF6'),
          os.path.join(TRUNK, 'MARMITESutilities')):
    if p not in sys.path:
        sys.path.insert(0, p)

mf6 = pytest.importorskip('marmites_mf6')


def test_the_fill_of_the_mobile_range_is_carried():
    ref = (np.full(4, 0.44), np.full(4, 0.45))
    carried = np.array([0.44, 0.445, 0.45, np.nan])
    thtr, thts = np.full(4, 0.275), np.full(4, 0.30)
    thti, n = mf6.carry_thti(carried, ref, thtr, thts, np.full(4, 0.275))
    assert n == 3
    assert thti == pytest.approx([0.275, 0.2875, 0.30, 0.275])


def test_without_a_reference_the_content_is_clipped():
    """The old behaviour, for a state that does not say what it was left
    with and no reference to stand in for it."""
    thti, n = mf6.carry_thti(np.array([0.445, 0.20, np.nan]), None,
                             np.full(3, 0.275), np.full(3, 0.30),
                             np.full(3, 0.28))
    assert n == 2 and thti == pytest.approx([0.30, 0.275, 0.28])


def test_a_reference_without_a_mobile_range_falls_back_to_clipping():
    thti, _n = mf6.carry_thti(np.array([0.29]), (np.array([0.3]),
                                                 np.array([0.3])),
                              np.array([0.275]), np.array([0.30]),
                              np.array([0.275]))
    assert thti == pytest.approx([0.29])


def test_the_same_water_contents_carry_unchanged():
    """The spin-up's cycles carry onto the model they came from."""
    rng = np.random.default_rng(3)
    thtr = rng.uniform(0.1, 0.3, (2, 1, 5))
    thts = thtr + rng.uniform(0.005, 0.05, (2, 1, 5))
    wc = thtr + rng.uniform(0, 1, (2, 1, 5)) * (thts - thtr)
    thti, _n = mf6.carry_thti(wc, (thtr, thts), thtr, thts, thtr.copy())
    assert np.allclose(thti, wc)


# --------------------------------------------------- the saved state says it

def _rl():
    spec = importlib.util.spec_from_file_location(
        '_mm_runner_carry', os.path.join(HERE, 'run_lamata_mf6.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_saved_state_keeps_the_water_contents_it_was_left_with(tmp_path):
    rl = _rl()
    nlay, nrow, ncol, ncell, nsl = 2, 3, 1, 3, 2
    shape = (nlay, nrow, ncol)
    b = SimpleNamespace(nlay=nlay, nrow=nrow, ncol=ncol,
                        thti_from_wc=lambda wc: np.full(shape, 0.445),
                        uzf_thtr_built=np.full(shape, 0.44),
                        uzf_thts_built=np.full(shape, 0.45))
    cpl = SimpleNamespace(uzf_wc_final=np.zeros(4), carry_out={})
    st = SimpleNamespace(Ssoil_ini=np.zeros((ncell, nsl)))
    ctx = SimpleNamespace(ncell=ncell, _nslmax=nsl)
    rl.save_run_state(str(tmp_path / 'hi'), b, cpl, st)
    got = rl.load_run_state(str(tmp_path / 'hi'), b, ctx)
    r0, s0 = got['uzf_ref']
    assert np.allclose(r0, 0.44) and np.allclose(s0, 0.45)
    # an older state, without them: no reference from the state itself
    b2 = SimpleNamespace(nlay=nlay, nrow=nrow, ncol=ncol,
                         thti_from_wc=lambda wc: np.full(shape, 0.445))
    rl.save_run_state(str(tmp_path / 'old'), b2, cpl, st)
    assert rl.load_run_state(str(tmp_path / 'old'), b2, ctx)['uzf_ref'] \
        is None


def test_the_driver_carries_against_the_state_or_the_reference_model():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert ("b.uzf_thti_carry_ref = (saved_state.get('uzf_ref')\n"
            "                                    or getattr(cMF, "
            "'uzf_carry_ref', None))") in src.replace('\r\n', '\n')
    assert ('b.uzf_thti_carry_ref = (b.uzf_thtr_built, b.uzf_thts_built)'
            in src)
