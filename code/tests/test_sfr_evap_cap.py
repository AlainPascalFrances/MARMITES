# -*- coding: utf-8 -*-
"""The cap on the stream's open-water evaporation, against MODFLOW 6 itself.

Analysis §8.23: MF6 takes a reach's evaporation from the depth its previous
solve stored, so a reach receiving less than E0 x w x L over a water table
below its bed flip-flops between evaporating all of it and leaking all of
it. sfr_fn's perturbed derivative straddles the two states, and the solver
cannot converge at a tight outer_dvclose. The coupler caps EVAP at
sfr.evap_inflow_fraction of the reach's inflow (_cap_sfr_evap, with the
upstream and mover inflow copied after every step by _keep_sfr_inflow).

This drives the valley of diag_sfr_evap_flipflop.py through libmf6 in the
coupler's order -- prepare_time_step, EVAP written, do_time_step (ATS),
finalize_time_step, inflow kept -- with the coupler's OWN two methods, at
outer_dvclose 0.001. A toy MF6 model, not La Mata.
"""
import importlib.util
import os
import re
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


coup = _load('marmites_coupler_cap', os.path.join(CODE, 'marmites_coupler.py'))
diag = _load('diag_sfr_evap_flipflop',
             os.path.join(HERE, 'diag_sfr_evap_flipflop.py'))


def _libmf6():
    import mm_paths
    if not os.path.isfile(mm_paths.LIBMF6) or not os.path.isfile(mm_paths.MF6_EXE):
        pytest.skip('MF6 (libmf6 and mf6) not where panel 0 puts it')
    pytest.importorskip('modflowapi')
    pytest.importorskip('flopy')
    return mm_paths.LIBMF6, mm_paths.MF6_EXE


class _Cpl(coup.MF6Coupler):
    """The coupler's stream-evaporation writer, bound to a real SFR."""

    def __init__(self, api, frac):
        pre = 'TOY/SFR'
        self.conv_fact = 1000.0
        self.Eo_zonesSP = np.array([[diag.E0 * 1000.0]])          # mm/d
        self.p_sfr_evap = api.get_value_ptr(pre + '/EVAP')
        self.nreaches = self.p_sfr_evap.size
        self.sfr_evap_zone = np.zeros(self.nreaches, dtype=int)
        self.nlakes = 0
        self.p_lak_evap = None
        self.p_sfr_inflow = api.get_value_ptr(pre + '/INFLOW')
        self.p_sfr_usflow = api.get_value_ptr(pre + '/USFLOW')
        self.p_sfr_qfrommvr = None
        self.sfr_wl = (np.array(api.get_value_ptr(pre + '/LENGTH'), float)
                       * np.array(api.get_value_ptr(pre + '/WIDTH'), float))
        self.sfr_evap_frac = frac
        self.sfr_up_prev = None
        self.sfr_evap_capped = self.sfr_evap_writes = 0


def _march(ws, libmf6, capped, frac=0.5):
    """Returns (failed attempts, worst EVAP x w L / inflow over the run)."""
    from modflowapi import ModflowApi
    api = ModflowApi(libmf6, working_directory=str(ws))
    api.initialize(os.path.join(str(ws), 'mfsim.nam'))
    worst = 0.0
    try:
        c = _Cpl(api, frac)
        while api.get_current_time() < api.get_end_time() - 1e-9:
            api.prepare_time_step(api.get_time_step())
            if capped:
                c._write_openwater_evap(0)
            else:
                c.p_sfr_evap[:] = diag.E0
            qin = np.array(c.p_sfr_inflow, float) + (
                0.0 if c.sfr_up_prev is None else c.sfr_up_prev)
            wet = qin > 0
            if wet.any():
                worst = max(worst, float(np.max(c.p_sfr_evap[wet] * c.sfr_wl[wet]
                                                / qin[wet])))
            api.do_time_step()
            api.finalize_time_step()
            c._keep_sfr_inflow()
    finally:
        api.finalize()
    lst = open(os.path.join(str(ws), 'mfsim.lst'), errors='replace').read()
    return len(re.findall(r'Solution 1 did not converge', lst)), worst


def test_the_capped_stream_converges_at_a_tight_tolerance(tmp_path):
    libmf6, exe = _libmf6()
    ws = tmp_path / 'capped'
    ws.mkdir()
    diag.build(str(ws), exe)
    diag.set_dvclose(str(ws), 0.001)
    nfail, worst = _march(ws, libmf6, capped=True)
    assert worst <= 0.5 + 1e-9, 'a reach evaporated more than half its inflow'
    assert nfail == 0, ('%d failed attempts at outer_dvclose 0.001 with the '
                        'cap' % nfail)


def test_with_the_streams_not_evaporating_it_converges(tmp_path):
    """sfr.evap_inflow_fraction = 0 (8.23.2): no evaporation, no
    flip-flop."""
    libmf6, exe = _libmf6()
    ws = tmp_path / 'zero'
    ws.mkdir()
    diag.build(str(ws), exe)
    diag.set_dvclose(str(ws), 0.001)
    nfail, worst = _march(ws, libmf6, capped=True, frac=0.0)
    assert worst == 0.0, 'a reach evaporated with the fraction at 0'
    assert nfail == 0, '%d failed attempts with no stream evaporation' % nfail


def test_without_the_cap_mf6_still_flip_flops(tmp_path):
    """Pins the reason for the cap. If this starts failing, MF6 no longer
    flip-flops on a nearly dry reach (analysis §8.23) and the cap may no
    longer be needed."""
    libmf6, exe = _libmf6()
    ws = tmp_path / 'uncapped'
    ws.mkdir()
    diag.build(str(ws), exe)
    diag.set_dvclose(str(ws), 0.001)
    nfail, worst = _march(ws, libmf6, capped=False)
    assert worst > 1.0, 'the toy no longer has reaches offered more than they get'
    assert nfail > 0, 'MF6 converged without the cap: see analysis §8.23'
