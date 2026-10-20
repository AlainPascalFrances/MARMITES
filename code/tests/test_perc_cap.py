# -*- coding: utf-8 -*-
"""The coupler caps each cell's soil percolation at what its land UZF object
can take (2026-10-18, analysis doc 8.26):

    cap = min(VKS, max(top - surfdep - h0, 0) * (thts - thtr) / day)

MMsoil keeps the rest, and what UZF still refuses runs off. Without the
cap the refused water came back into the soil from below the next day and
every rejecting cell flip-flopped day by day."""

import importlib.util
import os
import sys
import types

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_SOIL = os.path.join(HERE, '..', 'MARMITESsoil')
if os.path.abspath(_SOIL) not in sys.path:
    sys.path.insert(0, os.path.abspath(_SOIL))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CM = _load('t_coupler_mock_cap', os.path.join(HERE, 'test_coupler_mock.py'))
coup = CM.coup


def _self(vks, sy, top, surfdep=0.25, perlen=(1.0,)):
    n = len(vks)
    return types.SimpleNamespace(
        _uzf_vks=np.asarray(vks, float), _uzf_sy=np.asarray(sy, float),
        _uzf_top=np.asarray(top, float), _uzf_surfdep=surfdep, ncell=n,
        perlen=list(perlen), conv_fact=1000.0)


def test_the_cap_is_vks_or_the_room_above_the_water_table():
    s = _self(vks=[0.025, 0.025, 0.025, 0.002], sy=[0.01] * 4,
              top=[760.0] * 4)
    h = np.array([750.0, 759.5, 759.9, 750.0])
    cap = coup.MF6Coupler.perc_cap(s, h, 0)
    # deep water table: VKS (25 mm/d); 0.25 m of room x 0.01 = 2.5 mm;
    # above top - surfdep: no room; a low VKS wins over the room
    assert cap == pytest.approx([25.0, 2.5, 0.0, 2.0])


def test_a_head_that_is_not_one_leaves_vks():
    s = _self(vks=[0.01, 0.01, 0.01], sy=[0.01] * 3, top=[760.0] * 3)
    cap = coup.MF6Coupler.perc_cap(s, np.array([np.nan, -1e30, 1e30]), 0)
    assert cap == pytest.approx([10.0, 10.0, 10.0])


def test_the_room_is_spread_over_the_stress_period():
    s = _self(vks=[1.0], sy=[0.02], top=[760.0], perlen=(1.0, 5.0))
    assert coup.MF6Coupler.perc_cap(s, np.array([759.0]), 1) == \
        pytest.approx([0.02 * 0.75 / 5.0 * 1000.0])


def test_a_builder_without_uzf_properties_gives_no_cap():
    cpl, _api, _ctx = CM._setup(nper=1)
    assert cpl.perc_cap(np.zeros(cpl.ncell), 0) is None


def test_the_coupler_reads_the_builders_uzf_properties_and_tops():
    """A land object in layer 2 (layer 1 absent) takes the bottom of layer 1
    as its top."""
    cpl, api, ctx = CM._setup(nper=2)
    n = ctx.ncell
    b = cpl.mf6b
    nlay, nrow, ncol = 2, b.nrow, b.ncol
    b.top = np.full((nrow, ncol), 760.0)
    b.botm = np.stack([np.full((nrow, ncol), 740.0),
                       np.full((nrow, ncol), 700.0)])
    b.surf_cells = [(i, j, 1 if k % 2 else 0)
                    for k, (i, j, _l) in enumerate(b.surf_cells)]
    b.uzf_land_vks = np.full(n, 0.025)
    b.uzf_land_sy = np.full(n, 0.01)
    b.uzf_surfdep = 0.25
    cpl2 = coup.MF6Coupler(cpl.mm, ctx, cpl.state, b, conv_fact=1000.0)
    k = np.array([c[2] for c in b.surf_cells])
    assert np.allclose(cpl2._uzf_top[k == 0], 760.0)
    assert np.allclose(cpl2._uzf_top[k == 1], 740.0)
    cap = cpl2.perc_cap(np.full(n, 739.9), 0)
    assert np.allclose(cap[k == 1], 0.0)            # at its top: no room
    assert np.allclose(cap[k == 0], 25.0)           # 20 m down: VKS


def test_the_soil_step_receives_the_cap_every_day():
    cpl, api, ctx = CM._setup(nper=3)
    n, b = ctx.ncell, cpl.mf6b
    b.top = np.full((b.nrow, b.ncol), 760.0)
    b.uzf_land_vks = np.full(n, 0.025)
    b.uzf_land_sy = np.full(n, 0.01)
    b.uzf_surfdep = 0.25
    cpl2 = coup.MF6Coupler(cpl.mm, ctx, cpl.state, b, conv_fact=1000.0)
    seen = []
    step0 = cpl2.mm.step

    def spy(*a, **kw):
        seen.append(kw.get('perc_cap'))
        return step0(*a, **kw)
    cpl2.mm.step = spy
    cpl2.run(api)
    assert len(seen) == ctx.cMF.nper
    assert all(s is not None and s.shape == (n,) for s in seen)
    assert all(np.all((s >= 0.0) & (s <= 25.0)) for s in seen)


def test_a_full_unsaturated_zone_takes_nothing_more():
    """2026-10-19: the room is what the zone does not yet hold, thts - wc
    (WCNEW), not the empty zone's thts - thtr: a slow zone kept being fed
    to near saturation, and the water table rising into it stalled MF6's
    Newton loop (UZF's uz_rise has no derivative). Empty, half full, full,
    and a WCNEW outside [thtr, thts] held within it."""
    s = _self(vks=[1.0] * 4, sy=[0.02] * 4, top=[760.0] * 4)
    s._uzf_thtr = np.array([0.25] * 4)
    s.p_wcnew = np.array([0.25, 0.26, 0.27, 0.0])
    cap = coup.MF6Coupler.perc_cap(s, np.array([750.0] * 4), 0)
    room = 9.75                         # m between top - surfdep and h0
    assert cap == pytest.approx([room * 0.02 * 1000.0, room * 0.01 * 1000.0,
                                 0.0, room * 0.02 * 1000.0])
    # no water content (or no thtr): the empty zone, as before
    assert coup.uzf_room_fraction([0.02], None, [0.25], 1).tolist() == [0.02]
    assert coup.uzf_room_fraction([0.02], [0.26], None, 1).tolist() == [0.02]


def test_the_builder_gives_the_land_objects_thtr():
    import inspect
    src = inspect.getsource(coup.MF6Coupler.__init__)
    assert "'uzf_land_thtr'" in src
    with open(os.path.join(HERE, '..', 'ppMF6', 'marmites_mf6.py'),
              encoding='utf-8') as fh:
        assert 'self.uzf_land_thtr = ' in fh.read()


def test_the_wetting_front_saturation_limits_the_cap():
    """VKS * s^eps (2026-10-20): the flux whose Brooks-Corey front is s
    saturated -- at VKS a water table rising into the saturated front
    stalls MF6; a fixed fraction of VKS does not hold across eps (toy: 75 %
    gave 42 Newton iterations a day at eps 3.5, 94 at 14; s 0.92 39-42)."""
    s = _self(vks=[0.025] * 3, sy=[0.01] * 3, top=[760.0] * 3)
    s._uzf_eps = np.array([3.5, 7.0, 14.0])
    s.max_front_saturation = 0.92
    cap = coup.MF6Coupler.perc_cap(s, np.full(3, 700.0), 0)
    assert cap == pytest.approx(25.0 * 0.92 ** np.array([3.5, 7.0, 14.0]))
    assert cap[0] == pytest.approx(18.7, abs=0.05)    # ~75 % at eps 3.5
    # the room left still wins where it is smaller
    near = coup.MF6Coupler.perc_cap(s, np.full(3, 759.5), 0)
    assert near == pytest.approx([2.5, 2.5, min(2.5, 25.0 * 0.92 ** 14)])
    s.max_front_saturation = 1.0                       # up to VKS
    assert coup.MF6Coupler.perc_cap(s, np.full(3, 700.0), 0) == \
        pytest.approx([25.0] * 3)


def test_the_saturation_limit_is_a_uzf_field_passed_to_the_coupler():
    import inspect
    sys.path.insert(0, os.path.join(HERE, '..'))
    import marmites_config as mcfg
    c = mcfg.RunConfig.from_dict({})
    assert c.uzf.max_front_saturation == 0.92
    for bad in (0.0, 1.2):
        c.uzf.max_front_saturation = bad
        assert any('uzf.max_front_saturation must be in (0, 1]' in e
                   for e in c.problems())
    assert "'uzf_land_eps'" in inspect.getsource(coup.MF6Coupler.__init__)
    with open(os.path.join(HERE, '..', 'ppMF6', 'marmites_mf6.py'),
              encoding='utf-8') as fh:
        assert 'self.uzf_land_eps = ' in fh.read()
    with open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8') as fh:
        src = fh.read()
    assert 'uzf_front_saturation=float(cfg.uzf.max_front_saturation)' in src
    assert 'cpl.max_front_saturation = float(a.uzf_front_saturation)' in src
