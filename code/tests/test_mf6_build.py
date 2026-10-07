# -*- coding: utf-8 -*-
"""Phase-3a tests: clsMF6 builds a loadable MODFLOW 6 simulation from the
real La Mata configuration (no mf6 binary required -- write + flopy reload)."""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

flopy = pytest.importorskip('flopy')
import matplotlib  # noqa: E402
matplotlib.use('agg')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mf6mod = _load('marmites_mf6', os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'))


@pytest.fixture(scope='module')
def cmf():
    # La Mata's model description as the run builds it -- no parameter file
    # (lamata_model; it derives the outcrop layer too)
    import lamata_model
    c = lamata_model.lamata_cmf()
    # short synthetic time base: 3 daily SPs (build test only)
    c.nper, c.perlen, c.nstp = 3, [1, 1, 1], [1, 1, 1]
    return c


def _build(cmf, tmp, daily=True):
    top = np.asarray(cmf.elev, dtype=float)
    botm = np.asarray(cmf.botm, dtype=float)
    b = mf6mod.clsMF6(cmf, top=top, botm=botm, sim_ws=str(tmp), daily=daily)
    b.build()
    return b


def test_build_and_write(cmf, tmp_path):
    b = _build(cmf, tmp_path)
    assert b.ncell == int((cmf.outcropL > 0).sum()) == 1954
    b.write()
    name = cmf.modelname.lower()
    files = os.listdir(tmp_path)
    for needed in ('mfsim.nam', f'{name}.nam', f'{name}.dis', f'{name}.npf',
                   f'{name}.sto', f'{name}.evt', f'{name}.drn', f'{name}.uzf',
                   f'{name}.ic', f'{name}.oc'):
        assert needed in files, f'{needed} missing from {sorted(files)[:12]}...'
    # no GHB on La Mata (ghb_yn=0); no WEL: the ETg wells went with the
    # WEL route (2026-10-07), groundwater ET is the two EVT packages
    assert f'{name}.ghb' not in files
    assert f'{name}.wel' not in files
    assert sum(f.endswith('.evt') for f in files) == 2


def test_reload_roundtrip(cmf, tmp_path):
    b = _build(cmf, tmp_path)
    b.write()
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(tmp_path), verbosity_level=0)
    gwf = sim.get_model(cmf.modelname.lower())
    # grid + periods: steady SP + 3 daily
    assert sim.tdis.nper.get_data() == 4
    dis = gwf.get_package('dis')
    # the layer count is the panel's (layers.nlay), not a literal
    assert (dis.nlay.get_data(), dis.nrow.get_data(),
            dis.ncol.get_data()) == (cmf.nlay, 65, 60)
    # newton on
    assert gwf.name_file.newtonoptions.get_data() is not None
    # UZF: land cells first, count matches active columns
    uzf = gwf.get_package('uzf')
    pak = uzf.packagedata.get_data()
    assert uzf.nuzfcells.get_data() == b.nuzfcells >= b.ncell
    land = pak[:b.ncell]
    assert all(int(r[2]) == 1 for r in land), 'first ncell UZF objects must be landflag=1'
    assert all(int(r[2]) == 0 for r in pak[b.ncell:]), 'subsurface objects must be landflag=0'
    # EVT: two packages (Eg, Tg), one record per surface cell, inert as
    # written -- the coupler writes each day's curve
    for pname in ('evt_eg', 'evt_tg'):
        evt = gwf.get_package(pname)
        assert evt.maxbound.get_data() == b.ncell
    assert gwf.get_package('wel') is None
    # DRN present with the legacy count
    drn = gwf.get_package('drn')
    assert drn.maxbound.get_data() == len(cmf.layer_row_column_elevation_cond[0])


def test_vertical_columns_chain(cmf, tmp_path):
    b = _build(cmf, tmp_path)
    # every land cell with active cells below must chain via ivertcon
    pak = {int(r[0]): r for r in b.uzf_packagedata}
    n_chained = 0
    for n, (i, j, k) in enumerate(b.surf_cells[:200]):
        rec = pak[n]
        ivert = int(rec[3])
        below_active = any(b.idomain[kk, i, j] > 0 for kk in range(k + 1, b.nlay))
        if below_active:
            assert ivert >= b.ncell, (n, ivert)
            n_chained += 1
            # the child must sit in the same (i, j) column
            child = pak[ivert]
            assert (int(child[1][1]), int(child[1][2])) == (i, j)
        else:
            assert ivert == -1
    assert n_chained > 0


def test_k33_ratio_semantics(cmf, tmp_path):
    # layvka=1 on La Mata -> stored VKA is the hk/vk ratio -> k33 = hk/vka
    b = _build(cmf, tmp_path)
    hk = b._prop3d('hk_actual')
    vka = b._prop3d('vka_actual')
    k33 = b.gwf.get_package('npf').k33.get_data()
    L, i, j = 0, *b.surf_cells[0][:2]
    if vka[L, i, j] > 0:
        assert np.isclose(k33[L, i, j], hk[L, i, j] / vka[L, i, j], rtol=1e-6)
