# -*- coding: utf-8 -*-
"""Phase-4 tests: UZF6 input constraints.

The first real MODFLOW 6 run of La Mata failed with ~22,000 UZF validation
errors of two kinds -- THTR = 0 and EPSILON out of range. Both are genuine
NWT->MF6 semantic differences that the build stage can and must catch, so
they are asserted here rather than discovered by the solver.

MODFLOW 6 UZF rules covered:
    THTR > 0            (UZF1 tolerated 0 when specifythtr = 0)
    THTS > THTR
    THTR <= THTI <= THTS
    3.5 <= EPSILON <= 14.0   (UZF1 accepted any value)
    VKS > 0
"""
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

pytest.importorskip('flopy')
import matplotlib  # noqa: E402
matplotlib.use('agg')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mf6mod = _load('marmites_mf6_uzf', os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'))


@pytest.fixture(scope='module')
def cmf():
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s3L.ini')):
        pytest.skip('La Mata dataset not present')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s3L.ini',
                   xllcorner=739300.0, yllcorner=4553050.0)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 2, [1, 1], [1, 1]
    return c


def _build(cmf, tmp, grid='dis'):
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp), grid=grid)
    b.build()
    return b


# --------------------------------------------------------------------- #
# the real La Mata package must satisfy every UZF6 rule
# --------------------------------------------------------------------- #

def test_lamata_uzf_packagedata_satisfies_mf6_rules(cmf, tmp_path):
    b = _build(cmf, tmp_path)
    # packagedata: (iuzno, cellid, landflag, ivertcon, surfdep, vks,
    #               thtr, thts, thti, eps)
    for rec in b.uzf_packagedata:
        surfdep, vks, thtr, thts, thti, eps = (float(rec[4]), float(rec[5]),
                                               float(rec[6]), float(rec[7]),
                                               float(rec[8]), float(rec[9]))
        assert thtr > 0.0, 'THTR must be > 0 (UZF6)'
        assert thts > thtr, 'THTS must exceed THTR'
        assert thtr <= thti <= thts, 'THTI must lie in [THTR, THTS]'
        assert mf6mod.clsMF6.EPS_MIN <= eps <= mf6mod.clsMF6.EPS_MAX, \
            'EPSILON must be within the UZF6 range'
        assert vks > 0.0, 'VKS must be > 0'
        assert surfdep > 0.0, 'SURFDEP must be > 0'


def test_epsilon_clamped_and_reported(cmf, tmp_path):
    """La Mata's ini has EPSILON = 2.0 (valid in UZF1, invalid in UZF6):
    the builder must clamp it AND record that it did."""
    b = _build(cmf, tmp_path)
    assert b.eps_clamped is not None, 'clamping must be recorded, not silent'
    original, clamped = b.eps_clamped
    assert original == 2.0
    assert clamped == mf6mod.clsMF6.EPS_MIN == 3.5


def test_thtr_taken_from_ini_even_when_specifythtr_zero(cmf, tmp_path):
    """specifythtr=0 in the ini, but a THTR value is supplied on the same
    line; UZF6 needs it, so it must be used rather than defaulted to 0."""
    assert int(cmf.specifythtr) == 0
    ini_thtr = float(np.ravel(np.asarray(cmf.thtr, dtype=float))[0])
    assert ini_thtr > 0
    b = _build(cmf, tmp_path)
    assert np.isclose(float(b.uzf_packagedata[0][6]), ini_thtr)


# --------------------------------------------------------------------- #
# the validator itself
# --------------------------------------------------------------------- #

def _validator(cmf, tmp_path):
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float), sim_ws=str(tmp_path))
    return b


def test_validator_rejects_zero_thtr(cmf, tmp_path):
    b = _validator(cmf, tmp_path)
    with pytest.raises(mf6mod.MF6BuildError, match='THTR'):
        b._validate_uzf_params(0.0, 0.45, 0.15, 4.0)


def test_validator_rejects_thts_below_thtr(cmf, tmp_path):
    b = _validator(cmf, tmp_path)
    with pytest.raises(mf6mod.MF6BuildError, match='THTS'):
        b._validate_uzf_params(0.5, 0.45, 0.45, 4.0)


def test_validator_rejects_thti_out_of_range(cmf, tmp_path):
    b = _validator(cmf, tmp_path)
    with pytest.raises(mf6mod.MF6BuildError, match='THTI'):
        b._validate_uzf_params(0.05, 0.45, 0.9, 4.0)


def test_validator_passes_valid_values_untouched(cmf, tmp_path):
    b = _validator(cmf, tmp_path)
    out = b._validate_uzf_params(0.05, 0.45, 0.15, 4.0)
    assert out == (0.05, 0.45, 0.15, 4.0)
    assert b.eps_clamped is None


def test_gwseep_enabled_so_exfiltration_can_exist(cmf, tmp_path):
    """UZF6 computes groundwater discharge (= MARMITES Exf_g) ONLY with
    SIMULATE_GWSEEP. Without it the GWD array stays identically zero and
    exfiltration is structurally impossible -- which is exactly what happened
    in the first full La Mata run. Guard the option in the written file."""
    b = _build(cmf, tmp_path)
    assert b.gwseep is True
    b.write()
    uzf = os.path.join(str(tmp_path), '%s.uzf' % cmf.modelname.lower())
    with open(uzf) as f:
        head = f.read(4000).upper()
    assert 'SIMULATE_GWSEEP' in head, 'SIMULATE_GWSEEP missing -> Exf_g == 0'


def test_gwseep_can_be_disabled_explicitly(cmf, tmp_path):
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp_path), gwseep=False)
    b.build(); b.write()
    uzf = os.path.join(str(tmp_path), '%s.uzf' % cmf.modelname.lower())
    with open(uzf) as f:
        head = f.read(4000).upper()
    assert 'SIMULATE_GWSEEP' not in head


def test_validator_clamps_high_epsilon(cmf, tmp_path):
    b = _validator(cmf, tmp_path)
    _, _, _, eps = b._validate_uzf_params(0.05, 0.45, 0.15, 99.0)
    assert eps == mf6mod.clsMF6.EPS_MAX == 14.0
    assert b.eps_clamped == (99.0, 14.0)
