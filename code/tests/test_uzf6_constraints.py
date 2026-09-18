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
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s1L.ini')):
        pytest.skip('La Mata dataset not present')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s1L.ini',
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


def test_exactly_one_seepage_mechanism_exists(cmf, tmp_path):
    """Exfiltration (MARMITES Exf_g) must be structurally POSSIBLE -- the
    first full La Mata run produced none because neither mechanism was
    active -- and must not be provided twice, or the discharge is counted
    twice.

    Which mechanism is a separate question with its own default. Under
    'uzf' it is UZF6's SIMULATE_GWSEEP; under 'drn' (the default since a
    new catchment must not inherit the deprecated option) it is the
    drn_seep package, whose SIMVALS the coupler reads back into the soil
    column. So this pins the INVARIANT, in both settings.
    """
    for kind in ('uzf', 'drn'):
        ws = tmp_path / kind
        ws.mkdir()
        b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                          botm=np.asarray(cmf.botm, dtype=float),
                          sim_ws=str(ws))
        b.seep = kind
        assert b.gwseep is True
        b.build()
        b.write()
        with open(os.path.join(str(ws),
                               '%s.uzf' % cmf.modelname.lower())) as f:
            in_uzf = 'SIMULATE_GWSEEP' in f.read(4000).upper()
        by_drain = b.ndrnseep > 0
        assert in_uzf != by_drain, (
            'seep=%r gives %s in UZF and %s by drain: exfiltration is '
            'either impossible or counted twice'
            % (kind, in_uzf, by_drain))
        assert in_uzf if kind == 'uzf' else by_drain


def test_the_default_mechanism_is_the_drain(cmf, tmp_path):
    """Not SIMULATE_GWSEEP, which MODFLOW 6 deprecates and which switches
    discharge on and off discontinuously."""
    b = _build(cmf, tmp_path)
    assert b.seep == 'drn'
    b.write()
    with open(os.path.join(str(tmp_path),
                           '%s.uzf' % cmf.modelname.lower())) as f:
        assert 'SIMULATE_GWSEEP' not in f.read(4000).upper()
    assert b.ndrnseep > 0, 'no seepage face at all: Exf_g would be zero'


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


# --------------------------------------------------------------------- #
# the ET split: ETsoil (MM) + ETuzf (MF) + ETg (MM)
# --------------------------------------------------------------------- #

def _uzf_text(b, tmp):
    b.write()
    with open(os.path.join(str(tmp), '%s.uzf' % b.cMF.modelname.lower())) as f:
        return f.read().upper()


def test_uzf_et_is_off_until_it_is_asked_for(cmf, tmp_path):
    b = _build(cmf, tmp_path)
    assert 'SIMULATE_ET' not in _uzf_text(b, tmp_path)


def test_uzf_does_the_unsaturated_zone_and_never_the_groundwater(cmf,
                                                                 tmp_path):
    """THE ruling: total ET is ETsoil + ETuzf + ETg -- the soil column and
    the groundwater are MARMITES's, the unsaturated zone between them is
    UZF's. MODFLOW 6 supports exactly that: "et can be simulated in the uzf
    cell and not the gwf cell by omitting keywords linear_gwet and
    square_gwet". Asking for groundwater ET here would remove the same
    water MARMITES already removes through WEL."""
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp_path))
    b.uzf_et = True
    b.uzf_extdp = np.full((cmf.nlay, cmf.nrow, cmf.ncol), 2.5)
    b.build()
    txt = _uzf_text(b, tmp_path)
    assert 'SIMULATE_ET' in txt
    assert 'UNSAT_ETWC' in txt
    assert 'LINEAR_GWET' not in txt, 'MODFLOW would remove ETg a second time'
    assert 'SQUARE_GWET' not in txt, 'MODFLOW would remove ETg a second time'


def test_the_unsaturated_formulation_follows_the_configuration(cmf, tmp_path):
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp_path))
    b.uzf_et, b.uzf_et_form = True, 'etae'
    b.uzf_extdp = np.full((cmf.nlay, cmf.nrow, cmf.ncol), 2.5)
    b.build()
    txt = _uzf_text(b, tmp_path)
    assert 'UNSAT_ETAE' in txt and 'UNSAT_ETWC' not in txt


def test_the_pet_demand_starts_at_zero_for_the_coupler_to_write(cmf,
                                                                tmp_path):
    """PET is a daily quantity MARMITES computes, not a property of the
    model: the build leaves it at 0 and the coupler writes it each step.
    A non-zero constant here would be a demand nobody chose, applied every
    day of the run."""
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp_path))
    b.uzf_et = True
    b.uzf_extdp = np.full((cmf.nlay, cmf.nrow, cmf.ncol), 2.5)
    b.build()
    pet = {float(rec[2]) for rec in b.uzf_perioddata}
    assert pet == {0.0}, 'the build wrote a PET demand: %s' % sorted(pet)
    extdp = {round(float(rec[3]), 3) for rec in b.uzf_perioddata}
    assert extdp == {2.5}, extdp
