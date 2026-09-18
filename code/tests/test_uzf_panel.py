# -*- coding: utf-8 -*-
"""UZF as its own sub-panel, and where every UZF1 name went.

The parameter file carried twenty-two UZF names and MODFLOW 6 keeps
eight. The eight are asked; the rest are absent WITH A REASON, which is
recorded here so that "it is gone" cannot quietly become "we forgot it".
"""

import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
REF = os.path.join(CODE, 'configs', 'lamata.toml')
DS = os.path.abspath(os.path.join(CODE, '..', 'example', 'LaMata'))
PAGE = os.path.join(CODE, 'app', 'pages',
                    '4_Unsaturated_zone_and_groundwater.py')

for _p in (CODE, os.path.join(CODE, 'app'), os.path.join(CODE, 'ppMF6'),
           os.path.join(CODE, 'MARMITESutilities'),
           os.path.join(CODE, 'ppMF_FloPy')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


cfgmod = _load('marmites_config_u', os.path.join(CODE, 'marmites_config.py'))
props = _load('marmites_props_u', os.path.join(CODE, 'ppMF6',
                                               'marmites_props.py'))
from lib import schema                                        # noqa: E402


@pytest.fixture
def cfg():
    return cfgmod.load_run_config(REF)


# ------------------------------------------------------------ what is asked

def test_every_uzf_field_names_its_flopy_argument():
    for dotted in ('uzf.ntrailwaves', 'uzf.nwavesets', 'uzf.surfdep',
                   'uzf.eps', 'uzf.thtr', 'uzf.thts', 'uzf.thti', 'uzf.vks'):
        help_ = schema.describe(dotted)[2]
        assert 'flopy' in help_.lower() or 'packagedata' in help_, (
            '%s does not say what it becomes' % dotted)


def test_vks_scale_says_it_is_not_a_modflow_field():
    """It multiplies vks before writing, to offset a clamp. Presenting it
    beside eight real MF6 fields without saying so would be misleading."""
    help_ = schema.describe('uzf.vks_scale')[2]
    assert 'NOT a MODFLOW field' in help_


def test_iuzfopt_is_asked_as_what_it_means(cfg):
    """1 and 2 meant 'read a raster' and 'use the layer k33'. The numbers
    were never the question."""
    assert cfg.uzf.vks_from in ('layer', 'raster')
    assert schema.choices_for('uzf.vks_from') == ['layer', 'raster']


def test_the_legacy_names_are_all_accounted_for():
    """Every UZF1 name the modeller asked about appears in the table with
    an answer -- kept, or gone with the reason."""
    asked = ('SPECIFYTHTR', 'SPECIFYTHTI', 'NOSURFLEAK', 'nuztop', 'iuzfopt',
             'irunflg', 'ietflg', 'iuzfcb1', 'iuzfcb2', 'ntrail2', 'nsets',
             'nuzgag', 'surfdep', 'uzf_iuzfbnd', 'vks', 'eps', 'thts',
             'thti', 'iuzrow', 'iuzcol', 'iftunit', 'iuzopt', 'finf_user')
    table = ' | '.join('%s -> %s' % (a, b) for a, b in schema.UZF_LEGACY)
    for name in asked:
        assert name in table, '%s is not accounted for' % name


def test_the_things_modflow_6_dropped_are_not_config_fields():
    """A field for one of these would be a question with nowhere to go."""
    import dataclasses
    fields = {f.name for f in dataclasses.fields(cfgmod.Uzf)}
    for dead in ('nuztop', 'irunflg', 'ietflg', 'iuzfcb1', 'iuzfcb2',
                 'nuzgag', 'iuzrow', 'iuzcol', 'iftunit', 'iuzopt',
                 'specifythtr', 'specifythti', 'nosurfleak', 'iuzfbnd'):
        assert dead not in fields, '%s reached the panel' % dead


# ----------------------------------------------------- groundwater ET is out

def test_groundwater_et_in_modflow_is_gone(cfg):
    """MARMITES computes ETg and applies it through WEL. The old switch was
    forced false AND read by nothing, which is worse than absent."""
    import dataclasses
    fields = {f.name for f in dataclasses.fields(cfgmod.Et)}
    assert 'gwet_in_mf' not in fields
    assert 'et.gwet_in_mf' not in schema.FIELDS
    page = open(PAGE, encoding='utf-8').read()
    assert 'Groundwater ET is not asked' in page


def test_the_extinction_depth_follows_the_usual_rule(cfg):
    """A raster, a column of the vegetation layer (as CdL does it), or one
    value -- not a three-way enum that has to be kept in step with it."""
    assert schema.is_source(cfg.et.extdp)
    for producer in ('raster', 'layer', 'value'):
        assert hasattr(cfg.et.extdp, producer)
    assert not hasattr(cfg.et, 'extdp_source')
    assert not hasattr(cfg.et, 'extdp_default')
    help_ = schema.describe('et.extdp')[2]
    assert 'vegetation layer' in help_.lower()


def test_uzf_et_needs_an_extinction_depth(cfg):
    cfg.et.uzf_et = True
    cfg.et.extdp = cfgmod.VectorSource()
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'et.extdp' in str(e.value)


def test_the_panel_admits_uzf_et_is_not_wired(cfg):
    """simulate_et is hard-coded False in the build. Switching it on is a
    modelling decision -- UZF drying the unsaturated zone alongside MMsoil
    -- so the panel says so instead of implying the switch works."""
    page = open(PAGE, encoding='utf-8').read()
    assert 'simulate_et' in page and 'hard-coded' in page
    build = open(os.path.join(CODE, 'ppMF6', 'marmites_mf6.py'),
                 encoding='utf-8').read()
    assert 'simulate_et=False' in build, (
        'if simulate_et now follows the panel, this test and the panel '
        'warning should both be replaced')


# ------------------------------------------------- what the run receives

def test_the_uzf6_limits_are_refused_here_not_by_modflow(cfg):
    """Better a message on the panel than a failure after the run has
    started writing files."""
    for field, bad, needle in (('eps', 2.0, '3.5'),
                               ('thtr', 0.0, 'thtr'),
                               ('thti', 0.9, 'thti'),
                               ('ntrailwaves', 0, 'ntrailwaves')):
        c = cfgmod.load_run_config(REF)
        setattr(c.uzf, field, bad)
        with pytest.raises(cfgmod.ConfigError) as e:
            c.validate()
        assert needle in str(e.value), field


@pytest.mark.skipif(not os.path.isdir(os.path.join(DS, 'MF_ws')),
                    reason='the La Mata dataset is not present')
def test_the_panel_gives_uzf_what_the_parameter_file_gave_it(cfg):
    """Everything but eps is identical. eps differs BECAUSE the build
    clamps it: the file said 2.0, MODFLOW 6 forbids below 3.5, so the model
    has always run at 3.5 -- the panel says 3.5 instead of being corrected
    silently, and the built model is the same."""
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF

    cUTIL = MMutils.clsUTILITIES(verbose=0)
    cMF = ppMF.clsMF(cUTIL, MM_ws=DS, MM_ws_out=DS,
                     MF_ws=os.path.join(DS, 'MF_ws'),
                     MF_ini_fn='__inputMF_flopy_v3_2s1L.ini',
                     xllcorner=739300.0, yllcorner=4553050.0)

    def scalar(v):
        return float(np.ravel(np.asarray(v, dtype=float))[0])

    names = ('ntrail2', 'nsets', 'surfdep', 'thtr', 'thts', 'thti', 'iuzfopt')
    before = dict((n, scalar(getattr(cMF, n))) for n in names)
    eps_ini = scalar(cMF.eps)
    props.apply_uzf(cfg, cMF, DS, verbose=False)
    for n in names:
        assert abs(before[n] - scalar(getattr(cMF, n))) < 1e-12, n
    assert eps_ini == 2.0 and scalar(cMF.eps) == 3.5
    # ... and both land on 3.5 once the build has applied its clamp
    assert max(3.5, min(14.0, eps_ini)) == scalar(cMF.eps)


def test_the_run_applies_the_unsaturated_zone():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert 'props.apply_uzf' in src
