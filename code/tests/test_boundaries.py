# -*- coding: utf-8 -*-
"""GHB and DRN as front-end questions, one sub-panel per MF6 package.

These tabs ASK. They are not wired to the run yet, and the panel says so
rather than implying otherwise -- a test below pins that admission, so it
cannot quietly disappear before the wiring lands.
"""

import importlib.util
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
REF = os.path.join(CODE, 'configs', 'lamata.toml')
PAGE = os.path.join(CODE, 'app', 'pages',
                    '4_Unsaturated_zone_and_groundwater.py')

for _p in (CODE, os.path.join(CODE, 'app')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


cfgmod = _load('marmites_config_b', os.path.join(CODE, 'marmites_config.py'))
from lib import schema                                        # noqa: E402


@pytest.fixture
def cfg():
    return cfgmod.load_run_config(REF)


def test_both_packages_are_sections_of_their_own(cfg):
    """One sub-panel per MF6 package means one SECTION per package: a
    boundary buried in [layers] would be a boundary nobody finds."""
    assert 'ghb' in cfgmod._SECTIONS and 'drn' in cfgmod._SECTIONS
    assert cfg.ghb.enable is False        # La Mata runs without a GHB
    assert cfg.drn.enable is True         # and with six outlet drain cells


def test_every_field_names_its_flopy_argument():
    for dotted in ('ghb.enable', 'ghb.head', 'ghb.cond',
                   'drn.enable', 'drn.elevation', 'drn.cond'):
        help_ = schema.describe(dotted)[2]
        assert 'Modflow' in help_, '%s does not say what it becomes' % dotted


def test_the_layers_are_counted_from_one(cfg):
    """MODFLOW, the parameter file and every conversation about the model
    count layers from 1. A panel that counted from 0 would be the only
    thing in the room that did."""
    assert cfg.drn.layers == [1, 2]
    for bad, needle in (([1, 3], 'layer 3'), ([], 'empty'), ([1, 1], 'repeat')):
        cfg.drn.layers = bad
        with pytest.raises(cfgmod.ConfigError) as e:
            cfg.validate()
        assert needle in str(e.value)


def test_a_package_that_is_off_is_not_interrogated(cfg):
    """Half-filled answers for something that is not built are not errors:
    La Mata's GHB is off and its head/cond are still named."""
    assert cfg.ghb.enable is False
    cfg.ghb.head = cfgmod.VectorSource()
    cfg.ghb.cond = cfgmod.VectorSource()
    cfg.ghb.layers = []
    cfg.validate()          # must not raise


def test_a_package_that_is_on_must_be_answerable(cfg):
    cfg.ghb.enable = True
    cfg.ghb.head = cfgmod.VectorSource()
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'ghb.head' in str(e.value)


def test_a_conductance_of_zero_is_refused(cfg):
    """A boundary with no conductance is a boundary that does nothing,
    which is a way of spelling 'off' that the switch already spells."""
    cfg.drn.cond = cfgmod.VectorSource(value=0.0)
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'cond' in str(e.value)


def test_the_layer_base_is_a_switch_not_a_magic_number(cfg):
    """The parameter file wrote -1 in the raster and explained it in a
    comment. Every La Mata drain cell is -1, so this is the normal case and
    it deserves a name."""
    assert cfg.drn.at_layer_base is True
    help_ = schema.describe('drn.at_layer_base')[2]
    assert 'botm' in help_ and '-1' in help_


def test_the_seepage_face_is_not_confused_with_this_drain():
    """MARMITES builds TWO drain packages -- drn here, drn_seep from
    [seep] -- and a modeller who thought this tab configured the seepage
    face would turn the catchment outlet off looking for it."""
    page = open(PAGE, encoding='utf-8').read()
    assert 'drn_seep' in page and 'not the seepage face' in page.lower()
    help_ = schema.describe('drn.enable')[2]
    assert 'drn_seep' in help_


def test_the_panel_admits_these_are_not_wired_yet():
    """The layers tab earned the right to say 'every field reaches the
    run'. These two have not, and must not imply it."""
    lib = open(os.path.join(CODE, 'app', 'lib', 'panelui.py'),
               encoding='utf-8').read()
    assert 'Not read by a run yet' in lib, (
        'boundary_note no longer admits the packages are unwired -- if the '
        'wiring landed, this test should be replaced by one that checks the '
        'run reads [ghb] and [drn]')


def test_which_layers_gets_a_real_widget():
    """A bare list falls back to 'edit in the TOML tab', which is not an
    answer to 'which layers does this apply to'."""
    assert 'ghb.layers' in schema.LAYER_LIST
    assert 'drn.layers' in schema.LAYER_LIST
    lib = open(os.path.join(CODE, 'app', 'lib', 'panelui.py'),
               encoding='utf-8').read()
    assert 'def layer_box(' in lib
    assert 'schema.LAYER_LIST' in lib
