# -*- coding: utf-8 -*-
"""GHB and DRN as front-end questions, one sub-panel per MF6 package.

Both are wired: the run builds ModflowGwfghb and ModflowGwfdrn from
[ghb] and [drn], not from the parameter file. The acceptance test is the
same one the layer properties had -- the drain list the panel produces
must equal the one the parameter file produced, entry for entry.
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


def test_the_run_builds_both_packages_from_the_panel():
    """They were unwired for one commit and the tab said so. It says the
    opposite now, so this checks the opposite is true."""
    src = open(os.path.join(CODE, 'tests', 'run_lamata_mf6.py'),
               encoding='utf-8').read()
    assert 'props.apply_boundaries' in src
    # AFTER the properties: a drain taken at the base of its layer reads
    # botm, and botm is the panel thickness now
    assert (src.index('props.apply_layer_properties')
            < src.index('props.apply_boundaries'))
    lib = open(os.path.join(CODE, 'app', 'lib', 'panelui.py'),
               encoding='utf-8').read()
    assert 'Not read by a run yet' not in lib


def test_which_layers_gets_a_real_widget():
    """A bare list falls back to 'edit in the TOML tab', which is not an
    answer to 'which layers does this apply to'."""
    assert 'ghb.layers' in schema.LAYER_LIST
    assert 'drn.layers' in schema.LAYER_LIST
    lib = open(os.path.join(CODE, 'app', 'lib', 'panelui.py'),
               encoding='utf-8').read()
    assert 'def layer_box(' in lib
    assert 'schema.LAYER_LIST' in lib


# ------------------------------------- what the run actually builds from it

DS = os.path.abspath(os.path.join(CODE, '..', 'example', 'LaMata'))
for _p in (os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'MARMITESutilities'),
           os.path.join(CODE, 'ppMF_FloPy')):
    if _p not in sys.path:
        sys.path.insert(0, _p)
props = _load('marmites_props_b', os.path.join(CODE, 'ppMF6',
                                               'marmites_props.py'))


def _built(cfg):
    """A clsMF parsed from the parameter file, then given the panel's
    answers -- the properties first, because a drain at the base of its
    layer reads botm and botm is the panel thickness now."""
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    cUTIL = MMutils.clsUTILITIES(verbose=0)
    cMF = ppMF.clsMF(cUTIL, MM_ws=DS, MM_ws_out=DS,
                     MF_ws=os.path.join(DS, 'MF_ws'),
                     MF_ini_fn='__inputMF_flopy_v3_2s1L.ini',
                     xllcorner=739300.0, yllcorner=4553050.0)
    from_ini = [list(r) for r in cMF.layer_row_column_elevation_cond[0]]
    props.apply_layer_properties(cfg, cMF, DS, verbose=False)
    props.apply_boundaries(cfg, cMF, DS, verbose=False)
    return cMF, from_ini


@pytest.mark.skipif(not os.path.isdir(os.path.join(DS, 'MF_ws')),
                    reason='the La Mata dataset is not present')
def test_the_drains_reproduce_the_parameter_file_entry_for_entry(cfg):
    """THE acceptance test. Six outlet cells on each of two layers, each at
    the base of its own layer, with the layer's own conductance."""
    cMF, from_ini = _built(cfg)
    from_panel = [list(r) for r in cMF.layer_row_column_elevation_cond[0]]
    assert len(from_panel) == len(from_ini) == 12
    for a, b in zip(sorted(from_ini), sorted(from_panel)):
        assert a[:3] == b[:3]
        assert abs(a[3] - b[3]) < 1e-9, 'drain elevation moved'
        assert abs(a[4] - b[4]) < 1e-9, 'drain conductance changed'


@pytest.mark.skipif(not os.path.isdir(os.path.join(DS, 'MF_ws')),
                    reason='the La Mata dataset is not present')
def test_a_drain_at_the_layer_base_sits_just_above_botm(cfg):
    """botm + 0.01 m: on the bottom face the cell would be dry before the
    drain ever took water."""
    cMF, _ini = _built(cfg)
    for (l, i, j, elev, _cond) in cMF.layer_row_column_elevation_cond[0]:
        assert abs(elev - (cMF.botm[l, i, j] + 0.01)) < 1e-9


@pytest.mark.skipif(not os.path.isdir(os.path.join(DS, 'MF_ws')),
                    reason='the La Mata dataset is not present')
def test_the_ghb_is_built_where_the_head_raster_has_a_head(cfg):
    """La Mata's GHB is off, so there is no list from the parameter file to
    compare with: this checks against the rasters themselves -- six cells a
    layer, heads between 734 and 734.5, each layer's own conductance."""
    cfg.ghb.enable = True
    cMF, _ini = _built(cfg)
    spd = cMF.layer_row_column_head_cond[0]
    assert len(spd) == 12, 'six boundary cells on each of two layers'
    assert {int(r[0]) for r in spd} == {0, 1}
    for (_l, _i, _j, head, cond) in spd:
        assert 734.0 <= head <= 734.5
        assert cond > 0
    # each layer keeps its OWN conductance: 0.0101 and 0.7501
    per_layer = {}
    for (l, _i, _j, _h, cond) in spd:
        per_layer.setdefault(int(l), set()).add(round(float(cond), 4))
    assert per_layer[0] == {0.0101} and per_layer[1] == {0.7501}


@pytest.mark.skipif(not os.path.isdir(os.path.join(DS, 'MF_ws')),
                    reason='the La Mata dataset is not present')
def test_a_layer_left_out_carries_no_boundary(cfg):
    """`layers` is the answer to which layers, and leaving one out has to
    MEAN something -- otherwise it is a question with no effect."""
    cfg.drn.layers = [1]
    cMF, _ini = _built(cfg)
    got = cMF.layer_row_column_elevation_cond[0]
    assert len(got) == 6 and {int(r[0]) for r in got} == {0}


@pytest.mark.skipif(not os.path.isdir(os.path.join(DS, 'MF_ws')),
                    reason='the La Mata dataset is not present')
def test_a_package_switched_off_leaves_nothing_behind(cfg):
    """Off has to reach the run as ghb_yn/drn_yn = 0, or the package would
    be built from whatever the parameter file left in place."""
    cfg.drn.enable = False
    cMF, _ini = _built(cfg)
    assert cMF.drn_yn == 0
    assert cMF.ghb_yn == 0
