# -*- coding: utf-8 -*-
"""The layer properties come from the TOML, not from the MODFLOW ini.

The acceptance test is not "it runs": it is that the panel's answers
reproduce what ``__inputMF_flopy_v3_2s1L.ini`` produced, cell for cell.
A wiring that changed the model while claiming to move an input would be
worse than no wiring at all.
"""

import dataclasses
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
REF = os.path.join(CODE, 'configs', 'lamata.toml')
DS = os.path.abspath(os.path.join(CODE, '..', 'example', 'LaMata'))

for _p in (CODE, os.path.join(CODE, 'ppMF6'),
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


props = _load('marmites_props_t', os.path.join(CODE, 'ppMF6',
                                               'marmites_props.py'))
cfgmod = _load('marmites_config_lp', os.path.join(CODE, 'marmites_config.py'))


@pytest.fixture
def cfg():
    return cfgmod.load_run_config(REF)


# ------------------------------------------------------- the resolver alone

def test_a_value_is_uniform_over_every_layer():
    src = cfgmod.VectorSource(value=2.5)
    assert props.resolve_source(src, 3, DS, 'layers.k') == [2.5, 2.5, 2.5]


def test_an_unanswered_source_leaves_the_ini_in_charge():
    """None means 'not answered', which is NOT the same as an error: the
    parameter file still supplies it until the panel says otherwise."""
    assert props.resolve_source(cfgmod.VectorSource(), 2, DS, 'x') is None


def test_a_pattern_becomes_one_absolute_path_per_layer():
    src = cfgmod.VectorSource(raster='MF_ws/thick_l%d.asc')
    got = props.resolve_source(src, 2, DS, 'layers.thickness')
    assert [os.path.basename(p) for p in got] == ['thick_l1.asc',
                                                  'thick_l2.asc']
    assert all(os.path.isabs(p) for p in got), (
        'checkarray joins with MF_ws, so a relative name would be resolved '
        'against the wrong folder')


def test_one_raster_without_the_placeholder_serves_every_layer():
    src = cfgmod.VectorSource(raster='MF_ws/thick_l1.asc')
    got = props.resolve_source(src, 2, DS, 'layers.thickness')
    assert len(got) == 2 and got[0] == got[1]


def test_a_missing_raster_is_refused_by_name():
    src = cfgmod.VectorSource(raster='MF_ws/not_here_l%d.asc')
    with pytest.raises(props.PropertyError) as e:
        props.resolve_source(src, 2, DS, 'layers.k')
    assert 'not_here_l1.asc' in str(e.value)


def test_the_wrong_case_is_refused_even_though_windows_would_open_it(tmp_path):
    """Ss_l2.asc and ss_l2.asc are the same file here and two different
    files on Linux. A run that works on one machine only is a bug."""
    folder = tmp_path / 'MF_ws'
    folder.mkdir()
    (folder / 'k_l1.asc').write_text('x')
    (folder / 'k_l2.asc').write_text('x')
    src = cfgmod.VectorSource(raster='MF_ws/K_l%d.asc')
    with pytest.raises(props.PropertyError) as e:
        props.resolve_source(src, 2, str(tmp_path), 'layers.k')
    assert 'spelled' in str(e.value) and 'k_l1.asc' in str(e.value)


def test_a_shapefile_is_refused_rather_than_guessed_at():
    """Wrapping a polygon attribute onto the layers is the converter's job
    and it does not do it yet. Falling back to the ini silently is the
    confusion this module exists to end."""
    src = cfgmod.VectorSource(layer='geology.shp', column='K')
    with pytest.raises(props.PropertyError) as e:
        props.resolve_source(src, 2, DS, 'layers.k')
    assert 'not implemented' in str(e.value)


# --------------------------------------------- against the real parameter file

@pytest.mark.skipif(not os.path.isdir(os.path.join(DS, 'MF_ws')),
                    reason='the La Mata dataset is not present')
def test_the_panel_reproduces_the_ini_cell_for_cell(cfg):
    """THE acceptance test for this wiring."""
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF

    cUTIL = MMutils.clsUTILITIES(verbose=0)
    cMF = ppMF.clsMF(cUTIL, MM_ws=DS, MM_ws_out=DS,
                     MF_ws=os.path.join(DS, 'MF_ws'),
                     MF_ini_fn='__inputMF_flopy_v3_2s1L.ini',
                     xllcorner=739300.0, yllcorner=4553050.0)
    f2a = cMF.cPROCESS.float2array

    def snapshot():
        out = dict((n, np.array(f2a(getattr(cMF, n + '_actual')), dtype=float))
                   for n in ('hk', 'vka', 'ss', 'sy'))
        out['thick'] = np.array(cMF.thick, dtype=float)
        out['botm'] = np.ma.filled(np.asarray(cMF.botm), -9999.0).astype(float)
        return out

    from_ini = snapshot()
    assert cfg.layers.nlay == 2, 'this test is pinned to the 2-layer set'
    done = props.apply_layer_properties(cfg, cMF, DS, verbose=False)
    assert {d[0] for d in done} == {'ibound', 'thickness', 'k', 'k33',
                                    'ss', 'sy'}
    from_panel = snapshot()
    for name in ('thick', 'hk', 'vka', 'ss', 'sy', 'botm'):
        assert np.array_equal(from_ini[name], from_panel[name]), (
            '%s changed when it moved from the ini to the panel' % name)


def test_the_flags_are_fanned_out_over_the_layers(cfg):
    """The ini gave one integer per layer and every parameter set here
    repeats the same one, so the panel asks once."""
    class FakeProcess(object):
        def checkarray(self, v, **kw):
            return v

        def float2array(self, v):
            return v

    class FakeMF(object):
        nlay = 4
        cPROCESS = FakeProcess()

    cMF = FakeMF()
    # only the flags here: the sources are checked against the real dataset
    # above, and this model has four layers the rasters do not cover
    for name in ('ibound', 'thickness', 'k', 'k33', 'ss', 'sy'):
        setattr(cfg.layers, name, cfgmod.VectorSource())
    cfg.layers.convertible = True
    cfg.layers.k33_as_ratio = False
    props.apply_layer_properties(cfg, cMF, DS, verbose=False)
    assert list(cMF.laytyp) == [1, 1, 1, 1]
    assert list(cMF.layvka) == [0, 0, 0, 0]


def test_k33_as_a_ratio_must_be_positive(cfg):
    """It divides k on the way to k33, and below 1 it would say the aquifer
    conducts water more easily downwards than sideways."""
    cfg.layers.k33 = cfgmod.VectorSource(value=0.0)
    cfg.layers.k33_as_ratio = True
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'k33' in str(e.value)


def test_what_modflow_6_has_not_is_not_asked():
    """layavg, laywet, laycbd and hdry are MODFLOW-2005/NWT only, and the
    run sets hdry to None outright. A panel field for any of them would be
    a question with nowhere to go."""
    fields = {f.name for f in dataclasses.fields(cfgmod.Layers)}
    for dead in ('layavg', 'laywet', 'laycbd', 'hdry'):
        assert dead not in fields, '%s reached the panel' % dead
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert 'cMF.hdry = None' in src


def test_the_reference_config_names_a_pattern_not_one_layer(cfg):
    """`thick_l1.asc` for a 2-layer model gives layer 2 layer 1's numbers.
    For Ss that is 9.5e-05 where the confining unit wants 4e-07."""
    for name in ('thickness', 'k', 'ss', 'sy'):
        raster = getattr(cfg.layers, name).raster
        assert '%d' in raster, (
            'layers.%s = %r names a single layer; with nlay > 1 every layer '
            'would get it' % (name, raster))


def test_the_run_applies_the_properties_before_anything_reads_them():
    """A panel field the run ignores is decoration; one applied too late is
    worse, because half the model would have the old value."""
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert 'props.apply_layer_properties' in src
    assert src.index('props.apply_layer_properties') < src.index('conv_fact = ')


# ------------------------------------- ibound, and the geographic reference

@pytest.mark.skipif(not os.path.isdir(os.path.join(DS, 'MF_ws')),
                    reason='the La Mata dataset is not present')
def test_ibound_is_per_layer_and_reproduces_the_parameter_file(cfg):
    """A layer can pinch out INSIDE the catchment: La Mata's layer 1 is
    active in 1870 cells and layer 2 in 1954, layer 1 being a strict
    subset. The catchment polygon cannot know that -- thick_l1 is 20-35 m
    in the 84 cells where layer 1 is absent, so the thickness cannot
    express it either. Only this map can."""
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF

    rect, _n, _o = props.dataset_grid(DS)
    cMF = ppMF.clsMF(MMutils.clsUTILITIES(verbose=0), MM_ws=DS, MM_ws_out=DS,
                     MF_ws=os.path.join(DS, 'MF_ws'),
                     MF_ini_fn='__inputMF_flopy_v3_2s1L.ini', grid=rect)
    from_ini = np.abs(np.asarray(cMF.ibound, dtype=int))
    botm_ini = np.ma.filled(np.asarray(cMF.botm), -9999.0).astype(float)

    props.apply_layer_properties(cfg, cMF, DS, verbose=False)
    from_panel = np.abs(np.asarray(cMF.ibound, dtype=int))

    assert from_panel.shape == (cMF.nlay, cMF.nrow, cMF.ncol)
    assert np.array_equal(from_ini, from_panel)
    per_layer = [int((from_panel[L] != 0).sum()) for L in range(cMF.nlay)]
    assert per_layer == [1870, 1954], per_layer
    # the layers genuinely differ, which is the whole reason this is asked
    assert per_layer[0] != per_layer[1]
    # ... and botm is unchanged, though it multiplies by |ibound|
    assert np.array_equal(
        botm_ini, np.ma.filled(np.asarray(cMF.botm), -9999.0).astype(float))


@pytest.mark.skipif(not os.path.isdir(os.path.join(DS, 'MF_ws')),
                    reason='the La Mata dataset is not present')
def test_the_catchment_is_the_geographic_reference_not_the_ibound(cfg):
    """It checks that the active cells sit inside it; it does not decide
    which they are. The two are allowed to differ -- the polygon touches
    2092 cells and 1954 are active -- because a cell clipped at the
    boundary is a modelling choice."""
    import mm_paths
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF

    rect, _n, _o = props.dataset_grid(DS)
    cMF = ppMF.clsMF(MMutils.clsUTILITIES(verbose=0), MM_ws=DS, MM_ws_out=DS,
                     MF_ws=os.path.join(DS, 'MF_ws'),
                     MF_ini_fn='__inputMF_flopy_v3_2s1L.ini', grid=rect)
    rep = props.check_catchment(cfg, cMF, mm_paths.GIS, verbose=False)
    assert rep['outside'] == 0, 'active cells fall outside the catchment'
    assert rep['fraction_inside'] == 1.0
    assert rep['catchment'] > rep['active'], (
        'the polygon should touch more cells than the model activates')


@pytest.mark.skipif(not os.path.isdir(os.path.join(DS, 'MF_ws')),
                    reason='the La Mata dataset is not present')
def test_a_model_somewhere_else_entirely_is_refused(cfg):
    """The case that is never intentional: rasters and catchment in
    different coordinate systems. An edge disagreement is fine; a model
    mostly outside its own catchment is not."""
    import mm_paths
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF

    rect, _n, _o = props.dataset_grid(DS)
    moved = (rect[0] + 50000.0, rect[1] + 50000.0) + tuple(rect[2:])
    cMF = ppMF.clsMF(MMutils.clsUTILITIES(verbose=0), MM_ws=DS, MM_ws_out=DS,
                     MF_ws=os.path.join(DS, 'MF_ws'),
                     MF_ini_fn='__inputMF_flopy_v3_2s1L.ini', grid=moved)
    with pytest.raises(props.PropertyError) as e:
        props.check_catchment(cfg, cMF, mm_paths.GIS, verbose=False)
    assert 'coordinate system' in str(e.value)
