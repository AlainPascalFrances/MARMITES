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
                    '4_4_-_Unsaturated_zone_and_groundwater.py')

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
    """MARMITES computes ETg and the EVT packages take it. The old switch was
    forced false AND read by nothing, which is worse than absent."""
    import dataclasses
    fields = {f.name for f in dataclasses.fields(cfgmod.Et)}
    assert 'gwet_in_mf' not in fields
    assert 'et.gwet_in_mf' not in schema.FIELDS
    page = open(PAGE, encoding='utf-8').read()
    # still computed by MARMITES and never asked of UZF; the panel asks
    # the shape of the EVT curves (which package applied it -- WEL or
    # EVT, et.gw_route -- was asked 2026-10-05 to 2026-10-07)
    assert 'Groundwater ET is computed by MARMITES, not asked of UZF' in page
    assert 'schema.GW_ET_ROWS' in page


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


def test_the_extinction_depth_is_always_needed(cfg):
    """There is no "off" in which a missing extdp would be harmless: UZF
    always simulates unsaturated-zone ET, so the depth it stops at is
    always read."""
    cfg.et.extdp = cfgmod.VectorSource()
    with pytest.raises(cfgmod.ConfigError) as e:
        cfg.validate()
    assert 'et.extdp' in str(e.value)


def test_there_is_no_switch_for_uzf_et(cfg):
    """WP2's ruling: total ET has three sources and the deep unsaturated
    zone is one of them, so an off position could only ever produce a
    model that evaporates nothing from it. A forced-true switch read by
    nothing is worse than absent -- gwet_in_mf and wel_yn went the same
    way."""
    assert not hasattr(cfg.et, 'uzf_et')
    flat = [r for row in schema.UZF_ET_ROWS for r in row]
    assert 'et.uzf_et' not in flat, 'the switch is still on the panel'
    page = open(PAGE, encoding='utf-8').read()
    assert 'et.uzf_et' not in page, 'the panel still draws the switch'


def test_the_panel_says_uzf_et_is_wired_end_to_end(cfg):
    """WP2 closed the demand chain: the coupler writes the demand to PETMAX
    and reads UZF's ACTUAL ET back. "Half wired" was the truth until then
    and would now be the wrong one."""
    page = open(PAGE, encoding='utf-8').read()
    assert 'Half wired' not in page
    assert 'Wired end to end' in page and 'PETMAX' in page
    assert 'linear_gwet' in page and 'square_gwet' in page
    assert 'actual' in page and 'residual' in page
    build = open(os.path.join(CODE, 'ppMF6', 'marmites_mf6.py'),
                 encoding='utf-8').read()
    assert 'simulate_et=False' not in build, (
        'simulate_et is hard-coded off again')


# ------------------------------------------------- what the run receives

def test_the_uzf6_limits_are_refused_here_not_by_modflow(cfg):
    """Better a message on the panel than a failure after the run has
    started writing files."""
    # thtr and the thti range are the modeller's only with thtr_from =
    # source; with 'sy' thtr is derived and thti clipped at the build.
    for field, bad, needle, mode in (('eps', 2.0, '3.5', 'sy'),
                                     ('thtr', 0.0, 'thtr', 'source'),
                                     ('thti', 0.9, 'thti', 'source'),
                                     ('ntrailwaves', 0, 'ntrailwaves', 'sy'),
                                     ('thtr_from', 'raster', 'thtr_from',
                                      'sy')):
        c = cfgmod.load_run_config(REF)
        c.uzf.thtr_from = mode
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
    silently, and the built model is the same. (What the file gave is
    frozen in lamata_model.INI_UZF: a run reads no parameter file.)"""
    import lamata_model as LM
    cMF = LM.lamata_cmf(cfg, boundaries=False, uzf=False)

    def scalar(v):
        return float(np.ravel(np.asarray(v, dtype=float))[0])

    names = ('ntrail2', 'nsets', 'surfdep', 'thtr', 'thts', 'thti', 'iuzfopt')
    eps_ini = LM.INI_UZF['eps']
    props.apply_uzf(cfg, cMF, DS, verbose=False)
    for n in names:
        assert abs(LM.INI_UZF[n] - scalar(getattr(cMF, n))) < 1e-12, n
    assert eps_ini == 2.0 and scalar(cMF.eps) == 3.5
    # ... and both land on 3.5 once the build has applied its clamp
    assert max(3.5, min(14.0, eps_ini)) == scalar(cMF.eps)


def test_the_run_applies_the_unsaturated_zone():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'), encoding='utf-8').read()
    assert 'props.apply_uzf' in src


# ------------------------------- the Brooks-Corey parameters are per cell

@pytest.fixture(scope='module')
def built():
    """The real La Mata UZF package as clsMF6 builds it."""
    pytest.importorskip('flopy')
    import matplotlib
    matplotlib.use('agg')
    import lamata_model
    mf6 = _load('marmites_mf6_up', os.path.join(CODE, 'ppMF6',
                                                'marmites_mf6.py'))

    def make(tmp, **over):
        # La Mata's model description as the run builds it -- no
        # parameter file (lamata_model derives the outcrop layer too)
        c = lamata_model.lamata_cmf()
        c.nper, c.perlen, c.nstp = 2, [1, 1], [1, 1]
        for k, v in over.items():
            setattr(c, k, v)
        b = mf6.clsMF6(c, top=np.asarray(c.elev, dtype=float),
                       botm=np.asarray(c.botm, dtype=float),
                       sim_ws=str(tmp), grid='dis')
        b.build()
        return b
    return make


def test_a_uniform_answer_builds_what_the_scalar_built(built, tmp_path):
    """The parameter file gave one number for the catchment. Per-cell must
    reproduce it exactly when the answer is uniform, or this is a model
    change dressed up as a wiring change."""
    b = built(tmp_path)
    eps = {float(r[9]) for r in b.uzf_packagedata}
    thtr = {float(r[6]) for r in b.uzf_packagedata}
    assert eps == {3.5}, 'the clamped scalar no longer reaches every object'
    # NOT the ini's 0.05: with SPECIFYTHTR 0 UZF1 never read it and built
    # thtr = thts - Sy = 0.45 - 0.01, which is what the build does now.
    assert {round(t, 9) for t in thtr} == {0.44}


def test_a_per_layer_answer_actually_varies(built, tmp_path):
    """`layers` and rasters would be decoration if every UZF object still
    took the same number. One value per LAYER is the cheapest proof that
    the packagedata is indexed rather than broadcast."""
    b0 = built(tmp_path / 'probe')
    per_layer = [3.5 + 0.5 * k for k in range(b0.nlay)]
    b = built(tmp_path, eps=per_layer)
    # cellid -> layer is the first element on a DIS grid
    by_layer = {}
    for rec in b.uzf_packagedata:
        k = rec[1][0]
        by_layer.setdefault(int(k), set()).add(round(float(rec[9]), 6))
    assert len(by_layer) > 1, 'the model has objects in one layer only'
    for k, values in by_layer.items():
        assert values == {per_layer[k]}, (
            'layer %d got %s, expected %g' % (k, values, per_layer[k]))


def test_the_uzf6_rules_are_checked_on_every_cell(built, tmp_path):
    """A rule that only held for the catchment mean would let one bad cell
    reach MODFLOW, which is where the 22,000 validation errors came from."""
    mf6 = sys.modules['marmites_mf6_up']
    b0 = built(tmp_path / 'probe')
    bad = [0.0] + [0.05] * (b0.nlay - 1)      # layer 1 only
    with pytest.raises(mf6.MF6BuildError) as e:
        built(tmp_path, thtr=bad)
    assert 'THTR' in str(e.value)


def test_the_five_uzf_properties_take_the_usual_three_producers(cfg):
    """eps, thtr, thts, thti and vks -- raster, layer or value, like every
    other spatial input on these panels."""
    for name in ('eps', 'thtr', 'thts', 'thti', 'vks'):
        src = getattr(cfg.uzf, name)
        assert schema.is_source(src), 'uzf.%s is not a source' % name
        for producer in ('raster', 'layer', 'value'):
            assert hasattr(src, producer)


# ------------------------------------------- state that does not fit the grid
# The panel said "This run will not start." That was true while the guard
# raised a CONFIG ERROR, and false from the moment resolve_initial_heads began
# falling back to the land surface and saying so -- a panel threatening a
# refusal the run does not make is the same defect as a switch the run does
# not read.

def test_the_panel_does_not_threaten_a_refusal_the_run_never_makes():
    # WHAT THE PANEL SAYS, not what its comments explain: the comment
    # recording why the old wording was wrong necessarily quotes it.
    page = open(PAGE, encoding='utf-8').read()
    spoken = '\n'.join(
        ln for ln in page.split('with tab_init:', 1)[-1].splitlines()
        if not ln.strip().startswith('#'))
    assert 'This run will not start' not in spoken, (
        'the panel still claims the run refuses; resolve_initial_heads '
        'starts from the land surface instead')
    assert 'will not be used' in spoken, 'it no longer says what happens'
    assert 'land surface' in spoken


def test_the_panel_says_how_to_regenerate_the_spin_up():
    """The message asks for a spin-up on this mesh and nothing said how.
    A spin-up is a RUN, not a button, so the panel sets up the fields that
    make the next run one and names where to launch it."""
    page = open(PAGE, encoding='utf-8').read()
    assert 'mkspinup' in page, 'no way to set up a spin-up'
    for field in ('spinup.cycles', 'spinup.save_strt', 'spinup.strt_heads'):
        assert "park('%s'" % field in page, '%s is not set up' % field


def test_the_run_really_does_start_without_usable_state(tmp_path):
    """The claim the panel now makes, asked of the code that decides it."""
    import marmites_props as props

    cfg = cfgmod.RunConfig.from_dict({})
    cfg.grid.kind = 'voronoi'
    cfg.spinup.strt_heads = 'hi_spinup'      # no sidecar, no files
    kind, payload, why = props.resolve_initial_heads(cfg, str(tmp_path),
                                                     verbose=False)
    assert kind == 'dem', 'the run refuses instead of starting cold'
    assert payload == props.DEFAULT_STRT_DEM
    assert 'starts' in why


def test_a_tail_panel_is_named_not_numbered():
    """panel_name fell back to "panel 8", putting the number back into the
    prose the function exists to keep it out of."""
    assert schema.panel_name(7) == 'Run'
    assert schema.panel_name(8) == 'Results'
