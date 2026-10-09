# -*- coding: utf-8 -*-
"""Phase-2 tests for the TOML config module and the legacy-ini converter."""
import importlib.util
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.join(HERE, '..')
DATASET_INI = os.path.join(HERE, '..', '..', 'example', 'LaMata', '__inputMM_v3.ini')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cfgmod = _load('marmites_config', os.path.join(TRUNK, 'marmites_config.py'))


def test_convert_real_dataset_ini(tmp_path):
    if not os.path.exists(DATASET_INI):
        pytest.skip('DataSet_LaMata ini not present')
    toml_path = str(tmp_path / 'mm.toml')
    cfg = cfgmod.convert_ini_file(DATASET_INI, toml_path)
    # known values from the LaMata annotated ini (includes maxYearsTick fields)
    assert cfg.run_name.startswith('2s3L')
    assert cfg.verbose == 1
    assert cfg.iniMonthHydroYear == 10
    assert cfg.maxYearsTickTrimester == 10
    assert cfg.maxYearsTickSemester == 20
    assert cfg.plt_WB_unit == 'year'
    assert cfg.irr_yn is True
    assert cfg.gridIRR_fn == 'inputIRRzones.asc'
    assert cfg.xllcorner == 739300.0
    assert cfg.yllcorner == 4553050.0
    # a TOML file was written and reloads to an equivalent config
    assert os.path.exists(toml_path)
    cfg2 = cfgmod.load_mm_config(toml_path)
    assert cfg2.run_name == cfg.run_name
    assert cfg2.maxYearsTickSemester == cfg.maxYearsTickSemester
    assert cfg2.xllcorner == cfg.xllcorner


def test_roundtrip_toml(tmp_path):
    cfg = cfgmod.MMConfig(run_name='demo', iniMonthHydroYear=1, plt_WB_unit='day',
                          MMsoil_yn=-1, chunks=1)
    p = str(tmp_path / 'x.toml')
    cfgmod.dump_toml(cfg, p)
    back = cfgmod.load_mm_config(p)
    assert back.run_name == 'demo'
    assert back.iniMonthHydroYear == 1
    assert back.plt_WB_unit == 'day'
    assert back.MMsoil_yn == -1
    assert back.chunks == 1


def test_validation_rejects_bad_month():
    with pytest.raises(cfgmod.ConfigError):
        cfgmod.MMConfig(run_name='x', iniMonthHydroYear=13).validate()


def test_validation_rejects_bad_wb_unit():
    with pytest.raises(cfgmod.ConfigError):
        cfgmod.MMConfig(run_name='x', plt_WB_unit='week').validate()


def test_validation_rejects_empty_run_name():
    with pytest.raises(cfgmod.ConfigError):
        cfgmod.MMConfig(run_name='').validate()


def test_deprecated_fields_captured(tmp_path):
    if not os.path.exists(DATASET_INI):
        pytest.skip('DataSet_LaMata ini not present')
    cfg = cfgmod.convert_ini_file(DATASET_INI, str(tmp_path / 'm.toml'))
    # Picard-loop fields are parsed but quarantined as deprecated
    assert 'convcrit' in cfg.deprecated
    assert 'ccnum' in cfg.deprecated


# =========================================================================
# WP0 -- the run-configuration layer (RunConfig): defaults, unknown-key
# rejection, --set overrides, the guards, and provenance round-tripping.
# =========================================================================

REPO = os.path.join(HERE, '..', '..')
REF_CONFIG = os.path.join(REPO, 'code', 'configs', 'lamata.toml')


def test_empty_config_reproduces_todays_flag_defaults():
    """An empty file must behave exactly like the pre-WP0 defaults, so WP0
    stays a pure refactor -- with TWO deliberate exceptions, each of which
    a run can still ask for by name:

      grid.kind  'structured' -> 'voronoi'  (WP1d)
      seep.kind  'uzf' -> 'drn'             SIMULATE_GWSEEP is deprecated
                                            in MF6 (6.5.0); a new catchment
                                            must not inherit it in
                                            silence.
    """
    c = cfgmod.RunConfig.from_dict({})
    # the coupling is lagged and not asked (run.mode / run.relax retired
    # 2026-10-07)
    assert not hasattr(c.run, 'mode') and not hasattr(c.run, 'relax')
    assert c.run.nsp == 0                 # 0 = all stress periods
    assert c.run.daily is True
    assert c.run.ats is True
    assert c.run.max_discrepancy == 1.0
    assert c.grid.kind == 'voronoi'       # WP1d; was 'structured' == --grid dis
    assert c.layers.nlay == 6
    assert c.uzf.vks_scale == 1.0
    assert c.seep.kind == 'drn'       # WP1d; was 'uzf' == SIMULATE_GWSEEP
    assert c.seep.cond == 10000.0
    assert c.sfr.enable is False
    assert c.lak.enable is False
    assert c.lak.bedleak == 1e-3
    assert c.crr.enable is False
    assert c.spinup.cycles == 1
    assert c.spinup.tol == 0.05
    assert c.postproc.sankey_min_flux == 0.05
    assert c.postproc.sankey_full is True
    assert c.postproc.map_days == 6
    assert c.et.unsat_form == 'etwc'    # UZF ET is always simulated
    assert c.run_tag == '6lay'          # was 6lay_<mode> until 2026-10-07


def test_reference_config_loads_and_is_the_canonical_run():
    """The reference configuration still describes the La Mata run.

    WHAT IS ASSERTED HERE IS STRUCTURAL, never a working value. This file is
    EDITED THROUGH THE FRONT-END between test runs -- it is the modeller's
    live configuration, not a fixture -- so pinning a field they are meant to
    change turns ordinary work into a red suite. It has now happened three
    times: run.surface when MMsurf was unplugged, paths.libmf6 when the
    library was set, and spinup.strt_heads when state that did not fit the
    mesh was cleared to set up a spin-up. A test that needs a particular
    value must SET it -- see test_app_pages._force -- not find it.
    """
    if not os.path.exists(REF_CONFIG):
        pytest.skip('code/configs/lamata.toml not present')
    c = cfgmod.load_run_config(REF_CONFIG)
    assert c.layers.nlay == 2
    assert c.seep.kind == 'drn'
    # a working value -- 1e4 and 100 have both been run (2026-10-03)
    assert c.seep.cond > 0.0
    assert c.postproc.enable is True
    assert not hasattr(c.postproc, 'preproc'), 'one input-maps switch'
    # The spin-up fields NAME saved state; they are filled and cleared as the
    # grid changes, so that they are strings is the only durable claim.
    assert isinstance(c.spinup.strt_heads, str)
    assert isinstance(c.spinup.steady_means, str)


@pytest.mark.parametrize('bad', [
    {'run': {'moed': 'lagged'}},            # misspelt key
    {'grid': {'kind': 'structured', 'xx': 1}},
    {'sfr': {'enabel': True}},
])
def test_unknown_key_raises(bad):
    """A typo must never fall through to a default."""
    with pytest.raises(cfgmod.ConfigError) as e:
        cfgmod.RunConfig.from_dict(bad)
    assert 'unknown key' in str(e.value)


def test_unknown_section_raises():
    with pytest.raises(cfgmod.ConfigError) as e:
        cfgmod.RunConfig.from_dict({'runn': {}})
    assert 'unknown section' in str(e.value)


@pytest.mark.parametrize('bad,frag', [
    ({'et': {'evt_ramp': 0.0}}, 'evt_ramp'),
    ({'layers': {'nlay': 0}}, 'layers.nlay'),   # a count, so >= 1
    ({'seep': {'kind': 'drn', 'cond': 0.0}}, 'free-draining'),
    ({'grid': {'kind': 'nope'}}, 'grid.kind'),
    ({'crr': {'beta': 1.5}}, 'crr.beta'),
    ({'spinup': {'strt_dem': [1.0]}}, 'strt_dem'),
    ({'meta': {'config_version': 2}}, 'config_version'),
])
def test_validation_rejects_bad_values(bad, frag):
    with pytest.raises(cfgmod.ConfigError) as e:
        cfgmod.RunConfig.from_dict(bad)
    assert frag in str(e.value)


def test_groundwater_et_in_modflow_cannot_be_asked_for_at_all():
    """ETg is computed by MM and taken by the EVT packages; UZF must not
    remove it as well (cookbook WP2).

    This used to be a GUARD -- a field forced false. The field is gone, so
    what is asserted now is that it cannot be set at all: an unknown key
    is refused by name, and simulate_et is not fed by it either."""
    import dataclasses
    assert 'gwet_in_mf' not in {f.name
                                for f in dataclasses.fields(cfgmod.Et)}
    with pytest.raises(cfgmod.ConfigError) as e:
        cfgmod.RunConfig.from_dict({'et': {'gwet_in_mf': True}})
    assert 'gwet_in_mf' in str(e.value)


def test_the_coupled_model_requires_the_drain_seepage_face():
    """SIMULATE_GWSEEP is deprecated in MODFLOW 6 (6.5.0), and the drain
    reproduces it exactly; MARMITES reads the seepage back into the soil
    column every step, so it is held to the mechanism MODFLOW 6 maintains.
    MODFLOW alone may still use it."""
    with pytest.raises(cfgmod.ConfigError) as e:
        cfgmod.RunConfig.from_dict({'run': {'model': True},
                                    'seep': {'kind': 'uzf'}})
    assert 'seep.kind' in str(e.value)
    # ... and it is a coupling rule, not a ban
    c = cfgmod.RunConfig.from_dict({'run': {'model': False},
                                    'seep': {'kind': 'uzf'}})
    assert c.seep.kind == 'uzf'


def test_grid_dis_is_accepted_as_an_alias_for_structured():
    c = cfgmod.RunConfig.from_dict({'grid': {'kind': 'dis'}})
    assert c.grid_kind == 'structured'


@pytest.mark.parametrize('kind', list(cfgmod.GRID_KINDS))
def test_every_declared_grid_producer_exists(kind):
    """WP1c built them all, so the schema and the producers must agree. The
    guard stays wired for the next producer that gets declared before it
    works -- an unimplemented kind has to stop the run, not quietly do
    something else."""
    c = cfgmod.RunConfig.from_dict({'grid': {'kind': kind}})
    c.require_implemented_grid()          # must not raise
    import marmites_meshes
    assert kind in marmites_meshes._PRODUCERS


def test_an_unknown_grid_kind_is_rejected():
    """from_dict validates, so an unknown producer never reaches a run."""
    with pytest.raises(cfgmod.ConfigError) as e:
        cfgmod.RunConfig.from_dict({'grid': {'kind': 'hexagons'}})
    assert 'grid.kind' in str(e.value)


def test_set_overrides_coerce_to_the_replaced_type():
    c = cfgmod.RunConfig.from_dict({})
    c.apply_overrides(['run.nsp=365', 'sfr.enable=true', 'seep.cond=5000',
                       'meta.name=try', 'spinup.strt_dem=[0.9995, -2.0]'],
                      echo=False)
    assert c.run.nsp == 365 and isinstance(c.run.nsp, int)
    assert c.sfr.enable is True
    assert c.seep.cond == 5000.0
    assert c.meta.name == 'try'
    assert c.spinup.strt_dem == [0.9995, -2.0]


def test_set_override_rejects_unknown_key_and_bad_type():
    c = cfgmod.RunConfig.from_dict({})
    for bad in ('run.nope=1', 'nope.key=1', 'run=1'):
        with pytest.raises(cfgmod.ConfigError):
            c.apply_overrides([bad], echo=False)
    with pytest.raises(cfgmod.ConfigError):
        c.apply_overrides(['run.nsp=banana'], echo=False)
    with pytest.raises(cfgmod.ConfigError):
        c.apply_overrides(['run.ats=maybe'], echo=False)


def test_set_override_is_revalidated():
    """An override that produces an invalid configuration must raise, not
    quietly stand."""
    c = cfgmod.RunConfig.from_dict({})
    with pytest.raises(cfgmod.ConfigError):
        c.apply_overrides(['layers.nlay=0'], echo=False)


def test_config_hash_is_stable_and_sensitive():
    a = cfgmod.RunConfig.from_dict({})
    b = cfgmod.RunConfig.from_dict({})
    assert a.config_hash() == b.config_hash()
    b.apply_overrides(['layers.nlay=2'], echo=False)
    assert a.config_hash() != b.config_hash()
    # the hash must not depend on where the file came from
    c = cfgmod.RunConfig.from_dict({}, source_path='somewhere/else.toml')
    assert c.config_hash() == a.config_hash()


def test_resolved_config_round_trips(tmp_path):
    """The resolved configuration written into the run folder must load back
    to exactly the same configuration -- that is what makes it provenance."""
    if not os.path.exists(REF_CONFIG):
        pytest.skip('code/configs/lamata.toml not present')
    c = cfgmod.load_run_config(REF_CONFIG)
    c.apply_overrides(['run.nsp=30', 'sfr.enable=true'], echo=False)
    out = str(tmp_path / 'resolved_config.toml')
    c.write_toml(out)
    back = cfgmod.load_run_config(out)
    assert back.config_hash() == c.config_hash()
    assert back.run.nsp == 30 and back.sfr.enable is True


def test_param_source_accepts_a_bare_number_and_rejects_two_producers():
    c = cfgmod.RunConfig.from_dict({'sfr': {'manning': 0.04}})
    assert c.sfr.manning.producer() == 'value'
    assert c.sfr.manning.value == 0.04
    with pytest.raises(cfgmod.ConfigError):
        cfgmod.RunConfig.from_dict(
            {'sfr': {'rhk': {'value': 0.1, 'column': 'rhk'}}})
    with pytest.raises(cfgmod.ConfigError):
        cfgmod.RunConfig.from_dict({'sfr': {'rhk': {}}})


def test_drainage_width_producer_needs_a_and_b():
    """Drainage-scaled width (CdL's method) is the default width producer."""
    c = cfgmod.RunConfig.from_dict({})
    assert c.sfr.width.producer() == 'drainage'
    assert c.sfr.width.drainage == {'a': 0.5, 'b': 0.35}
    with pytest.raises(cfgmod.ConfigError):
        cfgmod.RunConfig.from_dict({'sfr': {'width': {'drainage': {'a': 0.5}}}})


def test_mm_paths_resolves_and_reports():
    import io
    paths = _load('mm_paths', os.path.join(TRUNK, 'mm_paths.py'))
    assert paths.REPO.exists()
    assert paths.dataset_dir('LaMata').name == 'LaMata'
    assert paths.dataset_dir('LaMata').parent.name == 'example'
    buf = io.StringIO()
    paths.report_paths('LaMata', stream=buf)
    text = buf.getvalue()
    assert 'REPO' in text and 'DATASET' in text and 'WS_ROOT' in text


def test_the_tools_are_found_wherever_the_layout_put_them(tmp_path,
                                                         monkeypatch):
    """2026-10-09: GRIDGEN and PEST++ were reported missing on the server,
    where MODFLOW_DIR holds win64/gridgen.exe and
    pestpp-5.2.27-win/bin/pestpp-ies.exe -- not the one layout mm_paths
    knew. A known layout first, then the file by name; the variable wins."""
    paths = _load('mm_paths', os.path.join(TRUNK, 'mm_paths.py'))
    for env in ('MM_GRIDGEN_EXE', 'MM_PESTPP_IES', 'MM_TRIANGLE_EXE'):
        monkeypatch.delenv(env, raising=False)
    root = tmp_path / 'MODFLOWandCo'
    for rel in ('win64/gridgen.exe', 'win64/triangle.exe',
                'pestpp-5.2.27-win/bin/pestpp-ies.exe',
                'pestpp-5.2.27-win/bin/pestpp-glm.exe'):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text('exe')
    assert paths._find_tool('GRIDGEN_EXE', root) == str(root / 'win64' /
                                                        'gridgen.exe')
    assert paths._find_tool('PESTPP_IES', root) == str(
        root / 'pestpp-5.2.27-win' / 'bin' / 'pestpp-ies.exe')
    assert paths._find_tool('TRIANGLE_EXE', root) == str(root / 'win64' /
                                                         'triangle.exe')
    # the x64 build of the known layout is preferred when both are there
    x64 = root / 'gridgen.1.0.02' / 'bin' / 'gridgen_x64.exe'
    x64.parent.mkdir(parents=True)
    x64.write_text('exe')
    assert paths._find_tool('GRIDGEN_EXE', root) == str(x64)
    monkeypatch.setenv('MM_GRIDGEN_EXE', 'D:/elsewhere/gridgen.exe')
    assert paths._find_tool('GRIDGEN_EXE', root) == 'D:/elsewhere/gridgen.exe'
    # nothing there: the known layout, which report_paths flags as missing
    assert paths._find_tool('PESTPP_IES', tmp_path / 'none') == str(
        tmp_path / 'none' / 'pestpp' / 'pestpp-ies.exe')


def test_the_run_and_its_figures_read_panel_zeros_dataset():
    """2026-10-09: plot_water_budget looked for the grid rasters in its own
    <repo>/example/LaMata while the run read panel 0's example_root, and
    the figures were skipped once the checkout's example folder was gone.
    The run's scripts take the dataset from mm_paths alone."""
    import re

    def code_of(name):
        src = open(os.path.join(HERE, name), encoding='utf-8').read()
        return '\n'.join(re.sub(r'#.*', '', ln) for ln in src.splitlines())

    for name in ('run_lamata_mf6.py', 'plot_water_budget.py',
                 'diagnose_coupling.py', 'make_quadtree_lamata.py'):
        code = code_of(name)
        assert not re.search(r"['\"]example['\"]\s*,\s*['\"]LaMata['\"]",
                             code), '%s builds the dataset path itself' % name
        assert "'E:'" not in code, '%s has a drive literal' % name
    drv = code_of('run_lamata_mf6.py')
    assert 'make_figures(a.ws, out_dir=a.out_dir, dataset_dir=DS)' in drv
    assert 'DS = str(mm_paths.dataset_dir(cfg.paths.case))' in drv


def test_a_retired_key_is_dropped_and_reported_not_refused():
    """Unknown keys raise, so DELETING a key would stop every existing file
    loading -- a hard refusal at launch over a setting that no longer does
    anything. et.uzf_et is the first: UZF always simulates unsaturated-zone
    ET, so there is no off position for it to be left in."""
    cfg = cfgmod.RunConfig.from_dict({'et': {'uzf_et': False,
                                           'unsat_form': 'etae'}})
    assert not hasattr(cfg.et, 'uzf_et')
    assert cfg.et.unsat_form == 'etae', 'the rest of the section was lost'
    said = ' '.join(cfg.migrated)
    assert 'et.uzf_et is gone' in said, cfg.migrated
    assert 'three sources' in said, 'the reason is not carried'


def test_a_genuinely_unknown_key_still_raises():
    """The retirement list must not become a hole a typo falls through."""
    with pytest.raises(cfgmod.ConfigError) as exc:
        cfgmod.RunConfig.from_dict({'et': {'uzf_ett': True}})
    assert 'uzf_ett' in str(exc.value)


def test_one_switch_for_the_input_maps():
    """The Plots panel showed "Input maps" twice: postproc.preproc ran the
    input stage and postproc.input_maps the parameter maps inside it, so
    on/off drew the general map alone. One switch now; an old file with
    preproc still loads, and says it was dropped."""
    cfg = cfgmod.RunConfig.from_dict({'postproc': {'preproc': True,
                                                 'input_maps': True}})
    assert not hasattr(cfg.postproc, 'preproc')
    assert cfg.postproc.input_maps is True
    assert 'postproc.preproc is gone' in ' '.join(cfg.migrated)
    import sys
    app = os.path.join(TRUNK, 'app')
    if app not in sys.path:
        sys.path.insert(0, app)
    from lib import schema
    labels = [schema.describe(d)[0] for d, _v in
              schema.fields_of(cfgmod.RunConfig.from_dict({}), 'postproc')]
    assert labels.count('Input maps') == 1, labels
    src = open(os.path.join(TRUNK, 'tests', 'run_lamata_mf6.py'),
               encoding='utf-8').read()
    assert 'preproc=cfg.postproc.input_maps' in src

def test_the_retired_coupling_and_route_keys_still_load():
    """2026-10-07: the coupling is lagged and groundwater ET is EVT's. A
    file that still says run.mode, run.relax or et.gw_route loads -- a hard
    refusal at launch over a setting that no longer does anything helps
    nobody -- with each drop reported, and a save writes the file without
    them."""
    c = cfgmod.RunConfig.from_dict({'run': {'mode': 'iterative', 'relax': 0.5},
                                    'et': {'gw_route': 'wel'}})
    said = ' '.join(c.migrated)
    for key in ('run.mode', 'run.relax', 'et.gw_route'):
        assert '%s is gone' % key in said, said
    assert c.run_tag == '%dlay' % c.layers.nlay
