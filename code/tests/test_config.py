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
    stays a pure refactor. WP1c is what flips grid.kind to voronoi."""
    c = cfgmod.RunConfig.from_dict({})
    assert c.run.mode == 'lagged'
    assert c.run.relax == 0.6
    assert c.run.nsp == 0                 # 0 = all stress periods
    assert c.run.daily is True
    assert c.run.ats is True
    assert c.run.max_discrepancy == 1.0
    assert c.grid.kind == 'structured'    # == the old --grid dis
    assert c.layers.nlay == 6
    assert c.layers.aggregate is False
    assert c.uzf.vks_scale == 1.0
    assert c.seep.kind == 'uzf'
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
    assert c.et.uzf_et is False           # WP2 turns this on
    assert c.run_tag == '6lay_lagged'


def test_reference_config_loads_and_is_the_canonical_run():
    if not os.path.exists(REF_CONFIG):
        pytest.skip('code/configs/lamata.toml not present')
    c = cfgmod.load_run_config(REF_CONFIG)
    # --nlay 2 --seep drn --strt-heads hi_spinup --steady-means hi_spinup
    # --preproc --postproc
    assert c.layers.nlay == 2
    assert c.seep.kind == 'drn'
    assert c.seep.cond == 10000.0
    assert c.spinup.strt_heads == 'hi_spinup'
    assert c.spinup.steady_means == 'hi_spinup'
    assert c.postproc.enable is True
    assert c.postproc.preproc is True


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
    ({'run': {'mode': 'nope'}}, 'run.mode'),
    ({'run': {'relax': 0.0}}, 'run.relax'),
    ({'layers': {'nlay': 3}}, 'layers.nlay'),
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


def test_gwet_guard_refuses_groundwater_et_in_modflow():
    """ETg is computed by MM and applied as a WEL sink; MODFLOW must not
    remove it as well (cookbook WP2)."""
    with pytest.raises(cfgmod.ConfigError) as e:
        cfgmod.RunConfig.from_dict({'et': {'gwet_in_mf': True}})
    assert 'gwet_in_mf' in str(e.value)


def test_grid_dis_is_accepted_as_an_alias_for_structured():
    c = cfgmod.RunConfig.from_dict({'grid': {'kind': 'dis'}})
    assert c.grid_kind == 'structured'


@pytest.mark.parametrize('kind', ['voronoi', 'quadtree'])
def test_unimplemented_grid_producer_fails_fast(kind):
    """Valid in the schema, but WP1c has not built the producer yet: the run
    must stop with a clear message rather than silently doing something else."""
    c = cfgmod.RunConfig.from_dict({'grid': {'kind': kind}})
    with pytest.raises(cfgmod.ConfigError) as e:
        c.require_implemented_grid()
    assert 'WP1c' in str(e.value)


def test_set_overrides_coerce_to_the_replaced_type():
    c = cfgmod.RunConfig.from_dict({})
    c.apply_overrides(['run.nsp=365', 'sfr.enable=true', 'seep.cond=5000',
                       'run.mode=iterative', 'spinup.strt_dem=[0.9995, -2.0]'],
                      echo=False)
    assert c.run.nsp == 365 and isinstance(c.run.nsp, int)
    assert c.sfr.enable is True
    assert c.seep.cond == 5000.0
    assert c.run.mode == 'iterative'
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
        c.apply_overrides(['layers.nlay=4'], echo=False)


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
