# -*- coding: utf-8 -*-
"""WP1b -- the app's lib layer, and WP1.7 the MMsurf parameter enumeration.

These are the parts of the front-end with substance, and they are deliberately
streamlit-free so they can be tested in the MODEL environment. The pages are
thin views over them; what is tested here is what could actually be wrong.
"""

import importlib.util
import json
import os
import sys
import time

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
REPO = os.path.abspath(os.path.join(CODE, '..'))
DS = os.path.join(REPO, 'example', 'LaMata')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


loaders = _load('app_loaders', os.path.join(CODE, 'app', 'lib', 'loaders.py'))
runlib = _load('app_runs', os.path.join(CODE, 'app', 'lib', 'runs.py'))
msc = _load('mmsurf_config', os.path.join(CODE, 'mmsurf_config.py'))

MMSURF_INI = os.path.join(DS, 'MMsurf_ws', '__inputMMsurf.ini')


# ------------------------------------------------------------------ loaders
def test_digest_changes_with_content():
    """The cache key must notice an edit, or the app shows stale data."""
    import tempfile
    p = os.path.join(tempfile.mkdtemp(), 'x.txt')
    with open(p, 'w') as fh:
        fh.write('a')
    d1 = loaders.digest(p)
    time.sleep(0.01)
    with open(p, 'w') as fh:
        fh.write('ab')
    assert loaders.digest(p) != d1


def test_read_asc_returns_the_real_grid_and_header():
    p = os.path.join(DS, 'MF_ws', 'elev.asc')
    if not os.path.exists(p):
        pytest.skip('elev.asc not present')
    import numpy as np
    arr, hdr = loaders.read_asc(p)
    assert (hdr['ncols'], hdr['nrows']) == (60, 65)
    assert hdr['cellsize'] == 50
    assert (hdr['xllcorner'], hdr['yllcorner']) == (739300, 4553050)
    assert arr.shape == (65, 60)
    # La Mata's DEM covers the whole rectangle, so this file has no nodata --
    # masking is tested on a synthetic grid below.
    assert 700 < np.nanmin(arr) < np.nanmax(arr) < 900


def test_read_asc_masks_nodata(tmp_path):
    """A nodata sentinel must become NaN, or it is plotted as an elevation of
    -9999 and silently wrecks every colour scale."""
    import numpy as np
    p = tmp_path / 'g.asc'
    p.write_text('ncols 2\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 50\n'
                 'NODATA_value -9999\n1.0 -9999\n2.0 3.0\n', encoding='utf-8')
    arr, hdr = loaders.read_asc(str(p))
    assert hdr['nodata_value'] == -9999
    assert np.isnan(arr[0, 1])
    assert np.nanmin(arr) == 1.0 and np.nanmax(arr) == 3.0


def test_inventory_reports_missing_files_rather_than_dropping_them():
    inv = loaders.inventory(DS)
    assert inv and all(len(items) for _g, items in inv)
    flat = {rel: ok for _g, items in inv for rel, _s, ok in items}
    assert flat.get('inputSOILthick.asc') is True
    assert any(k.startswith('inputObsHEADS_') for k in flat)   # glob expanded
    inv2 = loaders.inventory(os.path.join(DS, 'does_not_exist'))
    assert any(not ok for _g, items in inv2 for _r, _s, ok in items)


def test_read_table_separates_provenance_from_data():
    p = os.path.join(DS, 'inputPONDS.csv')
    if not os.path.exists(p):
        pytest.skip('inputPONDS.csv not generated')
    prov, header, rows = loaders.read_table(p)
    assert any('gis_to_dataset.py' in x for x in prov)
    assert header[:4] == ['fid', 'x', 'y', 'area_m2']
    assert len(rows) == 12
    assert all(not r[0].startswith('#') for r in rows)


def test_crs_note_names_the_unnamed_projection_case():
    """The La Mata DEM has an equivalent but unnamed PROJCS; reporting it as
    unknown would be misleading and reprojecting blindly would be wrong."""
    assert 'UNDECLARED' in loaders.crs_note(None)

    class _NoEpsg:
        def to_epsg(self):
            return None

    assert 'unnamed PROJCS' in loaders.crs_note(_NoEpsg())

    class _Epsg:
        def to_epsg(self):
            return 23029

    assert loaders.crs_note(_Epsg()) == 'EPSG:23029'


# --------------------------------------------------------------- run registry
def test_launch_is_detached_and_records_status(tmp_path):
    """A run must survive the page that started it, so it is a real detached
    process with a status file and a log."""
    cfgp = tmp_path / 'c.toml'
    cfgp.write_text('[meta]\nconfig_version = 1\n', encoding='utf-8')
    driver = tmp_path / 'fake_driver.py'
    driver.write_text('import sys, time\n'
                      'print("driver up", flush=True)\n'
                      'print(" ".join(sys.argv[1:]), flush=True)\n',
                      encoding='utf-8')
    runs = tmp_path / 'runs'
    rid, st = runlib.launch(str(cfgp), str(runs), overrides=['run.nsp=5'],
                            run_tag='unit', driver=str(driver),
                            python_exe=sys.executable, cwd=str(tmp_path))
    assert rid.endswith('_unit')
    assert st['state'] == 'running' and st['pid']
    assert '--set' in st['cmd'] and 'run.nsp=5' in st['cmd']
    for _ in range(100):
        if not runlib._pid_alive(st['pid']):
            break
        time.sleep(0.1)
    got = runlib.status(str(runs), rid)
    assert got['state'] == 'finished'
    assert got['outcome'] == 'completed'
    assert 'driver up' in runlib.log_tail(str(runs), rid)
    assert 'run.nsp=5' in runlib.log_tail(str(runs), rid)
    saved = json.loads((runs / rid / 'status.json').read_text(encoding='utf-8'))
    assert saved['run_id'] == rid


def test_launch_refuses_a_bad_override_and_a_missing_config(tmp_path):
    """The page must never build a command from free text."""
    cfgp = tmp_path / 'c.toml'
    cfgp.write_text('[meta]\nconfig_version = 1\n', encoding='utf-8')
    with pytest.raises(ValueError):
        runlib.launch(str(cfgp), str(tmp_path / 'r'), overrides=['--evil'])
    with pytest.raises(ValueError):
        runlib.launch(str(cfgp), str(tmp_path / 'r'), overrides=['no_equals'])
    with pytest.raises(FileNotFoundError):
        runlib.launch(str(tmp_path / 'nope.toml'), str(tmp_path / 'r'))


def test_failed_run_is_reported_as_failed(tmp_path):
    cfgp = tmp_path / 'c.toml'
    cfgp.write_text('[meta]\nconfig_version = 1\n', encoding='utf-8')
    driver = tmp_path / 'bad.py'
    driver.write_text('raise SystemExit("CONFIG ERROR: boom")\n', encoding='utf-8')
    runs = tmp_path / 'runs'
    rid, st = runlib.launch(str(cfgp), str(runs), driver=str(driver),
                            python_exe=sys.executable, cwd=str(tmp_path))
    for _ in range(100):
        if not runlib._pid_alive(st['pid']):
            break
        time.sleep(0.1)
    assert runlib.status(str(runs), rid)['outcome'] == 'failed'


def test_list_runs_is_newest_first(tmp_path):
    runs = tmp_path / 'runs'
    for name in ('20260101000001_a', '20260101000002_b'):
        d = runs / name
        d.mkdir(parents=True)
        (d / 'status.json').write_text(
            json.dumps({'run_id': name, 'state': 'finished'}), encoding='utf-8')
    got = [r['run_id'] for r in runlib.list_runs(str(runs))]
    assert got == ['20260101000002_b', '20260101000001_a']


# ------------------------------------------------------- WP1.7 MMsurf schema
def test_mmsurf_ini_parses_with_the_expected_counts():
    if not os.path.exists(MMSURF_INI):
        pytest.skip('MMsurf ini not present')
    ms = msc.load_mmsurf_ini(MMSURF_INI)
    assert ms.counts == {'NMETEO': 1, 'NVEG': 3, 'NCRP': 1, 'NFIELD': 1,
                         'NSOIL': 3}
    assert [n for n, _v in ms.veg] == ['grassMU', 'Qilex', 'Qpyr']
    assert [n for n, _v in ms.soil] == ['alluvium', 'regolith', 'outcrop']


def test_mmsurf_Zr_matches_what_the_model_actually_reads():
    """The parser must agree with what MMsoil consumes -- the check that makes
    the enumeration trustworthy rather than decorative.

    That used to mean ``__inputMMsurf4MMsoil.txt``, the file MMsurf wrote and
    the driver read back. WP1d removed it, because it was authoritative and
    disagreed with the ini; the configuration is now the single source, so
    this compares the ini against the CONFIGURATION instead.
    """
    if not os.path.exists(MMSURF_INI):
        pytest.skip('MMsurf ini not present')
    ref = os.path.join(REPO, 'code', 'configs', 'lamata.toml')
    if not os.path.exists(ref):
        pytest.skip('reference configuration not present')
    cfgmod = _load('cfg_for_app_lib', os.path.join(CODE, 'marmites_config.py'))
    ms = msc.load_mmsurf_ini(MMSURF_INI)
    cfg = cfgmod.load_run_config(ref)
    assert [n for n, _v in ms.veg] == [v.name for v in cfg.surface.vegetation]
    assert [v['Zr'] for _n, v in ms.veg] == \
        [v.root_depth for v in cfg.surface.vegetation]
    assert [n for n, _v in ms.soil] == [s.name for s in cfg.surface.soil]
    assert [v['por'] for _n, v in ms.soil] == \
        [s.porosity for s in cfg.surface.soil]


def test_every_mmsurf_parameter_carries_units_and_a_description():
    if not os.path.exists(MMSURF_INI):
        pytest.skip('MMsurf ini not present')
    params = msc.enumerate_parameters(msc.load_mmsurf_ini(MMSURF_INI))
    assert len(params) == 1 * 8 + 3 * 19 + 1 * 11 + 3 * 8      # 100
    for p in params:
        assert p.units and p.description, p.dotted
        assert isinstance(p.value, float)
    dotted = {p.dotted for p in params}
    assert 'veg.Qilex.Zr' in dotted and 'soil.alluvium.por' in dotted


def test_a_short_row_is_reported_against_a_named_field(tmp_path):
    """The ini is positional: a short line shifts everything after it. The
    parser must say which entry and which fields, not fail later somewhere
    unrelated."""
    p = tmp_path / 'bad.ini'
    p.write_text('#\n1\n41.0 6.1 798.0 0.0 1.0 3.0 6.0 6.0\n2\n'
                 'grass 0.01 0.5\n', encoding='utf-8')
    with pytest.raises(msc.MMsurfError) as e:
        msc.load_mmsurf_ini(str(p))
    assert 'veg entry 1' in str(e.value) and 'grass' in str(e.value)


# ------------------------------------------------------------- the view layer
def test_every_streamlit_page_compiles():
    """A page needs a Streamlit runtime to execute, so it cannot be run under
    pytest even though Streamlit is installed in this env -- so at least
    guarantee they parse: a syntax error in a page is otherwise only
    discovered on the server."""
    import glob
    import py_compile
    import tempfile
    pages = glob.glob(os.path.join(CODE, 'app', '**', '*.py'), recursive=True)
    assert pages, 'no app pages found'
    for f in pages:
        py_compile.compile(f, doraise=True, cfile=os.path.join(
            tempfile.mkdtemp(), 'x.pyc'))


def test_the_app_lib_layer_is_streamlit_free():
    """loaders.py and runs.py hold the substance and must stay testable in the
    MODEL environment, where Streamlit is not installed."""
    import re
    pat = re.compile(r'^\s*(?:import\s+streamlit|from\s+streamlit\b)', re.M)
    for name in ('loaders.py', 'runs.py'):
        p = os.path.join(CODE, 'app', 'lib', name)
        with open(p, encoding='utf-8') as fh:
            assert not pat.search(fh.read()), '%s imports streamlit' % name


def test_a_uniform_mesh_gets_one_area_class():
    """A structured grid has ONE cell area, so every quantile is the same
    number and the class edges collapse to a single value. That left no
    interval, no trace and a blank map -- the easiest mesh to draw was the
    one that did not get drawn."""
    import numpy as np

    got = loaders.area_classes(np.full(24, 2500.0))
    assert len(got) == 1
    lo, hi, sel = got[0]
    assert (lo, hi) == (2500.0, 2500.0)
    assert sel.all(), 'every cell has to fall in the single class'


def test_a_graded_mesh_gets_separated_classes():
    import numpy as np

    got = loaders.area_classes(np.linspace(10.0, 4000.0, 40))
    assert len(got) > 1
    assert sum(int(sel.sum()) for _lo, _hi, sel in got) == 40, \
        'every cell belongs to exactly one class'
    assert got[0][0] == 10.0 and got[-1][1] == 4000.0


def test_area_classes_of_nothing_is_nothing():
    assert loaders.area_classes([]) == []


def test_two_sizes_do_not_collapse_into_one_class():
    import numpy as np

    a = np.array([100.0] * 10 + [900.0] * 10)
    got = loaders.area_classes(a)
    assert len(got) >= 2
    assert sum(int(sel.sum()) for _lo, _hi, sel in got) == 20


# ------------------------------------------------------ the MODFLOW library
# paths.libmf6 is the BMI SHARED LIBRARY, not mf6.exe: the coupler steps
# MODFLOW one stress period at a time through the API, which the executable
# cannot do. Both sit in the same bin folder, so the two are easy to confuse
# and the confusion used to surface as "libmf6 not found" AFTER the build.

def _runner():
    return _load('_runner_libmf6',
                 os.path.join(CODE, 'tests', 'run_lamata_mf6.py'))


def test_a_folder_is_completed_to_the_library_inside_it(tmp_path):
    """The natural thing to paste is the bin directory."""
    r = _runner()
    binp = tmp_path / 'bin'
    binp.mkdir()
    lib = binp / 'libmf6.dll'
    lib.write_bytes(b'')
    got = r.resolve_libmf6(str(binp))
    assert os.path.normcase(got) == os.path.normcase(str(lib)), got
    assert r.check_libmf6(got) == got


def test_a_folder_without_a_library_says_so(tmp_path):
    r = _runner()
    empty = tmp_path / 'nothing'
    empty.mkdir()
    with pytest.raises(r.LibMF6Error) as exc:
        r.check_libmf6(r.resolve_libmf6(str(empty)))
    msg = str(exc.value)
    assert 'folder' in msg and 'libmf6.dll' in msg, msg


def test_the_executable_is_named_as_the_wrong_one(tmp_path):
    """mf6.exe is a whole simulation; the coupler needs the library."""
    r = _runner()
    binp = tmp_path / 'bin'
    binp.mkdir()
    (binp / 'mf6.exe').write_bytes(b'')
    (binp / 'libmf6.dll').write_bytes(b'')
    with pytest.raises(r.LibMF6Error) as exc:
        r.check_libmf6(r.resolve_libmf6(str(binp / 'mf6.exe')))
    msg = str(exc.value)
    assert 'EXECUTABLE' in msg and 'libmf6.dll' in msg, msg


def test_a_missing_library_points_at_auto(tmp_path):
    r = _runner()
    with pytest.raises(r.LibMF6Error) as exc:
        r.check_libmf6(r.resolve_libmf6(str(tmp_path / 'nope' / 'libmf6.dll')))
    assert 'auto' in str(exc.value)


def test_blank_stays_blank_and_is_not_an_error():
    """Blank means build and stop -- a choice, not a mistake."""
    r = _runner()
    assert r.resolve_libmf6('') == ''
    assert r.resolve_libmf6(None) == ''
    assert r.check_libmf6('') == ''


def test_quotes_from_a_pasted_path_are_stripped(tmp_path):
    r = _runner()
    lib = tmp_path / 'libmf6.dll'
    lib.write_bytes(b'')
    got = r.resolve_libmf6('"%s"' % lib)
    assert os.path.normcase(got) == os.path.normcase(str(lib)), got


def test_the_run_page_and_the_driver_use_one_resolver():
    """Two implementations of 'where is the library' would drift."""
    runs = _load('_runs_for_lib', os.path.join(CODE, 'app', 'lib', 'runs.py'))
    mod = runs.import_runner()
    assert hasattr(mod, 'resolve_libmf6') and hasattr(mod, 'check_libmf6')
    assert os.path.isfile(runs.driver_path())


# ------------------------------------------------------ the validation checks
# ONE answer serves the validation panel and the Launch button. These pin the
# levels down, because the level is what decides whether a run is blocked, and
# a check promoted from warning to error silently would stop every run.

def _checks():
    return _load('_checks_mod', os.path.join(CODE, 'app', 'lib', 'checks.py'))


def _cfg(**over):
    import marmites_config as mcfg
    cfg = mcfg.RunConfig.from_dict({})
    for dotted, val in over.items():
        sec, key = dotted.split('__')
        setattr(getattr(cfg, sec), key, val)
    return cfg


def test_a_level_that_does_not_exist_is_refused():
    chk = _checks()
    with pytest.raises(chk.CheckError):
        chk.Check('fatal', 'nope')


def test_worst_and_the_counts():
    chk = _checks()
    some = [chk.Check(chk.INFO, 'a'), chk.Check(chk.WARNING, 'b'),
            chk.Check(chk.ERROR, 'c')]
    assert chk.worst(some) == chk.ERROR
    assert chk.worst([]) == ''
    assert chk.worst(some[:2]) == chk.WARNING
    assert chk.count_by_level(some) == {chk.ERROR: 1, chk.WARNING: 1,
                                        chk.INFO: 1}


def test_a_switch_the_panels_contradict_is_an_error():
    """It is the one that made MMsurf run after being unplugged."""
    chk = _checks()
    got = list(chk.check_switches(_cfg(),
                                  unsaved=[('run.surface', False, True)]))
    assert len(got) == 1 and got[0].level == chk.ERROR
    assert 'run.surface' in got[0].title


def test_a_stress_period_cap_is_a_warning_not_an_error():
    """It does not stop a run; it makes the result mean something else."""
    chk = _checks()
    got = [c for c in chk.check_run_scope(_cfg(run__nsp=365))
           if c.key == 'run.nsp']
    assert len(got) == 1 and got[0].level == chk.WARNING
    assert '365' in got[0].title
    assert not [c for c in chk.check_run_scope(_cfg(run__nsp=0))
                if c.key == 'run.nsp']


def test_a_blank_library_is_a_note_and_not_a_problem():
    """Blank means build and stop -- a choice, not a mistake."""
    chk = _checks()
    got = list(chk.check_libmf6(_cfg(paths__libmf6='')))
    assert len(got) == 1 and got[0].level == chk.INFO


def test_a_library_that_is_not_there_is_an_error():
    chk = _checks()

    class _Runner(object):
        class LibMF6Error(Exception):
            pass

        @staticmethod
        def resolve_libmf6(given):
            return given

        @staticmethod
        def check_libmf6(path):
            raise RuntimeError('not there')

    got = list(chk.check_libmf6(_cfg(paths__libmf6='C:/nope/libmf6.dll'),
                                runner=_Runner))
    assert len(got) == 1 and got[0].level == chk.ERROR
    assert got[0].panel == chk.RUN_PANEL


def test_a_check_that_blows_up_becomes_an_error_and_not_a_crash():
    """The panel's whole job is to be reachable when something is wrong."""
    chk = _checks()

    def _boom(cfg, **kw):
        raise ZeroDivisionError('boom')
        yield                                        # pragma: no cover

    old = chk.CHECKS
    chk.CHECKS = (_boom,)
    try:
        got = chk.collect(_cfg())
    finally:
        chk.CHECKS = old
    assert len(got) == 1 and got[0].level == chk.ERROR
    assert '_boom' in got[0].title and 'boom' in got[0].title


def test_collect_puts_the_worst_first():
    chk = _checks()
    cfg = _cfg(run__nsp=10)
    got = chk.collect(cfg, unsaved=[('run.surface', False, True)])
    levels = [c.level for c in got]
    assert levels == sorted(levels, key=lambda l: {'error': 0, 'warning': 1,
                                                   'info': 2}[l])


def test_problems_returns_what_validate_raises():
    """The panel needs them one at a time; validate joins them."""
    import marmites_config as mcfg
    cfg = mcfg.RunConfig.from_dict({})
    assert cfg.problems() == []
    cfg.run.relax = 9.0
    errs = cfg.problems()
    assert errs and any('relax' in e for e in errs)
    with pytest.raises(mcfg.ConfigError) as exc:
        cfg.validate()
    for e in errs:
        assert e in str(exc.value)


def test_a_table_is_saved_through_the_same_save_as_the_fields():
    """Two ways to write the file would take turns discarding each other."""
    import shutil
    import marmites_config as mcfg
    editor = _load('_editor_mod', os.path.join(CODE, 'app', 'lib', 'editor.py'))
    src = os.path.join(CODE, 'configs', 'lamata.toml')
    tmp = os.path.join(CODE, 'configs', '_tbltest.toml')
    shutil.copy2(src, tmp)
    try:
        cfg = mcfg.load_run_config(tmp)
        rows = editor.table_rows(cfg, 'surface.vegetation')
        assert rows, 'no vegetation to edit'
        rows[0] = dict(rows[0])
        rows[0]['name'] = 'renamed_by_test'
        applied, _digest = editor.save(
            cfg, tmp, {'run.relax': 0.42},
            {'surface.vegetation': rows})
        assert any('run.relax' in a for a in applied), applied
        assert any('surface.vegetation' in a for a in applied), applied
        back = mcfg.load_run_config(tmp)
        assert back.run.relax == 0.42
        assert back.surface.vegetation[0].name == 'renamed_by_test'
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def test_an_unchanged_table_is_not_reported_as_a_change():
    import shutil
    import marmites_config as mcfg
    editor = _load('_editor_mod', os.path.join(CODE, 'app', 'lib', 'editor.py'))
    src = os.path.join(CODE, 'configs', 'lamata.toml')
    tmp = os.path.join(CODE, 'configs', '_tbltest2.toml')
    shutil.copy2(src, tmp)
    try:
        cfg = mcfg.load_run_config(tmp)
        rows = editor.table_rows(cfg, 'surface.vegetation')
        applied, _d = editor.save(cfg, tmp, {}, {'surface.vegetation': rows})
        assert applied == [], applied
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
