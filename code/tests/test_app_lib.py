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
    assert flat.get('inputSTREAMw.asc') is True      # renamed in WP1.2
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
    """The parser must agree with the handover file MMsoil consumes. This is
    the check that makes the enumeration trustworthy rather than decorative."""
    if not os.path.exists(MMSURF_INI):
        pytest.skip('MMsurf ini not present')
    ms = msc.load_mmsurf_ini(MMSURF_INI)
    zr = [v['Zr'] for _n, v in ms.veg]
    handover = os.path.join(DS, '__inputMMsurf4MMsoil.txt')
    if not os.path.exists(handover):
        pytest.skip('handover file not present')
    with open(handover, encoding='utf-8-sig') as fh:
        toks = [ln.split('#')[0].strip() for ln in fh
                if ln.split('#')[0].strip()]
    # the row after the vegetation-name row is Zr
    names = [n for n, _v in ms.veg]
    idx = next(i for i, t in enumerate(toks) if t.split() == names)
    assert [float(x) for x in toks[idx + 1].split()] == zr


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
    """The pages cannot be executed here (Streamlit lives in its own
    environment, on purpose), so at least guarantee they parse -- a syntax
    error in a page is otherwise only discovered on the server."""
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
