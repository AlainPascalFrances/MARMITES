# -*- coding: utf-8 -*-
"""WP1.6 / WP1b.10 -- the repository rules, enforced by the suite.

Two rules that have each already cost time, and that are easy to break by
accident:

  1. The repository holds code, docs, and STRICTLY the input files MM and MF
     read directly. No GIS, no run output, however small or convenient.
     (Shapefiles landed in the repo once; that is why this test exists.)

  2. The model must never import the Streamlit app or Streamlit itself, so it
     keeps running headless from Spyder and from a PEST worker on a machine
     where Streamlit is not installed.
"""
import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CODE = os.path.join(REPO, 'code')

# Extensions that are never model input: source cartography and run output.
FORBIDDEN_EXT = ('.shp', '.shx', '.dbf', '.sbn', '.sbx', '.prj', '.cpg',
                 '.tif', '.tiff', '.gpkg', '.qgz', '.mxd',
                 '.h5', '.hds', '.cbc', '.grb', '.lst', '.jcb')
# Folder names that only ever hold generated output.
FORBIDDEN_DIRS = ('_output', 'figures_nwt_comparison', 'postproc', 'preproc',
                  'master', 'workers', 'ies_output')
# The one exception: pre-processing figures a run writes are under _input/ in
# the WORKSPACE, never in the repo, so _input is forbidden here too -- except
# that name is also used by nothing in the repo, so it stays in the list.
FORBIDDEN_DIRS = FORBIDDEN_DIRS + ('_input',)

SKIP_DIRS = {'.git', '__pycache__', '.pytest_cache', '.ruff_cache', 'venv',
             '.idea'}


def _walk():
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        yield root, dirs, files


def test_no_gis_or_run_output_in_the_repository():
    """Rule 1. `.prj` and friends travel with a shapefile, so banning the
    sidecars is what actually stops one being committed."""
    offenders = []
    for root, _dirs, files in _walk():
        for f in files:
            if f.lower().endswith(FORBIDDEN_EXT):
                offenders.append(os.path.relpath(os.path.join(root, f), REPO))
    assert not offenders, (
        'These files belong in the WORKSPACE or the DATA_ROOT, not the repo '
        '(rule 1: code, docs, and strictly what MM/MF read directly):\n  '
        + '\n  '.join(sorted(offenders)[:40]))


def test_no_run_output_folders_in_the_repository():
    offenders = []
    for root, dirs, _files in _walk():
        for d in dirs:
            if d in FORBIDDEN_DIRS:
                offenders.append(os.path.relpath(os.path.join(root, d), REPO))
    assert not offenders, (
        'Run output must go to $MM_WS_ROOT, never into the repository:\n  '
        + '\n  '.join(sorted(offenders)))


def _python_files(where, exclude=()):
    for root, dirs, files in os.walk(where):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS and d not in exclude]
        for f in files:
            if f.endswith('.py'):
                yield os.path.join(root, f)


def test_the_model_never_imports_streamlit_or_the_app():
    """Rule 2 (WP1b.10). `code/app/` may import the model; the model may never
    import `code/app/` or streamlit."""
    pat = re.compile(r'^\s*(?:import\s+streamlit|from\s+streamlit\b|'
                     r'import\s+app\b|from\s+app\b|from\s+code\.app\b)',
                     re.MULTILINE)
    offenders = []
    for path in _python_files(CODE, exclude={'app'}):
        with open(path, encoding='utf-8', errors='replace') as fh:
            if pat.search(fh.read()):
                offenders.append(os.path.relpath(path, REPO))
    assert not offenders, (
        'The model must stay runnable where Streamlit is not installed '
        '(Spyder, a PEST worker). These modules import it:\n  '
        + '\n  '.join(sorted(offenders)))


def test_geospatial_libraries_stay_off_the_model_path():
    """Decision D3: geopandas and rasterio are imported by the WP1 converter
    and by the Streamlit app, and nowhere else. The one historical exception is
    the optional general-map figure, which guards its import."""
    # The indentation must be CAPTURED, not swallowed by the match: a guarded
    # import (inside a try, or in a function) is fine and is how the optional
    # general-map figure does it. Only a module-level import is the fault.
    pat = re.compile(r'^([ \t]*)(?:import|from)\s+'
                     r'(geopandas|rasterio|fiona|shapely)\b', re.MULTILINE)
    allowed = {os.path.join('code', 'tools', 'gis_to_dataset.py')}
    offenders = []
    # SFR_LAK_CRR holds the CdL model verbatim, for reference only. It is not
    # part of the MARMITES build and is never imported by it.
    for path in _python_files(CODE, exclude={'app', 'tests', 'SFR_LAK_CRR'}):
        rel = os.path.relpath(path, REPO)
        if rel.replace('/', os.sep) in allowed:
            continue
        with open(path, encoding='utf-8', errors='replace') as fh:
            text = fh.read()
        for m in pat.finditer(text):
            if len(m.group(1)) == 0:             # a TOP-LEVEL import is the fault
                offenders.append('%s: %s' % (rel, m.group(0).strip()))
    assert not offenders, (
        'Geospatial libraries must not be imported at module level on the '
        'model path (decision D3). Move it to code/tools/ or guard it:\n  '
        + '\n  '.join(sorted(offenders)))


@pytest.mark.parametrize('name', ['inputSTREAMw.asc', 'inputSTREAMhmax.asc',
                                  'inputPONDw.asc', 'inputPONDhmax.asc'])
def test_the_stream_rasters_are_retired(name):
    """WP1d retired them; WP1.2 had only renamed them.

    They were never a channel map: the values came from ``Soil_type.shp``,
    where ``PONDw`` is 1.5 m on the two alluvium polygons and 0 elsewhere, so
    what the model called "the stream network" was the alluvium footprint,
    with one width for the whole catchment. The network is the mapped
    hydrography now (``inputSTREAM.csv``, burned onto the grid at run time),
    and the surface reservoir they also fed went to SFR and LAK.
    """
    ds = os.path.join(REPO, 'example', 'LaMata')
    if not os.path.isdir(ds):
        pytest.skip('example/LaMata not present')
    assert not os.path.exists(os.path.join(ds, name)), (
        '%s is back -- WP1d retired it. The stream network comes from '
        'inputSTREAM.csv and the width from [sfr] width.' % name)


def test_converter_outputs_are_present_and_carry_provenance():
    """WP1.1. Every generated Tier-A file says where it came from."""
    ds = os.path.join(REPO, 'example', 'LaMata')
    for name in ('inputSTREAM.csv', 'inputSTREAM_param.csv', 'inputPONDS.csv',
                 'inputWATERSHED.csv'):
        p = os.path.join(ds, name)
        if not os.path.exists(p):
            pytest.skip('%s not generated yet' % name)
        with open(p, encoding='utf-8') as fh:
            head = ''.join(fh.readline() for _ in range(6))
        assert 'generated by code/tools/gis_to_dataset.py' in head, name
        assert '# source' in head and '# crs' in head, name
