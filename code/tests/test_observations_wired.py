# -*- coding: utf-8 -*-
"""The State variables panel names the files the post-processing reads.

The audit behind the cookbook's Appendix B found nothing on that panel read
at all: the point table, the head / soil-moisture / runoff prefixes and the
point layer's name column were literals in five places of
marmites_postprocess. They are one table now, filled from [obs].
"""

import os
import shutil
import sys
import types

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(CODE, '..', 'example', 'LaMata'))
for _p in (CODE, os.path.join(CODE, 'ppMF6')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

pp = pytest.importorskip('marmites_postprocess')


@pytest.fixture(autouse=True)
def _restore():
    saved = dict(pp.OBS)
    yield
    pp.OBS.clear()
    pp.OBS.update(saved)


def _cfg(**obs):
    base = dict(table='inputObs.txt', heads_prefix='inputObsHEADS',
                sm_prefix='inputObsSM', ro_prefix='inputObsRo',
                aet_prefix='', name_column='Name', layer='')
    base.update(obs)
    return types.SimpleNamespace(obs=types.SimpleNamespace(**base))


def test_the_defaults_are_the_la_mata_names_the_literals_were():
    assert pp.OBS == {'table': 'inputObs.txt', 'heads': 'inputObsHEADS',
                      'sm': 'inputObsSM', 'ro': 'inputObsRo',
                      'name_column': 'Name'}


@pytest.mark.skipif(not os.path.exists(os.path.join(DS, 'inputObs.txt')),
                    reason='the La Mata dataset is not present')
def test_renamed_files_are_found_under_the_names_the_panel_gives(tmp_path):
    """The proof that the panel is read: rename the table and the head files
    in a scratch dataset, tell the panel, and the readers find them."""
    shutil.copy2(os.path.join(DS, 'inputObs.txt'),
                 str(tmp_path / 'wells.txt'))
    heads = [f for f in os.listdir(DS) if f.startswith('inputObsHEADS_')]
    assert heads, 'La Mata has no head series to rename'
    for f in heads:
        shutil.copy2(os.path.join(DS, f),
                     str(tmp_path / f.replace('inputObsHEADS_', 'piezo_')))
    # The baseline: the La Mata names, on the La Mata dataset. Not the count
    # of head FILES -- three of them belong to points disabled with '#' in
    # the table, which the reader skips by design.
    pp.use_observations(_cfg())
    want = sorted(p['name'] for p in pp.obs_points(DS)
                  if pp.obs_series(DS, p['name']) is not None)
    assert want, 'La Mata reads no head series even under its own names'
    # ... and the same points and series under the names the panel gives
    pp.use_observations(_cfg(table='wells.txt', heads_prefix='piezo'))
    got = sorted(p['name'] for p in pp.obs_points(str(tmp_path))
                 if pp.obs_series(str(tmp_path), p['name']) is not None)
    assert got == want, 'renamed: %s, original: %s' % (got, want)


def test_the_old_literal_is_no_longer_what_is_read(tmp_path):
    """With the panel pointing elsewhere, the La Mata name is NOT used."""
    (tmp_path / 'inputObs.txt').write_text('P1 0 0 1\n', encoding='utf-8')
    pp.use_observations(_cfg(table='elsewhere.txt'))
    with pytest.raises((IOError, OSError)):
        pp.obs_points(str(tmp_path))


def test_no_reader_still_names_a_file_itself():
    src = open(os.path.join(CODE, 'ppMF6', 'marmites_postprocess.py'),
               encoding='utf-8').read()
    code = '\n'.join(ln.split('#')[0] for ln in src.splitlines())
    body = code[code.index('def use_observations'):]
    for lit in ("'inputObs.txt'", "'inputObsHEADS'", "'inputObsSM'",
                "'inputObsRo'", "'inputObsHEADS_'", "gx['Name']"):
        assert lit not in body, '%s is still a literal after OBS' % lit


def test_the_run_takes_them_from_the_panel_before_reading_any():
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'),
               encoding='utf-8').read()
    assert '_pp.use_observations(a.config)' in src
    assert (src.index('_pp.use_observations(a.config)')
            < src.index('resolve_obs_cells(cMF, ctx, DS)'))
