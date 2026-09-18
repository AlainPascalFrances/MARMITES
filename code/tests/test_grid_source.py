# -*- coding: utf-8 -*-
"""The grid the rasters stand on, not the one the parameter file claims.

nrow, ncol, delr, delc and the origin were declared twice -- in
``__inputMF_flopy_v3_*.ini`` and implicitly by every raster in the dataset
-- and nothing compared them. They agreed in La Mata by luck: all 42
rasters happen to carry the rectangle the ini names. A new catchment
would not be so lucky, and on a mesh the consequence is not a crash but a
model built on a displaced grid, because ``project_model`` uses these as
the SOURCE rectangle it projects from.
"""

import importlib.util
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
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


props = _load('marmites_props_g', os.path.join(CODE, 'ppMF6',
                                               'marmites_props.py'))


class FakeMF(object):
    """Only what check_grid reads."""

    def __init__(self, nrow=65, ncol=60, cell=50.0,
                 xll=739300.0, yll=4553050.0):
        self.nrow, self.ncol = nrow, ncol
        self.delr, self.delc = [cell] * ncol, [cell] * nrow
        self.xllcorner, self.yllcorner = xll, yll


needs_dataset = pytest.mark.skipif(
    not os.path.isdir(os.path.join(DS, 'MF_ws')),
    reason='the La Mata dataset is not present')


@needs_dataset
def test_the_rasters_declare_the_grid():
    rect, names, others = props.dataset_grid(DS)
    assert rect == (739300.0, 4553050.0, 65, 60, 50.0)
    assert len(names) > 1, 'one raster is not a consensus'
    assert not others, 'the dataset rasters disagree among themselves'


@needs_dataset
def test_agreement_is_checked_rather_than_assumed():
    assert props.check_grid(FakeMF(), DS, verbose=False) is not None


@needs_dataset
@pytest.mark.parametrize('kw,frag', [
    ({'xll': 739250.0}, 'xll'),        # the voronoi rectangle mismatch
    ({'yll': 4553000.0}, 'yll'),
    ({'nrow': 64}, 'nrow'),
    ({'ncol': 59}, 'ncol'),
    ({'cell': 25.0}, 'cell'),
])
def test_a_displaced_or_reshaped_grid_stops_the_run(kw, frag):
    """It CANNOT be corrected at this point -- clsMF has already sized and
    read every array with the parameter file's shape -- so the only honest
    thing left is to refuse."""
    with pytest.raises(props.GridMismatch) as e:
        props.check_grid(FakeMF(**kw), DS, verbose=False)
    assert frag in str(e.value)
    assert '42 raster(s) agree' in str(e.value), (
        'the message must say what the rasters actually declare')


def test_a_dataset_with_no_raster_has_nothing_to_disagree_with(tmp_path):
    """A brand new catchment. Refusing here would make it impossible to
    start one."""
    assert props.dataset_grid(str(tmp_path))[0] is None
    assert props.check_grid(FakeMF(), str(tmp_path), verbose=False) is None


def test_the_origin_is_no_longer_hard_coded_in_the_driver():
    """739300 / 4553050 were literals in the driver, repeated from the
    parameter file. One number, one place -- the raster headers."""
    src = open(os.path.join(HERE, 'run_lamata_mf6.py'),
               encoding='utf-8').read()
    assert 'xllcorner=739300.0' not in src
    assert 'props.dataset_grid(DS)' in src
    assert 'props.check_grid(' in src
