# -*- coding: utf-8 -*-
"""Cell areas, the same on both sides of the coupling (2026-09-24).

Every m3 <-> mm conversion between MMsoil and MODFLOW 6 divides by a cell
area. On La Mata's Voronoi mesh (0.0127 .. 4121 m2) they disagreed by up to
1.9 %: MARMITES computed its areas by the shoelace formula on absolute
coordinates (products ~3e12 m2, so a 0.0127 m2 cell lost 1.5 % to
rounding), and flopy wrote the DISV vertices with an 8-digit mantissa (a
northing of 4553208.97 m to the centimetre). UZF's ET read back on the
coupler's area came out above PET.
"""

import importlib.util
import os
import re
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'app')):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# a 0.1 m x 0.127 m cell at La Mata's coordinates: 0.0127 m2, the smallest
X0, Y0 = 741329.769123, 4553173.481234
CELL = [(X0, Y0), (X0 + 0.1, Y0), (X0 + 0.1, Y0 + 0.127), (X0, Y0 + 0.127)]


def _exact(pts):
    """The area of the polygon the doubles DESCRIBE, in exact rational
    arithmetic -- not 0.0127: X0 + 0.1 is itself rounded to the double."""
    from fractions import Fraction as F
    q = [(F(x), F(y)) for x, y in pts]
    return float(abs(sum(q[i][0] * q[i - 1][1] - q[i - 1][0] * q[i][1]
                         for i in range(len(q)))) / 2)


TRUE = _exact(CELL)


def _every_shoelace():
    import marmites_grid
    import marmites_mesh
    import marmites_vector
    ras = _load('marmites_rasterise_area',
                os.path.join(CODE, 'ppMF6', 'marmites_rasterise.py'))
    return {'marmites_grid.polygon_area': marmites_grid.polygon_area,
            'marmites_mesh._shoelace': marmites_mesh._shoelace,
            'marmites_vector._signed_area':
                lambda p: abs(marmites_vector._signed_area(p)),
            'marmites_rasterise._polygon_area': ras._polygon_area}


@pytest.mark.parametrize('name', sorted(_every_shoelace()))
def test_a_tiny_cell_far_from_the_origin_keeps_its_area(name):
    got = _every_shoelace()[name](CELL)
    assert got == pytest.approx(TRUE, rel=1e-12), name


def test_the_sign_still_says_the_orientation():
    import marmites_vector
    assert marmites_vector._signed_area(CELL) > 0.0
    assert marmites_vector._signed_area(CELL[::-1]) < 0.0


def test_the_panel_loader_too():
    from lib import loaders
    gp = {'vertices': [[k, x, y] for k, (x, y) in enumerate(CELL)],
          'cell2d': [[0, X0 + 0.05, Y0 + 0.06, 4, 0, 1, 2, 3]]}
    _polys, areas, _c = loaders.mesh_polygons(gp)
    assert areas[0] == pytest.approx(TRUE, rel=1e-12)


def test_no_absolute_shoelace_is_left():
    """One pattern, five copies, all fixed: a sixth would bring it back."""
    bad = re.compile(r'x, y = a\[:, 0\], a\[:, 1\]\s*\n\s*(return|polys)')
    for d, _dirs, files in os.walk(CODE):
        if '__pycache__' in d or os.sep + 'legacy' in d:
            continue
        for f in files:
            if f.endswith('.py'):
                s = open(os.path.join(d, f), encoding='utf-8',
                         errors='replace').read()
                assert not bad.search(s), os.path.join(d, f)


def test_modflow_reads_the_vertices_at_full_precision(tmp_path):
    """The DISV file written the way the build writes it holds the vertices
    to the last bit: 17 significant digits round-trip a double."""
    flopy = pytest.importorskip('flopy')
    src = open(os.path.join(CODE, 'ppMF6', 'marmites_mf6.py'),
               encoding='utf-8').read()
    assert 'sim.simulation_data.float_precision = 16' in src
    assert 'sim.simulation_data.float_characters = 24' in src
    sim = flopy.mf6.MFSimulation(sim_ws=str(tmp_path))
    sim.simulation_data.float_precision = 16
    sim.simulation_data.float_characters = 24
    flopy.mf6.ModflowTdis(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname='m')
    flopy.mf6.ModflowGwfdisv(
        gwf, nlay=1, ncpl=1, nvert=4, top=10.0, botm=[0.0],
        vertices=[[k, x, y] for k, (x, y) in enumerate(CELL)],
        cell2d=[[0, X0 + 0.05, Y0 + 0.06, 4, 0, 3, 2, 1]])
    sim.write_simulation(silent=True)
    txt = open(str(tmp_path / 'm.disv'), encoding='utf-8').read()
    block = txt[txt.index('BEGIN vertices'):txt.index('END vertices')]
    got = [tuple(float(v) for v in ln.split()[1:3])
           for ln in block.splitlines()[1:] if ln.strip()]
    assert got == CELL
