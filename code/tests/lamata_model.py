# -*- coding: utf-8 -*-
"""La Mata's model description for the tests -- from the configuration,
with no parameter file (2026-10-07).

The fixtures used to build ``clsMF`` by parsing La Mata's legacy MODFLOW
parameter file (in example/LaMata/MF_ws). A run reads no such file any
more, so the tests build the model the way the run does
(``marmites_props.model_from_config``), from La Mata's configuration with
its boundaries on the placement the parameter file had -- the drains and
the general-head boundary from their RASTERS. (The run's own configuration
places the drains by a line: tests/test_boundary_lines.py.)

``INI_*`` are what the parameter file produced, frozen from its last parse
on 2026-10-07, so the acceptance tests can still hold the panel to it
without the file.

Not a test module itself (no ``test_`` prefix): imported by them.
"""

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
REF = os.path.join(CODE, 'configs', 'lamata.toml')
DS = os.path.abspath(os.path.join(CODE, '..', 'example', 'LaMata'))
for _p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    if os.path.join(CODE, _p) not in sys.path:
        sys.path.insert(0, os.path.join(CODE, _p))


# ---- what the parameter file produced (frozen 2026-10-07) -----------------
# per layer: (sum, sum of squares, min, max) over the 65 x 60 dataset grid
INI_LAYERS = {
    'hk': [(236.60000000000005, 25.35000000000001, 0.05, 0.5),
           (236.60000000000005, 25.35000000000001, 0.05, 0.5)],
    'vka': [(7800.0, 15600.0, 2.0, 2.0),
            (7800.0, 15600.0, 2.0, 2.0)],
    'ss': [(0.3722000000000001, 3.7380000000000005e-05, 5e-05, 0.0004),
           (0.0014351999999999995, 5.615999999999998e-10, 1e-07, 4e-07)],
    'sy': [(39.0, 0.3900000000000001, 0.01, 0.01),
           (39.0, 0.3900000000000001, 0.01, 0.01)],
    'thick': [(84760.0, 1944800.0, 20.0, 45.0),
              (127140.0, 4375800.0, 30.0, 67.5)],
    'ibound': [(1870.0, 1870.0, 0.0, 1.0),
               (1954.0, 1954.0, 0.0, 1.0)],
}
# the twelve outlet drains: (layer, row, col, conductance)
INI_DRAINS = [(0, 2, 2, 0.035), (0, 3, 1, 0.035), (0, 4, 0, 0.035),
              (0, 5, 0, 0.035), (0, 6, 0, 0.035), (0, 7, 0, 0.035),
              (1, 2, 2, 0.025), (1, 3, 1, 0.025), (1, 4, 0, 0.025),
              (1, 5, 0, 0.025), (1, 6, 0, 0.025), (1, 7, 0, 0.025)]


# the unsaturated zone (its UZF1 block): eps 2.0 is below what MF6 allows,
# so the build has always clamped it to 3.5 -- the panel says 3.5
INI_UZF = {'ntrail2': 15.0, 'nsets': 500.0, 'surfdep': 0.25, 'thtr': 0.05,
           'thts': 0.45, 'thti': 0.15, 'iuzfopt': 2.0, 'eps': 2.0}


def layer_stats(arr):
    """(sum, sum of squares, min, max) per layer, as INI_LAYERS holds them."""
    a = np.asarray(arr, dtype=float)
    return [(float(np.sum(a[L])), float(np.sum(a[L] ** 2)),
             float(a[L].min()), float(a[L].max())) for L in range(a.shape[0])]


def same_stats(got, want, rel=1e-12):
    """Two layer_stats lists agree (to rounding of a different summation)."""
    return len(got) == len(want) and all(
        abs(g - w) <= rel * max(1.0, abs(w))
        for row_g, row_w in zip(got, want) for g, w in zip(row_g, row_w))


# ---- building it -----------------------------------------------------------

def available():
    """The dataset, its DEM and the configuration are all there."""
    import marmites_dem as mdem
    return (os.path.isdir(os.path.join(DS, 'MF_ws'))
            and os.path.exists(mdem.dem_path(DS)) and os.path.exists(REF))


def lamata_config():
    """configs/lamata.toml on the parameter file's placements: the drains
    by their rasters, not by a line, and the cold-start heads the file
    named (hi_topL1.asc for both layers)."""
    import marmites_config as mcfg
    c = mcfg.load_run_config(REF)
    c.drn.line = ''
    c.drn.cond = mcfg.VectorSource(raster='MF_ws/drn_cond_l%d.asc')
    c.ghb.line = ''
    c.layers.strt = mcfg.VectorSource(raster='MF_ws/hi_topL1.asc')
    return c


def lamata_cmf(cfg=None, boundaries=True, uzf=True, verbose=False):
    """La Mata's model description, as the run builds it before the forcing:
    the dataset grid, the land surface, the layers, the cold-start heads,
    the boundaries and the unsaturated zone -- plus the outcrop layer the
    driver derives. Skips the test when the dataset is not present."""
    if not available():
        pytest.skip('La Mata dataset (or its DEM / configuration) not present')
    import marmites_props as props
    cfg = cfg or lamata_config()
    c = props.model_from_config(cfg, DS, verbose=verbose)
    if boundaries:
        props.apply_boundaries(cfg, c, DS, verbose=verbose)
    if uzf:
        props.apply_uzf(cfg, c, DS, verbose=verbose)
    props.boundary_cell_counts(c)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    return c
