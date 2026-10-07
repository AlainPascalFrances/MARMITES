# -*- coding: utf-8 -*-
"""A run reads NO legacy MODFLOW parameter file (user, 2026-10-07).

``__inputMF_flopy_v3_*.ini`` used to be parsed first and overridden by the
panels section by section -- and still supplied the grid shape, the layer
count, the elevation raster, the initial heads, the cold-start recharge and
the raster names it loaded at parse time. A new catchment has no such file,
so the model description is built from the configuration and the dataset
alone: clsMF.from_config sets the scalars, marmites_props fills the arrays.
"""
import os
import re
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6',
           'app'):
    if os.path.join(CODE, _p) not in sys.path:
        sys.path.insert(0, os.path.join(CODE, _p))

import marmites_config as mcfg                                 # noqa: E402
import marmites_props as props                                 # noqa: E402

GRID = (739300.0, 4553050.0, 4, 3, 50.0)        # xll, yll, nrow, ncol, cell


def _cmf(tmp_path, nlay=2, recharge=2.0e-4):
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    ds = str(tmp_path)
    return ppMF.clsMF.from_config(
        MMutils.clsUTILITIES(verbose=0), MM_ws=ds, MM_ws_out=ds,
        MF_ws=os.path.join(ds, 'MF_ws'), grid=GRID, nlay=nlay,
        hnoflo=9999.999, modelname='toy', steady_recharge=recharge)


def test_the_model_description_needs_no_file(tmp_path):
    """An EMPTY folder: nothing to parse, and nothing is looked for."""
    c = _cmf(tmp_path, recharge=3e-4)
    assert c.MF_ini_fn is None
    assert (c.nlay, c.nrow, c.ncol) == (2, 4, 3)
    assert c.delr == [50.0] * 3 and c.delc == [50.0] * 4
    assert (c.xllcorner, c.yllcorner) == (739300.0, 4553050.0)
    assert c.hnoflo == 9999.999 and c.hdry is None
    assert (c.itmuni, c.lenuni) == (4, 2)
    assert c.perc_user == 3e-4                  # spinup.steady_recharge
    # what the figures read: every layer plotted, labelled by its number,
    # and the groundwater-ET terms drawn (wel_yn is MARMITESplot's old name)
    assert c.h_plt == [1, 1] and c.h_lbl == ['1', '2'] and c.wel_yn == 1
    assert list(c.drncells) == [0, 0] and c.ghb_yn == 0 and c.drn_yn == 0
    assert c.cPROCESS.nrow == 4 and c.cPROCESS.hnoflo == 9999.999


def test_the_driver_names_no_parameter_file():
    """The run's setup builds the model through props.model_from_config
    (clsMF.from_config) and opens no .ini: neither the constructor that
    parses one nor a file name -- in the driver, in that construction, or
    in the tests' shared fixture."""
    def code_of(path):
        src = open(path, encoding='utf-8').read()
        return '\n'.join(re.sub(r'#.*', '', ln) for ln in src.splitlines())

    drv = code_of(os.path.join(HERE, 'run_lamata_mf6.py'))
    prp = code_of(os.path.join(CODE, 'ppMF6', 'marmites_props.py'))
    assert 'props.model_from_config(' in drv
    assert 'ppMF.clsMF.from_config(' in prp
    for code in (drv, prp, code_of(os.path.join(HERE, 'lamata_model.py'))):
        assert not re.search(r'clsMF\(\s*(cUTIL|MMutils)', code), \
            'the parser is called'
        assert 'MF_ini_fn' not in code and '_2s1L' not in code
        assert not re.search(r"['\"][^'\"]*\.ini['\"]", code), 'an .ini is named'


def test_no_test_fixture_parses_the_parameter_file():
    """2026-10-07: the fixtures build La Mata the way the run does
    (tests/lamata_model.py); none may parse the legacy file again."""
    for name in sorted(os.listdir(HERE)):
        if not name.endswith('.py') or name == os.path.basename(__file__):
            continue
        src = open(os.path.join(HERE, name), encoding='utf-8',
                   errors='replace').read()
        code = '\n'.join(re.sub(r'#.*', '', ln) for ln in src.splitlines())
        assert 'MF_ini_fn=' not in code, '%s parses the parameter file' % name


def test_a_blank_layer_property_stops_the_run(tmp_path):
    cfg = mcfg.RunConfig.from_dict({'layers': {
        'nlay': 2, 'ibound': {'value': 1}, 'thickness': {'value': 10.0},
        'k33': {'value': 2.0}, 'ss': {'value': 1e-4}, 'sy': {'value': 0.1}}})
    c = _cmf(tmp_path)
    c.elev = np.full((4, 3), 800.0)
    with pytest.raises(props.PropertyError, match='layers.k is not answered'):
        props.apply_layer_properties(cfg, c, str(tmp_path), verbose=False,
                                     required=True)


def test_the_land_surface_needs_the_dem(tmp_path):
    cfg = mcfg.RunConfig.from_dict({'grid': {'dem': ''}})
    with pytest.raises(props.PropertyError, match=r'\[grid\] dem is blank'):
        props.land_surface(cfg, _cmf(tmp_path), str(tmp_path), verbose=False)
    cfg = mcfg.RunConfig.from_dict({'grid': {'dem': 'some_dem'}})
    with pytest.raises(props.PropertyError, match='not in the dataset'):
        props.land_surface(cfg, _cmf(tmp_path), str(tmp_path), verbose=False)


def test_an_active_cell_without_elevation_is_refused(tmp_path):
    c = _cmf(tmp_path)
    c.elev = np.full((4, 3), 800.0)
    c.elev[0, 0] = c.hnoflo
    c.ibound = np.ones((2, 4, 3), dtype=int)
    with pytest.raises(props.PropertyError, match='does not reach 1 active'):
        props.check_land_surface(c)
    c.ibound[:, 0, 0] = 0                       # inactive there: fine
    props.check_land_surface(c)


def test_cold_start_heads_from_layers_strt(tmp_path):
    cfg = mcfg.RunConfig.from_dict({'layers': {'strt': {'value': 795.5}}})
    c = _cmf(tmp_path)
    c.elev = np.full((4, 3), 800.0)
    what = props.initial_heads(cfg, c, str(tmp_path), verbose=False)
    assert 'layers.strt' in what
    assert c.strt.shape == (2, 4, 3) and np.all(c.strt == 795.5)


def test_cold_start_heads_from_the_dem_when_strt_is_blank(tmp_path):
    cfg = mcfg.RunConfig.from_dict({'spinup': {'strt_dem': [1.0, -3.0]}})
    c = _cmf(tmp_path)
    c.elev = np.full((4, 3), 800.0)
    c.elev[3, 2] = c.hnoflo                      # no elevation: stays hnoflo
    what = props.initial_heads(cfg, c, str(tmp_path), verbose=False)
    assert 'layers.strt is blank' in what
    assert np.all(c.strt[:, :3, :] == 797.0)
    assert np.all(c.strt[:, 3, 2] == c.hnoflo)


def test_where_a_run_starts(tmp_path):
    """Saved state if it fits; otherwise a COLD start from layers.strt, or
    the DEM rule when that is blank -- and the log says which."""
    blank = mcfg.RunConfig.from_dict({})
    assert props.resolve_initial_heads(blank, str(tmp_path),
                                       verbose=False)[0] == 'dem'
    strt = mcfg.RunConfig.from_dict({'layers': {
        'strt': {'raster': 'MF_ws/hi_%d.asc'}}})
    kind, payload, why = props.resolve_initial_heads(strt, str(tmp_path),
                                                     verbose=False)
    assert (kind, payload) == ('strt', None) and 'layers.strt' in why
    missing = mcfg.RunConfig.from_dict({
        'layers': {'strt': {'value': 790.0}},
        'spinup': {'strt_heads': 'no_such_state'}})
    kind, _p, why = props.resolve_initial_heads(missing, str(tmp_path),
                                                verbose=False)
    assert kind == 'strt' and 'no_such_state is named' in why


def test_the_boundary_counts_are_the_cells_built(tmp_path):
    c = _cmf(tmp_path)
    c.layer_row_column_elevation_cond = {0: [[0, 1, 1, 790.0, 1.0],
                                             [1, 1, 1, 760.0, 1.0],
                                             [1, 2, 1, 760.0, 1.0]]}
    c.layer_row_column_head_cond = {0: []}
    drn, ghb = props.boundary_cell_counts(c)
    assert list(drn) == [1, 2] and list(ghb) == [0, 0]


def test_the_cold_start_recharge_is_a_panel_field():
    assert mcfg.RunConfig.from_dict({}).spinup.steady_recharge == 2.0e-4
    with pytest.raises(mcfg.ConfigError, match='steady_recharge'):
        mcfg.RunConfig.from_dict({'spinup': {'steady_recharge': -1.0}})
    from lib import schema                     # the app's
    for key in ('spinup.steady_recharge', 'layers.strt'):
        assert key in schema.FIELDS, key


def test_la_mata_builds_its_description_without_the_file(tmp_path):
    """The real dataset, from its configuration: grid, land surface, every
    layer property and the cold-start heads -- with no parameter file
    involved (the one still lying in MF_ws is not read)."""
    import mm_paths
    cfg_fn = os.path.join(CODE, 'configs', 'lamata.toml')
    if not os.path.exists(cfg_fn):
        pytest.skip('no La Mata configuration')
    cfg = mcfg.load_run_config(cfg_fn)
    ds = str(mm_paths.dataset_dir(cfg.paths.case))
    rect = props.dataset_grid(ds)[0]
    import marmites_dem as mdem
    if rect is None or not cfg.grid.dem or not os.path.exists(
            mdem.dem_path(ds)):
        pytest.skip('La Mata dataset / DEM not present')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    c = ppMF.clsMF.from_config(
        MMutils.clsUTILITIES(verbose=0), MM_ws=ds, MM_ws_out=ds,
        MF_ws=os.path.join(ds, 'MF_ws'), grid=rect,
        nlay=int(cfg.layers.nlay), hnoflo=float(cfg.layers.hnoflo),
        modelname=cfg.meta.model_name(cfg.paths.case))
    props.land_surface(cfg, c, ds, verbose=False)
    props.apply_layer_properties(cfg, c, ds, verbose=False, required=True)
    props.check_land_surface(c)
    props.initial_heads(cfg, c, ds, verbose=False)
    active = (np.abs(np.asarray(c.ibound)) != 0).any(axis=0)
    assert c.botm.shape == (int(cfg.layers.nlay), c.nrow, c.ncol)
    assert np.all(np.asarray(c.elev)[active] > 500.0)   # La Mata ~730-800 m
    assert np.all(np.isfinite(c.strt))
