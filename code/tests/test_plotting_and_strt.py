# -*- coding: utf-8 -*-
"""Tests for the recovered MARMITESplot module and DEM-based initial heads.

Task 1: MARMITESplot_v3 must work again under Python 3.12 / matplotlib 3.10.
        The blocker was `interval_type` compared against bytes literals
        (b'linspace'), which stopped matching once callers passed str: the
        `ticks` variable was then never assigned and the colorbar raised
        UnboundLocalError. Both str and bytes must now be accepted.

Task 2: initial heads may come from a DEM regression, head = a*elev + b,
        the usual first guess for a sedimentary basin. Starting far above the
        water table forces MODFLOW to drain the excess, which shows up as slow
        convergence and rejected infiltration.
"""
import importlib.util
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))
for p in ('', 'MARMITESutilities', 'MARMITESutilities/MARMITESplot',
          'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

pytest.importorskip('flopy')
import matplotlib  # noqa: E402
matplotlib.use('agg')
import matplotlib.pyplot as plt  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P = _load('mmplot_v3', os.path.join(TRUNK, 'MARMITESutilities', 'MARMITESplot',
                                    'MARMITESplot_v3.py'))
mf6mod = _load('marmites_mf6_strt', os.path.join(TRUNK, 'ppMF6', 'marmites_mf6.py'))


# --------------------------------------------------------------------- #
# Task 1: the recovered plotting module
# --------------------------------------------------------------------- #

def _layer_args(nrow=6, ncol=8, nlay=2):
    V = np.random.default_rng(0).random((1, nlay, nrow, ncol))
    mask = np.zeros((nlay, nrow, ncol), dtype=bool)
    mask[:, 0, 0] = True
    return V, mask, nrow, ncol, nlay


@pytest.mark.parametrize('interval_type', ['linspace', b'linspace', 'arange', b'arange'])
def test_plotlayer_accepts_str_and_bytes(interval_type, tmp_path):
    """The regression that broke the module: bytes-vs-str interval_type."""
    V, mask, nrow, ncol, nlay = _layer_args()
    P.plotLAYER(days=[0], str_per=[0], Date='NA', JD='NA', ncol=ncol, nrow=nrow,
                nlay=nlay, nplot=nlay, V=V, cmap=plt.cm.viridis, CBlabel='x (mm)',
                msg='', plt_title='T', MM_ws=str(tmp_path),
                interval_type=interval_type, interval_num=5, interval_diff=0.2,
                contours=False, Vmax=[V.max()], Vmin=[V.min()], ntick=5,
                fmt='%.2f', points=None, mask=mask, hnoflo=-999.9)
    assert any(f.endswith('.png') for f in os.listdir(str(tmp_path)))


def test_plotlayer_accepts_bytes_title(tmp_path):
    """The legacy driver passed bytes titles too."""
    V, mask, nrow, ncol, nlay = _layer_args()
    P.plotLAYER(days=[0], str_per=[0], Date='NA', JD='NA', ncol=ncol, nrow=nrow,
                nlay=nlay, nplot=nlay, V=V, cmap=plt.cm.viridis, CBlabel=b'x',
                msg=b'', plt_title=b'BYTES_TITLE', MM_ws=str(tmp_path),
                interval_type='linspace', interval_num=5, contours=False,
                Vmax=[V.max()], Vmin=[V.min()], ntick=5, fmt='%.2f',
                mask=mask, hnoflo=-999.9)
    assert any('BYTES_TITLE' in f for f in os.listdir(str(tmp_path)))


def test_as_str_helper():
    assert P._as_str(b'abc') == 'abc'
    assert P._as_str('abc') == 'abc'
    assert P._as_str(None) is None


def test_scatter_to_grid_uses_explicit_shape():
    """Grid shape must be explicit: active cells alone under-report it."""
    cells_ij = np.array([[0, 0], [1, 2]])
    g = P.scatter_to_grid(cells_ij, [1.0, 2.0], grid_shape=(4, 5))
    assert g.shape == (4, 5)
    assert g[0, 0] == 1.0 and g[1, 2] == 2.0
    assert np.isnan(g[3, 4])


def test_plot_heads_figure(tmp_path):
    heads = 700.0 + np.random.default_rng(1).random((20, 6))
    ij = np.array([[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1]])
    fn = str(tmp_path / 'heads.png')
    P.plotHEADS(heads, ij, (4, 3), fn)
    assert os.path.getsize(fn) > 0


def test_plot_coupling_figure(tmp_path):
    it = np.random.default_rng(2).integers(1, 40, size=30)
    rej = np.random.default_rng(3).random((30, 5))
    exf = np.random.default_rng(4).random((30, 5))
    fn = str(tmp_path / 'coup.png')
    P.plotCOUPLING(it, rejinf=rej, exf=exf, plt_export_fn=fn)
    assert os.path.getsize(fn) > 0
    P.plotCOUPLING(it, plt_export_fn=str(tmp_path / 'c2.png'))   # iterations only
    assert os.path.getsize(str(tmp_path / 'c2.png')) > 0


# --------------------------------------------------------------------- #
# Task 2: DEM-regression initial heads
# --------------------------------------------------------------------- #

@pytest.fixture(scope='module')
def cmf():
    if not os.path.exists(os.path.join(DS, 'MF_ws', '__inputMF_flopy_v3_2s3L.ini')):
        pytest.skip('La Mata dataset not present')
    import MARMITESutilities as MMutils
    import ppMODFLOW_flopy_v3 as ppMF
    c = ppMF.clsMF(MMutils.clsUTILITIES(verbose=1), MM_ws=DS, MM_ws_out=DS,
                   MF_ws=os.path.join(DS, 'MF_ws'),
                   MF_ini_fn='__inputMF_flopy_v3_2s3L.ini',
                   xllcorner=739300.0, yllcorner=4553050.0)
    c.outcropL = np.zeros((c.nrow, c.ncol), dtype=int)
    for L in range(c.nlay):
        ib = (np.abs(np.asarray(c.ibound))[L] != 0)
        c.outcropL += ((c.outcropL == 0) & ib) * (L + 1)
    c.nper, c.perlen, c.nstp = 2, [1, 1], [1, 1]
    return c


def _build(cmf, tmp, strt_from_dem=None):
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp), strt_from_dem=strt_from_dem)
    b.build()
    return b


def test_default_uses_ini_arrays(cmf, tmp_path):
    b = _build(cmf, tmp_path)
    assert b.strt_from_dem is None
    assert np.allclose(b.initial_heads(), np.asarray(cmf.strt, dtype=float))


def test_dem_regression_applied(cmf, tmp_path):
    a, off = 0.9995, -2.0
    b = _build(cmf, tmp_path, strt_from_dem=(a, off))
    strt = b.initial_heads()
    dem = np.asarray(np.ma.filled(np.asarray(cmf.elev), np.nan), dtype=float)
    active = np.isfinite(dem) & (b.idomain[0] > 0)
    expect = a * dem + off
    # equal where the regression sits above the layer bottom
    got = strt[0][active]
    exp = expect[active]
    close = np.isclose(got, exp, atol=1e-6)
    assert close.mean() > 0.95, 'regression must set the heads on active cells'
    # every layer starts from the same DEM-derived surface (before clipping)
    assert strt.shape == (cmf.nlay, cmf.nrow, cmf.ncol)


def test_dem_heads_stay_above_layer_bottom(cmf, tmp_path):
    """Newton must not start from a dry cell."""
    b = _build(cmf, tmp_path, strt_from_dem=(0.9995, -2.0))
    strt = b.initial_heads()
    botm = np.asarray(b.botm, dtype=float)
    ok = strt >= botm - 1e-9
    assert ok.all(), 'initial heads must not start below the layer bottom'


def test_dem_heads_are_lower_than_ini_heads(cmf, tmp_path):
    """The point of the change: start closer to the water table, not above it."""
    b0 = _build(cmf, tmp_path / 'a')
    b1 = _build(cmf, tmp_path / 'b', strt_from_dem=(0.9995, -2.0))
    ini = np.asarray(b0.initial_heads(), dtype=float)[0]
    dem = np.asarray(b1.initial_heads(), dtype=float)[0]
    active = b1.idomain[0] > 0
    assert np.nanmean(dem[active]) < np.nanmean(ini[active])


def test_strt_written_into_ic_file(cmf, tmp_path):
    b = _build(cmf, tmp_path, strt_from_dem=(0.9995, -2.0))
    b.write()
    ic = os.path.join(str(tmp_path), '%s.ic' % cmf.modelname.lower())
    assert os.path.exists(ic) and os.path.getsize(ic) > 0


def test_heads_asc_roundtrip(cmf, tmp_path):
    """A saved equilibrated head field must load back unchanged on active cells,
    so a spin-up result can seed later runs via --strt-heads."""
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float), sim_ws=str(tmp_path))
    b.build()
    botm = np.asarray(cmf.botm, dtype=float)
    field = botm + 30.0
    pref = os.path.join(str(tmp_path), 'hi_test')
    paths = b.save_heads_asc(field, pref)
    assert len(paths) == cmf.nlay and all(os.path.exists(p) for p in paths)
    back = b.load_heads_asc(pref)
    active = np.asarray(b.idomain) > 0
    assert np.allclose(back[active], field[active], atol=1e-4)
    # inactive cells become nodata -> nan on read
    assert np.isnan(back[~active]).all()
    # and it drives the IC when fed back as strt_array
    b.strt_array = np.nan_to_num(back, nan=float(np.nanmin(back)))
    assert np.allclose(b.initial_heads()[active], field[active], atol=1e-3)


def test_load_heads_asc_missing_file_raises(cmf, tmp_path):
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float), sim_ws=str(tmp_path))
    b.build()
    with pytest.raises(Exception, match='not found'):
        b.load_heads_asc(os.path.join(str(tmp_path), 'does_not_exist'))


def test_uzf_vks_scale_multiplies_only_uzf_vks(cmf, tmp_path):
    """Raising the UZF vks to offset the EPSILON clamp must scale the UZF
    packagedata vks column but leave the aquifer NPF k33 untouched."""
    import flopy

    def _build(scale):
        b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                          botm=np.asarray(cmf.botm, dtype=float),
                          sim_ws=str(tmp_path / ('s%s' % scale)))
        b.uzf_vks_scale = scale
        b.build()
        return b

    b1, b3 = _build(1.0), _build(3.0)
    v1 = np.array([row[5] for row in b1.uzf_packagedata])   # vks column
    v3 = np.array([row[5] for row in b3.uzf_packagedata])
    assert np.allclose(v3, 3.0 * v1), 'UZF vks must scale by the factor'
    # NPF k33 (aquifer) must be identical regardless of the UZF scale
    k1 = np.asarray(b1.gwf.get_package('npf').k33.array)
    k3 = np.asarray(b3.gwf.get_package('npf').k33.array)
    assert np.allclose(k1, k3), 'aquifer k33 must not change with the UZF scale'


def test_strt_array_overrides_everything(cmf, tmp_path):
    """A spin-up cycle feeds the previous cycle's final heads back as the IC.
    An explicit strt_array must win over both the ini arrays and the DEM
    regression, clipped only to stay above the layer bottom."""
    b = mf6mod.clsMF6(cmf, top=np.asarray(cmf.elev, dtype=float),
                      botm=np.asarray(cmf.botm, dtype=float),
                      sim_ws=str(tmp_path), strt_from_dem=(0.9995, -2.0))
    botm = np.asarray(cmf.botm, dtype=float)
    want = botm + 25.0                       # a plausible equilibrated field
    b.strt_array = want
    b.build()
    got = b.initial_heads()
    active = b.idomain > 0
    assert np.allclose(got[active], want[active], atol=1e-6)
    # and it is clipped above the bottom even if the array dips below
    b.strt_array = botm - 5.0
    assert (b.initial_heads() >= botm - 1e-9).all()
