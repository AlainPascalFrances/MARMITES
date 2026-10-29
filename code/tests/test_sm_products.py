# -*- coding: utf-8 -*-
"""Satellite soil moisture against the model's saturation Se (user,
2026-10-29, branch swi-soil-moisture; ppMF6/marmites_smrs): the Copernicus
Soil Water Index (0.1 deg, 1 km), the Copernicus Surface Soil Moisture 1 km
and SMOS-BEC 1 km, read from the files the user downloads, against the
catchment-mean Se of the run's own outputs."""
import importlib.util
import os
import sys
from datetime import timedelta

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
for _p in (CODE, os.path.join(CODE, 'ppMF6'), os.path.join(CODE, 'app')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import marmites_daymap as DM                         # noqa: E402
import marmites_smrs as SM                           # noqa: E402
from lib import results as R                         # noqa: E402

_spec = importlib.util.spec_from_file_location(
    '_tdm_for_sm', os.path.join(HERE, 'test_day_maps.py'))
TDM = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(TDM)
run = TDM.run                                        # the small MF6 model
START = TDM.START


def test_the_exponential_filter():
    """The recursive SWI filter: the first value as it is, then K_n =
    K_{n-1} / (K_{n-1} + exp(-dt/T)); a day without a value carries the
    last one; a constant stays constant."""
    days = [START + timedelta(days=i) for i in range(5)]
    out = SM.exp_filter(days, np.array([0.2, np.nan, 0.6, 0.6, np.nan]), 5)
    k1 = 1.0 / (1.0 + np.exp(-2 / 5))
    s1 = 0.2 + k1 * 0.4
    k2 = k1 / (k1 + np.exp(-1 / 5))
    s2 = s1 + k2 * (0.6 - s1)
    assert out.tolist() == pytest.approx([0.2, 0.2, s1, s2, s2])
    assert np.allclose(SM.exp_filter(days, np.full(5, 0.3), 10), 0.3)


def test_the_scores():
    days = [START + timedelta(days=i) for i in range(40)]
    x = 50 + 10 * np.sin(np.arange(40) / 5.0)
    m = SM.metrics(days, x, x + 5)
    assert (m['n'], m['bias']) == (40, pytest.approx(5.0))
    assert m['r'] == pytest.approx(1.0) and m['ubrmse'] == pytest.approx(0.0)
    assert m['r_anom'] == pytest.approx(1.0)
    assert SM.metrics(days[:2], x[:2], x[:2]) == {'n': 2}


def _se_mean(day):
    """The toy run's catchment-mean Se (TDM._h5: columns of 100 + 300 mm
    with S_r 25 mm and S_sat 130 mm on cells 0, 1, 3, 5; the third day's
    storage halved)."""
    ss = np.array([40.0, 80.0, 120.0, 160.0]) * (0.5 if day == 2 else 1.0)
    return np.clip((ss - 25) / 105, 0, 1)


def _nc(path, days, fields, crs_epsg=23029):
    """A product file over the toy model, 10 m pixels on its cells (row 0
    the top), its CRS written as the file's crs_wkt."""
    import pandas as pd
    import xarray as xr
    from pyproj import CRS
    ds = xr.Dataset({k: (('time', 'y', 'x'), v) for k, v in fields.items()},
                    coords={'time': pd.to_datetime(days),
                            'y': [15.0, 5.0], 'x': [5.0, 15.0, 25.0]})
    ds['crs'] = xr.DataArray(0, attrs={'crs_wkt': CRS.from_epsg(crs_epsg).to_wkt()})
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ds.to_netcdf(path)


def _pixels(day, f):
    """Per pixel ``f(Se)`` on the soil columns, ``f(mean Se)`` on the two
    pixels without soil -- so the catchment mean is ``f(mean Se)``."""
    se = _se_mean(day)
    v = np.full(6, f(se.mean()))
    v[[0, 1, 3, 5]] = f(se)
    return v.reshape(2, 3)


def test_an_swi_against_the_catchment_mean_se(run, tmp_path):
    """SWI = 0.95 (100 Se) + 2 everywhere (above 100 % would be a flag):
    r 1, the bias and ubRMSE of that line on the three days; the T-value
    asked for, or the nearest; the figure in the run's sm_products, in the
    calibration tab."""
    days = [START + timedelta(days=i) for i in range(3)]
    swi = np.stack([_pixels(i, lambda s: 95 * s + 2) for i in range(3)])
    x = np.array([100 * _se_mean(i).mean() for i in range(3)])
    y = 0.95 * x + 2
    _nc(str(tmp_path / 'swi' / 'c_gls_SWI_test.nc'), days,
        {'SWI_010': swi, 'SWI_020': swi + 1, 'QFLAG_010': swi * 0})
    info = dict(DM.outputs_of(run['out']), h5=TDM._h5(run['ws'], tmp_path),
                soil=TDM.SOIL, crs='EPSG:23029')
    r = SM.compare(run['out'], 'clms_swi_0.1deg', str(tmp_path / 'swi'),
                   info=info, T=12)
    st = r['stats']
    assert st['T'] == 10, 'the nearest T-value'
    assert st['pixels'] == 6 and sum(st['pixel_shares']) == pytest.approx(1)
    f = st['fit']
    assert f['n'] == 3 and f['r'] == pytest.approx(1.0)
    d = (y - y.mean()) - (x - x.mean())
    assert f['bias'] == pytest.approx(np.mean(y - x))
    assert f['ubrmse'] == pytest.approx(np.sqrt(np.mean(d ** 2)))
    assert set(st['by_T']) == {'10', '20'}
    assert len(st['per_pixel']) == 6
    name = os.path.basename(r['png'])
    assert name == 'smrs_clms_swi_0.1deg_T10.png'
    assert os.path.dirname(r['png']) == os.path.join(run['out'], '_output',
                                                     SM.SUBDIR)
    assert R.classify('_output/sm_products', name) == 'calib'
    assert R.title('_output/sm_products', name) == (
        "Soil moisture: the model's Se against the Copernicus Soil Water "
        "Index 0.1 deg (SWI T = 10 days)")
    assert '_output/sm_products' in R.FOLDERS
    # the run's outputs stay its own: the comparison is not one of the
    # files that date the run
    assert DM.outputs_of(run['out'])['why'] == ''


def test_a_surface_product_is_filtered_then_scaled(run, tmp_path):
    """SMOS-BEC in m3/m3 = 0.1 + 0.2 Se: with a filter of no memory, both
    scaled to their own range, the series are one."""
    days = [START + timedelta(days=i) for i in range(3)]
    sm = np.stack([_pixels(i, lambda s: 0.1 + 0.2 * s) for i in range(3)])
    _nc(str(tmp_path / 'bec' / 'BEC_SM_L4_test.nc'), days, {'SM': sm})
    info = dict(DM.outputs_of(run['out']), h5=TDM._h5(run['ws'], tmp_path),
                soil=TDM.SOIL, crs='EPSG:23029')
    r = SM.compare(run['out'], 'smos_bec_1km', str(tmp_path / 'bec'),
                   info=info, T=1e-3)
    assert r['stats']['fit']['r'] == pytest.approx(1.0)
    assert 'scaled' in r['stats']['comparison']


def test_why_a_product_cannot_be_compared(run, tmp_path):
    days = [START - timedelta(days=100 + i) for i in range(3)]
    _nc(str(tmp_path / 'old' / 'swi.nc'), days,
        {'SWI_010': np.ones((3, 2, 3)) * 50})
    info = dict(DM.outputs_of(run['out']), h5=TDM._h5(run['ws'], tmp_path),
                soil=TDM.SOIL, crs='EPSG:23029')
    with pytest.raises(ValueError, match='no day of this run'):
        SM.compare(run['out'], 'clms_swi_0.1deg', str(tmp_path / 'old'),
                   info=info)
    with pytest.raises(ValueError, match='no file'):
        SM.compare(run['out'], 'clms_ssm_1km', str(tmp_path / 'none'),
                   info=info)


def test_the_products_found(tmp_path):
    d = tmp_path / 'soil_moisture' / 'clms_swi_1km'
    d.mkdir(parents=True)
    (d / 'a_20180615.nc').write_bytes(b'')
    assert SM.available(str(tmp_path)) == {'clms_swi_1km': 1}
    assert set(SM.PRODUCTS) == {'clms_swi_0.1deg', 'clms_swi_1km',
                                'clms_ssm_1km', 'smos_bec_1km'}
