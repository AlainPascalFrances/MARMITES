# -*- coding: utf-8 -*-
"""Readers for the Streamlit app.  WP1b, steps 1b.2-1b.4.

DELIBERATELY STREAMLIT-FREE. Every function here is a plain function that can be
unit-tested in the model environment; the pages wrap them in ``st.cache_data``
themselves. That keeps the testable substance out of the view layer, and it is
also what lets this module be exercised on a machine where Streamlit is not
installed.

Caching contract: ``digest(path)`` returns (size, mtime_ns), which is what the
pages key their cache on -- the same idea as the cbc digest cache in
marmites_postprocess. Without it Streamlit re-reads a 180 MB coupled HDF5 on
every widget click, because it reruns the whole script each interaction.
"""

import os
from pathlib import Path

__all__ = ['digest', 'read_asc', 'read_table', 'provenance_of',
           'inventory', 'read_vector', 'crs_note', 'TIER_A_GROUPS']

# La Mata: ED50 / UTM 29N. Web maps need WGS84.
MODEL_EPSG = 23029
WEB_EPSG = 4326


def digest(path):
    """(size, mtime_ns) -- the cache key for everything read from disk."""
    st = os.stat(path)
    return (st.st_size, st.st_mtime_ns)


def read_asc(path):
    """ESRI ASCII grid -> (2-D masked-by-nodata array, header dict)."""
    import numpy as np
    hdr = {}
    with open(path, encoding='utf-8-sig') as fh:
        for _ in range(6):
            k, v = fh.readline().split()
            hdr[k.lower()] = float(v)
    arr = np.loadtxt(path, skiprows=6)
    nodata = hdr.get('nodata_value', -9999.0)
    arr = np.where(np.isclose(arr, nodata), np.nan, arr)
    return arr, hdr


def provenance_of(path):
    """The leading ``#`` comment block of a generated file, as a list."""
    out = []
    with open(path, encoding='utf-8-sig') as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            out.append(line.rstrip('\n'))
    return out


def read_table(path):
    """A converter CSV -> (provenance lines, header list, rows).

    Kept dependency-free on purpose: these files are small and the app should
    not need pandas to show them.
    """
    import csv
    prov = provenance_of(path)
    with open(path, encoding='utf-8-sig') as fh:
        rows = [r for r in csv.reader(
            line for line in fh if not line.startswith('#'))]
    header = rows[0] if rows else []
    return prov, header, rows[1:]


# The Tier-A inventory, grouped the way example/<case>/README.md groups it, so
# the page and the README cannot drift apart.
TIER_A_GROUPS = [
    ('Run control & parameters', [
        '__inputMM_v3.ini', '__inputMMsurf4MMsoil.txt',
        'MF_ws/__inputMF_flopy_v3_2s1L.ini', 'MF_ws/__inputMF_flopy_v3_2s3L.ini',
        'MF_ws/inputSOILparam.txt', 'MMsurf_ws/__inputMMsurf.ini']),
    ('Time & forcing series', [
        'inputDATE.txt', 'inputZON*_d.txt', 'inputZON_*_stp.txt',
        'MMsurf_ws/__meteoTB.txt', 'MMsurf_ws/__IRR_TS.txt']),
    ('Spatial parameter rasters', [
        'inputMETEOzones.asc', 'inputSOILzones.asc', 'inputIRRzones.asc',
        'inputVEG1area.asc', 'inputVEG2area.asc', 'inputVEG3area.asc',
        'inputSOILthick.asc', 'inputSTREAMw.asc', 'inputSTREAMhmax.asc']),
    ('Aquifer grids', [
        'MF_ws/elev.asc', 'MF_ws/elev_sinkfil.ASC', 'MF_ws/ibound_l1.asc',
        'MF_ws/ibound_l2.asc', 'MF_ws/hk_l1.asc', 'MF_ws/hk_l2.asc',
        'MF_ws/thick_l1.asc', 'MF_ws/thick_l2.asc', 'MF_ws/Sy_l1.asc',
        'MF_ws/sy_l2.asc', 'MF_ws/Ss_l1.asc', 'MF_ws/ss_l2.asc',
        'MF_ws/uzf_iuzfbnd.asc', 'MF_ws/drn_cond_l1.asc', 'MF_ws/drn_elev_l1.asc',
        'MF_ws/ghb_cond_l1.asc', 'MF_ws/ghb_head_l1.asc']),
    ('Saved spin-up state', [
        'MF_ws/hi_spinup_l1.asc', 'MF_ws/hi_spinup_l2.asc',
        'MF_ws/hi_spinup_perc.asc', 'MF_ws/hi_spinup_etg.asc']),
    ('Observations', [
        'inputObs.txt', 'inputObsHEADS_*.txt', 'inputObsSM_*.txt',
        'inputObsRo_catchment.txt']),
    ('Generated from the cartography (WP1.1)', [
        'inputSTREAM.csv', 'inputSTREAM_param.csv', 'inputPONDS.csv',
        'inputWATERSHED.csv']),
]


def inventory(dataset_dir):
    """Resolve TIER_A_GROUPS against a case folder.

    Returns [(group, [(relpath, size_bytes, exists), ...])]. Patterns with a
    ``*`` are expanded; a listed file that is missing is reported as missing
    rather than silently dropped, because a silently missing input is exactly
    what this page exists to catch.
    """
    import glob
    root = Path(dataset_dir)
    out = []
    for group, patterns in TIER_A_GROUPS:
        items = []
        for pat in patterns:
            if '*' in pat:
                hits = sorted(glob.glob(str(root / pat)))
                for h in hits:
                    rel = os.path.relpath(h, str(root)).replace('\\', '/')
                    items.append((rel, os.path.getsize(h), True))
                if not hits:
                    items.append((pat, 0, False))
            else:
                p = root / pat
                items.append((pat, p.stat().st_size if p.exists() else 0,
                              p.exists()))
        out.append((group, items))
    return out


def crs_note(crs):
    """Describe a CRS, naming the unnamed-but-equivalent case explicitly.

    The La Mata DEM rasters carry an equivalent but UNNAMED PROJCS with no EPSG
    code, so reporting "unknown" would be misleading and reprojecting blindly
    would be wrong.
    """
    if crs is None:
        return 'UNDECLARED (no .prj) -- assumed EPSG:%d' % MODEL_EPSG
    try:
        code = crs.to_epsg()
    except Exception:
        code = None
    if code:
        return 'EPSG:%d' % code
    return ('unnamed PROJCS, no EPSG code -- treated as EPSG:%d'
            % MODEL_EPSG)


def read_vector(path, to_web=True):
    """Read a vector layer for MAPPING. Returns (gdf, source_crs_note).

    geopandas is imported HERE, in the app, never on the model path (D3).
    """
    import geopandas as gpd
    gdf = gpd.read_file(path)
    note = crs_note(gdf.crs)
    if gdf.crs is None:
        gdf = gdf.set_crs(epsg=MODEL_EPSG, allow_override=True)
    if to_web:
        gdf = gdf.to_crs(epsg=WEB_EPSG)
    return gdf, note


def vector_layers(gis_dir):
    """The Tier-B layers the Inputs page offers, and whether each is present."""
    names = [
        ('hydrography.shp', 'stream network', '#1f77b4'),
        ('lm_ponds.shp', 'ponds (charcas)', '#17becf'),
        ('Limite.shp', 'catchment boundary', '#d62728'),
        ('202109MonitPts.shp', 'monitoring points', '#2ca02c'),
        ('202109ObsPts.shp', 'observation points', '#9467bd'),
        ('Irr_Fields.shp', 'irrigated fields', '#bcbd22'),
        ('ECtower_footprint.shp', 'EC tower footprint', '#ff7f0e'),
    ]
    root = Path(gis_dir)
    return [(n, label, colour, (root / n).exists()) for n, label, colour in names]


# --------------------------------------------------------------------- #
# the mesh  (WP1c.7)
# --------------------------------------------------------------------- #

def mesh_cache_paths(ws_root, kind):
    """Where a built mesh is cached, per grid kind.

    The driver writes it under the mesh's own MF6 workspace, so the front-end
    can show the grid a run WOULD use without loading the model -- building it
    needs the .ini, the time discretisation and Triangle, which is a ~20 s
    round trip the page should not make on every rerun.
    """
    sub = {'structured': 'MF6_ws', 'disv': 'MF6_ws_disv'}.get(
        kind, 'MF6_ws_%s' % kind)
    d = Path(ws_root) / sub / '_mesh'
    return d / ('mesh_%s.json' % kind), d / ('mesh_%s.sig.json' % kind)


def read_mesh(ws_root, kind):
    """Cached DISV gridprops for one grid kind, or None if not built yet.

    Returns ``(gridprops, signature)``.
    """
    import json
    grid_fn, sig_fn = mesh_cache_paths(ws_root, kind)
    if not grid_fn.exists():
        return None, None
    with open(grid_fn, encoding='utf-8') as fh:
        gp = json.load(fh)
    sig = None
    if sig_fn.exists():
        try:
            with open(sig_fn, encoding='utf-8') as fh:
                sig = json.load(fh)
        except ValueError:
            pass
    return gp, sig


def mesh_polygons(gridprops):
    """``(polygons, areas, centres)`` for a DISV gridprops.

    Polygons are lists of (x, y) in model CRS -- TRUE cell outlines, not the
    display raster the figure suite uses, so the page shows the mesh as it
    actually is.
    """
    import numpy as np
    vxy = {int(v[0]): (float(v[1]), float(v[2]))
           for v in gridprops['vertices']}
    polys, areas, centres = [], [], []
    for rec in gridprops['cell2d']:
        n = int(rec[3])
        p = [vxy[int(iv)] for iv in rec[4:4 + n]]
        a = np.asarray(p, dtype=float)
        x, y = a[:, 0], a[:, 1]
        polys.append(p)
        areas.append(0.5 * abs(float(np.dot(x, np.roll(y, -1))
                                     - np.dot(y, np.roll(x, -1)))))
        centres.append((float(rec[1]), float(rec[2])))
    return polys, np.asarray(areas), np.asarray(centres)
