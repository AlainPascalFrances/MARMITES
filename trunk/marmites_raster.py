# -*- coding: utf-8 -*-
"""Grid-agnostic sampling of ESRI ASCII rasters onto a MODFLOW grid (Phase 4).

The legacy reader (`clsPROCESS.convASCIIraster2array`) requires the raster to
match the structured grid cell-for-cell -- it asserts identical nrow/ncol and
cellsize. That is fine for DIS, and it stays untouched.

For DISV (and for refined grids in general) the model cells no longer line up
with the raster, so values must be *sampled* at the cell centres. This module
does that with no third-party dependency (`rasterio`/`shapely` are optional in
the target environment and absent in CI):

  * `read_esri_ascii`  -- header + array, keeping the georeference;
  * `sample_points`    -- nearest-cell lookup at arbitrary (x, y);
  * `sample_to_cells`  -- values for a MARMITES cell list on any grid.

Zone rasters (soil, meteo, irrigation, vegetation ids) MUST use nearest
sampling -- interpolating integer class codes is meaningless. `sample_to_cells`
is nearest-only by design; add bilinear explicitly if a continuous field ever
needs it.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import os

import numpy as np

__all__ = ['EsriAscii', 'read_esri_ascii', 'sample_points', 'sample_to_cells']


class EsriAscii:
    """An ESRI ASCII grid with its georeference.

    Attributes: array (nrows, ncols), ncols, nrows, xllcorner, yllcorner,
    cellsize, nodata. Row 0 is the NORTHERNMOST row (ESRI convention, same
    as MODFLOW DIS).
    """

    def __init__(self, array, ncols, nrows, xllcorner, yllcorner, cellsize, nodata):
        self.array = np.asarray(array)
        self.ncols, self.nrows = int(ncols), int(nrows)
        self.xllcorner, self.yllcorner = float(xllcorner), float(yllcorner)
        self.cellsize = float(cellsize)
        self.nodata = nodata

    @property
    def bounds(self):
        """(xmin, ymin, xmax, ymax)."""
        return (self.xllcorner, self.yllcorner,
                self.xllcorner + self.ncols * self.cellsize,
                self.yllcorner + self.nrows * self.cellsize)

    def __repr__(self):
        return ('<EsriAscii %dx%d cellsize=%g origin=(%g, %g)>'
                % (self.nrows, self.ncols, self.cellsize,
                   self.xllcorner, self.yllcorner))


def read_esri_ascii(path, dtype=float):
    """Read an ESRI ASCII grid (6-line header, then values)."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    hdr = {}
    with open(path) as f:
        for _ in range(6):
            pos = f.tell()
            parts = f.readline().split()
            if len(parts) != 2:
                f.seek(pos)
                break
            try:
                hdr[parts[0].lower()] = float(parts[1])
            except ValueError:
                hdr[parts[0].lower()] = parts[1]
        data = np.loadtxt(f, dtype=float)
    ncols = int(hdr['ncols'])
    nrows = int(hdr['nrows'])
    cellsize = float(hdr['cellsize'])
    nodata = hdr.get('nodata_value', -9999)
    if 'xllcorner' in hdr:
        xll, yll = float(hdr['xllcorner']), float(hdr['yllcorner'])
    else:  # xllcenter/yllcenter variant
        xll = float(hdr['xllcenter']) - 0.5 * cellsize
        yll = float(hdr['yllcenter']) - 0.5 * cellsize
    arr = np.asarray(data, dtype=float).reshape(nrows, ncols)
    if dtype is int:
        arr = np.rint(arr).astype(int)
    return EsriAscii(arr, ncols, nrows, xll, yll, cellsize, nodata)


def sample_points(ras, x, y, outside=np.nan):
    """Nearest-cell sample of `ras` at points (x, y).

    Points outside the raster get `outside`. Row 0 is northernmost, so the
    row index grows as y decreases.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xmin, ymin, xmax, ymax = ras.bounds
    col = np.floor((x - xmin) / ras.cellsize).astype(int)
    row = np.floor((ymax - y) / ras.cellsize).astype(int)
    inside = (x >= xmin) & (x < xmax) & (y > ymin) & (y <= ymax)
    col = np.clip(col, 0, ras.ncols - 1)
    row = np.clip(row, 0, ras.nrows - 1)
    vals = ras.array[row, col].astype(float)
    if outside is not None:
        vals = np.where(inside, vals, outside)
    return vals


def sample_to_cells(path, cell_xy, dtype=float, hnoflo=None, nodata_to=None):
    """Sample a raster onto arbitrary model cells (DIS, DISV, refined...).

    Parameters
    ----------
    path : str            ESRI ASCII raster.
    cell_xy : (ncell, 2)  cell-centre coordinates in model CRS.
    dtype : float | int   int forces nearest + integer output (zone rasters).
    hnoflo : float        value written where the raster is NODATA
                          (defaults to the raster's own nodata).
    nodata_to : float     alias for hnoflo (kept explicit for callers).

    Returns (ncell,) array.
    """
    ras = read_esri_ascii(path, dtype=float)
    xy = np.asarray(cell_xy, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError('cell_xy must be (ncell, 2)')
    fill = hnoflo if hnoflo is not None else (nodata_to if nodata_to is not None
                                              else ras.nodata)
    vals = sample_points(ras, xy[:, 0], xy[:, 1], outside=fill)
    # map raster NODATA to the requested fill
    if ras.nodata is not None:
        vals = np.where(np.isclose(vals, float(ras.nodata)), fill, vals)
    if dtype is int:
        # an int raster cannot carry NaN: non-finite fills become the raster's
        # integer NODATA (callers test against it, as with the legacy reader)
        bad = ~np.isfinite(vals)
        if bad.any():
            nod = int(ras.nodata) if ras.nodata is not None else -9999
            vals = np.where(bad, nod, vals)
        return np.rint(vals).astype(int)
    return vals.astype(float)


def cell_centres_structured(delr, delc, xorigin=0.0, yorigin=0.0, i_arr=None, j_arr=None):
    """Cell-centre coordinates for a structured grid (helper for testing and
    for sampling a DIS model through the same code path as DISV)."""
    delr = np.asarray(delr, dtype=float)
    delc = np.asarray(delc, dtype=float)
    x_edge = np.concatenate(([0.0], np.cumsum(delr))) + float(xorigin)
    y_top = float(yorigin) + float(np.sum(delc))
    y_edge = y_top - np.concatenate(([0.0], np.cumsum(delc)))
    xc = 0.5 * (x_edge[:-1] + x_edge[1:])
    yc = 0.5 * (y_edge[:-1] + y_edge[1:])
    if i_arr is None or j_arr is None:
        nrow, ncol = delc.size, delr.size
        i_arr = np.repeat(np.arange(nrow), ncol)
        j_arr = np.tile(np.arange(ncol), nrow)
    i_arr = np.asarray(i_arr, dtype=int)
    j_arr = np.asarray(j_arr, dtype=int)
    return np.column_stack([xc[j_arr], yc[i_arr]])
