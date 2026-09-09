# -*- coding: utf-8 -*-
"""Grid-agnostic map plotting for MARMITES (Phase 4).

The legacy `MARMITESplot_v3.plotLAYER` draws maps with `imshow` on an
(nrow, ncol) image, which only exists for DIS. This helper uses flopy's
`PlotMapView`, which renders DIS and DISV identically, so the same call
produces a map whatever the grid is.

Scope (per the Phase-4 plan): a small helper for new outputs. Porting the
whole legacy plotting module is explicitly NOT part of Phase 4.

Typical use
-----------
    from marmites_plot import plot_cell_values
    plot_cell_values(modelgrid, cells, perc, 'Recharge (m/d)', 'rch.png')
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import numpy as np

__all__ = ['cells_to_grid_array', 'plot_cell_values']


def cells_to_grid_array(cells, values, ncpl, nodata=np.nan):
    """Scatter per-cell MARMITES values onto a full map-view array.

    cells  : MARMITES cell list [(cid, i, j, node), ...]
    values : (ncell,) values in cell order
    ncpl   : number of map cells of the target grid
    Returns (ncpl,) with `nodata` on inactive cells.

    Works for DIS and DISV alike because `node` is the flat map index in
    both (i*ncol+j == icell2d for the DIS-equivalent vertex grid).
    """
    values = np.asarray(values, dtype=float)
    if values.shape[0] != len(cells):
        raise ValueError('values must be one per active cell')
    out = np.full(int(ncpl), nodata, dtype=float)
    nodes = np.array([c[3] for c in cells], dtype=int)
    out[nodes] = values
    return out


def plot_cell_values(modelgrid, cells, values, label='', fname=None,
                     cmap='viridis', vmin=None, vmax=None, title=None,
                     contours=False, ax=None):
    """Map per-cell MARMITES values on any MODFLOW grid (DIS or DISV).

    modelgrid : flopy StructuredGrid or VertexGrid
    Returns the matplotlib Axes (saves to `fname` when given).
    """
    import matplotlib
    if fname is not None and matplotlib.get_backend().lower() != 'agg':
        matplotlib.use('agg')
    import matplotlib.pyplot as plt
    from flopy.plot import PlotMapView

    ncpl = int(np.atleast_1d(modelgrid.ncpl)[0]) if hasattr(modelgrid, 'ncpl') \
        else modelgrid.nrow * modelgrid.ncol
    arr = cells_to_grid_array(cells, values, ncpl)

    created = ax is None
    if created:
        _fig, ax = plt.subplots(figsize=(7, 7))
    pmv = PlotMapView(modelgrid=modelgrid, ax=ax, layer=0)
    plot_arr = np.ma.masked_invalid(
        arr.reshape(modelgrid.nrow, modelgrid.ncol)
        if getattr(modelgrid, 'grid_type', '') == 'structured' else arr)
    quad = pmv.plot_array(plot_arr, cmap=cmap, vmin=vmin, vmax=vmax)
    if contours:
        try:
            pmv.contour_array(plot_arr, colors='k', linewidths=0.4)
        except Exception:
            pass          # contouring is best-effort on irregular grids
    pmv.plot_grid(lw=0.15, color='0.6')
    cb = ax.get_figure().colorbar(quad, ax=ax, shrink=0.7)
    if label:
        cb.set_label(label)
    if title:
        ax.set_title(title, fontsize=10)
    ax.set_aspect('equal')
    if fname is not None:
        ax.get_figure().savefig(fname, dpi=140, bbox_inches='tight')
        if created:
            plt.close(ax.get_figure())
    return ax
