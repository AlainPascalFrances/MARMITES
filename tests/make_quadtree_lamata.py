# -*- coding: utf-8 -*-
"""Build a quadtree (refined) DISV grid for La Mata with GRIDGEN (Phase 4).

Refines around the DRN network (the stream/drain cells), which is where the
groundwater gradients and the seepage that MARMITES exchanges are steepest.

Requires the gridgen executable, e.g.
    C:\\00MODFLOW\\gridgen.1.0.02\\bin\\gridgen_x64.exe

Usage (Anaconda Prompt, env with flopy):
    python tests\\make_quadtree_lamata.py ^
        --gridgen C:\\00MODFLOW\\gridgen.1.0.02\\bin\\gridgen_x64.exe --level 2

Outputs into <ws-root>/MF6_ws_quadtree/:
    * gridgen scratch files and the DISV gridprops summary,
    * quadtree_grid.png  -- the refined grid with the refinement features,
    * gridprops.npz      -- vertices/cell2d/top/botm for reuse by the coupler.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..', 'trunk'))
WS_ROOT = os.environ.get('MARMITES_WS_ROOT', os.path.join('E:' + os.sep, '00code_ws', 'LaMata_MM-MF6'))
DS = os.path.abspath(os.path.join(HERE, '..', 'DataSet_LaMata'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

import matplotlib  # noqa: E402
matplotlib.use('agg')
import matplotlib.pyplot as plt  # noqa: E402
import MARMITESutilities as MMutils  # noqa: E402
import ppMODFLOW_flopy_v3 as ppMF  # noqa: E402
from marmites_gridgen import build_quadtree  # noqa: E402
from marmites_grid import VertexGeometry  # noqa: E402


def drain_points(cMF):
    """Cell-centre coordinates of the DRN cells (refinement features)."""
    if getattr(cMF, 'drn_yn', 0) != 1:
        return []
    delr = np.asarray(cMF.delr, dtype=float)
    delc = np.asarray(cMF.delc, dtype=float)
    x_edge = np.concatenate(([0.0], np.cumsum(delr))) + float(cMF.xllcorner)
    y_top = float(cMF.yllcorner) + float(np.sum(delc))
    y_edge = y_top - np.concatenate(([0.0], np.cumsum(delc)))
    xc = 0.5 * (x_edge[:-1] + x_edge[1:])
    yc = 0.5 * (y_edge[:-1] + y_edge[1:])
    seen, pts = set(), []
    for (_l, i, j, _e, _c) in cMF.layer_row_column_elevation_cond[0]:
        if (i, j) in seen:
            continue
        seen.add((i, j))
        pts.append([(float(xc[int(j)]), float(yc[int(i)]))])
    return pts


def obs_points(ds_dir):
    """Observation-point coordinates from inputObs.txt (name x y lay ...)."""
    fn = os.path.join(ds_dir, 'inputObs.txt')
    if not os.path.exists(fn):
        return []
    pts = []
    with open(fn) as f:
        first = f.readline().split()
        dc = first[0] if first else '#'
        for line in f:
            tok = line.split(dc)[0].split()
            if len(tok) >= 3:
                try:
                    pts.append([(float(tok[1]), float(tok[2]))])
                except ValueError:
                    continue
    return pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gridgen', required=True, help='path to gridgen_x64(.exe)')
    ap.add_argument('--level', type=int, default=2, help='refinement level')
    ap.add_argument('--features', choices=['drn', 'obs', 'both'], default='both',
                    help='refine around drains, observation points, or both')
    ap.add_argument('--ws', default=os.path.join(WS_ROOT, 'MF6_ws_quadtree'))
    a = ap.parse_args()

    cUTIL = MMutils.clsUTILITIES(verbose=1)
    cMF = ppMF.clsMF(cUTIL, MM_ws=DS, MM_ws_out=DS, MF_ws=os.path.join(DS, 'MF_ws'),
                     MF_ini_fn='__inputMF_flopy_v3_2s3L.ini',
                     xllcorner=739300.0, yllcorner=4553050.0)
    # aquifer top/botm as the model uses them (soil thickness subtracted later
    # in the coupled driver; for grid generation the raw surfaces are enough)
    cMF.top = np.asarray(cMF.elev, dtype=float)
    cMF.botm = np.asarray(cMF.botm, dtype=float)

    pts = []
    if a.features in ('drn', 'both'):
        d = drain_points(cMF)
        print('  %d drain cells' % len(d))
        pts += d
    if a.features in ('obs', 'both'):
        o = obs_points(DS)
        print('  %d observation points' % len(o))
        pts += o
    print('refinement: %d point features, level %d' % (len(pts), a.level))
    feats = [(pts, 'point', a.level)] if pts else []
    if not feats:
        print('WARNING: no refinement features -> the grid will equal the base grid')

    gp = build_quadtree(cMF, a.gridgen, a.ws, refine_features=feats, level=a.level)
    ncpl = int(gp['ncpl'])
    print('quadtree built: ncpl=%d (base grid was %d), nvert=%d, nlay=%d'
          % (ncpl, cMF.nrow * cMF.ncol, len(gp['vertices']), gp.get('nlay', cMF.nlay)))

    geom = VertexGeometry.from_vertices(gp['vertices'], gp['cell2d'],
                                        np.arange(ncpl), nlay=cMF.nlay)
    print('cell areas: min %.1f m2, max %.1f m2, mean %.1f m2'
          % (geom.area.min(), geom.area.max(), geom.area.mean()))

    np.savez_compressed(os.path.join(a.ws, 'gridprops.npz'),
                        vertices=np.array(gp['vertices'], dtype=object),
                        cell2d=np.array(gp['cell2d'], dtype=object),
                        ncpl=ncpl, area=geom.area, allow_pickle=True)

    # picture of the refined grid
    try:
        from flopy.discretization import VertexGrid
        from flopy.plot import PlotMapView
        vg = VertexGrid(vertices=gp['vertices'], cell2d=gp['cell2d'],
                        nlay=1, ncpl=ncpl)
        fig, ax = plt.subplots(figsize=(9, 9))
        pmv = PlotMapView(modelgrid=vg, ax=ax)
        pmv.plot_array(geom.area, cmap='viridis')
        pmv.plot_grid(lw=0.2, color='0.4')
        if pts:
            xy = np.array([p[0] for p in pts])
            ax.plot(xy[:, 0], xy[:, 1], 'r.', ms=1.5, label='DRN cells')
            ax.legend(loc='upper right', fontsize=8)
        ax.set_aspect('equal')
        ax.set_title('La Mata quadtree grid (level %d around drains)\n'
                     'ncpl=%d, cell area %.0f-%.0f m2'
                     % (a.level, ncpl, geom.area.min(), geom.area.max()), fontsize=10)
        png = os.path.join(a.ws, 'quadtree_grid.png')
        fig.savefig(png, dpi=150, bbox_inches='tight')
        print('grid map: %s' % png)
    except Exception as exc:      # plotting must never break generation
        print('WARNING: could not plot the grid (%r)' % (exc,))

    print('\nNext: the refined grid changes cell count and areas, so the '
          'MARMITES inputs must be resampled onto it\n'
          '(marmites_gridgen.resample_inputs + RefinedModel) before a coupled '
          'run on this grid.')


if __name__ == '__main__':
    main()
