"""
================================================================================
Groundwater Flow Model — Catchment "CdL"  (DRYAD / NbS)               v2 "fable"
FloPy 3.10 + MODFLOW 6.7 — Voronoi DISV grid, monthly transient 1981–2026
================================================================================
ARCHITECTURE (the validated design — see the referenced literature per block):

  GRID    Voronoi/DISV via Triangle (§3). Two designs, GRID_DESIGN=1|2:
            1 = SFRmaker style (Leaf et al. 2021): no stream refinement, the
                network is mapped onto the 100 m background (~3.3k cells, fast).
            2 = smooth graded 40 m stream corridor -> 100 m background
                (Daoud et al. 2022 transitions; ~6.6k cells, production).
          Ponds: ONE pond-scale cell each via centroid seeding (Voronoi cells
          cannot follow rim polygons — faces bisect them). Degenerate slivers
          are mass-conservatively merged at build time.
  LAYERS  nlay=3 = U1 alluvium / U2 terraces / U3 Plio-Miocene (one numerical
          layer per geologic unit, K = 25/40/1.5 m/d, all convertible). Unit
          presence from build_layers.py's own npz (never mixed-source!); absent
          units pinch to idomain=-1; present-but-thin cells are extended to
          MIN_ACTIVE_THK=4 m (feather-edge conditioning). Base at -35 m.
  SURFACE UZF (per-cell soil rasters; rain/ET forcing on the TOP ACTIVE layer
          of each column, top_act — never "layer 0") + DRN-SEEP smoothed
          seepage faces -> MVR + DRN outlets at the boundary.
  STREAMS SFR built directly (§8b DEM routing + §12): one reach per cell,
          Manning rectangular, drainage-scaled width (arbolate-sum equivalent),
          downstream-MONOTONIC bed elevations (SFRmaker rule).
  PONDS   19 charcas as LAK, EMBEDDEDV — THE DEFINITIVE CONNECTION TYPE (user
          decision 2026-07-04; a VERTICAL connection's bed = the CELL TOP, so a
          dug-in pond would be permanently "dry" — TM6-A55 p7-13/7-18). One
          embedded connection in the seeded pond cell; 4-col stage/vol/sarea/
          barea wedge table (graceful drying); surfdep-smoothed; equilibrium
          initial stage from the SS water table; one EXTERNAL Manning spill
          outlet per lake (flopy gotcha: external = lakeout -1, NOT 0).
          On-channel ponds: SFR reaches inside footprints are EXCISED and the
          stream is routed THROUGH the lake by MVR (inlet reach -> LAK;
          LAK spill -> downstream reach).
  CRR     Daoud et al. 2022 cascade-routing & reinfiltration (§13, CRR_ENABLE):
          rejected infiltration (UZF) + exfiltration (DRN-SEEP) are MVR-routed
          to downslope neighbours by MFD slope factors (Eq. 23, β=CRR_BETA);
          receiver priority: LAK pond > SFR reach > UZF cell (reinfiltration);
          topographic sinks evaporate. CRR_ENABLE=False = legacy nearest-reach.

RUN RECIPES (the guards enforce the order and stop misconfigurations):
  IC run     : SPINUP_MODE="ss"        (ponds auto-stripped; writes spinup_heads.npy)
  Production : SPINUP_MODE="transient" (all 19 ponds; needs the matching IC)
  Grid switch: set GRID_DESIGN + REBUILD_GRID=True, run with SPINUP_MODE="ss"
               (ONE run: rebuilds the grid, AUTO-runs build_layers.py on the
               ncpl mismatch, regenerates the IC), then REBUILD_GRID=False and
               run the transient. The grid cache carries a design tag; a
               design/cache mismatch stops with a clear error.

Units: metres, days.  CRS: EPSG:3763.  Workspace: E:\\00code_ws\\DRYAD\\CdL_model
History: v1 "opus" (Claude Opus 4.7, 2026-06) -> v2 "fable" (Claude Fable 5,
2026-07-20: coherency review — retired dead ends removed: excavated & on-top
VERTICAL lake modes, LAK-seepage->UZF leg, RCH runoff-boost fallback).
Author: alain.frances@lneg.pt
Written by Claude Fable 5
================================================================================
"""

import os
import pickle
from pathlib import Path
from datetime import datetime
from collections import deque, defaultdict

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from shapely.geometry import Polygon, MultiPolygon, LineString, Point
from shapely.ops import unary_union, linemerge
from shapely.validation import make_valid
from shapely.strtree import STRtree
from scipy.spatial import cKDTree

import flopy
from flopy.discretization import VertexGrid
from flopy.utils.voronoi import VoronoiGrid
from flopy.utils.triangle import Triangle
from flopy.utils.gridintersect import GridIntersect

# -----------------------------------------------------------------------------
# 0. PATHS & GLOBAL SETTINGS
# -----------------------------------------------------------------------------
MF6_EXE = r"C:\00MODFLOW\mf6.7.0_win64\bin\mf6.exe"
TRIANGLE_EXE = r"C:\00MODFLOW\win64\triangle.exe"

WATERSHED_GPKG = r"E:/zzCloud/OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia/DRYAD/GIS/dryad_modelo_NbS.gpkg"
WATERSHED_LAYER = "watershed_cdl_fixed"
# Streams now live in the GeoPackage (user-fixed layer); routing is derived in section 8b.
STREAMS_LAYER = "streams_cdl"
OBS_PTS_LAYER = "obs_points_cdl"

# DEM-based SFR routing parameters (section 8b)
BRIDGE_TOL = 20.0   # m: snap dangling endpoints onto nearest line within this gap (T-junctions)
SNAP_TOL   = 1.0    # m: endpoints closer than this collapse to a single node

DEM_TIF = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\GIS\Geodatabase_LIDAR_DGT\Geodatabase_CdL\dem_cdl.tif"

P_CSV  = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\WP3_modeling\3.1.1\modelo_numerico\p_month_198101_202605.csv"
ET0_CSV = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\WP3_modeling\3.1.1\modelo_numerico\et0_month_198101_202605.csv"

# Soil-hydraulic rasters (DRYAD\dados\solos) -> per-cell UZF parameters (section 4b).
SOLOS_DIR     = r"E:\zzCloud\OneDrive - LNEG - Laboratorio Nacional de Energia e Geologia\DRYAD\dados\solos"
SOIL_VKS_TIF  = SOLOS_DIR + r"\ks_cdl.tif"     # saturated K           -> UZF vks
SOIL_THTS_TIF = SOLOS_DIR + r"\ths_cdl.tif"    # saturated water cont. -> thts
SOIL_THTR_TIF = SOLOS_DIR + r"\wp_cdl.tif"     # wilting point         -> thtr (residual)
SOIL_THTI_TIF = SOLOS_DIR + r"\fc_cdl.tif"     # field capacity        -> thti (initial)
# Unit conversions: ths/wp/fc rasters are in PERCENT (-> *0.01); ks is in cm/d
# (-> *0.01 = m/d), per the data provider.
SOIL_WC_TO_FRAC = 0.01
SOIL_VKS_TO_MD  = 0.01     # ks_cdl.tif is cm/d

WORKSPACE = Path(r"E:\00code_ws\DRYAD\CdL_model")
WORKSPACE.mkdir(exist_ok=True, parents=True)
# Pre-processing outputs (model map, parameter maps) -> timestamped subfolder
# (time the script starts).  The stamp is written to disk so post-processing lands
# in a matching postproc\<stamp>\ folder for the same run.
RUN_STAMP   = datetime.now().strftime("%Y%m%d%H%M")
PREPROC_DIR = WORKSPACE / "preproc" / RUN_STAMP
PREPROC_DIR.mkdir(parents=True, exist_ok=True)
(WORKSPACE / "last_run_stamp.txt").write_text(RUN_STAMP)

# --- SPIN-UP MODE / INITIAL HEADS (hybrid equilibration) ---------------------
# "transient": all periods transient — the lakes converge this way (a steady-state
#   solve won't close with lakes), but a uniform initial water table does NOT fully
#   equilibrate the deep/upland cells in a 1-yr spin-up (they drain over the run).
# "ss": steady-state period 0 (equilibrates the water table) — converges only
#   WITHOUT lakes.
# HYBRID workflow for the lake model:
#   (A) POND_SUBSET=[]  + SPINUP_MODE="ss"        -> converges, writes spinup_heads.npy
#   (B) POND_SUBSET=[…] + SPINUP_MODE="transient" + INIT_HEADS_FILE=SPINUP_HEADS_OUT
#       -> the lake run starts from the equilibrium water table (no upland drainage).
SPINUP_MODE      = "transient"          # ◀ THE ONLY FLAG TO FLIP between run types:
                                 #   "ss"        = no-lake SS+full-period run; regenerates spinup_heads.npy (the IC). Ponds OFF (auto).
                                 #   "transient" = the all-19 embedded-LAK run from the SS IC. Ponds ON (auto).
                                 # Do NOT touch the POND_SUBSET guard below — it follows this flag by itself.
SPINUP_HEADS_OUT = WORKSPACE / "spinup_heads.npy"   # SS-run equilibrium heads written here
INIT_HEADS_FILE  = None if SPINUP_MODE == "ss" else SPINUP_HEADS_OUT   # ss cold-starts; transient starts from the SS IC

MODEL_NAME = "cdl_gwf"
SIM_START = datetime(1981, 1, 1)
SIM_END   = datetime(2026, 5, 1)

# --- RUN LENGTH CONTROL ------------------------------------------------------
# How much of the period to simulate.  The 12-month spin-up (replays year 1) is
# ALWAYS prepended and is NOT counted here:
#   RUN_YEARS  < 0  -> full period (SIM_START .. SIM_END)
#   RUN_YEARS  > 0  -> short test run of this many years AFTER the spin-up
#                      (e.g. 5 -> 1 spin-up year + 5 test years = 6 years total)
# NOTE: changing this changes nper, so all period-aware files must be rewritten
#       -> run with WRITE_ALL=True (REBUILD_GRID can stay False; grid is unchanged).
RUN_YEARS = -1          # full period (>0 = that many years after the 12-month spin-up, for quick tests)
# Auto-run postprocess_cdl.py after a SUCCESSFUL run (a non-convergence raises before the hook).
RUN_POSTPROC = True   # set False to skip post-processing during quick solver/BC iterations
# -----------------------------------------------------------------------------

# ── GRID DESIGN (user 2026-07-04): TWO options, switch with GRID_DESIGN ───────────────────────
#  1 = SFRMAKER methodology (Leaf et al. 2021): NO stream refinement — the network is MAPPED
#      onto the CELL_FAR background (one reach per cell; robustness = reach attributes: DEM
#      elevations + downstream-monotonic smoothing + drainage-scaled width, not small cells).
#      COARSE & FAST (ncpl≈3295, full period ≈15 min) → the TEST grid.
#  2 = SMOOTH GRADED CORRIDOR (the grid the 45-yr coupled run validated, pre-"Daoud dimensions"):
#      40 m cells along the streams grading through Daoud-style transition zones to the 100 m
#      background (ncpl≈6609) → the PRODUCTION grid.
#  Ponds are centroid-seeded one-cell in BOTH. ⚠ SWITCHING DESIGNS = the full grid-rebuild chain
#  (REBUILD_GRID=True run → build_layers.py → SS run for the IC); the ncpl assert + the IC
#  fail-fast guards enforce the order, so a forgotten step stops with a clear error.
GRID_DESIGN = 2
if GRID_DESIGN == 1:        # option 1 — SFRmaker (Leaf 2021), coarse/fast
    STREAM_REFINE, CELL_NEAR_STREAM, CELL_FAR, STREAM_BUFFER = False, 20.0, 100.0, 60.0
    _TRANS_LEVELS_DESIGN = [20.0, 40.0, 70.0]           # pond grading only (no stream corridor)
else:                       # option 2 — graded 40 m corridor (the validated production grid)
    STREAM_REFINE, CELL_NEAR_STREAM, CELL_FAR, STREAM_BUFFER = True, 40.0, 100.0, 60.0
    _TRANS_LEVELS_DESIGN = [10.0, 20.0, 40.0, 70.0]     # corridor level (40) + transitions

# --- SMOOTH GRID TRANSITIONS (after Daoud et al. 2022) -----------------------
# Grade the Voronoi cell size GRADUALLY (CVFD requirement: the cell-centre line must
# ~⊥-bisect the shared edge — badly violated by an abrupt 2.5 m pond next to 40-150 m
# cells).  Concentric buffer-ring zones step the target size outward in ~2x increments.
SMOOTH_TRANSITIONS = True    # Daoud-style smooth grid grading via a CONTINUOUS sizing field
                             # (verified in preview_transition_grid.py).  Takes effect only on a grid REBUILD
                             # (set REBUILD_GRID=True -> then re-run build_layers.py + regenerate the SS IC).
TRANS_GRADE  = 0.5           # cell-size gradient: cell width grows ~GRADE m per 1 m of distance from a feature
TRANS_LEVELS = _TRANS_LEVELS_DESIGN   # set by GRID_DESIGN above: pond grading always; with STREAM_REFINE=True
                                      #   the levels >= CELL_NEAR_STREAM also grade the stream corridor outward.

# --- CATTLE PONDS (charcas) -> LAK package -----------------------------------
# Design 🔒 (user decision 2026-07-04, definitive): every pond is an EMBEDDEDV lake in its
# centroid-seeded pond cell — tune WITHIN this design (bedleak, table, surfdep, solver), never
# the connection type. (Excavated and on-top VERTICAL modes are retired dead ends: excavation
# leaves an unconditionable partial-saturation cell; a VERTICAL connection's bed sits at the
# CELL TOP so a dug-in pond would never seep — TM6-A55 p7-13/7-18.)
POND_LAYER  = "ponds_cdl"
POND_SUBSET = [4,5,11,15,17,18,0,1,2,3,6,7,8,9,10,12,13,14,16]   # all 19 ponds (FID order:
                      #   6 in-contact valley ponds first, then the perched upland ponds)
# ─── DO NOT EDIT THIS GUARD — it switches the ponds ON/OFF automatically. ─────────────────────
# "ss"        -> ponds STRIPPED (an SS run WITH lakes cannot converge: no storage to damp the
#                lake balance). "transient" -> ponds ACTIVE (the list above, as-is).
# To go from the IC run to the lake run change ONLY SPINUP_MODE — never this guard.
if SPINUP_MODE == "ss":
    POND_SUBSET = []
# --- CRR: cascade-routing & reinfiltration (Daoud et al. 2022, Eq. 23) -------------------------
# Rejected infiltration (UZF provider) + groundwater exfiltration (DRN-SEEP provider) are routed
# to the DOWNSLOPE NEIGHBOUR cells by MVR FACTORs α_ij = β·S_ij/ΣS_ij (S_ij = slope to neighbour j,
# negative slopes → 0). Receivers: the neighbour's LAK (pond), its SFR reach, or its land-surface
# UZF cell (REINFILTRATION — the leg standard MF6 lacks; 8.4% of P in Sardon). A cell with no
# downslope receiver is a SINK: its unrouted water leaves as evaporation (Daoud's convention).
# CRR_ENABLE=False restores the old blanket nearest-reach routing (for comparison runs).
CRR_ENABLE = True
CRR_BETA   = 1.0                 # β flow-partitioning factor (Daoud calibrates 0.8–1; 1 = route everything downslope)
POND_DEPTH  = 4.0     # m, pond depth below the local ground: lake bottom = mean(footprint DEM) − POND_DEPTH
                      #   (drives the stage/volume/area table + the equilibrium initial stage)
VORONOI_LAYERS = r"E:\00code_ws\DRYAD\CdL_model\conceptual\layers\voronoi_layers.npz"  # build_layers.py
LAK_BEDLEAK = 1e-3    # 1/d, lakebed leakance — the physical clay-lined-charca value (El-Zehairy 2018:
                      #   7e-4..1.5e-3/d; Turawa clay-lined 0.0007–0.0015/d). USER-CONFIRMED behaviour
                      #   2026-07-04: low seepage + never-dry ponds is the EXPECTED CdL physics.
                      #   THE calibration lever for the NbS recharge question.
LAK_SURFDEP = 0.5     # m, LAK option (TM6-A55): smooths the connection wetted area over [top, top+surfdep]
                      #   — "enhances convergence"; the lever that first stabilised the perched ponds.
LAK_MAXITER  = 200    # LAK-internal Newton cap (MF6 default 100): lets the stage settle in hard periods.
LAK_STAGECHG = 1e-4   # m, LAK stage closure tolerance (default 1e-5): 10x looser kills the tiny-tolerance
                      #   overshoot on small ponds; mass-balance-safe at 1e-4 m on a ~4 m pond.

# Reuse the cached Voronoi grid (WORKSPACE/voronoi_grid.pkl) instead of re-running
# the slow Triangle + Voronoi build every time.  Set True after you change the
# watershed / streams / CELL_* / STREAM_BUFFER — it will rebuild and re-cache.
REBUILD_GRID = False   # Design-2 grid cached at ncpl=6609 (40 m graded corridor -> 100 m background, centroid-seeded
                       # ponds; rebuilt 2026-07-20; layers npz at 6609). Set True only to rebuild the mesh — §5 then
                       # AUTO-runs build_layers.py on the ncpl mismatch (warning printed) and the run continues.

# Fast iteration when tuning solver / BC params (NOT UZF, forcing, grid, or the
# set of packages): WRITE_ALL=False skips the big UZF build+write (nper x ncpl)
# and rewrites ONLY the packages in WRITE_THESE, reusing every other input file
# already on disk.  Requires a prior full run (WRITE_ALL=True) so the rest exist.
WRITE_ALL   = True    # mode change (SS+lakes, STO steady_state, ATS periods, print_option) -> full write
# packages to refresh in fast mode; e.g. ['ims', 'npf', 'sto', 'ic', 'DRN-OUT', 'SFR', 'MVR', 'oc']
WRITE_THESE = ["lak"]   # (only used when WRITE_ALL=False) — LAK_BEDLEAK lowered, only connectiondata changed

# -----------------------------------------------------------------------------
# CRS HANDLING — detect & reproject all geographic inputs to ETRS89 / PT-TM06
# -----------------------------------------------------------------------------
from pyproj import CRS

TARGET_CRS = CRS.from_epsg(3763)   # ETRS89 / Portugal TM06 — metres

def report_crs(name, gdf_or_crs):
    """Print a one-line CRS summary for a GeoDataFrame or CRS object."""
    crs = gdf_or_crs.crs if hasattr(gdf_or_crs, "crs") else gdf_or_crs
    if crs is None:
        print(f"   [{name}] CRS = UNDEFINED  !!")
        return None
    crs = CRS.from_user_input(crs)
    unit = crs.axis_info[0].unit_name if crs.axis_info else "?"
    epsg = crs.to_epsg()
    print(f"   [{name}] CRS = {crs.name}  (EPSG:{epsg}, units: {unit})")
    return crs

def ensure_projected(gdf, name, target=TARGET_CRS):
    """
    Verify CRS is defined, report it, and reproject to `target` if needed.
    Returns the (possibly reprojected) GeoDataFrame.
    """
    src = report_crs(name, gdf)
    if src is None:
        raise RuntimeError(
            f"Layer '{name}' has no CRS defined. "
            f"Set it explicitly with gdf.set_crs(<epsg>, inplace=True) before reprojecting."
        )
    if src.to_epsg() == target.to_epsg():
        print(f"   [{name}] already in target CRS — no reprojection needed.")
        return gdf
    print(f"   [{name}] reprojecting -> EPSG:{target.to_epsg()} ({target.name})")
    return gdf.to_crs(target)

def check_raster_crs(path, name, target=TARGET_CRS):
    """Report a raster's CRS and warn if it differs from target."""
    with rasterio.open(path) as src:
        rcrs = src.crs
    report_crs(name, rcrs)
    if rcrs is None:
        print(f"   [{name}] !! raster has no CRS — assuming target.")
        return
    if CRS.from_user_input(rcrs).to_epsg() != target.to_epsg():
        print(
            f"   [{name}] !! raster CRS differs from target EPSG:{target.to_epsg()}.\n"
            f"      Sampling will be done in raster CRS — reproject points before sampling,\n"
            f"      OR pre-process the DEM with gdalwarp / rasterio.warp to EPSG:{target.to_epsg()}."
        )

# -----------------------------------------------------------------------------
# 1. LOAD & CLEAN GEOMETRIES  (with CRS detection + reprojection)
# -----------------------------------------------------------------------------
print(">> Loading & validating watershed polygon …")
ws = gpd.read_file(WATERSHED_GPKG, layer=WATERSHED_LAYER)
ws = ensure_projected(ws, "watershed")
ws["geometry"] = ws.geometry.apply(lambda g: make_valid(g) if not g.is_valid else g)
ws["geometry"] = ws.geometry.buffer(0)
ws_poly = unary_union(ws.geometry.values)
if isinstance(ws_poly, MultiPolygon):
    ws_poly = max(ws_poly.geoms, key=lambda p: p.area)
print(f"   watershed area = {ws_poly.area:,.0f} m²,  "
      f"perimeter = {ws_poly.length:,.0f} m")

print(">> Loading & validating stream network …")
streams = gpd.read_file(WATERSHED_GPKG, layer=STREAMS_LAYER)
streams = ensure_projected(streams, "streams", target=CRS.from_user_input(ws.crs))
streams["geometry"] = streams.geometry.apply(lambda g: make_valid(g))
streams = streams[streams.geometry.type.isin(["LineString", "MultiLineString"])]
merged = linemerge(unary_union(streams.geometry.values))
stream_lines = [merged] if merged.geom_type == "LineString" else list(merged.geoms)

def dedup(line):
    coords = list(dict.fromkeys(line.coords))
    return LineString(coords) if len(coords) > 1 else None
stream_lines = [dedup(l) for l in stream_lines if l is not None]
stream_lines = [l for l in stream_lines if l is not None]
stream_union = unary_union(stream_lines)

# DEM CRS check (the DEM is sampled later — warn if it's not in the same CRS)
print(">> Checking DEM CRS …")
check_raster_crs(DEM_TIF, "DEM")

# -----------------------------------------------------------------------------
# 2. FIND CATCHMENT OUTLET (stream ∩ watershed boundary)
# -----------------------------------------------------------------------------
print(">> Locating catchment outlet(s) …")
boundary = ws_poly.boundary
# Outlets = stream endpoints that reach the watershed boundary.  Do NOT use
# stream∩boundary: endpoints that only *graze* the boundary aren't returned by
# GEOS, which dropped the two low top-left outlets (the DRN then landed on the
# single high crossing instead of the real low ones).
OUTLET_BND_TOL = 5.0   # m: a stream endpoint within this of the boundary is an outlet
_eps = list(dict.fromkeys([p for ln in stream_lines for p in (ln.coords[0], ln.coords[-1])]))
_bd = np.array([Point(p).distance(boundary) for p in _eps])
outlet_xys = [_eps[i] for i in range(len(_eps)) if _bd[i] <= OUTLET_BND_TOL]
if not outlet_xys:
    outlet_xys = [_eps[int(np.argmin(_bd))]]          # fallback: closest endpoint
with rasterio.open(DEM_TIF) as src:
    _oe = [float(next(src.sample([p]))[0]) for p in outlet_xys]
outlet_xy = outlet_xys[int(np.argmin(_oe))]           # lowest outlet (routing hint + map marker)
print(f"   {len(outlet_xys)} boundary outlet(s); lowest @ "
      f"({outlet_xy[0]:.0f}, {outlet_xy[1]:.0f}) elev={min(_oe):.2f} m")

# -----------------------------------------------------------------------------
# 2b. LOAD CATTLE PONDS (charcas) — the subset to model with LAK
# -----------------------------------------------------------------------------
print(">> Loading cattle ponds (ponds_cdl) …")

def _largest_poly(g):
    """Return the largest Polygon part of a (Multi)Polygon."""
    if isinstance(g, MultiPolygon):
        return max(g.geoms, key=lambda p: p.area)
    return g

ponds_all = gpd.read_file(WATERSHED_GPKG, layer=POND_LAYER)
ponds_all = ensure_projected(ponds_all, POND_LAYER, target=CRS.from_user_input(ws.crs))
ponds_all["geometry"] = ponds_all.geometry.apply(lambda g: make_valid(g))

# ALL ponds clipped to watershed — used as Triangle constraints regardless of POND_SUBSET
pond_polys_all, pond_fids_all = [], []
for _fid in range(len(ponds_all)):
    _g = _largest_poly(ponds_all.geometry.iloc[_fid]).intersection(ws_poly)
    if not _g.is_empty:
        pond_polys_all.append(_largest_poly(_g))
        pond_fids_all.append(_fid)
print(f"   {len(pond_polys_all)} pond(s) found in watershed (all will be meshed into the grid).")

# POND_SUBSET ponds — active LAK lakes
pond_polys, pond_meta = [], []
for fid in POND_SUBSET:
    if fid < 0 or fid >= len(ponds_all):
        print(f"   !! pond FID {fid} out of range (0..{len(ponds_all) - 1}); skipped.")
        continue
    g = _largest_poly(ponds_all.geometry.iloc[fid]).intersection(ws_poly)
    if g.is_empty:
        print(f"   !! pond FID {fid} outside the watershed; skipped.")
        continue
    g = _largest_poly(g)
    d2s = g.distance(stream_union)
    pond_polys.append(g)
    pond_meta.append({"fid": fid, "area": g.area, "dist_stream": d2s, "online": d2s <= 1.0})
    print(f"   pond FID {fid}: area={g.area:6.0f} m²  dist_stream={d2s:5.1f} m  "
          f"{'ON-channel' if d2s <= 1.0 else 'off-channel'}")
print(f"   {len(pond_polys)} pond(s) selected for LAK.")

# -----------------------------------------------------------------------------
# 3. BUILD VORONOI GRID (refined along streams)  -- ROBUST VERSION
# -----------------------------------------------------------------------------
print(">> Triangulating domain with refinement along streams …")

tri_ws = WORKSPACE / "tri"
tri_ws.mkdir(exist_ok=True)
# clean any previous Triangle run
for f in tri_ws.glob("_triangle.*"):
    f.unlink()

# ---- 3a. Prepare the outer watershed polygon
# Densify the outer boundary so Triangle can satisfy the area constraint
# along the boundary without inserting huge skinny triangles
ws_dense = ws_poly.exterior
# resample at ~CELL_FAR spacing
n_outer = max(int(np.ceil(ws_dense.length / CELL_FAR)), 50)
domain_xy = [ws_dense.interpolate(i / n_outer, normalized=True).coords[0]
             for i in range(n_outer)]

# ---- 3b. Prepare stream refinement zones
# Sanity-check CRS units before buffering
crs_units = ws.crs.axis_info[0].unit_name.lower() if ws.crs and ws.crs.axis_info else "unknown"
if "metre" not in crs_units and "meter" not in crs_units:
    raise RuntimeError(
        f"Watershed CRS units are '{crs_units}', not metres. "
        f"Reproject to a projected CRS (e.g. EPSG:3763 for PT-TM06/ETRS89) before buffering."
    )

print(f"   watershed area = {ws_poly.area:,.0f} m²  "
      f"(perimeter = {ws_poly.length:,.0f} m)")

# Shrink inward by STREAM_BUFFER (not 2x — that was overly conservative)
ws_inner = ws_poly.buffer(-STREAM_BUFFER)

# Handle MultiPolygon / GeometryCollection results
if ws_inner.is_empty:
    raise RuntimeError(
        f"Inner watershed is empty after buffer(-{STREAM_BUFFER}). "
        f"Watershed area={ws_poly.area:.1f} m². Reduce STREAM_BUFFER."
    )
if isinstance(ws_inner, MultiPolygon):
    ws_inner = max(ws_inner.geoms, key=lambda p: p.area)
elif ws_inner.geom_type == "GeometryCollection":
    polys = [g for g in ws_inner.geoms if g.geom_type in ("Polygon", "MultiPolygon")]
    if not polys:
        raise RuntimeError("Inner watershed collapsed to non-polygonal geometry.")
    ws_inner = max(polys, key=lambda p: p.area)

# clip stream union to the inner watershed
streams_clipped = stream_union.intersection(ws_inner)
if streams_clipped.is_empty:
    # fallback: clip to the watershed itself, then shrink the buffer instead
    print("   !! No streams in inner watershed — using full watershed and shrinking buffer.")
    streams_clipped = stream_union.intersection(ws_poly)
    if streams_clipped.is_empty:
        raise RuntimeError("Streams do not intersect the watershed at all — check CRS / layers.")

# single dissolved buffer => no overlap, no self-intersect
stream_buf = streams_clipped.buffer(STREAM_BUFFER, cap_style=1, join_style=1)  # round: mitre joins spiked and broke Triangle
stream_buf = unary_union(stream_buf)
if isinstance(stream_buf, MultiPolygon):
    stream_buf_polys = list(stream_buf.geoms)
else:
    stream_buf_polys = [stream_buf]

# ---- 3c. Build Triangle PSLG
# POND MESHING = CENTROID SEEDING (user 2026-07-04, "one cell per pond"): one Voronoi generator
# node at each pond centroid + a helper ring -> the centroid's cell is pond-SIZED, pond-CENTRED
# and COVERS the footprint. (Rim-polygon constraints CANNOT do this: Voronoi faces bisect the rim,
# so rim-based cells always straddle it half-in/half-out — proven on the 5922 rim-mesh attempt.)
# Everything else (streams, smooth Daoud transitions) is untouched.
POND_RING_FACTOR = 2.5   # ring nodes at this × pond-radius -> the centroid's Voronoi cell shrinks to ~pond size
POND_RING_N      = 6     #   (without the ring the centroid cell takes the local background size)
_seed = []
for _pg in pond_polys_all:
    _c = _pg.representative_point(); _r = float(np.sqrt(_pg.area / np.pi))
    _seed.append((_c.x, _c.y))
    for _k in range(POND_RING_N):
        _th = 2.0 * np.pi * _k / POND_RING_N
        _px, _py = _c.x + POND_RING_FACTOR * _r * np.cos(_th), _c.y + POND_RING_FACTOR * _r * np.sin(_th)
        if ws_poly.contains(Point(_px, _py)):                     # keep ring nodes inside the watershed
            _seed.append((_px, _py))
_pond_seed_pts = np.array(_seed) if _seed else None
print(f"   pond seeds: {len(pond_polys_all)} centroid(s) + rings ({POND_RING_N}@{POND_RING_FACTOR}×r) = {len(_seed)} nodes")
tri = Triangle(
    model_ws=str(tri_ws),
    angle=20,   # 30 sits at Triangle's stability edge -> "ran out of precision"
    exe_name=TRIANGLE_EXE,
    nodes=_pond_seed_pts,     # one generator node per pond -> one pond-scale Voronoi cell each
    additional_args=["-j"],   # remove unused vertices; helps stability
)

# Outer polygon = watershed boundary, with coarse target area
tri.add_polygon(domain_xy)

# Coarse region marker: a point well inside the watershed but OUTSIDE any
# stream buffer.  Use the watershed centroid; if it falls inside a buffer,
# walk outward until it doesn't.
coarse_pt = ws_poly.representative_point()
if any(p.contains(coarse_pt) for p in stream_buf_polys):
    # find a safer interior point: difference of watershed and stream buffers
    safe_area = ws_poly.difference(unary_union(stream_buf_polys))
    if safe_area.is_empty:
        raise RuntimeError("Entire watershed is inside stream buffer.")
    if safe_area.geom_type == "MultiPolygon":
        safe_area = max(safe_area.geoms, key=lambda p: p.area)
    coarse_pt = safe_area.representative_point()

tri.add_region((coarse_pt.x, coarse_pt.y),
               attribute=0, maximum_area=CELL_FAR ** 2)

# Stream corridor at CELL_NEAR_STREAM (one region per disjoint buffer piece) — added only
# when NOT smoothing; the smooth sizing field below puts the 40 m corridor at level 40 and
# grades it outward, so adding it here too would double the boundary.
if STREAM_REFINE and not SMOOTH_TRANSITIONS:
    for buf in stream_buf_polys:
        if buf.is_empty or not buf.is_valid:
            continue
        ring = buf.exterior
        n_ring = max(int(np.ceil(ring.length / CELL_NEAR_STREAM)), 8)
        coords = [ring.interpolate(i / n_ring, normalized=True).coords[0] for i in range(n_ring)]
        tri.add_polygon(coords)
        rp = buf.representative_point()
        tri.add_region((rp.x, rp.y), attribute=1, maximum_area=CELL_NEAR_STREAM ** 2)

# Pond refinement: NONE here — ponds are meshed by the centroid seeding above (one pond-scale cell
# each). No rim polygons (Voronoi cells cannot follow them) and no interior area constraints.
# The smooth sizing field below still grades outward from POND_RIM_SIZE around the footprints.

# Smooth transition zones (after Daoud et al. 2022) — CONTINUOUS sizing field so the cell
# size grades smoothly from the 2.5 m ponds / 40 m stream corridor out to CELL_FAR, honouring
# the CVFD requirement.  Two things make Triangle accept the PSLG (both verified in
# preview_transition_grid.py): (1) the buffer distance per level grows CONTINUOUSLY -> the
# level sets STRICTLY nest, so no two level boundaries coincide; (2) each feature is clipped
# BEFORE buffering (_clip_then_buffer) so its zone reaches the watershed edge only tangentially
# (no shared clip edge between levels).  The pond cores (2.5 m) are added by the pond loop above.
POND_RIM_SIZE = 10.0   # m, base size of the pond transition grading (~the seeded pond-cell scale)
if SMOOTH_TRANSITIONS:
    _ponds_u = unary_union(pond_polys_all)   # coarse pond rims (loop above); the sizing field grades them outward
    def _clip_then_buffer(geom, B):
        inside = ws_poly.buffer(-B)
        if inside.is_empty:
            return None
        g = geom.intersection(inside)
        return None if g.is_empty else g.buffer(B)
    _prev = _ponds_u
    for _L in TRANS_LEVELS:
        _parts = []
        if _ponds_u is not None and _L > POND_RIM_SIZE:   # strict >: a zero-distance buffer would re-add
            _gp = _clip_then_buffer(_ponds_u, (_L - POND_RIM_SIZE) / TRANS_GRADE)   # the rim segments coincidentally
            if _gp is not None:                                                    # (Triangle "topological inconsistency")
                _parts.append(_gp)
        if STREAM_REFINE and _L >= CELL_NEAR_STREAM:      # SFRmaker mode: NO stream-corridor refinement
            _gs = _clip_then_buffer(streams_clipped, STREAM_BUFFER + (_L - CELL_NEAR_STREAM) / TRANS_GRADE)
            if _gs is not None:
                _parts.append(_gs)
        if not _parts:
            continue
        _G = unary_union(_parts)
        if _G.is_empty:
            continue
        for _poly in (_G.geoms if isinstance(_G, MultiPolygon) else [_G]):
            for _ring in [_poly.exterior] + list(_poly.interiors):
                _n = max(int(np.ceil(_ring.length / _L)), 8)
                tri.add_polygon([_ring.interpolate(i / _n, normalized=True).coords[0] for i in range(_n)])
        _annulus = _G.difference(_prev) if _prev is not None else _G
        for _ap in (_annulus.geoms if isinstance(_annulus, MultiPolygon) else [_annulus]):
            if _ap.is_empty or _ap.area < _L ** 2:
                continue
            _sp = _ap.representative_point()
            tri.add_region((_sp.x, _sp.y), attribute=int(_L), maximum_area=_L ** 2)
        _prev = _G
    print(f"   smooth transition zones: levels {TRANS_LEVELS} m (GRADE={TRANS_GRADE})")

# ---- Sliver cleanup: merge degenerate cells into a neighbour (user 2026-06-29) ----
SLIVER_MERGE_AREA = 1.0   # m², a cell smaller than this is a degenerate sliver -> mass-conservatively DISSOLVED into the
#   non-sliver neighbour sharing the longest edge (better than idomain deactivation: no domain holes). Matches the
#   validated criterion width < SLIVER_WMIN=1 m ⇒ area < 1 m²; §5's sliver_set is the residual guard (normally empty).
def _merge_slivers(gp_disv, min_area):
    """Dissolve each sliver-cluster (cells with area < min_area) into the adjacent non-sliver cell sharing the longest edge,
    then rebuild the DISV gridprops (flopy CVFD utils resolve shared/hanging vertices)."""
    from flopy.utils.cvfdutil import to_cvfd, get_disv_gridprops
    _vxy = {v[0]: (v[1], v[2]) for v in gp_disv["vertices"]}
    _polys = [Polygon([_vxy[iv] for iv in c[4:4 + c[3]]]).buffer(0) for c in gp_disv["cell2d"]]
    _n = len(_polys); _ar = np.array([p.area for p in _polys]); _sl = set(np.where(_ar < min_area)[0])
    if not _sl:
        print("   sliver-merge: no slivers found"); return gp_disv
    _tree = STRtree(_polys)
    def _nb(i):
        out = []
        for j in _tree.query(_polys[i]):
            j = int(j)
            if j != i and getattr(_polys[i].boundary.intersection(_polys[j].boundary), "length", 0.0) > 1e-6:
                out.append((j, _polys[i].boundary.intersection(_polys[j].boundary).length))
        return out
    _seen, _clusters = set(), []                                  # connected components of slivers (sliver-sliver adjacency)
    for s in _sl:
        if s in _seen:
            continue
        _st, _grp = [s], []
        while _st:
            u = _st.pop()
            if u in _seen:
                continue
            _seen.add(u); _grp.append(u)
            for v, _ in _nb(u):
                if v in _sl and v not in _seen:
                    _st.append(v)
        _clusters.append(_grp)
    _assign = {}                                                  # each cluster -> non-sliver neighbour with the longest shared edge
    for _grp in _clusters:
        _best, _bl = None, 0.0
        for m in _grp:
            for j, L in _nb(m):
                if j not in _sl and L > _bl:
                    _bl, _best = L, j
        if _best is not None:
            _assign.setdefault(_best, []).extend(_grp)
    _merged = {s for sl in _assign.values() for s in sl}
    _new = [(unary_union([_polys[i]] + [_polys[s] for s in _assign[i]]).buffer(0) if i in _assign else _polys[i])
            for i in range(_n) if i not in _merged]
    _vd = {k: list((p if p.geom_type == "Polygon" else max(p.geoms, key=lambda g: g.area)).exterior.coords)
           for k, p in enumerate(_new)}
    _verts, _iverts = to_cvfd(_vd, skip_hanging_node_check=True)
    _gp = get_disv_gridprops(np.array(_verts), _iverts)
    print(f"   sliver-merge: {len(_merged)} cell(s) < {min_area} m² in {len(_clusters)} cluster(s) -> "
          f"merged into {len(_assign)} neighbour(s); ncpl {_n} -> {_gp['ncpl']}")
    return _gp

# ---- Cache the slow Triangle + Voronoi build (gate behind REBUILD_GRID) ----
GRID_CACHE = WORKSPACE / "voronoi_grid.pkl"
GRID_TAG   = WORKSPACE / "voronoi_grid.design.txt"   # records which GRID_DESIGN built the cache
if (not REBUILD_GRID) and GRID_CACHE.exists():
    # FAIL FAST on a design/cache mismatch (2026-07-20: the user set GRID_DESIGN=2 and ran, but the
    # cached design-1 grid won silently — GRID_DESIGN only affects the TRIANGLE BUILD, so switching
    # designs REQUIRES a rebuild; without this check the run proceeds on the wrong grid unnoticed).
    _cached_design = int(GRID_TAG.read_text().strip()) if GRID_TAG.exists() else None
    if _cached_design is not None and _cached_design != GRID_DESIGN:
        raise RuntimeError(
            f"CONFIG ERROR: GRID_DESIGN={GRID_DESIGN} but the cached grid ({GRID_CACHE.name}) was built "
            f"with design {_cached_design}. To switch designs: set REBUILD_GRID=True and run with "
            "SPINUP_MODE='ss' (ONE run rebuilds the grid, auto-runs build_layers.py and regenerates "
            "the IC), then set REBUILD_GRID=False and run the transient.")
    print(f"   loading cached Voronoi grid ({GRID_CACHE.name}, design {_cached_design}; "
          f"set REBUILD_GRID=True to rebuild) …")
    with open(GRID_CACHE, "rb") as _f:
        gridprops_vg, gridprops_disv = pickle.load(_f)
else:
    # Run Triangle WITH verbose output so we see any future error immediately
    print("   running triangle.exe …")
    try:
        tri.build(verbose=True)
    except Exception as e:
        # dump the .poly file path so the user can inspect it
        poly_file = tri_ws / "_triangle.poly"
        print(f"!! Triangle failed.  Inspect the PSLG here: {poly_file}")
        print(e)
        raise

    print(">> Building Voronoi tessellation …")
    vor = VoronoiGrid(tri)
    gridprops_disv = vor.get_disv_gridprops()         # for ModflowGwfdisv (has 'nvert')
    gridprops_disv = _merge_slivers(gridprops_disv, SLIVER_MERGE_AREA)   # dissolve degenerate transition-ring slivers
    # VertexGrid gridprops derived from the (sliver-cleaned) disv gridprops — same vertices/cell2d, ncpl
    gridprops_vg = {"vertices": gridprops_disv["vertices"], "cell2d": gridprops_disv["cell2d"], "ncpl": gridprops_disv["ncpl"]}

    with open(GRID_CACHE, "wb") as _f:
        pickle.dump((gridprops_vg, gridprops_disv), _f)
    GRID_TAG.write_text(str(GRID_DESIGN))            # tag the cache with the design that built it
    print(f"   cached Voronoi grid -> {GRID_CACHE.name} (design {GRID_DESIGN})")

vgrid = VertexGrid(**gridprops_vg, nlay=1)
ncpl = vgrid.ncpl
print(f"   Voronoi cells per layer: {ncpl}")

xc = np.array([vgrid.xcellcenters[i] for i in range(ncpl)])
yc = np.array([vgrid.ycellcenters[i] for i in range(ncpl)])

# --- MESH QUALITY: flag degenerate "sliver" cells (Voronoi duals of near-degenerate triangles) -------------
# The Daoud smooth-transition rebuild creates a cluster of ~zero-area cells (min width ~0.03 m at x=-57732) that
# breaks GW convergence — a limit cycle stuck just above outer_dvclose at the cluster.  Flag cells thinner than
# SLIVER_WMIN; §5 deactivates them (idomain=0 in every layer) so they leave the solve.  They are <1 m and away
# from ponds/streams (transition-ring geometry), so the domain change is negligible.  (User 2026-06-27: "fix
# slivers only" — keep the nlay=5 thin-CONN scheme so it can be judged on a clean grid.)
SLIVER_WMIN = 1.0   # m, equivalent width sqrt(area) below which a Voronoi cell is treated as a degenerate sliver
_cell_area = np.empty(ncpl)
for _i in range(ncpl):
    _cv = np.asarray(vgrid.get_cell_vertices(_i)); _cx, _cy = _cv[:, 0], _cv[:, 1]
    _cell_area[_i] = 0.5 * abs(np.dot(_cx, np.roll(_cy, 1)) - np.dot(_cy, np.roll(_cx, 1)))
sliver_set = set(np.where(np.sqrt(_cell_area) < SLIVER_WMIN)[0].tolist())
print(f"   mesh quality: min cell width {np.sqrt(_cell_area).min():.3f} m; "
      f"deactivating {len(sliver_set)} sliver cell(s) < {SLIVER_WMIN} m (degenerate transition-ring cells)")

# Build GridIntersect robustly across FloPy versions
try:
    ix = GridIntersect(vgrid, method="vertex")
except TypeError:
    # older FloPy: no `method` kwarg
    try:
        ix = GridIntersect(vgrid, rtree=True)
    except TypeError:
        ix = GridIntersect(vgrid)

# -----------------------------------------------------------------------------
# 4. CLIP DEM ONTO VORONOI POLYGONS  (mean elevation per cell)
# -----------------------------------------------------------------------------
print(">> Sampling DEM at Voronoi cell centroids …")
with rasterio.open(DEM_TIF) as src:
    samples = list(src.sample(list(zip(xc, yc))))
    top_elev = np.array([s[0] for s in samples], dtype=float)
    nodata = src.nodata
if nodata is not None:
    top_elev[top_elev == nodata] = np.nan
# fill any nan with mean
top_elev = np.where(np.isnan(top_elev), np.nanmean(top_elev), top_elev)

# -----------------------------------------------------------------------------
# 4b. SAMPLE SOIL-HYDRAULIC RASTERS -> per-cell UZF parameters
# -----------------------------------------------------------------------------
print(">> Sampling soil-hydraulic rasters at cell centroids …")

def _sample_raster(path, xy, scale=1.0):
    """Sample a raster at (x,y) points; nodata -> mean; multiply by `scale`."""
    with rasterio.open(path) as src:
        vals = np.array([v[0] for v in src.sample(xy)], dtype=float)
        nd = src.nodata
    if nd is not None:
        vals[vals == nd] = np.nan
    vals[np.abs(vals) > 1e30] = np.nan
    vals = vals * scale
    return np.where(np.isnan(vals), np.nanmean(vals), vals)

_xy = list(zip(xc, yc))
vks_cell  = _sample_raster(SOIL_VKS_TIF,  _xy, SOIL_VKS_TO_MD)     # m/d
thts_cell = _sample_raster(SOIL_THTS_TIF, _xy, SOIL_WC_TO_FRAC)    # -
thtr_cell = _sample_raster(SOIL_THTR_TIF, _xy, SOIL_WC_TO_FRAC)    # -
thti_cell = _sample_raster(SOIL_THTI_TIF, _xy, SOIL_WC_TO_FRAC)    # -
# enforce MF6 UZF ordering thtr < thti < thts (clip noisy/edge pixels)
thts_cell = np.clip(thts_cell, 0.05, 0.95)
thtr_cell = np.clip(thtr_cell, 0.01, thts_cell - 0.02)
thti_cell = np.clip(thti_cell, thtr_cell + 0.01, thts_cell - 0.01)
print(f"   vks  {vks_cell.min():.3g}–{vks_cell.max():.3g} m/d (mean {vks_cell.mean():.3g})  "
      f"[ks cm/d × {SOIL_VKS_TO_MD}]")
print(f"   thts {thts_cell.min():.3f}–{thts_cell.max():.3f}  |  "
      f"thti {thti_cell.min():.3f}–{thti_cell.max():.3f}  |  "
      f"thtr {thtr_cell.min():.3f}–{thtr_cell.max():.3f}")

# -----------------------------------------------------------------------------
# 5. LAYER GEOMETRY  (nlay=3 — ONE numerical layer per geologic unit, Daoud-style)
#    L0 = U1 alluvium (K=25) / L1 = U2 terraces (K=40) / L2 = U3 Plio-Miocene (K=1.5, base -35 m).
#    Unit bottoms from build_layers.py -> voronoi_layers.npz (top + bU1/bU2/bU3 at these DISV
#    centroids). Where a unit is ABSENT (presence tested against the npz's OWN top — never the
#    mixed-source model-DEM thickness, the phantom-U1 lesson) the layer pinches to idomain=-1
#    vertical pass-through, so the TOP ACTIVE layer at each cell IS the outcropping unit.
#    Present-but-thin cells are extended down to MIN_ACTIVE_THK=4 m (feather-edge conditioning).
#    The topmost validated fit: full 45-yr no-lake RMSE 2.5 m / bias +0.5 m.
# -----------------------------------------------------------------------------
_vl = np.load(VORONOI_LAYERS)
if int(_vl["ncpl"]) != ncpl:
    # Layers npz from a DIFFERENT grid (normal right after a grid rebuild) -> AUTO-HEAL
    # (user 2026-07-20): warn, run build_layers.py for this grid, reload, continue.
    print(f"\n   !! voronoi_layers.npz is for a {int(_vl['ncpl'])}-cell grid but this grid has "
          f"{ncpl} cells\n   !! -> auto-running build_layers.py to regenerate it (takes a minute) …")
    import subprocess, sys
    _bl = Path(__file__).resolve().parent / "build_layers.py"
    _rc = subprocess.run([sys.executable, str(_bl)],
                         env=dict(os.environ, MPLBACKEND="Agg"), check=False).returncode
    if _rc != 0:
        raise RuntimeError(f"build_layers.py failed (exit {_rc}) — fix it and re-run.")
    _vl = np.load(VORONOI_LAYERS)
    if int(_vl["ncpl"]) != ncpl:
        raise RuntimeError(
            f"build_layers.py wrote ncpl={int(_vl['ncpl'])} but the grid has {ncpl} cells — "
            "it read a different voronoi_grid.pkl? Check WORKSPACE paths in both scripts.")
    print(f"   auto-heal OK: voronoi_layers.npz regenerated at ncpl={ncpl}\n")
ub = np.asarray(_vl["botm"], float)                        # [bot_U1, bot_U2, bot_U3(-35)] at the DISV centroids
AQ_BASE = -35.0                                            # m a.s.l., formation/model base (matches build_layers.py U3 base)
MIN_FLOOR = 0.1                                            # m, min layer thickness (keeps botm strictly decreasing)
MIN_THK   = 0.5                                            # m, a unit-layer thinner than this -> idomain=-1 vertical PASS-THROUGH
KH_ALLU, KH_TERR, KH_FORM = 25.0, 40.0, 1.5              # m/d (U1 alluvium / U2 terraces repr.16-77 / U3 formation 1-2) — the CALIBRATED
                                                          #   values (no-lake SS RMSE ~2.1-2.6 m). (The x0.1 reduction was a lake-convergence
                                                          #   experiment — it lowered the LAK blow-up but failed AND pushed heads +5 m; reverted.)
nlay = 3                                                  # U1 + U2 + U3 — one numerical layer per geologic unit
top = top_elev.copy()
# FEATHER-EDGE FIX (2026-07-03, first nlay=3 run): where a unit is PRESENT but thinner than
# MIN_ACTIVE_THK, EXTEND its base down to MIN_ACTIVE_THK.  Without this, the unit pinch-out
# boundaries leave 0.5-1.3 m thick HIGH-K convertible TOP cells carrying DRN-SEEP (cond 1e4)
# + UZF (surfdep 0.5 ~ the whole cell) -> drown/dry limit cycle, SP5 non-convergence (the six
# failing cells all had top-active thk 0.58-1.34 m).  4 m = the proven HOST-layer thickness of
# the nlay=4/5 scheme, so feather edges locally reproduce the configuration that converged 45 yr.
MIN_ACTIVE_THK = 4.0
# UNIT PRESENCE from build_layers' OWN sampling (npz top vs npz bottoms — internally consistent).
# GOTCHA (2026-07-03, "phantom U1"): the model samples the DEM itself (top_elev) while the npz has
# its own top; they differ by ±0.3 m. Testing presence with the MIXED thickness (top_elev − ub)
# faked 0.5–4 m of U1 on 60 cells in the middle of the U2/U3 outcrop, which the feather-edge
# extension then activated as 4 m of K=25. Presence MUST come from one consistent source.
_top_npz  = np.asarray(_vl["top"], float)
_present1 = (_top_npz - ub[0]) > MIN_THK                              # U1 really exists here
_present2 = (ub[0] - ub[1])    > MIN_THK                              # U2 really exists here
botm = np.zeros((nlay, ncpl))
botm[0] = np.where(_present1, np.minimum(ub[0], top - MIN_FLOOR), top - MIN_FLOOR)   # absent -> pinched
_t0 = top - botm[0]
_grow0 = _present1 & (_t0 < MIN_ACTIVE_THK)                           # feather edge of a REAL unit -> extend to 4 m
botm[0] = np.where(_grow0, top - MIN_ACTIVE_THK, botm[0])
botm[1] = np.where(_present2, np.minimum(ub[1], botm[0] - MIN_FLOOR), botm[0] - MIN_FLOOR)
_t1 = botm[0] - botm[1]
_grow1 = _present2 & (_t1 < MIN_ACTIVE_THK)
botm[1] = np.where(_grow1, botm[0] - MIN_ACTIVE_THK, botm[1])
botm[2] = np.minimum(np.minimum(ub[2], AQ_BASE), botm[1] - MIN_FLOOR) # U3 layer base (flat -35)
thk = np.vstack([top - botm[0], botm[0] - botm[1], botm[1] - botm[2]])
print(f"   unit presence (npz-consistent): U1 {int(_present1.sum())}, U2 {int(_present2.sum())} of {ncpl}; "
      f"feather-edge floor: U1 {int(_grow0.sum())} + U2 {int(_grow1.sum())} thin cell(s) extended to {MIN_ACTIVE_THK} m")

# K = the unit K per layer (no point-in-polygon HOST needed: the TOP ACTIVE layer at each cell IS the outcropping unit, via
# the pinch-out below). The _hydro sjoin is kept only for the alluvium/terraces percentages in the summary + parameter plots.
_hydro = gpd.read_file(WATERSHED_GPKG, layer="dryad_modelo_nbs__gc_35a_cdl_hydrostrat").to_crs(TARGET_CRS)
_hydro["Codigo"] = _hydro["Codigo"].astype(str)
_cellpts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(xc, yc), crs=TARGET_CRS)
_sj = gpd.sjoin(_cellpts, _hydro[["Codigo", "geometry"]], how="left", predicate="within")
_code = _sj[~_sj.index.duplicated(keep="first")].reindex(range(ncpl))["Codigo"].fillna("3").to_numpy()
allu = _code == "1"; terr = _code == "2"
k_array = np.vstack([np.full(ncpl, KH_ALLU), np.full(ncpl, KH_TERR), np.full(ncpl, KH_FORM)])  # U1=25 / U2=40 / U3=1.5
k33 = 0.1 * k_array                                        # Kh/Kv = 10
icelltype = np.ones((nlay, ncpl), dtype=int)               # all convertible

# base idomain: U3 (basal, L2) always active; U1 (L0) and U2 (L1) active where the unit is substantial, else -1 PASS-THROUGH
# -> the topmost ACTIVE layer at each cell is the outcropping unit (the embedded lake connects to it, §5b).
idomain_unit = np.ones((nlay, ncpl), dtype=int)
for L in (0, 1):
    idomain_unit[L] = np.where(thk[L] > MIN_THK, 1, -1)
for _s in sliver_set:                                      # (sliver_set is empty now — slivers are MERGED at grid build — kept as a guard)
    idomain_unit[:, _s] = 0
print(f"   nlay={nlay} = U1 + U2 + U3 (one layer per unit; no HOST/CONN); U1 active {int((idomain_unit[0]==1).sum())}, "
      f"U2 active {int((idomain_unit[1]==1).sum())} of {ncpl} (rest -1 pass-through); "
      f"alluvium {100*np.mean(allu):.0f}%/terraces {100*np.mean(terr):.0f}% outcrop; base {AQ_BASE}")

# -----------------------------------------------------------------------------
# 5b. POND -> LAKE CELLS  (EMBEDDEDV geometry, built by hand)
#     flopy's get_lak_connections() does not support embedded lakes on a
#     VertexGrid, so the geometry is explicit: each pond = ONE EMBEDDEDV
#     connection in its centroid-seeded cell (§3), lake bottom = mean footprint
#     DEM − POND_DEPTH, bathymetry from the 4-col wedge table (§12b), initial
#     stage = the SS water table clamped to [bottom+0.1, ground].
# -----------------------------------------------------------------------------
def pond_cells_of(pg):
    """Cells OWNED by a pond = cells whose CENTER lies inside the footprint (majority-inside
    proxy — the raw intersect set also catches rim-straddling neighbours), with a nearest-cell
    fallback for a pond smaller than its (centroid-seeded) cell.  Used by §5b (LAK cells) and
    by the §14 pond-zoom pages (planned cells in an SS build)."""
    try:
        _res = ix.intersect(pg, geo_dataframe=False)
    except TypeError:
        _res = ix.intersect(pg)
    _cand = sorted({int(c) for c in _res["cellids"]})
    _nds = [nd for nd in _cand if pg.contains(Point(float(xc[nd]), float(yc[nd])))]
    if not _nds:
        _c = pg.representative_point()
        _nds = [int(_cand[int(np.argmin([(xc[nd] - _c.x) ** 2 + (yc[nd] - _c.y) ** 2 for nd in _cand]))])] \
            if _cand else [int(np.argmin((xc - _c.x) ** 2 + (yc - _c.y) ** 2))]
    return _nds

pond_nodes = []
nlakes = len(pond_polys)
idomain_lak = idomain_unit.copy()    # = the unit pattern (1 active / -1 pass-through). NO lake mode excavates anymore.

# TOP ACTIVE layer per cell = the land-surface cell of each column (unit-layer scheme: the layer of
# the OUTCROPPING unit; -1 pass-through cells skipped).  EVERY surface package (DRN-SEEP cellids,
# UZF landflag/surfdep, UZF finf/pet forcing, UZF->SFR runoff MVR) AND the LAK bed connections MUST
# key on this — NOT on layer 0: with nlay=3, layer 0 = U1 alluvium exists on only ~5% of the
# catchment (the old nlay=4/5 HOST layer was active everywhere, whence the layer-0 shortcuts).
top_act = np.full(ncpl, -1, dtype=int)
for _lay in range(nlay):
    _m = (idomain_lak[_lay] == 1) & (top_act < 0)
    top_act[_m] = _lay
print("   top-active (land-surface) layer: " +
      " | ".join(f"L{k}: {int((top_act == k).sum())}" for k in range(nlay)) +
      f" | inactive: {int((top_act < 0).sum())}")
lak_conn, lak_pkg, lak_strt, lak_geom = [], [], [], []   # lak_geom: (lbot, gmean, area) per lake -> stage-vol-area table
if nlakes:
    # FAIL FAST on the two configurations that can never converge (both actually happened —
    # a 13 s guaranteed failure each; make them clear errors instead):
    if SPINUP_MODE == "ss":
        raise RuntimeError(
            f"CONFIG ERROR: SPINUP_MODE='ss' but {nlakes} lake(s) are active. An SS run WITH lakes "
            "cannot converge (no storage to damp the lake balance). The POND_SUBSET guard (§ flags) "
            "must read `if SPINUP_MODE == \"ss\": POND_SUBSET = []` — restore it "
            "and flip ONLY SPINUP_MODE to switch run types.")
    print(">> Mapping ponds to lake cells (manual LAK geometry) …")
    # SS water table per cell (from the IC) — used to start each lake at its LOCAL
    # equilibrium stage (see the stage clamp below), so neither perched nor in-contact
    # ponds get a startup shock at TS=1.
    _pond_wt = None
    if INIT_HEADS_FILE is not None and Path(INIT_HEADS_FILE).exists():
        _hic = np.load(INIT_HEADS_FILE)
        if _hic.size == nlay * ncpl:
            _hic = np.where(np.abs(_hic.reshape(nlay, ncpl)) < 1e29, _hic.reshape(nlay, ncpl), np.nan)
            _pond_wt = np.nanmax(_hic, axis=0)        # water table = highest saturated head in each column
    if _pond_wt is None:
        raise RuntimeError(
            f"CONFIG ERROR: a transient LAKE run needs the SS water table for the equilibrium lake "
            f"stages, but {SPINUP_HEADS_OUT.name} is missing or stale (needs size {nlay}x{ncpl}). "
            "Run once with SPINUP_MODE='ss' first to regenerate it — starting the lakes blind at "
            "their spill inverts diverges at SP1 TS1.")
    for L, pg in enumerate(pond_polys):
        nodes = pond_cells_of(pg)                       # LAK cells = centroid-inside (+ nearest fallback)
        pond_nodes.append(nodes)
        gmean = float(np.mean(top[nodes]))
        lbot = gmean - POND_DEPTH
        _wt = float(np.nanmean(_pond_wt[nodes]))       # guaranteed loaded by the fail-fast above
        _perched = _wt < lbot                          # (reporting only; the stage clamp handles both)
        # EMBEDDEDV in ONE host cell (the footprint cell nearest the centroid — usually THE seeded
        # pond cell). The host stays the normal UZF surface cell that converges in the no-lake
        # model. NLAKECONN=1; belev/telev ignored for EMBEDDEDV; connlen/connwidth > 0 required.
        # The lake's real geometry (incl. the dug-in bottom at lbot) lives in its 4-col table.
        _cx, _cy = float(np.mean(xc[nodes])), float(np.mean(yc[nodes]))
        _rep = int(nodes[int(np.argmin((xc[nodes] - _cx) ** 2 + (yc[nodes] - _cy) ** 2))])
        _klay = 0                                   # connect at the TOP ACTIVE layer = the outcropping unit
        while _klay < nlay - 1 and idomain_lak[_klay, _rep] != 1:
            _klay += 1
        _clen = 0.5 * POND_DEPTH                    # m, node->lake distance (≈ half the pond depth)
        _cwid = float(np.sqrt(max(pond_meta[L]["area"], 1.0)))   # m, representative connection width (>0 required)
        lak_conn.append([L, 0, (_klay, _rep), "EMBEDDEDV", LAK_BEDLEAK, 0.0, 0.0, _clen, _cwid])  # iconn 0-based!
        # Start each lake at its LOCAL EQUILIBRIUM stage = clip(SS water table, lake bottom+0.1, ground): perched ponds start
        # ~empty, in-contact start ~full — avoids a TS=1 startup shock.
        strt = float(np.clip(_wt, lbot + 0.1, gmean))
        lak_strt.append(strt)
        lak_geom.append((lbot, gmean, float(pond_meta[L]["area"])))   # for the stage-volume-(area+barea) table
        lak_pkg.append([L, strt, 1, f"pond{pond_meta[L]['fid']}"])
        print(f"   lake {L} (FID {pond_meta[L]['fid']}): {len(nodes)} footprint cell(s) -> EMBEDDEDV in 1 cell; "
              f"bottom {lbot:.2f} m, WT {_wt:.2f} m, init stage {strt:.2f} m"
              f"  [{'PERCHED' if _perched else 'contact'}, {'ON' if pond_meta[L]['online'] else 'off'}-channel]")
    print(f"   {nlakes} lake(s); {len(lak_conn)} connection(s); EMBEDDEDV (no excavation)")

# node -> lake index: drives the SFR excision + the SFR↔LAK MVR hand-offs (the stream must pass
# THROUGH the lake, never bypass it through a footprint cell) and the CRR pond-capture receiver.
# The UZF/DRN-seep EXCLUSION set stays EMPTY — the 45-yr-converged configuration: footprint cells
# keep UZF + the DRN-seep pin (the stabiliser at in-contact ponds); the small rain double-count
# (LAK rainfall on sarea + UZF finf on the same cells) is accepted.
node_lake = {nd: L for L, nds in enumerate(pond_nodes) for nd in nds}
lake_node_set = set()

# -----------------------------------------------------------------------------
# 6. STRESS-PERIOD TIME DISCRETISATION  (monthly steps + 1-yr spin-up)
# -----------------------------------------------------------------------------
print(">> Building monthly stress periods (with 1-yr spin-up) …")
months = pd.date_range(SIM_START, SIM_END, freq="MS")
# Spin-up: prepend first 12 months (replays year 1; always present, not counted in RUN_YEARS)
spinup = months[:12]
# Optional short test period: keep only the first RUN_YEARS years of the main
# series (see RUN_YEARS at top).  RUN_YEARS < 0 -> full period.
if RUN_YEARS is not None and RUN_YEARS > 0:
    months_main = months[:RUN_YEARS * 12]
    print(f"   TEST RUN: {RUN_YEARS} yr after spin-up -> "
          f"{len(spinup)} + {len(months_main)} = {len(spinup) + len(months_main)} stress periods")
else:
    months_main = months
    print(f"   FULL RUN: {len(spinup)} spin-up + {len(months_main)} = "
          f"{len(spinup) + len(months_main)} stress periods")
all_months = spinup.append(months_main)
# PERLEN = number of days in each month (always positive).
# Do NOT difference consecutive all_months: the 12-month spin-up repeats the
# same calendar months as the start of the main series, so the difference goes
# negative at the join -> "PERLEN cannot be less than 0.0 for stress period 12".
perlen = [int(m.days_in_month) for m in all_months]
nper = len(perlen)
nstp = [1] * nper        # daily? -> 1 step per SP means monthly stepping
tsmult = [1.0] * nper
period_data = list(zip(perlen, nstp, tsmult))
# Adaptive time stepping: try the full month first; if a period fails to converge
# (e.g. the wet near-stream discharge cells in some transient months), MF6 shrinks
# the step (dt / dtfailadj), retries, then grows it back (x dtadj) — only where
# needed, so converging months stay 1-step.
# record: (iperiod[1-based], dt0, dtmin, dtmax, dtadj, dtfailadj)
# ATS on all transient periods; excluded from steady-state period 0 (SPINUP_MODE="ss",
# where MF6 warns ATS may misbehave).  On a failure: shrink dt (/dtfailadj), retry, grow.
_ats0 = 1 if SPINUP_MODE == "ss" else 0
ATS_DTMIN = 1.0    # d, min ATS sub-step — floored at 1 d (was 0.01): sub-day steps are pointless for a monthly
                   #    model and just grind on periods that fail anyway; if dt=1 d won't converge, fail fast.
ats_records = [(i, perlen[i], ATS_DTMIN, perlen[i], 2.0, 5.0) for i in range(_ats0, nper)]  # (iperiod0, dt0, dtmin, dtmax, dtadj, dtfailadj)
print(f"   total stress periods: {nper}")

# -----------------------------------------------------------------------------
# 7. READ FORCINGS (P and ET) and convert mm -> m/day per stress period
# -----------------------------------------------------------------------------
print(">> Reading rainfall and ET0 time series …")
p_df  = pd.read_csv(P_CSV,  sep=r"\s+|,|;", engine="python", parse_dates=["Date"] if False else None)
et0_df  = pd.read_csv(ET0_CSV,  sep=r"\s+|,|;", engine="python", parse_dates=["Date"] if False else None)
# Robust parse:
#p_df  = pd.read_csv(P_CSV,  delim_whitespace=True)
#et0_df = pd.read_csv(ET0_CSV)

def parse_series(df, value_col="MONTHLY_TOTAL"):
    # find date column
    date_col = [c for c in df.columns if "date" in c.lower()][0]
    df[date_col] = pd.to_datetime(df[date_col], format="%m/%d/%Y", errors="coerce")
    df = df[[date_col, value_col]].rename(columns={date_col: "date", value_col: "val"})
    df = df.dropna().set_index("date").sort_index()
    return df

p_df  = parse_series(p_df)
et0_df = parse_series(et0_df)

# Build per-stress-period rate (m/day): monthly_total_mm / 1000 / days_in_month
def monthly_rate(series, month_start, days):
    if month_start in series.index:
        return (series.loc[month_start, "val"] / 1000.0) / days
    # fallback: use month-of-year average
    moy = series[series.index.month == month_start.month]["val"].mean()
    return (moy / 1000.0) / days

# Spin-up uses first year duplicated
rain_rate = []
et0_rate   = []
for i, m in enumerate(all_months[:nper]):
    days = perlen[i]
    if i < 12:                                            # spin-up year
        ref = SIM_START + pd.DateOffset(months=i)
    else:
        ref = m
    rain_rate.append(monthly_rate(p_df,  pd.Timestamp(ref.year, ref.month, 1), days))
    et0_rate.append( monthly_rate(et0_df, pd.Timestamp(ref.year, ref.month, 1), days))

# WT-CALIB: the SS period (period 0 = SIM_START = Jan 1981) is a FREAK-DRY month -> ~0 recharge -> an
# artificially deep steady-state water table. Drive the SS spin-up instead with a representative
# recharge = 20% of the mean annual precipitation (user-chosen recharge coefficient).  HARDCODED from
# the data: mean annual P = 508.8 mm/yr (p_month CSV, 45 complete years 1981-2025) x 0.20 = 101.8 mm/yr.
SS_RECHARGE_MMYR = 101.8     # = 0.20 x mean-annual P (hardcoded; see E:\tmp_claude\recharge_20pct.py)
if SPINUP_MODE == "ss":
    rain_rate[0] = SS_RECHARGE_MMYR / 1000.0 / 365.0
    et0_rate[0]  = 0.0
    print(f"   [WT-CALIB] SS recharge = 20% of mean annual P = {SS_RECHARGE_MMYR:.1f} mm/yr (hardcoded)")

# -----------------------------------------------------------------------------
# 9. IDENTIFY OUTLET CELL  ->  DRAIN
# -----------------------------------------------------------------------------

def find_cell_for_point(xy, ix_obj, xc, yc, nudge=0.01):
    """Find the Voronoi cell containing point xy. Robust to boundary cases."""
    pt = Point(xy)

    # Attempt 1: GridIntersect, silence the deprecation, allow boundary hits
    try:
        res = ix_obj.intersect(
            pt,
            return_all_intersections=True,
            geo_dataframe=False,
        )
        if len(res) > 0:
            return int(res.cellids[0])
    except TypeError:
        # older flopy without these kwargs
        res = ix_obj.intersect(pt)
        if len(res) > 0:
            return int(res.cellids[0])

    # Attempt 2: nudge the point slightly inward (toward the mesh centroid)
    cx, cy = float(np.mean(xc)), float(np.mean(yc))
    dx, dy = cx - xy[0], cy - xy[1]
    norm = (dx * dx + dy * dy) ** 0.5
    if norm > 0:
        nx = xy[0] + nudge * dx / norm
        ny = xy[1] + nudge * dy / norm
        try:
            res = ix_obj.intersect(
                Point(nx, ny),
                return_all_intersections=True,
                geo_dataframe=False,
            )
            if len(res) > 0:
                return int(res.cellids[0])
        except TypeError:
            res = ix_obj.intersect(Point(nx, ny))
            if len(res) > 0:
                return int(res.cellids[0])

    # Attempt 3 (always works): nearest cell centroid
    # On a Voronoi mesh, the nearest generator point IS the containing cell.
    tree = cKDTree(np.column_stack([xc, yc]))
    _, idx = tree.query([xy[0], xy[1]])
    print(f"   !! GridIntersect found no cell for outlet {xy}; "
          f"falling back to nearest centroid (cell {int(idx)}).")
    return int(idx)

outlet_cell = find_cell_for_point(outlet_xy, ix, xc, yc)
print(f"   outlet cell = {outlet_cell}  "
      f"(centroid: {xc[outlet_cell]:.2f}, {yc[outlet_cell]:.2f}, "
      f"top: {top[outlet_cell]:.2f} m)")

drn_cond = 100.0
# DRN at EVERY boundary outlet (not just the one GEOS intersection returns), in
# ALL layers — a full-column drain at each outlet.  Each layer's elevation is the
# outlet datum (top-0.5) clamped >= that layer's cell bottom (MF6 requires drain
# elev >= cell bottom; the 0.1 m top layer sits above top-0.5).
outlet_cells = sorted({find_cell_for_point(xy, ix, xc, yc) for xy in outlet_xys})
# A boundary outlet can land on a deactivated sliver cell -> the DRN-outlet record would be dropped (idomain check)
# and that outlet would not drain.  Snap any sliver outlet to its nearest ACTIVE cell so every outlet keeps a drain.
if sliver_set:
    _act = np.array([i for i in range(ncpl) if i not in sliver_set])
    _acttree = cKDTree(np.c_[xc[_act], yc[_act]])
    outlet_cells = sorted({(int(_act[_acttree.query([xc[c], yc[c]])[1]]) if c in sliver_set else c)
                           for c in outlet_cells})
drn_records = []
for c in outlet_cells:
    d_elev = top[c] - 0.5
    for lay in range(nlay):
        if idomain_lak[lay, c] != 1:          # skip -1 pass-through / 0 lake cells (no BC on inactive cells)
            continue
        cb = botm[lay, c]
        drn_records.append([(lay, c), d_elev if d_elev > cb else cb + 0.01, drn_cond])
drn_spd = {0: drn_records}
print(f"   DRN at {len(outlet_cells)} outlet cell(s) × {nlay} layers = {len(drn_records)} records")

# -----------------------------------------------------------------------------
# 10. BUILD MF6 SIMULATION
# -----------------------------------------------------------------------------
print(">> Assembling MODFLOW 6 simulation …")
sim = flopy.mf6.MFSimulation(
    sim_name=MODEL_NAME, exe_name=MF6_EXE, sim_ws=str(WORKSPACE), version="mf6"
)

tdis = flopy.mf6.ModflowTdis(
    sim, time_units="DAYS", nper=nper, perioddata=period_data,
    ats_perioddata={"maxats": len(ats_records), "perioddata": ats_records},
)
ims = flopy.mf6.ModflowIms(
    sim, complexity="COMPLEX", linear_acceleration="BICGSTAB",
    print_option="SUMMARY",        # SUMMARY: one line/time-step (ALL bloated the listing to 50 MB + slowed I/O)
    outer_maximum=500, inner_maximum=100,    # transient + ATS: let ATS subdivide hard periods rather than grind
    outer_dvclose=3e-2, inner_dvclose=1e-4,   # 3 cm head tolerance: the lakes lift the WT to the
    #   surface over an area near the ponds -> a seepage limit-cycle (~5 mm at TS1, ~2 cm once the
    #   lake exchange builds) that won't cross 1e-3 (blow-up is fixed; 3 cm << accuracy for a ~50 m
    #   aquifer w/ metre-scale seasonal swing). Tighten via gentler seepage smoothing in a follow-up.
)

gwf = flopy.mf6.ModflowGwf(
    sim, modelname=MODEL_NAME, save_flows=True, newtonoptions="NEWTON UNDER_RELAXATION"
)

disv = flopy.mf6.ModflowGwfdisv(
    gwf, length_units="METERS", nlay=nlay, ncpl=ncpl,
    nvert=gridprops_disv["nvert"],
    vertices=gridprops_disv["vertices"],
    cell2d=gridprops_disv["cell2d"],
    top=top, botm=botm,
    idomain=idomain_lak,   # 1 active / -1 vertical pass-through (absent unit); nothing is excavated
)

# NPF
# k_array / k33 were built per-cell (zoned by hydrostrat unit) in the layers block above
npf = flopy.mf6.ModflowGwfnpf(
    gwf, icelltype=icelltype, k=k_array, k33=k33, save_specific_discharge=True,
    alternative_cell_averaging="amt-hmk",   # UNTRIED LEVER (after MF6 ex-gwf-sfr-p01b): arithmetic-mean saturated THICKNESS /
    #   harmonic-mean K for interface conductance.  MF6-recommended for water-table flow with a strong K contrast (here K=40
    #   terraces abutting K=1.5 formation) — keeps inter-cell conductance well-behaved as the perched-pond connection cells
    #   dewater, vs the default harmonic mean (collapses conductance to ~0 there, a conditioning trap).
)

# Storage
sto = flopy.mf6.ModflowGwfsto(
    gwf, iconvert=icelltype, ss=1e-5, sy=0.15,
    # SPINUP_MODE: "ss" -> steady-state period 0 (equilibrates; no-lake only);
    # "transient" -> all transient (lakes converge this way).  See hybrid notes (§0).
    steady_state=({0: True} if SPINUP_MODE == "ss" else None),
    transient=({i: True for i in range(1, nper)} if SPINUP_MODE == "ss"
               else {i: True for i in range(nper)}),
)

# Initial heads: start the water table ~2 m BELOW land surface, not at it.
# strt = top everywhere starts the aquifer fully saturated, so step 1 sheds huge
# storage out through UZF seepage / SFR (the SP1 non-convergence).  The
# steady-state first period (STO) then relaxes this to equilibrium.
INIT_WT_DEPTH = 1.0   # m below land surface (used only when INIT_HEADS_FILE is None)
_ic_arr = None
if INIT_HEADS_FILE is not None and Path(INIT_HEADS_FILE).exists():
    _ic_arr = np.load(INIT_HEADS_FILE)
    if _ic_arr.size != nlay * ncpl:
        # Stale IC left over from a DIFFERENT grid (e.g. after refining the Voronoi mesh):
        # never reshape-crash — fall back to a cold start and tell the user to regenerate it.
        print(f"   !! {Path(INIT_HEADS_FILE).name} has {_ic_arr.size} values but this grid needs "
              f"{nlay*ncpl} (nlay={nlay} x ncpl={ncpl}). STALE IC from a different grid — cold-starting "
              f"instead. Run once with SPINUP_MODE='ss' to regenerate it for this grid.")
        _ic_arr = None
if _ic_arr is not None:
    strt = _ic_arr.reshape(nlay, ncpl)                    # hybrid: equilibrium heads from the no-lake SS run
    print(f"   initial heads loaded from {Path(INIT_HEADS_FILE).name} (hybrid spin-up)")
else:
    if SPINUP_MODE != "ss":
        raise RuntimeError(
            f"CONFIG ERROR: transient run but no usable IC ({SPINUP_HEADS_OUT.name} missing or stale). "
            "A transient run must start from the SS equilibrium heads — run once with SPINUP_MODE='ss' "
            "to regenerate them. (A cold start is only legitimate for the 'ss' run itself.)")
    strt = np.tile((top - INIT_WT_DEPTH).reshape(1, ncpl), (nlay, 1))
    print(f"   initial heads = cold start (top - {INIT_WT_DEPTH:.1f} m)")
ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

# DRN at outlet
drn = flopy.mf6.ModflowGwfdrn(
    gwf, maxbound=len(drn_records), stress_period_data=drn_spd, pname="DRN-OUT",
    mover=True,
)

# Surface-seepage drain — replaces UZF simulate_gwseep (deprecated + non-smooth,
# which caused single-cell limit cycles).  A DRN at land surface in the top layer
# of every cell, with cubic smoothing over DDRN so groundwater discharge ramps in
# gradually instead of switching on/off.  mover=True -> routed to SFR (baseflow).
DRN_SEEP_COND = 10000.0    # m2/d per cell (seepage-face conductance). REVERTED to 10000 (2026-06-25): dropping it to
                           # 1000 BACKFIRED — the firm high-COND pin holds the WT at the surface (stable); the softer
                           # 1000 let it float -> SP2 ground for 1645 s CPU (the COND=10000 run reached SP24). Keep 10000.
DRN_SEEP_DDRN = 3.5       # m, cubic smoothing depth (widened: damps the 2-pond seepage toggle at 2422/2480)
# land-surface seepage drain on every NON-lake cell (lakes handle their own surface), placed in the
# TOP ACTIVE layer of the column (top_act; layer 0 is pass-through wherever U1 is absent)
seep_cells = [icell for icell in range(ncpl)
              if icell not in lake_node_set and icell not in sliver_set and top_act[icell] >= 0]
seep_records = [[(int(top_act[icell]), icell), float(top[icell]), DRN_SEEP_COND, DRN_SEEP_DDRN]
                for icell in seep_cells]
drnseep_id = {icell: i for i, icell in enumerate(seep_cells)}   # cell -> DRN-SEEP boundary index (for MVR)
drn_seep = flopy.mf6.ModflowGwfdrn(
    gwf, maxbound=len(seep_records), stress_period_data={0: seep_records},
    auxiliary=["ddrn"], auxdepthname="ddrn", mover=True, pname="DRN-SEEP",
)

# -----------------------------------------------------------------------------
# 10c. HEAD OBSERVATIONS  (time series at obs_points from the geopackage)
# -----------------------------------------------------------------------------
# Continuous OBS output -> WORKSPACE/<MODEL_NAME>.obs.head.csv: one column per
# (point, layer), one row per time step.  obs_points has no screen depths, so we
# observe head in EVERY layer at each point's cell (post-processing picks/plots).
print(">> Building head observations at obs_points …")
obs_gdf = gpd.read_file(WATERSHED_GPKG, layer=OBS_PTS_LAYER)
obs_gdf = ensure_projected(obs_gdf, OBS_PTS_LAYER)
_namecol = "Name" if "Name" in obs_gdf.columns else obs_gdf.columns[0]
obs_cells = {}                                    # point name -> node (0-based)
for _i, _row in obs_gdf.iterrows():
    _g = _row.geometry
    if _g is None or _g.is_empty:
        continue
    _nm = str(_row[_namecol]).strip().replace(" ", "_")
    obs_cells[_nm] = find_cell_for_point((_g.x, _g.y), ix, xc, yc)
OBS_HEAD_CSV = f"{MODEL_NAME}.obs.head.csv"
obs_records = [(f"{nm}_L{lay + 1}", "head", (lay, node))
               for nm, node in obs_cells.items() for lay in range(nlay)
               if idomain_lak[lay, node] == 1]                # skip -1 pass-through / 0 lake cells (not observable)
head_obs = flopy.mf6.ModflowUtlobs(
    gwf, pname="head_obs", print_input=False,
    continuous={OBS_HEAD_CSV: obs_records},
)
print(f"   {len(obs_cells)} obs point(s) -> {len(obs_records)} head observation(s) "
      f"(active layers only) -> {OBS_HEAD_CSV}")

# (The old §10c "recharge-fallback" pond runoff-boost is RETIRED: CRR (§13) routes each pond's
#  upslope runoff into its LAK receiver physically, replacing the C_RUNOFF×catchment finf boost.)

# -----------------------------------------------------------------------------
# 11. UZF PACKAGE  (sandy-clay defaults)
# -----------------------------------------------------------------------------
print(">> Building UZF package …")
# thtr/thts/thti/vks are PER-CELL from the soil rasters (section 4b); eps
# (Brooks-Corey exponent) has no raster -> kept uniform.
eps = 4.0
# Two passes: first assign iuzno to every active UZF cell, then build packagedata
# so the vertical connection (ivertcon) points to the correct cell below.
# (lake_node_set is EMPTY in the embedded design — footprint cells keep their UZF.)
uzf_map = {}        # (lay, cell) -> iuzno   (idomain<=0 pass-through/inactive cells excluded)
iuzno = 0
for lay in range(nlay):
    for icell in range(ncpl):
        if icell in lake_node_set:
            continue
        if idomain_lak[lay, icell] <= 0:                 # skip -1 pass-through (pinched unit-layer) cells
            continue
        uzf_map[(lay, icell)] = iuzno
        iuzno += 1
nuzfcells = iuzno
uzf_pkdata = []
for lay in range(nlay):
    for icell in range(ncpl):
        if (lay, icell) not in uzf_map:
            continue
        landflag = 1 if lay == top_act[icell] else 0     # land surface = the TOP ACTIVE layer of the column (not layer 0)
        ivertcon = -1                                    # connect to the next ACTIVE UZF cell below (hop over pass-through layers)
        for _lb in range(lay + 1, nlay):
            if (_lb, icell) in uzf_map:
                ivertcon = uzf_map[(_lb, icell)]; break
        # surfdep smooths the UZF groundwater-seepage discharge as the water table
        # crosses land surface (must stay < cell thickness).
        surfdep = 0.5 if landflag else 0.05
        uzf_pkdata.append([
            uzf_map[(lay, icell)], (lay, icell), landflag, ivertcon, surfdep,
            vks_cell[icell], thtr_cell[icell], thts_cell[icell], thti_cell[icell], eps,
        ])

# Period data + UZF package — the big nper x ncpl block.  Build/attach it only
# for a full write; in fast-rerun mode reuse the cdl_gwf.uzf already on disk.
# (uzf_map above is always built so MVR can still reference the UZF cells.)
if WRITE_ALL:
    # Period data: apply P (finf) and ET on top layer only
    uzf_perioddata = {}
    for kper in range(nper):
        rows = []
        finf = max(rain_rate[kper], 0.0)
        pet  = max(et0_rate[kper], 0.0)
        for icell in range(ncpl):
            _ks = int(top_act[icell])            # forcing goes on the LAND-SURFACE UZF cell = top active layer
            if _ks < 0 or (_ks, icell) not in uzf_map:   # fully-inactive column / lake column -> no UZF forcing
                continue
            iuz = uzf_map[(_ks, icell)]
            rows.append([iuz, finf, pet, 2.0, float(thtr_cell[icell]), 0.5, 2.0, 1.0])
            # iuzno finf pet extdp extwc(=thtr) ha hroot rootact
        uzf_perioddata[kper] = rows

    uzf = flopy.mf6.ModflowGwfuzf(
        gwf, simulate_et=True, linear_gwet=True,   # simulate_gwseep removed -> smoothed DRN-SEEP instead
        unsat_etwc=True, ntrailwaves=7, nwavesets=40,
        mover=True, save_flows=True,
        budget_filerecord=f"{MODEL_NAME}.uzf.bud",   # unsaturated-zone compartment budget
        nuzfcells=nuzfcells, packagedata=uzf_pkdata, perioddata=uzf_perioddata,
        pname="UZF",
    )
else:
    print("   (fast rerun: skipping UZF build/write; reusing cdl_gwf.uzf on disk)")

# (The LAK package is built in section 12b, AFTER SFR, so its outlets can be wired
#  to the stream reaches the on-channel ponds hand off to — see MVR coupling.)

# -----------------------------------------------------------------------------
# 8b. DEM-BASED STREAM ROUTING  (build SFR ID / toID from the MDT)
# -----------------------------------------------------------------------------
print(">> Building DEM-based stream routing (ID/toID) …")


class _DSU:
    """Union-find for node clustering and connected components."""
    def __init__(self, n):
        self.p = list(range(n)); self.r = [0] * n
    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]; x = self.p[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.r[ra] < self.r[rb]:
            ra, rb = rb, ra
        self.p[rb] = ra
        if self.r[ra] == self.r[rb]:
            self.r[ra] += 1


def _line_parts(geom):
    if geom is None or geom.is_empty:
        return []
    if geom.geom_type == "LineString":
        return [geom]
    if geom.geom_type in ("MultiLineString", "GeometryCollection"):
        out = []
        for p in geom.geoms:
            out += _line_parts(p)
        return out
    return []


def _clean_lines(lines):
    out = []
    for ln in lines:
        c = list(dict.fromkeys(ln.coords))
        if len(c) >= 2 and LineString(c).length > 0:
            out.append(LineString(c))
    return out


def _bridge_gaps(segs, tol):
    """Snap each dangling endpoint onto the nearest OTHER line within `tol`
    (a connector), then re-node. Fixes T-junction gaps so tributaries connect."""
    tree = STRtree(segs)
    conns = []
    for k, s in enumerate(segs):
        for p in (Point(s.coords[0]), Point(s.coords[-1])):
            best_d, best_q = tol + 1.0, None
            for j in tree.query(p.buffer(tol)):
                j = int(j)
                if j == k:
                    continue
                d = p.distance(segs[j])
                if 0.0 < d < best_d:
                    best_d, best_q = d, segs[j].interpolate(segs[j].project(p))
            if best_q is not None and best_d <= tol:
                conns.append(LineString([(p.x, p.y), (best_q.x, best_q.y)]))
    if not conns:
        return segs, 0
    return _clean_lines(_line_parts(linemerge(unary_union(segs + conns)))), len(conns)


def _build_routing_gdf(in_lines, ws_poly, dem_tif, model_crs,
                       bridge_tol, snap_tol, outlet_xy=None):
    """Clip to watershed, node + bridge gaps, then derive ID/toID routing with
    flow direction from the DEM.  Returns a GeoDataFrame for sfrmaker."""
    # clip to the watershed so routing stays in-basin and node elevs stay in the DEM
    clipped = []
    for ln in in_lines:
        clipped += _line_parts(ln.intersection(ws_poly))
    if not clipped:
        clipped = list(in_lines)
    segs = _clean_lines(_line_parts(linemerge(unary_union(clipped))))
    nconn = 0
    if bridge_tol > 0:
        segs, nconn = _bridge_gaps(segs, bridge_tol)
    N = len(segs)
    if N == 0:
        raise RuntimeError("routing: no usable stream segments after cleaning.")

    # cluster endpoints into shared nodes
    endpts = np.array([pt for ln in segs for pt in (ln.coords[0], ln.coords[-1])])
    dsu = _DSU(len(endpts))
    for i, j in cKDTree(endpts).query_pairs(r=snap_tol):
        dsu.union(i, j)
    roots = [dsu.find(i) for i in range(len(endpts))]
    remap = {r: k for k, r in enumerate(sorted(set(roots)))}
    nid = [remap[r] for r in roots]
    nnode = len(remap)
    node_a = [nid[2 * k] for k in range(N)]
    node_b = [nid[2 * k + 1] for k in range(N)]
    node_xy = np.zeros((nnode, 2)); cntn = np.zeros(nnode)
    for e in range(len(endpts)):
        node_xy[nid[e]] += endpts[e]; cntn[nid[e]] += 1
    node_xy /= cntn[:, None]
    segs_at = defaultdict(list)
    for k in range(N):
        segs_at[node_a[k]].append(k); segs_at[node_b[k]].append(k)
    degree = {n: len(s) for n, s in segs_at.items()}

    # connected components
    dsu_c = _DSU(N)
    for n, sl in segs_at.items():
        for s in sl[1:]:
            dsu_c.union(sl[0], s)
    comp_id = {}; comp_of = [0] * N
    for k in range(N):
        comp_of[k] = comp_id.setdefault(dsu_c.find(k), len(comp_id))
    comps = defaultdict(list)
    for k in range(N):
        comps[comp_of[k]].append(k)

    # DEM elevation at each node (src.sample is fine in the activated env)
    with rasterio.open(dem_tif) as src:
        nod = src.nodata
        node_elev = np.array([v[0] for v in src.sample(
            [(float(x), float(y)) for x, y in node_xy])], dtype=float)
    if nod is not None:
        node_elev[node_elev == nod] = np.nan
    node_elev[np.abs(node_elev) > 1e30] = np.nan
    if np.isnan(node_elev).any():
        node_elev = np.where(np.isnan(node_elev), np.nanmean(node_elev), node_elev)

    forced = None
    if outlet_xy is not None:
        d2 = (node_xy[:, 0] - outlet_xy[0]) ** 2 + (node_xy[:, 1] - outlet_xy[1]) ** 2
        forced = int(np.argmin(d2))

    def pick_outlet(cn):
        if forced is not None and forced in cn:
            return forced
        leaves = [n for n in cn if degree[n] == 1]
        cand = leaves if leaves else list(cn)
        return cand[int(np.nanargmin([node_elev[n] for n in cand]))]

    # propagate routing + flow orientation UPSTREAM from each outlet (BFS)
    DOWN = [None] * N; dnode = [None] * N; vis = [False] * N
    other = lambda s, n: node_b[s] if n == node_a[s] else node_a[s]
    for sl in comps.values():
        cn = {node_a[s] for s in sl} | {node_b[s] for s in sl}
        out_node = pick_outlet(cn)
        q = deque()
        for s in segs_at[out_node]:
            if vis[s]:
                continue
            DOWN[s], dnode[s], vis[s] = -1, out_node, True
            q.append((other(s, out_node), s))
        while q:
            node, dseg = q.popleft()
            for s in segs_at[node]:
                if vis[s]:
                    continue
                DOWN[s], dnode[s], vis[s] = dseg, node, True
                q.append((other(s, node), s))
    for k in range(N):
        if not vis[k]:
            DOWN[k], dnode[k] = -1, node_b[k]

    # oriented geometry (upstream -> downstream) + attribute table
    geoms, up_e, dn_e = [], [], []
    for k in range(N):
        ln = segs[k]
        if dnode[k] == node_a[k]:
            ln = LineString(list(ln.coords)[::-1])
        geoms.append(ln)
        un = node_a[k] if dnode[k] == node_b[k] else node_b[k]
        up_e.append(round(float(node_elev[un]), 3))
        dn_e.append(round(float(node_elev[dnode[k]]), 3))
    toID = [0 if DOWN[k] == -1 else DOWN[k] + 1 for k in range(N)]
    gdf = gpd.GeoDataFrame(
        {"ID": np.arange(1, N + 1, dtype=int), "toID": np.array(toID, dtype=int),
         "up_elev": up_e, "dn_elev": dn_e, "comp": comp_of,
         "is_outlet": (np.array(toID) == 0).astype(int)},
        geometry=geoms, crs=model_crs)
    print(f"   routing: {N} segments, {len(comps)} network(s), "
          f"{int((gdf['toID'] == 0).sum())} outlet(s), {nconn} bridged gap(s)")
    return gdf


routed = _build_routing_gdf(stream_lines, ws_poly, DEM_TIF, TARGET_CRS,
                            BRIDGE_TOL, SNAP_TOL, outlet_xy)
routed.to_file(WORKSPACE / "streams_cdl_routed.shp")   # QA artifact (optional)

# -----------------------------------------------------------------------------
# 12. SFR  (built directly with FloPy on the DISV / Voronoi grid)
#     sfrmaker 0.13.2 cannot create MF6 SFR on an unstructured grid (its SFR
#     creation is hardcoded to structured row/col grids), so we intersect the
#     DEM-routed streams (section 8b) with the Voronoi grid via GridIntersect
#     and assemble ModflowGwfsfr ourselves, with cellid = (layer, node).
# -----------------------------------------------------------------------------
print(">> Building SFR directly on the Voronoi grid …")

SFR_LAYER  = 0       # top candidate layer; each reach descends to the first layer deep enough
SFR_WIDTH  = 1.0     # m,   MIN (headwater) channel width
SFR_W_MAX  = 8.0     # m,   MAX (outlet) width — reach width scales SFR_WIDTH->SFR_W_MAX by sqrt(drainage), so the
                     #      outlet trunk is not a 1 m pipe carrying the whole-catchment discharge (the reach-828 divergence)
SFR_RBTH   = 0.5     # m,   streambed thickness
SFR_RHK    = 0.1     # m/d, streambed K (was 1.0: streams leaked ~34,000 m3/d into GW, feeding the seepage-loop oscillation)
SFR_MAN    = 0.035   # Manning's n
SFR_MINSLP = 1e-3    # min reach gradient (0.1%, realistic lowland-stream slope). RAISED 1e-4 -> 1e-3 (2026-06-25):
                     # the near-flat outlet reach 828 sat exactly on the old 1e-4 floor -> ill-conditioned stage rating
                     # -> SFR-828-inflow Newton DIVERGED (residual 2.26e9) at nlay=4. 0.1% conditions the flat reaches.
SFR_MINLEN = 0.1     # m,   floor on reach length

ids   = routed["ID"].to_numpy()
toids = routed["toID"].to_numpy()
geoms = routed.geometry.to_numpy()
up_e  = routed["up_elev"].to_numpy()
dn_e  = routed["dn_elev"].to_numpy()
id2line = {int(i): k for k, i in enumerate(ids)}
incoming = {k: [] for k in range(len(routed))}          # line -> upstream line indices
for k in range(len(routed)):
    t = int(toids[k])
    if t != 0 and t in id2line:
        incoming[id2line[t]].append(k)

# process headwater lines first so global reach 0 is never a downstream target
# (a downstream reach is written as -rno; for rno 0 the sign would be lost)
order = sorted(range(len(routed)), key=lambda k: len(incoming[k]))

reach_cell, reach_len, reach_rtp, reach_grd = [], [], [], []
line_reaches = {k: [] for k in range(len(routed))}
no_hits = 0
for k in order:
    ln = geoms[k]
    L = ln.length
    try:                                       # flopy >= 3.10 warns unless geo_dataframe set
        res = ix.intersect(ln, geo_dataframe=False)
    except TypeError:
        res = ix.intersect(ln)
    cells = list(res["cellids"])
    lens = list(res["lengths"])
    if len(cells) == 0:
        no_hits += 1
        continue
    # order the per-cell pieces from upstream to downstream along the line
    try:
        shapes = list(res["ixshapes"])
        pos = [ln.project(shapes[i].centroid) for i in range(len(cells))]
    except Exception:
        pos = [ln.project(Point(xc[int(cells[i])], yc[int(cells[i])])) for i in range(len(cells))]
    seq = sorted(range(len(cells)), key=lambda i: pos[i])
    grd = max((float(up_e[k]) - float(dn_e[k])) / L, SFR_MINSLP) if L > 0 else SFR_MINSLP
    cum = 0.0
    for i in seq:
        rl = max(float(lens[i]), SFR_MINLEN)
        mid = (cum + rl / 2.0) / L if L > 0 else 0.5
        cum += rl
        rtp = float(up_e[k]) + (float(dn_e[k]) - float(up_e[k])) * min(max(mid, 0.0), 1.0)
        line_reaches[k].append(len(reach_cell))
        reach_cell.append(int(cells[i]))
        reach_len.append(rl)
        reach_rtp.append(rtp)
        reach_grd.append(grd)

nreaches = len(reach_cell)
if nreaches == 0:
    raise RuntimeError("SFR: no reaches — routed streams did not intersect the grid.")
if no_hits:
    print(f"   !! {no_hits} routed segment(s) did not intersect the grid (skipped).")

first_reach = {k: line_reaches[k][0]  for k in range(len(routed)) if line_reaches[k]}
last_reach  = {k: line_reaches[k][-1] for k in range(len(routed)) if line_reaches[k]}

# ---- per-reach channel width scaled by DRAINAGE (accumulated upstream channel length, topological accumulation) ----
#      width = SFR_WIDTH + (SFR_W_MAX-SFR_WIDTH)*(d/dmax)**2 : the **2 exponent keeps mid/headwater reaches near
#      SFR_WIDTH (the uniform 1 m that converged to SP24) and widens ONLY the outlet TRUNK — so fixing the outlet
#      (reach 828) does NOT redistribute flow into small reaches (that whack-a-moled the divergence to reach 77).
_indeg = {k: len(incoming[k]) for k in range(len(routed))}
_q = [k for k in range(len(routed)) if _indeg[k] == 0]; _topo = []
while _q:                                                  # Kahn topological order (headwaters -> outlet)
    _k = _q.pop(); _topo.append(_k)
    _t = id2line.get(int(toids[_k]), None)
    if _t is not None and _t in _indeg:
        _indeg[_t] -= 1
        if _indeg[_t] == 0: _q.append(_t)
reach_drain = [0.0] * nreaches
for _k in _topo:
    _acc = sum(reach_drain[last_reach[u]] for u in incoming[_k] if u in last_reach)   # inflow from upstream lines
    for _r in line_reaches.get(_k, []):
        _acc += reach_len[_r]; reach_drain[_r] = _acc
_dmax = max(reach_drain) if any(reach_drain) else 1.0
reach_wid = [SFR_WIDTH + (SFR_W_MAX - SFR_WIDTH) * (d / _dmax) ** 2.0 for d in reach_drain]
print(f"   SFR width scaled by drainage: {min(reach_wid):.1f}–{max(reach_wid):.1f} m "
      f"(outlet trunk drains {_dmax/1000:.1f} km of channel)")

# ---- SFRmaker-style DOWNSTREAM-MONOTONIC bed elevations (Leaf et al. 2021, fig. 2) ------------
# The per-segment endpoint interpolation can produce beds that RISE downstream (DEM noise at the
# segment endpoints; coarser cells make it worse). Enforce a running minimum along the routing
# topology: no reach bed is ever higher than the lowest bed encountered upstream of it — exactly
# SFRmaker's smoothing rule. (Uses the Kahn topological line order computed for the widths.)
_line_out_elev = {}
_nmono = 0
for _k in _topo:
    _run = min([_line_out_elev[u] for u in incoming[_k] if u in _line_out_elev], default=np.inf)
    for _r in line_reaches.get(_k, []):
        if reach_rtp[_r] > _run:
            reach_rtp[_r] = _run
            _nmono += 1
        else:
            _run = reach_rtp[_r]
    _line_out_elev[_k] = _run
if _nmono:
    print(f"   SFRmaker monotonic bed smoothing: {_nmono} reach elevation(s) clamped (were rising downstream)")

# colocated reaches (several stream branches meeting in one coarse cell) — SFRmaker consolidates
# them for performance in huge models; at CdL's scale MF6 handles them natively, so just report.
_percell = {}
for _c in reach_cell:
    _percell[_c] = _percell.get(_c, 0) + 1
_ncoloc = sum(1 for v in _percell.values() if v > 1)
print(f"   colocated reaches: {_ncoloc} cell(s) host >1 reach (confluences; kept — MF6 routes them natively)")

packagedata, connectiondata = [], []
for k in range(len(routed)):
    LR = line_reaches[k]
    tline = id2line.get(int(toids[k]), None)
    for p, r in enumerate(LR):
        ups, dns = [], []
        if p == 0:
            ups += [last_reach[u] for u in incoming[k] if u in last_reach]
        else:
            ups.append(LR[p - 1])
        if p == len(LR) - 1:
            if tline is not None and tline in first_reach:
                dns.append(first_reach[tline])
        else:
            dns.append(LR[p + 1])
        conns = [int(u) for u in ups] + [-int(d) for d in dns]   # downstream negated
        connectiondata.append([r] + conns)
        # Place the reach in the topmost layer whose cell bottom is below the
        # streambed bottom (rtp - rbth).  MF6 requires rtp-rbth > cell bottom, and
        # the 0.1 m top layer is far too thin to contain a 0.5 m streambed, so a
        # fixed layer 0 fails for nearly every reach.
        klay = SFR_LAYER
        bed_bot = reach_rtp[r] - SFR_RBTH
        # descend past layers too shallow for the streambed AND past excavated
        # (lake) cells, so no reach is placed in a deactivated pond cell.
        while klay < nlay - 1 and (botm[klay, reach_cell[r]] >= bed_bot
                                   or idomain_lak[klay, reach_cell[r]] <= 0):
            klay += 1
        packagedata.append([
            r, (klay, reach_cell[r]), reach_len[r], reach_wid[r],
            reach_grd[r], reach_rtp[r], SFR_RBTH, SFR_RHK, SFR_MAN,
            len(conns), 1.0, 0,
        ])
packagedata.sort(key=lambda x: x[0])
connectiondata.sort(key=lambda x: x[0])

# --- Excise SFR reaches that fall in lake (excavated) cells --------------------
# On-channel ponds: the stream passes through cells that are now lakes.  Such a
# reach would compete with the lake for the same GW cell (destabilising the solve)
# and sit beneath a deactivated cell.  Remove those reaches and renumber; flow is
# reconnected through the lakes via MVR (sfr_to_lake / lake_to_sfr) below.
sfr_to_lake, lake_to_sfr = [], []        # (reach_new, lake) ; (lake, reach_new)
if nlakes:
    lake_reaches = {r for r in range(nreaches) if reach_cell[r] in node_lake}
    if lake_reaches:
        conn_of = {row[0]: row[1:] for row in connectiondata}
        # detect stream<->lake hand-offs from the ORIGINAL connectivity
        for r in range(nreaches):
            if r in lake_reaches:
                continue
            for c in conn_of.get(r, []):
                if c < 0 and (-c) in lake_reaches:        # r flows into a lake reach
                    sfr_to_lake.append((r, node_lake[reach_cell[-c]]))
                elif c > 0 and c in lake_reaches:         # a lake reach flows into r
                    lake_to_sfr.append((node_lake[reach_cell[c]], r))
        downstream_targets = {(-c) for row in connectiondata for c in row[1:] if c < 0}
        keep = [r for r in range(nreaches) if r not in lake_reaches]
        if keep and keep[0] in downstream_targets:        # keep a headwater as new reach 0 (avoid -0)
            hw = [r for r in keep if r not in downstream_targets]
            if hw:
                keep.remove(hw[0]); keep.insert(0, hw[0])
        newno = {old: i for i, old in enumerate(keep)}
        pkg_of = {row[0]: row for row in packagedata}
        new_pkg, new_conn = [], []
        for old in keep:
            nn = newno[old]
            conns = [(1 if c >= 0 else -1) * newno[abs(c)]
                     for c in conn_of[old] if abs(c) not in lake_reaches]
            new_conn.append([nn] + conns)
            row = list(pkg_of[old]); row[0] = nn; row[9] = len(conns)
            new_pkg.append(row)
        # remap the MVR hand-off reach numbers to the new numbering
        sfr_to_lake = [(newno[r], L) for r, L in sfr_to_lake]
        lake_to_sfr = [(L, newno[r]) for L, r in lake_to_sfr]
        packagedata, connectiondata = new_pkg, new_conn
        reach_cell = [reach_cell[old] for old in keep]     # reindex by new reach number
        nreaches = len(keep)
        print(f"   excised {len(lake_reaches)} reach(es) in lake cells -> {nreaches} reaches; "
              f"MVR hand-offs: {len(sfr_to_lake)} in, {len(lake_to_sfr)} out")

sfr = flopy.mf6.ModflowGwfsfr(
    gwf, save_flows=True, mover=True, pname="SFR",
    unit_conversion=86400.0,             # meters + DAYS (Manning's coeff is m^1/3/s)
    budget_filerecord=f"{MODEL_NAME}.sfr.bud",   # surface-water (stream) compartment budget
    nreaches=nreaches,
    packagedata=packagedata,
    connectiondata=connectiondata,
)

# reach -> its cell centroid, for the MVR mapping below
reach_xy = np.array([[xc[c], yc[c]] for c in reach_cell])
print(f"   SFR reaches: {nreaches}  (from {len(routed)} routed segments)")

# -----------------------------------------------------------------------------
# 12b. LAK PACKAGE  (cattle ponds as lakes; built after SFR for the MVR coupling)
#   Vertical bed seepage to the thick aquifer below + P/ET on the lake surface.
#   Every lake gets a MANNING spill outlet at ~rim (bounds the stage); for the
#   on-channel ponds, MVR (section 13) routes upstream SFR -> lake and the lake
#   outlet -> downstream SFR.  Stage + budget saved per pond for post-processing.
# -----------------------------------------------------------------------------
if nlakes:
    print(">> Building LAK package …")
    lak_period = {}
    for kper in range(nper):
        recs = []
        for L in range(nlakes):
            recs.append([L, "RAINFALL", max(rain_rate[kper], 0.0)])
            recs.append([L, "EVAPORATION", max(et0_rate[kper], 0.0)])
        lak_period[kper] = recs
    # ONE MANNING spill outlet per lake (outletno = L), invert ~0.5 m below the rim, EXTERNAL.
    # §13 MVR routes it to the downstream SFR reach for on-channel ponds (lake_to_sfr);
    # off-channel ponds spill out of the model.
    #   record: outletno, lakein, lakeout, couttype, invert, width, rough, slope
    # ⚠ flopy 0-based GOTCHA (found 2026-07-04, the pond-4 TS1 blow-up): lakeout is a LAKE-NUMBER
    # field, so flopy ADDS 1 on write. "External" must be flopy -1 (-> MF6 0). lakeout=0 means
    # LAKE 1 = pond 4 -> every spill in the catchment poured into pond 4 and detonated its cell.
    LAK_EXTERNAL = -1
    lak_outlets = []
    for L in range(nlakes):
        rim = float(np.mean(top[pond_nodes[L]])) - 0.5
        lak_outlets.append([L, L, LAK_EXTERNAL, "MANNING", rim, 1.0, 0.03, 1e-3])

    # ---- STAGE-VOLUME-AREA tables (one per lake) so a lake DRIES GRACEFULLY -------
    # The excavated pond has a FLAT bottom -> the wetted surface area jumps 0->full as
    # a step at the bottom, so a bone-dry perched lake degenerates (SP24: LAK residual
    # stuck at 2.7e18 while GW is converged to 1e-13).  An assumed wedge bathymetry
    # (sarea: ~0 at the deepest point -> full footprint at the rim) replaces that step
    # with a SMOOTH sarea->0, so storage/rainfall/evaporation vanish smoothly when dry.
    def _lake_svt(lbot, gmean, area, nw=8, headroom=2.0):
        D = max(gmean - lbot, 0.1)
        smin = max(area * 0.01, 1.0)                  # tiny nonzero area at the bottom (avoids exact-zero degeneracy)
        stages = [lbot + D * i / nw for i in range(nw + 1)] + [gmean + headroom]
        rows, vol, ps, pa = [], 0.0, None, None
        for s in stages:
            a = area if s >= gmean else smin + (area - smin) * (s - lbot) / D
            if ps is not None:
                vol += 0.5 * (a + pa) * (s - ps)      # trapezoidal volume = integral of sarea(stage)
            rows.append((round(s, 4), round(vol, 6), round(a, 4)))
            ps, pa = s, a
        return rows

    lak_tables = []
    for L in range(nlakes):
        _lb, _gm, _ar = lak_geom[L]
        _tab = _lake_svt(_lb, _gm, _ar)                          # rows = (stage, volume, sarea)
        _tab = [(s, v, a, a) for (s, v, a) in _tab]              # EMBEDDEDV needs ncol=4: barea = sarea (wetted bed = exchange area)
        flopy.mf6.ModflowUtllaktab(
            gwf, nrow=len(_tab), ncol=4, table=_tab,
            filename=f"{MODEL_NAME}.lak{L + 1}.tab", pname=f"laktab_{L + 1}",
        )
        lak_tables.append([L, f"{MODEL_NAME}.lak{L + 1}.tab"])   # flopy adds the TAB6/FILEIN keywords itself
        print(f"   lake {L} (FID {pond_meta[L]['fid']}) stage-vol-sarea-barea table: "
              f"{len(_tab)} rows, area {_ar:.0f} m², stage {_tab[0][0]:.2f}–{_tab[-1][0]:.2f} m")

    lak = flopy.mf6.ModflowGwflak(
        gwf, pname="LAK", boundnames=True,
        print_stage=True, print_flows=True, save_flows=True,
        budget_filerecord=f"{MODEL_NAME}.lak.bud",
        stage_filerecord=f"{MODEL_NAME}.lak.stage",
        mover=True,                        # coupled to SFR via MVR
        length_conversion=1.0, time_conversion=86400.0,   # FID 12 OVERFILL FIX: the MANNING outlet eqn
        #   needs the m-days unit conversion (like SFR's unit_conversion=86400). Without it the outlet
        #   discharged ~86400x too little per day, so the spilling on-channel pond (FID 12) backed up
        #   metres above its rim. With it, the pond spills freely at the invert (stays at rim level).
        surfdep=LAK_SURFDEP,               # PHASE-0: smooth vertical connection wetted area -> convergence (MF6 TM6-A55)
        maximum_iterations=LAK_MAXITER,    # LAK-internal Newton cap (default 100) — untried lever for the deep-U3 overflow
        maximum_stage_change=LAK_STAGECHG, # LAK stage closure tol (default 1e-5); looser -> no overshoot into float overflow
        nlakes=nlakes, noutlets=len(lak_outlets), ntables=nlakes,
        packagedata=lak_pkg, connectiondata=lak_conn,
        outlets=lak_outlets, perioddata=lak_period,
        tables=lak_tables,                 # stage-volume-area per lake -> graceful drying (no bone-dry degeneracy)
    )
    print(f"   LAK: {nlakes} lake(s), {len(lak_conn)} connection(s), "
          f"{len(lak_outlets)} outlet(s), P/ET forcing on {nper} periods")

# -----------------------------------------------------------------------------
# 13. MVR  (UZF rejected infiltration + GW-seep -> nearest SFR reach)
# -----------------------------------------------------------------------------
print(">> Building MVR package …")

mvr_pkgs = [["UZF"], ["SFR"], ["DRN-OUT"], ["DRN-SEEP"]]
# LAK joins MVR only when it actually provides/receives movers: the SFR<->LAK hand-offs of the
# on-channel ponds, and/or CRR runoff captured by a pond (checked in the CRR block below).
_lak_in_mvr = bool(nlakes) and (bool(sfr_to_lake) or bool(lake_to_sfr))
if _lak_in_mvr:
    mvr_pkgs.append(["LAK"])

# nearest SFR reach to each cell — vectorized (used by the legacy blanket routing)
_, nearest_reach = cKDTree(reach_xy).query(np.column_stack([xc, yc]))

if CRR_ENABLE:
    # ---- CRR (Daoud et al. 2022, Eq. 23): cascade-routing & reinfiltration --------------------
    # Each cell's rejected infiltration (UZF provider) and groundwater exfiltration (DRN-SEEP
    # provider) are split among its DOWNSLOPE shared-face neighbours with MFD slope weights
    # α_ij = β·S_ij/ΣS_ij. Receiver at neighbour j: its LAK pond (runoff capture — replaces the
    # old §10c C_RUNOFF finf-boost), else its SFR reach (direct runoff), else its land-surface
    # UZF cell (REINFILTRATION). Cells with no downslope receiver are sinks (runoff evaporates).
    from collections import defaultdict as _dd
    _e2c = _dd(list)                                   # shared Voronoi face -> the (max 2) cells on it
    for _c2 in gridprops_vg["cell2d"]:
        _icc = int(_c2[0]); _nv = int(_c2[3]); _ivl = [int(v) for v in _c2[4:4 + _nv]]
        for _a, _b in zip(_ivl, _ivl[1:] + _ivl[:1]):
            _e2c[(_a, _b) if _a < _b else (_b, _a)].append(_icc)
    _nbrs = _dd(set)
    for _cl in _e2c.values():
        if len(_cl) == 2:
            _nbrs[_cl[0]].add(_cl[1]); _nbrs[_cl[1]].add(_cl[0])
    _cell_reach = {}
    for _rr, _cc in enumerate(reach_cell):
        _cell_reach.setdefault(_cc, _rr)               # dominant (first) reach in a cell
    period_records = []
    _nsink = 0; _crr_to_lak = False
    for _i in range(ncpl):
        _prov_uzf = uzf_map.get((int(top_act[_i]), _i)) if top_act[_i] >= 0 else None
        _prov_drn = drnseep_id.get(_i)
        if _prov_uzf is None and _prov_drn is None:
            continue
        _down = []
        for _j in _nbrs[_i]:
            _dist = float(np.hypot(xc[_i] - xc[_j], yc[_i] - yc[_j])) or 1.0
            _S = float(top[_i] - top[_j]) / _dist                       # land-surface slope i -> j
            if _S <= 0.0:
                continue                                                # uphill/flat -> no flow (Eq. 23)
            if _j in node_lake:
                _rcv = ("LAK", int(node_lake[_j]))                      # pond captures the runoff
            elif _j in _cell_reach:
                _rcv = ("SFR", int(_cell_reach[_j]))                    # direct runoff to the stream
            elif top_act[_j] >= 0 and (int(top_act[_j]), _j) in uzf_map:
                _rcv = ("UZF", int(uzf_map[(int(top_act[_j]), _j)]))    # downslope REINFILTRATION
            else:
                continue
            _down.append((_S, _rcv))
        if not _down:
            _nsink += 1                                                 # topographic sink -> evaporates
            continue
        _Ssum = sum(s for s, _ in _down)
        for _S, _rcv in _down:
            _f = CRR_BETA * _S / _Ssum
            if _rcv[0] == "LAK":
                _crr_to_lak = True
            if _prov_uzf is not None:
                period_records.append(["UZF", _prov_uzf, _rcv[0], _rcv[1], "FACTOR", _f])
            if _prov_drn is not None:
                period_records.append(["DRN-SEEP", _prov_drn, _rcv[0], _rcv[1], "FACTOR", _f])
    if _crr_to_lak and not _lak_in_mvr:
        mvr_pkgs.append(["LAK"]); _lak_in_mvr = True
    print(f"   CRR: {len(period_records)} downslope mover record(s) (β={CRR_BETA}); "
          f"{_nsink} sink cell(s) evaporate their runoff")
else:
    # LEGACY blanket routing: UZF rejected infiltration AND the surface-seepage drain -> the
    # NEAREST reach (pre-CRR behaviour, kept for comparison runs).
    period_records = [
        ["UZF", uzf_map[(int(top_act[icell]), icell)], "SFR", int(nearest_reach[icell]), "FACTOR", 1.0]
        for icell in range(ncpl) if top_act[icell] >= 0 and (int(top_act[icell]), icell) in uzf_map
    ] + [
        ["DRN-SEEP", drnseep_id[icell], "SFR", int(nearest_reach[icell]), "FACTOR", 1.0]
        for icell in seep_cells
    ]
# stream <-> pond coupling (on-channel ponds): upstream reach -> lake (full outflow),
# and lake outlet -> downstream reach (split if a lake feeds several reaches).
_recv = {}
for _L, _r in lake_to_sfr:
    _recv[_L] = _recv.get(_L, 0) + 1
for _r, _L in sfr_to_lake:                                   # (1) upstream SFR reach -> LAK (full inflow)
    period_records.append(["SFR", int(_r), "LAK", int(_L), "FACTOR", 1.0])
for _L, _r in lake_to_sfr:                                   # (2) lake SPILL outlet -> downstream SFR.
    # A LAK provider's MVR id is its OUTLET number = L (one spill outlet per lake, §12b).
    period_records.append(["LAK", int(_L), "SFR", int(_r), "FACTOR", 1.0 / _recv[_L]])

# Movers are time-invariant, so specify them ONCE: MF6 reuses the previous
# period's MVR list for any period without a PERIOD block (verified with a
# 2-period test). This makes the .mvr file ncpl records instead of nper*ncpl.
mvr_perioddata = {0: period_records}

mvr = flopy.mf6.ModflowGwfmvr(
    gwf,
    maxmvr=len(period_records),
    maxpackages=len(mvr_pkgs),
    packages=mvr_pkgs,
    perioddata=mvr_perioddata,
    pname="MVR",
)
# -----------------------------------------------------------------------------
# 14. OUTPUT CONTROL
# -----------------------------------------------------------------------------
oc = flopy.mf6.ModflowGwfoc(
    gwf,
    head_filerecord=f"{MODEL_NAME}.hds",
    budget_filerecord=f"{MODEL_NAME}.cbc",
    saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    printrecord=[("BUDGET", "LAST")],
)

# -----------------------------------------------------------------------------
# 14b. MAP: Voronoi grid, streams, and boundary conditions (DRN + SFR)
# -----------------------------------------------------------------------------
print(">> Plotting model map (grid, streams, boundary conditions) …")
try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as _MplPoly, Patch
    from matplotlib.collections import PatchCollection
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(10, 11))
    pmv = flopy.plot.PlotMapView(modelgrid=vgrid, ax=ax)
    pmv.plot_grid(lw=0.2, color="0.8")

    # SFR reach cells (shaded)
    _patches = [_MplPoly(vgrid.get_cell_vertices(nd)) for nd in sorted(set(reach_cell))]
    ax.add_collection(PatchCollection(_patches, facecolor="deepskyblue",
                                      alpha=0.45, edgecolor="none", zorder=2))

    # cattle ponds (LAK) — the lake-owned cells (node_lake; lake_node_set is empty by design)
    if nlakes:
        _lk = [_MplPoly(vgrid.get_cell_vertices(nd)) for nd in sorted(node_lake)]
        ax.add_collection(PatchCollection(_lk, facecolor="navy", alpha=0.65,
                                          edgecolor="none", zorder=3))
        gpd.GeoSeries(pond_polys, crs=TARGET_CRS).boundary.plot(
            ax=ax, color="navy", lw=1.0, zorder=4)

    # ALL pond footprints (outline only), drawn regardless of LAK so the map shows
    # where ponds sit vs SFR — an on-channel pond's stream-crossed cells ARE SFR
    # (deepskyblue) even when it is not activated as a lake.
    if pond_polys_all:
        gpd.GeoSeries(pond_polys_all, crs=TARGET_CRS).boundary.plot(
            ax=ax, color="darkorange", lw=0.8, zorder=4)
        for _pp, _pf in zip(pond_polys_all, pond_fids_all):     # pond id label (number only, user 2026-07-04)
            ax.annotate(f"{_pf}", (_pp.centroid.x, _pp.centroid.y), color="darkorange", fontsize=7.5,
                        fontweight="bold", ha="center", va="bottom", xytext=(0, 3),
                        textcoords="offset points", zorder=7)

    # routed stream network, colored by independent network (comp)
    routed.plot(ax=ax, column="comp", categorical=True, cmap="tab10", lw=1.3, zorder=3)

    # watershed boundary
    gpd.GeoSeries([ws_poly], crs=TARGET_CRS).boundary.plot(ax=ax, color="k", lw=1.2, zorder=4)

    # DRN outlet cells (section 9) — all boundary outlets
    ax.scatter([xc[c] for c in outlet_cells], [yc[c] for c in outlet_cells],
               marker="v", s=170, c="red", edgecolors="k", linewidths=0.8, zorder=6)

    # SFR network outlets (downstream end of toID == 0 segments)
    _out = routed[routed["toID"] == 0]
    if len(_out):
        ax.scatter([g.coords[-1][0] for g in _out.geometry],
                   [g.coords[-1][1] for g in _out.geometry],
                   marker="*", s=160, c="gold", edgecolors="k", linewidths=0.6, zorder=6)

    # observation points (section 10c)
    for _nm, _nd in obs_cells.items():
        ax.scatter(xc[_nd], yc[_nd], marker="s", s=55, c="magenta",
                   edgecolors="k", linewidths=0.8, zorder=7)
        ax.annotate(_nm, (xc[_nd], yc[_nd]), textcoords="offset points",
                    xytext=(5, 4), fontsize=8, fontweight="bold", zorder=7)

    _handles = [
        Patch(facecolor="deepskyblue", alpha=0.45, label="SFR reach cells"),
        Line2D([0], [0], color="tab:blue", lw=1.3, label="stream network"),
        Line2D([0], [0], color="k", lw=1.2, label="watershed boundary"),
        Line2D([0], [0], marker="v", color="w", markerfacecolor="red",
               markeredgecolor="k", markersize=12, label="DRN outlet"),
        Line2D([0], [0], marker="*", color="w", markerfacecolor="gold",
               markeredgecolor="k", markersize=14, label="SFR outlets"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="magenta",
               markeredgecolor="k", markersize=10, label="obs points"),
        Line2D([0], [0], color="darkorange", lw=0.8, label="pond footprints (all)"),
        Patch(facecolor="navy", alpha=0.65, label="ponds (LAK)"),
    ]
    ax.legend(handles=_handles, loc="upper right", fontsize=8)
    ax.set_title(f"CdL model — {ncpl} Voronoi cells, {nreaches} SFR reaches")
    ax.set_xlabel("X (m, EPSG:3763)"); ax.set_ylabel("Y (m)")
    ax.set_aspect("equal")
    _map_png = PREPROC_DIR / "model_map.png"
    fig.savefig(_map_png, dpi=150, bbox_inches="tight")
    print(f"   wrote {_map_png}")
    if plt.isinteractive():     # show inline in Spyder; skip in headless/CLI runs (plt.show() blocks there)
        plt.show()
except Exception as _e:
    print(f"   (map skipped: {_e!r})")

# -----------------------------------------------------------------------------
# 14b-bis. POND ZOOM PAGES (grid + SFR cells + LAK cells per pond) -> PREPROC_DIR
#   Plots the IN-MEMORY build: reach_cell is post-excision and the LAK cells come
#   straight from §5b -> immune to the stale-written-file artifact that hit the
#   standalone diag_sfr_ponds_zoom.py (2026-07-04). In an SS build (no lakes) the
#   PLANNED lake cells (same §5b selection rule) are shown hatched instead.
# -----------------------------------------------------------------------------
print(">> Plotting pond zoom pages …")
try:
    _sfr_cells = set(reach_cell)
    if nlakes:
        _zoom = [(pond_meta[L]["fid"], pond_polys[L], pond_nodes[L]) for L in range(nlakes)]
        _planned = False
    else:                        # SS build: show the PLANNED lake cells for every pond
        _zoom = [(_f, _p, pond_cells_of(_p)) for _f, _p in zip(pond_fids_all, pond_polys_all)]
        _planned = True
    _PPP = 6                     # ponds per page (2 x 3 panels)
    _npag = int(np.ceil(len(_zoom) / _PPP))
    for _pgno in range(_npag):
        _chunk = _zoom[_pgno * _PPP:(_pgno + 1) * _PPP]
        fig, axes = plt.subplots(2, 3, figsize=(16, 11))
        for ax in axes.ravel()[len(_chunk):]:
            ax.axis("off")
        for ax, (_fid, _g, _nds) in zip(axes.ravel(), _chunk):
            _x0, _y0, _x1, _y1 = _g.bounds
            _x0 -= 25; _y0 -= 25; _x1 += 25; _y1 += 25
            _win = [i for i in range(ncpl) if _x0 <= xc[i] <= _x1 and _y0 <= yc[i] <= _y1]
            ax.add_collection(PatchCollection([_MplPoly(vgrid.get_cell_vertices(i)) for i in _win],
                                              facecolor="none", edgecolor="0.75", lw=0.4))
            _sw = [i for i in _win if i in _sfr_cells]
            if _sw:
                ax.add_collection(PatchCollection([_MplPoly(vgrid.get_cell_vertices(i)) for i in _sw],
                                                  facecolor="deepskyblue", alpha=0.55, edgecolor="none"))
            if _nds:
                ax.add_collection(PatchCollection([_MplPoly(vgrid.get_cell_vertices(i)) for i in _nds],
                                                  facecolor="orange", alpha=0.75, edgecolor="crimson",
                                                  lw=1.2, hatch="//" if _planned else None))
            _xs, _ys = _g.exterior.xy
            ax.plot(_xs, _ys, color="navy", lw=1.8)
            _sx = stream_union.intersection(_g.buffer(25))
            for _ln in ([_sx] if _sx.geom_type == "LineString" else getattr(_sx, "geoms", [])):
                if _ln.geom_type == "LineString":
                    _lx, _ly = _ln.xy; ax.plot(_lx, _ly, color="tab:blue", lw=1.2)
            _nsl = len(set(_nds) & _sfr_cells)          # SFR reaches still inside LAK cells (transient target: 0)
            ax.set_title(f"pond {_fid}: {len(_nds)} LAK cell(s), SFR∩LAK = {_nsl}"
                         + ("" if _planned else " (expect 0)"), fontsize=10)
            ax.set_xlim(_x0, _x1); ax.set_ylim(_y0, _y1); ax.set_aspect("equal")
            ax.tick_params(labelsize=7)
        fig.legend(handles=[
            Patch(facecolor="deepskyblue", alpha=0.55, label="SFR reach cells (this build)"),
            Patch(facecolor="orange", alpha=0.75, edgecolor="crimson", hatch="//" if _planned else None,
                  label="PLANNED LAK cells (SS build)" if _planned else "LAK cells (this build)"),
            Patch(facecolor="none", edgecolor="navy", label="pond footprint"),
        ], loc="upper center", ncol=3, fontsize=10)
        fig.suptitle(f"Ponds — LAK vs SFR cells (page {_pgno + 1}/{_npag})", y=0.955)
        if _planned:
            fig.text(0.5, 0.935, "SS build (no LAK package): orange = PLANNED lake cells; "
                     "the transient build excises the SFR inside them", ha="center",
                     color="darkorange", fontsize=11, fontweight="bold")
        fig.savefig(PREPROC_DIR / f"ponds_zoom_p{_pgno + 1}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    print(f"   wrote {_npag} ponds_zoom_p*.png -> {PREPROC_DIR}")
except Exception as _e:
    print(f"   (pond zoom skipped: {_e!r})")

# -----------------------------------------------------------------------------
# 14c. PRE-PROCESSING PARAMETER MAPS  (soil + aquifer) -> PREPROC_DIR
# -----------------------------------------------------------------------------
print(">> Plotting parameter maps (soil + aquifer) …")
try:
    import matplotlib.pyplot as plt
    # soil (UZF) parameters from the rasters
    _soil = [("vks (m/d)", vks_cell, "viridis"), ("thts (sat)", thts_cell, "YlGnBu"),
             ("thti (field cap.)", thti_cell, "YlGnBu"), ("thtr (wilting pt)", thtr_cell, "YlOrBr")]
    fig, axes = plt.subplots(2, 2, figsize=(13, 12))
    for ax, (ttl, arr, cmap) in zip(axes.ravel(), _soil):
        pmv = flopy.plot.PlotMapView(modelgrid=vgrid, ax=ax)
        ca = pmv.plot_array(arr, cmap=cmap)
        gpd.GeoSeries([ws_poly], crs=TARGET_CRS).boundary.plot(ax=ax, color="k", lw=0.8)
        fig.colorbar(ca, ax=ax, shrink=0.7)
        ax.set_title(f"Soil — {ttl}"); ax.set_aspect("equal")
    fig.suptitle("CdL — soil-hydraulic (UZF) parameters sampled from rasters")
    fig.tight_layout()
    fig.savefig(PREPROC_DIR / "soil_params.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # aquifer (NPF) Kh per layer — ONE shared K ramp + a single colorbar (K33 = 0.1*Kh).  The 3 MF6
    # NUMERICAL layers are NOT the 3 geologic units; the titles spell out the relationship.
    _klab = [f"MF6 Layer 1 — U1 alluvium (K={KH_ALLU:.0f})",
             f"MF6 Layer 2 — U2 terraces (K={KH_TERR:.0f})",
             f"MF6 Layer 3 — U3 Plio-Miocene (K={KH_FORM})"]
    fig, axes = plt.subplots(1, nlay, figsize=(6 * nlay, 7), squeeze=False)
    _ca = None
    for L in range(nlay):
        ax = axes[0, L]; pmv = flopy.plot.PlotMapView(modelgrid=vgrid, ax=ax)
        _ca = pmv.plot_array(k_array[L], cmap="viridis", vmin=KH_FORM, vmax=KH_TERR)
        gpd.GeoSeries([ws_poly], crs=TARGET_CRS).boundary.plot(ax=ax, color="k", lw=0.8)
        ax.set_aspect("equal"); ax.set_title(_klab[L] if L < len(_klab) else f"MF6 Layer {L + 1}", fontsize=10)
    _cb = fig.colorbar(_ca, ax=axes.ravel().tolist(), shrink=0.85, pad=0.02)
    _cb.set_label("Kh (m/d) — shared ramp")
    _cb.set_ticks([KH_FORM, KH_ALLU, KH_TERR])
    _cb.set_ticklabels([f"{KH_FORM} (U3 formation)", f"{KH_ALLU} (U1 alluvium)", f"{KH_TERR} (U2 terraces)"])
    fig.suptitle("CdL — horizontal K per MF6 layer (K33 = 0.1*Kh). ONE numerical layer per unit: "
                 "HOST (top 4 m, K = outcropping unit) + U1 (25) + U2 (40) + U3 (1.5).",
                 fontsize=11)
    fig.savefig(PREPROC_DIR / "aquifer_params.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # hydrostratigraphic units — top / bottom / thickness per cell (from build_layers.py -> voronoi_layers.npz)
    _ut, _ub = _vl["top"], _vl["botm"]          # top + botm[3] = U1/U2/U3 boundaries
    _units = [("U1 alluvium", _ut, _ub[0]), ("U2 terraces", _ub[0], _ub[1]), ("U3 Plio-Miocene", _ub[1], _ub[2])]
    fig, axes = plt.subplots(3, 3, figsize=(16, 15), squeeze=False)
    for _r, (_lab, _t, _b) in enumerate(_units):
        for _c, (_arr, _ttl, _cm, _un) in enumerate([(_t, "TOP", "viridis", "m a.s.l."),
                                                      (_b, "BOTTOM", "viridis", "m a.s.l."),
                                                      (_t - _b, "THICKNESS", "YlGnBu", "m")]):
            ax = axes[_r, _c]
            pmv = flopy.plot.PlotMapView(modelgrid=vgrid, ax=ax)
            ca = pmv.plot_array(np.asarray(_arr, float), cmap=_cm)
            gpd.GeoSeries([ws_poly], crs=TARGET_CRS).boundary.plot(ax=ax, color="k", lw=0.8)
            fig.colorbar(ca, ax=ax, shrink=0.7); ax.set_aspect("equal")
            ax.set_title(f"{_lab} — {_ttl} ({_un})", fontsize=9)
    fig.suptitle("CdL — hydrostratigraphic units: top / bottom / thickness (sampled on the Voronoi grid)")
    fig.tight_layout()
    fig.savefig(PREPROC_DIR / "layer_surfaces.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # idomain (ibound) per MF6 layer + UZF surface boundary (iuzfbnd at the top-active layer)
    from matplotlib.colors import ListedColormap as _LC, BoundaryNorm as _BN
    _iuzfbnd = np.array([1 if top_act[_ic] >= 0 and (int(top_act[_ic]), _ic) in uzf_map else 0
                         for _ic in range(ncpl)])
    _bc = _LC(["#9ecae1", "crimson", "0.85"]); _bn = _BN([-1.5, -0.5, 0.5, 1.5], _bc.N)  # -1 pass-through / 0 lake / 1 active
    _pan = [(idomain_lak[L], f"idomain — MF6 Layer {L + 1}") for L in range(nlay)] + \
           [(_iuzfbnd, "uzf_iuzfbnd (surface)")]
    fig, axes = plt.subplots(1, len(_pan), figsize=(5.4 * len(_pan), 6.8), squeeze=False)
    _im = None
    for ax, (arr, ttl) in zip(axes[0], _pan):
        arr = np.asarray(arr); pmv = flopy.plot.PlotMapView(modelgrid=vgrid, ax=ax)
        _im = pmv.plot_array(arr, cmap=_bc, norm=_bn)
        gpd.GeoSeries([ws_poly], crs=TARGET_CRS).boundary.plot(ax=ax, color="k", lw=0.8)
        _inact = np.where(arr == 0)[0]                            # lake cells (tiny) -> mark them visibly
        ax.scatter(xc[_inact], yc[_inact], s=12, c="crimson", edgecolors="none", zorder=5)
        for _pp, _pm in zip(pond_polys, pond_meta):
            gpd.GeoSeries([_pp], crs=TARGET_CRS).boundary.plot(ax=ax, color="navy", lw=1.3, zorder=6)
            ax.annotate(f"{_pm['fid']}", (_pp.centroid.x, _pp.centroid.y), color="navy", fontsize=7.5,
                        fontweight="bold", ha="center", va="bottom", xytext=(0, 4),
                        textcoords="offset points", zorder=7)
        ax.set_aspect("equal"); ax.set_title(ttl, fontsize=10)
    _cb = fig.colorbar(_im, ax=axes.ravel().tolist(), shrink=0.8, pad=0.02, ticks=[0, 1])
    _cb.set_ticklabels(["0 = inactive / no UZF (lake)", "1 = active / UZF"])
    fig.suptitle("CdL — idomain (ibound) per MF6 layer + UZF surface boundary "
                 "(-1 pass-through where a unit is absent; no cells are excavated)", fontsize=11)
    fig.savefig(PREPROC_DIR / "ibound_uzf_map.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"   wrote soil_params.png, aquifer_params.png, layer_surfaces.png, ibound_uzf_map.png -> {PREPROC_DIR}")
except Exception as _e:
    print(f"   (parameter maps skipped: {_e!r})")

# -----------------------------------------------------------------------------
# 15. WRITE & RUN
# -----------------------------------------------------------------------------
if WRITE_ALL:
    print(">> Writing all MODFLOW 6 input files …")
    sim.write_simulation()
else:
    print(f">> Fast write: refreshing only {WRITE_THESE} (reusing the rest on disk) …")
    for _name in WRITE_THESE:
        _pkg = sim.get_package(_name)
        if _pkg is None:
            try:
                _pkg = gwf.get_package(_name)
            except Exception:
                _pkg = None
        if _pkg is None:
            raise RuntimeError(
                f"WRITE_THESE: package '{_name}' not found — use WRITE_ALL=True for a full "
                f"write, or fix the name (e.g. 'ims','npf','sto','ic','DRN-OUT','SFR','MVR','oc').")
        _pkg.write()
        print(f"   refreshed {_name}")

print(">> Running MODFLOW 6 …")
success, buff = sim.run_simulation(silent=False)
if not success:
    # AUTO-FORENSICS on divergence: parse the listing + map WHERE it failed
    # -> diag\divergence_map.png + divergence_report.csv (EPSG:3763 coords for GIS).
    try:
        import subprocess, sys
        _diag = Path(__file__).resolve().parent / "diag_divergence.py"
        print(f">> Divergence — running forensics ({_diag.name}) …")
        subprocess.run([sys.executable, str(_diag), str(WORKSPACE / "mfsim.lst")], check=False)
        print(f"   -> {WORKSPACE / 'diag' / 'divergence_map.png'}  (+ divergence_report.csv)")
    except Exception as _fe:
        print(f"   (auto-forensics step skipped: {_fe})")
    raise RuntimeError("MODFLOW 6 did not converge — see the listing file + diag\\divergence_map.png")
print(">> Simulation finished successfully.")

# -----------------------------------------------------------------------------
# 16. QUICK POST-PROCESS
# -----------------------------------------------------------------------------
print(">> Post-processing: heads & budget summary …")
head_obj = gwf.output.head()
times = head_obj.get_times()
head_last = head_obj.get_data(totim=times[-1])

# Hybrid spin-up: save the steady-state (period 0) equilibrium heads so a following
# transient lake run can load them as initial heads (INIT_HEADS_FILE=SPINUP_HEADS_OUT).
if SPINUP_MODE == "ss":
    np.save(SPINUP_HEADS_OUT, np.asarray(head_obj.get_data(kstpkper=(0, 0))).reshape(nlay, ncpl))
    print(f"   saved SS equilibrium heads -> {SPINUP_HEADS_OUT}")
    # Auto-run the SS observed-vs-computed piezometer scatter (user request, 2026-06-27) — it reads the
    # spinup_heads.npy just written + the qualitative piezo targets -> preproc/<stamp>/ss_obs_vs_computed_heads.png.
    try:
        import subprocess, sys
        _ssoc = Path(__file__).resolve().parent / "ss_obs_vs_computed.py"
        print(f">> SS calibration check: running {_ssoc.name} …")
        subprocess.run([sys.executable, str(_ssoc)], check=False, env=dict(os.environ, MPLBACKEND="Agg"))
    except Exception as _se:
        print(f"   (ss_obs_vs_computed skipped: {_se})")

bud_obj = gwf.output.budget()
records = bud_obj.get_unique_record_names()
print("   Budget record names:", records)

# Save a CSV of mean head per layer for QA
mean_h = [np.nanmean(np.where(np.abs(head_last[lay]) < 1e10, head_last[lay], np.nan))
          for lay in range(nlay)]
pd.DataFrame({"layer": np.arange(1, nlay + 1), "mean_head_m": mean_h})\
  .to_csv(WORKSPACE / "mean_head_per_layer.csv", index=False)

print(">> Done.  Model workspace:", WORKSPACE.resolve())

# -----------------------------------------------------------------------------
# 17. AUTO POST-PROCESSING  (only on a SUCCESSFUL run — a non-convergence raises above)
# -----------------------------------------------------------------------------
# Runs postprocess_cdl.py as a separate process (same interpreter), reading
# last_run_stamp.txt -> postproc\<RUN_STAMP>\.  Mirrors the diag_divergence hook;
# MPLBACKEND=Agg forces headless figure generation so it never blocks on a GUI.
if RUN_POSTPROC:
    try:
        import subprocess, sys
        _pp = Path(__file__).resolve().parent / "postprocess_cdl.py"
        print(f">> Auto post-processing: running {_pp.name}  (-> postproc\\{RUN_STAMP}) …")
        _rc = subprocess.run([sys.executable, str(_pp)], check=False,
                             env=dict(os.environ, MPLBACKEND="Agg")).returncode
        if _rc == 0:
            print(f">> Post-processing complete -> {WORKSPACE / 'postproc' / RUN_STAMP}")
        else:
            print(f"   !! post-processing exited with code {_rc} — see its output above.")
    except Exception as _pe:
        print(f"   (auto post-processing skipped: {_pe})")