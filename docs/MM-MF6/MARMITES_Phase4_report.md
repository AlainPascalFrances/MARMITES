# MARMITES — Phase 4 completion report

**Date:** 2026-07-21
**Model:** Claude Opus 4.8
**Scope:** Phase 4 of `MARMITES_code_review.md` §5/§6 — unstructured-grid
support via **DISV**: cell-geometry abstraction, DISV model builder, coupler
node mapping, grid-agnostic raster sampling, and grid-agnostic plotting.
DISU remains out of scope (locked decision).

## Result

All five work items of `MARMITES_Phase4_kickoff.md` are complete and the four
"definition of done" criteria are met. The test suite grew from **30 to 65
tests, all passing**; ruff error classes clean on every Phase-4 module.

MARMITES no longer contains a single `delr`/`delc` index in its physics or
coupling paths — the spatial dimension is fully abstracted.

---

## 1. `CellGeometry` (`trunk/marmites_grid.py`)

One provider supplies the only two spatial quantities the soil model and the
coupler need per cell: **`area`** (volumetric-flux conversion) and **`width`**
(the ponding/open-water shape geometry).

- `StructuredGeometry` — `area = delc[i]*delr[j]`, `width = delr[j]`,
  reproducing the legacy formulas exactly.
- `VertexGeometry` — area by shoelace on the cell polygon (no `shapely`
  dependency, which is absent in this environment), `width = sqrt(area)`.
  For square cells the two agree identically; for rectangles the areas still
  agree exactly and only the width differs — **documented approximation**,
  affecting only the ponding shape factor.
- `disv_from_structured()` — builds MODFLOW `vertices`/`cell2d` from a DIS
  grid, `icell2d = i*ncol + j`, north-first row ordering, clockwise vertices.
- `geometry_for()` — factory used by the drivers; rejects `disu` explicitly.

Wired into the three marked call sites (`Ssurf_max`/`Eosurf_max` in
`_cell_step`, the exfiltration conversion in `runMMsoil`, `MF6Coupler.area`).
`build_context(..., geom=None)` defaults to the structured provider, so every
pre-existing call site is untouched — the refactor is behaviour-preserving
(the full inherited suite stayed green throughout).

**One latent bug fixed in passing:** the exfiltration conversion used
`delr[0]*delc[0]` — a single uniform cell area for the whole grid, with a
legacy comment admitting it "assumes there is no grid refinement". It is now
per-cell. Identical on uniform grids (so La Mata results are unchanged);
correct on refined and DISV grids, where the old code would have been wrong.

## 2. DISV builder (`trunk/ppMF6/marmites_mf6.py`)

`clsMF6(..., grid='dis'|'disv')` — one class, two grids:

- `ModflowGwfdisv` with generated (or supplied) `vertices`/`cell2d`;
  `idomain`, `top`, `botm`, `k`, `k33`, `ss`, `sy`, `strt` reshaped
  `(nlay, nrow, ncol) → (nlay, ncpl)` by a single `_griddata()` helper.
- `_cellid()` returns `(lay, row, col)` or `(lay, icell2d)`, so WEL, DRN, GHB
  and UZF packagedata are written correctly for either grid **from the same
  code** — in particular the UZF **column-chaining logic (land-first
  ordering + `ivertcon`) is completely shared**, as the kickoff required.

## 3. Coupler node mapping (`trunk/marmites_coupler.py`)

`_bind()` now computes user nodes as `k*ncpl + icell2d` for DISV and
`k*nrow*ncol + i*ncol + j` for DIS, keeping the `NODEUSER` reduced-grid path
for both. The coupler also takes its cell areas from the shared geometry, so
the soil model and the exchange layer cannot disagree about them.

## 4. Raster sampling (`trunk/marmites_raster.py`)

The legacy reader asserts the raster matches the structured grid cell-for-cell
— unusable for DISV. The new module samples any ESRI ASCII raster at arbitrary
cell centres (nearest-neighbour, no `rasterio`/`shapely` needed):
`read_esri_ascii`, `sample_points`, `sample_to_cells`,
`cell_centres_structured`. Zone rasters stay integer by construction
(nearest only — interpolating class codes is meaningless). The legacy reader
is untouched and still used for DIS.

## 5. Plotting (`trunk/marmites_plot.py`)

`plot_cell_values(modelgrid, cells, values, ...)` renders a MARMITES per-cell
field on **either** grid through `flopy.plot.PlotMapView` — the same call for
DIS and DISV. Scoped as planned: a helper for new outputs; the legacy
`MARMITESplot_v3` imshow path is deliberately not ported.

## 6. Runner

`tests/run_lamata_mf6.py --grid dis|disv` (default `dis`; DISV writes to
`DataSet_LaMata/MF6_ws_disv`). Verified end-to-end on La Mata: the DISV
simulation builds with **3,900 cell2d, 3,721 vertices, 11,472 UZF objects,
1,954 wells** — the same model counts as DIS.

---

## Verification (65 tests)

New in Phase 4 (35 tests):

- **`test_grid_geometry.py` (7)** — shoelace correctness; `StructuredGeometry`
  reproduces the legacy formulas cell-by-cell; DISV-from-DIS vertices/areas;
  DIS≡DISV for square cells; documented width divergence for rectangles;
  DISU rejected.
- **`test_disv_equivalence.py` (5)** — **the primary acceptance test**:
  running the soil model through the DISV geometry reproduces the DIS results
  *exactly* (`np.array_equal` on `MM`, `MM_S`, `perc`, `etg` for every stress
  period), including the carry-over state and the exfiltration path (the most
  geometry-sensitive exchange), with mass balance still closing to <1e-6.
- **`test_mf6_disv_build.py` (5)** — DISV builds, writes and reloads through
  flopy on the real La Mata model; cellids are 2-tuples using `icell2d`;
  DIS and DISV describe the same aquifer (idomain/k/top mapped through
  `i*ncol+j`); UZF chaining preserved; polygon areas equal the DIS areas.
- **`test_raster_sampling.py` (10)** — header/bounds, nearest sampling,
  orientation (row 0 = north), NODATA/outside handling, integer zone safety,
  and **sampling the real La Mata rasters at DIS cell centres reproduces the
  legacy reader exactly** (soil zones, meteo zones, soil thickness, pond width).
- **`test_coupler_mock.py` (+3)** — DISV node mapping reads the correct `X`
  entries (unique per-node tagging so a wrong mapping cannot pass); DIS and
  DISV mappings agree on the equivalent grid; coupler area comes from the
  geometry.
- **`test_plot_helper.py` (5)** — scatter logic, and one call rendering both
  grid types.

## 7. Quadtree refinement via GRIDGEN (`trunk/marmites_gridgen.py`)

Added once the gridgen executable became available
(`C:\00MODFLOW\gridgen.1.0.02\bin\gridgen_x64.exe`).

- `build_quadtree()` wraps flopy's `Gridgen`: base `StructuredGrid` →
  refinement features → `get_gridprops_disv()`.
- `resample_inputs()` samples the MARMITES rasters onto the refined cell
  centres (via `marmites_raster`, nearest for zones).
- `RefinedModel` ties gridprops + active nodes + inputs + `VertexGeometry`.
- `clsMF6(..., cell_nodes=...)` accepts an explicit icell2d per surface cell,
  so `_cellid` is correct on a refined grid where `icell2d != i*ncol + j`.

**The refined-grid input layout.** A quadtree cell has no (row, column), but
the kernel indexes inputs as `grid[i, j]`. Rather than rewrite the kernel (and
break the frozen-oracle comparison), Phase 4 uses the degenerate layout from
§5 of the design: cells are `(cid, cid, 0, icell2d)` and inputs are
`(ncell, 1)` column vectors, so every `grid[i, j]` resolves to `grid[cid, 0]`.
`flux()` keeps its exact signature and the legacy oracle tests keep comparing
like for like.

`tests/test_refined_layout.py` (5 tests) proves this is sound **without the
gridgen binary**: the same physical model fed through the 2-D layout and
through the degenerate layout produces bit-identical results (all fluxes, all
SPs, carry-over state, with and without exfiltration), plus a guard test that
perturbing one cell's column-vector input changes only that cell.

`tests/make_quadtree_lamata.py` generates the La Mata quadtree, refining
around the DRN network (6 cells) and/or the 11 observation points
(`--features drn|obs|both`), writes `gridprops.npz` and a map of the refined
grid. Verified here up to the binary call (feature extraction and extents
checked against the model domain).

## 8. Coupled run achieved, and rejected infiltration closed

The first successful coupled MARMITES-MF6 run of La Mata (60 daily SPs,
lagged mode) required two fixes found by bisection against the live library
(`tests/diagnose_coupling.py`, `--level 0..5`):

- **UZF6 input rules.** `THTR = 0` (legal in UZF1 when `specifythtr = 0`) is
  rejected by UZF6; the value supplied in the ini is now always used.
  `EPSILON = 2.0` is outside the UZF6 range 3.5-14.0 and is clamped to 3.5
  with a loud warning -- a genuine physics difference, recorded in
  `clsMF6.eps_clamped`, that the modeller should settle explicitly.
- **WEL rates must be written to `Q`, never to `BOUND`.** MF6 6.7 exposes
  both. Writing `BOUND` (a defensive "write both" I had added) corrupts MF6's
  memory manager: the write appears to succeed and the process then dies
  inside the next `prepare_time_step()` with no message at all. Bisection
  isolated it precisely -- FINF fine, FINF+Q fine, FINF+Q+BOUND fatal.
  `BOUND` is now bound only when `Q` is absent, with a sentinel test.

Diagnostics on the successful run then exposed a genuine **water-balance
leak**: UZF rejected ~17% of the percolation MARMITES delivered
(`INFILTRATION 977` vs `REJ-INF -164`), and that water simply vanished.

**Resolved, following the published formulation.** Rejected infiltration
enters the **bottom soil layer**, exactly like groundwater exfiltration, and
is carried upward by the existing saturation-excess cascade: it fills the
soil column from below, spills into the surface store as `Rexf[0]`, and only
the excess above `Ssurf_max` becomes runoff.

This follows Francés & Lubczynski (2023) rather than a direct surface
injection. Equation 1,

    dSsurf/dt = Pe + Exf_g1 - I - Eow - Ro

admits only one subsurface->surface input, `Exf_g1` -- exfiltration arriving
*from the topmost soil layer* -- and section 2.3 states that water the
subsurface cannot accept "eventually creat[es] saturation-excess overland
flow (also known as Dunnian flow) if all the soil layers turn saturated".
An earlier implementation that injected the water straight into `Ssurf` was
rejected for this reason: it bypassed the soil column and produced runoff
without first saturating the profile.

A useful consequence of routing through the same pathway: when the subsurface
refuses water, the bottom layer's percolation is suppressed for that period
(the existing `elif EXF_ini == 0.0` guard), so MARMITES stops pushing water
into a zone that cannot take it.

`flux()` gained `REJINF_ini` (default 0.0, so the uncoupled path is
bit-identical), `step()` gained `rejinf_cell`, and the coupler reads UZF
`REJINF` each period and reports the rejected fraction. Eleven tests in
`tests/test_rejected_infiltration.py` pin the behaviour -- storage below
capacity, `Exf_g1` on saturation, Dunnian runoff when soil *and* surface are
full, equivalence with the exfiltration pathway, suppressed bottom-layer
percolation, and closure of both the Eq.-1 surface balance and the soil-column
balance.

**Open question for the modeller.** The original code carried the UZF1 budget
term `HORT+DUNN` with the comment *"should be 0"* -- i.e. the design intent
was that the unsaturated zone always accepts MARMITES' percolation. A ~17%
rejection in the MF6 run is therefore itself a signal worth investigating,
most plausibly the `EPSILON` clamp (2.0 -> 3.5, which changes unsaturated
relative permeability) or the `VKS` taken from the layer `k33`. The routing
above makes the water balance correct either way, but reducing the rejection
at source may be the better fix.

### First-run results (60 daily SPs, summer 2008, lagged)

| Quantity | Value |
|---|---|
| Outer iterations | mean 5.3, max 95 |
| Heads | 734-1132 m, mean drawdown 10.6 m over 60 days |
| Percolation | mean 3.35e-4 m/d (~122 mm/yr equivalent rate) |
| ETg | mean 7.88e-4 m/d (~288 mm/yr equivalent rate) |
| Exfiltration | 0 -- confirmed against MF6's own UZF budget (GWF ~1e-16) |

Exfiltration being zero is physically correct for these summer days, but it
means the exfiltration feedback is still unexercised; the wet season of a
full-length run will test it. The rapid drawdown (ETg drawing ~4148 L^3/T
against near-zero recharge, supplied from specific yield) should be compared
against the published NWT results before the magnitudes are trusted.

## Ready to run (binaries now unblocked)

Everything above is verified in-sandbox by build/reload round-trips, a mocked
API and bit-identity tests. The binaries are only needed for the *execution*
step, which now runs on the Windows machine. Suggested order, from the
Anaconda Prompt with the `flopy` env active, in `E:\tmp_claude_marmites\MARMITES`:

```bat
:: 1. smoke test: 60 SPs, structured grid, lagged coupling
python tests\run_lamata_mf6.py --libmf6 C:\00MODFLOW\mf6.7.0_win64\bin\libmf6.dll ^
       --mode lagged --nsp 60

:: 2. same on the DISV (vertex) grid -- should agree closely with (1)
python tests\run_lamata_mf6.py --libmf6 C:\00MODFLOW\mf6.7.0_win64\bin\libmf6.dll ^
       --mode lagged --grid disv --nsp 60

:: 3. iterative coupling (removes the one-SP lag)
python tests\run_lamata_mf6.py --libmf6 C:\00MODFLOW\mf6.7.0_win64\bin\libmf6.dll ^
       --mode iterative --relax 0.6 --nsp 60

:: 4. full length once the smoke tests pass (drop --nsp)

:: 5. quadtree grid with GRIDGEN
python tests\make_quadtree_lamata.py ^
       --gridgen C:\00MODFLOW\gridgen.1.0.02\bin\gridgen_x64.exe --level 2 --features both
```

Results: `DataSet_LaMata/MF6_ws[_disv]/_coupled_<mode>.h5` (heads, exf, perc,
ETg, outer-iteration counts per SP); the quadtree writes `gridprops.npz` and
`quadtree_grid.png`.

What to check, and what is expected:
- **DIS vs DISV (1 vs 2)** should agree closely but not bit-exactly: MARMITES'
  soil results are proven identical, while MF6 assembles a vertex grid's
  cell connections differently, so the groundwater solution can differ at
  solver tolerance.
- **lagged vs iterative (1 vs 3)** should converge together on daily SPs;
  if outer iterations hit the cap at exfiltration cells, lower `--relax`
  (0.6 → 0.4).
- If a `CouplingError: cannot resolve <var>` appears, MF6 6.7 has moved a
  memory address: send the message plus `mf6.get_input_var_names()` and the
  fix in `MF6Coupler._bind` is a one-liner.

## What Phase 4 deliberately did not do

- **Porting all legacy plotting** (explicitly out of scope).
- **DISU** (locked decision).
- **A coupled run on the refined quadtree grid.** The grid generation, the
  input resampling and the refined input layout are all implemented and
  tested; what remains is the data-preparation step of rebuilding the layered
  `idomain`/outcrop and the aquifer properties on the refined cells, which is
  best done once you have looked at the generated grid and decided on the
  refinement strategy. The DIS and DISV-from-DIS coupled runs need none of it.

## Suggested next steps

1. Run the five commands above; confirm the water balance against the
   published La Mata results.
2. Inspect `quadtree_grid.png`, settle the refinement strategy, then resample
   the full input set onto it (`marmites_gridgen.resample_inputs` /
   `RefinedModel`) for a coupled run on the refined grid.
3. Fold the MF6 + grid options into `marmites_config.py` (TOML) now that the
   CLI surface has settled, and do the driver/export split deferred since
   Phase 2 — both are natural once the first coupled runs are green.
