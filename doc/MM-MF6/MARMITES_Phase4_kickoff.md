# MARMITES — Phase 4 kickoff (handoff to Opus)

**Written:** 2026-07-21, end of Phase 3 (Fable 5).
**Goal of Phase 4:** unstructured-grid support via **DISV** (vertex grid) —
grid generation, input sampling, geometry plumbing, DISV variant of the MF6
builder, plotting. Phase 3 delivered the MF6/DIS backend and the API coupler;
Phase 4 generalizes the spatial dimension. **DISU is out of scope** (locked
decision — the outcrop-layer concept requires layered grids).

Read first, in order: `MARMITES_code_review.md` (§5 is the DISV design, §6 the
plan), `MARMITES_Phase3_kickoff.md` (interface contract + environment traps —
they all still apply), `MARMITES_Phase3_report.md` (what exists now). Then this.

## State you inherit

- 30/30 tests green (`python3 -m pytest tests/ -q --ignore=tests/validate_lamata.py`).
- MF6/DIS backend: `trunk/ppMF6/marmites_mf6.py` (`clsMF6`), API coupler
  `trunk/marmites_coupler.py` (lagged + iterative, mock-tested), runner
  `tests/run_lamata_mf6.py`.
- **Pending outside your control:** the first real coupled binary run on the
  user's machine is blocked by antivirus (libmf6.dll load → WinError 87); IT
  unblocking is awaited. Do NOT let Phase 4 depend on it. If the user reports
  a `CouplingError: cannot resolve <var>` or oscillation when the DLL is
  freed, that is Phase-3 debugging: fix the address candidates in
  `MF6Coupler._bind` / lower `relax`; the mock tests define correct behavior.

## The one contract that makes DISV cheap

Everything spatial flows through the **cell list**: MMsoil iterates
`ctx.cells = [(cid, i, j, node), ...]` and the coupler requires
`clsMF6.surf_cells` to be in the same order (asserted at construction).
For DISV the *only* semantic change is: a map-view cell is an `icell2d`
instead of `(i, j)`. Strategy (§5, locked): keep the cell list as the single
source of truth and make `(i, j)` degenerate — for DISV build cells as
`(cid, icell2d, 0, icell2d)`-style with a `CellGeometry` provider supplying
what today comes from `delr/delc` indexing.

## Work items (suggested order, each with its own tests)

1. **`CellGeometry` abstraction.** New small module (e.g.
   `trunk/marmites_grid.py`): per-cell `area(cid)`, `width(cid)` (for the
   Ssurf/Eosurf geometry), plus `ncpl`, `nlay`, `top/botm(cid)`. Two
   implementations: `StructuredGeometry(delr, delc, i_arr, j_arr)`
   (reproduces today's `delr[j]*delc[i]` exactly — regression-tested against
   current behavior) and `VertexGeometry(flopy VertexGrid)` (area from
   `modelgrid.geo_dataframe` or `get_cell_vertices`; width = sqrt(area) —
   document this approximation for the Ssurf shape factor).
   Then replace the three marked call sites:
   - `MARMITESsoil_v3._cell_step`: `Ssurf_max`, `Eosurf_max` (search
     `TODO Phase 4`), and the exf area conversion in `runMMsoil`;
   - `marmites_coupler.MF6Coupler.__init__`: `self.area`.
   Guardrail: after this refactor the entire existing suite must still pass
   bit-identically (StructuredGeometry is a pure refactor).

2. **DISV variant of `clsMF6`.** `ModflowGwfdisv` (vertices/cell2d from a
   flopy `VertexGrid`), `cellid = (lay, icell2d)` for WEL/DRN/GHB/UZF
   packagedata, `idomain (nlay, ncpl)`. Keep one class with a `grid`
   parameter or a sibling class — prefer whichever keeps the UZF
   column-chaining code (land-first ordering + `ivertcon`) shared; that logic
   is grid-independent once cells are 1-D. The coupler needs a matching
   node-mapping change in `_bind`: user node = `k*ncpl + icell2d` (DISV)
   instead of `k*nrow*ncol + i*ncol + j`; keep the NODEUSER reduced-grid path.

3. **Grid generation for La Mata.** Quadtree via flopy's `Gridgen` wrapper
   needs the gridgen executable (same download channel as mf6 — likely
   blocked in the sandbox; check first). If unavailable, generate the DISV
   demonstration grid WITHOUT gridgen: a rectangular vertex grid is directly
   constructible from the DIS geometry (each (i,j) → one cell2d) — that is
   the correctness test anyway (DISV-from-DIS must reproduce DIS results);
   true quadtree refinement can then be produced on the user's machine.
   **Primary acceptance test of Phase 4: La Mata rebuilt as DISV-from-DIS
   gives the same MMsoil outputs as the DIS run** (same cells, same areas).

4. **Input sampling replacing ESRI-ASCII-only readers.** Add a path that
   samples rasters onto an arbitrary modelgrid: `flopy.utils.Raster
   .resample_to_grid` (needs `rasterio` — check installability in the
   sandbox; if blocked, implement nearest-neighbour sampling from the ASC
   grids directly — they are plain arrays with a georeference, and the La
   Mata case is axis-aligned). Zones (int grids) must use nearest/mode, not
   interpolation. Keep `convASCIIraster2array` untouched for DIS.

5. **Plotting.** Replace the imshow-based `plotLAYER` path with
   `flopy.plot.PlotMapView` for maps (works for DIS and DISV alike). Scope
   this pragmatically: a new small plotting helper used by new outputs is
   enough; porting all of `MARMITESplot_v3` is NOT required for Phase 4
   acceptance.

## Do NOT touch

- `tests/legacy/MARMITESsoil_v3_legacy.py` (frozen oracle) and the physics of
  `flux()`/`_cell_step` (geometry call sites only).
- The coupler exchange semantics (mock tests define them). If DISV needs a
  different pointer shape, extend `_bind`, don't fork the exchange logic.
- The 30 green tests. Add, don't break; run the suite after every work item.

## Environment reminders (all verified the hard way)

- Sandbox Python 3.10 (target 3.12): no `tomllib` (tomli installed), strip
  every `readFile` token, use `os.path.abspath`, never read the multi-GB h5
  files whole (slice `[0:n]`), bash timeout 45 s, git commits blocked on the
  mount (files persist; user snapshots locally).
- flopy 3.10.0 + modflowapi 0.2.0 installed. **No MF6/gridgen binaries in the
  sandbox** (proxy blocks GitHub CDN + conda) — anything needing execution of
  MODFLOW tools goes through `tests/run_lamata_mf6.py`-style entry points for
  the user's machine (mf6 at `C:\00MODFLOW\mf6.7.0_win64\bin\`, pending AV
  unblock).
- User context: Windows, Miniconda env `flopy` (needs admin rights to install
  packages), runs via Spyder or Anaconda Prompt — prefer prompting them to
  use the Anaconda Prompt; Spyder run-configs mangle CLI args.

## Definition of done for Phase 4

1. `CellGeometry` in place, suite green, DIS results unchanged.
2. `clsMF6` builds a loadable DISV simulation (flopy round-trip test, like
   `test_mf6_build.py`).
3. DISV-from-DIS equivalence test passes (same MMsoil outputs as DIS on a
   truncated La Mata run through the mock or file path).
4. Runner gains `--grid disv`; docs updated; Phase-4 report written in the
   style of the Phase 1-3 reports.
