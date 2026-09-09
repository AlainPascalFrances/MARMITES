# MARMITES — Phase 2 completion report

**Date:** 2026-07-21
**Model:** Claude Opus 4.8
**Scope:** Phase 2 of `MARMITES_code_review.md` §6 — internal refactor: the grid-agnostic `(nper, ncpl)` data model, the per-stress-period `step()` interface (the Phase-3 API hinge), TOML configuration replacing positional ini parsing, and shared-module extraction.

---

## Result

The soil model is now grid-agnostic and exposes the per-SP `step()` interface that the MODFLOW 6 API driver will call in Phase 3. Positional ini parsing is replaced by a validated TOML schema with a legacy converter. The test suite grew from 9 to **20 tests, all passing**; core modules are clean under ruff's error classes.

---

## Phase 2 — grid-agnostic cell model + `step()` (the Phase-3 hinge)

`MARMITESsoil_v3.py` was restructured around a 1-D **active-cell list** instead of nested `(nrow, ncol)` loops:

- `build_cell_list(cMF)` returns active cells as `(cid, i, j, node)`, where `node = i*ncol + j`. For a structured grid this is the flat index; for a DISV grid (Phase 4) `node` becomes the `icell2d` and the `(i, j)` fields collapse onto it. This is the concrete `ncpl` model — a structured grid is the special case `ncell = count(outcropL > 0)`.
- `init_state(ctx)` holds the inter-SP carry-over (`Ssoil_ini`, `Ssurf_ini`) as `(ncell,)` / `(ncell, nslmax)` arrays — the state is now indexed by cell id, not by row/col.
- `step(ctx, n, tstart_MF, heads_cell, exf_cell, state)` advances one stress period over all cells and returns per-cell `MM`, `MM_S`, `perc`, `etg` arrays, mutating `state`. **This is the interface the Phase-3 API driver calls**: heads and exfiltration come in as `(ncell,)` arrays (from MODFLOW memory), percolation and ETg go out as `(ncell,)` arrays (to UZF `finf` and WEL `Q`). No HDF5 needed on that path.
- `runMMsoil(...)` keeps its old signature (driver-compatible) but is now a thin loop: it reads structured heads/exfiltration from the MODFLOW HDF5, gathers them to per-cell arrays, calls `step()`, and scatters the per-cell results back to the structured `(nper, nrow, ncol)` datasets the export/plot code consumes.

The per-cell physics kernel (`flux()`) is unchanged — the refactor is purely structural, which the tests confirm.

### On full SIMD vectorization — an explicit engineering decision

The review’s Phase-2 line item said “vectorize the inner loop.” I iterate the **1-D cell list** (the grid-agnostic flatten that unblocks DISV and the API), but I did **not** rewrite `flux()` to operate on arrays of cells simultaneously. Reasoning: the per-cell kernel is deeply sequential (inter-layer percolation cascade, per-vegetation PT consumption, upward saturation cascade, root-ordered Tg with head feedback), so true across-cell SIMD is a high-risk rewrite of the core physics — and the dominant cost was already removed in Phase 1 (Decimal elimination). The cell-list loop delivers the structural goal (grid independence, the `step()` seam) at zero physics risk. If profiling on a large DISV grid later shows the Python cell loop dominates, across-cell vectorization can be added behind the same `step()` interface without touching callers. Flagged, not silently skipped.

## Phase 2 — TOML configuration (`marmites_config.py`)

Replaces the fragile `l += 1` positional parsing of the MM ini (where one missing line silently shifted every later parameter):

- `MMConfig` dataclass with typed fields, defaults, and a `validate()` method (checks hydrologic-year month, `plt_WB_unit`, `MMsoil_yn ∈ {-1,0,1}`, irrigation consistency, etc.).
- `load_mm_config(path)` reads TOML (stdlib `tomllib` on 3.12, `tomli` fallback for older test interpreters) into a validated `MMConfig`.
- `convert_ini_file(ini, toml)` / `legacy_ini_to_config(values)` migrate existing datasets automatically.
- Picard-loop fields (`convcrit`, `convcritmax`, `ccnum`, `MF_yn`, `MF_lastrun`), obsolete since Phase 1, are parsed but quarantined into a `deprecated` sub-table and ignored.

Converting the real `DataSet_WRR/__inputMM_v3.ini` surfaced a concrete instance of the very fragility this fixes: the driver’s parser expects `maxYearsTickTrimester/Semester`, but that dataset’s ini does not contain them — a positional drift that would misalign every subsequent field. The converter tolerates their absence (they are cosmetic tick controls with defaults); the finding is documented in the code.

## Phase 2 — shared-module extraction (`marmites_indices.py`)

The `INDEX_MM` / `INDEX_MM_SOIL` output-column maps were a literal block in the driver and re-declared in tests. They now live in one module imported by the driver and the tests — one source of truth.

### On the driver/export physical split — deliberately deferred

The review placed a `config / driver / export` split in Phase 2. I extracted **config** and **indices** (stable, reused, testable). I deliberately **did not** physically split the ~1,800-line export/plotting tail of `startMARMITES_v3.py`, for two reasons: (1) it cannot be executed end-to-end yet (no MF6 backend, no legacy NWT runtime), so a blind extraction of tightly-coupled plotting code could introduce regressions no test would catch; (2) Phase 3 rewrites the driver core for the API coupling, so splitting it now is partly wasted effort. Recommendation: perform the driver/export split as part of the Phase-3 rewrite, when the code is runnable and the API driver is being written anyway.

---

## Verification

- `py_compile` on all core modules + the two new modules: OK.
- **20 tests pass** (`tests/`): 9 flux regression (legacy-vs-new), 5 `runMMsoil`/`step` integration (mass-balance closure, inactive-cell masking, and `step()` reproducing `runMMsoil` percolation/ETg exactly), 6 config (real-dataset conversion, TOML round-trip, validation).
- `ruff --select F,E9,E722,E711,E721` on core: **All checks passed**.

### A physics finding worth recording

The integration tests showed the **soil-column mass balance closes to ~1e-14 for any SP length**, but the **surface mass balance closes only for daily (`perlen=1`) stress periods**. This is inherent to MARMITES’ SP-averaging (flux() computes one representative day then rate-averages), and it is exactly why the review chose **daily stress periods as the Phase-3 default**. The tests encode both facts: soil-MB closure is asserted for mixed multi-day SPs; full (soil+surface) closure is asserted only for daily SPs.

---

## Handed to Phase 3

- `step(ctx, n, tstart_MF, heads_cell, exf_cell, state)` is ready to be driven by `modflowapi`: supply `heads_cell`/`exf_cell` from MODFLOW memory each SP (lagged mode) or each outer iteration (iterative mode), take `perc`/`etg` back to UZF6/WEL. `build_context()` / `init_state()` set up the static context and carry-over once.
- Dry-cell handling still uses the `hdry`/`hnoflo` sentinels inside `_cell_step` (marked with a NOTE); Phase 3 replaces them with `idomain` + `h < botm`.
- Cell geometry (`Ssurf_max`, `Eosurf_max`, exf area conversion) still uses `delr/delc` (marked TODO Phase 4); DISV supplies per-node area then.
- `marmites_config.py` is the config contract; extend with an MF6-package section in Phase 3 (the NWT ini schema was intentionally not modelled).
- Style-class ruff findings (UP031 printf, B007/B905) remain across the driver/plot tail — cosmetic, safe to sweep alongside the Phase-3 driver rewrite.

## For the next model

Read `MARMITES_code_review.md` (contract), then `tests/` (behavioral spec: 20 tests must stay green). `tests/legacy/MARMITESsoil_v3_legacy.py` is the frozen v0.3 oracle — do not modify. The Phase-3 coupling driver is the subtle part (iteration-embedded exchange, under-relaxation, dry-cell semantics under MF6 Newton); build it against `step()` and validate on La Mata.
