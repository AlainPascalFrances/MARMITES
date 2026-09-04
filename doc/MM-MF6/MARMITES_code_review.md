# MARMITES — Critical code review and MF6/unstructured-grid migration analysis

**Scope:** `E:\tmp_claude_marmites\MARMITES\trunk` (~18,100 lines of Python), reviewed against Francés & Lubczynski (2023, *Front. Water* 5:1055934), the MF6 documentation set (TM 6-A55, TM 6-A57, mf6io 6.x), the MODFLOW API paper (Hughes et al. 2022), the FloPy unstructured-workflows paper (Hughes et al. 2024) and GRIDGEN OFR 2014-1109.

**Date:** 2026-07-21

---

## 1. What the code is today

| Module | Lines | Role |
|---|---|---|
| `startMARMITES_v3.py` | 3,331 | Monolithic driver script: reads ini files, orchestrates MMsurf → MF init → MM↔MF Picard loop → export/plots |
| `ppMF_FloPy/ppMODFLOW_flopy_v3.py` | 1,503 | `clsMF`: parses a custom MODFLOW ini, builds MF2005/NWT packages with `flopy.modflow`, runs MF, harvests heads + cbc into HDF5 |
| `MARMITESsoil/MARMITESsoil_v3.py` | 731 | Cell-by-cell soil-zone water balance (interception → Ssurf/Ro → Esoil/Tsoil → Rp), Eg (Shah et al. 2007), phenomenological Tg |
| `MARMITESutilities/MARMITESprocess_v3.py` | 781 | ESRI ASCII import, soil/obs/time-series parsing, SP→daily post-processing |
| `MARMITESsurf/` | ~1,940 | FAO56-PM PET/PE/interception preprocessing |
| `MARMITESplot/` | ~3,400 | Custom plotting (imshow-based maps, time series, Sankey) |
| Rest (`MM_XLSstuff`, `pyEARTH1D`, `Tg*.py`, `*_BAK*`) | ~5,000 | One-off scripts, abandoned experiments, wxPython GUI, backups under version control |

Coupling scheme (matches the paper): MM computes percolation (`perc`) and ETg per cell; these are passed to MODFLOW as UZF1 `finf` and per-cell negative WEL; MF returns heads and UZF `SURFACE LEAKAGE` (exfiltration), which feed the next MM iteration; loop repeats until average/max head change < convergence criteria. Exchange is entirely file-based (HDF5), and each iteration reruns the full transient MF model.

The scientific core is sound and published. The implementation, however, is a 2010–2013 research code that has been minimally touched since: it targets MODFLOW-NWT/MF2005 through an old `flopy.modflow` API, assumes a regular structured grid everywhere, and carries substantial dead weight.

---

## 2. Bugs found (correctness, before any migration)

These are live defects in code paths you rely on:

1. **`MARMITESsoil_v3.py:192` — `Rexf_tmp /= perlen` on a Python *list*.** `Rexf_tmp` is built as a list of `Decimal` (lines 170–172). `list /= int` raises `TypeError`. This executes whenever `EXF_ini > 0`, i.e. the groundwater-exfiltration feedback — a headline feature of the paper — crashes (or is being silently absorbed by an enclosing bare `except`). Must be an array, or the division applied per element.
2. **`MARMITESsoil_v3.py:530` — `if CROP_tmp != None:`** where `CROP_tmp` can be a numpy array → `ValueError: truth value of an array is ambiguous`. Same pattern at line 548 (`if LAIveg_tmp[v] > 1.0E-5:` on a time-series row). The irrigation branch cannot run on modern numpy.
3. **`ppMODFLOW_flopy_v3.py:1467, 1484` — `for i in range(self.row)`**: attribute is `self.nrow`; the fallback path for large arrays raises `AttributeError`.
4. **`ppMODFLOW_flopy_v3.py:314` — `self.versionsys.exit()`**: garbled line, crashes if ever reached.
5. **`startMARMITES_v3.py:22` — `if mpl.get_backend != 'agg':`** compares the *function object* to a string; always True. Harmless here, but symptomatic.
6. **`MARMITESsoil_v3.py:113–131` — `global rp_tmp` in `perc()`**: if neither branch assigns (e.g. `s_tmp` is NaN, all comparisons False), the function silently returns the *previous cell's* percolation. `evp()` has the mirrored problem (`UnboundLocalError`). NaNs propagate as wrong numbers instead of errors.
7. **`wel_dum` early-`break` logic (`ppMODFLOW_flopy_v3.py:1149–1203`)**: when the first active cell has zero flux, the triple `break` leaves all remaining stress periods with a single dummy well; intent (dummy WEL to keep the package alive) is obscured and fragile.
8. **78 bare `except:` clauses in `startMARMITES_v3.py` alone** (17 in `ppMF`, 12 in `MARMITESprocess`). They convert every genuine bug (including items 1–3 above) into a generic "Type error in the input file" message. This is the single most damaging pattern in the codebase — it has probably been hiding defects for a decade.
9. **`global` state across functions** (`MARMITESprocess_v3.py:104,218,706`; `MARMITESsoil_v3.py:379`): read-order dependencies invisible to the caller; `convASCIIraster2array` leaves `NODATA_value`, `fin`, etc. as module globals.
10. Cosmetic but telling: unreachable `del` after `return` (several functions), duplicated line `Ssurf_ini_tmp_array = ...` (`MARMITESsoil_v3.py:401–402`), `raise sys.exit()` idiom, `*_BAK*.py` and `TgNEW_v2..v5.py` under version control.

---

## 3. Python 3.12 compatibility assessment

The core chain (`startMARMITES_v3` → `ppMF` → `MMsoil` → `MMproc` → plots) has already been mechanically converted to Python-3 syntax and will *parse* under 3.12. It will not run reliably. Findings:

**Hard blockers**
- `MARMITESutilities/MM_XLSstuff/*`, `MFarray2ascarray.py`, and others still use Python-2 `print` statements → `SyntaxError`. Recommendation: delete or archive the whole `MM_XLSstuff` and `pyEARTH1D` (wxPython 2.x GUI) trees rather than porting them.
- `np.float`, `np.int` (removed in numpy ≥ 1.24) in several utility scripts.
- Elementwise `!= None` array comparisons (28 occurrences across the core) — deprecated/broken semantics; must become `is not None`.
- `flopypth = os.path.dirname(sys.executable) + r'\Lib\site-packages\flopy'` + `sys.path.append` (`ppMF:23–25`): remove; it can shadow the installed flopy and is Windows-only.

**The Decimal problem.** 63 `Decimal(...).quantize(...)` call sites, most of them inside the innermost per-cell, per-stress-period loop of `MMsoil.flux()`. Three separate objections:

1. *Performance*: `Decimal` arithmetic is 10–100× slower than `float64`, in a loop executed `nrow × ncol × nper` times. This, together with the pure-Python double loop, is why runs take hours.
2. *Correctness under Python 3*: `Decimal * float` raises `TypeError` (Python 2 silently compared/coerced in some paths). The code survives only through defensive `float(...)`/`Decimal(...)` wrappers applied inconsistently — e.g. `MARMITESsoil_v3.py:248` mixes `Decimal`, `np.float32` and `float` in one expression. Any refactor will keep tripping over this.
3. *It doesn't buy anything*: quantizing to 1e-5 mm does not remove water-balance error; it just rounds it. Tolerance-based comparison (`np.isclose`) on `float64` is the right tool.

Recommendation: eliminate `Decimal` entirely in phase 1. This is the highest-leverage single change in the port.

**Soft issues**: bytes keys (`b'iP'`, `b'SURFACE LEAKAGE'`) as HDF5/dict indices — a Python-2 str/bytes leftover; `mpl.dates.datestr2num(datetime.today().isoformat())` idioms (work, but replace with `datetime`/`pandas`); `ffmpeg.exe` hard-coded (Windows-only animation); `os.path.expanduser('~'), 'Desktop'` + `startMM_fn.txt` bootstrap (replace with a CLI argument, e.g. `python -m marmites run path/to/config.toml`).

**Environment to target**: Python 3.12, numpy ≥ 2.0, flopy ≥ 3.9, h5py ≥ 3.11, matplotlib ≥ 3.9. None of the current core survives numpy 2.x untouched, so pinning and a regression test come first (§6).

---

## 4. MODFLOW-NWT → MODFLOW 6

### 4.1 Where the NWT coupling actually lives

The MF-specific surface is well contained — essentially all in `ppMODFLOW_flopy_v3.py` plus the cbc-record names used in `startMARMITES_v3.py:948–991` and `MARMITESprocess_v3.procMF()`. That is good news: the soil physics does not need to change.

| Current (NWT/flopy.modflow) | MF6 equivalent | Impact |
|---|---|---|
| `Modflow(version='mfnwt')` | `MFSimulation` + `ModflowGwf` (Newton via `newtonoptions`) | rewrite of model construction |
| `ModflowDis` (nper/perlen in DIS) | `ModflowTdis` (simulation level) + `ModflowGwfdis`/`disv`/`disu` | time moved out of DIS |
| `ModflowBas` (`ibound`, `strt`, `hnoflo`) | `idomain` (in DIS) + `ModflowGwfic` | **no ibound<0 constant heads** — needs CHD package; `hnoflo`/`hdry` sentinel values disappear (see below) |
| `ModflowUpw`/`Lpf` (`laytyp`, `hk`, `vka`, `ss`, `sy`, `iphdry`) | `ModflowGwfnpf` (`icelltype`, `k`, `k33`) + `ModflowGwfsto` (`ss`, `sy`, steady/transient flags per period) | storage/steady-state control moves to STO |
| `ModflowNwt` / `ModflowPcg` | `ModflowIms` | straightforward |
| `ModflowWel` + `options=['SPECIFY 0.05 iunitramp']` | `ModflowGwfwel` with `auto_flow_reduce` | the whole `iunitramp`/`flopy.pakbase` hack (`ppMF:1295–1305`) goes away |
| `ModflowDrn`, `ModflowGhb` (`[l,i,j,...]`) | same packages, `cellid` tuples: `(lay,row,col)` DIS / `(lay,icell2d)` DISV / `(node)` DISU | mechanical, but must be grid-agnostic |
| `ModflowUzf1` | `ModflowGwfuzf` (UZF6) | biggest package change, see 4.2 |
| `ModflowOc` (`ihedfm`, spd dict) | `ModflowGwfoc` (`saverecord`) | simple |
| `flopy.utils.HeadFile` | unchanged (or `gwf.output.head()`) | trivial |
| `CellBudgetFile`: `'FLOW RIGHT FACE'`, `'FLOW FRONT FACE'`, `'FLOW LOWER FACE'` | single `'FLOW-JA-FACE'` on the connection list | **structural change**: face flows are no longer per-axis arrays; `procMF` FRF/FFF/FLF datasets must be redesigned (or dropped — check whether anything downstream truly needs face flows rather than budgets) |
| cbc `'STORAGE'` | `'STO-SS'` + `'STO-SY'` (two records) | update index logic |
| `itmuni`/`lenuni` integer codes + `conv_fact` from `lenuni` | `time_units`/`length_units` strings in TDIS/DIS | update the `conv_fact` block (`startMARMITES_v3.py:506–527`) |

**hnoflo/hdry deserve emphasis.** MM currently detects dry cells by `abs(h - hdry) < 1e-5` (`MARMITESsoil_v3.py:556`) and masks arrays with `hnoflo ± 0.09`. Under MF6 with Newton, dry cells stay active with heads below cell bottom; with standard formulation they can convert. There is no `hdry` sentinel. The dry-cell logic must be rewritten as `h < botm[cell]` tests, and the inactive-cell masking driven by `idomain`, not magic values. This touches MMsoil, the plotting, and every masked-array operation.

### 4.2 UZF1 → UZF6, or drop UZF?

MM uses UZF1 for exactly three things: (a) deliver `finf` (MM percolation) to the water table with unsaturated-zone delay, (b) return `SURFACE LEAKAGE` (exfiltration) to MM, (c) `UZF RECHARGE` for reporting Rg. GW-ET inside UZF is optional (`ietflg`), and runoff routing (`irunflg`) is unused.

If you keep UZF: UZF6 differs materially — per-cell UZF objects (`packagedata` with `landflag`, `ivertcon` for vertical stacking, `surfdep`, `vks`, `thtr/thts/thti/eps`), period data for `finf`, and budget records renamed: `UZF-GWRCH` (recharge), `UZF-GWD` (groundwater discharge to land surface = your SURFACE LEAKAGE), `UZF-GWET`. The `nuzgag` gage machinery is replaced by UZF observations (OBS6 utility) — the whole `row_col_iftunit_iuzopt` / `uzf_filenames` block (`ppMF:1087–1109`) is rewritten in ~10 lines of `obs` dicts. `HORT+DUNN` no longer exists as a record; rejected infiltration appears as `UZF-INF`-vs-applied difference or via the MVR budget if routed.

**DECIDED (2026-07-21): UZF6 always.** The RCH+DRN alternative was considered and rejected: the unsaturated travel-time delay between soil bottom and water table must be preserved, because MARMITES will be applied in areas with deep water tables (La Mata is only the pilot). Consequence: the thts/thtr/thti/eps initialization logic (`ppMF.array_ini:543–597`) must be ported and cleaned, `ivertcon` vertical stacking configured for multi-layer UZF columns, and the `landflag` cells aligned with `outcropL`.

### 4.3 Coupling architecture

**DECIDED (2026-07-21): remove the Picard scheme; couple through the MODFLOW 6 API (`modflowapi`/`xmipy`). MODFLOW-NWT support is deleted, not ported — no initial NWT run remains in the workflow.**

The driver becomes a single march over stress periods with two coupling modes sharing one `step()` interface and one pointer-exchange layer:

- **`lagged` mode (baseline, implement first).** For each SP *n*, MM computes the soil water balance using heads and UZF exfiltration from the end of SP *n−1*, writes `finf` (UZF6) and ETg (WEL) into MF6 memory, then advances MF6 through SP *n*. Simple, robust, ideal for debugging. Lag error bounded by SP length.
- **`iterative` mode (default once validated).** MM is re-evaluated inside MF6's outer-iteration loop: the XMI sequence `prepare_solve` → repeated `solve()` (one Newton outer iteration per call, `iteration_start` callback in `modflowapi.run_simulation`) → `finalize_solve` allows recomputing `finf`/ETg each outer iteration from the *current head iterate* of SP *n*. MM then acts as a head-dependent boundary package inside MF6's nonlinear solve — the iMOD MetaSWAP–MF6 coupler architecture. This eliminates the SP−1 lag entirely; Eg/Tg/exfiltration become consistent with end-of-SP heads, which matters most for deep-water-table applications where ETg is strongly head-sensitive. Requirements: the vectorized MM step of Phase 2 (MM runs once per outer iteration, 5–20× per SP — prohibitive with the current per-cell Decimal loop, negligible when vectorized) and under-relaxation of the exchanged fluxes (damping ~0.5–0.7 on ΔFINF/ΔQ between iterations) to prevent wet/dry oscillation at exfiltration cells, alongside MF6 Newton backtracking.

Common implications:

- Results will differ from the published Picard runs — La Mata validation becomes a plausibility/tolerance comparison, not a regression match.
- Exchange is pointer-based: `get_value_ptr` on UZF `FINF`, WEL `Q` (bound with `auto_flow_reduce`), heads `X`, and the UZF groundwater-discharge/rejected-infiltration arrays. The HDF5 shuttle of 4-D arrays between MM and MF disappears; HDF5 remains only as MARMITES' own output store.
- **Time discretization can be simplified.** The rainfall-based SP aggregation (`ppMFtime`, ~380 lines) existed to shorten MODFLOW runs repeated by the Picard loop. With a single API-driven march, daily stress periods become the sensible default — removing the flux-averaging approximation, most of the aggregation bookkeeping, and (in lagged mode) reducing the lag to one day. Keep `perlenmax` aggregation as an option for very long simulations.
- The first-SP steady-state trick (`dum_sssp1`) survives as an ordinary initial steady-state SP advanced through the API before the transient march.
- MF6 package files are written once at startup (flopy.mf6); only memory values change per SP/iteration.

---

## 5. Unstructured grid readiness

Blunt assessment: the code is structurally welded to a regular grid. `(nrow, ncol)` indexing appears ~200× across the core; cell area is `delr[j]*delc[i]` (with `reggrid==1` shortcuts assuming uniform spacing, e.g. `exf4MM ... /(delr[0]*delc[0])`, `MARMITESsoil:392`); *all* spatial input is ESRI ASCII rasters whose cellsize is asserted equal to the MF grid (`MARMITESprocess:129`); Ssurf geometry uses `delr[j]` directly (`MARMITESsoil:513–514`); observation wells are located by `pp_indexfromcoordinates` row/col arithmetic; every map plot is an `nrow×ncol` image.

Two decisions make this tractable:

1. **Target DISV, not DISU.** MARMITES' central spatial concept is the *outcropping layer* (`outcropL`): each map-view cell has one soil column sitting on the uppermost active layer. DISV preserves the layer × cell2d structure, so `outcropL`, `iuzfbnd`-style tops, and the soil-column-per-surface-cell concept map one-to-one onto `icell2d`. Pure DISU abandons layers and would force a redesign of the outcrop logic for zero benefit — quadtree/Voronoi grids are exactly what DISV + GRIDGEN/triangle/voronoi (FloPy `Gridgen`, `VertexGrid`) produce. Support DIS and DISV; treat DISU as out of scope unless you have a concrete nested-grid use case.
2. **Refactor MM's data model from `(nper, nrow, ncol, ...)` to `(nper, ncpl, ...)`.** A structured grid is then the special case `ncpl = nrow*ncol`. This is the same refactor needed to vectorize the per-cell double loop (`for i in range(nrow): for j in range(ncol):` → operations over 1-D node arrays), so the unstructured requirement and the performance problem have a common solution. Cell geometry (area, x, y, top) comes from `flopy.discretization.StructuredGrid`/`VertexGrid` (`modelgrid.geo_dataframe`, `cell_thickness`, etc.) instead of `delr/delc`.

Consequences elsewhere: spatial inputs move from ASCII rasters to raster/shapefile sampling onto the model grid (`flopy.utils.GridIntersect`, or `flopy.utils.Raster.resample_to_grid`); plotting moves to `flopy.plot.PlotMapView` (works identically for DIS and DISV, replacing most of `MARMITESplot_v3.plotLAYER`); HDF5 layouts become `(nper, nlay, ncpl)`.

---

## 6. Agreed action plan (updated 2026-07-21)

Decisions locked in: (1) fix bugs 1–10 of §2; (2) remove `Decimal` (§3); (3) **UZF6 always** — no RCH+DRN variant (§4.2); (4) **remove the Picard scheme, couple via the MODFLOW 6 API** — `lagged` mode first (MM at SP *n* uses heads from SP *n−1*), `iterative` mode (MM inside MF6's outer-iteration loop) as the optimized default (§4.3); (5) **DISV**, not DISU (§5); (6) **MODFLOW-NWT support is deleted**, no NWT baseline run — regression scope is the MM soil module itself.

Ordered to keep the model verifiable at every step; each phase gates the next.

**Phase 0 — Safety net (small effort, do first).**
Put the trunk in git (if not already); archive `MM_XLSstuff`, `pyEARTH1D`, `*_BAK*`, `TgNEW_v*` out of the package; create `pyproject.toml` with pinned deps (Python 3.12, numpy ≥ 2, flopy ≥ 3.9, h5py, modflowapi/xmipy); add `ruff` + `pytest`. Since the NWT path is deleted, the regression reference is MM-only: build soil-column test cases (single cell, a few soil types, prescribed heads/exfiltration series covering shallow/deep WT, EXF>0, saturation) and record `flux()` outputs as the reference for Phases 1–2. La Mata remains the end-to-end plausibility benchmark from Phase 3 on.

**Phase 1 — Bug fixes + Decimal removal + Python 3.12 port (actions 1 and 2).**
Fix §2 bugs 1–10; remove all `Decimal` usage (float64 + `np.isclose` tolerances); replace bare `except:` with typed exceptions (`MarmitesError`); kill `global`s; bytes→str keys; remove the flopy sys.path hack; CLI entry point instead of Desktop `startMM_fn.txt`. **Delete the NWT/mf2005 code paths** (`version != 'mfnwt'` branches, `ModflowNwt/Upw/Lpf/Pcg` construction, `runMF`, the initial-MF-run block in the driver, and the MM↔MF convergence-loop scaffolding — `h_diff` bookkeeping, `plt_ConvLoop`, `dum_sssp1` array-shuffling). Acceptance: soil-column tests match Phase-0 reference within float tolerance. Expect a large speedup from Decimal removal alone.

**Phase 2 — Internal refactor.**
Replace the positional-line ini parsing (the `l += 1` pattern, ~400 lines, where one missing line silently shifts every subsequent parameter) with TOML + validation (dataclasses or pydantic). Split `startMARMITES_v3.py` into `config / driver / export` modules. Flatten MM arrays to `(nper, ncpl)` and vectorize `runMMsoil`'s inner loop — this is also the DISV enabler. Critically for Phase 3: restructure the MM driver so the soil model exposes a clean per-SP step interface (`step(n, heads, exf) -> (finf, etg, state)`) instead of the monolithic all-SP `runMMsoil`. Acceptance: same regression test.

**Phase 3 — MF6 backend with UZF6, API coupling (actions 3 and 4; structured DIS).**
Three sub-steps:
- *3a — model construction.* New `clsMF6` on `flopy.mf6` per §4.1: TDIS/DIS/IC/NPF/STO/IMS (Newton), WEL with `auto_flow_reduce`, DRN/GHB with grid-agnostic `cellid`, UZF6 (`packagedata` with `landflag`/`ivertcon` aligned to `outcropL`, OBS6 instead of gages). New dry/inactive logic: `idomain` masking, `h < botm` dryness, no `hnoflo`/`hdry` sentinels. Daily stress periods as default (SP aggregation optional). Validate the generated files by one conventional full run.
- *3b — coupling driver, `lagged` mode.* Per-SP march through `modflowapi`: steady initial SP; then for each transient SP *n*, MM `step()` with heads/UZF-GWD from SP *n−1*, write `FINF` and WEL `Q` via memory pointers, advance MF6. Post-processing reads MF6 budget names (`UZF-GWRCH`, `UZF-GWD`, `UZF-GWET`, `STO-SS`/`STO-SY`, `FLOW-JA-FACE`).
- *3c — coupling driver, `iterative` mode.* MM `step()` re-evaluated at each MF6 outer iteration (`iteration_start` callback / repeated `solve()`), with under-relaxation on the exchanged fluxes. ~50 lines on top of 3b given the shared interface; becomes the default after validation against 3b (differences should shrink as SPs shorten).
Acceptance: La Mata heads RMSE and water-balance components within documented tolerance of the published results; document residual differences attributable to UZF6 and the coupling scheme; lagged vs iterative comparison reported.

**Phase 4 — DISV / unstructured (action 5).**
Grid-agnostic geometry from `flopy.discretization.VertexGrid` (area, top, thickness per `icell2d`); GridIntersect/`Raster.resample_to_grid` input sampling replacing ESRI ASCII; `PlotMapView` plotting replacing the imshow-based `plotLAYER`; La Mata rebuilt on a quadtree DISV grid (GRIDGEN) as the demonstration case, compared against the DIS run.

---

## 7. Summary judgment

The science in MARMITES is coherent and the MF-facing code is more contained than typical for research codes of this vintage — the migration is very feasible. The three real obstacles are not MODFLOW 6 itself but: (1) the `Decimal`-laden, per-cell pure-Python soil loop; (2) error handling that hides defects (§2 items 1–3 are almost certainly biting today, silently); (3) the positional ini-file parsing and `(nrow, ncol)` assumptions welded through the driver. Fix those in Phases 1–2 and the MF6/DISV work (Phases 3–4) becomes a bounded, testable rewrite of one class plus the cbc post-processing.

All structuring design decisions are now taken (§6): DISV for unstructured support, UZF6 always for the recharge/exfiltration interface (unsaturated travel-time delay is required for deep-water-table applications), and per-SP MODFLOW 6 API coupling replacing the Picard loop.
