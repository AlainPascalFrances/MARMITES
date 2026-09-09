# MARMITES / MODFLOW 6 — handoff & next-steps

**THIS FILE (in the repo, `doc/MM-MF6/`) IS THE CANONICAL HANDOFF.** A copy also
exists under the old scratch tree `E:\tmp_claude_marmites\` — it is ABANDONED as
of 2026-09-07; do not read or update it.

Full technical history is in `MARMITES_SFR_LAK_CRR_analysis.md` (read §8.x, esp.
8.15.x for the recharge-coupling bug and 8.14 for the drawdown diagnosis).

**Dev repo:** `E:\00code\MARMITES`, branch `MM-MF6`
(github.com/AlainPascalFrances/MARMITES). Run everything from here.
Dataset (inputs, tracked): `DataSet_LaMata\`
MF6 workspace (outputs): `$MARMITES_WS_ROOT\MF6_ws` — OUTSIDE the repo (see below)
libmf6: `C:\00MODFLOW\mf6.7.0_win64\bin\libmf6.dll`

### CURRENT STATE — READ THIS FIRST (2026-09-09)

**The conversion and the whole figure suite are DONE.** MARMITES runs coupled
to MODFLOW 6 through BMI, is mass-conserving, matches the MODFLOW-NWT
reference, and the native `MARMITESplot_v3` suite is driven entirely off the
MF6 output. `pytest tests -q` = **192 passed, 7 skipped, 0 failed**.

**The next step is §4: run SFR + LAK coupled and validate, then CRR.** Nothing
in §3 remains to do; it is kept as the record of how the plotting was built.

Redraw every figure from a run already on disk, no MODFLOW, ~1 min:
```
python tests\run_lamata_mf6.py --postproc-only --preproc --nlay 2 --mode lagged --run-tag <tag>
```

Two standing rules, both learned the hard way:
1. **The repo holds code, docs, and STRICTLY the files MM and MF read
   directly.** Nothing else — no GIS, no run output, however small or
   convenient. `.gitignore` enforces it; see the split below.
2. **Never `git add -A` while on `master`** (it carries only upstream's
   minimal `.gitignore`).

---

### Repository vs workspace (restructured 2026-09-07)

The repo holds **code, input data and docs only**. All run output goes to a
workspace outside it, set by `--ws-root` or `$MARMITES_WS_ROOT`
(default `E:\00code_ws\LaMata_MM-MF6`):

```
E:\00code\MARMITES\                     REPO
  trunk\  tests\  doc\
  DataSet_LaMata\                       INPUTS ONLY (20 MB) -- strictly the
                                        files MM and MF read; no GIS, no output
    *.txt *.asc  GIS\  MF_ws\ (grids + .ini)  MMsurf_ws\ (MMsurf inputs)

E:\00code_ws\LaMata_MM-MF6\             WORKSPACE (never in git)
  MF6_ws\                               MODFLOW 6 model + output
  MMsurf_ws\                            MMsurf output
  out_<YYYYMMDDHHMM>_<tag>\             MM results:
    _input\                             input parameter maps (IN_*) + general map
    _output\                            all result figures + CSVs
    figures_nwt_comparison\             the 01-07 MF-NWT vs MF6 figures
```

Every map in `_input` and `_output` carries MODFLOW row/column indices on the
top and right and projected coordinates (km) on the bottom and left, drawn by
`MARMITESplot_v3.add_real_coord_axes` -- shared by the native `plotLAYER`
pages (1-based cell centres, `frame='centre1'`) and the plain `imshow` overlay
(0-based indices, `frame='index0'`).

Saved run state (`hi_spinup_*`) is written to the workspace; reading prefers
the workspace copy and falls back to the baseline committed in
`DataSet_LaMata\MF_ws`, so `--strt-heads hi_spinup` works on a fresh clone.
Use `--run-tag` to label a run's results folder.

Branch note: `MM-MF6` was cut at `c849cab`, before master's pyEARTH1D commit
`c766171`. Merging MM-MF6 -> master will raise **7 modify/delete conflicts** on
`trunk/pyEARTH1D/*` (master modified those files; the Phase-1 cleanup deleted
them here). Decide then whether pyEARTH1D stays on master.

### Legacy MODFLOW-NWT reference data (outside the repo)

**Authoritative original — the PhD / paper dataset:**
`E:\00code_ws\LaMata_new_PhD_artigo_2s3L`  (6.1 GB)

| Item | Size | What it is |
|---|---|---|
| `_h5_MM.h5` | 1.3 GB | **NWT reference MM output** — what `tests/plot_water_budget.py::load_reference()` reads for the NWT-vs-MF6 comparison |
| `MF_ws\_h5_MF.h5` | 3.1 GB | NWT reference aquifer output |
| `MF_ws\` | 4.4 GB, 68 files | the full MODFLOW-NWT model (bas6/dis/upw/uzf/wel/drn, `__inputMF_flopy_v3_2s3L.ini`, ASCII grids, .chk) |
| `MMsurf_ws\` | 71 files | MMsurf inputs/outputs (meteo, LAI/veg, irrigation series) — the provenance of the `inputZON_*` series the MF6 pipeline consumes |
| `out_2023*_2s3L_..._PAPER_newTg\` | ~57-61 MB each | the 12 published paper runs |

Both `_h5_*.h5` are dated **2023-01-12** and **cannot be regenerated**: the
MODFLOW-NWT build path was removed from the code base in Phase 1. Treat this
folder as read-only archival provenance — it is also the only source for the
MMsurf side, which is py3-ported but wired into nothing.

Duplicates of the two `.h5` that had accumulated under
`E:\tmp_claude_marmites\MARMITES\DataSet_LaMata\` were verified byte-for-byte
identical to the originals above (size, mtime and full `cmp`) and **deleted on
2026-09-07** to reclaim 1.5 GB. Do not re-copy them; point at the original.

To restore the NWT-vs-MF6 comparison figures, either copy
`LaMata_new_PhD_artigo_2s3L\_h5_MM.h5` into the run's `DataSet_LaMata\`
(untracked, .gitignored) or point `DS` in `plot_water_budget.py` at that folder.
Without it the script prints "reference not loaded" and emits new-run figures only.
---

## 0-ter. COMMAND-LINE FLAGS of `tests/run_lamata_mf6.py`

The single entry point. Everything below is a flag of that script; the
canonical spun-up run is in section 0.

**Where things go**
| flag | default | meaning |
|---|---|---|
| `--ws-root DIR` | `$MARMITES_WS_ROOT`, else `E:\00code_ws\LaMata_MM-MF6` | root of ALL run output, outside the repo |
| `--ws DIR` | `<ws-root>\MF6_ws` | the MODFLOW 6 workspace itself |
| `--run-tag TAG` | `<nlay>lay_<mode>` | names this run's results folder `out_<YYYYMMDDHHMM>_<TAG>` |

**Model set-up**
| flag | default | meaning |
|---|---|---|
| `--libmf6 PATH` | none | path to `libmf6.dll`; without it the run stops after writing the MF6 files |
| `--nlay {2,6}` | 6 | parameter set: 2 reads `_2s1L.ini` directly, 6 reads `_2s3L.ini` |
| `--grid {dis,disv}` | `dis` | structured or vertex grid |
| `--mode {lagged,iterative}` | `lagged` | coupling scheme |
| `--relax F` | 0.6 | relaxation for `iterative` |
| `--nsp N` | all | truncate to N stress periods (for quick tests) |
| `--aggregated` | off | stress periods instead of daily |
| `--aggregate` | off | derive the 2-layer model from the 6-layer one (comparison only) |
| `--seep {uzf,drn}` | `uzf` | seepage face; `drn` is the validated choice |
| `--seep-cond F` | 10000 | DRN-seep conductance (must be free-draining) |
| `--uzf-vks-scale F` | 1.0 | scale UZF vertical K |
| `--sfr` / `--sfr-rhk F` | off / 0.1 | stream network and its bed K |
| `--lak [SHP]` / `--lak-bedleak F` | off / 1e-3 | ponds as lakes and their bed leakance |
| `--no-ats` | ATS on | disable adaptive time stepping |
| `--strt-dem A B` | none | initial heads from a DEM regression |

**Spin-up and saved state** (state is written to the workspace; reads fall back
to the baseline committed in `DataSet_LaMata\MF_ws`)
| flag | default | meaning |
|---|---|---|
| `--spinup N` / `--spinup-tol M` | 1 / 0.05 | repeat the run N times to equilibrate |
| `--strt-heads PREFIX` | none | start from saved heads, e.g. `hi_spinup` |
| `--steady-means PREFIX` | none | drive the steady SP0 with saved mean recharge/ETg |
| `--save-strt [PREFIX]` / `--save-means [PREFIX]` | auto after a spin-up | save that state for reuse |

**Post-processing** (all output lands in `<ws-root>\out_<stamp>_<tag>\`)
| flag | default | meaning |
|---|---|---|
| `--postproc` | off | the full native figure suite + the 01-07 water-budget figures |
| `--preproc` | off | input maps |
| `--gis-ws DIR` | `MARMITES_GIS_WS`, else `E:/00code_ws/LAMATA_new/GIS` (the shapefiles stay in the WORKSPACE -- nothing in MM or MF reads one, so none belongs in the repo) | GIS layers for the general map (`_input/IN_000_general_map.png`); the figure is skipped when the workspace or geopandas is missing |
| `--postproc-only` | off | re-draw from a run already on disk: reads the coupled HDF5 and the MF6 output, builds nothing and runs nothing (implies `--postproc`; combine with `--preproc` for the input maps) |
| `--sankey-min-flux MM` | 0.05 | hide flows below this on the CORE Sankey |
| `--no-sankey-full` | full on | skip the all-flux Sankey |
| `--map-days N` | 6 | head maps on N evenly spaced days (0 = time mean only) |

**Diagnostics**
| flag | meaning |
|---|---|
| `--build-only` | write the MF6 files and stop |
| `--standalone MF6EXE` | build then run `mf6` directly (no API) — isolates model faults from coupling faults |
| `--probe` | list the MF6 memory variables and exit |
| `--max-discrepancy F` / `--allow-bad-budget` | budget guard threshold / do not fail on it |

**Environment variables**
| variable | meaning |
|---|---|
| `MARMITES_WS_ROOT` | default for `--ws-root` |
| `MARMITES_NWT_REF` | the legacy `_h5_MM.h5` used for the MODFLOW-NWT comparison panels |

---

## 0-bis. OPEN ISSUES (recap 2026-09-09)

Ordered by what blocks "post-processing reproduces MARMITESplot_v3 output".

### A. Post-processing vs the native MARMITESplot_v3 suite — COMPLETE ✅
| Native routine | State |
|---|---|
| `plotTIMESERIES_CATCH` | ❌ NOT produced -- its panels came out empty and the per-point series cover the same fluxes (removed 2026-09-09) |
| `plotWBsankey` catchment | ✅ works, core + full, MM-side MB ~0 % |
| `plotWBsankey` per obs point | ✅ **FIXED 2026-09-07** — all 11 points render |
| `plotLAYER` | ✅ **CORRECTED 2026-09-07** — 17 flux maps, legacy conventions, obs points overlaid |
| map layout | ✅ **UNIFIED 2026-09-09** — `MMmap_*` and `GWmap_*` come from ONE function (`_native_result_maps`) with one `plotLAYER` call, so they share a sheet layout. They had drifted: the MM maps passed `nlay=1`, which puts plotLAYER in single-column mode and gave them a full-width panel against the aquifer maps' half-width one |
| map axes | ✅ **2026-09-09** — every map carries MODFLOW row/column indices top and right, projected coordinates (km) bottom and left, from one shared helper `MARMITESplot_v3.add_real_coord_axes` (`frame='centre1'` for the plotLAYER mesh, `'index0'` for a plain imshow) |
| input maps | ✅ **2026-09-09** — 26 native `IN_*` parameter-field pages (geometry, aquifer properties, DRN/GHB, UZF soil parameters, ibound, MM zoning, vegetation areas). The duplicate `aq_*`/`mm_*` imshow set is gone |
| general map | ✅ **NEW 2026-09-09** — `_input/IN_000_general_map.png`, the site map rebuilt from the workspace GIS layers after `GIS/LaMata_MM_MF_202109.png`. The only figure needing geopandas; skipped cleanly without it |
| `obs_heads.png` | ✅ **2026-09-09** — all 11 monitoring points (was 4: the others have no `inputObsHEADS_*` file, and the filter dropped them), blue tones |
| `plotTIMESERIES` (obs soil column) | ✅ **DONE 2026-09-07** — 11/11 obs points |
| `plotTIMESERIES_flxGW` | ✅ **DONE 2026-09-07** — 11/11 obs points |
| `plotCALIBCRIT` (RMSE/RSR/NSE/R) | ✅ **DONE 2026-09-07** — 4 criteria; heads at 4 pts, SM at 2 |
| full `plotLAYER` set + time selection | ✅ **DONE 2026-09-07** — 7 per-layer aquifer maps + head time series |

**Per-point Sankey — diagnosed.** All 11 points raise, from matplotlib:
`ValueError: The connection cannot be made, which may occur if the magnitude
of flow 1 of diagram 3 is less than the specified tolerance`.

Diagram 3 is the `MF UZF` block; its flow 1 is `-Rg_1` (recharge to layer 1),
which the next block (`MF layer 1`) connects to via `prior=3, connect=(1,0)`.
At a SINGLE cell Rg_1 is frequently ~0 over a given hydrological year, and
matplotlib refuses to connect a flow below its `tolerance` (default 1e-6).
The native code only guards LABELS (via `treshold`); it never guards the
CONNECTING flow, because at catchment scale that term is never ~0.

This also explains why `_obs_C1_WBsankey_0whole.png` exists alone: plotWBsankey
saves inside the hydro-year loop, so C1's k=0 page was written before a later
year raised. C1's per-year pages are missing too - it is not a success.

Fix options (a native change, in the spirit of "modify, do not reimplement"):
expose matplotlib's `tolerance` on `plotWBsankey`, and/or floor the connecting
flow to a small non-zero epsilon so the chain can always be built. Must be
validated at several obs cells, not just one.

Stage 2 is now cheap: the `mm_obs` / `mms_obs` capture it needs already exists.

### B. HDF5 vs reading MF6 output — SETTLED 2026-09-08 🟡

Decision: **the aquifer side reads the MF6 output directly; the MM side keeps
its (small) HDF5.** That split is forced, not a preference — the MM soil and
surface fluxes are computed in Python by MARMITES and exist in NO MODFLOW file,
so no amount of cbc reading can recover them.

Done: every aquifer term the figures need (recharge, storage, seepage, drainage,
groundwater ET, vertical exchange, heads) now comes from the cell budget and the
`.hds`, via `_aquifer_pass` / `_aquifer_map_pass` with their cached digests. No
legacy HDF5 is in any figure path.

Still open, and deliberately so: `_coupled_<mode>.h5` is ~182 MB for 1949 SPs,
96 % of it six per-cell x per-SP arrays (heads, perc, etg, exf, rejinf, runoff,
29 MB each). Slimming it would leave ~8 MB. The measured facts, so the decision
can be made on evidence rather than memory:

* the FIGURES take ETg and Ro from the `wb_ts` / `wb_map` aggregates, **not**
  from the `etg` / `runoff` arrays — so those two plots are safe either way;
* `heads`, `rejinf` and `exf` are recoverable from the `.hds` and the cbc
  (`REJ-INF`, and the seepage face now that paknam2 separates it);
* `perc` is read only for `.shape[0]`;
* `load_new()` loads five of them eagerly, so any removal must relax that;
* `--save-means` uses `res['etg']` IN MEMORY during a run, so it is unaffected.

Cost of keeping them: ~174 MB per run, untracked, in the workspace. Cheap to
revisit; nothing depends on doing it.

### C. MODFLOW-NWT comparison — RESTORED 2026-09-07 ✅
`plot_water_budget.load_reference()` reads `DS/_h5_MM.h5`, but `DS` is now the
repo's inputs-only `DataSet_LaMata`, which no longer holds it. The function
returns None, so every comparison figure quietly degrades to "new run only".
Fix: point it at the archive
`E:\00code_ws\LaMata_new_PhD_artigo_2s3L\_h5_MM.h5` (configurable), and cache a
small digest (catchment-mean series + time-mean maps, a few MB) so the 1.3 GB
file is read once rather than every time. The legacy `_h5_MF.h5` (3.1 GB,
aquifer side) is not used at all yet and is the reference for comparing the
MF6 aquifer terms.

### D. Reading the cbc — RESOLVED 2026-09-07 ✅

Benchmarked 2026-09-07 on the 1950-SP run (cbc 862 MB, uzf.cbc 773 MB,
coupled h5 191 MB):

| Operation | Time |
|---|---|
| HDF5 open | 3 ms |
| HDF5 read `wb_ts` (0.4 MB) | **3 ms** |
| HDF5 read `heads` (29 MB) | 19 ms (~1.5 GB/s) |
| CBC open + build index | **1.9-4.4 s** (irreducible floor) |
| CBC `get_data(text=...)` one call, ALL 1950 SPs | 3.6 s for 6 record types |
| `get_data(..., full3D=True)`, single SP | ~2 ms |
| `get_structured_faceflows(grb_file=...)`, single SP | **10.6 ms** |
| FLF for all SPs via a precomputed JA index | **0.02 s** |

**Answer to "can the cbc match the HDF5?" - no, and it does not need to.**
The gap is fundamental: the cbc is 1.6 GB of per-cell, per-package, per-SP data
that must be parsed record by record, against a 0.4 MB pre-aggregated
contiguous slab. Cold read is ~250x slower. But the useful figure is seconds,
not minutes, and it is a one-time cost next to an 11.5-minute model run.

**The 9.3 minutes was self-inflicted, not inherent to the cbc.** Three faults in
`_aquifer_layer_fluxes`, worth ~100x together:
1. `get_data` called once per stress period instead of once per record type.
2. `get_structured_faceflows(grb_file=...)` **re-parses the .grb on every
   call** - 10.6 ms x 1950 = 17.6 s. Precomputing each cell's downward JA
   connection index once (0.01 s) reduces the whole FLF extraction to 0.02 s.
3. The entire extraction repeated per target (catchment + 11 obs = x12) instead
   of reducing for every target inside one pass.

**Design adopted:** one pass over the stress periods, reducing for ALL targets
simultaneously, FLF via the precomputed JA index, and the small result cached as
a digest (a few hundred KB) keyed on the cbc's size+mtime. First post-processing
pays the pass; every re-plot afterwards is HDF5-speed.

**Two flopy defects found (flopy 3.10.0; verified still present on master /
3.11 on 2026-09-07):**
- `get_structured_faceflows(ia=..., ja=...)` **without** `grb_file` raises
  `UnboundLocalError`/`NameError`: the body does `for n in range(grb.nodes)`
  while `grb` is only bound inside `if grb_file is not None`. PR #1968 added the
  `nlay/nrow/ncol` alternative but left that reference, so the documented
  ia/ja path has never worked. No open issue - worth reporting upstream.
  Workaround: pass `grb_file` (slow) or bypass the helper, as we do.
- Sign convention: flopy applies `flows[face][n] = -1 * flowja[i]`. A
  hand-rolled JA extraction must negate to match; pin it with a test.

### Resolved on 2026-09-07

**Per-point Sankey.** matplotlib imposes two OPPOSING conditions on the flow
that connects one block to the previous: it must be ABOVE `tolerance` to get
an arrow angle, and the two sides must cancel to WITHIN `tolerance`. No
tolerance satisfies both when the connecting flow is ~0, which is why lowering
it merely moved the error. `plotWBsankey` now floors the five terms that carry
a connection (Pe, I, Rp, Rg_L, FLF_L) to `eps_connect` (1e-4) -- the same
number on both sides, so cancellation stays exact. All 11 points render, 8
pages each; MM-side closure at P0 is MMsurf 0.1 %, MMsoil 0.0 %, MF-UZF -0.0 %.

**Cell-budget reading.** 558 s -> 16.4 s cold, 0.017 s from the cached digest.
One pass reduces every record for every target; FLF comes from a precomputed
JA index instead of get_structured_faceflows; each record series is read in a
single call with a MemoryError fallback to per-stress-period reads.
Correctness fix found on the way: MF6 writes BOTH drain packages under the
text `DRN`, so records are now addressed by (text, paknam2). The 1954-cell
seepage face (~-2958 m3/d in layer 1) had been missed entirely.

**NWT comparison.** Now reads the archive copy; first quantitative check of
the conversion: forcing (P/Pe/Ei) identical, ETsoil within 0.7 %, and the
differences confined to the surface split (I +45.0, Ro -20.6, Rp +36.4 mm/y),
consistent with UZF6 EPSILON 3.5 vs UZF1 2.0.

### E. Carried over
- **SFR / LAK** build and reload but have never been run coupled or validated.
- **CRR** (Daoud cascade routing + reinfiltration) not started.
- **MARMITESsurf** is py3-ported but wired into nothing; its inputs are now in
  the repo, so re-running it is possible again but unexercised.
- Merging `MM-MF6` -> `master` raises 7 modify/delete conflicts on
  `trunk/pyEARTH1D/*`; `master` carries only upstream's minimal `.gitignore`.
- The water-table drawdown remains a **calibration** matter, out of scope
  (see §2.8).
- **flopy `get_structured_faceflows(ia=, ja=)` is still broken upstream** —
  PR #1968 added the parameters but left `for n in range(grb.nodes)`, so it
  raises. We carry `_ja_down_index` in `marmites_postprocess.py` instead. No
  issue is filed upstream. Remember this if SFR/LAK work needs face flows.
- **Six tests cannot run without an `mf6` binary on PATH.** Their fixture
  (`tests/test_postprocess.py::tiny_run`) was repaired on 2026-09-09 (its UZF
  `packagedata` was a field short of the MF6 schema, so all six errored the
  moment `mf6` was reachable), but the repair itself has never been executed.
- Cosmetic, pre-existing: unused locals/imports in `run_lamata_mf6.py`
  (`NCROP`, `check`, a re-imported `flopy`) and two unused imports in
  `marmites_postprocess.py`.

---

## 0. STATUS SNAPSHOT (what works today)

- Python-3.12 / MODFLOW-6 / UZF6 conversion is COMPLETE and behaves like the
  MODFLOW-NWT reference (the conversion goal). MM-side fluxes match; the coupled
  run is mass-conserving and converges.
- The recharge coupling is FIXED (was the big hidden bug — see §2 below).
- Two-layer model (`--nlay 2`, reads `_2s1L.ini` directly), DRN-seep seepage
  (`--seep drn`, cond 10000), ATS, spin-up with reusable `hi_spinup` heads +
  means, DEM-regression initial heads, and a pre/post figure suite all work.
- SFR network, LAK (EMBEDDEDV ponds), MVR stream-through-ponds all BUILD and
  reload; not yet run/validated coupled (see §4).
- The **native figure suite is complete** (§0-bis A): per-point and
  catchment Sankeys, per-point time series, calibration criteria, the MM flux
  maps and the per-layer aquifer maps, the input parameter maps and the site's
  general map. One `--postproc-only --preproc` emits the lot with no skips.
- Test suite: **192 passed, 7 skipped, 0 failed**. Six of the skips need an
  `mf6` binary on PATH (they exercise `run_postproc`/`run_preproc` against a
  tiny real model); the seventh needs `pyshp`. Worth running locally with
  `E:\00code\moflowapi\emsdatasets\bin` on PATH after touching pre/post.

Canonical run command (recharge now couples correctly):
```
python tests\run_lamata_mf6.py --libmf6 C:\00MODFLOW\mf6.7.0_win64\bin\libmf6.dll ^
  --mode lagged --nlay 2 --seep drn --strt-dem 0.9995 -2.0 --spinup 5 --preproc --postproc
```
Reuse a saved equilibrium (skip spin-up):
```
python tests\run_lamata_mf6.py --libmf6 ... --nlay 2 --seep drn ^
  --strt-heads hi_spinup --steady-means hi_spinup --preproc --postproc
```

---

## 1. VERIFIED 2026-07-25 — the unverified edits are sound ✅

`python -m pytest tests -q` (flopy env) = **192 passed, 7 skipped, 0 failed**
as of 2026-09-09 (190 when this section was written).
The wedged-shell edits all compile and pass:
- `trunk/marmites_coupler.py`: `_bind_first` now accepts a reachable,
  correctly-sized pointer even if not in `get_input_var_names()` (SINF fix);
  binds `SINF` before `FINF`; writes fluxes AFTER `prepare_solve` via a
  `write_cb` threaded through `_advance`/`_one_step`; SFR runoff folded into the
  same callback; `check_solution` has an infiltration-fidelity guard.
- `tests/plot_water_budget.py`: body refactored into `make_figures(ws, mode,
  no_reference, verbose)`; `main()` calls it.
- `tests/run_lamata_mf6.py`: `--postproc` calls `pwb.make_figures(a.ws, mode=)`.

The initial run showed 3 failures, all in `tests/test_coupler_mock.py` — STALE
tests, not a coupler bug. They simulated a "missing" MF6 variable by only
dropping its name from `get_input_var_names()`, while the mock's
`get_value_ptr` still returned a valid array for every leaf. Under the SINF fix
(the var list is NOT authoritative) those "missing" vars correctly bind, so the
mock no longer reproduced absence. Real MF6 signals absence by RAISING in
`get_value_ptr` (that is exactly why SINF — absent from the list, resolvable —
binds and a truly-absent var does not). Fix (mock only, coupler untouched):
`_BadApi.get_value_ptr` now raises for the leaf each mode hides
(`no_finf`→SINF/FINF, `no_gwd`→GWD, `no_q`→Q). All 21 mock tests pass.

Pending confirmation from the user's Spyder run (§6 step 2): one short coupled
run should log `coupler: UZF infiltration bound to LAMATAMM/UZF/SINF` and the
aquifer-balance recharge should VARY (not a flat 977 m3/d).

---

## 2. CRITICAL LESSONS (do not relearn these the hard way)

1. **UZF infiltration = SINF, not FINF.** MF6's operative infiltration array in
   the memory manager is `<MODEL>/UZF/SINF`. `FINF` resolves to a valid but
   non-operative pointer. Binding FINF made UZF apply a constant `perc_user`
   (977 m3/d) for the whole project while looking wired.
2. **Write API inputs AFTER `prepare_solve`.** UZF re-derives SINF from its
   period data during BOTH `prepare_time_step` and `prepare_solve`. A value
   written before `prepare_solve` is overwritten. Proven with `tests/diag_sinf.py`
   (position A=fails, B/C=stick). This applies to every API-written input
   (UZF SINF, WEL Q, SFR INFLOW).
3. **`get_input_var_names()` is not authoritative** — advanced-package vars
   (UZF SINF) are reachable via `get_value_ptr` but absent from that list.
   Validate a bound pointer by SIZE, not list membership.
4. **The infiltration-fidelity guard** in `check_solution` fails the run if the
   written percolation varies but UZF INFILTRATION is constant — keep it; it is
   what caught the bug.
5. **A steady-state SP ignores the initial-head file** (solves dh/dt=0). To seed
   near equilibrium, drive SP0 with mean recharge/ETg (`--steady-means`), not
   just carried heads.
6. **max_outer must equal the IMS OUTER_MAXIMUM** or ATS never sees a failed
   step and cannot retry with a smaller dt.
7. **Files can be stale/mismatched.** Always confirm `.hds`/`.lst`/`.uzf.cbc`
   and `_coupled_*.h5` are from the SAME run (timestamps, nper) before drawing
   conclusions. A 1-year `--nsp 365` test's cbc mixed with a full-run h5 wasted
   a diagnosis.
8. **The remaining water-table drawdown is a CALIBRATION matter, not a bug.**
   Recharge (~67 mm/yr) < discharge (~97): the aquifer bleeds storage and the
   table sinks below the (stable) observed heads. `--uzf-vks-scale` does NOT fix
   it (recharge is not UZF-throttled once SINF is bound). The NWT reference
   drains too. This is La Mata hydrology (recharge/discharge/BCs/soil params),
   the modeller's domain — out of scope for the conversion.

---

## 3. PLOTTING — COMPLETE 2026-09-09 ✅ (kept as the build record)

**Nothing here remains to do.** All three stages shipped and were reviewed
figure by figure with Alain. What `--postproc-only --preproc` now emits, into
`out_<stamp>_<tag>\`:

| folder | content |
|---|---|
| `_input\` | `IN_000_general_map.png` + 26 native `IN_*` parameter-field pages |
| `_output\` | per-point time series (3 pages x 11 points), catchment + per-point Sankeys (whole period, and per hydrological year for the catchment), calibration criteria (NSE/RMSE/RSR/r), 17 `MMmap_*` flux maps, 15 `GWmap_*` aquifer maps incl. a head time selection, `obs_heads.png`, listing/UZF/SFR budget figures, and the head / depth / storage grids as CSV |
| `figures_nwt_comparison\` | the 01-07 MODFLOW-NWT vs MF6 figures |

Decisions that stuck: the native functions were MODIFIED, never reimplemented;
the aquifer side is read straight from the MF6 `.cbc`/`.hds` (see B) with a
per-workspace digest cache; the MM side still comes from the coupled HDF5.

The rest of this section is the original plan and the stage-by-stage record.
Read it for WHY something is the way it is, not for what to do next.

Goal (user): MODIFY `trunk/MARMITESutilities/MARMITESplot/MARMITESplot_v3.py`
to consume MF6 output and produce ALL its native figures — do NOT reimplement/
imitate. Decisions already taken:
- **Feed data via the legacy MARMITES HDF5 layout** (emit the full arrays the
  old driver read), so the native functions run almost unchanged.
- **Deliver in stages, review each** before moving on.

### Data contract (from the old driver `trunk/startMARMITES_v3.py`, post-proc
section ~lines 1550-2740). The native functions consume:
- `flx` / `flxCatch_lst` / `flxObs_lst`: a LIST indexed by flux position, each
  element a per-stress-period time series (catchment-mean, or at an obs cell).
- `flxIndex_lst`: dict mapping exact flux-name -> position. Names MUST match or
  the native functions raise KeyError. Required names include:
  MM soil: `iP iEi iPe iPE iPT iRo iEow iI iSsurf idSsurf idSsoil iperc idSu
  iETg iEg iTg iETsoil`; per soil layer: `iEsoil iTsoil iRsoil iExf_1`; per
  MF layer L: `iEg_L iTg_L iRg_L idSg_L iFRF_L iFFF_L iFLF_L iEXFg_L iWEL_L
  iDRN_L iGHB_L iCH_L`.
- `cMF` attributes used: `inputDate` (matplotlib date numbers), `Mlay`, `Mnlay`,
  `nlay`, `nrow`, `ncol`, `wel_yn`, `drn_yn`, `ghb_yn`, `drncells`, `ghbcells`,
  `ibound`, `perlen`, `hnoflo`.
- `ncell_MM` (active-cell count per layer), `HYindex`/`indexTime` (hydro-year
  boundary indices), `year_lst`, `iniMonthHydroYear`, `DATE`.
- Observations: `inputObsHEADS_{C1..C5,P0,W1}.txt` (date, head), `inputObsSM_*`
  (soil moisture), `inputObsRo_catchment.txt` (runoff); `inputObs.txt` lists
  point name/X/Y/lay (lines starting with # or ## are disabled).

### Native functions to drive (all in MARMITESplot_v3.py):
- `plotTIMESERIES_CATCH(cMF, flx, flxLbl, fn, title, hmax, hmin,
  iniMonthHydroYear, date_ini, date_end, flxIndex_lst, obs_catch=..., ...)`
  — ALREADY WORKS from `wb_ts` (see `native_suite` in
  `trunk/ppMF6/marmites_postprocess.py`, which reassembles the combined MM+soil
  flux array). Returns calibration stats.
- `plotWBsankey(path, DATE, flx, flxIndex, fn, indexTime, year_lst, cMF,
  ncell_MM, obspt, fntitle, ibound4Sankey, ...)` — the Sankey. Needs the FULL
  per-layer flux list above. Splits per-MF-layer terms onto units via
  `cMF.Mlay`. Reads `flx[flxIndex['iX']][i:indexend]` per hydro-year.
- `plotTIMESERIES(cMF, i, j, flx, flxLbl, flxIndex, Sm, Sr, fn, suptitle, title,
  clr, hmax, ...)` — per-obs-cell soil-column time series with obs overlay.
- `plotTIMESERIES_flxGW(...)` — per-obs-cell groundwater fluxes.
- `plotCALIBCRIT(calibcritSM, ..., calibcritHEADS, ..., fn, title, calibcrit,
  ...)` — RMSE/RSR/NSE/R for heads and soil moisture vs obs.
- `plotLAYER(days, str_per, Date, JD, ncol, nrow, nlay, nplot, V, cmap, CBlabel,
  msg, plt_title, MM_ws, mask=..., hnoflo=..., pref_plt_title=..., cMF=...)`
  — per-layer maps; V shape (ndays, nlay, nrow, ncol). ALREADY partly wired
  (native_layer_maps / _native_flux_maps in marmites_postprocess.py). `cmap`
  must be a colormap OBJECT (matplotlib.colormaps['viridis']), not a string.

### Legacy HDF5 to emit from the coupled run (basis for all stages):
`_h5_MF.h5` datasets, each `(nper, nlay, nrow, ncol)`, from the MF6 hds+cbc:
  `heads_d`, `RCH_d` (=UZF-GWRCH), `DRN_d` (drn + drn_seep), `GHB_d`, `STO_d`,
  `WEL_d`, `EXF_d`, and `FLF_d` (layer-to-layer, from the GWF FLOW-JA-FACE or
  UZF/lower-layer exchange). `FRF_d`/`FFF_d` were zeroed by the old driver — OK
  to write zeros.
`_h5_MM.h5` datasets: `MM` per-cell/per-SP `(nper, ncell, nidx)` and `MM_S`
  `(nper, ncell, nsl, nidx_s)`. The coupler currently AGGREGATES these
  (`wb_ts`, `wb_map`); it must additionally stream the full arrays to H5 during
  the run (write per-SP to avoid holding ~0.7-1.5 GB in memory; use gzip).
  Per-layer groundwater ET split (`iEg_L`, `iTg_L`) and recharge-per-layer
  (`iRg_L`) must be derived — check how MMsoil partitions ETg across layers.

### Stage plan (tasks #42-45 in the task list):
- **Stage 1** (#42, #43): CATCHMENT Sankey DONE + VALIDATED 2026-07-25. ✅
  `plotWBsankey` now takes a `treshold=0.05` keyword (was hardcoded 5E-2);
  `marmites_postprocess._native_sankey` + `_aquifer_layer_fluxes` +
  `_hydro_year_index` assemble `flxCatch_lst`/`flxIndex_lst` (MM from
  `wb_ts`/`wb_ts_soil`, aquifer per-layer from the cbc: UZF-GWRCH→Rg,
  STO-SS/SY→dSg, WEL→Eg/Tg split by catchment ratio, DRN→DRN, FLOW-JA-FACE→FLF
  via get_structured_faceflows; FRF/FFF/GHB/CH=0; ETg drawn as Eg/Tg so WEL not
  double-counted). `native_suite` emits a decluttered CORE diagram + optional
  FULL (`treshold=0`); runner flags `--sankey-min-flux`, `--no-sankey-full`.
  Depth conv: `conv_fact/total_surface_area` (same basis as wb_ts).
  VALIDATION: MM-side MB closes (MMsurf 0.0 / MMsoil -1.1 / MFUZF -0.1 %); MF6's
  own layer-2 budget closes to net -0.00 with the extracted terms at a settled
  SP. CAVEAT on the 30-day cold-start test run: absolute values inflated
  (Ro>P) by the SP1 pulse, and MFL2 whole-period closure reads 157% only because
  FLF oscillates sign per-SP and averages to ~0 while recharge is one-signed --
  a spun-up `--seep drn` run's hydro-year SUMS will close. plotTIMESERIES_CATCH
  still worked then (native_wb_catchment.png; dropped 2026-09-09, see A).
  Suite 190 pass / 7 skip / 0 fail.
  PER-POINT Sankeys (flxObs_lst): DONE 2026-07-25. ✅ The coupler now captures
  the FULL per-SP MM flux vectors at the obs cells (`MF6Coupler(obs_idx=,
  obs_names=)` -> `res['mm_obs']` (nper,nobs,nidx), `mms_obs`, `obs_ij`,
  `obs_names`; auto-saved to `_coupled_*.h5`; tiny vs the whole grid).
  `resolve_obs_cells` reuses `cPROCESS.inputObs` to map inputObs.txt points ->
  MM cell-list positions (runner resolves once when `--postproc`, passes to the
  coupler each spin-up cycle). Postproc `_native_sankey_obs` builds a per-point
  `flx` (MM from that cell's captured series, aquifer from that single cell via
  `_aquifer_layer_fluxes(sel_ij=[(i,j)])`) and renders one Sankey per point
  (`_obs_<NAME>_WBsankey_*.png`). VALIDATED offline: obs resolution real, per-
  cell aquifer extraction genuinely cell-specific (P0 Rg1=1502/FLF=39.1 vs
  catchment Rg1=2766/FLF=0.7), renders coherently; MM side needs a real re-run
  to populate `mm_obs` (offline test broadcast wb_ts as a placeholder). Coupler
  mock tests added (obs capture + back-compat). **TO GET REAL per-point values:
  re-run with `--postproc` (obs capture auto-on); the 09:04 h5 predates the
  capture so has no mm_obs.** RELOAD in Spyder: run_lamata_mf6.py,
  marmites_coupler.py, marmites_postprocess.py, MARMITESplot_v3.py.
  NEXT: Stage 2 (per-obs-cell plotTIMESERIES + plotCALIBCRIT) reuses the SAME
  mm_obs/mms_obs capture.
- **Stage 2** (#44): stream full `MM`/`MM_S` to `_h5_MM.h5`; drive
  `plotTIMESERIES` + `plotTIMESERIES_flxGW` at each obs cell; `plotCALIBCRIT`
  from `inputObsHEADS_*`/`inputObsSM_*`.
- **Stage 3** (#45): drive `plotLAYER` over the full flux set + time selection
  from the legacy H5 (heads, recharge, exf, ETg, storage, per layer).

### Build-and-test discipline (the reason we restarted):
Implement each stage against a REAL short run (`--nsp 30`) and view the PNGs
with the Read tool before declaring it done. Do not ship untested index-matching
code. `diag_sinf.py` is a template for isolating BMI issues.

---

## 4. NEXT BIG STAGE: SFR / LAK / CRR

Design + decisions are in `MARMITES_SFR_LAK_CRR_analysis.md` §7-8. Status:

DONE (builds + reloads; NOT yet run coupled/validated):
- **Two-layer aggregation** — superseded: `--nlay 2` reads `_2s1L.ini` directly
  (`marmites_layers.py` kept for `--aggregate` comparison only).
- **DRN-seep** seepage (`trunk/ppMF6/marmites_mf6.py`, `seep='drn'`, cond 10000
  — NB user's original "10" was wrong; a seepage face must be free-draining).
- **SFR** network from `inputPONDw.asc` via priority-flood routing
  (`trunk/ppMF6/marmites_sfr.py`); outlet reaches replace the 6 outlet DRN
  cells; `EVAPORATION=0` (MARMITES does E_ow); runoff delivered as reach INFLOW
  by the coupler (now via the post-prepare_solve callback).
- **LAK** — 12 EMBEDDEDV lakes from `GIS/lm_ponds.shp`
  (`trunk/ppMF6/marmites_lak.py`); ponds are all sub-grid so one lake per host
  cell with a stage-volume table; assigned by centroid.
- **MVR** routes the stream through the 11 on-channel ponds (inlet reach->LAK,
  LAK outlet->downstream reach).
- Flags: `--sfr`, `--lak`, `--seep drn`, `--sfr-rhk`, `--lak-bedleak`.

REMAINING:
- **Run SFR+LAK coupled and validate.** First do it now that recharge couples
  (SINF fix). Command adds `--sfr --lak` to the canonical run. Watch: SFR INFLOW
  binds (advanced-package var — the SINF lesson applies), MVR balance, LAK
  stages sane, no non-convergence/discrepancy from the guard.
- **Where the source geometry lives, and a naming trap.** The network is to
  come from `hydrography.shp` (SFR) and `lm_ponds.shp` (LAK), both in the
  WORKSPACE GIS `E:\00code_ws\LAMATA_new\GIS` — never copy them into the
  repo (see the rules at the top). Alain has flagged a legacy misnomer that
  will bite: **`inputPONDw.asc` and `inputPONDhmax.asc` hold the STREAM
  network, not ponds** (hence `gridSsurfw` = "stream width"). Those two
  rasters, and the `IN_007_gridSsurfhmax` / `IN_008_gridSsurfw` maps built
  from them, lose their job once SFR reads the shapefile.
- **Decide the dependency question first.** Reading a shapefile on the MODEL
  path puts **geopandas** (or pyshp — note pyshp is NOT installed) into the
  build, where today only the one optional general map needs it and the rest
  of the module is deliberately free of geospatial dependencies. Either accept
  it, guard it, or convert the shapefiles to a plain intermediate once.
- **Wanted with it:** an input map of the LAK configuration. `_fig_general_map`
  already draws ponds and hydrography, so it is the natural home — a second
  mode there rather than a new function.
- **CRR** (task #32) — NOT STARTED. Daoud et al. 2022 cascade routing &
  reinfiltration, `CRR_BETA=1.0`. Must run in PYTHON in the coupler (MVR cannot
  reach MARMITES' soil column). MFD weights alpha_ij = beta * S_ij/sum(S_ij).
  Receivers in priority LAK > SFR > downslope MARMITES soil; topographic sinks
  evaporate. REQUIRES a descending-topographic cell ordering, which changes the
  Phase-2 cell-list order contract — plan that carefully. Reference:
  `trunk/SFR_LAK_CRR/cdl_gwf_model_fable_v2.py` (user's CdL model).
- SFR INFLOW / LAK / MVR writes must all go through the post-prepare_solve
  callback (same BMI timing rule as SINF).

---

## 5. FILE MAP (the pieces)

- `trunk/ppMF6/marmites_mf6.py` — clsMF6: builds the MF6 sim (DIS/DISV), UZF6,
  DRN-seep, SFR, LAK, MVR, ATS, initial heads (DEM/array/ini), save/load heads
  asc, layer handling.
- `trunk/marmites_coupler.py` — MF6Coupler: per-SP API coupling (lagged/
  iterative), SINF/Q/INFLOW writes after prepare_solve, ATS sub-stepping,
  check_solution guard, spin-up means, wb_ts/wb_map aggregates.
- `trunk/ppMF6/marmites_sfr.py`, `marmites_lak.py`, `marmites_layers.py`.
- `trunk/ppMF6/marmites_postprocess.py` — all pre/post. `run_preproc` ->
  `_input\`, `run_postproc` -> `_output\`, `native_suite` drives the native
  figures. Key pieces: `_native_result_maps` (MMmap_* AND GWmap_*, one draw),
  `_native_input_maps` (the IN_* set), `_fig_general_map` (the GIS map),
  `_native_sankey` / `_native_sankey_obs`, `_native_obs_timeseries`,
  `_native_calibcrit`, `_aquifer_pass` / `_aquifer_map_pass` (cached cbc
  digests), `_ja_down_index` (our replacement for flopy's broken
  `get_structured_faceflows`).
- `tests/run_lamata_mf6.py` — the driver/CLI (all flags).
- `tests/plot_water_budget.py` — 01-07 figures incl. NWT-vs-MF6, now callable
  via make_figures().
- `tests/diag_sinf.py` — BMI write-position diagnostic (template).
- `trunk/MARMITESutilities/MARMITESplot/MARMITESplot_v3.py` — the native
  suite, fully wired. Modified, not reimplemented. Additions worth knowing:
  `add_real_coord_axes` (the shared projected-coordinate axes) and
  `_nice_tick` / `_series_colour`.
- `trunk/startMARMITES_v3.py` — the OLD NWT driver; the reference for the
  post-proc assembly to port (do not run it; read it).
- `trunk/MM_MF6_conversion/` — checkpoint snapshot (earlier milestone).

---

## 6. SUGGESTED ORDER AFTER RESTART
1. `pytest tests -q` — DONE (§1): **192 pass / 7 skip / 0 fail**. ✅
2. One short coupled run `--nsp 30 --nlay 2 --seep drn` — DONE. ✅
3. Plotting Stages 1, 2, 3 — **DONE 2026-09-09**. ✅ The whole native suite
   runs off the MF6 output; `--postproc-only --preproc` redraws everything in
   ~1 min. Reviewed with Alain and corrected over several passes (colormaps,
   axes, Sankeys, obs points, map layout, the general map).
4. **Run + validate SFR+LAK coupled.**  ← NEXT UP: the surface network moves
   from the `inputPONDw/PONDhmax` rasters to **SFR** driven by
   `GIS/hydrography.shp`, and the ponds to **LAK** driven by
   `GIS/lm_ponds.shp` (both in the workspace GIS, not the repo). Alain has
   flagged the current file names as a legacy misnomer: `inputPONDw.asc` and
   `inputPONDhmax.asc` hold the STREAM network, not ponds. Also wanted then:
   an input map of the LAK configuration.
5. CRR implementation.
(Calibration of the water-table deficit is the modeller's separate task.)
