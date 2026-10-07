# MARMITES → MODFLOW 6, branch `MM-MF6_SFR_LAK_CRR` — instructions for Claude

This file carries the working rules and the project state from Alain's local
sessions (whose memory does not travel) to every session on this branch,
cloud ones included. **Read it whole before acting. Keep it current:** at the
end of each block of work, update §6 (state) and §7 (next steps) and commit
it with the work.

Owner: Alain P. Francés (the user). Last update: 2026-10-07, the §7.1
analysis of run 20261006164429 (analysis notes §8.19), on top of `b483d0e`.

---------------------------------------------------------------------------

## 1. The project

MARMITES (distributed soil-water balance, Francés & Lubczynski 2023) moved from
MODFLOW-NWT to **MODFLOW 6 driven through BMI/XMI** (`modflowapi`, libmf6
6.7.0, flopy 3.10.0). The current phase merges the CdL model's approach into
MARMITES. The plan is `docs/MARMITES_x_CdL_merge_cookbook.docx`, with these work packages:
WP0 config file · WP1 data tiers · WP1b Streamlit UI · WP1c meshes · WP2
three-source ET · WP3 SFR · WP4 LAK · WP5 CRR · WP6 post-processing + obs
exports · **WP7 PEST++-IES** · WP8 hygiene.

Case study: **La Mata** (`example/LaMata/`), daily, 365-day spin-up cycles
on a Voronoi mesh. CdL cannot run yet: it has monthly forcing only, and MM
needs daily forcing.

### Layout
```
code/                    all Python
  tests/run_lamata_mf6.py   THE DRIVER (yes, in tests/): --config --set k=v
                            --run-tag --probe --postproc-only
  marmites_coupler.py       MF6Coupler: MMsoil <-> MF6 per stress period
  marmites_config.py        RunConfig schema + validation of configs/*.toml
  configs/lamata.toml       La Mata run configuration (USER-OWNED, see §2.4)
  mm_paths.py               machine paths: MM_* env var > configs/paths.local.toml
                            (git-ignored) > Windows defaults
  MARMITESsoil/MARMITESsoil_v3.py   MMsoil (flux(), _cell_step(), step())
  ppMF6/marmites_mf6.py     clsMF6: builds the MF6 simulation (UZF, SFR, LAK,
                            MVR, DRN, GHB, WEL, EVT, ATS, IMS)
  ppMF6/marmites_crr.py     WP5 cascade routing (CascadeNetwork, Daoud 2022 Eq. 23)
  ppMF6/marmites_evt.py     Eg/Tg as MF6 EVT curves (et.gw_route = 'evt')
  ppMF6/marmites_postprocess.py   figures, Sankey, lake/obs series, obs exports
  ppMF6/marmites_topology.py      shared-face adjacency on DISV meshes
  marmites_meshes.py / marmites_mesh.py   mesh producers / projection
  marmites_indices.py       MM flux-vector indices (iRunon 31, iEcrr 32, iReinf 33 ...)
  app/                      Streamlit front-end: streamlit run code/app/Home.py
                            pages 1 Grid .. 7 Run (Validation + Run tabs) .. 8 Results
                            app/lib/schema.py = panel help texts
  tools/gis_to_dataset.py   converter shapefiles -> dataset tables (run by hand / Launch)
docs/                    cookbook .docx, docs/MM-MF6/*.md reports
                         (MARMITES_SFR_LAK_CRR_analysis.md = running analysis notes;
                          MARMITES_NEXT_STEPS.md is STALE since 2026-09-09)
example/LaMata/          Tier-A inputs MM and MF read directly — nothing else
```
**Repo = code + input data that MM/MF read DIRECTLY + docs.** Never run output,
never shapefiles. The user asked to be reminded whenever that gets mixed up.

---------------------------------------------------------------------------

## 2. Working rules (the user's, binding)

1. **Never launch La Mata model runs** (spin-up, transient, 1-year). Alain runs
   them himself from the front-end (Run panel). Implement, verify by tests and
   by **build-only / initialize-only** checks in a scratch workspace
   (`--set run.build_only=true --set paths.ws=<scratch>`), then hand off with
   what to set on which panel and what to expect. Never write into the
   user's run workspace. Move any `out_*buildcheck` folder out of WS_ROOT.
   Small toy MF6 models for probing the API are fine.
2. **The front-end owns the model configuration.** Every model choice is a
   config field that a panel asks, plus the code that honours it. Both land
   together, with help text in `app/lib/schema.py`, validation in
   `marmites_config.py`, and the key read by the run. `tests/test_wiring.py`
   (NOT_WIRED list) enforces this. Never hard-wire a choice in the driver.
   A control that changes nothing is worse than none.
3. **Total ET can never exceed PET** (check_solution fails a run on any
   cell-period with ET > PE + PT).
4. **Git.** Stage your own files explicitly. Never `git add -A`.
   **Never commit** `code/configs/lamata.toml`, `docs/MM-MF6/*.txt` or
   `example/LaMata/*`: they are the user's working copies. End every commit
   message with the attribution line the harness gives (2026-10-06:
   `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`). The user
   pushes; do not push unless asked. In a cloud session, commits live only
   in the container until pushed, so ask before ending a block of work.
   **Cloud sessions (user decision 2026-10-06):** Claude pushes to the
   SESSION branch only (never to `MM-MF6_SFR_LAK_CRR`); the user then
   fast-forwards `MM-MF6_SFR_LAK_CRR` from it in his own terminal.
5. When code that the Streamlit app imports changes (schema, config, panels,
   meshes), **remind the user to restart Streamlit**. A model run starts a fresh
   process and picks up model code without a restart.
6. When fixes start to turn into whack-a-mole, step back. Name the structural
   cause and propose a **literature-grounded redesign** (Alain co-authored
   Daoud et al. 2022 and El-Zehairy et al.; he reaches for those).
7. Analyse runs from their logs and outputs. When an earlier statement turns
   out wrong, say so explicitly. Report test failures with their output.
8. While the user's model is running on the same machine, avoid heavy CPU:
   run targeted tests, not the full suite.
9. Code style: match the surrounding code. Explanatory docstrings, dated
   findings in comments where a trap was hit, small functions, no new
   dependencies without asking.

---------------------------------------------------------------------------

## 3. Environments

### Alain's Windows machine (where runs happen)
- Python env `C:\miniconda3\envs\flopy` (py3.12, flopy 3.10.0, streamlit,
  geopandas, rasterio, pyshp, scipy ...). **Always run it activated** (by
  prefix). A bare `python.exe` crashes with 0xC06D007F because MKL's
  delay-loaded DLLs are off PATH. A second, EMPTY env also named `flopy`
  exists under the user profile.
- Launchers + every helper/scratch script: `E:\tmp_claude\helpers`
  (`pt.bat` = pytest in the activated env, `run_py.bat <script>`). Never use
  AppData\Local\Temp, which gets wiped.
- Worktree `E:\00code\MM-MF6_SFR_LAK_CRR`; `E:\00code\MARMITES` = branch MM-MF6.
- Run workspace `E:\00code_ws\LaMata_MM-MF6\MF6_ws_voronoi` (runs, saved
  states `hi_voronoi_lamata_*` = `_l1/_l2/_etg/_perc.asc`, `_state.npz`,
  `.scope.json`); results in `E:\00code_ws\LaMata_MM-MF6\out_<stamp>_<tag>`.
- GIS sources `E:\00code_ws\LAMATA_new\GIS`. MF6 6.7.0 + its Fortran source
  at `C:\00MODFLOW\mf6.7.0_win64` (read the source before trusting a budget line).
- Legacy NWT reference `E:\00code_ws\LaMata_new_PhD_artigo_2s3L\_h5_MM.h5`
  (cannot be regenerated).

### The Windows server (runs from 2026-10-06 on)
- **THE local repo is `X:\3p1p1\MARMITES`** (user, 2026-10-07): update the
  code, docs and everything else there, and commit there. It sits on a
  network share, and git refuses it as "dubious ownership": pass
  `-c safe.directory=*` on each command rather than changing the global
  config.
- `C:\00code\MM-MF6_SFR_LAK_CRR` is another checkout. Never edit it. It is
  read-only for a non-elevated session anyway (BUILTIN\Users: RX). The
  2026-10-06 runs were launched from it (`status.json` cwd), at `84d487f`,
  behind `X:`. Check which checkout a run used before attributing a result
  to a commit.
- Python env `C:\Users\su-alain.frances\AppData\Local\miniconda3\envs\mf6models`
  (the runs use it). It has NO pytest: pytest sits in `X:\tmp_claude\pylibs`
  (pip --target), so the env is untouched. Helpers in `X:\tmp_claude\helpers`:
  `run_py.bat <script>` runs a script in the activated env, `pt.bat <args>`
  runs pytest from the X: repo root, `buildcheck.bat` builds La Mata from
  the X: code into `X:\tmp_claude\scratch_ws` (the user's inputs, copied
  mesh cache and saved state) and initializes MF6 (`--probe`), never a run.
- Machine paths: `C:\00code\...\code\configs\paths.local.toml` (the X:
  checkout has none, so its defaults are E:\). MF6 6.7.0 + source:
  `C:\sw\MODFLOWandCo\mf6.7.0_win64`.
- Test baseline on this server (2026-10-07), `pt.bat` (which gives the
  server's read-only data paths, as the X: checkout has no paths.local.toml):
  all pass but **4 environmental failures** in test_app_pages -- Home.py and
  the two panel-0 tests (Streamlit refuses a page on a network path: run
  them from C:), and the pond-footprint map (no mesh cache, because WS_ROOT
  is deliberately not given to tests). ~8 min.
- `X:\3p1p1\MF6models\LaMata\` holds `DATA_ROOT`, `NWT_REF` and `WS_ROOT`.
  Run workspace `WS_ROOT\MF6_ws_voronoi`, results `WS_ROOT\out_<stamp>_<tag>`,
  per-run `run.log` + `status.json` in `WS_ROOT\runs\<run_id>\`.
- Copied from E:\ (2026-10-06): the runs\ logs, the saved states and
  `out_202610052201_2lay_drn_spinup` (CSVs + its `mfsim.lst`). Not copied:
  the other `out_*` folders, and any run's binary outputs (each run
  overwrites `lamata.hds/.cbc` and `_coupled_lagged.h5`). No GIS workspace,
  so the general map is skipped.

### A cloud session
- None of the E:\ / C:\ resources exist: no run outputs, saved states, GIS,
  NWT reference or mesh cache. Tests that need them skip or fail for that
  reason, not because of a regression. To analyse a run, ask the user to
  share the log and the relevant output files.
- Setup (Linux): `pip install numpy scipy pandas matplotlib h5py flopy==3.10.0
  modflowapi pytest` (plus `shapely pyproj geopandas rasterio pyshp
  streamlit triangle` for the mesh/app/tool tests). MF6 binaries: `python -m
  flopy.utils.get_modflow --repo modflow6 --release-id 6.7.0 <dir>`.
- **Portability gap:** `mm_paths.py` builds Windows names:
  `LIBMF6 = <MM_MODFLOW_DIR>/mf6.7.0_win64/bin/libmf6.dll`, also `mf6.exe`
  and `win64/triangle.exe`. Workaround: create that tree with symlinks
  (`libmf6.dll -> libmf6.so`, `mf6.exe -> mf6`) and set `MM_MODFLOW_DIR`.
  Making mm_paths platform-aware is a reasonable small fix; ask first.
- `code/configs/lamata.toml` in the repo is the user's current run
  configuration: he committed it with the converter output in
  `example/LaMata/` (`a4f7b27`, 2026-10-06). He commits these working copies
  himself at handoff points, so between handoffs his local copies may be
  ahead of the repo. Main settings (2026-10-06), as changes from the
  previous committed copy: nsp 365, run.ats_dtmin 0.0002; mesh cell_near_stream 20, stream_buffer 60,
  trans_levels [20, 60], cell_pond 40; soil.irr_infiltration "column";
  drn line "lm_outlet_drn.shp", cond_per length, value 0.00084; seep cond 100,
  ddrn 0.125, base 0, cond_from value; et.gw_route "evt", evt_nseg 8,
  evt_ramp 0.1; crr.enable true; sfr.enable true, depth 0.5, rbth 0.2
  (rhk 0.1); lak.enable true; spinup strt_heads "hi_voronoi_lamata_rbth02",
  save_strt/save_means "hi_voronoi_lamata_evt"; postproc.input_maps true;
  [solver] complex, outer_dvclose 0.025, outer_maximum 100, inner_dvclose
  0.001, inner_rclose 0.01, cell_averaging "amt-hmk".
- Raw meteo spreadsheets (`_meteoSARDON_TB_200709_201011.xlsx`,
  `_200709_201309.xlsx`) are source data that no code reads (MMsurf reads
  `__meteoTB.txt`). At the user's request they live in
  `E:\00code_ws\LAMATA_new\GIS` and are no longer tracked (2026-10-06), and
  .gitignore no longer lets `.xlsx` into `MMsurf_ws`.

Test baseline on Windows: **1251 passed / 9 skipped** (full suite,
2026-10-06, at 9ab2991, ~7 min). In the cloud, expect
the skips and failures of §3 for missing local data or binaries. Compare
against a run of the same tests at the parent commit before calling
something a regression.

---------------------------------------------------------------------------

## 4. Model design decisions in force (do not undo silently)

- Coupling: **lagged daily, the only mode** (user, 2026-10-07: `run.mode`
  and `run.relax` retired). MMsoil runs first (soil column on top of MF6),
  then MF6 does one stress period through **`do_time_step`**, so MF6's own ATS
  retry runs. The iterative mode drove the solve by hand (a failed step was
  accepted, never retried), and the head feedback it existed for is EVT's
  now. Fluxes are written to arrays MF6 re-applies on every retry:
  UZF `SINF_PVAR` (recharge), `PET_PVAR` (UZF demand). EVT/SFR/LAK period
  values are NOT reloaded on a retry; only PERIOD blocks reset them. A split
  period's rates are booked as **sub-step means** (SUBSTEP_RATES). Results
  file: `MF6Coupler.RESULTS_H5` = `_coupled_lagged.h5` (name kept).
- **NO LEGACY PARAMETER FILE** (user, 2026-10-07: "nothing can come from
  the legacy ini files" -- a new catchment has none). A run builds its model
  description with `clsMF.from_config` (scalars from the TOML) and
  marmites_props (arrays): grid = the dataset rasters' grid; nlay, hnoflo,
  model name from the panel; land surface = `[grid] dem` wrapped onto the
  dataset grid (`land_surface`, REQUIRED); every layer property REQUIRED
  (`apply_layer_properties(required=True)`); cold-start heads =
  `layers.strt` (blank: elevation x a + b, `spinup.strt_dem`); cold steady
  recharge = `spinup.steady_recharge`; DRN/GHB counts from the cells built.
  The dataset is `mm_paths.dataset_dir(paths.case)` (it was the driver's own
  `<repo>/example/LaMata`). The parsing constructor `clsMF(..., MF_ini_fn)`
  survives only for the legacy NWT scripts and some test fixtures.
- MF6 top = **soil base** (top = land − soil thickness). The land surface is
  `clsMF6._land_surface()` (cMF.elev). UZF extdp = roots below the soil.
- **Groundwater ET = EVT, the only coupled path** (user, 2026-10-07:
  `et.gw_route` retired, the WEL route and its ETg wells removed). Two EVT
  packages `evt_eg`/`evt_tg`, one curve per cell per day, evaluated by MF6 at
  the solved head (marmites_evt); the coupler sets `ctx.gw_evt` so MMsoil
  hands over the potential. Each EVT curve starts at the start-of-day head
  with MM's rate and is flat above, so EVT can only fall within the day. The
  UZF demand is capped at PETuzf − EVT max (`uzf_demand(petuzf, etg_cap)`),
  which keeps ET ≤ PET exact. EVT arrays (SURFACE/RATE/DEPTH/PXDP/PETM)
  persist across do_time_step retries. A steady first SP takes the mean ETg
  as a FLAT EVT curve down to the surface cell's bottom (`flat_curve`, all on
  EVT_EG). MMsoil's own Eg/Tg with drawdown remain for the uncoupled path.
  **WEL is left for real pumping** (not built yet; post-processing still
  reads a WEL term for runs made before 2026-10-07).
- Seepage face = DRN_SEEP at the soil base (not UZF GWSEEP: deprecated in MF6
  6.5+, and validate() refuses seep.kind='uzf' under coupling). **No seepage
  drain in stream cells or in any pond-footprint cell.**
- SFR: stream total depth = soil depth + channel depth + rbth (bed top =
  aquifer top − depth); downstream-monotonic beds. A panel VALUE for
  manning/rhk/rbth wins over the converted table (`sfr_param_fixed`).
- LAK (La Mata): EMBEDDEDV, footprint = cell centres inside the pond, host =
  centroid cell, rim = LAND surface, bed = aquifer top − pond depth, stream
  excised from footprints with MVR in/out, MVR factor 1.0; stage carried
  across cycles and in `_state.npz`.
- CRR (WP5): receivers LAK > SFR > soil, α_ij = β·S_ij/ΣS_ij over face
  neighbours, 1−β evaporates each hop, sinks evaporate|route, run-on
  infiltrates top-down (Eq. 1b). Irrigation water infiltrates the whole column
  (`soil.irr_infiltration = column`) in irrigated cells (inputIRRzones.asc).
- Spin-up: periodic cycles carrying heads, UZF water content, soil state and
  previous-day exchanges; each passed cycle saved as `<save>_lastcycle`;
  converged when mean |dWT| < tol.
- Post-processing reads binary outputs **per stress period**
  (`period_steps`/`period_end_kk`: state = last record, rate = time-weighted
  mean). Since ATS there are several records per period.
- Obs exports (WP7's contract): `obs_heads/obs_sm/obs_sfr/obs_et.csv`, long
  format, one row per stress period, stable `obsnme`, written on every run to
  `_output/` and `<ws>/obs_exports/`.

---------------------------------------------------------------------------

## 5. Traps already hit (each cost hours)

- MF6 UZF **routes FINF but reports SINF**, and uzf_ad resets both from
  SINF_PVAR. PET is reset from PETMAX on every iteration, so write the PVAR
  arrays, after prepare_solve.
- In the UZF cbc, `node2` is the UZF object id, not the GWF cell.
- A DISV model writes `<name>.disv.grb` (not `.dis.grb`). On a mesh,
  `delc×delr` is a placeholder (1 m²); use `model_cell_area`. flopy `get_ts` on DISV
  .hds wants `(lay, 0, icell2d)`. Mesh projection uses the (ncpl, 1) convention:
  nrow := ncpl, ncol := 1.
- `np.ma.masked_values(x, v)` overwrites data under an existing mask; use
  `_mask_sentinel`.
- flopy VoronoiGrid duplicates vertex ids at one coordinate: match vertices by
  coordinate.
- flopy writes DISV at 8 digits by default; use float_precision 16, else MF6 and MM
  cell areas differ by up to 2%.
- Dry-cell test in MMsoil uses `botm_l0` = bottom of the **outcrop layer**
  (fixed 6b9defa). It was layer 1's bottom everywhere, which was wrong where layer 1
  pinches out.
- `_pad` in marmites_evt: a Tg curve can have no inner breakpoint (a root
  tip within the ramp of the head); that is a straight line, not an error.
- `budget_uzf/sfr/lak.csv` (`package_budget`, max_samples=120) average an
  even subsample of the ATS records. That is off by up to 6 % (UZF GWF), with
  storage terms that should be zero, and two runs sample different days.
  Never compare runs on these files. Use `budget_terms.csv` (from the listing),
  the run log, or `package_budget(..., max_samples=None)` (2026-10-07, notes §8.19).
- **Machine paths have ONE definition, `mm_paths`** (MM_* env > configs/
  paths.local.toml > defaults, which are the OLD machine's E:\ paths). A
  literal default anywhere else is a bug: marmites_postprocess had its own
  GIS folder, E:/00code_ws/LAMATA_new/GIS, so on the server the general map
  was skipped while panel 0 pointed at the GIS (fixed 2026-10-07). On the
  X: checkout (no paths.local.toml) mm_paths falls back to E:\.
- The Streamlit app keeps old modules in memory: after code changes, restart it,
  or the grid/schema is stale.
- AppTest page tests time out while the user's model runs; re-run them alone.
- `test_aquifer_balance::test_a_real_list_file_reads` fails while a run is
  writing lamata.lst (environmental).

---------------------------------------------------------------------------

## 6. State (2026-10-07)

Done: WP0, WP1, WP1b, WP1c (meshes), WP2 (three-source ET), WP3 (SFR), WP4
(LAK), WP5 (CRR), WP6.3 (soil moisture at depth figure), WP6.4 (obs
exports), EVT route; 2026-10-07: the WEL route for Eg/Tg and the iterative
coupling removed (EVT and lagged are the only paths; §4), and **no legacy
parameter file read by a run** (§4). The run of 2026-10-07 (WEL/iterative
removal) was launched by the user before the parameter-file change. The
parameter-file change is verified by tests and by a build + MF6 initialize
from a copy of the La Mata dataset with every .ini deleted: the MF6 input it
writes equals the parameter-file build to 5e-13 (initial heads differ only
in inactive cells, UZF extdp by 3e-7 m).

Recent runs (La Mata, 4566-cell Voronoi mesh: 20 m stream corridor ratio 2,
40 m pond cells; 2 layers; 1-year spin-up cycles):
- LAK run 20261004183453: converged, MB 0.05%. C1–C3 are pond piezometers
  (hcorr RMSE 0.35/0.48/0.40 m); ponds full all year.
- rbth 0.2 run 20261005081712 (WEL route): no material change vs 0.5; the user
  keeps 0.2. Dropped by the user: the irrigated-field groundwater mound, and the
  gauge record check.
- EVT run 20261005220110 vs that WEL run: in layer-1 cells (96% of the area)
  ETg is identical (23.9 vs 24.0 mm/yr), so the EVT route is validated. The
  catchment gap (24.4 vs 32.0) was entirely the botm_l0 bug in 181
  layer-2-outcrop cells near streams. The WEL route pumped ~175 mm/yr per cell
  from a water table 3.4 m deep; the EVT route reserved that rate (withholding
  it from UZF) but took little. Fixed in 6b9defa. EVT run: 530/615 extra
  sub-steps vs 986/989; outlet 91 vs 84 mm/yr, gauge bias −28% vs −34%.
  Under EVT, hcorr = start-of-day head (MM no longer draws down): calibrate on
  raw MF6 heads. Obs exports had the expected row counts on the real run.
- EVT-fix run 20261006164429 (with 6b9defa; notes §8.19):
  - The fix works: MMsoil now reads the real layer-2 head in the 181 cells.
  - Catchment unchanged against the EVT run (ETg 24.4, every GWF term within
    0.05 mm/yr, heads within 2 cm). That is expected: EVT already
    evaluated the curve at the real head, and the bug only inflated the
    reserve withheld from UZF, which is limited by water there.
  - L2 cells in soil zones 1-2 behave like L1 cells at the same depth
    (L2/L1 0.85-0.93 at 2-5 m). Zone 3 (132 of the 181, Shah 'sand') has
    no Eg and 3.9 mm/yr of Tg (thin soil, fewer trees).
  - The WEL gap is the bug in those cells: L2 contribution 8.1 vs 0.4 mm/yr,
    L1 23.9 vs 24.0.
  - Convergence unchanged: the same 18 SPs fail, and the failures sit at
    SFR cells on storm days, not at the ETg cells.
- Gap to the NWT reference is still large (EXFg ~13 vs 67, Ro ~31 vs 83 mm/yr).
  Sy 0.01 / K 0.05 is a calibration matter (WP7).

---------------------------------------------------------------------------

## 7. Next steps (in order)

1. **User checks the runs of 2026-10-07** (WEL/iterative removal; then the
   parameter-file removal). Expect the EVT-fix run's results to the digit
   but the run tag (`<nlay>lay` by default), no WEL in the budget, the three
   retired keys reported until the config is saved from a panel, and a log
   that names no parameter file. On the panel: `layers.strt` (initial heads
   of a cold start; La Mata's old file gave MF_ws/hi_topL1.asc for both
   layers) and `spinup.steady_recharge` are new. Restart Streamlit first.
2. **Test fixtures off the parameter file** (small): tests that build
   `clsMF` from `__inputMF_flopy_v3_2s1L.ini` as a fixture (test_mf6_build,
   test_lak_ponds, ...) should use `clsMF.from_config`; then the parser and
   the file can go with the legacy NWT scripts.
3. **WEL for real boreholes/extraction** (when the user wants it): a `[wel]`
   section on a panel, plain MF6 input with the coupler hands-off; WEL =
   pumping in the balance, maps and Sankey -- relabel the post-processing's
   WEL term ('ET (groundwater)' today, for old WEL-route runs).
4. **Exact package budgets in post-processing** (small, before 6.6).
   `package_budget` should read every record (time-weighted, about 40 s per
   cbc on La Mata), or the per-period means. `layer_storage_change` should
   be time-weighted too. See the §5 trap.
5. Validate WP6.3 figures (`sm_depth_<pt>.png`) on real output.
6. WP6 remainder: 6.2 (pond volume panels, MVR accounting, water-balance
   graphs), 6.5 (calibcrit groups for streamflow and ET), 6.6 (Results page:
   run picker, run-to-run comparison).
7. **WP7 PEST++-IES** on the obs exports. Lessons from the CdL calibration:
   draw the prior ensemble from the geostatistical structure (pyEMU
   `pf.draw` → `prior_pe.jcb`, `ies_parameter_ensemble`). A diagonal
   bounds-only prior gave spatially white pilot points, checkerboard K, 58% of
   parameters on bounds, ensemble collapse and a diverging forward model. Use
   `ies_autoadaloc`, ~150 realisations, and check posterior Moran's I and
   bound-hitting. The forward run must complete (physical-plausibility gate).
   Runs at that scale belong on a server.
8. WP8 hygiene.
