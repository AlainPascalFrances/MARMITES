# MARMITES — Phase 3 kickoff (handoff to the next model)

**Written:** 2026-07-21, at the end of Phase 2.
**Goal of Phase 3:** replace the removed MODFLOW-NWT backend with a **MODFLOW 6**
model built via `flopy.mf6`, and couple it to MARMITES through the **MODFLOW 6
API** (`modflowapi`/`xmipy`) instead of the old file-based Picard loop.

Read these first, in order: `MARMITES_code_review.md` (the contract — decisions
are in §4.1/§4.2/§4.3/§5/§6), then `MARMITES_Phase1_report.md`,
`MARMITES_Phase2_report.md`, `MARMITES_LaMata_validation.md`. Then this file.

---

## Locked decisions (do not relitigate)

1. **UZF6 always** — no RCH+DRN variant. The unsaturated travel-time delay must
   be preserved (deep-water-table applications). §4.2.
2. **API coupling, no Picard loop.** Two modes sharing one interface:
   `lagged` (MM at SP *n* uses heads from SP *n−1*) first, then `iterative`
   (MM re-evaluated inside MF6's outer-iteration loop, with flux under-relaxation)
   as the default once validated. §4.3.
3. **DISV**, not DISU, for the unstructured target (that's Phase 4; Phase 3 is DIS).
4. **Daily stress periods** as the default (the surface mass balance only closes
   for daily SPs — see validation). Keep `perlenmax` aggregation optional.
5. **MODFLOW-NWT is deleted, not ported.** Don't reintroduce it.

## The `step()` interface — this is the coupling seam

`MARMITESsoil_v3.clsMMsoil` already exposes the grid-agnostic per-SP interface the
API driver must call. Do **not** change its physics; build the MF6 driver around it.

```
cells = mm.build_cell_list(cMF)          # [(cid, i, j, node), ...] active cells; node=i*ncol+j (=icell2d for DISV)
ctx   = mm.build_context(cMF, cells, ...) # static inputs bundled once
state = mm.init_state(ctx)                # carry-over: Ssoil_ini (ncell,nslmax), Ssurf_ini (ncell,)
out   = mm.step(ctx, n, tstart_MF, heads_cell, exf_cell, state)   # advances SP n, mutates state
# out = {'MM':(ncell,24), 'MM_S':(ncell,nslmax,9), 'perc':(ncell,), 'etg':(ncell,)}
```

- **`heads_cell` (ncell,)** — representative groundwater head **[m]** per active
  cell for this SP, from MODFLOW. In `lagged` mode: end of SP *n−1*. In
  `iterative` mode: current outer-iteration head iterate of SP *n*.
- **`exf_cell` (ncell,)** — exfiltration **into** the soil column **[mm/d]**,
  **sign: positive = seepage up into the soil** (the file driver passes
  `-exf_MF_ini`, i.e. it negates the MODFLOW surface-leakage sign; preserve this).
- **`perc` (ncell,)** — gross recharge to the water table **[m/d]**
  (`Rp_bottom[mm/d] / conv_fact`). Feed to **UZF6 `finf`**.
- **`etg` (ncell,)** — groundwater ET **[m/d]** (`ETg[mm/d] / conv_fact`). Feed to
  **WEL** as `Q = -etg_cell * cell_area` [m³/d] (as the old WEL path did).

`conv_fact = {1:304.8, 2:1000.0, 3:10.0}[cMF.lenuni]` (La Mata: lenuni=2 → 1000).

The file-based `runMMsoil` is a thin loop over `step()` and stays as a reference /
non-API fallback. The API driver replaces the loop body: read heads/exf from MF6
memory instead of the HDF5, scatter `perc`/`etg` to MF6 memory instead of scattering
to the structured HDF5.

## Phase-3 task order (from §6)

- **3a — MF6 model construction (`clsMF6`, DIS).** `flopy.mf6`: TDIS (daily),
  DIS, IC, NPF, STO (Newton via `newtonoptions`), IMS, WEL (`auto_flow_reduce`),
  DRN/GHB with grid-agnostic `cellid`, UZF6. See the NWT→MF6 package mapping table
  in **§4.1 of the code review** — it is the build spec. Validate the generated
  files with one conventional full MF6 run before wiring the API.
- **3b — API driver, `lagged` mode.** `modflowapi`/`xmipy`: init, advance the
  steady initial SP, then per SP: `mm.step(...)` with heads/`UZF-GWD` from SP *n−1*
  → write `FINF` and WEL `Q` via `get_value_ptr` → advance MF6. Persist MM outputs
  to HDF5 (the export/plot code still expects the structured datasets).
- **3c — API driver, `iterative` mode.** Re-evaluate `step()` at each outer
  iteration (`modflowapi` `iteration_start` callback / repeated `solve()`), with
  under-relaxation (damping ~0.5–0.7 on ΔFINF/ΔQ between iterations) to prevent
  wet/dry oscillation at exfiltration cells. ~50 lines on top of 3b.

## The subtle parts (where reasoning depth matters)

- **Dry-cell semantics.** `_cell_step` currently detects dry cells by
  `abs(h - cMF.hdry) < 1e-5` and substitutes `botm_l0*1000`. **MF6 has no `hdry`
  sentinel.** Replace with: inactive driven by `idomain`, dryness by `h < botm[cell]`.
  Search `MARMITESsoil_v3.py` for the `NOTE MF6 (Phase 3)` marker — that's the exact
  spot. `hnoflo` masking (`hnoflo`/±0.09) must likewise be driven by `idomain`, not
  magic values (La Mata `hnoflo = 9999.999`, positive — don't assume negative).
- **UZF6 budget names.** `SURFACE LEAKAGE`→`UZF-GWD` (discharge to land surface =
  your exfiltration), `UZF RECHARGE`→`UZF-GWRCH`, GW-ET→`UZF-GWET`. Face flows:
  `FLOW RIGHT/FRONT/LOWER FACE`→ single `FLOW-JA-FACE`. Storage `STORAGE`→
  `STO-SS`+`STO-SY`. The export/`procMF` code and the `startMARMITES_v3.py` cbc-name
  index block need updating for these.
- **UZF6 setup.** Per-cell `packagedata` with `landflag`/`ivertcon` (vertical
  stacking for multi-layer columns) aligned to `cMF.outcropL`; OBS6 replaces the
  `nuzgag`/`row_col_iftunit_iuzopt` gage machinery.
- **Cell geometry.** `_cell_step` still uses `delr[j]*delc[i]` for area
  (`Ssurf_max`, `Eosurf_max`, exf conversion) — marked `TODO Phase 4`. Fine for DIS;
  DISV supplies per-node area from `flopy.discretization.VertexGrid`.
- **Config.** Extend `marmites_config.py` with an MF6 section (the NWT ini schema
  was intentionally not modelled). Normalize path separators (La Mata's
  `SOILparam_fn` uses a Windows backslash).

## Validation — you have a real regression baseline

`tests/validate_lamata.py` stands up the full La Mata model and runs the refactored
MMsoil against real MODFLOW heads. Use it as the template for the API-run comparison.

- The refactor is already proven faithful: identical inputs → identical outputs,
  MB closes to 1e-14. See `MARMITES_LaMata_validation.md`.
- **Important caveat when comparing to `_h5_MM.h5`:** the reference output is **one
  Picard iteration out of sync** with `_h5_MF.h5` (~0.09 m in heads). So a
  whole-field match is not expected against the reference; validate instead that
  (a) MB closes, (b) the water-balance components are physically sensible, and
  (c) `lagged` vs `iterative` converge as SPs shorten. La Mata is a plausibility
  benchmark, not a bit-match, once the coupling scheme changes (this was the
  agreed acceptance in §6 Phase 3).

## Environment traps I hit (save yourself the debugging)

- **Sandbox is Python 3.10; target is 3.12.** `tomllib` is absent on 3.10 —
  `marmites_config` falls back to `tomli` (already installed). Code targets 3.12.
- **`clsUTILITIES.readFile` does NOT strip tokens** — every caller must `.strip()`.
  (This bit the validation harness: an unstripped `inputDATE.txt\n` "didn't exist".)
- **Don't pass paths containing `..`** to the ppMF/readFile layer — `os.path.exists`
  behaved inconsistently; use `os.path.abspath`.
- **Bash has a 45 s timeout and the reference `_h5_MM.h5` is 1.4 GB / `_h5_MF.h5`
  is 3.3 GB.** Never read a full dataset; slice `[0:nd]` (a few days). A full-array
  read wedged the shell during Phase 2.
- **git commits past the baseline fail** on this mount (lock-file permission). All
  changes save to the files regardless; snapshot with git on the real machine.
- `matplotlib`, `h5py`, `numpy` 2.x are available; `flopy`/`modflowapi` are **not
  yet installed** in the sandbox — `pip install --break-system-packages flopy
  modflowapi` when you start 3a.

## Do NOT touch

- `tests/legacy/MARMITESsoil_v3_legacy.py` — the frozen v0.3 Decimal oracle.
- The 20 existing tests must stay green through Phase 3 (`python3 -m pytest tests/ -q`).
- `MARMITESsoil_v3.flux()` / `_cell_step` physics — only the dry-cell/`hnoflo`
  handling changes for MF6, nothing else.

## Current test status

20 passed: 9 flux regression (new vs legacy), 5 runMMsoil/step integration,
6 config. Plus `tests/validate_lamata.py` (needs the La Mata `_h5_MF.h5`, present).
