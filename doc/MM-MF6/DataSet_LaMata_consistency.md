# DataSet_LaMata — consistency check against the code

**Date:** 2026-07-21

## 1. Code fix made

The Phase-2 config converter was aligned to the **deleted DataSet_WRR** ini,
which lacked the `maxYearsTickTrimester/Semester` fields. LaMata's ini
**does** include them (matching the driver's canonical field order), so the
converter was misaligned by two fields on LaMata. Fixed: both fields
restored to `marmites_config._LEGACY_ORDER`. LaMata's `__…__inputMM_v3.ini`
now converts correctly (irrigation block, coordinates, all fields verified),
and the config test suite was repointed from WRR to LaMata. **20/20 tests pass.**

This is the positional-parsing fragility the TOML migration exists to kill:
two datasets, two different field counts, silent drift. `MMConfig.validate()`
now guards against it (a misaligned parse trips the `plt_WB_unit`/month checks).

## 2. Internal consistency — GOOD

- All rasters (MM_asc and MF_asc) are **60 cols × 65 rows, cellsize 50 m**,
  origin (739300, 4553050) — consistent with each other and with both ini files.
- MF model: **6 layers**, `itmuni=4` (days), `lenuni=2` (m), `nper=10`,
  `version=mfnwt`. Daily driving series `inputDATE.txt` has ~1,949 days from
  2008-05-31 (~5.3 years).
- The MMsurf-generated stress-period series (`inputZON_*_stp.txt`) and
  `inputDATE.txt` are present, so MMsurf has already been run (consistent with
  `MMsurf_yn=0` in the ini) — MMsoil can proceed without re-running MMsurf.

## 3. Layout mismatch — this is an ARCHIVED OUTPUT, not a runnable input tree

The `__<timestamp>__` filename prefixes and the `MM_asc/` + `MF_asc/` folders
are exactly the driver's **output-copy convention**. The code, however, expects
a different input layout:

| Code expects (relative to `MM_ws`) | Present in dataset | Status |
|---|---|---|
| grids at `MM_ws/` root (`inputSOILzones.asc`, …) | `MM_asc/…` | moved |
| `MM_ws/MF_ws/` with MF rasters + MF ini + `inputSOILparam.txt` | `MF_asc/…`, ini at root (timestamped) | moved / renamed |
| `MM_ws/MMsurf_ws/` (`__inputMMsurf.ini`, `__meteoTB.txt`, `__IRR_TS.txt`) | **absent** | missing (but MMsurf outputs present) |
| `SOILparam_fn = MF_ws\inputSOILparam.txt` | only `__…__inputSOILparam.txt` at root | renamed |
| `_h5_MF.h5` (heads/exfiltration for MMsoil) | **absent** | needs an MF run |

To run MMsoil from this dataset, a small restructuring is needed:
1. `MM_ws/` = a folder holding the contents of `MM_asc/` (grids + MMsurf outputs + obs).
2. `MM_ws/MF_ws/` = contents of `MF_asc/` + the MF ini + `inputSOILparam.txt`
   (strip the `__timestamp__` prefixes; fix `SOILparam_fn` path).
3. `MM_ws/MMsurf_ws/` only needed if re-running MMsurf; since its outputs are
   already present and `MMsurf_yn=0`, it can stay empty.
4. `_h5_MF.h5` must come from a MODFLOW run — see §4.

## 4. MODFLOW dependency — expected Phase-3 gap

The MF ini is MODFLOW-**NWT** (`version=mfnwt`, hardcoded
`C:\00MODFLOW\MODFLOW-NWT_1.1.4\…NWT_64.exe`, 122 tokens of NWT/UZF1 package
config). Since Phase 1 removed `runMF`, `ppMODFLOW_flopy_v3` now only *parses*
this ini and computes the stress-period discretization (`ppMFtime`) — it no
longer builds or runs MODFLOW. So:
- The MF ini is consistent enough to drive time discretization today.
- Its package/exe content is NWT-era and is **superseded by the MF6
  configuration in Phase 3**; `_h5_MF.h5` will be produced by the MF6/API
  backend, not by this NWT ini.
- Practically, MMsoil cannot run end-to-end on LaMata until the Phase-3 MF6
  backend exists (or a legacy NWT run supplies `_h5_MF.h5`). This is the
  documented Phase-1/2 gap, not a dataset defect.

## 5. Portability note

`SOILparam_fn` (and `inputFile_TSirr_fn`) use Windows backslash separators
(`MF_ws\inputSOILparam.txt`). On non-Windows targets `os.path.join` will not
split these. The Phase-3 config should normalize path separators when it takes
over MF configuration.

## Bottom line

The dataset is **internally consistent** (grids, dimensions, time series all
agree) and the **converter now parses it correctly**. It is packaged as an
archived output snapshot, so it needs light restructuring into the expected
`MM_ws / MF_ws` layout to be runnable, and — like everything else — it waits on
the Phase-3 MF6 backend to produce `_h5_MF.h5`. I can do the restructuring into
a clean runnable input tree on request.
