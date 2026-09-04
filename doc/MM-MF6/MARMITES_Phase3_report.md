# MARMITES — Phase 3 completion report (code side)

**Date:** 2026-07-21
**Model:** Claude Fable 5
**Scope:** Phase 3 of `MARMITES_code_review.md` §6 — MODFLOW 6 backend (3a),
API coupling driver in lagged (3b) and iterative (3c) modes, MF6 dry-cell
semantics. All decisions of the kickoff honored: UZF6 always, no Picard loop,
daily stress periods, DIS now / DISV in Phase 4, NWT not reintroduced.

## Status: all Phase-3 code is written and tested. One step remains on your machine — the coupled binary run — because the sandbox cannot download MF6 executables (GitHub CDN and conda are blocked by its proxy). Everything up to the binary call is verified: the MF6 simulation builds and round-trips through flopy on the real La Mata model, and both coupling modes are fully exercised against a mocked MF6 API.

---

## What was built

### 1. `trunk/ppMF6/marmites_mf6.py` — `clsMF6` (Phase 3a)

Builds the MF6 simulation with `flopy.mf6` per the §4.1 mapping:

- **TDIS**: steady initial SP (the `dum_sssp1` convention) + transient SPs;
  daily by default (§4.3). On La Mata, daily discretization gives exactly
  **1949 SPs / 1949 days** (obtained by `perlenmax=1` through the existing
  `ppMFtime` — no new code path).
- **DIS**: 6×65×60, `idomain` from |ibound| (no ibound<0 on La Mata; a guard
  raises if constant heads appear — CHD would then be needed).
- **NPF**: `icelltype = laytyp`; `layvka=1` semantics honored: stored VKA is
  the hk/vk **ratio**, so `k33 = hk / vka` (unit-tested).
- **STO**: `iconvert=laytyp`, steady SP0, transient after.
- **IMS**: COMPLEX (from the NWT `options`), `outer_dvclose = headtol` (0.05),
  Newton with under-relaxation on the GWF model.
- **WEL**: one well per active surface cell (1,954), q=0 initially,
  **AUTO_FLOW_REDUCE 0.05** — replaces the NWT `SPECIFY 0.05 iunitramp` hack.
- **DRN**: the legacy 5-tuple list mapped to `cellid` form (count asserted).
  GHB skipped (La Mata `ghb_yn=0`), guard present for datasets that have it.
- **UZF6**: **vertical columns of UZF objects** — the first 1,954 objects are
  the land-surface cells (landflag=1) in MARMITES cell order, deeper active
  cells chain via `ivertcon` (landflag=0). Total **11,472 UZF objects = every
  active cell**, i.e. the full unsaturated column NWT-UZF1 handled internally
  is reproduced explicitly. `iuzfopt=2` honored: vertical K taken from layer
  k33. ET off in UZF (MARMITES does ET). Budget file wired for
  `UZF-GWRCH`/`UZF-GWD` terms.

### 2. `trunk/marmites_coupler.py` — `MF6Coupler` (Phases 3b/3c)

One exchange layer, two modes, API object injectable (that's what makes it
testable without the binary):

- **Exchange (pointers, no files):** MM→MF6: `perc [m/d] → UZF FINF`
  (land cells; deeper set 0), `ETg → WEL Q = −ETg·area [m³/d]`.
  MF6→MM: heads from `X` (with `NODEUSER` reduced-node mapping when idomain
  reduces the grid; identity fallback), exfiltration from `UZF GWD [m³/d]`
  → `/area·conv_fact` → **mm/d positive into the soil** (sign test included).
- **`lagged`**: heads/exf from the END of SP n−1 → `step()` → write fluxes →
  advance MF6 one SP.
- **`iterative`**: MM re-evaluated at every outer iteration from the current
  head iterate; under-relaxation `relax·new + (1−relax)·prev` on both fluxes;
  the soil state is snapshot-restored so exactly **one** state advance
  survives per SP regardless of iteration count.
- Steady SP0 driven with `perc_user` and zero ETg before the transient march.

### 3. MF6 dry-cell semantics in MMsoil (the `NOTE MF6` marker)

`cMF.hdry = None` now switches `_cell_step` from the NWT sentinel test to the
MF6 rule **dry ⇔ h < botm of layer 1** (fallback head = botm, as before).
Backward compatible: all legacy-comparison tests still pass with `hdry=1e30`.

### 4. `tests/run_lamata_mf6.py` — the entry point for your machine

Replicates the driver setup (same code as the validated harness), builds the
MF6 workspace `DataSet_LaMata/MF6_ws/`, and runs the coupled model when given
`--libmf6`. Verified here up to the binary call.

---

## Verification (all in-sandbox)

- **30/30 tests pass** (was 20): +4 MF6-build tests on real La Mata
  (write + flopy reload round-trip, UZF land-first ordering and `ivertcon`
  chaining into the same column, WEL/DRN counts, k33 ratio semantics) and
  +6 coupler tests against a mock API:
  - lagged mode demonstrably consumes SP n−1 heads (deterministic head
    evolution in the fake);
  - FINF/WEL pointer contents equal `step()` outputs (and Q = −ETg·area);
  - iterative mode: N evaluations per SP, relaxation trace, and **exactly one
    state advance per SP** (iterative(relax=1, static heads) ≡ lagged);
  - exfiltration sign/scaling: 25 m³/d on a 100×100 cell → +2.5 mm/d into soil.
- ruff error classes: clean. All modules compile.

## What remains — one command on your machine

```bat
:: once: binaries (in the MARMITES python env)
pip install flopy modflowapi
get-modflow :flopy            :: downloads mf6.exe + libmf6.dll

:: smoke test, ~60 SPs, lagged
python tests\run_lamata_mf6.py --libmf6 <path-to>\libmf6.dll --mode lagged --nsp 60

:: full runs
python tests\run_lamata_mf6.py --libmf6 <path>\libmf6.dll --mode lagged
python tests\run_lamata_mf6.py --libmf6 <path>\libmf6.dll --mode iterative --relax 0.6
```

Results land in `DataSet_LaMata/MF6_ws/_coupled_<mode>.h5` (heads, exf, perc,
ETg, outer-iteration counts per SP). Acceptance per §6: mass-balance closure,
water-balance components within tolerance of the published La Mata results
(plausibility, not bit-match — coupling scheme and UZF generation changed),
and lagged vs iterative converging as SPs are daily.

Two things to watch on the first real run (flagged in advance):
1. **MF6 variable addresses.** The coupler resolves `X`, `FINF`, `GWD`,
   `BOUND`, `NODEUSER` defensively, but exact address paths can vary by MF6
   version. If it raises `CouplingError: cannot resolve...`, send me the error
   text — `mf6.get_input_var_names()` output makes the fix a one-liner.
2. **Iterative-mode oscillation.** If outer iterations hit the cap at
   exfiltration cells, lower `--relax` (0.5 → 0.4). The lagged run is the
   robust reference either way.

## Deferred / notes

- The export/plot pipeline still reads the legacy `_h5_MF.h5` layout; wiring
  the coupler results into it (and the cbc-name updates `STO-SS/SY`,
  `FLOW-JA-FACE`) is the natural next slice after the first successful
  coupled run, together with the driver/export split deferred from Phase 2.
- `marmites_config` MF6 section: the runner takes CLI flags for now; fold
  into TOML once the run settles.
- The daily run of `ppMFtime` overwrites the `inputZON_*_stp.txt` files in
  `DataSet_LaMata/` with daily series (regenerated on every run; the
  aggregated versions are rebuilt automatically if you run in aggregated mode).
- Phase 4 (DISV) hooks are in place: cell-order contract (`surf_cells` ==
  `build_cell_list`), grid-agnostic cellids, and the geometry TODOs marked.
