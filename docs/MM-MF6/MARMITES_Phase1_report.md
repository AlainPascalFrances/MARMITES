# MARMITES — Phase 1 completion report

**Date:** 2026-07-21
**Model:** Claude Opus 4.8
**Scope:** Phase 0 (safety net) + Phase 1 (bug fixes, Decimal removal, Python 3.12 port, strip NWT + Picard loop) of the agreed plan in `MARMITES_code_review.md` §6.

---

## Result

All Phase-0 and Phase-1 tasks are complete. Core modules compile and import under Python 3.12 / numpy 2.x; the soil-flux regression suite (9 tests) passes; the MODFLOW-NWT model construction, the initial MF run, and the MM↔MF Picard convergence loop are removed. The code is now a clean base for the Phase-3 MF6/API backend.

Net change: **+1,107 / −3,755 lines** across the trunk (16 files).

---

## Phase 0 — safety net

- Baseline committed to git (`10dc0dc`) before any change. *(Note: subsequent commits are blocked by a filesystem-permission quirk on the mounted drive in this environment; all changes are saved to the files themselves. Re-run `git add -A && git commit` on your own machine to snapshot Phase 1.)*
- Legacy/dead code moved to `MARMITES/_archive/`: `MM_XLSstuff/`, `pyEARTH1D/` (wxPython GUI), `MARMITESsoil_v3_BAK20221216.py`, `Tg.py`, `Tg4Maciek.py`, `TgNEW_v*.py` (8 files).
- `pyproject.toml` created: Python ≥ 3.12, numpy ≥ 2, matplotlib ≥ 3.9, h5py, flopy ≥ 3.9, modflowapi; ruff + pytest configured.
- Regression reference built as **MM soil-column unit tests** (no NWT baseline, per your decision to delete NWT): `tests/test_soil_flux.py`, with the frozen legacy soil module preserved at `tests/legacy/MARMITESsoil_v3_legacy.py`.

## Phase 1 — bug fixes (§2 of the review)

All ten fixed:

1. `Rexf_tmp /= perlen` on a list → now a numpy array; the **exfiltration feedback path (EXF>0) runs** (was crashing). Covered by `test_new_exf_mass_balance`.
2. `if CROP_tmp != None` / array truth-value ambiguity in the irrigation branch → scalar crop id + `is not None`.
3. `for i in range(self.row)` typo → removed with `runMF` (the only site).
4. `self.versionsys.exit()` garbled line → removed.
5. `mpl.get_backend != 'agg'` → `mpl.get_backend() != 'agg'`.
6. `global rp_tmp` / silent stale-value returns in `perc()`/`evp()` → pure static methods `_perc`/`_evp` with `np.isfinite` guards raising `MarmitesSoilError`.
7. `wel_dum` early-break well logic → removed with `runMF`.
8. Bare `except:` → `MarmitesError` introduced; all active bare excepts in the run path converted to typed/`except Exception:` (the input-parsing and observation ones narrowed to `TypeError/ValueError/OSError/IndexError`).
9. Module `global`s (`MARMITESprocess`, `MARMITESsoil`) → removed; state passed explicitly. Fixed a latent bug where `verifObs` returned the previous point's `obs_yn` when a file was missing.
10. Unreachable `del`-after-`return`, duplicated lines, `raise sys.exit()` → cleaned.

Plus a pre-existing latent bug found during import testing: `mpl.dates` used without `import matplotlib.dates`, evaluated in an import-time default argument — fixed in `MARMITESutilities`.

## Phase 1 — Decimal removal (§3)

All 63 `Decimal(...).quantize(...)` sites removed across `MARMITESsoil_v3.py`, `MARMITESprocess_v3.py`, `startMARMITES_v3.py`. Arithmetic is now float64 with `np.isclose` tolerances. `MARMITESsoil_v3.flux()` was rewritten around numpy arrays (the inner per-cell loop no longer allocates Decimals) — this alone should give a large speedup and is the groundwork for the Phase-2 vectorization.

## Phase 1 — Python 3.12 / numpy 2.x

- `!= None` / `== None` → `is`/`is not None` (elementwise-array hazard) across the core.
- bytes dict/HDF5 keys (`b'iP'`, …) → str, in driver and `MARMITESplot_v3.py` (264 sites); cbc **record names kept as bytes** (flopy returns bytes) — only MARMITES' own index keys were converted.
- `np.float`/`np.int` aliases: none remained in the run path (were only in archived utilities).
- flopy `sys.path` hack removed; UTF-8 BOM stripped from the driver; `plt.register_cmap` → `matplotlib.colormaps.register`.

## Phase 1 — strip NWT + Picard (§6, action 6)

- `startMARMITES_v3.py`: removed the initial MODFLOW run block, the whole `while` convergence loop, `h_diff*` bookkeeping, the ConvLoop and HEADSmaxdiff plots. MMsoil now runs **once**, reading heads/exfiltration from an existing `_h5_MF.h5` (a clearly-marked temporary bridge until Phase 3 supplies heads via the API).
- `ppMODFLOW_flopy_v3.py`: `runMF()` (≈440 lines of flopy.modflow NWT/mf2005/UZF1/WEL construction and the h5 harvest) removed; the mf2005-vs-mfnwt error branches removed. The class now only parses the MF ini, imports arrays, and computes stress-period discretization (`ppMFtime`).

---

## Verification

- `python3 -m py_compile` on all six core modules: OK.
- Runtime import of `MARMITESsoil_v3`, `MARMITESutilities`, `MARMITESprocess_v3`, `ppMODFLOW_flopy_v3` under numpy 2.2: OK.
- `ruff check --select F,E9,E722,E711,E721` on the five core computation modules: **All checks passed** (remaining ruff findings are style-class: UP031 printf-format, B007/B905 — deferred to Phase 2).
- `pytest tests/` : **9 passed**. Tests confirm legacy-vs-new agreement on the working paths (dry, rain/percolation, runoff, shallow-WT Eg/Tg, deep-WT no-Eg), that the legacy EXF path raised `TypeError` while the new one closes mass balance to < 1e-6 mm, and that the column saturates correctly under large exfiltration.

---

## Known limitations handed to Phase 2/3

- MMsoil still requires a pre-existing `_h5_MF.h5`; end-to-end running resumes only when the Phase-3 MF6 backend exists. This is the intended, documented gap.
- `runMMsoil` is still the whole-run `(nper, nrow, ncol)` double loop. Phase 2 must (a) flatten to `(nper, ncpl)`, (b) vectorize the inner loop, and (c) expose a per-SP `step(n, heads, exf) -> (finf, etg, state)` interface — the hinge for the API coupling.
- Dry-cell logic still uses the `hdry`/`hnoflo` sentinels (`h_MF_ini - cMF.hdry`); Phase 3 replaces these with `idomain` + `h < botm` under MF6.
- Style-class ruff findings (printf formatting, zip-strict, loop vars) remain; cosmetic, safe to sweep in Phase 2.
- `MARMITESplot_v3.py` was updated only for the str-key change; its `nrow×ncol` imshow maps are a Phase-4 (DISV) concern.

## For the next model (Phase 2)

Start by reading `MARMITES_code_review.md` (the contract) and `tests/test_soil_flux.py` (the behavioral spec). The regression tests must keep passing through the Phase-2 refactor — they are the guardrail that the vectorization preserves the soil physics. `tests/legacy/MARMITESsoil_v3_legacy.py` is the frozen v0.3 oracle; do not modify it.
