# MARMITES — branch `MM-MF6_SFR_LAK_CRR`

MARMITES is a transient, distributed model of the land surface and the soil zone (recharge,
ET partitioning and ET sourcing), coupled to **MODFLOW 6** through the BMI/API. The model
itself is described in [`docs/MARMITES_overview.md`](docs/MARMITES_overview.md).

This branch adds the surface-water side — **SFR** (stream routing), **LAK** (ponds) and
**CRR** (cascade routing with reinfiltration) — together with a **three-source total ET**, a
**configuration file** in place of the command-line flags, a **Streamlit** interface, and a
**PEST++-IES** calibration.

- The plan: [`docs/MARMITES_x_CdL_merge_cookbook.docx`](docs/MARMITES_x_CdL_merge_cookbook.docx)
  (work packages WP0–WP8 plus WP1b).
- The running state: [`docs/MM-MF6/MARMITES_NEXT_STEPS.md`](docs/MM-MF6/MARMITES_NEXT_STEPS.md).
- The design record: [`docs/MM-MF6/MARMITES_SFR_LAK_CRR_analysis.md`](docs/MM-MF6/MARMITES_SFR_LAK_CRR_analysis.md).

## Layout

```
docs/       documentation
  MM-MF6/     the MODFLOW 6 conversion record: handoff, design analysis, phase reports
  MARMITES_x_CdL_merge_cookbook.docx     THE PLAN for this branch
  MARMITES_overview.md                   what the model is and does
code/       all the code
  MARMITESsoil/ MARMITESsurf/ MARMITESutilities/   the MARMITES model itself
  ppMF6/        the MODFLOW 6 side: build, SFR, LAK, layers, pre/post-processing
  ppMF_FloPy/   the legacy MODFLOW pre-processor (kept to read, not to run)
  SFR_LAK_CRR/  the CdL reference model — reference only, not part of the build
  tests/        the test suite and the run drivers (run_lamata_mf6.py)
  tools/        one-off utilities: GIS -> dataset conversion, the PEST chain    (WP1, WP7)
  app/          the Streamlit interface                                          (WP1b)
  configs/      run configuration files, one per model configuration             (WP0)
example/    the input data MM and MF read directly — one folder per case study
  LaMata/     La Mata (Spain) — the development and test case
  CdL/        Casa de Lobos (Portugal) — the production case      (to be populated)
  WRR/        the small WRR dataset a few tests use
```

## The repository / workspace rule

**The repository holds code, docs and strictly the input files MM and MF read directly.**
No GIS, no run output, however small or convenient. Everything a run produces goes to a
workspace outside the repository, set by `--ws-root` or `$MARMITES_WS_ROOT`
(default `E:\00code_ws\LaMata_MM-MF6`):

```
<WS_ROOT>/MF6_ws/                 the MODFLOW 6 model and its output
<WS_ROOT>/MMsurf_ws/              MMsurf output
<WS_ROOT>/out_<stamp>_<tag>/      MM results: _input/ (parameter maps) + _output/ (figures, CSV)
```

Raw geospatial sources (shapefiles, DEM) live outside the repository as well, under
`$MARMITES_DATA_ROOT` (default `E:\00code_ws\LAMATA_new`). A one-off converter turns them
into the small plain-text tables under `example/<case>/` that the model actually reads, which
is what keeps geopandas and rasterio off the model path.

## Running

```
python code/tests/run_lamata_mf6.py --libmf6 C:\00MODFLOW\mf6.7.0_win64\bin\libmf6.dll ^
  --mode lagged --nlay 2 --seep drn --strt-heads hi_spinup --steady-means hi_spinup ^
  --preproc --postproc
```

Redraw every figure from a run already on disk (no MODFLOW, about a minute):

```
python code/tests/run_lamata_mf6.py --postproc-only --preproc --nlay 2 --mode lagged --run-tag <tag>
```

WP0 of the cookbook replaces those flags with `--config code/configs/<case>.toml`.

## Tests

```
python -m pytest code/tests -q
```

Six tests need an `mf6` binary on PATH and one needs `pyshp`. A few plotting tests can crash
natively in a headless shell with some matplotlib/flopy combinations — a pre-existing
environment issue, not a code fault; run the suite from the same environment as the model.

## Relation to the other branches

Cut from `MM-MF6` at the point where the MODFLOW 6 conversion was complete and validated
(192 passed / 7 skipped / 0 failed, 2026-09-09). Relative to that branch this one is
**reorganised, not rewritten**:

| was | is |
|---|---|
| `trunk/` | `code/` |
| `tests/` | `code/tests/` |
| `doc/` | `docs/` |
| `DataSet_LaMata/` | `example/LaMata/` |
| `DataSet_WRR/` | `example/WRR/` |

Dropped here and kept on `MM-MF6` and in the history: the committed Python-2 era `venv/`,
the JetBrains and PyScripter IDE settings, and the vestigial root `__init__.py` (a
flopy-derived stub that nothing imported and that made the repository root masquerade as a
package).
