# WP7 — PEST++-IES on La Mata: setup

*Draft, 2026-10-11, with the user's decisions of the same day. Open points
are marked **OPEN**.* It follows the cookbook's WP7 (§7.1–7.9, design §8)
and corrects it where the model has moved on: the run time, the
observations that exist, the vegetation cover and the ET data.

Baseline, the full-period run `20261011113811_2lay_fullperiod` (1949 days,
2008-05-31 → 2013-09-30):

| target | fit now |
|---|---|
| heads P0 (near the outlet drains) | RMSE 1.21 m, r 0.58 |
| heads C1, C2, C3 | RMSE 0.69 / 0.55 / 0.54 m, r 0.78–0.85 |
| soil moisture P0, horizons 1–2 | RMSE 2.0 / 2.2 % |
| soil moisture SM, horizons 1–2 | RMSE 3.3 / 6.1 % |
| outlet streamflow | 194 mm/yr; NSE −5.0, **bias +33 %** |
| UZF | rejects 30 % of the percolation |

---

## Decisions (2026-10-11)

| # | question | decision |
|---|---|---|
| 1 | vegetation cover | grass from the **legacy fractional cover** `inputVEG1area.asc` (`[[soil.veg_class]] cover_raster`); the trees from the crown polygons. Not the outcrop cut. |
| 2 | "amplitude ×10" | multiplier **0.1–10** |
| 3 | Ss | **one value per layer** |
| 4 | DRN | **one conductance per layer** |
| 5 | additions | **UZF eps and vks, CRR β**: yes. Q. ilex root depth and soil thickness: no, measured |
| 6 | observations | **C4, C5 and W1 used**. The eddy-covariance ET is **validation** only. ET for calibration from a satellite product, as on CdL (see §1, **OPEN**) |
| 7 | weights | heads **35 %**, soil moisture **35 %**, streamflow **20 %**, ET **10 %** |
| 8 | speed | profiled, then **numba** for the groundwater-ET curves and the soil step: **~45 min a full run, was 2 h 20** (§4) |
| 9 | workers | **this server only** |
| 10 | workspace | asked on **Home** (`pest_ws`), with the size warning (CdL's PEST folder is ~40 GB) |

---

## 0. The vegetation cover

Interception was 29.3 mm/yr in MF6 against 16.2 in the NWT reference over
the same 500 days. The cause was the grass cover:

| | legacy NWT (`inputVEG*area.asc`) | crown polygons (`lm_veg.shp`) |
|---|---|---|
| grass | **25.9 %** | **91.2 %** |
| Q. ilex | 7.6 % | 7.8 % |
| Q. pyrenaica | 0.8 % | 0.8 % |
| bare soil / rock | 65.8 % | 0.2 % |

The polygon layer is the tree crowns inside grass polygons that fill
everything else. The legacy grass raster is a **fractional cover**: on
average 39 % in alluvium cells, 24 % on regolith, 5.5 % on the outcrops, and
32 % even where there are no trees. Cutting the outcrops (4.5 % of the
catchment) would only have taken grass to 87 %.

**Done:** `[[soil.veg_class]]` has a `cover_raster` column (Soil panel,
vegetation table). The grass row names `inputVEG1area.asc`; its cover is put
onto the mesh by exact area and kept within what the crowns leave. The
trees still come from the polygons, and they match the legacy rasters to
0.02 points.

**Summer wilting is accounted for.** A vegetation type counts only while its
LAI > 1e-5; grass has `lai_dry = 0`, so from late June to early October
(112 days a year) its area is bare soil: no interception, and soil
evaporation instead of transpiration. The rain in that season is 72 mm/yr
(16 % of P); its interception would have been 3.7 mm/yr.

---

## 1. What we calibrate against

Counts are in the calibration window of §3 (2009-10-01 → 2013-09-30).

| group | data | obs | in the forward run | share of phi |
|---|---|---|---|---|
| `hds_p0` | P0 head, near-daily | 496 | head at the period end | 35 % for all heads: P0 half, the quarterly wells half |
| `hds_c` | C1–C5 quarterly to 2012-07; W1 | 5 × 12 + 2 | head at the period end | (in the 35 %) |
| `sm_p0`, `sm_sm` | θ at P0 and SM, 2 horizons, daily | 2 × 1267, 2 × 959 | θ at the period end | 35 %, the two probes half each |
| `q_log` | outlet streamflow, daily | 725 | log10(Q + 0.01 mm/d) | 15 % |
| `q_year` | outlet volume per hydrological year | 4 | constrains the water balance (bias +33 % now) | 5 % |
| `et_rs` | satellite actual ET, catchment and per cover type | per composite | the period mean of total ET (Ei + Eow + ETsoil + ETuzf + ETg) | 10 % |
| `et_ec` | eddy-covariance ET at `EC`, daily, 2009-08-01 → 2010-10-01 | 266 | total ET at the EC cell | **0: validation**, reported never fitted |

- C4, C5 and W1 are commented out in `inputObs.txt`; the PEST setup adds
  them to its own copy of the observation table, so the model writes their
  simulated heads.
- The eddy-covariance series is in
  `DATA_ROOT\MMsurf_ws\Eddy_MM_202110.xlsx`, sheet `compil4python`, column
  `ET_EC`: 266 days, mean 1.25 mm/d.
- Measurement noise is set apart from the weights (`standard_deviation`:
  heads 0.1 m, θ 0.01, log flow 0.1, ET 20 %), so rebalancing the shares does
  not change the noise the IES realisations get.

### The satellite ET — **OPEN**

- **The Copernicus product named** (Evapotranspiration 2025–present, 300 m,
  10-daily, Sentinel-3) **starts in 2025**. It cannot overlap La Mata's
  2008–2013, and the product page lists no older version.
- **CdL did not use Copernicus:** its ET group is MODIS **MOD16A2GF**
  (500 m, 8-day, 2000–present), zonal means per cover type in
  `CdL_pest\mod16_aet_zonal.csv`, from the tiles in
  `CdL_model\RS\MODIS4CDL` (`CopernicusET4CDL` is empty).
- CdL's tiles are **h17v05** (30–40°N). La Mata (41.1°N, 6.1°W) is in
  **h17v04**, so it needs its own download: 6 years × 46 composites ≈ 276
  files, ~3.3 GB as whole tiles, or a few MB as a catchment subset through
  NASA AppEEARS.
- Either way the download needs a **NASA Earthdata login**, which only you
  can use. You can download the files (or the AppEEARS subset) yourself, or
  set up the login on this machine (an Earthdata `.netrc` or token) for a
  download script to use.
- At 500 m the catchment is ~20 pixels: catchment and per-cover-type means
  per 8-day composite, as on CdL, not pixel by pixel.

---

## 2. Parameters

Bounds are on the multiplier of today's value unless stated; log-uniform
priors, standard deviation = range/4 in log space.

| group | what | style | bounds | count |
|---|---|---|---|---|
| `k_l1`, `k_l2` | horizontal K (today 3 zones: 0.05 / 0.1 / 0.5 m/d) | pilot points on the existing zones | ×0.1–10 | ~30 per layer |
| `sy_l1`, `sy_l2` | specific yield (today 0.01) | pilot points | ×0.1–10 (Sy 0.001–0.1) | ~30 per layer |
| `ss_l1`, `ss_l2` | specific storage (~1e-4 L1, ~4e-7 L2) | one per layer | ×0.1–10 | 2 |
| `drn_l1`, `drn_l2` | outlet drain conductance (today 1e-4 m²/d per m; legacy 7e-4 / 5e-4) | one per layer | ×0.1–100 | 2 |
| `sfc_*` | field capacity as its position f = (sfc − sr)/(smax − sr), per horizon | additive on f | f ± 0.10: sfc stays inside (sr, smax) | 5 |
| `ks_*` | soil saturated K per horizon (1.5–5 mm/d) | multiplier | ×0.1–10 | 5 |
| `sfr_rhk` | streambed K (0.1 m/d; `rbth` fixed, fully correlated) | multiplier | ×0.01–10 | 1 |
| `lak_bedleak` | pond-bed leakance (1e-3 1/d) | multiplier | ×0.1–10 | 1 |
| `uzf_eps` | Brooks-Corey exponent (3.5, MF6's minimum) | value | 3.5–8 | 1 |
| `uzf_vks` | UZF vks multiplier (UZF rejects 30 % of the percolation now) | multiplier | ×0.2–5 | 1 |
| `crr_beta` | share of the runoff that moves on in the cascade | value | 0.8–1.0 | 1 |

About 140 parameters. Ensemble methods cost per realisation, not per
parameter. With 7 head wells, K and Sy away from them stay near their prior,
and IES reports that as uncertainty.

**Pilot points on the Voronoi mesh:** ~400 m apart (about 30 inside the
catchment per layer), at mesh-cell centres. pyEMU 1.4 takes an unstructured
grid as a dictionary of cell centres; this is checked on a toy model first.
Exponential variogram, range 1200–1600 m (3–4 × the spacing, as on CdL). The
prior ensemble is drawn from that structure, never from bounds alone: on CdL
the bounds-only prior gave the checkerboard K and the collapse.

**Fixed:**
- measured: Q. ilex root depth, soil thickness;
- numerical: `seep.cond`, which must stay free-draining;
- no information: `sfr.rbth` (duplicates `rhk`), `sfr.manning`, `k33` (no
  layer-2 head data), `canopy_mm` (applied in MMsurf), the Tg coefficients
  (already calibrated).

---

## 3. Time window and initial state

- Every realisation starts from the same saved state, the end of the
  full-period run (`hi_voronoi_lamata_evt_drn1`).
- A parameter set moves the equilibrium, so **the first hydrological year is
  warm-up** (2008-05-31 → 2009-09-30, weight 0) and the objective covers
  **2009-10-01 → 2013-09-30**, four hydrological years.
- The forward run simulates all 1949 days.

---

## 4. Run time and numba

**Before:** the full-period run took about **4 s per day, 2 h 20 min**,
twelve times the cookbook's 12 minutes.

**Profile** (30 days in a scratch workspace; the profiler roughly doubles
Python's own time):

| part of a day | profiled | share |
|---|---|---|
| groundwater-ET curves, built per cell in Python (`_evt_curves`) | 3.5 s | 43 % |
| soil column + CRR cascade, per cell in Python (MMsoil) | 3.1 s | 39 % |
| MODFLOW 6 (compiled) | 1.4 s | 17 % |

**Done, both compiled with numba:**
- `marmites_evt.curves` builds every cell's Eg and Tg curves at once:
  0.023 s a day instead of 1.6 s (70×).
- `marmites_soil_nb` runs the soil step and the cascade in one pass. The
  Python code (`clsMMsoil.step`) stays the reference: it runs without numba,
  with `MM_SOIL_PYTHON` set, or on the uncoupled NWT path.
- The tests hold the compiled code to the reference to 1e-12, cell by cell:
  every output, the cascade, the carried soil state, the tallies and the
  warnings, on toys built to hit every branch.
- On La Mata (30 days) the results match the old code to round-off: heads
  within 8e-8 m, every flux within 6e-8 relative, identical solver
  iterations.

**After: about 1.35 s per day**, so a full run takes **~45 min** (the first
run of a session adds ~15 s of compilation). MF6 is now most of the day; its
solver settings are the next lever, if needed. At 6 workers (8 physical
cores, 2 left to you):

| step | runs | at ~0.8 h per run |
|---|---|---|
| base run (`noptmax = 0`) | 1 | < 1 h |
| pilot: 30 realisations × (prior + 2 iterations) | ~95 | ~13 h, overnight |
| production: 100 realisations × (prior + 3 iterations) | ~400 | ~2–2.5 days |

How much 6 simultaneous runs slow each other down is measured on the pilot.
Numba caches the compiled code next to the module, so the setup compiles it
once before the workers start.

---

## 5. Build steps

Everything in the PEST workspace set on Home (`pest_ws`; `org/`,
`template/`, `master/`, `workers/`), on a copy of the configuration, never in
WS_ROOT and never on `code/configs/lamata.toml`. **Each worker needs its own
copy of the dataset folder** (27 MB): a run rewrites the aggregated forcing
files (`inputZON_*_stp.txt`) there.

1. **Model hooks** (code, with tests):
   - a `[pest]` multiplier directory: one file per property and layer, one
     value per line per mesh cell, applied to K and Sy after they are
     resampled onto the mesh;
   - a scalar table (`name, value`) applied to the configuration, including
     list items (`soil.horizon[i].sfc`, `.ks`), which `--set` cannot address;
   - the ET exports: total ET per cover type over the catchment, for the
     `et_rs` group.
2. **`code/tools/pest_build_lamata.py`**: `PstFrom`, pilot points at the
   mesh-cell centres, the geostatistical prior ensemble, the observation
   groups and weights of §1 from the base-run exports.
3. **`code/tools/pest_forward_run.py`**: apply the parameters → write the
   configuration → run `run_lamata_mf6.py` (MMsurf off, figures off) →
   transform (log flow, annual volumes, ET means per composite) → write the
   observation files. Guarded by `multiprocessing.freeze_support()`.
4. **Base run** (`noptmax = 0`): its phi must match the residuals computed by
   hand.
5. **Pilot IES**, then production.
6. **`[pest] use_pest_params`** (declared, not wired yet): a normal run with
   the calibrated set.

**pestpp-ies settings:** `ies_num_reals` 30 (pilot) / 100, `noptmax` 2 / 3,
`ies_autoadaloc true`, `ies_bad_phi_sigma 2.5`, `ies_lambda_mults 0.1,1,10`,
`overdue_giveup_fac 1.5`, `save_binary true`, 6 PANTHER workers. After each
iteration, check:
- the share of parameters on their bounds;
- the posterior Moran's I of the K fields;
- the failed realisations (< 30 %).
