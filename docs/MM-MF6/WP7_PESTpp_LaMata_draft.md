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
| 6 | observations | **C4, C5 and W1 used**. The eddy-covariance ET is **validation** only. ET for calibration from **MODIS MOD16A2GF**, as on CdL; downloaded to 2010-11 so far (§1) |
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

### The satellite ET: MODIS MOD16A2GF

- **MODIS MOD16A2GF** (500 m, 8-day), as on CdL (whose ET group is
  `CdL_pest\mod16_aet_zonal.csv`, zonal means per cover type). The
  Copernicus 300 m product starts in 2025 and cannot cover 2008–2013.
- Downloaded by the user into
  `DATA_ROOT\RS\MODIS4LaMata\MOD16A2GF_061-20261011_201206`: tile
  **h17v04**, the right one for La Mata (CdL's is h17v05).
- **The download stops in late November 2010:** 119 composites, 2008 DOY 121
  → 2010 DOY 329. The model runs to 2013-09-30, so about 135 composites are
  missing (December 2010 – September 2013), most of the calibration window.
  **OPEN:** complete the order, or calibrate ET on 2009-10 → 2010-11 only.
- The remote-sensing folder is a machine path, asked on Home (`rs`, default
  `DATA_ROOT\RS`), one sub-folder per product and case.
- At 500 m the catchment is ~20 pixels: catchment and per-cover-type means
  per 8-day composite, as on CdL, not pixel by pixel; the 8-day sums are
  compared with the model's total ET over the same 8 days.

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

## 5. The interface, as built (2026-10-11)

**Model side** (`code/marmites_pest.py`, no pyemu on the model path):
- `[pest] multipliers` names a folder of multiplier files, applied on the
  mesh after the 50 m rasters were projected: `k_l<L>`, `sy_l<L>` per cell,
  `ss_l<L>`, `drn_l<L>` one value per layer.
- K carries k33 and UZF's vks with it; Sy carries UZF's thtr.
- `set_param` sets a configuration value by name, list items included
  (`soil.horizon[2].ks`).

**Tools** (`code/tools/`):

| script | does |
|---|---|
| `pest_modis_et.py` | MOD16A2GF composites → the catchment's area-weighted ET per composite (35 pixels, 4.835 km²; fills left out; HDF4 read through rasterio's GDAL) |
| `pest_lamata_obs.py` | the targets (groups of §1) from a run's exports + MODIS + EC; their collection from any run; the group balancing |
| `pest_lamata_forward.py` | the forward run in a worker: scalars into the config, the coupled model, `pest_obs.csv` (a failed run writes the 1e10 sentinel, never a missing file) |
| `pest_lamata_setup.py` | `org` (code snapshot, config copy, dataset copy with C4/C5/W1 on, caches + saved state, parameter files, ET records) → `base` (the base run, the targets, the weights) → `pst` (PstFrom: pilot points on the mesh, the scalars, the observations, the prior ensemble) |
| `pest_lamata_run.py` | pestpp-ies with local PANTHER agents (at most 6: 2 cores stay with the user) |

**The PEST workspace** (`X:\3p1p1\MF6models\LaMata\PEST`, set on Home):
`code/` (the snapshot), `org/`, `template/`, `master/`, `workers/`.
Every realisation starts from `hi_voronoi_lamata_evt_drn1`, the end state of
the run `2lay_fullperiod_numba` (fractional grass).

**The template:**
- 215 parameters: 49 pilot points for each of K and Sy in both layers
  (400 m apart, exponential variogram, range 1500 m, log), plus 19
  scalars.
- Each scalar sits in its own `pars/<name>.csv`: PstFrom parameterises a
  file once in direct style.
- The prior ensemble is drawn from the geostatistical structure.

**Checked on a 20-day model** (scratch workspace):
- setup → base run → template;
- pestpp-ies `noptmax 0` through an agent reproduces the setup's phi exactly
  (8500);
- a prior Monte Carlo of 4 realisations on 2 agents: 4 runs, 0 failed,
  phi 8205–11371 around the base's 8500.

**Run it** (after the setup):

    python code/tools/pest_lamata_run.py --workers 6 --noptmax -1   # prior MC
    python code/tools/pest_lamata_run.py --workers 6 --noptmax 2    # pilot IES

**Settings in the control file:** `ies_num_reals 30`,
`ies_autoadaloc true`, `ies_bad_phi_sigma 2.5`,
`ies_lambda_mults 0.1,1,10`, `overdue_giveup_fac 1.5`, `save_binary true`.
After each iteration, check:
- the share of parameters on their bounds;
- the posterior Moran's I of the K fields;
- the failed realisations.

**Starting on a bound:** `uzf_eps` (3.5, MF6's minimum) and `crr_beta`
(1.0) start on a bound, so half of their prior draws sit on it.

---

## 6. The prior Monte Carlo (night of 2026-10-11/12)

`pest_lamata_run.py --workers 6 --noptmax -1`: 29 prior realisations plus
the base, 6 agents, **5.6 h**. A run took 53–64 min with six at once, 40
alone. Report: `PEST\master\report_0\` (`pest_lamata_report.py`).

- **28 ran, 2 failed: exactly the two realisations whose streambed K was
  clipped at its upper bound, 1.0 m/d.** Run 21 was killed as overdue
  (95 min against a 64-min average); run 11 failed after 82 min. Every
  realisation with rhk ≤ 0.65 m/d ran in normal time.
- **Phi:** base 10000 (the balanced shares); prior 8981–30416, mean 14162.
  Several realisations already beat the base.

| group | prior phi min–max | |
|---|---|---|
| hds_c (C1–C5, W1) | 1278–21206 | very sensitive to K and Sy |
| hds_p0 | 652–3157 | sensitive |
| sm_sm, sm_p0 | 1284–3729, 1638–2674 | moderate |
| q_year | 320–795 | moderate |
| q_log | 1403–1825 | weak |
| **et_rs (MODIS)** | **944–1125** | **almost insensitive** |

- **On a bound in the prior:** `uzf_eps` 64 % (it starts at MF6's minimum
  3.5), `crr_beta` 57 % (it starts at 1.0), `drn` 5 %.

### What the Monte Carlo shows

1. **The C1–C5 heads are not measurements.** Each well repeats exactly
   base + (0, +1.5, +2.0, +1.0 m) in Oct/Jan/Apr/Jul, identically every
   year from 2007 to 2012, and W1 is a constant 779.7 m. As drafted they
   carry half of the heads' share of phi (17.5 %) as if they were data.
2. **The ET group cannot be fitted with these parameters.** The model
   follows MODIS's seasons but misses the spring peaks (MODIS 17–27 mm per
   composite, the model 10–15), and the prior barely moves it. What
   controls ET (transpiration demand, phenology, the Tg coefficients,
   rooting) is fixed.
3. **`sfr_rhk` at 1.0 m/d makes the model slow, not wrong.** Realisation
   11 rerun alone (scratch) completes: 64 min against 40, 63 linear
   iterations a day against 39, 2 ATS sub-steps. Six at once, that passes
   PEST's overdue limit. Its phi, 10892, is no better than the base.
4. **`overdue_giveup_fac 1.5`** kills slow but valid runs when six run at
   once.

### Proposed before the pilot IES (for the user's decision)

- **C1–C5, W1:** one target per well, its mean level (level information,
  not a time series), or drop them; the heads' share stays with P0.
- **ET:** add the parameters that move it, or keep the 10 % and accept a
  structural misfit, or compare seasonal sums rather than 8-day
  composites.
- **`sfr_rhk`:** upper bound 0.5 m/d.
- **`uzf_eps` starts at 4.5, `crr_beta` at 0.95:** inside their ranges, so
  the prior no longer piles on a bound.
- **`overdue_giveup_fac 2.0`.**
- **The pilot:** `--noptmax 2`, restarting from this Monte Carlo's
  observation ensemble (no prior runs again): ~42 runs an iteration (30
  realisations + the lambda subset), ~14 h for two iterations at 6 workers.

### Decided (user, 2026-10-12) and launched

- **C1–C5, W1 are soft data:** each weighs 25 % of a measured P0 head
  (weight 0.48 against 1.94; noise 0.4 m against 0.1). P0 takes the heads'
  whole 35 %. The soft group's initial phi is then 6.5: the base already
  sits within about half a metre of those levels.
- **ET stays at 10 %, effectively validation;** no vegetation parameter
  varies.
- **`sfr_rhk` ≤ 0.5 m/d; `uzf_eps` starts at 4.5, `crr_beta` at 0.95;
  `overdue_giveup_fac` 2.0; `lambda_scale_fac` 1.0** (three trial upgrades
  an iteration).
- **The pilot was launched 2026-10-12 10:38** (`--noptmax 2`, 6 agents). It
  draws a FRESH prior of 30: restarting from the Monte Carlo would have kept
  that prior's eps/beta draws on their bounds and three realisations above
  the new rhk cap. Expected about 20 h: 31 prior runs, then per iteration
  12 trial runs on the 4-realisation subset plus 26 updates.
- The Monte Carlo's results are kept in `PEST\master_prior_mc`.
