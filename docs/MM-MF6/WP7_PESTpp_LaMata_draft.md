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
| `et_pt` | measured ET at the observation points (the EC tower; lysimeters when added), daily; EC 2009-08-01 → 2010-10-01 | 266 | total ET at the point's cell | **0: validation**, reported never fitted |

- C4, C5 and W1 are commented out in `inputObs.txt`; the PEST setup adds
  them to its own copy of the observation table, so the model writes their
  simulated heads.
- The measured ET at the points is read by the model run itself
  (2026-10-15), like the heads and the soil moisture: the dataset file
  `<obs.aet_prefix>_<point>.txt` (State variables panel; `inputObsAET`), one
  per point of `inputObs.txt`, which gives the coordinates. The run writes it
  beside the simulated ET in `obs_et.csv` (`obs_et_mmd`), and the targets take
  it from there for every point that has one. For EC this is
  `example\LaMata\inputObsAET_EC.txt`: 266 days, mean 1.25 mm/d, the
  values of the original `Eddy_MM_202110.xlsx` (sheet `compil4python`,
  column `ET_EC`). The user's spreadsheet copy is
  `DATA_ROOT\MMsurf_ws\inputObsAET_EC.txt`/`.xlsx`. A lysimeter is
  added the same way: a point in `inputObs.txt` and its series file.
  - The production run of 2026-10-13 still has the older group `et_ec`
    (read from the Excel sheet at setup): report it with the snapshot's
    `PEST\code\tools\pest_lamata_report.py`.
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
| `pest_lamata_obs.py` | the targets (groups of §1) from a run's exports + MODIS; their collection from any run; the group balancing |
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

---

## 7. The pilot IES (2026-10-12 10:38 → 10-13 06:28)

A fresh prior of 30, then 2 iterations; 6 agents, 10 from 11:25. 19.8 h,
0 failed runs. Report: `PEST\master\report_2\`.

**Phi** (mean of the ensemble; the base was 10000): prior 10788, iteration
1 6679, iteration 2 **6394** (−36 %).

| group | base | iteration 2 |
|---|---|---|
| hds_p0 | 3500 | **604** |
| sm_sm | 1750 | 1333 |
| q_log | 1500 | 1337 |
| q_year | 500 | 334 |
| sm_p0 | 1750 | 1747 |
| et_rs | 1000 | 1027 |
| hds_c (soft) | 6.5 | 13.5 |

**What the calibration moved** (ensemble means):
- K layer 1 ×0.19, with 31 % of the ensemble on the ×0.1 bound;
- Sy layer 2 ×2.1; UZF vks ×2.7; uzf_eps 3.54 (83 % on MF6's 3.5);
- DRN ×0.45 (27 % on ×0.1); streambed K 0.05 m/d; pond-bed leakance
  0.002 1/d; soil Ks ×1.2.

**The spread** fell by 60–87 % per group, and the phi spread to 27 (0.4 %
of the mean): the collapse of a 30-member ensemble against 215
parameters. The pilot shows the direction; the uncertainty needs the 100
of the production.

**10 agents were slower than 6.** Runs took 122 min against 55–65 with 6
(8 physical cores, the network drive): 0.08 runs a minute against 0.10.
Six to eight agents for the production.

### What the pilot shows

1. **P0's head cannot rise above the soil base.** The MF6 top is the base
   of the soil column (~738.2 m at P0) and the seepage drain removes any
   rise; the observed winter heads reach 739.6 (the land surface). The
   corrected head (`sim_hcorr_m`) stays equal to MF6's: MMsoil adds a
   horizon only when it is EXACTLY full while groundwater seeps up, and
   P0's lower horizon sits near, not at, porosity all winter (as observed:
   θ2 ≈ 0.19). The peaks above the soil base are unreachable, and fitting
   them may be what drives K layer 1 onto its bound.
2. **Summer low flows are too high:** the gauge drops to ~0 in the 2010
   and 2011 summers, while the model keeps ~0.05 mm/d of baseflow
   (GW→SFR).
3. **ET:** unchanged, as expected (10 %, validation-like).

### Proposed before the production run (for the user)

- **P0 above the soil base:** compare those heads with the water table
  including the saturated soil (a tolerance on "full" in the corrected
  head), or make them an inequality ("at least the soil base"; PEST++
  `g_` groups).
- **K layer 1:** keep the ×0.1 bound (0.005 m/d at the lowest) until
  item 1 is settled; it is likely compensating.
- **Production:** 100 realisations, 3 iterations, 6–8 agents — about 3
  days.

---

## 8. The production run (launched 2026-10-13 22:33)

The user's decisions after the pilot (2026-10-13):

1. **Heads against MMsoil's corrected head** (`sim_hcorr_m`, MF6's head plus
   the saturated soil above the soil base) at every point. It equals MF6's
   head where the soil is not saturated.
2. **No summer flow at La Mata.**
   - Where it came from, in a calibrated pilot run: the summer baseflow
     (~230 m³/d of gains) came mostly from reaches on layer 2, whose
     summer head stood 1.7–2.3 m above the bed under the low-K layer 1 the
     pilot chose to lift P0's head; plus layer-1 reaches near the outlet
     (0.1–0.3 m).
   - Stream evaporation was tried as the lever and dropped: with the
     calibrated parameters it slows the model to 8–120 s a day (against
     1.4), even capped at 0.5.
   - Instead, a new group **`q_dry`**: the model's flow against zero on
     the gauge's dry days (≤ 0.005 mm/d), in mm/d; 33 days in the window.
     Flow keeps 20 %: q_log 10, q_dry 5, q_year 5.
   - The corrected heads should also remove the pressure that pushed
     K layer 1 down.
3. **Bounds widened where the pilot hit them:** K ×0.01–30, Sy ×0.03–30,
   DRN ×0.01–100, field capacity ±0.2. `uzf_eps` stays ≥ 3.5, MF6's
   minimum.

**Run:** 100 realisations, 3 iterations, 6 agents (measured faster than
10). Expected about 3 days: 101 prior runs, then per iteration 12 trial
runs plus 96 updates, at ~1 h a run with six at once.
8955 targets, 5955 weighted.

The pilot's results are kept in `PEST\master_pilot`, the first Monte
Carlo's in `PEST\master_prior_mc`.

**Interrupted and resumed.** A server reboot (2026-10-14 ~18:52) stopped
the production run during its prior: 88 of 100 runs completed, 3 failed,
9 unfinished. Kept as `PEST\master_prod_prior`.

Resumed 2026-10-15 15:29 (user's go) with
`pest_lamata_run.py --restart-from master_prod_prior`. It reads the 88
completed runs back from the run storage and starts at iteration 1:
- 2 realisations dropped as bad (phi above the 2.5-sigma threshold), 86
  continue; the best starting phi is 6745;
- per iteration 24 trial runs (3 lambdas × 8) plus 78 updates;
- about 2.5–3 days at ~75 min a run with six at once.

**Stopped in iteration 1; the Tg curve fixed; warm start (2026-10-16).**
Iteration 1's runs took 130-460 min, not ~75. Runs past twice the
average were killed: 12 of the first 43, the slow ones. The cause was
MF6 stalling where summer water tables reached the 15 m root tips (the
Tg curve's ramp shrank to a step there; analysis doc §8.25). The user
stopped the run (kept as `PEST\master_iter1_stopped`).

- **Iteration 1's update kept.** Its 78 parameter sets (lambda 0.1) were
  read from the run storage into `PEST\iter1_warmstart\iter1_par.jcb`.
  The 34 runs that completed had a mean phi of 6882 (best 6052) with the
  old curve; the prior's mean was 11224 (best 6745).
- **Tg ramp anchored at the root tip** (99669a1). One realisation: 149 min
  instead of more than 8 h.
- **Setup again** with the fixed code (`org`, base run 40 min, `pst`). The
  measured ET at the points comes from the runs' exports (group `et_pt`,
  EC, weight 0). Overdue limit 4x the running mean. The previous setup is
  kept as `*_prod_20261013`.
- **Warm start, launched 17:55** (user's choice):
  `pest_lamata_run.py --workers 20 --noptmax 2 --start-ensemble
  iter1_warmstart\iter1_par.jcb`. PEST++ first runs the 78 sets with the
  fixed model, then makes 2 iterations: the 2nd and 3rd updates of the
  plan. 20 agents: the server has 32 cores since the 2026-10-15 reboot.

**The warm start stopped; a new calibration (2026-10-17/18).**

The warm start (iteration 1: mean phi 6899 → 6536, best 5840) was stopped
by the user, and kept as `PEST\master_warm_stopped`. Its two failures,
total ET above PET, led to a coupling artefact present in every run: UZF's
rejected infiltration flip-flopping day by day (analysis doc 8.26). The
calibration had been feeding it, by lowering layer-1 Sy and raising soil
Ks. Three changes before the new calibration:

1. **The soil percolates what UZF can take; the rest runs off** (4ad1d5b).
2. **A carried UZF state keeps its fill, not its content.** A state saved
   with θr/θs 0.44/0.45, restarting a realisation with other water
   contents, now carries (wc − θr0)/(θs0 − θr0) onto the new range,
   instead of a clipped content that started UZF full.
3. **Hydraulic parameters from the user's table**
   (`lamata_sardon_hydraulic_parameters.xlsx`). PEST only: the reference
   model keeps Sy 0.01, θs 0.45, K 0.05 m/d.

| Layer | Sy: start (bounds) | Sr: start (bounds) | n = θs = Sy + Sr | K [m/d]: start (bounds) |
|---|---|---|---|---|
| L1 saprolite | 2.5 % (0.5–5) | 27.5 % (4.5–30) | 30 % (5–35) | 0.05 (2.3e-5–1.12) |
| L2 fissured | 1 % (0.1–2) | 4 % (0.9–9) | 5 % (1–11) | 0.5 (2.3e-4–11.2) |

How each is set:
- **Sy:** pilot points as before, on the 0.01 rasters, with the table's
  start and bounds per layer.
- **Sr:** one value per layer (`sr_l<L>`). UZF's θs = Sy + Sr in every
  cell, so θr = Sr.
- **Widened maxima** (user): saprolite Sr 25 → 30 %, fissured Sy 1 → 2 %,
  so the table's medians sit inside their bounds.
- **K bounds:** they span every value gathered, from 2.3e-5 (the table's
  minimum) to 11.2 m/d (Legchenko et al. 2004). They are split so layer
  1's lie a decade below layer 2's (user). Layer 2 starts a decade above
  layer 1.
- **K order:** the model holds layer-1 K at most 0.1 × layer-2 K, cell by
  cell (`k_l1_max_ratio.txt`). The prior is drawn to respect it.

