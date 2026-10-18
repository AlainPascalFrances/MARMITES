# PEST++ in the front-end: plan

Proposed 2026-10-17, approved by the user the same day with the changes
noted below. Built in phases; each phase has its own tests and commit, and
the app is restarted through `code/tools/launch_app.py` after it.

## 1. Where it goes

Two pages after Run, before Results:

| Page | Purpose |
|---|---|
| **8 — PEST++** | Everything that defines a calibration, in tabs. Its choices are config fields, edited and saved like the other panels' |
| **9 — PEST++ progress** | Live monitor plus results. Refreshes itself; can stay open for days |
| **10 — Results** | The model run results (was page 8) |

**User decisions.**
- PEST before Results.
- The PEST settings live in the case config's `[pest]` section.

As everywhere else in the app (CLAUDE.md rule 2), every PEST choice is a
config field: help text in `app/lib/schema.py`, validation in
`marmites_config.py`, read by the PEST tools, covered by the wiring test.
The PEST tools read `cfg.pest` instead of their constants (calibration
start, pilot-point spacing, group shares, noise, bounds...). Today's values
become the defaults, so a setup made from the panel equals today's.

## 2. Page 8 — PEST++ (tabs)

**① Base and workspace**
- The PEST workspace stays on Home (user). Page 8 shows its path, read-only.
- The **run solution** to calibrate from: a finished run from the run
  registry. Its configuration becomes the model PEST runs; its saved end
  state is the start state of every realisation.
- Today a run only points at the live `lamata.toml`. So every launch
  also saves a copy of its effective configuration (with overrides) in
  `runs/<run_id>/`. For an older run the page falls back to the current
  configuration plus a saved state chosen by hand, with a warning.
- The calibration window: the first fitted day; the hydrological year
  before it is the warm-up.

**② Parameters**, one row per family: on/off, lower, upper, start, log.

| Family | Kind |
|---|---|
| K, Sy, Ss per layer | pilot points (spacing, variogram range); Ss with K and Sy (user, 2026-10-17; it was one value per layer) |
| DRN conductance per layer | one value per layer |
| Soil field capacity per horizon | scalar, **bounded by the wilting point (residual) and the porosity** (user): the range can never leave that interval |
| Soil Ks per horizon | scalar |
| SFR bed K, LAK bed leakance, UZF eps and vks, CRR beta | scalar |

The prior's spread is shown before launching, with the pilot points on the
mesh.

**③ Observations and weights**
- The groups and their share of phi (now heads 35 / soil moisture 35 /
  flow 20 / ET 10).
- The soft-data fraction (C wells), the noise per group, the dry-day
  threshold.
- The points used (C4, C5, W1 switched on).
- The ET sources: MODIS from the Home RS folder; the point ET (EC,
  lysimeters) as validation.
- A summary of the targets per group in the window.

**④ Ensemble and PEST++ options**
- Realisations, iterations (`noptmax`), lambda multipliers, bad-phi
  sigma, overdue limit.
- Agents, up to the cap (20).
- The start: the prior, a **warm start** from a previous run's iteration,
  or a **restart** of an interrupted run.

**⑤ Set up and launch** (user: launching from the app, yes)
- **Set up**: the template from the chosen base, with its base run and
  that run's log.
- **Launch**, with a confirmation: hours to days on many cores.
- **Add agents**, **Stop at the end of this iteration**, **Stop now**.
- Detached processes through the run registry (`app/lib/runs.py`), never
  inside Streamlit, so closing the browser stops nothing. One PEST job at a
  time.

## 3. Page 9 — PEST++ progress (refreshes every ~30 s)

**Status.** The stage: setting up, base run, initial ensemble, iteration k
(lambda test, then the remaining realisations), finished, stopped. Elapsed
time and an estimated end.

**Per iteration.** Runs planned, completed, running, failed, killed as
overdue; mean and longest run time.

**Agents.** Each agent's current realisation and elapsed time.

**Phi.** Per iteration (every realisation, the mean and the best), per
group, and the lambda tests.

**Failures.** Realisation and cause, read from the agent's model log (e.g.
"total ET above PET").

**Parameters.** Spread per iteration; the share on a bound.

**Results** at each iteration's end and at the run's end:
- the iteration report (`pest_lamata_report.py`): ensemble envelopes
  against the observations per group;
- the posterior against the prior; K, Sy, Ss multiplier maps on the mesh;
- the EC ET validation.

Then **"Write the calibrated parameters as a model configuration"**: the
best or the median realisation as a configuration plus its multiplier set,
named after the input model with **`_PEST`** appended (user), e.g.
`lamata.toml` -> `lamata_PEST.toml`. The Run page runs it with full
figures.

## 4. Phases

| # | Content |
|---|---|
| 1 | Progress page on the current `PEST\master` (no config change; Results renumbered 10) |
| 2 | `[pest]` config fields; the PEST tools read them; page 8 tabs ① to ④ |
| 3 | A configuration copy saved with every run; the run solution as the base |
| 4 | Set up, launch, add agents and stop from the app (tab ⑤) |
| 5 | The PEST results, and the calibrated parameters as `<model>_PEST` |

## 5. Phase 2, as built (2026-10-18)

**The `[pest]` section** (`marmites_config.Pest`, validated in
`RunConfig._pest_problems`, help texts in `app/lib/schema.py`):

| Tab | Fields |
|---|---|
| ① Base | `base_state` (a saved state of the run workspace), `calib_start` (blank: the first hydrological year after a whole year of record, 2009-10-01 on La Mata) |
| ② Parameters | `pp_space`, `vario_range`, `k_order_decades`; the table `[[pest.par]]`: family, layer (or horizon, or 0), on, start, lower, upper, log |
| ③ Observations | the table `[[pest.group]]`: share of phi, noise, relative noise (every group once); `soft_fraction`, `points_on`, `q_dry`, `q_eps` |
| ④ Ensemble | `nreal`, `noptmax`, `lambda_mults`, `bad_phi_sigma`, `overdue_giveup_fac`, `agents`, `start` (prior, warm, restart), `start_from` |

- **Units by family** (`PEST_FAMILIES`): K in m/d, Sy and Sr as water
  contents, the field capacity as a moisture, Ks in mm/d; Ss, the drains
  and UZF's vks as multipliers. The set-up turns K and Sy into multipliers
  of the reference rasters' median (uniform on La Mata, so exact).
- **An empty `par` means the defaults** for the configuration
  (`pest_lamata_setup.default_pars`: the calibration of 2026-10-18). The
  page shows them; saving an edit writes the table. A reset button puts
  the defaults back.
- **The field capacity cannot leave its interval** (user): a row's range
  must lie strictly between the horizon's wilting point `sr` and porosity
  `smax`, and PEST++ calibrates its place in that interval.
- **Ss on pilot points** with K and Sy (user); it was one value per layer.
- Also refused at save: a start outside its range, a log range from 0, a
  layer or horizon that does not exist, Sy + Sr up to 1, a K order no pair
  can hold, shares that do not sum to 1, a soft group with a share, a
  warm or restart start without `start_from`.
- **Retired:** `pest.use_pest_params`, `pest.params_npz` (read by nothing).

**The tools read it.** `pest_lamata_setup.py --config <file> [--step]`
(no more `--state`, `--nreal`, `--calib-start`); `pest_lamata_run.py
--config <file>` reads the launch settings at launch (no more `--workers`,
`--noptmax`, `--nreal`, `--start-ensemble`, `--restart-from`); `--join N`
stays. The set-up writes the resolved parameter table into its copy of the
configuration, so every later step reads the same rows.

**Not in phase 2:** the prior's spread with the pilot points on the mesh
(drawn at set-up, to show with phase 4), the run solution as the base
(phase 3), set-up and launch from the page (phase 4).
