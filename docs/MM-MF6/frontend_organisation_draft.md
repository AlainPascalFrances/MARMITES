# Front-end organisation — DRAFT for review

Draft 1, 2026-09-10. **This is for you to correct, not to approve.** Every
parameter list below was extracted mechanically from the ini files and the
code; every "live / superseded / dead" verdict comes from grepping the modules
that actually run. What I could not derive — grouping, which parameters a
modeller should see, sensible ranges — is what I need back.

Order follows your brief: **1 surface → 2 soil → 3 MF → 4 plotting**, each with
a master on/off switch.

---

## 1  Three findings that shape the design

### 1.1  The MF ini is not obsolete as a *file* — 34 of its parameters are live

It is read on every run (`run_lamata_mf6.py:178`, and the run prints
`parameter set: __inputMF_flopy_v3_2s1L.ini`). What is obsolete is the
**MF-NWT half of its content**: 72 of 106 named parameters are referenced
nowhere in the MF6 path.

The 34 that are live have **no home anywhere else** — `marmites_config.py`
holds run control (modes, spin-up, package switches, post-processing) and
contains no `hk`, `ss`, `sy`, `ibound`, `elev`, `thick`, `thts`, `thti`,
`eps`. So the aquifer parameterisation still comes only from that ini.

> **The goal of this remodelling, stated plainly:** the front-end is what lets
> all three ini files be retired. Migrate the live parameters into the
> configuration, drop the dead ones, delete the files.

### 1.2  `__inputMMsurf4MMsoil.txt` is a copy, and it currently wins

You are right that it is a reorganisation of what is already in the two ini
files. Worse than redundant: **it is authoritative today**. MMsoil reads `Zr`,
`kTg_min`, `kTg_max`, `kT_f`, `kT_s` from the handover, not from
`__inputMMsurf.ini` — so with MMsurf not running, editing those in the ini
changes nothing.

Proposal, per your instruction: **the handover disappears.** The front-end
holds one value per parameter; MMsurf and MMsoil both read it from the
configuration. Nothing is written to be read back.

### 1.3  The raster filenames in the MM ini are ignored — the driver hardcodes them

`gridMETEO_fn`, `gridSOIL_fn`, `gridSOILthick_fn` and the rest are in the ini,
but `setup_lamata` hardcodes `'inputMETEOzones.asc'`, `'inputSOILzones.asc'`,
`'inputSOILthick.asc'`, `'inputSTREAMhmax.asc'`, `'inputSTREAMw.asc'`,
`'inputIRRzones.asc'`. So those ini entries are inert. They should become
front-end fields that the code actually honours.

---

## 2  The four groups and their master switches

| # | Group | Switch | What OFF means |
|---|---|---|---|
| 1 | **Surface** (MMsurf) | `MARMsurf_yn` | MMsurf is not run. The *Time & forcing series* must already exist in the workspace, and the page says so with their dates and row counts. |
| 2 | **Soil** (MMsoil) | `MMsoil_yn` | No soil water balance. Only meaningful with MF off too, or as a forcing-only run. Legacy allowed `-1` = run MMsoil once, for calibration — keep? |
| 3 | **MODFLOW 6** | `MF_yn` | Write the MF6 input files and stop — today's `run.build_only`. |
| 4 | **Plotting** | `plot_yn` | No figures. Today's `postproc.enable` / `postproc.preproc`. |

Each group also carries its own sub-switches (per-package for MF, per-figure
family for plotting).

**What the code must do.** The switches cannot be decoration:

- `MARMsurf_yn = 1` → the driver calls `MARMITESsurf.startMARMITESsurface.MMsurf()`
  before the coupled loop, with the parameters from the configuration, and
  produces the forcing series into the **workspace** (not the repository).
- `MARMsurf_yn = 0` → it checks the forcing series exist and are long enough
  for the requested stress periods, and **fails clearly** if not, rather than
  running on a stale file.
- The state guard already built for the grid (WP0.6) is the pattern: a
  sidecar records which parameters produced the forcing, and a mismatch stops
  the run instead of silently reusing it.

---

## 3  Group 1 — SURFACE  (MMsurf)  ·  switch `MARMsurf_yn`

Source: `MMsurf_ws/__inputMMsurf.ini` (100 parameters) + its time-series
inputs. Labels are the ini's own names; the explanation column is the ini's
own comment, lightly rephrased.

### 1.1  Meteo station — 8 fields × NMETEO (currently 1)

| field | units | meaning |
|---|---|---|
| `phi` | deg | latitude of the station (> 0 = northern hemisphere) |
| `Lm` | deg | longitude, west of Greenwich |
| `Z` | m | altitude above sea level |
| `Lz` | deg | longitude of the centre of the local time zone |
| `FC` | h | time shift if the site is not in its nominal time zone |
| `DTS` | h | data time shift, for data not acquired at standard clock time |
| `z_m` | m | height of the wind-speed measurement |
| `z_h` | m | height of the humidity measurement |

### 1.2  Vegetation — 19 fields × NVEG (currently 3: grassMU, Qilex, Qpyr)

`h_d`, `h_w` (m, plant height dry/wet) · `S_w` (mm, canopy storage) ·
`C_leaf_star` (m/s, max leaf conductance) · `LAI_d`, `LAI_w` (m²/m²) ·
`f_s_vd`, `f_s_vw` (shelter factor) · `alfa_vd`, `alfa_vw` (albedo) ·
`J_vd`, `J_vw` (julian day, season start) · `TRANS_vdw`, `TRANS_vwd` (d,
transition length) · **`Zr`** (m, max root depth) · `kTg_min`, `kTg_max`,
`kT_f`, `kT_s` (transpiration sourcing factors).

> The last five are the ones duplicated into the handover today (§1.2), and
> `Zr` is also the natural source for the WP2 UZF extinction depth.

### 1.3  Crops — 11 fields × NCRP (currently 1), plus NFIELD

`h_c`, `S_w_c`, `C_leaf_star_c`, `LAI_c`, `f_s_c`, `alfa_c`, `Zr_c`,
`kTg_min_c`, `kTg_max_c`, `kT_f_c`, `kT_s_c` — the crop counterparts of the
vegetation set.

### 1.4  Surface soil — 8 fields × NSOIL (currently 3: alluvium, regolith, outcrop)

`por`, `fc` (m³/m³, surface 1 cm porosity and field capacity) · `alfa_sd`,
`alfa_sw` (albedo) · `J_sd`, `J_sw` (julian day) · `TRANS_sdw`, `TRANS_swd` (d).

> Note these are the **surface** soil properties MMsurf needs for evaporation.
> They are NOT the soil-column properties of group 2 — different file,
> different meaning, same word. Worth distinct labels in the UI.

### 1.5  Time-series inputs to MMsurf

| file | what it is |
|---|---|
| `__meteoTB.txt` | the meteorological record MMsurf converts into PET and rainfall |
| `__IRR_TS.txt` | irrigation time series per zone (only if `irr_yn = 1`) |
| `__inputFIELD1_crop_schedule.txt` | crop schedule per field |

### 1.6  *Time & forcing series* — MMsurf's OUTPUT

These are what your "Time & forcing series" group contains. They are produced
when `MARMsurf_yn = 1` and consumed as-is when it is 0:

`inputDATE.txt` · `inputZONRF_veg_d.txt` (rainfall) ·
`inputZONTF_veg_d.txt` (throughfall) · `inputZONPT_veg_d.txt` (potential
transpiration) · `inputZONLAI_veg_d.txt` · `inputZONPE_d.txt` (potential
evaporation) · `inputZONEo_d.txt` (open-water evaporation) — plus the
irrigation trio `inputZONRF_irr_d`, `inputZONTF_irr_d`, `inputZONPT_irr_d`
and `inputZONcrop_irr_d`.

**Open question A:** these currently live in `example/LaMata/`, i.e. in the
repository. If MMsurf generates them they become run OUTPUT and belong in the
workspace. Do you want them moved, or kept committed as a reproducible
baseline?

---

## 4  Group 2 — SOIL  (MMsoil)  ·  switch `MMsoil_yn`

### 2.1  Soil column parameters — `MF_ws/inputSOILparam.txt`

Per soil zone (3), per soil layer (`nsl` = 2, 2, 1):
`st` (soil type, text) · `slprop` (layer thickness proportion) ·
`Smax` (max soil moisture) · `Sfc` (field capacity) · `Sr` (residual) ·
`Si` (initial) · `Ks` (m/d, saturated conductivity).

### 2.2  Zone and property rasters (the filenames the driver currently hardcodes)

| raster | what it maps |
|---|---|
| `inputMETEOzones.asc` | meteo zone per cell |
| `inputSOILzones.asc` | soil zone per cell |
| `inputSOILthick.asc` | soil column thickness |
| `inputIRRzones.asc` | irrigation zone (if `irr_yn`) |
| `inputSTREAMhmax.asc` | max stream/surface water height |
| `inputSTREAMw.asc` | channel width |

### 2.3  Vegetation area fractions

`inputZONVEGarea_*` — the per-cell share of each vegetation type.

### 2.4  Observations

`inputObs.txt` (points: name, x, y, layer) · `inputObsHEADS_*.txt` ·
`inputObsSM_*.txt` · `inputObsRo_*.txt`, plus `rmseHEADSmax`, `rmseSMmax`.

**Open question B:** `rmseHEADSmax` / `rmseSMmax` were the legacy convergence
criteria for the MM↔MF Picard loop, which no longer exists. Do they still
mean anything, or do they go to §7?

---

## 5  Group 3 — MODFLOW 6  ·  switch `MF_yn`

The 34 live parameters of the MF ini, regrouped. **Everything else in that
file goes to §7.**

### 3.1  Grid and discretisation
`nlay`, `nrow`, `ncol`, `delr`, `delc`, `laycbd`, `lenuni` — plus the WP1c
`[grid]` block (`kind`, `resample`, the `[grid.voronoi]` sizes) which already
supersedes `reggrid`.

### 3.2  Geometry rasters
`elev_fn` (`elev_sinkfil.asc`) · `thick_fn` (`thick_l1/l2.asc`) ·
`ibound_fn` (`ibound_l1/l2.asc`) · `strt_fn` (`hi_topL1.asc`) · `hnoflo`.

### 3.3  Aquifer properties (per layer)
`hk_fn` (`hk_l1/l2.asc`) · `vka_fn` · `ss_fn` (`Ss_l1/l2.asc`) ·
`sy_fn` (`Sy_l1/l2.asc`) · `laytyp`, `layavg`, `layvka`, `laywet`.

### 3.4  UZF
`uzf_yn` · `iuzfopt` · `ntrail2` · `nsets` · `surfdep` · `vks` · `eps` ·
`thts` · `thti` — plus WP0's `[uzf] vks_scale` and WP2's `[et]` block.

### 3.5  Boundary packages
`wel_yn` · `drn_yn` + `drn_elev`, `drn_cond` · `ghb_yn` + `ghb_head`,
`ghb_cond` — plus the WP0 `[seep]` block (`kind`, `cond`) and the
`[sfr]`/`[lak]`/`[crr]` blocks.

### 3.6  Water-balance aggregation
`Mnlay`, `Mlay` (hydrogeological layers from MODFLOW layers) · `h_plt`,
`h_lbl`.

### 3.7  Already in TOML, keep there
`[run]` mode/relax/nsp/daily/ats · `[spinup]` · `[paths]`.

---

## 6  Group 4 — PLOTTING  ·  switch `plot_yn`

Today: `postproc.enable`, `postproc.preproc`, `postproc.only`,
`postproc.map_days`, `postproc.sankey_*`. From the MM ini, the ones that
describe *what to draw* rather than *how the model runs*:

`plt_out_obs` (per-point series) · `WBsankey_yn` (Sankey) ·
`plt_WB_unit` (`year` | `day`) · `iniMonthHydroYear` (hydrological year start
— **this one is live**, `marmites_postprocess.py` reads it) · `plt_input`
(input maps) · `MMsurf_plot`.

Sub-switches worth having, matching the figure families that exist:
input maps · result maps · Sankey (catchment + per point) · obs time series ·
water-budget figures · NWT comparison.

---

## 7  Superseded and dead — the list you asked for

Your intuition was right for the MM ini and **wrong for the MF ini**: most of
what goes there is not plotting, it is MF-NWT solver machinery.

### 7.1  MM ini — superseded by the TOML (plotting/run control)
`run_name` → `meta.name` · `plt_out`, `plt_freq`, `nrangeMM`, `nrangeMF`,
`ctrsMM`, `ctrsMF`, `ntick`, `animation`, `animation_freq`,
`maxYearsTickTrimester`, `maxYearsTickSemester` → `[postproc]` ·
`chunks` (HDF5 compression) → gone, the coupled writer does not use it.

### 7.2  MM ini — dead, the Picard loop they belonged to was removed in Phase 1
`MF_yn`* , `MF_lastrun`, `convcrit`, `convcritmax`, `ccnum`.
(*`MF_yn` returns as your group-3 master switch, but with a new meaning:
build-only versus run.)

### 7.3  MM ini — inert because the driver hardcodes the value (§1.3)
`gridMETEO_fn`, `gridSOIL_fn`, `gridSOILthick_fn`, `gridIRR_fn`,
`inputFile_PAR_fn`, `inputFile_TS_fn`, `inputFile_TSirr_fn`,
`outputFILE_fn`, `MMsurf_fn`. **These should come back as real fields.**

### 7.4  MF ini — dead, MF-NWT only (57 parameters)
Solver blocks: `ext_pcg`, `hclose`, `rclose`, `ext_nwt`, `HEADTOL`,
`FLUXTOL`, `MAXITEROUT`, `THICKFACT`, `LINMETH`, `IPRNWT`, `IBOTAV`,
`OPTIONS` — MF6 uses IMS, configured in `marmites_mf6.py`.
Flow packages: `ext_lpf`, `ilpfcb`, `nplpf`, `ext_upw`, `iupwcb`, `npupw`,
`iphdry`, `chani`, `storagecoefficient`, `constantcv`, `thickstrt`,
`nocvcorrection`, `novfc` — MF6 uses NPF/STO.
Output control: `ext_oc`, `ihedfm`, `iddnfm`, `ext_cbc`, `ext_heads`,
`ext_ddn`, `MFout_yn`.
File extensions: `namefile_ext`, `ext_dis`, `ext_bas`, `ext_uzf`, `ext_wel`,
`ext_ghb`, `ext_drn`, `ext_rch`.
UZF1-only: `nuztop`, `irunflg`, `ietflg`, `NUZF2`, `NUZF3`, `NUZF4`,
`EXTDP`, `EXTWC`, `iuzfcb1`, `iuzfcb2`, `nuzgag`, `iuzrow`, `iuzcol`,
`iftunit`, `iuzopt`.
Recharge package: `rch_yn`, `ext_rch`, `nrchop`, `rch_user` — MM supplies
recharge through UZF; RCH is never used.
Also `exe_name`, `version` (MF6 comes from `mm_paths.MF6_EXE`/`LIBMF6`),
`itmuni`, `reggrid`, `dum_sssp1` (now `[spinup]`), `finf_user`, `wel_user`.

---

## 8  Questions for you

1. **Open question A** — do the *Time & forcing series* move to the workspace
   once MMsurf generates them, or stay committed as a baseline?
2. **Open question B** — do `rmseHEADSmax` / `rmseSMmax` still mean anything?
3. `MMsoil_yn = -1` (run MMsoil once, for calibration) — keep, or is that now
   PEST's job?
4. Should the four master switches live in one `[run]` block
   (`run.surface`, `run.soil`, `run.mf`, `run.plot`) or in each group's own
   block? I lean to one block, so a run's shape is visible in one place.
5. Anything in §7 you want kept that I have proposed to delete.
