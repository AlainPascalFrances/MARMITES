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

### 1.4 Rastering or griding from input vector info
We will change the paradigm in order to ease the life of the users. The information will be provided as scalars or shape file (points, lines or polygons, with or without columns with parametric values. This vectorial information will be later converted into model input by wrapping it or projecting it on the grid (structured, voronoi, quadtree). The grid itself will be defined in a group, indicated the resolution, refinemnet, that will be specific to one of the 3 grid type and will be defined inside the catchment boundary, that is one of the main input that the front-end must require (polygon that will be in the main CRS of the project, must be projected, metric units) 
---

## 2  The four groups and their master switches

| # | Group | Switch | What OFF means |
|---|---|---|---|
| 1 | **Surface** (MMsurf) | `MARMsurf_yn` | MMsurf is not run. The *Time & forcing series* must already exist in the workspace, and the page says so with their dates and row counts. If the *Time & forcing series* does not exist and MMsurf is off, a warning message must appear and code stops.|
| 2 | **Soil** (MMsoil) | `MMsoil_yn` | No soil water balance. This option on/off will be removed since MMsoil and MF are run together now. So the legacy option that allowed `-1` = run MMsoil once, for calibration, will be also removed. |
| 3 | **MODFLOW 6** | `MF_yn` | Write the MF6 input files and stop — today's `run.build_only`. Option to be removed also, the MF run is compulsory, together with MMsoil, as stated in the previous line.|
| 4 | **Plotting** | `plot_yn` | No figures. Today's `postproc.enable` / `postproc.preproc`. |

So only groups 1 and 4 carries a sub-switches (per-figure family for plotting).

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

### 1.6  Zone and property rasters (the filenames the driver currently hardcodes)

| raster | what it maps |
|---|---|
| `inputMETEOzones.asc` | meteo zone per cell |
| `inputIRRzones.asc` | irrigation zone (if `irr_yn`) |
The inputMETEOzones.asc will be eliminated and produced using the grid characteristics as Thiessen polygon (derived from the grid defined in 1.4). If there is only one station, the meteo zone will be the same for the whole catchment. 

### 1.7  *Time & forcing series* — MMsurf's OUTPUT (to be eliminated from the front-end)

These are what your "Time & forcing series" group contains. They are produced
when `MARMsurf_yn = 1` and consumed as-is when it is 0. So they will not appear in the front-end.

`inputDATE.txt` · `inputZONRF_veg_d.txt` (rainfall) ·
`inputZONTF_veg_d.txt` (throughfall) · `inputZONPT_veg_d.txt` (potential
transpiration) · `inputZONLAI_veg_d.txt` · `inputZONPE_d.txt` (potential
evaporation) · `inputZONEo_d.txt` (open-water evaporation) — plus the
irrigation trio `inputZONRF_irr_d`, `inputZONTF_irr_d`, `inputZONPT_irr_d`
and `inputZONcrop_irr_d`.

These currently live in `example/LaMata/`, i.e. in the
repository. As MMsurf generates them, they are a run OUTPUT and belong in the
workspace. So they will be deleted from the repo. They do not appear in the front-end.

IMPORTANT NOTE:
In groups 1.1 to 1.4, the list of parameters are repeated as a function of NMETEO, NVEG, NCRP and NSOIL. So the front-end should put one column for each entry.
In relation to 1.5, verify in the code how should be the format please and make proposal how the data must be organized (is it one txt file with several columns, or the data are repeated sequencially? Instructions about the data format should be inserted in the front-end  
---

## 4  Group 2 — SOIL  (MMsoil)

### 2.1  Soil column parameters — `MF_ws/inputSOILparam.txt`

Per soil zone (3), per soil layer (`nsl` = 2, 2, 1):
`st` (soil type, text) · `slprop` (layer thickness proportion) ·
`Smax` (max soil moisture) · `Sfc` (field capacity) · `Sr` (residual) ·
`Si` (initial) · `Ks` (m/d, saturated conductivity).

### 2.2  Zone and property rasters (the filenames the driver currently hardcodes)

| raster | what it maps |
|---|---|
The following raster will be elimintaed. Instead, a vector, polygon layer is required, thta must contain the soil zone code that must correpond with the soil zone code of MF_ws/inputSOILparam.txt. The soil polygon layer may have a column called thick_m, in which the thickness will be uniform for each polygon of the layer. Note that if the inputSOILthick.asc exists, it will be the default value. If none of them exist raise an error.
| `inputSOILzones.asc` | soil zone per cell |
The following raster can exist and will be projected on the grid defined in 1.4. It must be coincident with the catchment polygon.
| `inputSOILthick.asc` | soil column thickness |
I think that the two following rasters must be eliminated. They will be produced using the shape file of hydrography that will have a column with spatial variation of these parameters (a value for each segment) or they will be fixed for the whole catchment. To produce them from the shape file, the script should use the hydrography layer (currently: E:\00code_ws\LAMATA_new\GIS\hydrography.shp)   
| `inputSTREAMhmax.asc` | max stream/surface water height |
| `inputSTREAMw.asc` | channel width |

### 2.3  Vegetation area fractions

`inputZONVEGarea_*` — the per-cell share of each vegetation type.

### 2.4  Observations

`inputObs.txt` (points: name, x, y, layer) · `inputObsHEADS_*.txt` ·
`inputObsSM_*.txt` · `inputObsRo_*.txt`
NOTE: `rmseHEADSmax` and `rmseSMmax` eliminated

---

## 5  Group 3 — MODFLOW 6

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
`iuzfopt` · `ntrail2` · `nsets` · `surfdep` · `vks` · `eps` ·
`thts` · `thti` — plus WP0's `[uzf] vks_scale` and WP2's `[et]` block.
NOTE: `uzf_yn` removed since it is compulsory

### 3.5  Boundary packages
`drn_yn` + `drn_elev`, `drn_cond` · `ghb_yn` + `ghb_head`,
`ghb_cond` — plus the WP0 `[seep]` block (`kind`, `cond`) and the
`[sfr]`/`[lak]`/`[crr]` blocks.
NOTE: 'wel_yn' removed since it is compulsory to compute ETG

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
YES, definitively. They are also going out from the front-end, as stated in point 1.7
 
2. **Open question B** — do `rmseHEADSmax` / `rmseSMmax` still mean anything?
NO, remove

3. `MMsoil_yn = -1` (run MMsoil once, for calibration) — keep, or is that now
   PEST's job?
NO, remove

4. Should the four master switches live in one `[run]` block
   (`run.surface`, `run.soil`, `run.mf`, `run.plot`) or in each group's own
   block? I lean to one block, so a run's shape is visible in one place.
run.surface in one block, run.soil and run.mf togerher in one block, run.plot in one block

5. Anything in §7 you want kept that I have proposed to delete.
NO, it is ok. Just verify thta the time series plots must be started at the used-defined starting hydrologic year.  
