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

## 2A  Panel 1 — THE GRID   *(added 2026-09-12, for review)*

> **Panel numbering.** This draft was written when the order was
> *1 surface → 2 soil → 3 MF → 4 plotting*. The panels are now
> **0 Overview · 1 Grid · 2 Surface · 3 Soil + MF · 4 Plots**, so group 1
> below is panel 2, groups 2 and 3 are panel 3, group 4 is panel 4. This
> section is the one that was missing: panel 1.

Panel 1 is answered **before every other panel** because everything else is
wrapped onto what it produces. The soil zones, the vegetation cover, the
stream network and the observation points are vector layers in the project
CRS; they are projected onto this grid at run time, so changing the grid does
**not** mean re-making any of them.

### 2A.1  Permanent block — asked once, common to every grid kind

Always visible, at the top, above the kind selector.

| field | TOML key | units | default | what it does |
|---|---|---|---|---|
| Catchment polygon | `grid.boundary` | — | `lm_lim.shp` | A **projected, metric** polygon in `DATA_ROOT/GIS`. Defines the active domain, the mesh boundary and the model rectangle. The shapefile stays in the GIS folder; it never enters the repository. |
| Project CRS | `grid.crs_epsg` | EPSG | `23029` | The CRS every other layer is assumed to be in, and the one the meteo station latitude/longitude is converted **into**. `0` = take it from the layer .prj of the Catchment polygon. |
| Grid kind | `grid.kind` | — | `voronoi` | **Combo box**: `structured` · `disv` · `voronoi` · `quadtree`. Chooses which sub-panel below opens. Voronoi is the default: in the schema, and in `lamata.toml` (D1 answer, (c)). |
| Wrapping rule | `grid.resample` | — | `auto` | How a value is carried onto a cell: `auto` = area-weighted mean for continuous fields, largest-overlap class for zone layers; `centre` = the value under the cell centre. Applies to every kind, so it belongs here. |

**Read-outs, not parameters** (computed from the polygon, shown beside the
fields so a wrong file is caught immediately): polygon area [km²], extent
[m × m], the CRS the .prj actually declares, and the number of cells the
current settings would give. An extent below 100 m is refused: that is
degrees, not metres.

### 2A.2  Sub-panel — `structured`

The rectangle, divided. The legacy DIS grid: kept permanently as the regression anchor and the MODFLOW-NWT comparison, no longer the starting point.

| field | TOML key | units | default | what it does |
|---|---|---|---|---|
| Cell size | `grid.cell_size` | m | `50` | Side of the square cell. |
| Buffer | `grid.buffer` | m | `0` | Extend the rectangle beyond the polygon bounding box. |
| Reproduce an existing grid | `grid.override.enable` | — | `false` | Escape hatch: take the origin and shape from the four fields below instead of deriving them from the polygon. |
| Origin X | `grid.override.xllcorner` | m | `739300` | Lower-left corner, only when the override is on. |
| Origin Y | `grid.override.yllcorner` | m | `4553050` | idem |
| Rows | `grid.override.nrow` | count | `65` | idem |
| Columns | `grid.override.ncol` | count | `60` | idem |

Derived and shown, not typed: rows, columns and origin when the override is
**off** — they come from the polygon, snapped **down** to a multiple of the
cell size. On `lm_lim.shp` at 50 m that is **63 × 60 from 739250, 4553100**,
against the legacy **65 × 60 from 739300, 4553050**.

### 2A.3  Sub-panel — `disv`

The same geometry, re-expressed as polygons. Nothing about the cells
changes, which is exactly what makes it the control for the unstructured
path. **Same fields as `structured`** — cell size, buffer, override.

### 2A.4  Sub-panel — `voronoi`

| field | TOML key | units | default | what it does |
|---|---|---|---|---|
| Background cell size | `grid.voronoi.cell_far` | m | `100` | Side of the **equivalent square**, so 100 aims at 10 000 m². |
| Refine along the streams | `grid.voronoi.stream_refine` | — | `true` | `false` = the SFRmaker approach, network mapped onto the background. True is the default.|
| Cell size near the stream | `grid.voronoi.cell_near_stream` | m | `40` | Ignored (as the other refinement options below) when the refinement is off, so the box is blocked while the option is not activated. |
| Stream corridor width | `grid.voronoi.stream_buffer` | m | `60` | Half-width of the refined corridor. |
| Maximum size ratio between bands | `grid.voronoi.grade_ratio` | — | `1.5` | How fast the cells may grow away from the corridor, and so how many bands there are. Held in (1, 3] — see D2. |
| Transition bands **(derived)** | `grid.voronoi.trans_levels` | m | `[15, 30, 45, 60]` | **Read-only.** The buffer DISTANCES from the centreline at which the size steps up, recomputed on every save from the corridor and the ratio. |
| Seed a cell per pond | `grid.voronoi.seed_ponds` | — | `true` | One pond-scale cell per pond, centroid-seeded (the CdL lesson: a cell can never follow a constraint polygon — seed a generator instead). |
| Buffer | `grid.buffer` | m | `0` | Still applies: it sets the domain rectangle that is triangulated. `grid.cell_size` does **not** appear here (D1 of 2A.7): the rectangle is snapped to `cell_far`. |

The build **prints the mean cell area it actually achieved** — the triangle
maximum-area conversion is calibrated, not exact. Read that line, not the
setting.

### 2A.5  Sub-panel — `quadtree`

| field | TOML key | units | default | what it does |
|---|---|---|---|---|
| Background cell size | `grid.cell_size` | m | `50` | The GRIDGEN base grid. |
| Buffer | `grid.buffer` | m | `0` | As above. |
| Refine along the streams | `grid.quadtree.refine_streams` | — | `true` | Off = an unrefined quadtree, which is geometrically the base grid. The producer says so rather than reporting success. |
| Refinement levels | `grid.quadtree.refine_level` | count | `2` | GRIDGEN halves a cell per level, so 2 on 50 m gives **12.5 m** along the streams. |

Needs the `gridgen` executable, and `pyshp` for the refinement features.
Both are environment facts, so the producer reports them, not this panel.

### 2A.6  The two buttons

**Create grid** — an *experiment*. Builds the mesh from the settings
**currently on screen**, reports cells / mean area / equivalent side /
signature, draws it, and keeps it in a list of attempts so two kinds or two
cell sizes can be compared side by side. It writes nothing the model reads,
and it can be pressed as often as wanted.

**Select this grid for the model** — the *commitment*. Writes the settings
that produced the selected attempt into the TOML, promotes its cached mesh
to the one a run will pick up, and stamps the configuration hash. From then
on the panel states which grid is selected and how it differs, if at all,
from the attempt on screen.

### 2A.7  Doubts, before I build this

1. **`cell_size` means two different things.** For `structured`, `disv` and
   `quadtree` it is the cell. For `voronoi` it is **not used for the cells at
   all** — `cell_far` is — it only snaps the bounding rectangle. Proposal:
   drop `cell_size` from the voronoi sub-panel entirely and snap the
   rectangle to `cell_far`. One number, one meaning.
2. **`grid.rebuild`** is a stored parameter today (re-mesh, ignore the
   cache). With *Create grid* as an explicit action it has no reason to be a
   saved setting. Proposal: delete it from the TOML; the button is the
   rebuild.
3. **The override block** is offered under `structured` and `disv` above, but
   `model_rectangle` honours it for `quadtree` too. Show it there as well, or
   restrict it to the two kinds the legacy comparison actually needs?
4. **Where does selected live?** Either `grid.*` in the TOML *is* the
   selection (simple, one source of truth, but then an experiment must be
   saved before it can be re-opened), or a small `[grid.selected]` sidecar
   records the signature and the panel warns when `grid.*` has drifted from
   it. I lean to the first.

### 2A.8  Second round of doubts (2026-09-12, after your edits)

Doubts 1, 2 and 4 of §2A.7 are settled as proposed; doubt 3 is settled as
**the override appears under `structured` and `disv` only**. Four things your
edits raise.

**D1 — which kind is the default? The section now says both.**
The table row says `structured`, your note on the same row says *Default is
Voronoi*, and §2A.2 still opens with "This is the default". In the code,
`marmites_config.Grid.kind` defaults to `structured`, with the comment that
it stays so until the WP1c.8 validation ladder clears the Voronoi mesh, and
`lamata.toml` sets it explicitly anyway. Three different things could be
meant:

  a. the **panel pre-selects** voronoi for a new project, the dataclass
     default stays `structured` (nothing about existing runs changes);
  b. the **dataclass default** becomes voronoi as well — it only affects a
     configuration that omits the key, so not La Mata;
  c. **`lamata.toml` itself** switches to voronoi, which changes the grid
     every committed comparison was run on.

I read your note as (a) + (b). Confirm, and say explicitly if you also want
(c).
ANSWER: I want (c).

**D2 — `trans_levels` are DISTANCES from the stream, not cell sizes.**
My label was wrong, and so is the one in `schema.py`. What the producer
actually does (`marmites_meshes._add_stream_regions`): it takes
`trans_levels ∪ {stream_buffer}` as **buffer radii in metres**, sorts them,
and gives the innermost band `cell_near_stream`, the outermost `cell_far`,
linearly in between. Two consequences:

  - With today's defaults the list is `[10, 20, 40, 70]` and
    `stream_buffer = 60`, so the bands are 10, 20, 40, 60, **70** m. The
    outermost is 70, which is *beyond* the declared corridor — the corridor
    really ends at 70 m and `stream_buffer` is just one more band. That alone
    justifies deriving the list rather than typing it.
  - To recompute it I need a rule. **Proposal**, which keeps the producer's
    existing linear size law and bounds how fast the cells may grow: let *r*
    be the largest acceptable size ratio between neighbouring bands
    (default **1.5**, the usual mesh-grading rule of thumb). Then

        n = ceil( (cell_far - cell_near_stream) / (cell_near_stream * (r - 1)) )
        band k is at distance  stream_buffer * k / n,   k = 1 .. n

    On your defaults (40 → 100 m, corridor 60 m, r = 1.5) that is **n = 3**:
    bands at 20, 40, 60 m carrying 40, 60, 80, 100 m cells. Accept r = 1.5 as
    the hidden constant, or do you want *r* exposed as a field, or the band
    count *n* fixed by hand instead?

  - And once it is derived: does `trans_levels` leave the TOML entirely (like
    `grid.rebuild`, doubt 2) with the bands printed in the build log, or stay
    as a read-only echo of what was computed? I lean to removing it — but if
    you ever want to hand-tune a band, it has to stay.
    
ANSWER: Keep r=1.5 as default constant, but can be changed by user, but with bounds (do you think that  >1 and < 5 is ok?). trans_levels stay as read-only of what is computed in the background.

**D3 — the override is hidden for voronoi and quadtree, but the code still
honours it.** `model_rectangle` applies `[grid.override]` whatever the kind.
So enabling it under `structured`, then switching to `voronoi`, silently
overrides the domain from a box that is no longer on screen. **Proposal**:
`validate()` refuses `override.enable = true` unless the kind is `structured`
or `disv`, with a message saying so. The alternative — the builder ignores it
and warns — leaves a live setting that does nothing, which is the thing this
remodelling exists to remove.
ANSWER: ok

**D4 — blocked boxes keep their values.** When *Refine along the streams* is
off, `cell_near_stream`, `stream_buffer` and the band count are shown greyed
out, **not cleared**, and their values stay in the TOML — so turning the
refinement back on restores what you had. Say if you would rather they be
dropped from the file while the option is off.
ANSWER: if the option no refinement is not selected, the voronoi should show no refinement, so parameters `cell_near_stream`, `stream_buffer` and the band count are shown grey and reset to nothing. When the option is activated, these parameters are comouted as defines previously /default values), ad can be changed by the user.
### 2A.9  Built, 2026-09-12 — and one correction to D2

All of §2A.7 and §2A.8 is implemented. What a reviewer should check against
the code:

| decision | where it lives now |
|---|---|
| voronoi is the default | `marmites_config.Grid.kind`, and `kind = "voronoi"` in `lamata.toml` |
| no `rebuild` key | gone from `[grid]`; `build_mesh(..., force=True)`, which *Create grid* passes |
| voronoi rectangle snapped to `cell_far` | `marmites_meshes.rectangle_cell_size` |
| override refused on a mesh | `RunConfig.validate()` |
| bands derived | `GridVoronoi.bands()` / `.refresh()`, called from `validate()` |
| corridor cleared when refinement is off | `GridVoronoi.refresh()` |
| per-kind sub-panels | `schema.GRID_PERMANENT` / `GRID_SUBPANEL` / `GRID_GATED` / `GRID_DERIVED` |
| Create grid ≠ Select this grid | `code/app/pages/1_Grid.py` |

**The correction.** §2A.8 D2 said 3 bands at 20, 40 and 60 m. That was wrong
by one: the producer spreads the sizes **across** the bands — *n* bands carry
*n* sizes, the innermost `cell_near_stream` and the outermost `cell_far` — so
bounding the ratio takes one more band than it takes intervals. The rule as
built is

    n = ceil( (cell_far - cell_near_stream) / (cell_near_stream * (r - 1)) ) + 1
    band k at  stream_buffer * k / n,   k = 1 .. n

and on the shipped settings (40 → 100 m, 60 m corridor, r = 1.5) that gives
**4 bands at 15, 30, 45 and 60 m, carrying 40, 60, 80 and 100 m cells** —
steepest step exactly 1.5, and the outermost band exactly at the corridor
edge. The rule and every answer are unchanged; only the arithmetic in that
one worked example was off.

**The bound on `r`.** You asked whether > 1 and < 5 is right. Not quite, in
both directions, so it is held in **(1, 3]**:

- it cannot be *> 1* alone, because the band count diverges as *r* → 1: at
  1.01 on a 40 → 100 m corridor the rule asks for 150 bands. There is a
  ceiling of 12, and rather than cap silently — which would mean the mesh
  does not grade the way the file says — a ratio needing more than 12 is
  **refused, naming the number it wanted**.
- 5 is too generous at the top: a 5× jump between neighbouring cells is the
  distortion the grading exists to remove. 3 is already permissive; 1.2–2 is
  the range worth using.

**One consequence of (c) you will meet on the next run.** `hi_spinup_l1.asc`
and `hi_spinup_l2.asc` in `example/LaMata/MF_ws` are 65 × 60 arrays produced
on the structured grid, and there is no scope sidecar for them anywhere. With
`kind = "voronoi"` the state guard stops the run rather than hand MODFLOW an
array of the wrong length for an 830-cell mesh:

    CONFIG ERROR: spinup.strt_heads = 'hi_spinup' has no scope sidecar, so it
    was produced on the structured grid, and grid.kind = 'voronoi' needs state
    on its own mesh.

That is the guard doing its job. Either clear `spinup.strt_heads` and
`spinup.steady_means`, or run a spin-up on the mesh and save it. Panel 1 now
says so **at selection time**, so it is a message when you choose the grid
rather than a stop at the start of the next run.

---

### 2A.10  The rectangle check — a stopgap, 2026-09-13

**The problem, in workflow order.**

1. Panel 1 derives the model rectangle from the catchment polygon. The
   polygon's western edge is at x = 739293.4 and the origin snaps *outward*
   to a whole cell, so the rectangle starts at **739250** — the same for all
   four grid kinds.
2. A run then assembles the *structured* model from the rasters committed in
   the dataset (`inputSOILzones.asc`, `inputVEG*area.asc`, and on the MF side
   `elev.asc`, `thick_l*.asc`, `hk_l*.asc`, `ibound_l*.asc`). Every one of
   them declares **739300** — one 50 m column further east, clipping 6.6 m
   off the catchment. That is a historical rectangle frozen into the exports.
3. `marmites_mesh.project_model` resamples that assembly onto the grid. Both
   are in absolute UTM, so they must stand on the same ground. They do not:
   the grid hangs 50 m west over nothing. There `top` and each `botm` are
   averaged over different subsets of source cells and can cross, and the
   projection stops with *"the cell top at or below its bottom after
   resampling"* — a message about layers, for a fault in rectangles.

**It is not a Voronoi problem.** The rectangle is the same for every kind.
`structured` and `disv` merely have an escape hatch — `grid.override`, which
reproduces the legacy grid and is refused on a mesh by D3 — so they can be
pinned back onto the rasters and a mesh cannot.

**Where it really belongs.** Those rasters are *grid-dependent derived data*
still sitting in the dataset. Under the two-tier rule they should come out of
the converter, from the cartography, onto whatever rectangle this panel
produces. That is step 2 above, it belongs to **panel 3 (Model)**, and it is
on stand-by until panels 1 and 2 are finished. The DEM is already the
exception that proves it: `[grid] dem` is kept at its own 5 m resolution and
wrapped onto the cells at run time, so it covers either rectangle and is
deliberately excluded from this check.

**The stopgap.** `marmites_meshes.dataset_rectangle()` reads the rectangle
the dataset's rasters declare — grouped, largest group wins — and
`rectangle_check()` compares it with the derived one:

| status | meaning |
|---|---|
| `ok` | the grid stands entirely on the rasters |
| `overhang` | it reaches past them — the model build will fail, with the distance per side |
| `shifted` | covered, but off their lattice, so every cell resamples from fractions of four |
| `none` | the dataset holds no raster yet — a new catchment has nothing to disagree with |
| `error` | the rectangle itself cannot be derived (no polygon, no override) |

Panel 1 draws it **above *Create grid***, because it is about the rectangle
the button would use. On a structured or disv grid it also offers *Pin the
grid to the rasters' rectangle*, which fills the override from the rasters'
own header; on a mesh it is not offered, because there the answer is panel 3.

On La Mata the check also found two rasters on a **third** rectangle entirely
— `MF_ws/vka_l1_old.asc` and `vka_l2_old.asc`, 69 × 72 @ 40 m at 739325 —
which is why disagreeing rasters are reported rather than out-voted.

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
