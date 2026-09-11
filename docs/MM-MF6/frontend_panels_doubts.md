# Panels 0-4 — what I verified, and what I need before implementing

Reply to your edits of `frontend_organisation_draft.md`, 2026-09-11.
Nothing implemented yet. Part A answers the two things you asked me to
check; part B is what I found that changes the plan; part C is the list of
doubts I want settled first.

---

## A  Answers to what you asked

### A.1  §1.5 — the format of the MMsurf time-series inputs

Read from `startMARMITESsurface.py` (lines 420-560). All three are read with
`np.loadtxt(..., skiprows=1)`, so **the header line is compulsory but never
parsed** — column order alone carries the meaning.

**`__meteoTB.txt` — ONE file, columns side by side, six per meteo zone.**

```
Date  Time | P Ta RHa Pa u_z_m Rs | P Ta RHa Pa u_z_m Rs | ...
           |  <--- zone 1 --->    |  <--- zone 2 --->    |
```

Indexing is `dataMETEOTS[:, 2 + n*6 .. 7 + n*6]`, so for NMETEO = 2 the file
has 14 columns. Not sequential blocks, not one file per station.

| col | name | units |
|---|---|---|
| 1 | Date | `YYYY-MM-DD` (parsed by `datestr2num`) |
| 2 | Time | `HH:MM` — **hourly data, the whole PM/FAO56 chain assumes it** |
| 3 | `RF_mm` (P) | mm, rainfall over the hour |
| 4 | `AirT` (Ta) | °C, at height `z_h` |
| 5 | `RH` (RHa) | %, at height `z_h` |
| 6 | `Pa` | kPa |
| 7 | `WS` (u_z_m) | m/s, at height `z_m` |
| 8 | `S_in` (Rs) | MJ·m⁻²·h⁻¹ |

Two traps worth surfacing in the UI: only `DTS[0]` is applied — **the data
time shift of zone 1 is used for every zone** — and nothing checks that the
columns are complete, so a missing station silently shifts all the ones
after it.

**`__IRR_TS.txt` — ONE file, one column per FIELD.**

```
Date  Time | FIELD1 | FIELD2 | ...        (read as [:, 2+n], n = 0..NFIELD-1)
```

The dates are **never parsed** — only columns 2+ are read — so the file is
assumed row-by-row aligned with `__meteoTB.txt` and nothing verifies it.
The committed example even uses a different date format (`31/05/2008` against
the meteo file's `2008-05-31`) and no one notices. The front-end should
validate row count and first/last date against the meteo file, and refuse
otherwise.

**`__inputFIELD<n>_crop_schedule.txt` — ONE FILE PER FIELD, name hardcoded.**

`'__inputFIELD%d_crop_schedule.txt' % (f+1)`, in `MMsurf_ws`, so NFIELD = 3
needs `__inputFIELD1..3_crop_schedule.txt` and the name is not configurable.
Five columns, tab separated, one row per growing season:

| col | name | meaning |
|---|---|---|
| 1 | `StartDate` | sowing |
| 2 | `EndDate` | harvest |
| 3 | `GrowingDuration` | d, ramp from bare to full crop |
| 4 | `WiltingDuration` | d, ramp back down |
| 5 | `CROP` | crop index, 1-based into the NCRP list |

Validated: start < end, and season *i*-1 must end before season *i* starts.

**Proposal for the front-end.** Keep the wide-table layout — it is what the
reader expects — but present it as an editable grid (one tab per station,
columns named), write the file from the UI, and check on entry: hourly
spacing, no gaps, row alignment across the three files, `CROP` within
`1..NCRP`, and schedule coverage against the meteo record. The panel carries
a short "expected format" block with the column tables above.

### A.2  Question 5 — do the time-series plots start on the hydrologic year?

**No. Your setting is ignored on the MF6 path.** `iniMonthHydroYear` is read
from the MM ini into `MMConfig`, but `MMConfig` is never constructed on a
run, and nothing ever sets the attribute on `cMF`. Every call site is

```python
int(getattr(cMF, 'iniMonthHydroYear', 10))      # 3 sites in marmites_postprocess.py
```

so the value is **hardcoded to October** by the default. La Mata's ini also
says 10, which is why this has never shown. Change it to 9 today and nothing
moves. It becomes a real field in panel 4 and has to be threaded through — it
drives the x-axis major locator in ~20 plots plus the Sankey year index.

---

## B  Findings that change the plan

### B.1  The committed MMsurf is OLDER than the forcing files it supposedly wrote

`startMARMITESsurface.py` writes `inputZONP_veg_d.txt` / `inputZONP_irr_d.txt`
and labels them `P zones` in the handover. The committed handover and the
committed forcing say `inputZONRF_veg_d.txt` / `RF zones`. It also never
writes `inputZONRFe_*`, yet two such files are committed — and
`inputZONRFe_veg_d.txt` has 4869 data lines, which is not a multiple of the
1949-day record, i.e. it is **truncated garbage**.

The numbers differ too: the ini says `kTg_min = 0.04440913`, `kT_s =
21.57149438`; the handover says `0.044694325` and `21.98151201`.

So "recover MMsurf" is not only re-wiring. The module in the repo cannot
reproduce the forcing the model runs on today, and we will not know it is
right until it does. I propose making that the acceptance test: **run the
recovered MMsurf on `__meteoTB.txt` and require it to reproduce the committed
`inputZON*_d.txt` to within rounding**, then delete them from the repo.

### B.2  `Pe` means throughfall, and `kT_s` is stored inverted

Two naming traps to fix in the labels rather than in the code:

- the driver reads the handover positionally and binds slot 2 (`TF`) to the
  variable `Pe_veg_fn`. So everywhere downstream `Pe` = **throughfall**, not
  effective rainfall. `RFe` is not used at all.
- the ini's 19th vegetation field is commented `kT_1/s` and holds 21.57; the
  driver does `kT_s = 1.0/float(x)`. So the file stores **1/s** and the model
  wants **s**. The front-end must say which one the user is typing. I suggest
  exposing `s` (0 < s < 1, as the comment says) and inverting on write.

### B.3  NVEG is silently +1

MMsurf prepends `grassFAO56` with hardcoded FAO-56 parameters and returns
`NVEG+1`; the handover writes `NVEG-1` to undo it. So "3 vegetation types" in
the UI means four internally, and vegetation index 0 is the reference grass.
The area rasters `inputVEG<v>area.asc` are 1-based over the **user's** types.
Worth stating in the panel, not hiding.

### B.4  MMsurf writes its output into the repository

Every `.out`, `.int` and `.png` goes to `pathMMsurf` =
`example/LaMata/MMsurf_ws/`, and the `paramTS` helper writes a PNG per
parameter per vegetation type on every call. That is run output inside the
repo, against the rule. It moves to `WS_ROOT/<case>/MMsurf_out/` when MMsurf
is recovered.

### B.5  The vector layers you want already exist, with the right attributes

`E:\00code_ws\LAMATA_new\GIS\Soil_type.shp` — 27 polygons, ED50 / UTM 29N:

```
SoilType | SoilCode | ibound_l1 | ibound_l2 | iuzfbnd | PONDhmax | PONDw | SOILthick
```

That is *exactly* the layer §2.2 describes, and it also carries the MF ibound
and the UZF boundary — so part of group 3 vectorises for free when we get
there. Note `PONDhmax` / `PONDw` live **on the soil polygons**, not on the
stream lines, which is the opposite of what §2.2 proposes (doubt C.2).

Vegetation: `lm_veg.shp`, 15586 crown polygons with a `Species` column — the
area fraction per cell is an area-weighted overlay of those against the grid,
the same clipping machinery `marmites_mesh.py` already has.

Observations: `202109ObsPts.shp` / `202109MonitPts.shp` carry
`Name, X, Y, lay, hi, h0, RC, STO, NameReal, onMap` — `inputObs.txt` as a
point layer, already.

Irrigation: `Irr_Fields.shp`, 3 polygons, attribute `Id` ∈ {0, 1} — not a
field index, so it cannot map to `__inputFIELD<n>_crop_schedule.txt` as it
stands (doubt C.6).

### B.6  Panel 1 inverts today's dependency, and that breaks the regression anchor

Today the grid is **declared**: `nrow = 65`, `ncol = 60`, `delr = delc = 50`
in the MF ini, `xllcorner/yllcorner` hardcoded in `run_lamata_mf6.py`, and
every `.asc` must match those dimensions or the reader fails. Panel 1 makes
the grid **derived**: catchment polygon + cell size → origin, extent, shape.

Derived from `Limite.shp`, a 50 m grid does not land on the same origin, so
the committed rasters, the spin-up state `hi_spinup_*.asc` and the WP0
byte-identical regression all stop matching. That is not a reason to avoid
the change — it is the right design — but it needs a decision (doubt C.1).

---

## C  Doubts — please settle these before I start

**C.1  Do we keep a way to reproduce the legacy 65 × 60 @ 50 m grid exactly?**
My proposal: panel 1 derives the grid from the polygon by default, but keeps
an explicit override block (`origin`, `nrow`, `ncol`) that reproduces today's
grid byte-for-byte, so the WP0 regression and rung (a) of the validation
ladder stay usable. Without it we lose every baseline at once. Keep or drop?

**C.2  §2.2 sentence is cut off:** *"To produce them from the shape file, the
script should use"* — use what? And related: `PONDhmax` / `PONDw` are today
attributes of the **soil polygons**, not the stream lines. You propose moving
them to the hydrography layer, one value per segment. Do you want me to (a)
add those two columns to the hydrography table in the converter, (b) keep
reading them from `Soil_type.shp`, or (c) make the source configurable per
parameter?

**C.3  Is MMsoil's `Ssurfw` the same quantity as SFR's channel width?**
Today one raster, `inputSTREAMw.asc`, feeds both `gridSsurfw` (the MMsoil
surface store) and `b.sfr_pondw` (the SFR network builder). SFR already has
its own `sfr.width` with a drainage law `w = a·A^b`. Once width is a
per-segment attribute, do these become one parameter, or stay two that may
disagree?

**C.4  Soil thickness precedence reads backwards to me.** You wrote that the
polygon may carry `thick_m`, and then *"if the inputSOILthick.asc exists, it
will be the default value"* — i.e. the raster wins over the attribute. Is
that intended? I would rather not decide by file existence at all: make it an
explicit source (`value` | `column` | `raster`), the same `ParamSource`
pattern WP0 already uses for SFR, and error if it is unset. Agreed?

**C.5  Meteo Voronoi needs station coordinates that do not exist yet.**
`__inputMMsurf.ini` holds `phi` / `Lm` — latitude and longitude **in
degrees**, used for the solar geometry of Penman-Monteith — not projected
coordinates. To build Thiessen polygons I need each station's x, y in the
project CRS. Add two fields per station in panel 2, or a stations point
shapefile? (With NMETEO = 1 it does not matter; I would still require it so
NMETEO = 2 does not surprise us later.)

**C.6  Irrigation fields.** `Irr_Fields.shp` has only `Id` ∈ {0,1} over three
polygons, so it cannot carry the field index that selects
`__inputFIELD<n>_crop_schedule.txt`. Shall I require a `field_id` column
(1..NFIELD, 0 = not irrigated) and have the front-end refuse a layer without
it?

**C.7  Do the shapefiles enter the repository?** The rule is code + docs +
strictly what MM/MF read. If panels 1-3 take vector layers as the primary
input, I would keep WP1's two tiers: shapefiles stay in `DATA_ROOT/GIS`
(outside the repo), the converter writes grid-independent tables into
`example/<case>/`, and **only those tables are committed**. The front-end
browses the GIS folder and drives the converter. Confirm that is what you
mean by "provided as vectorial layers", because the other reading — commit
the shapefiles — contradicts the rule.

**C.8  Scope of this task.** I read your list as: build panels 0, 1, 2, 4 and
the MMsoil half of panel 3; vectorise soil, vegetation, meteo zones,
irrigation and the stream attributes; leave the MF rasters (`hk_l*.asc`,
`ibound_l*.asc`, `thick_l*.asc`, `elev_sinkfil.asc`, `Ss` / `Sy`) as rasters
until you define the MF vector scheme. Correct?

**C.9  Panel 3 has no master switch left.** You removed `MMsoil_yn` and
`MF_yn` as redundant (MMsoil and MF always run together). So panel 3 is
unconditional and only panels 2 and 4 carry switches. Fine — or do you want a
single "run model" switch there anyway, so a surface-only run is possible?

**C.10  One new module, or extend `marmites_mesh.py`?** The vector→grid
wrapping is a genuinely new operation: polygon overlay (majority / area
fraction / area-weighted mean), line burning with per-segment attributes, and
point location. `marmites_mesh.py` already has the clipping and the cell
lookup, but it resamples *rasters onto a mesh*. I propose a sibling
`marmites_vector.py` with the same `how=` vocabulary, reusing the
projection's cell index. Any objection?

---

## D  Your answers, and the layer inventory  (2026-09-11)

Settled: **C.2** hydrography carries the stream attributes · **C.1**
`lm_lim.shp` is the catchment · **C.5** convert lat/long to the project CRS ·
**C.7** shapefiles stay in `DATA_ROOT/GIS`, never in the repo · **C.4** the
raster beats the polygon attribute.

Still open: **C.3** (is MMsoil `Ssurfw` the same as SFR `width`), **C.8**
(scope stops before the MF rasters), **C.9** (does panel 3 need a switch at
all), **C.10** (a new `marmites_vector.py`).

### D.1  What I found, per input

| front-end input | layer in `DATA_ROOT/GIS` | state |
|---|---|---|
| catchment boundary (panel 1) | `lm_lim.shp` | **OK** — 1 polygon, 4.844 km², ED50 / UTM 29N |
| stream network geometry | `hydrography.shp` | OK as geometry — 97 lines, attributes missing (D.2.1) |
| ponds | `lm_ponds.shp` | OK — 12 polygons |
| soil zones | `Soil_type.shp` `SoilCode` 1/2/3 | **OK** — matches `inputSOILparam.txt` order (alluvium, regolith, outcrop) |
| soil thickness | raster `inputSOILthick.asc`, fallback `Soil_type.SOILthick` | OK |
| ibound l1/l2, iuzfbnd (MF, later) | `Soil_type.shp` | OK |
| vegetation crowns | `lm_veg.shp` | needs a mapping (D.2.4) |
| observation points | `202109ObsPts.shp` + `202109MonitPts.shp` | needs a decision (D.2.5) |
| irrigation fields | `Irr_Fields.shp` | needs a field index (D.2.3) |
| meteo station position | none | needs a decision (D.2.2) |

**`lm_lim.shp` also corrects WP1.** `inputWATERSHED.csv` is currently built
from `Limite.shp`, which is 19.159 km² — nearly four times the catchment, and
the reason `_clip_to_grid` had to exist. `lm_lim.shp` is 4.844 km², against
4.885 km² for the 1954 active cells of today's model: the same domain. I will
switch the converter over.

Its bbox is x 739293..742223, y 4553110..4556240 — **inside** today's model
rectangle (739300..742300, 4553050..4556300), and 7 m west of its western
edge. So a grid derived from the polygon is not today's grid, and unless you
say otherwise I will keep the explicit origin/shape override from C.1 so the
WP0 regression and rung (a) of the ladder stay reproducible.

### D.2  What is missing — please say where these live

**D.2.1  Stream width and max height, per segment.** `hydrography.shp`
carries only `GRID_CODE` (1 or 2). Today the values come from
`Soil_type.shp`, where they are uniform per soil class:

```
Alluvium (n=2)   PONDhmax 1.0   PONDw 1.5   SOILthick 1.50
Regolith (n=2)   PONDhmax 0.0   PONDw 0.0   SOILthick 0.75
Outcrop  (n=23)  PONDhmax 0.0   PONDw 0.0   SOILthick 0.05
```

So `inputSTREAMw.asc` is really **the alluvium footprint**, not the stream
lines, and what MMsoil treats as its surface-water network is that footprint.
Moving these onto the segments is a genuine improvement, but it will change
the numbers. Is there another hydrography layer that carries width, or shall
I add `w_m` and `hmax_m` columns to `hydrography.shp` (with one constant for
the whole network as the fallback, which is what the data says today)?

**D.2.2  The meteo station position.** Converting the ini values as you asked
works, but the result is outside the catchment:

| source | result | inside the model rectangle |
|---|---|---|
| `phi 41.045`, `Lm 6.16 W` from WGS84 | x 738823, y 4547854 | **no — 5.2 km south** |
| the same from ED50 (EPSG:4230) | x 738710, y 4547718 | **no** |
| model centre, back-converted | lat 41.1058, lon 6.1339 W | — |

The ini values are 6.8 km south and 2.2 km west of the catchment centre, so
they are not the station's real position — they only have to be roughly right
for the solar geometry of Penman-Monteith. The `EC` point of
`202109MonitPts.shp` — x 739625, y 4555925, the eddy-covariance tower — is
inside. Which one is the meteorological station? At NMETEO = 1 the Thiessen
polygon is the whole catchment either way, so this only bites at NMETEO > 1;
I would rather store the right thing now.

**D.2.3  Irrigation field index.** `Irr_Fields.shp` has 3 polygons and one
attribute, `Id` in {0, 1}, so nothing selects
`__inputFIELD<n>_crop_schedule.txt` or the matching `__IRR_TS.txt` column.
Add a `field_id` column (1..NFIELD, 0 = not irrigated)?

**D.2.4  Vegetation species to VEG index.** `lm_veg.shp` holds 15586 crown
polygons with `Species`: **`i` 14464, `p` 921, blank 201**. I read that as
`i` to VEG2 (Qilex), `p` to VEG3 (Qpyr), and VEG1 (grassMU) as the remainder
of each cell — which is what the name `inputVEG1areaNoGRASS.asc` implies. Two
things to confirm: that mapping, and what the 201 blank ones are — drop them,
or treat them as ilex?

**D.2.5  Observation points: two layers, and no active flag.**
`202109ObsPts.shp` (O1, O2, I1, G1, G2) and `202109MonitPts.shp` (P0, SM,
C1..C5, EC) are 13 points together, carrying `Name, X, Y, lay, hi, h0, RC,
STO, NameReal, onMap`. `inputObs.txt` lists the same 13 but **comments some
out** (`W1`, `C4`, `C5`, `H3`, `H4`, `C6`, `I2`, `G3` — several of which are
not in the shapefiles at all). Do I merge the two layers and add an `active`
column, or is `onMap` already that flag?

---

## E  The overlay, validated against the rasters it replaces  (2026-09-11)

`code/marmites_vector.py` is written and tested (28 unit tests). The real
test is whether it reproduces the ASCII rasters, which were made in ArcGIS
from these same layers. Run on today's 65 × 60 @ 50 m grid:

| what | vector source | vector | raster | verdict |
|---|---|---|---:|---|
| Qilex area | `lm_veg.Species = 'i'` | mean **3.847 %** | 3.847 % | **exact** (max diff 0.014 pp) |
| Qpyr area | `lm_veg.Species = 'p'` | mean **0.374 %** | 0.374 % | **exact** (max diff 0.005 pp) |
| grass area | `lm_veg.Species = ''` | mean 45.221 % | 12.688 % | **wrong source** (E.1) |
| soil zone | `Soil_type.SoilCode`, majority | 1831/1962 cells agree | | E.2 |
| soil zone | `SOIL_map.SoilClas`, majority | 1342/1962 cells agree | | worse; `Soil_type` is the layer |
| soil thickness | `Soil_type.SOILthick`, area mean | mean 0.789 m | 0.612 m | confirms C.4 — the raster is a real thickness map, not per-zone constants |
| streams | `hydrography.shp` burned | 350 cells, 14 536 m | 244 cells | 234 shared, 116 vector-only |

The two crown classes coming out **exact** against an independent ArcGIS
computation is the evidence that the clipping, the ring handling and the
area weighting are right. Everything else below is about which layer to
point at, not about the machinery.

### E.1  grassMU has no vector source that I can find

Neither candidate matches the 12.688 % of `inputVEG1area.asc`:

| layer | blank-`Species` coverage |
|---|---|
| `lm_veg.shp` (15586 crowns) | 45.22 % |
| `lm_veg_raster.shp` (22471, already cut by the 50 m cells) | **95.78 %** |

`lm_veg_raster.shp` is `lm_veg` intersected with the model grid — its
`CODveg` is 0 / 1 / 2 for blank / `i` / `p`, and its ilex and pyr totals
(375 098 m² and 36 457 m²) reproduce the rasters to three digits. So in both
layers "blank" means **not a tree crown**, which is 96 % of the catchment,
not the 12.7 % of grass. Grass was mapped somewhere else.

Candidates in the GIS folder I have not been able to open as vectors:
`ilexveg1area`, `pyrveg2area`, `lm_veg_area.tab`,
`lm_veg_raster_PROCESSING.xlsx`. Note their numbering is the other way round
from the model's (`ilexVEG1area` against `inputVEG2area.asc`), which is its
own trap. **Where does the grass map live?**

### E.2  Soil zone: 93 % agreement, and the gap is alluvium

`Soil_type.shp` by majority vote against `inputSOILzones.asc`:

```
zone 1 alluvium   vector  167 cells   raster  258 cells
zone 2 regolith   vector 1840 cells   raster 1620 cells
zone 3 outcrop    vector   86 cells   raster   84 cells
```

The alluvium is a ribbon along the streams, mostly narrower than a 50 m
cell, so it rarely wins a majority vote — while the raster gives it 645 000 m²,
more than the polygon's own 573 644 m². This is the WP1c.3 trade-off in a new
place: majority is faithful per cell and loses minority classes in aggregate.
Since the grid is changing anyway (panel 1) exact reproduction is not
required, but a third fewer alluvium cells changes the soil water balance, so
I want you to see it rather than discover it later.

---

## F  MMsurf recovered, and proven faithful  (2026-09-11)

Finding B.1 said the MMsurf in the repository could not reproduce the forcing
it supposedly wrote. It can now, and the proof is exact.

**Run 1 — the recovered MMsurf, driven by the configuration.** Rainfall, LAI
and the crop schedule come out identical; the residual is 0.009 % on the mean
and confined to the three quantities that pass through Penman-Monteith:

| file | mean new | mean old | max diff |
|---|---:|---:|---:|
| `inputZONRF_veg_d.txt` | 1.201939 | 1.201939 | **0** |
| `inputZONLAI_veg_d.txt` | 3.516259 | 3.516259 | **0** |
| `inputZONRF_irr_d.txt` | 1.807574 | 1.807574 | **0** |
| `inputZONcrop_irr_d.txt` | 0.617753 | 0.617753 | **0** |
| `inputZONTF_veg_d.txt` | 0.880196 | 0.880193 | 1.2e-04 |
| `inputZONPT_veg_d.txt` | 0.514251 | 0.514252 | 1.8e-02 |
| `inputZONPE_d.txt` | 1.577433 | 1.577571 | 1.2e-01 |
| `inputZONEo_d.txt` | 1.761934 | 1.762088 | 1.3e-01 |

**Run 2 — the same, with the ini's station position pinned back** (41.045 N /
6.16 W instead of the 41.117 N / 6.149 W that 739508, 4555882 projects to).
Every one of the ten files reproduces **exactly, max diff 0**.

So the residual is the station move and nothing else: net radiation depends
on the station's latitude and longitude, and the station was 6.8 km south of
the catchment. The recovered module is faithful to the one that produced the
committed dataset, and the only behavioural change is a correction.

This also closes B.2's second half in passing. The writer labels its output
TF (throughfall) while `PET_P_INTER` computes it as `RFe_veg_d = RF - I` --
effective rainfall. They are the same quantity under two names, so the
driver's `Pe` is not a misnomer after all; only the file name is inconsistent
with the variable that fills it.
