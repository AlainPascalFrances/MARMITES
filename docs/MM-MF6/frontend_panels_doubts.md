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
