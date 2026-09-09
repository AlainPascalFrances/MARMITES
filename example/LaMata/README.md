# La Mata — Tier A inputs

Everything MM and MF read **directly**, and nothing else. No GIS, no run output,
however small or convenient — `code/tests/test_repo_hygiene.py` enforces it.

- **Raw cartography** (shapefiles, the 5 m DEM) lives outside the repository in
  `$MM_DATA_ROOT/GIS` (default `E:\00code_ws\LAMATA_new\GIS`).
- **Run output** goes to `$MM_WS_ROOT` (default `E:\00code_ws\LaMata_MM-MF6`).
- The `*.csv` files below are **generated** from the cartography by
  `python code/tools/gis_to_dataset.py --case LaMata`, and each carries a
  provenance header naming its source, that file's size and mtime, and its CRS.

Grid: 60 columns × 65 rows, 50 m cells, lower-left corner (739300, 4553050),
**EPSG:23029** (ED50 / UTM 29N). Note the DEM rasters in the GIS folder carry an
equivalent but *unnamed* PROJCS with no EPSG code — the converter reports that
explicitly rather than guessing silently.

## Run control and parameters

| File | Content | Read by |
|---|---|---|
| `__inputMM_v3.ini` | the legacy MM run-control block. **Superseded by `code/configs/lamata.toml`** (WP0); kept because the MMsurf and observation file names still come from it | `ppMODFLOW_flopy_v3`, `startMARMITES_v3` |
| `MF_ws/__inputMF_flopy_v3_2s1L.ini` | the **2-layer** MODFLOW parameter set — authoritative, hand-set, not a derived aggregation | `clsMF`, `setup_lamata` |
| `MF_ws/__inputMF_flopy_v3_2s3L.ini` | the 6-layer parameter set. **Written for MODFLOW-NWT; obsolete for MF6** and re-expressed in the `[layers]` / `[uzf]` / `[seep]` config sections | as above |
| `MF_ws/inputSOILparam.txt` | per soil zone: `nsl`, and per soil layer `st`, `slprop`, `Smax`, `Sfc`, `Sr`, `Si`, `Ks`. The file the WP7 `soil_*` parameter group calibrates | MMsoil |
| `MMsurf_ws/__inputMMsurf.ini` | MMsurf parameters: meteo-station geometry, vegetation (incl. `Zr` rooting depth and the `kTg_*` transpiration-sourcing factors), crops, soils, open water | MMsurf |
| `__inputMMsurf4MMsoil.txt` | the MMsurf → MMsoil handover: counts, series file names, `Zr`, `kTg_min/max`, `kT_f`, `kT_s` | `setup_lamata` |

## Time and forcing series

| File | Content |
|---|---|
| `inputDATE.txt` | the stress-period calendar |
| `inputZON*_d.txt` | **daily** per-zone series: `P`, `Pe`, `RF`, `RFe`, `TF`, `PT`, `PE`, `E0`, `Eo`, `LAI`, crop — for vegetation (`_veg`) and irrigation (`_irr`) |
| `inputZON_*_stp.txt` | the same series aggregated per stress period |
| `MMsurf_ws/__meteoTB.txt` | the raw meteorological table MMsurf reads |
| `MMsurf_ws/__IRR_TS.txt`, `__inputFIELD1_crop_schedule.txt` | irrigation series and crop schedule |

## Spatial parameter rasters (ESRI ASCII, model grid)

| File | Content |
|---|---|
| `inputMETEOzones.asc`, `inputSOILzones.asc`, `inputIRRzones.asc` | integer zone maps |
| `inputVEG{1,2,3}area.asc` | fractional vegetation cover per type |
| `inputSOILthick.asc` | soil-column thickness (m) |
| `inputSTREAMw.asc` | **channel width (m)** — feeds `gridSsurfw` |
| `inputSTREAMhmax.asc` | **channel depth (m)** — feeds `gridSsurfhmax` |

> **Renamed in WP1.2.** These two were `inputPONDw.asc` / `inputPONDhmax.asc`, a
> legacy misnomer: they hold the **stream network**, not ponds — which is why the
> arrays they feed have always been called `gridSsurfw` / `gridSsurfhmax`.
> `mm_paths.resolve_input()` still accepts the old names for one release.

## Aquifer grids (`MF_ws/`, ESRI ASCII)

| File | Content |
|---|---|
| `elev.asc` | land-surface elevation |
| `elev_sinkfil.ASC` | **sink-filled** DEM on the model grid — the DEM the SFR routing and CRR slopes use |
| `ibound_l*.asc` | active-cell mask per layer |
| `hk_l*.asc`, `thick_l*.asc` | horizontal K and thickness per layer |
| `Sy_l*.asc`, `Ss_l*.asc` | specific yield and storage |
| `uzf_iuzfbnd.asc` | the UZF footprint |
| `drn_cond_l*.asc`, `drn_elev_l*.asc` | drain conductance and elevation |
| `ghb_cond_l*.asc`, `ghb_head_l*.asc` | general-head boundary |
| `hi_spinup_l*.asc`, `hi_spinup_perc.asc`, `hi_spinup_etg.asc` | the committed **baseline** spin-up state, so `spinup.strt_heads` works on a fresh clone. A newer state in the workspace takes precedence |
| `*_old.asc`, `*_uniform.asc`, `hi_lastrun.asc`, `hi_topL1.asc` | earlier parameter variants, kept for reference |

## Observations

| File | Content |
|---|---|
| `inputObs.txt` | the observation points: name, X, Y, layer (a leading `#` disables one) |
| `inputObsHEADS_{C1..C5,P0,W1}.txt` | measured heads (date, m) |
| `inputObsSM_P0.txt`, `inputObsSM_SM.txt` | measured **soil moisture per probe depth** — the WP7 `sm` calibration group |
| `inputObsRo_catchment.txt` | catchment runoff / streamflow |
| `inputObs4GIS.txt` | the same points for mapping |

## Generated from the cartography (WP1.1)

| File | Content | Source |
|---|---|---|
| `inputSTREAM.csv` | stream geometry as `seg_id, seq, x, y` — **grid-independent**, in projected coordinates, so it survives the WP1c grid change | `GIS/hydrography.shp` (97 segments) |
| `inputSTREAM_param.csv` | per-segment `grid_code`, `length_m`, and the resolved `width_m`, `manning`, `rhk`, `rbth`. A value of `drainage:a=..,b=..` means the **model** computes `w = a·A^b` from its own flow accumulation — drainage area is a property of the routing, not of the shapefile | same |
| `inputPONDS.csv` | `fid`, centroid, `area_m2`, `perimeter_m`, and `dem_mean_m` / `dem_min_m` sampled over each footprint. The model sets rim = `dem_mean` and bottom = rim − pond depth | `GIS/lm_ponds.shp` (12 ponds, 341–2036 m²) |
| `inputWATERSHED.csv` | the catchment boundary as `ring_id, seq, x, y` — the Voronoi domain for WP1c | `GIS/Limite.shp` |

Regenerate with:

```
python code/tools/gis_to_dataset.py --case LaMata --dry-run   # preview
python code/tools/gis_to_dataset.py --case LaMata             # write
```

The converter is the only place besides the Streamlit app that imports geopandas
and rasterio; the model path stays free of them (decision D3).
