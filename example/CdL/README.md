# Casa de Lobos (CdL) — the production case study

Empty for now. This folder will hold the input data MM and MF read directly for CdL, in the
same spirit as `example/LaMata/`.

The MODFLOW 6 side already exists in the `MF6models` repository
(github.com/AlainPascalFrances/MF6models, `CdL/`): a Voronoi DISV grid (ncpl 6609, 3 layers),
monthly 1981-2026, with UZF / SFR / LAK / MVR / CRR and a pyEMU + pestpp-ies calibration.

What is missing before MARMITES can run here, and what the cookbook plans:

- **Daily meteorological forcing.** CdL currently carries MONTHLY precipitation and ET0 only
  (`p_month_*.csv`, `et0_month_*.csv`). MMsurf needs hourly-to-daily meteorology to produce
  the daily P / Pe / RF / RFe / TF / PT / PE / E0 / LAI series MMsoil consumes. This is the
  single blocking item.
- **Soil and vegetation zoning** in the MARMITES sense. The raw material exists — the soil
  hydraulic rasters (`ks`, `ths`, `wp`, `fc`) and the COS-2025 land cover — but it has to be
  turned into MM soil zones, vegetation areas and an `inputSOILparam.txt`.
- **A grid ingestion path** for the Voronoi DISV grid: MARMITES builds its DIS/DISV grid from
  ASCII rasters, whereas CdL's grid comes from `voronoi_grid.pkl` + `voronoi_layers.npz`.
- **Observations** mapped to the MARMITES formats: heads, soil moisture per depth,
  streamflow, actual ET.
