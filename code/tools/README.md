# One-off tools  (WP1, WP7)

Utilities that are run **by hand**, never by the model, so their dependencies stay off the
model path:

- `gis_to_dataset.py`  — reads the shapefiles and rasters in `$MARMITES_DATA_ROOT` and writes
  the small plain-text tables the model reads into `example/<case>/`. The only place
  geopandas and rasterio are imported besides the Streamlit app.
- the PEST chain — `pest_prep_mm.py`, `build_pst_mm.py`, `forward_run_mm.py`, `run_ies_mm.py`,
  `postproc_ies_mm.py`, `extract_pest_optimised_mm.py`.

See the cookbook, WP1 and WP7.
