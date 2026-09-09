# Streamlit interface  (WP1b)

A browser front-end to browse the input data (shapefiles included), edit and
validate a run configuration, launch a run locally or on the server, follow it,
and present the output.

```
streamlit run code/app/Home.py
```

Streamlit is installed in the **`flopy`** env on this machine. That was checked
before doing it: the conda solve touches nothing in the model stack — no
`numpy`, `flopy`, `rasterio`, `h5py`, `gdal`, `matplotlib`, `geopandas`,
`shapely` or `pyproj` is superseded, downgraded or changed; only
`ca-certificates`, `certifi` and `openssl` are updated. `environment-ui.yml` is
the recipe for a fresh machine, or a server where you would rather keep the
front-end separate.

Serving from the server: run it **there** and reach it over the LAN —
`streamlit run code/app/Home.py --server.address 0.0.0.0 --server.port 8501`.
The PEST workers already run there and the app only needs to see `$MM_WS_ROOT`.

## Layout

```
Home.py                 entry point, case selector, machine banner, recent runs
pages/1_Inputs.py       Tier A inventory + Tier B cartography on a map + the converter
pages/2_Configuration.py  schema-driven editor; MMsurf parameters pre-filled (WP1.7)
pages/3_Run.py          launch (detached) and follow the log
pages/4_Results.py      the figures a run wrote, with a run picker
lib/loaders.py          cached readers, keyed on (path, size, mtime)
lib/runs.py             the run registry: launch, status, log tail
```

## Two rules this app lives by

**It is a viewer and a launcher, never part of the model path.** `code/app/` may
import the model; the model may never import `code/app/` or `streamlit`.
`code/tests/test_repo_hygiene.py` asserts it, so the model keeps running headless
from Spyder and from a PEST worker where Streamlit is not installed.

**No model run happens inside a page.** Streamlit reruns the whole script on
every widget interaction, and a coupled run is ~11.5 minutes (an IES run is
hours). `lib/runs.py` launches a **detached** process that writes
`<runs_dir>/<run_id>/run.log` and `status.json`; the page polls them. A run
therefore survives closing the browser.

## Why the substance is in `lib/`

`loaders.py` and `runs.py` import no Streamlit, so they are unit-tested in
`code/tests/test_app_lib.py` like any other module. The pages are thin views
over them.

## Three Streamlit traps this app already hit

Recorded because each one renders *something* plausible rather than failing
loudly, so none of them shows up as an error:

1. **`st.button` for anything that must persist.** A button is `True` only on
   the rerun its own click causes, so the map vanished the moment any checkbox
   was touched. Use `st.toggle` (or session state) for state.
2. **`st.text_area(value=..., key=...)`.** With a `key`, session state becomes
   authoritative after the first render, so the log froze at whatever it was
   when the page first drew it — empty, for a run just started. The log is
   `st.code` now.
3. **`folium.Map()` with no location** centres on (0, 0) and opens on the whole
   world. Compute the layer bounds and `fit_bounds` them.
