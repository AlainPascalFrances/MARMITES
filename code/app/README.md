# Streamlit interface  (WP1b)

A browser front-end to browse the input data (shapefiles included), edit and
validate a run configuration, launch a run locally or on the server, follow it,
and present the output.

```
conda env create -f environment-ui.yml      # once -- a SEPARATE environment
conda activate mm-ui
streamlit run code/app/Home.py
```

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

`loaders.py` and `runs.py` import no Streamlit, so they are unit-tested in the
**model** environment (`code/tests/test_app_lib.py`) — which is also the only way
they could be tested at all, given Streamlit deliberately lives elsewhere. The
pages are thin views over them.
