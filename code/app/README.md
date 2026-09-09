# Streamlit interface  (WP1b)

A browser front-end to browse the input data (shapefiles included), edit and validate a run
configuration, launch a run locally or on the server, follow it, and present the output.

    streamlit run code/app/Home.py

It is a **viewer and a launcher**. `code/app/` may import the model; the model must never
import `code/app/` or `streamlit` — the test suite asserts this, so the model keeps running
headless from Spyder and from a PEST worker where Streamlit is not installed.

Planned pages: `1_Inputs`, `2_Configuration`, `3_Run`, `4_Results`, `5_Calibration`.
See the cookbook, WP1b.
