# -*- coding: utf-8 -*-
"""Build the MARMITES x CdL merge cookbook (.docx) with python-docx.

The cookbook is a DERIVED document: edit this script, re-run it, and commit
both. Editing the .docx by hand means the next regeneration silently discards
the change, so the script lives beside the document it produces.

    python docs/make_cookbook.py

Writes next to itself by default; set MM_COOKBOOK_OUT to write elsewhere.
"""
import datetime
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.section import WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

import os as _os

# Written beside this script, so `python docs/make_cookbook.py` refreshes the
# committed copy in place. Set MM_COOKBOOK_OUT to write somewhere else.
OUT = _os.environ.get(
    'MM_COOKBOOK_OUT',
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                  'MARMITES_x_CdL_merge_cookbook.docx'))

ACCENT = RGBColor(0x1F, 0x4E, 0x79)
GREY = RGBColor(0x59, 0x59, 0x59)

doc = Document()

# ---------------------------------------------------------------- page setup
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)      # A4
for m in ("left_margin", "right_margin"):
    setattr(sec, m, Cm(2.2))
sec.top_margin, sec.bottom_margin = Cm(2.0), Cm(2.0)

st = doc.styles["Normal"]
st.font.name = "Calibri"
st.font.size = Pt(10)
st.paragraph_format.space_after = Pt(6)
st.paragraph_format.line_spacing = 1.05

for name, size, col in (("Heading 1", 17, ACCENT), ("Heading 2", 13.5, ACCENT),
                        ("Heading 3", 11.5, ACCENT)):
    s = doc.styles[name]
    s.font.name = "Calibri"
    s.font.size = Pt(size)
    s.font.color.rgb = col
    s.font.bold = True
    s.paragraph_format.space_before = Pt(14 if name == "Heading 1" else 10)
    s.paragraph_format.space_after = Pt(5)
    s.paragraph_format.keep_with_next = True


# ---------------------------------------------------------------- helpers
def shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    sh = OxmlElement("w:shd")
    sh.set(qn("w:val"), "clear")
    sh.set(qn("w:color"), "auto")
    sh.set(qn("w:fill"), hexcolor)
    tcPr.append(sh)


def h(level, text):
    p = doc.add_heading(text, level=level)
    return p


def para(text="", bold=False, italic=False, size=None, space_after=None,
         align=None, color=None):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold, r.italic = bold, italic
    if size:
        r.font.size = Pt(size)
    if color:
        r.font.color.rgb = color
    if space_after is not None:
        p.paragraph_format.space_after = Pt(space_after)
    if align:
        p.alignment = align
    return p


def rich(parts, style=None, space_after=None):
    """parts = list of (text, 'b'|'i'|'c'|'') tuples."""
    p = doc.add_paragraph(style=style)
    for txt, kind in parts:
        if kind in ("", "b"):
            _emit(p, txt, bold=("b" in kind))
            continue
        r = p.add_run(txt)
        if "b" in kind:
            r.bold = True
        if "i" in kind:
            r.italic = True
        if "c" in kind:
            r.font.name = "Consolas"
            r.font.size = Pt(9)
    if space_after is not None:
        p.paragraph_format.space_after = Pt(space_after)
    return p


def bullet(text, level=0):
    p = doc.add_paragraph(style="List Bullet" if level == 0 else "List Bullet 2")
    _mixed(p, text)
    p.paragraph_format.space_after = Pt(3)
    return p


def numlist(text, level=0):
    p = doc.add_paragraph(style="List Number" if level == 0 else "List Number 2")
    _mixed(p, text)
    p.paragraph_format.space_after = Pt(3)
    return p


def _emit(p, text, size=None, bold=False):
    """Inline markup: **bold**, `code`, and `code` nested inside **bold**."""
    import re
    for tok in re.split(r"(\*\*.+?\*\*)", text):
        if not tok:
            continue
        if tok.startswith("**") and tok.endswith("**"):
            _emit_code(p, tok[2:-2], size, True)
        else:
            _emit_code(p, tok, size, bold)


def _emit_code(p, text, size, bold):
    import re
    for tok in re.split(r"(`.+?`)", text):
        if not tok:
            continue
        mono = tok.startswith("`") and tok.endswith("`")
        r = p.add_run(tok[1:-1] if mono else tok)
        r.bold = bold
        if mono:
            r.font.name = "Consolas"
            r.font.size = Pt((size or 10) - 0.7)
        elif size:
            r.font.size = Pt(size)


def _mixed(p, text):
    _emit(p, text)


def body(text):
    p = doc.add_paragraph()
    _mixed(p, text)
    return p


def code(lines, size=8.5):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.space_before = Pt(4)
    pf.space_after = Pt(8)
    pf.left_indent = Cm(0.4)
    pf.line_spacing = 1.0
    for i, ln in enumerate(lines):
        r = p.add_run(ln)
        r.font.name = "Consolas"
        r.font.size = Pt(size)
        if i < len(lines) - 1:
            r.add_break()
    # light box
    pPr = p._p.get_or_add_pPr()
    bdr = OxmlElement("w:pBdr")
    for side in ("top", "left", "bottom", "right"):
        e = OxmlElement("w:" + side)
        e.set(qn("w:val"), "single")
        e.set(qn("w:sz"), "4")
        e.set(qn("w:space"), "4")
        e.set(qn("w:color"), "BFBFBF")
        bdr.append(e)
    pPr.append(bdr)
    sh = OxmlElement("w:shd")
    sh.set(qn("w:val"), "clear")
    sh.set(qn("w:fill"), "F5F5F5")
    pPr.append(sh)
    return p


def table(headers, rows, widths, font=8.5, header_fill="1F4E79"):
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    hdr = t.rows[0].cells
    for i, hh in enumerate(headers):
        hdr[i].text = ""
        p = hdr[i].paragraphs[0]
        p.paragraph_format.space_after = Pt(2)
        r = p.add_run(hh)
        r.bold = True
        r.font.size = Pt(font)
        r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        shade(hdr[i], header_fill)
    for ri, row in enumerate(rows):
        cells = t.add_row().cells
        for i, val in enumerate(row):
            cells[i].text = ""
            p = cells[i].paragraphs[0]
            p.paragraph_format.space_after = Pt(2)
            _mixed_sized(p, str(val), font)
            if ri % 2 == 1:
                shade(cells[i], "F2F6FA")
    for row in t.rows:
        for i, w in enumerate(widths):
            row.cells[i].width = Cm(w)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)
    return t


def plain_table(rows, widths, font=9):
    t = doc.add_table(rows=0, cols=len(widths))
    t.style = "Table Grid"
    t.autofit = False
    for ri, row in enumerate(rows):
        cells = t.add_row().cells
        for i, val in enumerate(row):
            cells[i].text = ""
            pp = cells[i].paragraphs[0]
            pp.paragraph_format.space_after = Pt(2)
            _mixed_sized(pp, str(val), font)
            if i == 0:
                for r in pp.runs:
                    r.bold = True
            if ri % 2 == 1:
                shade(cells[i], "F2F6FA")
    for row in t.rows:
        for i, w in enumerate(widths):
            row.cells[i].width = Cm(w)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)
    return t


def _mixed_sized(p, text, size):
    _emit(p, text, size=size)


def callout(label, text, fill="FFF4E5"):
    t = doc.add_table(rows=1, cols=1)
    t.style = "Table Grid"
    c = t.rows[0].cells[0]
    c.text = ""
    p = c.paragraphs[0]
    r = p.add_run(label + "  ")
    r.bold = True
    r.font.size = Pt(9)
    _mixed_sized(p, text, 9)
    shade(c, fill)
    c.width = Cm(16.6)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)
    return t


def toc_field():
    p = doc.add_paragraph()
    r = p.add_run()
    fld = OxmlElement("w:fldChar"); fld.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText"); instr.set(qn("xml:space"), "preserve")
    instr.text = r'TOC \o "1-2" \h \z \u'
    sep = OxmlElement("w:fldChar"); sep.set(qn("w:fldCharType"), "separate")
    tx = OxmlElement("w:t"); tx.text = "Right-click > Update Field to build the table of contents."
    end = OxmlElement("w:fldChar"); end.set(qn("w:fldCharType"), "end")
    for e in (fld, instr, sep, tx, end):
        r._r.append(e)
    return p


def footer_pagenum():
    ftr = doc.sections[0].footer
    p = ftr.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("MARMITES x CdL merge cookbook  |  ")
    r.font.size = Pt(8); r.font.color.rgb = GREY
    r2 = p.add_run()
    f1 = OxmlElement("w:fldChar"); f1.set(qn("w:fldCharType"), "begin")
    it = OxmlElement("w:instrText"); it.set(qn("xml:space"), "preserve"); it.text = "PAGE"
    f2 = OxmlElement("w:fldChar"); f2.set(qn("w:fldCharType"), "end")
    for e in (f1, it, f2):
        r2._r.append(e)
    r2.font.size = Pt(8); r2.font.color.rgb = GREY


def pagebreak():
    doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)


footer_pagenum()

# =============================================================== TITLE PAGE
p = para("MARMITES × CdL", bold=True, size=30, align=WD_ALIGN_PARAGRAPH.LEFT,
         color=ACCENT, space_after=0)
para("Merge implementation cookbook", bold=True, size=19, color=ACCENT, space_after=10)
para("Bringing the CdL modelling approach (SFR / LAK / CRR / MVR, config-driven build, "
     "PEST++-IES calibration) into MARMITES–MODFLOW 6, with a three-source total ET and "
     "a soil-moisture-aware calibration.", size=11.5, italic=True, space_after=16)

plain_table(
    [["Prepared for", "Alain Pascal Francés — LNEG, Portugal"],
     ["Date", datetime.date.today().isoformat()],
     ["Source A (model approach)", "`MF6models` — `CdL/code/cdl_gwf_model_fable_v2.py` (v2 “fable”), "
                                   "`postprocess_cdl.py`, the PEST chain. "
                                   "Local `E:\\00code\\MF6models`; remote github.com/AlainPascalFrances/MF6models"],
     ["Source B (host code base)", "`MARMITES`, branch `MM-MF6`. Local `E:\\00code\\MARMITES`; "
                                   "remote github.com/AlainPascalFrances/MARMITES"],
     ["Baseline state", "**Updated 2026-09-10.** WP0, WP1, WP1b (stage 1) and WP1c are DONE and committed on "
                        "branch `MM-MF6_SFR_LAK_CRR`; `pytest code/tests -q` = **353 passed / 6 skipped**. "
                        "A coupled La Mata run completes on a **Voronoi mesh**, mass-conserving, with SFR "
                        "routing over the shared-face topology. Voronoi is NOT yet the default. WP2 onwards "
                        "not started. **Chapter 11 says what to verify and how.**"],
     ["Companion documents", "`doc/MM-MF6/MARMITES_NEXT_STEPS.md` (canonical running state), "
                             "`doc/MM-MF6/MARMITES_SFR_LAK_CRR_analysis.md` §7–8 (design record)"]],
    [4.0, 12.6], font=9.5)

para("")
callout("How to read this:",
        "this is the plan; **chapter 11, at the end, is the status report -- what has actually been built, and the command that checks each claim.** "
        "the work is cut into twelve work packages -- the ten numbered **WP0 -> WP9**, plus **WP1b** (the Streamlit interface) and "
        "**WP1c** (one grid path, Voronoi by default). Each one states its objective, the files it "
        "touches, a numbered step list, the configuration keys it introduces, and the acceptance criteria that "
        "must be green before the next package starts. Sections 6–8 are the reference material the packages "
        "point at: the input-data inventory, the configuration reference and the PEST-IES design.",
        fill="EAF1F8")

pagebreak()
h(1, "Contents")
toc_field()


pagebreak()

# =============================================================== 1
h(1, "1  Purpose, scope and working method")

h(2, "1.1  What this cookbook delivers")
body("Four things were asked for, and they are interleaved rather than sequential — so the plan orders them "
     "by dependency, not by the order they were requested:")
numlist("**SFR, LAK and CRR in MARMITES**, plus an **ET calibration with PEST++-IES** whose observations are "
        "heads (MF), soil moisture at several depths (MM), streamflow (SFR) and total ET (MM+MF). → WP3, WP4, WP5, WP7.")
numlist("**Total ET recomposed from three sources**: ETsoil from MM, ETuzf from MF (the deep unsaturated zone), "
        "ETg from MM — which means MODFLOW's own groundwater ET must stay switched off. → WP2.")
numlist("**A configuration file replacing the command-line tags**, in the style of the CdL script, plus an "
        "inventory of the input data MM and the SFR/LAK approach need, split between the repository and an "
        "outside-the-repo data root. → WP0, WP1, §6, §7.")
numlist("**A Streamlit web interface** to browse the input data — the shapefiles included — to launch the model "
        "locally or on the server, and to present the output. → WP1b.")

body("Two further packages were added because leaving them out would force rework: **WP6** (post-processing and "
     "the observation exports, without which WP7 has nothing to calibrate against) and **WP8** (tests, "
     "documentation and the `MM-MF6` → `master` merge hygiene).")
body("Two more follow from the decisions of 2026-09-09 (§1.4 and §3.3): **WP1c** — one grid path, **Voronoi by "
     "default**, so La Mata and CdL share a single discretisation and the surface packages are built once — and "
     "**WP9**, the port to CdL. Eleven packages in all.")

h(2, "1.2  Dependency order — why WP0 and WP1 come first")
body("WP2 to WP5 each introduce a handful of new settings and a handful of new input files. Adding them to the "
     "current `argparse` flag list and to an ad-hoc dataset layout, and only then refactoring to a configuration "
     "file, means writing the same settings twice and re-testing them twice. Doing the configuration and the "
     "data layout first costs about two days and is repaid by every package after it.")
body("**WP1b (the Streamlit interface) sits on both of them** — it edits the configuration file WP0 defines and "
     "browses the data WP1 sorts, so it cannot sensibly be built before either. It is then delivered in three "
     "stages rather than in one, because its output side has to follow the figures WP3–WP7 produce.")

body("**WP1c comes before WP2–WP5, not after.** SFR routing, LAK host cells and the CRR neighbour graph are all "
     "grid-topology dependent: build them on the raster and they are built twice.")

code([
    "WP0 config file    ─┬─► WP1b Streamlit UI (stage 1)",
    "WP1 data inventory  ├─► WP1c ONE GRID PATH — Voronoi by default",
    "                    │        │",
    "                    │        ├─► WP2  ET (ETsoil + ETuzf + ETg)",
    "                    │        ├─► WP3  SFR coupled + validated",
    "                    │        │        └─► WP4  LAK coupled + validated",
    "                    │        │                 └─► WP5  CRR (python cascade)",
    "                    └────────┴───────────────────► WP6  post-proc + obs exports",
    "                                                    │      └─► WP1b stage 2",
    "                                                    └─► WP7  PEST-IES",
    "                                                           ├─► WP1b stage 3",
    "                                                           └─► WP9  port to CdL",
    "",
    "                    WP8  tests / docs / merge   (runs alongside)",
])

h(2, "1.3  Working method (unchanged from the conversion phase)")
bullet("**Implement, then hand off.** The transient and multi-year runs are launched by Alain in Spyder. "
       "Every hand-off must name the files to reload in Spyder before the run, and the exact configuration to set.")
bullet("**One real short run per step.** Every package is validated against a real run "
       "(`run.nsp = 365`, one hydrological year) and the figures are looked at before the step is declared done. "
       "This is the discipline that ended the last restart; keep it.")
bullet("**Modify the native functions, never reimplement them.** That rule produced the complete figure suite; "
       "it applies equally to `MARMITESplot_v3.py` and to the CdL post-processing being ported in WP6.")
bullet("**`MARMITES_NEXT_STEPS.md` stays the canonical running state.** This cookbook is the plan; that file "
       "records what actually happened. Update its “CURRENT STATE” block at the end of every work package.")

h(2, "1.4  Study area: develop on La Mata, deliver on CdL")
body("**Decided 2026-09-09.** The run-time difference — about 10 minutes at La Mata against about 30 at CdL — is "
     "the least important factor. The decisive one is data.")
callout("The blocker:",
        "**CdL's forcing is monthly** (`p_month_198101_202605.csv`, `et0_month_198101_202605.csv`). MARMITES is a "
        "daily model: MMsurf needs hourly meteorology — rainfall, wind speed, air humidity, temperature, incoming "
        "solar radiation — to produce the daily P / Pe / RF / RFe / TF / PT / PE / E0 / LAI series MMsoil "
        "consumes. **MM cannot run on CdL today at all.** It is blocked on data, not merely slower.")
body("Three further gaps, none of them small: CdL has no soil or vegetation zoning in the MARMITES sense (the raw "
     "material exists — the `ks` / `ths` / `wp` / `fc` rasters and COS-2025 — but must be converted); MARMITES "
     "builds its grid from ASCII rasters whereas CdL's comes from `voronoi_grid.pkl` + `voronoi_layers.npz` "
     "(WP1c removes this one); and CdL's observations are not in MARMITES formats.")
body("The calibration argument also runs the other way. La Mata being uncalibrated is not a drawback here — "
     "**calibration is the deliverable** — and La Mata carries the better observation set for the ET and "
     "soil-moisture calibration this project is about: real head records at 11 points, real soil-moisture probes "
     "at several depths, catchment runoff, and an eddy-covariance tower. CdL's counterparts are **synthetic** "
     "piezometers and MODIS AET. (Worth confirming on the CdL side: per the project record the last IES "
     "overfitted to a checkerboard K field and the DRN-datum change invalidated that parameter set, so CdL has "
     "calibration machinery rather than a current trusted calibration.)")
body("**So: WP0–WP7 on La Mata, including a real calibration; CdL enters through WP9.** The configuration file "
     "and the `example/<case>/` layout make the port data and a grid reader, not a rewrite. One rule protects "
     "that: nothing case-specific in `code/` — if it names La Mata, it belongs in a configuration file or in "
     "`example/LaMata/`.")


# =============================================================== 2
h(1, "2  Where the two code bases stand today")

h(2, "2.1  Side-by-side")
table(
    ["Aspect", "MARMITES (branch MM-MF6)", "CdL (MF6models)"],
    [["Purpose", "Physically-based soil-water balance (MMsurf + MMsoil) coupled to MF6 through the BMI/API; "
                 "La Mata catchment, daily, 1949 stress periods",
      "MF6 groundwater-flow model of the Casa de Lobos catchment; monthly 1981–2026"],
     ["Grid", "Structured DIS (and DISV-from-DIS) inherited from the NWT model; 50 m cells. An **unstructured "
              "path already exists but is unwired** — `VertexGeometry`, `RefinedModel`, `resample_inputs` — see "
              "WP1c, which makes **Voronoi the default**",
      "Voronoi/DISV via Triangle; ncpl 6609 (design 2) or 3295 (design 1)"],
     ["Layers", "`--nlay 2` reads `_2s1L.ini` directly (6-layer set still available)",
      "nlay = 3 — one numerical layer per geologic unit (U1/U2/U3), Daoud-style"],
     ["Soil / unsaturated", "**MMsoil** does the soil column in Python; UZF6 is only the percolation zone "
                            "**below** the soil, with `simulate_et=False`",
      "**UZF is the soil zone** — per-cell soil rasters, `simulate_et` + `linear_gwet` on"],
     ["Seepage face", "DRN-SEEP (`seep='drn'`, cond 10000) — validated; exfiltration returned to the MM soil column",
      "DRN-SEEP with cubic smoothing over DDRN; routed on to SFR by MVR"],
     ["Streams", "`marmites_sfr.py` — priority-flood routing from `inputPONDw.asc`; **builds and reloads, "
                 "never run coupled**", "SFR built directly on the DISV grid; DEM routing, Manning rectangular, "
                                        "drainage-scaled width, downstream-monotonic bed (SFRmaker rule)"],
     ["Ponds", "`marmites_lak.py` — 12 EMBEDDEDV lakes from `lm_ponds.shp`; **builds, never run coupled**",
      "19 charcas as EMBEDDEDV LAK, centroid-seeded cells, wedge stage/vol table, Manning spill outlet"],
     ["CRR", "**Not started.** Design decided (analysis doc R2): must run in Python, receiver = MM soil column",
      "Implemented as static MVR movers, MFD weights α=β·S/ΣS, receivers LAK > SFR > UZF"],
     ["ET", "ETsoil + ETg both in MM; ETg applied to MF6 as a WEL sink; **no UZF ET at all**",
      "All ET in UZF (unsaturated + groundwater), extinction depth = COS rooting depth"],
     ["Coupling", "BMI/API through libmf6, `lagged` or `iterative`, per stress period",
      "None — a standalone MF6 run"],
     ["Configuration", "~40 `argparse` flags on `code/tests/run_lamata_mf6.py`; a legacy positional `.ini`; "
                       "a TOML schema exists (`code/marmites_config.py`) but only covers the MM run-control block",
      "`config.py` (three machine paths) + a documented module-level settings block at the top of the script"],
     ["Calibration", "None. Drawdown deficit flagged as a calibration matter, deferred",
      "pyEMU + pestpp-ies: pilot points, geostatistical prior, Tikhonov smoothness, "
      "heads / MODIS AET / SFR streamflow"],
     ["Post-processing", "Native `MARMITESplot_v3` suite, complete: Sankeys, per-point series, calibration "
                         "criteria, MM and GW maps, input parameter maps, general map",
      "`postprocess_cdl.py`: water balance in mm/yr, streamflow hydrographs, head and AET fits, "
      "lake stage/volume, MVR accounting, per-pond panels"],
     ["Repo rule", "Repo = code + docs + **strictly** the files MM/MF read; all output to "
                   "`$MARMITES_WS_ROOT` outside the repo", "Repo carries `input_data/` and a ready-to-run "
                                                            "`input_model/`; large rasters documented, not committed"]],
    [2.5, 7.2, 6.9], font=8)

h(2, "2.2  What already exists on the MARMITES side (do not rebuild it)")
bullet("`code/ppMF6/marmites_sfr.py` — priority-flood stream network, downstream-monotonic bed elevations, "
       "`EVAPORATION = 0`, outlet reaches replacing the six outlet DRN cells.")
bullet("`code/ppMF6/marmites_lak.py` — EMBEDDEDV lakes, wedge stage–volume–area table, host-cell assignment by "
       "centroid, on-channel / inlet / outlet reach attributes.")
bullet("`code/ppMF6/marmites_mf6.py` — MVR wiring for the stream running through the on-channel ponds; "
       "`_add_sfr_package`, `_add_lak_package`, `_add_mvr_package`.")
bullet("`code/marmites_coupler.py` — the SFR **INFLOW** pointer is already bound and the MM runoff of channel "
       "cells is already injected as reach inflow, through the post-`prepare_solve` callback.")
bullet("`code/marmites_config.py` — a validated TOML schema with an `MMConfig` dataclass and a legacy-ini "
       "converter. WP0 extends this rather than inventing a new mechanism.")

h(2, "2.3  What is genuinely missing")
bullet("A coupled **run** of SFR and LAK — the packages have never been exercised through the API.")
bullet("Any ET in the deep unsaturated zone (UZF), and therefore any three-source ET total.")
bullet("The CRR cascade, and the descending-topographic cell ordering it needs.")
bullet("A LAK↔MM exchange: lake inflow written from MM, lake stage read back so MM's open-water evaporation "
       "uses the real wetted area.")
bullet("Streamflow, lake-stage and total-ET observation exports — the inputs a PEST run consumes.")
bullet("Everything on the calibration side.")


# =============================================================== 3
h(1, "3  Target architecture")

h(2, "3.1  The merged picture")
code([
    "  MARMITES (python, per stress period)                     MODFLOW 6 (libmf6, BMI)",
    "  ────────────────────────────────────                     ─────────────────────────",
    "  MMsurf   P, PE, PT, LAI, interception  ─┐",
    "                                          │",
    "  MMsoil   per soil layer:                │   perc (m/d)        ┌──►  UZF   SINF",
    "            Esoil, Tsoil, Ssoil, Rsoil    ├───────────────────► │      (deep unsat zone only)",
    "            Ssurf, I, Ro, Eow             │   PETres (m/d)      └──►  UZF   PET   ── ETuzf",
    "            ETg = Eg + Tg                 │   ETg (m3/d)        ────►  WEL   Q",
    "                                          │",
    "  CRR      MFD cascade over the DEM       │   Ro_reach (m3/d)   ────►  SFR   INFLOW",
    "            β·S_ij/ΣS_ij, sinks evaporate ├───────────────────►",
    "            receivers: LAK > SFR > soil   │   Ro_lake (m3/d)    ────►  LAK   RUNOFF",
    "                                          │",
    "                            ◄─────────────┘   heads, exfiltration (DRN-SEEP),",
    "                                              rejected infiltration, LAK stage",
    "",
    "  TOTAL ET  =  Ei (MM)  +  Eow (MM)  +  ETsoil (MM)  +  ETuzf (MF/UZF)  +  ETg (MM, via WEL)",
    "                                                         ▲",
    "                                       UZF groundwater ET is OFF — no LINEAR_GWET / SQUARE_GWET",
])

h(2, "3.2  What crosses the API each stress period (extended)")
table(["Direction", "Quantity", "Target / source", "Status"],
      [["MM → MF6", "percolation (m/d)", "UZF `SINF` (never `FINF`)", "exists"],
       ["MM → MF6", "ETg (m³/d)", "WEL `Q`, AUTO_FLOW_REDUCE", "exists"],
       ["MM → MF6", "runoff of channel cells (m³/d)", "SFR `INFLOW`", "exists, never run coupled"],
       ["**MM → MF6**", "**residual PET after MM's own ET (m/d)**", "**UZF `PET`**", "**new — WP2**"],
       ["**MM → MF6**", "**runoff captured by a pond (m³/d)**", "**LAK `RUNOFF`**", "**new — WP4**"],
       ["MF6 → MM", "heads (m)", "`X`, per surface cell", "exists"],
       ["MF6 → MM", "exfiltration (m³/d)", "DRN-SEEP flows (or UZF `GWD`)", "exists"],
       ["MF6 → MM", "rejected infiltration (m³/d)", "UZF `REJ-INF`", "exists"],
       ["**MF6 → MM**", "**ETuzf actually taken (m³/d)**", "**UZF ET budget term**", "**new — WP2**"],
       ["**MF6 → MM**", "**lake stage (m) → wetted area for E_ow**", "**LAK `STAGE`**", "**new — WP4**"]],
      [2.2, 5.2, 5.2, 4.0], font=8.5)

callout("Timing rule (do not relearn it):",
        "every API-written input — `SINF`, `Q`, `INFLOW`, and now `PET` and LAK `RUNOFF` — must be written "
        "**after** `prepare_solve`, through the existing `write_cb` callback. UZF re-derives its arrays from the "
        "period data during both `prepare_time_step` and `prepare_solve`, so anything written earlier is silently "
        "overwritten. This is the bug that made UZF apply a constant 977 m³/d for the whole project.")

h(2, "3.3  Three decisions taken here, with their reasons")
h(3, "D1 — CRR runs in Python, not as MVR movers")
body("The CdL model routes the cascade with static MVR movers whose reinfiltration receiver is a UZF cell. "
     "In CdL that is correct, because UZF **is** the soil zone. In MARMITES, UZF is only the percolation zone below "
     "the soil column, so an MVR reinfiltration would inject water beneath `Ssoil` and bypass `Esoil`, `Tsoil` and "
     "the Eq. 1b infiltration limit — defeating the ET partitioning that is the point of the model. The cascade "
     "therefore runs in the coupler, reusing Daoud's Eq. 23 weights and the LAK > SFR > soil receiver priority; "
     "only the reinfiltration target differs. (Analysis doc, revision R2.)")

h(3, "D2 — LAK connection type is EMBEDDEDV, definitively")
body("The ponds are sub-grid (341–2036 m² against a 2500 m² cell; two contain no cell centre), so there is nothing "
     "to excavate. A VERTICAL connection puts the lakebed at the cell **top**, which would leave a dug-in pond "
     "permanently dry. This matches the locked CdL decision; tune inside the design (bed leakance, table, surfdep, "
     "solver settings), never the connection type.")

h(3, "D3 — geospatial libraries stay off the model path")
body("Reading `hydrography.shp` and `lm_ponds.shp` during the build would put geopandas (or pyshp, which is not "
     "installed) into the model path, where today only one optional figure needs it. Instead, a one-off converter "
     "(`code/tools/gis_to_dataset.py`, WP1) turns the shapefiles and rasters into small plain-text tables that live in "
     "the repository dataset and that MM/MF read directly. This answers the open dependency question in the "
     "handoff document and satisfies the repository rule at the same time.")

h(3, "D4 — two infiltration returns, in opposite directions; do not merge them")
body("Two distinct fluxes are easy to conflate, and conflating them would quietly break the soil budget. Both "
     "exist, and they chain:")
table(["", "Rejected infiltration returning to the soil", "CRR cascade reinfiltration"],
      [["What it is", "MF6/UZF refuses infiltration, or groundwater exfiltrates; the water is given back to the "
                      "**same** cell's soil column",
        "Surface runoff `Ro` routed from an **upslope** cell to a downslope one"],
       ["Enters at", "**the BOTTOM** soil layer, propagating upward", "**the TOP** soil layer"],
       ["Governed by", "saturation from below: the column fills, emerges as `Exf_g1` and becomes Dunnian runoff",
        "the infiltration law, Eq. 1b: `I = min[Ssurf , D_soil,1 (φ1 − θ1)]`"],
       ["Status", "**already coded** in MMsoil — the documented MARMITES behaviour (analysis doc §5). WP5 must "
                  "not disturb it", "**WP5**, new; D1 above is about this leg only"]],
      [2.2, 7.2, 7.2], font=8.5)
body("The chain runs: MF6 rejects infiltration → it returns bottom-up into that cell's column → it emerges as "
     "`Ro` → CRR routes that `Ro` downslope → it reinfiltrates top-down into the receiving column, or reaches a "
     "pond or a reach. Where D1 and this agree, and it is the point that matters: the receiver is always the "
     "**MARMITES soil column**, never a UZF cell.")

h(3, "D5 — one grid path, Voronoi by default")
body("La Mata moves off its 50 m structured raster onto a Voronoi/DISV grid built the way CdL builds its own, so "
     "both case studies share one discretisation and the surface packages are written once. This is **WP1c**, "
     "placed before WP2 because SFR routing, LAK host cells and the CRR neighbour graph all depend on grid "
     "topology. It also removes D2's sub-grid problem: with pond-scale seeded cells a pond cell really is a pond, "
     "so WP4's open-water bypass becomes exact rather than approximate. The structured producer is **kept "
     "permanently**, as the regression anchor and for the MODFLOW-NWT comparison.")


# =============================================================== 4
h(1, "4  Standing rules and known traps")
body("These are the things that have already cost time. Every step below assumes them.")

table(["#", "Rule / trap", "Why it matters"],
      [["1", "The repository holds code, docs and **strictly** the files MM and MF read directly. All run output "
             "goes to `$MARMITES_WS_ROOT`.", "No GIS, no output, however small. `.gitignore` enforces it."],
       ["2", "Never `git add -A` while on `master`.", "`master` carries only upstream's minimal `.gitignore`."],
       ["3", "UZF infiltration is `SINF`, not `FINF`.", "`FINF` resolves to a valid but non-operative pointer."],
       ["4", "Write every API input **after** `prepare_solve`.", "See §3.2. Applies to the two new writes as well."],
       ["5", "`get_input_var_names()` is not authoritative.", "Advanced-package variables (UZF `SINF`) are "
             "reachable but unlisted. Validate a bound pointer by **size**."],
       ["6", "Keep the infiltration-fidelity guard in `check_solution`.", "It is what caught trap 3. WP2 adds a "
             "sibling guard for groundwater ET."],
       ["7", "A steady-state period ignores the initial-head file.", "Seed it with mean recharge/ETg "
             "(`--steady-means` / `spinup.steady_means`), not carried heads."],
       ["8", "`max_outer` must equal the IMS `OUTER_MAXIMUM`.", "Otherwise ATS never sees a failed step and "
             "cannot retry with a smaller dt."],
       ["9", "Confirm `.hds` / `.lst` / `.uzf.cbc` and `_coupled_*.h5` come from the **same** run.",
         "A one-year cbc mixed with a full-run h5 has already wasted a diagnosis."],
       ["10", "flopy: an EXTERNAL lake outlet is `lakeout = -1`, **not** 0.",
         "The CdL transient failed on exactly this."],
       ["11", "flopy `get_structured_faceflows(ia=, ja=)` is broken upstream (3.10 and 3.11).",
         "`_ja_down_index` in `marmites_postprocess.py` is the replacement. Remember it if SFR/LAK "
         "post-processing needs face flows."],
       ["12", "`inputPONDw.asc` / `inputPONDhmax.asc` hold the **stream** network, not ponds.",
         "A legacy misnomer. WP1 renames them; the `IN_007` / `IN_008` maps built from them retire with SFR."],
       ["13", "A configuration silently overridden by a cache is a lost run.",
         "The CdL grid-design cache did this. WP0 stamps the resolved configuration into the run folder and "
         "fails fast on a mismatch."],
       ["14", "Remind Alain to reload the edited `.py` files in Spyder before flipping any setting.",
         "Every hand-off, without exception."],
       ["15", "**The two drain families have different datums.** OUTLET drains sit at the **bottom of the "
              "layer** (`botm[lay, cell]`, the CdL rule). **DRN-SEEP sits at the LAND SURFACE** and must stay "
              "there.", "DRN-SEEP represents groundwater emerging where the water table reaches the ground; move "
              "its datum to the layer bottom and the seepage face stops working — silently, since the package "
              "still builds and the run still converges. Confirmed 2026-09-09."]],
      [0.8, 7.0, 8.8], font=8.5)


# =============================================================== 5  WORK PACKAGES
h(1, "5  The work packages")

# ---------------- WP0
h(2, "WP0 — A configuration file replaces the command-line tags   [DONE]")
rich([("Objective. ", "b"),
      ("One human-readable file per model configuration, in the CdL spirit: machine paths in one place, "
       "model and run settings in a documented, validated block. The ~40 argparse flags become "
       "configuration keys; the CLI keeps only a `--config` argument and an explicit override switch.", "")])

h(3, "Design")
body("Two layers, mirroring CdL but reusing what MARMITES already has:")
table(["Layer", "File", "Holds", "CdL equivalent"],
      [["Machine paths", "`code/mm_paths.py`", "`REPO`, `DATA_ROOT`, `WS_ROOT`, `LIBMF6`, `MF6_EXE`, "
        "`PESTPP_IES`, `PYTHON_EXE` — each overridable by an `MM_*` environment variable",
        "`CdL/code/config.py`"],
       ["Model + run settings", "`code/configs/lamata.toml` (and siblings per configuration)",
        "everything that is a flag today, grouped in sections, plus the new SFR/LAK/CRR/ET/PEST keys",
        "the module-level settings block at the top of `cdl_gwf_model_fable_v2.py`"]],
      [3.2, 3.6, 7.4, 2.6], font=8.5)

body("TOML rather than a Python module, because `code/marmites_config.py` already parses TOML into a validated "
     "`MMConfig` dataclass with a legacy-ini converter and a test file (`code/tests/test_config.py`). Extending that "
     "schema costs a fraction of writing a second mechanism, and it keeps the configuration data — so it can be "
     "copied into a run folder, diffed between runs, and templated by pyEMU in WP7. Machine paths stay in Python "
     "so environment variables and `Path` arithmetic work exactly as in CdL.")

h(3, "Steps")
table(["#", "Action", "Files", "Done when"],
      [["0.1", "Extend the `MMConfig` dataclass with the new sections: `[paths] [run] [grid] [layers] [uzf] "
               "[seep] [et] [sfr] [lak] [crr] [spinup] [postproc] [pest] [meta]`. Give every key a default equal "
               "to today's flag default, so an empty file reproduces today's behaviour.",
        "`code/marmites_config.py`", "`load_mm_config()` round-trips a full file; unknown keys **raise**"],
       ["0.2", "Write `code/mm_paths.py` with the seven roots and their `MM_*` environment-variable overrides, "
               "plus a `report_paths()` that prints what resolved.",
        "new file", "importing it on a clean machine prints every root and flags the missing ones"],
       ["0.3", "Author `code/configs/lamata.toml` — the full annotated reference configuration (§7), reproducing the "
               "canonical spun-up run exactly.",
        "new file", "the run below is byte-identical to the flag-driven one"],
       ["0.4", "Rewrite the CLI: `run_lamata_mf6.py --config code/configs/lamata.toml [--set section.key=value ...] "
               "[--run-tag TAG]`. Keep `--set` for one-off overrides; delete the rest of the flags. Every override "
               "is echoed at startup.",
        "`code/tests/run_lamata_mf6.py`", "`--help` fits on a screen; no behaviour depends on an unprinted default"],
       ["0.5", "**Provenance**: copy the resolved configuration (after `--set`) into the run folder as "
               "`_input/resolved_config.toml`, and write its hash into the coupled HDF5 attributes.",
        "`run_lamata_mf6.py`, `marmites_coupler.py`", "every run folder says exactly what produced it"],
       ["0.6", "**Fail-fast guards**, in the spirit of the CdL design-tag sidecar: a cached grid / layer set / "
               "spin-up state carries the configuration hash it was built from; a mismatch stops with a message "
               "naming the offending key. No silent override, ever.",
        "`marmites_mf6.py`, `run_lamata_mf6.py`", "changing `layers.nlay` against a stale cache raises"],
       ["0.7", "Migration helper `code/tools/flags_to_toml.py`: takes an old command line, prints the equivalent TOML.",
        "new file", "the two canonical commands in `MARMITES_NEXT_STEPS.md` §0 convert cleanly"],
       ["0.8", "Extend `code/tests/test_config.py`: defaults, round-trip, unknown-key rejection, `--set` override, "
               "guard trip.", "`code/tests/test_config.py`", "suite green"]],
      [0.8, 8.0, 3.6, 4.2], font=8)

h(3, "Acceptance")
bullet("`python code/tests/run_lamata_mf6.py --config code/configs/lamata.toml` reproduces the current canonical run "
       "(same MF6 input files, same water balance).")
bullet("No model behaviour is reachable except through the configuration file or an echoed `--set`.")
bullet("Full flag → key map published (Appendix A) and the old commands in the handoff document rewritten.")

# ---------------- WP1
h(2, "WP1 — Input-data inventory and the DATA_ROOT split   [DONE]")
rich([("Objective. ", "b"),
      ("Name every input MM and the SFR/LAK/CRR approach need; put the raw geospatial sources in a data root "
       "outside the repository; and generate, once, the small plain-text intermediates that MM and MF actually "
       "read — which are the only new files allowed into the repository.", "")])

body("The full inventory is §6. The mechanism is a single converter run by hand, never by the model:")
code([
    "  DATA_ROOT (outside the repo)                     REPO  example/LaMata/  (read directly by MM/MF)",
    "  ────────────────────────────                     ──────────────────────────────────────────────",
    "  GIS/hydrography.shp        ─┐",
    "  GIS/lm_ponds.shp            │   tools/gis_to_dataset.py     inputSTREAM.csv     (reach table)",
    "  GIS/Limite.shp              ├──────────────────────────►    inputPONDS.csv      (pond table)",
    "  GIS/demfill, GIS/lm_dem     │   geopandas + rasterio,       inputSTREAMw.asc    (renamed)",
    "  GIS/202109MonitPts.shp     ─┘   run ONCE, by hand           inputSTREAMhmax.asc (renamed)",
    "                                                              inputDEMfill.asc    (CRR + SFR routing)",
])

h(3, "Steps")
table(["#", "Action", "Files", "Done when"],
      [["1.1", "Write `code/tools/gis_to_dataset.py`: reads the shapefiles/rasters named in `[paths] data_root`, "
               "writes the six intermediates above, prints a provenance header (source file, mtime, CRS, feature "
               "count) into each output.",
        "new file", "re-running it is idempotent; outputs carry their provenance"],
       ["1.2", "Rename the misnomer rasters `inputPONDw.asc` → `inputSTREAMw.asc` and `inputPONDhmax.asc` → "
               "`inputSTREAMhmax.asc`, with a compatibility shim that reads either name for one release.",
        "`example/LaMata/`, `marmites_sfr.py`, `marmites_postprocess.py`",
        "no code refers to a “pond” file that holds a stream"],
       ["1.3", "Retire the `IN_007_gridSsurfhmax` / `IN_008_gridSsurfw` input maps once SFR reads the reach table; "
               "replace them with an SFR-network map and a LAK-configuration map (a second mode of "
               "`_fig_general_map`, as flagged in the handoff).",
        "`marmites_postprocess.py`", "the input map set describes what the model actually uses"],
       ["1.4", "Point `[paths] nwt_reference` at the archive `E:\\00code_ws\\LaMata_new_PhD_artigo_2s3L\\_h5_MM.h5` "
               "instead of the repository dataset; keep the digest cache.",
        "`code/tests/plot_water_budget.py`", "the NWT comparison figures draw without a 1.3 GB copy in the repo"],
       ["1.5", "Write `example/LaMata/README.md` — the Tier A table of §6.1, one line per file: what reads it, "
               "what produced it, units.", "new file", "a fresh clone can be understood without archaeology"],
       ["1.6", "Add a `pytest` check that no file matching the excluded patterns (`*.shp`, `*.tif`, `*.h5`, "
               "output folders) exists under the repository.",
        "`code/tests/test_repo_hygiene.py`", "rule 1 is enforced by the suite, not by memory"],
       ["1.7", "**Enumerate every parameter of the three legacy ini files** (§6.4) into the WP0 schema — that is "
               "the authoritative definition of “the inputs MM needs”. `__inputMM_v3.ini` and "
               "`__inputMMsurf.ini` migrate parameter for parameter; `__inputMF_flopy_v3_2s3L.ini` is **obsolete** "
               "(MODFLOW-NWT) and is re-expressed for MF6, not carried across.",
        "`code/marmites_config.py`, `code/tools/flags_to_toml.py`",
        "every ini parameter has a schema key or is explicitly recorded as retired"],
       ["1.8", "Hand the enumeration to WP1b so the **front-end presents all of it already filled** from the "
               "current case study — nothing left for the user to discover.",
        "`code/app/pages/2_Configuration.py`", "opening the app on `example/LaMata` shows every parameter "
        "populated"]],
      [0.8, 8.0, 3.6, 4.2], font=8)

h(3, "Acceptance")
bullet("A fresh clone plus `DATA_ROOT` reproduces a run; nothing geospatial is imported on the model path.")
bullet("Every new repository file is one MM or MF reads directly, and is listed in the dataset README.")


# ---------------- WP1b
h(2, "WP1b — Streamlit web interface   [DONE, stage 1]")
rich([("Objective. ", "b"),
      ("A browser front-end that browses the input data — the shapefiles included — edits and validates the run "
       "configuration, launches a run on this machine or on the server, follows it while it runs, and presents "
       "the output. It is a **viewer and a launcher**. It never becomes part of the model path.", "")])

h(3, "Why it belongs here, and why it is not one delivery")
body("It sits on WP0 and WP1: it edits the configuration file WP0 defines and browses the data WP1 sorts into "
     "tiers. It is also the natural home for the geospatial libraries — decision **D3** pushes geopandas and "
     "rasterio off the model path, and this is the one place they are welcome, alongside the WP1 converter. "
     "But its output side can only show figures that exist, so it is delivered in three stages:")
table(["Stage", "Delivered with", "Content"],
      [["1", "WP1", "Inputs page (Tier A tables + a real map of the Tier B shapefiles and rasters), "
                    "Configuration page, Run page (launch + follow), a first Results page over the figures the "
                    "native suite already writes"],
       ["2", "WP6", "Results grows with the new material: streamflow hydrographs, per-pond lake stage and volume, "
                    "CRR routing, the ET split, soil moisture at depth against the probes"],
       ["3", "WP7", "Calibration page: phi progress, parameter and observation ensembles, per-group fits"]],
      [1.4, 2.6, 12.6], font=8.5)

h(3, "Layout")
code([
    "  app/                                REPO -- code only, as everything else in it",
    "    Home.py                           entry point; the run registry and the machine banner",
    "    pages/1_Inputs.py                 Tier A tables + Tier B shapefile / raster map",
    "    pages/2_Configuration.py          read, edit and VALIDATE configs/*.toml",
    "    pages/3_Run.py                    launch (local | server) and follow the log",
    "    pages/4_Results.py                figures, water balance, hydrographs, obs fits",
    "    pages/5_Calibration.py            pestpp-ies output                       (stage 3)",
    "    lib/loaders.py                    cached readers, keyed on (path, mtime, size)",
    "    lib/runs.py                       the run registry: launch, status, log tail",
    "  .streamlit/config.toml              theme, server address and port, upload cap",
])

h(3, "Steps")
table(["#", "Action", "Files", "Done when"],
      [["1b.1", "Scaffold the multi-page app and pin its dependencies in a **separate optional group** "
                "(`streamlit`, `streamlit-folium` or `pydeck`, `geopandas`, `rasterio`) so the model environment "
                "is untouched.", "`code/app/`, `environment-ui.yml`",
        "`streamlit run code/app/Home.py` opens; the model env still runs headless without any of it"],
       ["1b.2", "**Inputs page.** Tier A: each file as a table or a raster preview, with the units and the "
                "“read by” column of §6.1. Tier B: the shapefiles (hydrography, ponds, boundary, monitoring "
                "points, irrigated fields, EC tower) and the DEM on one map, layer toggles, click for attributes.",
        "`code/app/pages/1_Inputs.py`, `code/app/lib/loaders.py`",
        "every file listed in §6 is reachable from the page"],
       ["1b.3", "**CRS handling.** La Mata's projected CRS must be detected and reprojected to EPSG:4326 for a web "
                "map — reuse the CdL `report_crs` / reprojection pattern rather than assuming a CRS. Show the "
                "detected CRS in the sidebar.", "`code/app/lib/loaders.py`",
        "layers overlay correctly; a layer with a missing `.prj` is reported, not silently misplaced"],
       ["1b.4", "**Caching.** `@st.cache_data` on every reader, keyed on path + mtime + size (the same digest "
                "idea as the cbc cache); `@st.cache_resource` for the grid and open HDF5 handles. Without this, "
                "Streamlit re-reads a 180 MB coupled HDF5 on **every** widget click — it reruns the whole script "
                "each interaction.", "`code/app/lib/loaders.py`", "a widget click is instant on a full run"],
       ["1b.5", "**Configuration page** built from the `MMConfig` schema, not hand-written widgets, so new keys "
                "appear by themselves. Save refuses on a validation failure — the same fail-fast rule as WP0.6. "
                "Writes to `code/configs/`; never mutates a run's `resolved_config.toml`.",
        "`code/app/pages/2_Configuration.py`", "an invalid value cannot be saved; a new schema key appears unaided"],
       ["1b.6", "**Run page — local or server.** Launch the entry point as a **detached** subprocess writing "
                "`<ws_root>/runs/<stamp>/run.log` and `status.json`; the page polls the log tail. Never run the "
                "model in-line: a coupled run is ~11.5 min and an IES run is hours, and an in-line call would "
                "block the app and die on a browser refresh.", "`code/app/lib/runs.py`, `code/app/pages/3_Run.py`",
        "a run survives closing and reopening the browser; the log follows live"],
       ["1b.7", "**Server mode.** Two honest options; take the first: run Streamlit **on** the server "
                "(`streamlit run code/app/Home.py --server.address 0.0.0.0 --server.port 8501`) and reach it over the "
                "LAN — the PEST workers already run there, and the app only needs to see `$MARMITES_WS_ROOT`. "
                "The alternative, submitting over SSH from a local app, adds a credential path for no gain.",
        "`.streamlit/config.toml`, docs", "the same app serves both machines with one setting"],
       ["1b.8", "**Results page — show, do not redraw.** Present the PNGs `run_postproc` already wrote in "
                "`out_<stamp>_<tag>/_output/`, with a run picker and a figure filter. Add interactivity only "
                "where it buys something: a time slider on the head maps, a per-pond stage selector, an "
                "observation-point picker. Same rule as the native suite — one source of truth for a figure.",
        "`code/app/pages/4_Results.py`", "any archived run can be opened and compared without re-running anything"],
       ["1b.9", "**Run the WP1 converter from the app**: an “update dataset from GIS” action over "
                "`code/tools/gis_to_dataset.py`, with a preview of what would change before it writes. Makes the "
                "shapefile → plain-text step visible and repeatable instead of folklore.",
        "`code/app/pages/1_Inputs.py`", "the converter runs from the UI and reports its diff"],
       ["1b.10", "**Keep it off the model path.** `code/app/` may import `code/`; `code/` may never import `code/app/` or "
                 "`streamlit`. Add the assertion to the hygiene test. Bind to localhost by default. The launch "
                 "action takes a **validated configuration path only** — never a free-text command box, since the "
                 "page executes a process.", "`code/tests/test_repo_hygiene.py`",
        "the test fails if any `code/` module imports streamlit"]],
      [1.3, 7.5, 3.6, 4.2], font=8)

h(3, "Acceptance")
bullet("`streamlit run code/app/Home.py` on this machine and on the server serves the same app from one configuration "
       "setting.")
bullet("Every input in §6 is inspectable from the browser, shapefiles on a map with their attributes.")
bullet("A run is launched, followed and reviewed without touching a terminal — and the model still runs headless "
       "from Spyder and from a PEST worker with Streamlit not installed at all.")

callout("The one thing that would spoil this:",
        "the app quietly becoming a dependency of the model — a loader that only exists in `code/app/lib`, a "
        "configuration default that only the UI applies. The model must stay fully runnable from Spyder and from "
        "a PEST worker on a machine where Streamlit is not installed. Step 1b.10 is what keeps that true, and it "
        "is enforced by a test rather than by intent.")


# ---------------- WP1c
h(2, "WP1c — One grid path: Voronoi by default   [DONE, 1c.1-1c.8]")
rich([("Objective. ", "b"),
      ("Move La Mata off its 50 m structured raster onto a Voronoi/DISV grid built the way CdL builds its own, so "
       "the two case studies share one discretisation and every surface package below is written once. The "
       "structured grid stays available as a producer, permanently.", "")])

callout("Status 2026-09-10 - all eight steps implemented.",
        "A coupled La Mata run completes on a Voronoi mesh (ncpl 989, 469 active cells, cumulative "
        "discrepancy 0.000 % with SFR). The order of 1c.1/1c.2 was swapped and 1c.7 grew; see §11.5. "
        "**Voronoi is not yet the default** - `[grid] kind` stays `\"structured\"` until rung (c) of the "
        "ladder is judged on an equilibrated multi-year run. Modules built: `code/marmites_meshes.py` "
        "(producers), `code/marmites_mesh.py` (projection + resampling), "
        "`code/ppMF6/marmites_topology.py` (shared faces), `code/ppMF6/marmites_rasterise.py` (display "
        "adapter), `code/tools/validate_grid.py` (the ladder), `code/app/pages/5_Grid.py` (the mesh view).")

h(3, "What already exists — more than expected")
body("A survey of the code before planning this found the unstructured machinery largely present, and one piece "
     "of it genuinely clever:")
table(["Layer", "State", "Note"],
      [["`VertexGeometry`", "**generic already**", "accepts arbitrary `vertices`/`cell2d` or a flopy `VertexGrid`. "
        "`disv_from_structured` is one **producer**, not an assumption"],
       ["Cell list on an unstructured grid", "**exists**", "`marmites_gridgen.refined_cell_list()` returns "
        "`(cid, cid, 0, icell2d)` — setting **`i := cid`, `j := 0`** so every legacy `array[i, j]` access resolves "
        "against `(ncell, 1)` column vectors. This is what lets raster-era MM code run on arbitrary cells"],
       ["Raster → cell resampling", "**exists**", "`resample_inputs` + `marmites_raster.sample_to_cells`"],
       ["`RefinedModel`", "**exists, tested**", "ties gridprops + active nodes + inputs + geometry together "
        "(`test_refined_layout.py`)"],
       ["MF6 DISV build", "**exists**", "`--grid disv`"],
       ["Coupler", "**ready**", "carries `grid_kind` and `ncpl`; cell area comes from `geom`"],
       ["Driver wiring", "**missing**", "`RefinedModel` appears only in tests and a standalone quadtree script — "
        "it was never connected to the run"],
       ["SFR routing", "**needs work**", "priority-flood over the fixed 8-neighbour raster stencil `_NB8`"],
       ["CRR neighbours", "**easier on Voronoi**", "CdL's shared-face `_e2c` edge→cells map is the reference"],
       ["**Post-processing + figures**", "**raster-shaped**", "~104 `nrow`/`ncol` references across "
        "`marmites_postprocess.py` and `MARMITESplot_v3.py`, no vertex branch. `plotLAYER` takes `V` shaped "
        "`(ndays, nlay, nrow, ncol)`. **This is the real cost of the package**"]],
      [3.4, 2.6, 10.6], font=8)

h(3, "Design — the grid becomes a pluggable producer")
code([
    "  one function -> DISV gridprops + active_nodes, selected by  [grid] kind",
    "",
    "    structured   today's raster re-expressed as DISV (disv_from_structured).",
    "                 KEPT PERMANENTLY: the regression anchor and the MODFLOW-NWT comparison.",
    "    voronoi      Triangle + flopy.utils.voronoi.VoronoiGrid, seeded as CdL seeds it:",
    "                 graded stream corridor, centroid-seeded pond cells, background size.",
    "                 CdL's section 3 transfers almost verbatim.          <-- THE DEFAULT",
    "    quadtree     the existing gridgen path, kept because it already works.",
    "",
    "  everything downstream consumes RefinedModel -> .cells  .geom  .inputs",
])

h(3, "Steps")
table(["#", "Action", "Files", "Done when"],
      [["1c.1", "Wire `RefinedModel` into the driver — the missing link. `[grid] kind` selects the producer; the "
                "cached grid carries the configuration hash (WP0.6) so a stale mesh cannot be reused silently.",
        "`code/tests/run_lamata_mf6.py`, `marmites_gridgen.py`",
        "a run completes on a `quadtree` grid, proving the path end to end before Voronoi is added"],
       ["1c.2", "Add the **voronoi** producer: Triangle + `VoronoiGrid`, seeded from the watershed boundary, the "
                "stream lines and the pond polygons (the Tier A tables of WP1). Port CdL §3, including the "
                "Daoud-style graded transitions and the centroid pond seeding.",
        "`code/ppMF6/marmites_voronoi.py` (new)", "`triangle.exe` is at `C:\\00MODFLOW\\win64\\triangle.exe`; the "
        "mesh builds and its ncpl is reported"],
       ["1c.3", "**Upgrade the resampling.** `sample_to_cells` samples at the cell **centre**; a 100 m Voronoi cell "
                "over a 50 m raster then throws away most of the data. Make it **area-weighted** for continuous "
                "fields and **majority-vote** for zone rasters, keeping centre sampling as an option.",
        "`code/marmites_raster.py`", "resampled means match the raster means within a stated tolerance"],
       ["1c.4", "Resample every La Mata input onto the mesh: soil zones, vegetation areas, soil thickness, ibound, "
                "elevation, hk, thickness, Sy/Ss, DRN and GHB.",
        "`marmites_gridgen.resample_inputs`", "the `IN_*` input maps redraw on the new grid and look right"],
       ["1c.5", "**Shared-face neighbour topology** from the `cell2d` edge lists (CdL's `_e2c` map). This replaces "
                "the fixed 8-neighbour stencil `_NB8` in the SFR routing and is what CRR consumes in WP5.",
        "`code/ppMF6/marmites_topology.py` (new), `marmites_sfr.py`",
        "every interior face has exactly two cells; the graph is connected"],
       ["1c.6", "Resolve observation points by **nearest centroid / point-in-polygon** instead of `(i, j)`.",
        "`marmites_postprocess.resolve_obs_cells`", "all 11 points resolve, and to the right cells"],
       ["1c.7", "**Plotting adapter — do not rewrite the figure suite.** Rasterise the per-cell result vectors "
                "onto a display raster and hand the native functions the `(ndays, nlay, nrow, ncol)` array they "
                "already expect. A thin adapter against ~104 raster references, and it keeps the whole validated "
                "suite including `add_real_coord_axes`. True-polygon maps come free in the Streamlit app "
                "(`PlotMapView`) and can be added natively later.",
        "`code/ppMF6/marmites_rasterise.py` (new)", "every figure redraws with no change to `MARMITESplot_v3.py`"],
       ["1c.8", "**Validation ladder**, extending `test_disv_equivalence.py`: (a) DIS vs DISV-from-DIS — must be "
                "identical, the existing Phase-4 baseline; (b) DISV-from-DIS vs a uniform-size Voronoi — water "
                "balance within a stated tolerance; (c) the production Voronoi against the validated 45-year DIS "
                "run — every difference attributed to discretisation, not to a bug.",
        "`code/tests/test_voronoi_grid.py` (new)", "the three rungs pass and (c) is written up"]],
      [1.3, 7.5, 3.6, 4.2], font=8)

h(3, "Acceptance")
bullet("A full coupled La Mata run completes on a Voronoi grid, mass-conserving, with the whole figure suite "
       "drawn and no change to `MARMITESplot_v3.py`.")
bullet("`[grid] kind = \"structured\"` still reproduces the pre-WP1c run exactly.")

callout("Accept this deliberately:",
        "La Mata's validated baseline — the RMSE, the saved `hi_spinup_*.asc` state, the MODFLOW-NWT comparison — "
        "is all on the structured grid. Voronoi invalidates the saved state, and the NWT comparison becomes a "
        "cross-grid comparison rather than like-for-like. That is precisely why the structured producer is kept "
        "rather than deleted: the regression stays runnable on the grid it was established on.")


# ---------------- WP2
h(2, "WP2 — Total ET from three sources")
rich([("Objective. ", "b"),
      ("Total ET = Ei + Eow + ETsoil (MM) + ETuzf (MF/UZF) + ETg (MM). Today MM supplies ETsoil and ETg and UZF "
       "does nothing (`simulate_et=False`), so the deep unsaturated zone evaporates nothing. Turn UZF's "
       "unsaturated ET on, keep its groundwater ET off, and make the demand accounting airtight so nothing is "
       "counted twice.", "")])

h(3, "2a  The MF6 mechanism — verified, not assumed")
body("MF6's UZF package lets ET be simulated in the UZF cell **without** simulating it in the GWF cell: switch "
     "`SIMULATE_ET` on, choose the unsaturated formulation (`UNSAT_ETWC` — the water-content form — or "
     "`UNSAT_ETAE`), and simply omit both `LINEAR_GWET` and `SQUARE_GWET`. Omitting them is what suppresses "
     "groundwater ET. This was confirmed against the flopy 3.10 UZF package definition shipped in the pinned "
     "environment, not from memory. The resulting call differs from CdL's in exactly one keyword:")
code([
    "  CdL (all ET in UZF)                    MARMITES (unsaturated ET only)",
    "  ───────────────────                    ──────────────────────────────",
    "  simulate_et = True                     simulate_et  = True",
    "  linear_gwet = True   ◄── ETg in MF     unsat_etwc   = True",
    "  unsat_etwc  = True                     # linear_gwet / square_gwet  OMITTED  ◄── ETg stays in MM",
])

h(3, "2b  The demand chain")
body("PET must be spent once. The proposed order per stress period, which follows the MARMITES formulation and "
     "keeps the residual explicit:")
table(["Order", "Consumer", "Demand it sees", "Where computed"],
      [["1", "Interception `Ei`", "PE / PT partition from MMsurf", "MM (MMsurf)"],
       ["2", "Open water `Eow`", "PE over the wetted area (SFR reaches, LAK ponds)", "MM (MMsoil)"],
       ["3", "Soil `ETsoil` = Σ(Esoil, Tsoil) over soil layers", "what remains of PE and PT after 1–2",
         "MM (MMsoil)"],
       ["4", "**Deep unsaturated `ETuzf`**", "**the residual PET after 1–3, written to UZF `PET` per land-surface "
             "object; `extdp` = the rooting depth from land surface; `extwc` = θr**", "**MF (UZF) — new**"],
       ["5", "Groundwater `ETg` = Eg + Tg", "the residual after 1–4 — computed from **ETuzf ACTUAL**, not from "
             "the PET written at step 4 (see below) — then limited by depth to water table",
         "MM (MMsoil), applied as WEL `Q`"]],
      [1.2, 4.4, 7.0, 4.0], font=8.5)

callout("The detail that makes step 5 correct — actual uptake, not demand:",
        "what we write to UZF is a **demand**; what UZF takes is bounded by water content — with `UNSAT_ETWC` it "
        "extracts at the PET rate only while θ > `extwc`, and nothing below it. In a semi-arid catchment the deep "
        "unsaturated zone is dry for much of the year, so actual uptake is often a small fraction of the demand. "
        "If step 5 subtracted the **written demand**, ETg would be starved precisely when the deep zone is dry — "
        "which is exactly when the deep-rooted oaks draw on the water table instead. The totals would still look "
        "plausible while the **sourcing** — the thing this model exists to compute — would be wrong. So: "
        "`PET for ETg = PET_total − Ei − Eow − ETsoil − ETuzf_actual`, read back from the UZF budget. This is what "
        "makes step 2.5 load-bearing rather than a reporting convenience.", fill="FDE9E7")

body("**Two consequences worth stating at the call site.** First, **simultaneity**: ETg reaches MF6 as a WEL sink "
     "inside the same solve in which UZF computes its own ET, so within a stress period the two are simultaneous, "
     "not sequential. Strict sequential depletion is therefore exact only in `iterative` mode, where each outer "
     "iteration can use the current ETuzf iterate; `lagged` mode uses the previous stress period's actual. "
     "Second, **depth overlap**: UZF's `extdp` is measured from land surface and its unsaturated column "
     "conceptually overlaps the depth range the MM soil column already occupies. What prevents double counting is "
     "**the demand chain, not the geometry** — UZF can only take what remains of PET. So `extdp` stays the full "
     "rooting depth and must **not** be shrunk to “below the soil column”; that would suppress ETuzf for no "
     "reason.")

callout("Ordering caveat, state it in the code:",
        "in `lagged` mode MM cannot know step 4's outcome before MF6 solves, so ETg at step 5 uses the "
        "**previous** stress period's ETuzf actual. At daily stress periods that lag is small, but it is an "
        "approximation and must be documented at the call site and covered by a test, not discovered later.")

h(3, "Steps")
table(["#", "Action", "Files", "Done when"],
      [["2.1", "Add `[et]` configuration: `uzf_et` (bool), `unsat_form` (`etwc`|`etae`), `extdp_source` "
               "(`uniform`|`veg_zone`|`raster`), `extdp_default`, `extwc_source`, `gwet_in_mf` (must be `false`, "
               "guarded).", "`marmites_config.py`, `code/configs/lamata.toml`", "schema accepts and validates"],
       ["2.2", "Build the UZF package with `simulate_et=True`, the chosen unsaturated form, and **no** "
               "`linear_gwet` / `square_gwet`. Populate `extdp`, `extwc`, `ha`, `hroot`, `rootact` per land-surface "
               "object. Take extinction depth per vegetation zone (`inputVEG*area.asc`), as CdL does per COS class.",
        "`code/ppMF6/marmites_mf6.py`", "the written `.uzf` shows `SIMULATE_ET` and no GWET keyword"],
       ["2.3", "Bind the UZF `PET` pointer in the coupler (same `_bind_first` size validation as `SINF`) and write "
               "the residual demand in the post-`prepare_solve` callback.",
        "`code/marmites_coupler.py`", "the log line `coupler: UZF PET bound to …` appears; PET varies per SP"],
       ["2.4", "Compute the residual demand in MMsoil: expose `PET_residual` after step 3, and accept an "
               "`etuzf_actual_prev` argument for step 5 — the **actual** uptake of the previous stress period, "
               "never the demand written at step 4.",
        "`code/MARMITESsoil/MARMITESsoil_v3.py`", "unit test on a single column reproduces the chain by hand, "
        "including the dry-deep-zone case where actual << demand"],
       ["2.5", "**Load-bearing, not reporting:** read ETuzf **actual** back from the UZF budget and feed it to "
               "step 5, as well as into `wb_ts` / `wb_map`. The code today harvests only `UZF-GWRCH` and "
               "`REJ-INF`, addressed by `(text, paknam2)`; add the unsaturated-ET record and **confirm its exact "
               "record text on the first run** — the same discipline that uncovered the missing 1954-cell seepage "
               "face when both drain packages wrote under the text `DRN`.",
        "`marmites_coupler.py`, `marmites_postprocess.py`",
        "ETuzf actual reaches MMsoil and the water balance; the record text is pinned by a test"],
       ["2.5b", "Emit a **runtime** PET-balance line per stress period — "
                "`PET_total − Σ(Ei, Eow, ETsoil, ETuzf, ETg)` — so starvation or over-extraction shows up in the "
                "log while a run is going, not in a test afterwards.",
        "`marmites_coupler.py`", "the line appears and stays at ~0 through a one-year run"],
       ["2.6", "**GWET guard**, the sibling of the infiltration-fidelity guard: fail the run if any groundwater-ET "
               "term in the UZF budget is non-zero. Cheap, and it makes double counting impossible to ship.",
        "`marmites_coupler.check_solution`", "a deliberately mis-set `linear_gwet=True` fails the run"],
       ["2.7", "New flux indices `iETuzf` and `iETtot`; add them to `INDEX_MM`, to `flxIndex_lst`, to the Sankey "
               "(a new UZF-ET arm) and to the catchment/obs time series.",
        "`code/marmites_indices.py`, `MARMITESplot_v3.py`, `marmites_postprocess.py`",
        "the Sankey closes with the new arm; MM-side mass balance still ~0 %"],
       ["2.8", "**PET-closure test**: per cell per stress period, "
               "`Ei + Eow + ETsoil + ETuzf + ETg ≤ PE + PT + tol`, and the run-total ET is reported split by source.",
        "`code/tests/test_et_partition.py` (new)", "test green on a one-year run"]],
      [0.8, 8.0, 3.6, 4.2], font=8)

h(3, "Acceptance")
bullet("One-year coupled run: total ET splits into the five sources, the split sums to the demand, and the "
       "UZF groundwater-ET budget term is exactly zero.")
bullet("Catchment-scale ET is compared against the MODFLOW-NWT reference and against the eddy-covariance record; "
       "any change from the current 0.7 % ETsoil agreement is explained by the new ETuzf arm, not by a shift in "
       "ETsoil.")


# ---------------- WP3
h(2, "WP3 — SFR: run coupled, then validate")
rich([("Objective. ", "b"),
      ("Take the existing, never-executed SFR build to a validated coupled run, sourcing the network from the "
       "hydrography shapefile rather than from the misnamed raster, and adopting the CdL reach rules.", "")])

callout("Attribution, so the docs say the right thing:",
        "the CdL SFR build is **not** Daoud methodology. Daoud et al. 2022 is the source for the grid transitions "
        "(WP1c) and for CRR (WP5). CdL's SFR reads the stream lines from its GeoPackage, derives the ID/toID "
        "routing **from the DEM**, and gives each reach **SFRmaker-style attributes** (Leaf et al. 2021): one "
        "reach per cell, Manning rectangular, drainage-scaled width, downstream-monotonic bed. That is the method "
        "adopted here.")

h(3, "3a  Reach parameters: one source pattern, four producers")
body("CdL computes width from drainage area and keeps Manning, streambed K and thickness uniform. This project "
     "generalises that: **every reach parameter declares where its value comes from**, and the same pattern "
     "serves LAK and CRR.")
table(["Producer", "Meaning", "Typical use"],
      [["`{value: 0.035}`", "one number for the whole network", "Manning n, streambed thickness — the CdL default"],
       ["`{column: \"rhk\"}`", "from an attribute column of the stream shapefile", "a surveyed or zoned streambed K"],
       ["`{raster: \"...\"}`", "sampled from a raster", "a parameter that varies continuously"],
       ["`{drainage: {a: .., b: ..}}`", "**computed from contributing drainage area**, `w = a·A^b` "
        "(the arbolate-sum equivalent)", "**channel width — CdL's method, kept as an option and the default "
        "for width**"]],
      [3.6, 6.6, 6.4], font=8.5)
callout("Where this is resolved — the point that keeps D3 intact:",
        "the **converter** (WP1.1) reads the shapefile, evaluates whichever producer each parameter declares, and "
        "bakes the per-reach number into `inputSTREAM.csv`. The **Streamlit page** (WP1b) is where the mapping is "
        "chosen and the converter is run. The model itself never opens a shapefile.")

h(3, "Steps")
table(["#", "Action", "Files", "Done when"],
      [["3.1", "Source the network from `inputSTREAM.csv` (WP1) instead of `inputSTREAMw.asc`. Routing is derived "
               "**from the DEM** on the WP1c shared-face topology, not on the 8-neighbour raster stencil.",
        "`marmites_sfr.py`", "reach count and connectivity match the shapefile"],
       ["3.2", "Adopt the CdL reach rules: one reach per cell, Manning rectangular section, "
               "**downstream-monotonic** bed elevations. `EVAPORATION = 0` stays — MM does E_ow.",
        "`marmites_sfr.py`", "no reach bed sits above the lowest bed upstream of it"],
       ["3.2b", "Implement the **parameter-source pattern** of §3a, with `drainage` (`w = a·A^b`) kept as the "
                "default producer for channel width, and `value` for Manning, streambed K and thickness. "
                "Resolution happens in the converter, never in the model.",
        "`code/tools/gis_to_dataset.py`, `marmites_sfr.py`, `marmites_config.py`",
        "a width from drainage area and a width from a shapefile column both produce a valid network"],
       ["3.3", "**Only the DRN cell that actually becomes the stream outlet is replaced** by an outlet reach — "
               "the others are kept. CdL removes exactly one node's DRN-OUT record (the one that becomes the SFR "
               "inlet) and leaves the second outlet its drain. Replacing them all would delete legitimate "
               "boundary drainage; keeping the replaced one would drain the same water twice.",
        "`marmites_mf6.py`", "one drain record removed, the rest intact; outlet discharge appears once in the "
        "budget"],
       ["3.3b", "**Outlet-drain elevation = the BOTTOM of the layer**, as in CdL — `elev = botm[lay, cell]`, not "
                "a fixed depth below the cell top (CdL: `drnw_records.append([(lay, c), float(botm[lay, c]), "
                "cond])`). **DRN-SEEP is a different package and stays at the land surface**: it represents "
                "seepage where the water table reaches the ground, so its datum must not move.",
        "`marmites_mf6.py`", "outlet-drain elevations equal the layer bottoms; DRN-SEEP unchanged"],
       ["3.4", "Add an SFR continuous-observation block (inlet and outlet, per reach of interest) writing "
               "`<model>.obs.sfr.csv`, exactly as CdL does — WP7 consumes it.",
        "`marmites_mf6.py`", "the csv exists after a run, one column per observation"],
       ["3.5", "**First coupled run** — one hydrological year, `sfr.enable = true`, `seep = drn`, from the saved "
               "spin-up state. Watch: `INFLOW` binds (advanced-package variable — trap 5); the written reach inflow "
               "equals the MM runoff of channel cells; MVR balance; no non-convergence; discrepancy under the guard.",
        "hand-off to Spyder", "0 failed time steps, mass balance < 0.1 %"],
       ["3.6", "Post-processing: streamflow hydrograph in mm/yr and m³/d at the outlet, against "
               "`inputObsRo_catchment.txt`; SFR budget figure.",
        "`marmites_postprocess.py`", "hydrograph drawn, obs overlaid"],
       ["3.7", "Full-period run and comparison of the catchment water balance with and without SFR.",
        "hand-off", "the difference is explained term by term"]],
      [0.8, 8.0, 3.6, 4.2], font=8)

callout("Known physical limitation, document it once:",
        "an SFR reach exchanges **directly** with the aquifer — MF6 has no unsaturated zone beneath streams "
        "(MODFLOW-NWT's SFR2+UZF1 did; MF6 did not carry it over), and MVR cannot restore it because streambed "
        "leakage is an intrinsic package↔GWF connection, not water the package offers. In La Mata the channel "
        "follows the valley floor where the water table is shallow, so the bypass is acceptable — but say so in "
        "the model description rather than leaving it implicit.")

# ---------------- WP4
h(2, "WP4 — LAK: ponds coupled, with the MM open-water exchange")
rich([("Objective. ", "b"),
      ("Run the twelve EMBEDDEDV ponds coupled, wire the two missing MM↔LAK exchanges, and give the pre-processing "
       "the LAK configuration map that was asked for.", "")])

h(3, "4a  Open water bypasses the MM soil column")
body("**User decision, 2026-09-09:** wherever there is a LAK or SFR cell, the surface stands directly on UZF — "
     "the MM soil column is bypassed, so no ETsoil and no ETg there, only `Eow`, evaporation from open water.")
body("On the **Voronoi grid of WP1c this is applied wholesale, and it is exact**, because the mesh is built so "
     "that the cell **is** the feature: every pond is centroid-seeded into its own pond-scale cell, so a LAK cell "
     "is a pond and nothing else. This is the case CdL was designed around, and it is what makes the decision "
     "clean rather than approximate — one more reason the grid change comes first.")
body("Stream cells are the single place where the mesh cannot follow the feature: a channel is metres wide while "
     "a corridor cell is tens of metres, and meshing to channel width is not practical at catchment scale. The "
     "descriptor therefore carries a share, which is 1 for a pond cell by construction and the channel's share "
     "of the cell for a reach cell:")
code([
    "  per cell:   f_lake + f_stream + f_soil = 1",
    "",
    "    LAK cell   pond-scale, centroid-seeded    f_lake   = 1        wholesale bypass, exact",
    "    SFR cell   corridor cell, channel in it   f_stream = w*L / A  the channel's share",
    "    elsewhere                                 f_soil   = 1        the MM soil column, unchanged",
    "",
    "  over f_lake + f_stream :  no ETsoil, no ETg, only Eow;  water goes to LAK / SFR",
    "  over f_soil            :  the MM soil column, exactly as it behaves today",
])
bullet("With EMBEDDEDV the lake exchanges **directly** with GWF, so the UZF column beneath a pond receives only "
       "the soil fraction's percolation — which is zero in a pond cell. Worth a comment at the call site.")
bullet("“No ETg under open water” is not a loss of physics: groundwater discharge to the pond happens through the "
       "lakebed leakage term instead. Same flux, different name.")
bullet("On the **structured** producer — kept only for regression — ponds are sub-grid and `f_lake` would fall "
       "well below 1. The descriptor covers that case automatically; it is simply no longer the production one.")
callout("WP4 and WP5 must share ONE mechanism:",
        "both perturb the same thing — which cells have a soil column, and in what proportion. Build a single "
        "**per-cell surface descriptor** (`f_soil`, `f_lake`, `f_stream`, summing to 1) once at build time, with "
        "the cell-list order untouched. Two independent surgeries on the Phase-2 cell-list contract is how that "
        "contract gets broken.")

h(3, "Steps")
table(["#", "Action", "Files", "Done when"],
      [["4.1", "Source the ponds from `inputPONDS.csv` (WP1). Carry the CdL solver settings that stabilised the "
               "perched ponds: `surfdep` smoothing, `maxiter = 200`, `stagechg = 1e-4`; bed leakance from "
               "`[lak] bedleak` (1e-3 /d, the clay-lined charca value).",
        "`marmites_lak.py`, `marmites_mf6.py`", "package written, values traceable to the configuration"],
       ["4.2", "One EXTERNAL Manning spill outlet per lake — **`lakeout = -1`, not 0** (trap 10). Equilibrium "
               "initial stage from the steady-state water table.",
        "`marmites_lak.py`", "no lake starts dry; spill routes to the downstream reach"],
       ["4.3", "Keep the MVR wiring for on-channel ponds: inlet reach → LAK, LAK spill → downstream reach.",
        "`marmites_mf6.py`", "`diag`-style mover listing shows both legs"],
       ["4.4", "**New:** bind LAK `RUNOFF` and write the MM runoff captured by each pond, in the post-"
               "`prepare_solve` callback.", "`marmites_coupler.py`", "pond inflow varies with MM runoff"],
       ["4.5", "**New:** read LAK `STAGE` back and give MMsoil the wetted area, so open-water evaporation uses the "
               "real pond area rather than a fixed footprint. LAK's own evaporation stays at zero.",
        "`marmites_coupler.py`, `MARMITESsoil_v3.py`", "E_ow tracks stage; no double evaporation"],
       ["4.6", "Build the **per-cell surface descriptor** of §4a — `f_soil`, `f_lake`, `f_stream` — once at build "
               "time, on the cell list as ordered. MMsoil applies its column to `f_soil` only; `Eow` acts on "
               "`f_lake + f_stream`; ETsoil and ETg are suppressed over the open fraction. Shared with WP5.",
        "`marmites_grid.py` / the context builder, `MARMITESsoil_v3.py`",
        "fractions sum to 1 per cell; a pond-scale Voronoi cell gives `f_lake` ≈ 1; the catchment water balance "
        "still closes"],
       ["4.7", "LAK configuration input map — a second mode of `_fig_general_map`, plus per-pond zoom pages in the "
               "CdL style (grid cells, SFR cells, LAK cells, footprint, stream).",
        "`marmites_postprocess.py`", "`_input/IN_00x_lak_config.png` and the zoom pages exist"],
       ["4.8", "**Coupled run** — one year with `lak.enable = true` on top of WP3, then the full period.",
        "hand-off to Spyder", "stages physically coherent, LAK budget closes, no non-convergence"]],
      [0.8, 8.0, 3.6, 4.2], font=8)

callout("Expect low seepage and never-dry ponds.",
        "That was the CdL outcome at a bed leakance of 1e-3 /d and Alain confirmed it as the expected physics of "
        "clay-bottomed charcas. Do not treat it as a bug; it is the calibration lever, and in this project it "
        "becomes a PEST parameter (WP7).")


# ---------------- WP5
h(2, "WP5 — CRR: the cascade-routing and reinfiltration module")
rich([("Objective. ", "b"),
      ("Implement Daoud et al. (2022) Eq. 23 cascade routing with reinfiltration, in Python inside the coupler, "
       "with receiver priority LAK > SFR > downslope MARMITES soil column.", "")])

h(3, "The one structural risk")
body("CRR makes runoff non-local: an upslope cell must be solved before the cell that receives its runoff, so "
     "cells must be visited in **descending topographic order**. The cell list is in build order, and both the "
     "coupler and the DISV mapping depend on that order. The change must therefore be introduced as a "
     "**permutation computed once at build time** and stored on the context (`ctx.crr_order`), never as a re-sort "
     "inside the stress-period loop. The cell list itself does not move. It shares the per-cell surface "
     "descriptor built in WP4.6 — one mechanism, not two.")

h(3, "Steps")
table(["#", "Action", "Files", "Done when"],
      [["5.1", "Build the MFD weights once at build time on the **WP1c shared-face topology**: "
               "`α_ij = β · S_ij / Σ S_ij` over downslope neighbours; negative slopes drop out. Slopes from the "
               "sink-filled DEM sampled per cell.",
        "`code/ppMF6/marmites_crr.py` (new), `marmites_topology.py`", "weights sum to β per source cell"],
       ["5.2", "Compute and store `ctx.crr_order` — the descending-elevation permutation. Assert it is a "
               "permutation of the cell list and that every cell's receivers appear after it.",
        "`marmites_grid.py` / the context builder", "assertion holds on La Mata"],
       ["5.3", "Receiver resolution, priority **LAK > SFR > downslope soil column**; a cell with no downslope "
               "receiver is a topographic sink and its unrouted water evaporates (Daoud's convention). Count and "
               "report the sinks.", "`marmites_crr.py`", "sink count printed, matches the DEM"],
       ["5.4", "Reinfiltration enters the **topmost** MM soil layer, limited by Eq. 1b, "
               "`I = min[Ssurf , D_soil,1 (φ1 − θ1)]` — the same law as direct infiltration; the excess continues "
               "downslope. **Do not confuse this with the rejected-infiltration return**, which enters at the "
               "bottom and propagates upward and is already coded — see D4.",
        "`MARMITESsoil_v3.py`", "no cell exceeds its infiltration capacity; the bottom-up return still behaves "
        "exactly as before (regression test)"],
       ["5.5", "Deliver the residual reaching channel and pond cells to SFR `INFLOW` and LAK `RUNOFF` "
               "(WP3, WP4 wiring).", "`marmites_coupler.py`", "delivered volume equals the routed residual"],
       ["5.6", "Configuration: `[crr] enable`, `beta` (Daoud calibrates 0.8–1.0), `sink_behaviour`, "
               "`order_source`. `enable = false` restores today's behaviour exactly.",
        "`marmites_config.py`", "off-state is bit-identical to WP4"],
       ["5.7", "**Mass-balance test on a synthetic DEM** (a 3×3 and a 5×5 with a known sink and a known outlet), "
               "then on La Mata: routed in = routed out + reinfiltrated + evaporated at sinks.",
        "`code/tests/test_crr.py` (new)", "closure to machine precision on the synthetic case"],
       ["5.8", "Coupled run with CRR on; compare the outlet hydrograph and the catchment balance against WP4.",
        "hand-off", "the difference is attributable to reinfiltration, and quantified"]],
      [0.8, 8.0, 3.6, 4.2], font=8)

callout("Why not MVR:",
        "MVR moves water between MF6 **packages**, and MARMITES is not one. Routing the cascade through MVR would "
        "make the reinfiltration receiver a UZF cell — beneath the soil column — silently bypassing `Ssoil`, "
        "`Esoil`, `Tsoil` and Eq. 1b. Same physics and same parameters as CdL; only the receiver differs.")

# ---------------- WP6
h(2, "WP6 — Post-processing and the observation exports")
rich([("Objective. ", "b"),
      ("Extend the native suite to the new fluxes, port the CdL figures that have no MARMITES equivalent, and emit "
       "the four observation files WP7 consumes.", "")])

table(["#", "Action", "Files", "Done when"],
      [["6.1", "Native suite: add ETuzf and total ET to the Sankeys, the catchment and per-point time series and "
               "the flux maps; add SFR and LAK arms.",
        "`MARMITESplot_v3.py`, `marmites_postprocess.py`", "one `--postproc-only --preproc` emits the lot, no skips"],
       ["6.2", "Port from `postprocess_cdl.py`: streamflow hydrographs in mm/yr, lake stage and volume panels "
               "per pond, MVR accounting (per-mover budget), the water-balance graphs.",
        "`marmites_postprocess.py`", "per-pond panels drawn for all ponds"],
       ["6.3", "**Soil-moisture-at-depth figure and export** — `MM_S[:, :, iSsoil_pc_s]` per soil layer at each "
               "observation point, against `inputObsSM_*`. This is the new calibration target.",
        "`marmites_postprocess.py`", "one page per observation point, one curve per depth"],
       ["6.4", "**Observation exports** (the WP7 contract), written on every run: `obs_heads.csv`, `obs_sm.csv`, "
               "`obs_sfr.csv`, `obs_et.csv` — long format, period-end rows, deterministic row count.",
        "`marmites_postprocess.py`", "row counts equal nper regardless of ATS sub-stepping"],
       ["6.5", "Extend `plotCALIBCRIT` to the new groups (streamflow, ET) alongside heads and soil moisture.",
        "`MARMITESplot_v3.py`", "NSE / RMSE / RSR / r reported per group"],
       ["6.6", "**WP1b stage 2** — surface all of the above on the Results page: streamflow, per-pond stage and "
               "volume, CRR routing, the ET split, soil moisture at depth. Show the written figures; add a run "
               "picker and a run-to-run comparison.",
        "`code/app/pages/4_Results.py`", "a run is reviewed end to end from the browser"]],
      [0.8, 8.0, 3.6, 4.2], font=8)

callout("ATS determinism — the trap CdL already hit:",
        "adaptive time stepping gives a parameter-dependent number of rows in a time-series output, which breaks "
        "the positional PEST instruction files. Collapse every exported series to **stress-period-end rows** "
        "(CdL's `obs_collapse.py`) in the same way in the reference run and in every forward run.")


# ---------------- WP7
h(2, "WP7 — PEST++-IES calibration")
rich([("Objective. ", "b"),
      ("Calibrate the coupled MM+MF6 model against four observation groups — heads, soil moisture at depth, "
       "streamflow, total ET — with pyEMU and pestpp-ies, reusing the CdL calibration architecture and, more "
       "important, its hard-won failure modes.", "")])
body("The design is in §8. The steps here are the build order.")

table(["#", "Action", "Files", "Done when"],
      [["7.1", "Externalise the model inputs so pyEMU can template them: MF6 side via "
               "`sim.set_all_data_external()`; MM side by writing the soil/vegetation parameter maps and "
               "`inputSOILparam.txt` as flat, one-value-per-line files (pyEMU reads with `np.loadtxt`, which needs "
               "a rectangular file — CdL hit exactly this on wrapped DISV arrays).",
        "`code/tools/pest_prep_mm.py` (new)", "`org/` builds and MF6 runs from it unchanged"],
       ["7.2", "Build the PEST interface with `PstFrom`: pilot points for Kh per layer with a geostatistical "
               "structure, constants/multipliers for the boundary conditions, and the MM soil / vegetation / ET "
               "parameters (§8.1).", "`code/tools/build_pst_mm.py` (new)", "`mm.pst` builds; parameter count reported"],
       ["7.3", "Register the four observation groups from the WP6 exports, with the weighting scheme of §8.2.",
        "same", "`pst.obs_groups` = heads, sm, sfr, et"],
       ["7.4", "**Forward run** — the piece that differs most from CdL. `forward_run_mm.py` must: apply the "
               "multipliers; run the **coupled** MM+MF6 driver through libmf6 (not `mf6.exe`); collapse the "
               "outputs to period-end rows; write the four observation csvs. Guard it with "
               "`multiprocessing.freeze_support()` — without it the apply step silently kills the worker.",
        "`code/tools/forward_run_mm.py` (new)", "a base run (`noptmax = 0`) reproduces the reference phi"],
       ["7.5", "Draw the **geostatistical prior ensemble** (`prior_pe.jcb`) and build the first-order Tikhonov "
               "smoothness equations (`prior_cov.jcb`). Not optional — see the callout.",
        "`code/tools/build_pst_mm.py`", "prior drawn; PI equation count in the low thousands"],
       ["7.6", "Prior Monte-Carlo (`noptmax = -1`) to check the machinery and see the prior spread, then the base "
               "run, then the IES iterations.", "`code/tools/run_ies_mm.py` (new)", "prior MC completes; failure rate < 30 %"],
       ["7.7", "Post-process: phi progress, parameter and observation ensembles, per-group fits, and the "
               "soil-moisture-at-depth fit.", "`code/tools/postproc_ies_mm.py` (new)", "figures written to `ies_output/`"],
       ["7.8", "Extract the minimum-phi realisation into `pest_optimised.npz` and add a "
               "`[pest] use_pest_params` switch so a normal run can use the calibrated set.",
        "`code/tools/extract_pest_optimised_mm.py`, `marmites_config.py`", "a calibrated run reproduces the IES base"],
       ["7.9", "**WP1b stage 3** — the Calibration page: phi progress, parameter and observation ensembles, "
               "per-group fits, and the launch button for `run_ies_mm.py` on the server (detached, like any other "
               "run, but measured in hours).", "`code/app/pages/5_Calibration.py`",
        "an IES run is started and followed from the browser"]],
      [0.8, 8.0, 3.6, 4.2], font=8)

callout("Carry the CdL scars into the first build, not the second:",
        "the first CdL IES overfitted to a **checkerboard** K field — pestpp draws its prior from a diagonal, "
        "bounds-only covariance, so pilot points are spatially white. The fix package was: a geostatistical prior "
        "ensemble, first-order Tikhonov smoothness at `ies_reg_factor ≈ 0.5`, moderated bounds (wide enough to "
        "explore, tame enough that most realisations converge — ultra-wide bounds bought 60 % failures and a "
        "biased survivor set), `ies_bad_phi_sigma = 2.5`, and an absolute `overdue_giveup_minutes` cap. "
        "Build all of it in from the start.", fill="FDE9E7")

# ---------------- WP8
h(2, "WP8 — Tests, documentation and merge hygiene")
table(["#", "Action", "Done when"],
      [["8.1", "New test modules: `test_config.py` (extended), `test_et_partition.py`, `test_sfr_coupled.py`, "
               "`test_lak_coupled.py`, `test_crr.py`, `test_repo_hygiene.py` (repo rule **and** the "
               "“`code/` never imports streamlit” assertion).", "suite green, no new skips"],
       ["8.2", "Run the six `mf6`-binary tests with `E:\\00code\\moflowapi\\emsdatasets\\bin` on PATH — their "
               "fixture was repaired on 2026-09-09 but **has never been executed**.", "6 skips become 6 passes"],
       ["8.3", "Rewrite `MARMITES_NEXT_STEPS.md` §0-ter (flags → configuration) and §4 (SFR/LAK/CRR status); "
               "record the WP4.6 and WP2 ordering decisions in the analysis document.", "docs match the code"],
       ["8.4", "Clean the pre-existing cosmetics: unused locals and imports in `run_lamata_mf6.py` "
               "(`NCROP`, `check`, the re-imported `flopy`) and two in `marmites_postprocess.py`.", "lint clean"],
       ["8.5", "Plan the `MM-MF6` → `master` merge: **7 modify/delete conflicts on `trunk/pyEARTH1D/*`** (that path exists on `master`) "
               "(master modified them; the Phase-1 cleanup deleted them here). Decide whether pyEARTH1D stays on "
               "master before merging, not during.", "decision recorded; merge rehearsed on a scratch branch"],
       ["8.6", "Report the flopy `get_structured_faceflows(ia=, ja=)` defect upstream — it has been verified on "
               "3.10 and 3.11 and no issue is filed.", "issue opened"]],
      [0.8, 9.8, 6.0], font=8.5)


# ---------------- WP9
h(2, "WP9 — Port to CdL, the production case study")
rich([("Objective. ", "b"),
      ("Make Casa de Lobos the second case study of the same code, populating `example/CdL/` rather than writing "
       "a second model. WP1c has already removed the grid obstacle; what remains is mostly data.", "")])

callout("Step 9.1 is a data question, not a coding one — and it gates the rest:",
        "CdL carries **monthly** precipitation and ET0. MMsurf needs **hourly** meteorology (rainfall, wind speed, "
        "air relative humidity, temperature, incoming solar radiation) to produce the daily series MMsoil "
        "consumes. Until that record is found or assembled, MARMITES cannot run on CdL at any grid resolution. "
        "Answer this before scheduling the rest of the package.")

table(["#", "Action", "Done when"],
      [["9.1", "**Source the daily/hourly meteorology.** Options in order of preference: an existing station "
               "record for the catchment; a national network series (IPMA) disaggregated; or a reanalysis product "
               "downscaled. Whatever is chosen, MMsurf must be able to run from it and produce the "
               "`inputZON*_d.txt` series.", "MMsurf runs for CdL and writes the daily series"],
       ["9.2", "**Soil and vegetation zoning in the MARMITES sense.** Derive MM soil zones from the CdL soil "
               "hydraulic rasters (`ks`, `ths`, `wp`, `fc`) and vegetation areas from COS-2025 — the same rasters "
               "the CdL model already samples per cell. Write `inputSOILparam.txt` and the `.asc` zone maps.",
        "the `IN_*` parameter maps draw for CdL and are physically sensible"],
       ["9.3", "**Grid**: reuse CdL's own Voronoi mesh through the WP1c `voronoi` producer, or regenerate it from "
               "the same seeds. The `voronoi_layers.npz` layering (U1/U2/U3) maps onto the MM layer handling.",
        "the mesh loads and `ncpl` matches the CdL model"],
       ["9.4", "**Observations** into the MARMITES formats: `inputObs.txt` plus the per-point head, soil-moisture "
               "and streamflow files. Note that CdL's piezometers are synthetic/regionalised, so weight them "
               "accordingly in WP7's scheme.", "the observation exports populate for CdL"],
       ["9.5", "**SFR / LAK / CRR from the CdL geometry** — the streams and the 19 charcas, through the same "
               "converter and the same parameter-source pattern as La Mata. Nothing case-specific in `code/`.",
        "a CdL run builds with the surface packages on"],
       ["9.6", "**Run, validate, then calibrate**, reusing the WP7 chain unchanged. Compare against the existing "
               "MF6-only CdL model as a cross-check: the aquifer side should be recognisable.",
        "a coupled CdL run converges and its water balance is defensible"]],
      [0.9, 9.3, 6.4], font=8.5)

body("Estimate: dominated by 9.1 and 9.2, which are data-assembly tasks rather than programming ones. "
     "The code work (9.3, 9.5) should be small precisely because WP0, WP1 and WP1c were done first.")


# =============================================================== 6
h(1, "6  Input-data inventory")
body("Two tiers. **Tier A** lives in the repository because MM or MF reads it directly; **Tier B** lives in "
     "`DATA_ROOT` outside the repository and is consumed only by the one-off converter, by the pre-processing "
     "figures, or by the calibration.")

h(2, "6.1  Tier A — in the repository (example/<case>/, e.g. example/LaMata/)")
table(["Group", "File(s)", "Content / units", "Read by"],
      [["Run control", "`__inputMM_v3.ini` → `code/configs/lamata.toml`", "run and model settings (WP0)",
        "`run_lamata_mf6.py`"],
       ["Time", "`inputDATE.txt`", "the stress-period calendar", "MM, post-processing"],
       ["Meteo / vegetation series", "`inputZON*_d.txt`, `inputZON_*_stp.txt` (P, Pe, RF, RFe, TF, PT, PE, E0, "
                                     "Eo, LAI, crop)", "per zone, daily and per stress period, mm/d",
        "MMsurf → MMsoil"],
       ["Zoning maps", "`inputMETEOzones.asc`, `inputIRRzones.asc`, `inputSOILzones.asc`, "
                       "`inputVEG{1,2,3}area.asc`", "integer zone / fractional area rasters", "MM"],
       ["Soil parameters", "`inputSOILthick.asc`, `inputSOILparam.txt`", "soil thickness (m); per-zone θs, θr, "
                           "θfc, Ks, layer split", "MMsoil"],
       ["Aquifer grids", "`MF_ws/elev.asc`, `elev_sinkfil.ASC`, `ibound_l*.asc`, `hk_l*.asc`, `thick_l*.asc`, "
                         "`Sy_l*.asc`, `Ss_l*.asc`, `uzf_iuzfbnd.asc`", "per-layer aquifer properties",
        "`marmites_mf6.py`"],
       ["Boundary conditions", "`MF_ws/drn_cond_l*.asc`, `drn_elev_l*.asc`, `ghb_cond_l*.asc`, `ghb_head_l*.asc`",
        "DRN and GHB geometry and conductance", "`marmites_mf6.py`"],
       ["Saved state", "`MF_ws/hi_spinup_l*.asc`, `hi_spinup_perc.asc`, `hi_spinup_etg.asc`",
        "spin-up heads and mean recharge / ETg (baseline fallback)", "`marmites_mf6.py`, coupler"],
       ["MF6 parameter sets", "`MF_ws/__inputMF_flopy_v3_2s1L.ini`, `_2s3L.ini`", "2- and 6-layer parameter sets",
        "`marmites_mf6.py`"],
       ["Observations — heads", "`inputObs.txt`, `inputObsHEADS_{C1..C5,P0,W1}.txt`",
        "point name / X / Y / layer; date, head (m)", "post-processing, PEST"],
       ["Observations — soil moisture", "`inputObsSM_P0.txt`, `inputObsSM_SM.txt`",
        "date, θ (%) **per probe depth** — the WP7 `sm` group", "post-processing, PEST"],
       ["Observations — runoff", "`inputObsRo_catchment.txt`", "catchment runoff / streamflow",
        "post-processing, PEST"],
       ["**NEW — stream**", "**`inputSTREAM.csv`**", "**reach table: id, i/j (or node), toID, length, slope, "
                            "width, bed top, Manning n, streambed K and thickness**", "**`marmites_sfr.py`**"],
       ["**NEW — stream rasters**", "**`inputSTREAMw.asc`, `inputSTREAMhmax.asc`** (renamed from `inputPONDw` / "
                                    "`inputPONDhmax`)", "**channel width (m) and depth (m) per cell**",
        "**`marmites_sfr.py`, MM E_ow**"],
       ["**NEW — ponds**", "**`inputPONDS.csv`**", "**pond table: fid, centroid X/Y, area (m²), rim elevation, "
                           "bottom elevation, host cell, on-channel flag, inlet/outlet reach**",
        "**`marmites_lak.py`**"],
       ["**NEW — DEM**", "**`inputDEMfill.asc`**", "**sink-filled DEM (m) — CRR slopes and SFR routing**",
        "**`marmites_crr.py`, `marmites_sfr.py`**"],
       ["**NEW — grid seeds**", "**`inputWATERSHED.csv`**", "**catchment boundary polygon — the Voronoi domain "
        "(with the stream and pond tables above, these are the three mesh seeds)**",
        "**`marmites_voronoi.py` (WP1c)**"],
       ["**NEW — ET**", "**`inputEXTDP.asc`** (or a per-vegetation-zone table)",
        "**ET extinction / rooting depth (m) per cell — the UZF `extdp` column**", "**`marmites_mf6.py`**"],
       ["**NEW — ET observations**", "**`inputObsET_catchment.txt`** (or per zone)",
        "**observed actual ET (mm/d or mm/month) — the WP7 `et` group**", "**post-processing, PEST**"]],
      [2.6, 4.4, 5.6, 3.0], font=7.5)

h(2, "6.2  Tier B — outside the repository (DATA_ROOT)")
table(["Group", "Source (current location)", "Consumer", "Why it stays out"],
      [["Stream network", "`E:\\00code_ws\\LAMATA_new\\GIS\\hydrography.shp` (also `Streams.shp`, "
                          "`streams_dem.shp`, `stream_buffer.shp`)", "`code/tools/gis_to_dataset.py` → `inputSTREAM.csv`",
        "shapefile; needs geopandas"],
       ["Ponds", "`GIS\\lm_ponds.shp`", "→ `inputPONDS.csv`", "same"],
       ["Catchment boundary", "`GIS\\Limite.shp`", "the general map, the converter", "same"],
       ["DEM", "`GIS\\lm_dem`, `GIS\\demfill`, `GIS\\lm_streamord`", "→ `inputDEMfill.asc`; hillshade figures",
        "raster; needs rasterio"],
       ["Monitoring / observation points", "`GIS\\202109MonitPts.shp`, `202109ObsPts.shp`, `PilotsPoints.shp`",
        "the general map; PEST pilot-point layout", "shapefile"],
       ["Irrigated fields, EC tower", "`GIS\\Irr_Fields.shp`, `ECtower_footprint.shp`",
        "the general map; the ET observation footprint", "shapefile"],
       ["MMsurf inputs and outputs", "`E:\\00code_ws\\LaMata_MM-MF6\\MMsurf_ws`, and the provenance copy in "
                                     "`LaMata_new_PhD_artigo_2s3L\\MMsurf_ws`",
        "MMsurf (py3-ported, currently wired into nothing)", "run output / provenance"],
       ["MODFLOW-NWT reference", "`E:\\00code_ws\\LaMata_new_PhD_artigo_2s3L\\_h5_MM.h5` (1.3 GB), "
                                 "`MF_ws\\_h5_MF.h5` (3.1 GB)",
        "`plot_water_budget.load_reference()` — the NWT-vs-MF6 panels", "1.3 + 3.1 GB, cannot be regenerated"],
       ["All run output", "`$MARMITES_WS_ROOT` = `E:\\00code_ws\\LaMata_MM-MF6` "
                          "(`MF6_ws/`, `out_<stamp>_<tag>/`)", "everything", "rule 1"],
       ["PEST workspace", "`$MARMITES_WS_ROOT\\MM_pest` (`org/`, `template/`, `master/`, `workers/`)",
        "pyEMU, pestpp-ies", "generated; large; per machine"],
       ["Binaries", "`C:\\00MODFLOW\\mf6.7.0_win64\\bin\\libmf6.dll` and `mf6.exe`; the PEST++ suite",
        "the coupler, the calibration", "machine-specific — `mm_paths.py`"]],
      [2.6, 5.6, 4.6, 2.8], font=7.5)
body("Every shapefile and raster in this table is also a layer on the **WP1b Inputs page** — that page and the "
     "WP1 converter are the only two consumers of Tier B, and the only two places geopandas and rasterio are "
     "imported.")

h(2, "6.3  Data that does not exist yet and must be assembled")
table(["Item", "What is needed", "Candidate source"],
      [["Observed actual ET", "a catchment or zonal actual-ET series to weigh the `et` group against",
        "the eddy-covariance tower already in the GIS (`ECtower_footprint.shp`) — its flux record; failing that, "
        "MODIS MOD16A2GF zonal aggregation, exactly as CdL does (`download_mod16.py`, `make_mod16_zonal.py`)"],
       ["Streamflow at the outlet", "a gauged or regionalised series for the `sfr` group",
        "`inputObsRo_catchment.txt` if it is a discharge; otherwise the CdL donor-gauge regionalisation method"],
       ["Soil-moisture probe depths", "the depth of each probe, to map it onto an MM soil layer",
        "`inputObsSM_*` headers plus the field notes — must be pinned before WP7.3"],
       ["Streambed conductance", "reach-level streambed K and thickness",
        "not in the dataset today (the six DRN cells carry conductance 0.025–0.035). Start from the CdL value and "
        "make it a PEST parameter"],
       ["Pond bathymetry", "pond depth for the stage–volume table",
        "`POND_DEPTH` fallback (CdL uses 4.0 m below local ground); refine from survey if it exists"]],
      [3.0, 6.0, 7.6], font=8)

h(2, "6.4  The legacy ini files — the authoritative parameter list")
body("The tables above name the data **files**. The **parameters** are defined by the three ini files in "
     "`example/LaMata/`, and those are what the WP0 schema must cover key for key and the WP1b front-end must "
     "present already filled. Nothing may be quietly dropped in the migration.")
table(["File", "Blocks it defines", "Fate"],
      [["`__inputMM_v3.ini` (122 lines)",
        "**run control** (`run_name`, verbosity, `MMsurf_yn`, `MMsoil_yn`, `MF_yn`, `MF_lastrun`); "
        "**convergence loop** (`convcrit`, `convcritmax`, `ccnum`); **plotting** (`plt_out`, `plt_freq`, "
        "`nrangeMM/MF`, `ctrsMM/MF`, `ntick`, animation, `plt_out_obs`, `WBsankey_yn`, `plt_WB_unit`, "
        "`iniMonthHydroYear`, tick-grid years, `plt_input`); **workspaces and file names** (`MMsurf_ws`, "
        "`inputFile_PAR_fn`, `inputFile_TS_fn`, `MF_ws`, `MF_ini_fn`, `outputFILE_fn`, `MMsurf_fn`); "
        "**irrigation package** (`irr_yn`, `inputFile_TSirr_fn`, `gridIRR_fn`); **grid origin** (`xllcorner`, "
        "`yllcorner`); **the raster names** (`gridMETEO_fn`, `gridSOIL_fn`, `gridSOILthick_fn`, `gridPONDhmax`, "
        "`gridPONDw`); **soil parameters** (`SOILparam_fn`); **observations** (`inputObs_fn`, `inputObsHEADS_fn`, "
        "`inputObsSM_fn`, `inputObsRo_fn`, `rmseHEADSmax`, `rmseSMmax`); HDF5 compression",
        "**migrates key for key** into `[run] [postproc] [paths] [irr] [grid] [obs]`. `gridPONDhmax` / "
        "`gridPONDw` are the misnamed **stream** rasters (trap 12) and retire with SFR"],
       ["`__inputMMsurf.ini` (99 lines)",
        "**meteo zones** (`NMETEO`; per zone `phi`, `Lm`, `Z`, `Lz`, `FC`, `DTS`, `z_m`, `z_h`); "
        "**vegetation** (`NVEG`; per type `h_d`, `h_w`, `S_w`, `C_leaf_star`, `LAI_d`, `LAI_w`, `f_s_vd`, "
        "`f_s_vw`, `alfa_vd`, `alfa_vw`, `J_vd`, `J_vw`, `TRANS_vdw`, `TRANS_vwd`, **`Zr` max root depth**, and "
        "the transpiration-sourcing factors `kTg_min`, `kTg_max`, `kT_f`, `kT_1/s`); **crops** (`NCRP`, "
        "`NFIELD`, per crop `h_c`, `S_w_c`, `C_leaf_star_c`, `LAI_c`, `f_s_c`, `alfa_c`, `Zr_c`, its `kT*_c`); "
        "**soils** (`NSOIL`; per type `por`, `fc`, `alfa_sd`, `alfa_sw`, `J_sd`, `J_sw`, `TRANS_sdw`, "
        "`TRANS_swd`); **open water** (`alfa_w`)",
        "**migrates key for key** into `[mmsurf.meteo] [mmsurf.veg] [mmsurf.crop] [mmsurf.soil] "
        "[mmsurf.openwater]`. `Zr` is also the natural source for the WP2 `[et] extdp_source = \"veg_zone\"` "
        "rooting depth"],
       ["`__inputMF_flopy_v3_2s3L.ini`",
        "the MODFLOW-**NWT** model: `nlay/nrow/ncol/nper`, `itmuni`, `lenuni`, `delr`, `delc`, `elev_fn`, "
        "`ibound_fn`, `hnoflo`, the LPF/UPW blocks (`laytyp`, `layavg`, `chani`, `layvka`, `laywet`, `vka_fn`, "
        "`hdry`, `iphdry`…), OC and budget extensions",
        "**OBSOLETE — do not migrate.** Written for MF-NWT. Its *content* (layer geometry, aquifer properties, "
        "boundary conditions) is re-expressed for MF6 in `[layers] [uzf] [seep]` and the `MF_ws/*.asc` grids; the "
        "NWT-only solver and package keys are recorded as retired, not carried across"]],
      [2.6, 8.4, 5.6], font=7.5)
body("Also parameterised, and outside the ini files: **`MF_ws/inputSOILparam.txt`** — per soil zone, the number "
     "of soil layers `nsl`, and per layer the soil type `st`, thickness proportion `slprop`, `Smax`, `Sfc`, "
     "`Sr`, `Si` and `Ks`. That is the file the WP7 `soil_*` parameter group calibrates against, and it must be "
     "externalised in WP7.1.")


# =============================================================== 7
h(1, "7  Configuration reference")

h(2, "7.1  code/mm_paths.py — the machine layer")
code([
    '"""The only machine-specific paths. Edit here, or set the MM_* environment variables."""',
    "import os",
    "from pathlib import Path",
    "",
    "REPO       = Path(__file__).resolve().parents[1]                     # the MARMITES checkout",
    "DATA_ROOT  = Path(os.environ.get('MM_DATA_ROOT',  r'E:\\00code_ws\\LAMATA_new'))",
    "WS_ROOT    = Path(os.environ.get('MM_WS_ROOT',    r'E:\\00code_ws\\LaMata_MM-MF6'))",
    "NWT_REF    = Path(os.environ.get('MM_NWT_REF',    r'E:\\00code_ws\\LaMata_new_PhD_artigo_2s3L\\_h5_MM.h5'))",
    "MODFLOW_DIR= Path(os.environ.get('MM_MODFLOW_DIR',r'C:\\00MODFLOW'))",
    "",
    "LIBMF6     = str(MODFLOW_DIR / 'mf6.7.0_win64' / 'bin' / 'libmf6.dll')",
    "MF6_EXE    = str(MODFLOW_DIR / 'mf6.7.0_win64' / 'bin' / 'mf6.exe')",
    "PESTPP_IES = str(MODFLOW_DIR / 'pestpp' / 'bin' / 'pestpp-ies.exe')",
    "PYTHON_EXE = os.environ.get('MM_PYTHON_EXE', r'C:\\miniconda3\\envs\\flopy\\python.exe')",
    "",
    "CASE       = 'LaMata'                     # the case study",
    "DATASET    = REPO / 'example' / CASE      # Tier A — read directly by MM/MF",
    "GIS        = DATA_ROOT / 'GIS'            # Tier B — converter input only",
])
callout("Environment note:",
        "the flopy environment must be invoked **by prefix** — `C:\\miniconda3\\envs\\flopy\\python.exe` — because "
        "two identically-named environments exist on this machine and activating by name has crashed rasterio. "
        "`PYTHON_EXE` above is what the PEST workers will be given, so it must be the real one.")

h(2, "7.2  code/configs/lamata.toml — the model and run layer (skeleton)")
code([
    "[meta]",
    'config_version = 1',
    'name           = \"lamata_2lay_sfr_lak_crr\"    # -> out_<stamp>_<name>',
    "",
    "[paths]",
    'case      = \"LaMata\"                 # -> example/LaMata/, the Tier A inputs',
    'ws        = \"\"                       # blank -> WS_ROOT/MF6_ws',
    'gis_ws    = \"\"                       # blank -> DATA_ROOT/GIS (figures only, never the model path)',
    "",
    "[run]",
    'mode        = \"lagged\"      # lagged | iterative        (was --mode)',
    "relax       = 0.6           # iterative under-relaxation (was --relax)",
    "nsp         = 0             # 0 = all stress periods     (was --nsp)",
    "daily       = true          # false = aggregated         (was --aggregated)",
    "ats         = true          #                            (was --no-ats)",
    "build_only  = false         #                            (was --build-only)",
    "max_discrepancy = 1.0       #                            (was --max-discrepancy)",
    "",
    "[grid]                      # ---- WP1c ----",
    'kind = \"voronoi\"            # voronoi (DEFAULT) | structured | quadtree',
    "rebuild = false             # true = re-mesh; the cache carries the config hash",
    "  [grid.voronoi]",
    "  cell_far        = 100.0   # m, background cell size",
    "  cell_near_stream = 40.0   # m, stream corridor",
    "  stream_buffer   = 60.0    # m",
    "  stream_refine   = true    # false = SFRmaker style, network mapped onto the background",
    "  trans_levels    = [10.0, 20.0, 40.0, 70.0]   # Daoud-style graded transitions",
    "  seed_ponds      = true    # one pond-scale cell per pond, centroid-seeded",
    "",
    "[layers]",
    "nlay      = 2               # 2 reads _2s1L.ini directly (was --nlay)",
    "aggregate = false           # 6->2 derivation, comparison only (was --aggregate)",
    "",
    "[uzf]",
    "vks_scale = 1.0             #                            (was --uzf-vks-scale)",
    "",
    "[seep]",
    'kind = \"drn\"                # uzf | drn                 (was --seep)',
    "cond = 10000.0              # must be free-draining      (was --seep-cond)",
    "",
    "[et]                        # ---- WP2, all new ----",
    "uzf_et       = true         # UZF SIMULATE_ET on",
    'unsat_form   = \"etwc\"       # etwc | etae',
    "gwet_in_mf   = false        # MUST stay false: ETg comes from MM (guarded)",
    'extdp_source = \"veg_zone\"   # uniform | veg_zone | raster',
    "extdp_default = 2.0         # m, fallback rooting depth",
    "",
    "[sfr]                       # ---- WP3 ----",
    "enable = true               #                            (was --sfr)",
    "min_slope = 1e-4",
    'source = \"inputSTREAM.csv\"  # table from tools/gis_to_dataset.py',
    "monotonic_bed = true        # SFRmaker rule (Leaf et al. 2021)",
    "# parameter-source pattern (section 3a): each parameter says WHERE its value",
    "# comes from -- value | column | raster | drainage. Resolved by the CONVERTER,",
    "# never by the model, so no shapefile is ever on the model path.",
    "  [sfr.width]",
    "  drainage = { a = 0.5, b = 0.35 }  # w = a*A^b, A = drainage area -- CdL's method, THE DEFAULT",
    "  # column = \"width_m\"              # or: an attribute column of the stream shapefile",
    "  # value  = 2.0                    # or: uniform over the network",
    "  [sfr.manning]",
    "  value = 0.035",
    "  [sfr.rhk]",
    "  value = 0.1               # m/d streambed K            (was --sfr-rhk)",
    "  [sfr.rbth]",
    "  value = 0.5               # m streambed thickness",
    "",
    "[lak]                       # ---- WP4 ----",
    "enable   = true             #                            (was --lak)",
    "bedleak  = 1e-3             # 1/d                        (was --lak-bedleak)",
    "surfdep  = 0.05             # m",
    "maxiter  = 200              # LAK Newton cap",
    "stagechg = 1e-4             # m",
    'source   = \"inputPONDS.csv\"',
    "",
    "[crr]                       # ---- WP5 ----",
    "enable = true",
    "beta   = 1.0                # Daoud Eq. 23 partitioning (calibrated 0.8-1.0)",
    'sinks  = \"evaporate\"',
    'dem    = \"inputDEMfill.asc\"',
    "",
    "[spinup]",
    "cycles     = 1              #                            (was --spinup)",
    "tol        = 0.05           #                            (was --spinup-tol)",
    'strt_heads = \"hi_spinup\"    #                            (was --strt-heads)',
    'steady_means = \"hi_spinup\"  #                            (was --steady-means)',
    "strt_dem   = []             # [a, b] regression          (was --strt-dem)",
    "",
    "[postproc]",
    "enable = true               #                            (was --postproc)",
    "preproc = true              #                            (was --preproc)",
    "only    = false             #                            (was --postproc-only)",
    "sankey_min_flux = 0.05      #                            (was --sankey-min-flux)",
    "sankey_full = true          #                            (was --no-sankey-full)",
    "map_days = 6                #                            (was --map-days)",
    "",
    "[pest]                      # ---- WP7 ----",
    "use_pest_params = false",
    'params_npz      = \"pest_optimised.npz\"',
    "",
    "[ui]                        # ---- WP1b, read ONLY by app/, never by the model ----",
    'execution = \"local\"         # local | server  -- where the Run page launches the model',
    'runs_dir  = \"\"              # blank -> WS_ROOT/runs   (log + status.json per launch)',
    "port      = 8501",
    'address   = \"localhost\"     # \"0.0.0.0\" when serving from the server over the LAN',
    "poll_secs = 3               # log-tail refresh interval on the Run page",
])

h(2, "7.3  Rules the configuration layer must enforce")
bullet("**Unknown keys raise.** A typo must not fall through to a default.")
bullet("**The resolved configuration is written into the run folder** and hashed into the coupled HDF5.")
bullet("**Caches carry the hash they were built from**; a mismatch stops the run naming the offending key. "
       "The CdL grid-design cache silently overrode a user setting once — that is the whole reason for this rule.")
bullet("**`--set section.key=value` is the only override**, and every override is echoed at startup.")


# =============================================================== 8
h(1, "8  PEST-IES calibration design")

h(2, "8.1  Parameters")
table(["Group", "Parameter", "Style / bounds", "Notes"],
      [["`kh{1,2}`", "horizontal K per layer, pilot points", "multiplier, 0.2 – 5",
        "geostatistical structure per layer; range ≈ 4 × pilot-point spacing. Kv follows Kh at a fixed ratio in "
        "the forward run"],
       ["`sy`, `ss`", "specific yield / storage per layer", "multiplier, 0.2 – 5", "one constant per layer to start"],
       ["`drn_cond`", "outlet drain conductance", "multiplier, 0.05 – 20", "split by layer if the two-unit "
        "behaviour repeats what CdL saw"],
       ["`drnseep_cond`", "surface seepage-drain conductance", "multiplier, 0.05 – 20",
        "the baseflow lever; keep it free-draining"],
       ["`ghb_cond`, `ghb_head`", "lateral boundary", "multiplier 0.05 – 20 / additive ±5 m", "if GHB is retained"],
       ["**`sfr_rhk`**", "**streambed K**", "**multiplier, 0.05 – 20**", "**no field value exists — see §6.3**"],
       ["**`lak_bedleak`**", "**lakebed leakance**", "**multiplier, 0.1 – 10**",
        "**the NbS recharge lever; physical range 7e-4 – 1.5e-3 /d**"],
       ["**`crr_beta`**", "**cascade partitioning β**", "**value, 0.5 – 1.0**", "**Daoud calibrates 0.8 – 1.0**"],
       ["**`soil_*`**", "**per soil zone and layer: θs, θr, θfc, Ks, soil thickness**",
        "**multiplier, moderate bounds**", "**what the soil-moisture group actually informs**"],
       ["**`et_extdp`**", "**rooting / extinction depth per vegetation zone**", "**multiplier, 0.5 – 2**",
        "**splits ETuzf against ETg — the ET group's main lever**"],
       ["**`et_partition`**", "**PE / PT split, interception coefficient**", "**multiplier, narrow**",
        "**only if the ET residuals demand it; correlated with `et_extdp`**"]],
      [2.2, 4.6, 3.4, 6.0], font=8)

h(2, "8.2  Observations")
table(["Group", "Source file (WP6 export)", "Content", "Weighting"],
      [["`head`", "`obs_heads.csv`", "11 monitoring points, period-end heads",
        "1/σ with σ from the measurement record; then balanced against the other groups"],
       ["**`sm`**", "**`obs_sm.csv`**", "**θ per probe depth at the points that have probes — from "
                    "`MM_S[:,:,iSsoil_pc_s]` against `inputObsSM_*`**",
        "**the new group; weight so it carries a share comparable to `head`**"],
       ["`sfr`", "`obs_sfr.csv`", "outlet streamflow (and inlet, if one is prescribed), period-end",
        "relative weighting on the log of discharge, so low flows are not drowned by peaks"],
       ["**`et`**", "**`obs_et.csv`**", "**total ET = Ei + Eow + ETsoil + ETuzf + ETg, catchment or zonal, "
                    "monthly**", "**against the EC tower or a MODIS zonal aggregation (§6.3)**"]],
      [1.8, 3.6, 6.4, 4.4], font=8)

body("Start with **group-balanced** weights — each group contributing an equal share of the initial phi — and "
     "adjust only with a reason. This is where the model is most easily fooled: heads alone will happily be fitted "
     "by a checkerboard K field, and it takes streamflow, soil moisture and ET together to constrain the "
     "partitioning that MARMITES exists to compute.")

h(2, "8.3  The forward run")
code([
    "  forward_run_mm.py  (executed in every PANTHER worker directory)",
    "  ───────────────────────────────────────────────────────────────",
    "   1. pyemu.helpers.apply_list_and_array_pars(...)        # multipliers -> MF6 + MM input files",
    "   2. derive Kv = KV_RATIO * Kh for every layer file",
    "   3. run the COUPLED driver: MM + MF6 through libmf6     # NOT mf6.exe -- this is the CdL difference",
    "   4. collapse every time-series output to period-end rows (ATS determinism)",
    "   5. write obs_heads.csv / obs_sm.csv / obs_sfr.csv / obs_et.csv",
    "  guarded by  if __name__ == '__main__': multiprocessing.freeze_support()",
])

h(2, "8.4  Run-time budget — decide this before building the interface")
body("A full La Mata run is 1949 daily stress periods and takes roughly 11–12 minutes coupled. IES cost is "
     "`num_reals × (noptmax + 1)` runs:")
table(["Configuration", "Runs", "At 12 workers, 12 min/run", "Verdict"],
      [["150 reals × 3 iterations, full period", "600", "≈ 10 hours", "feasible on a server, overnight"],
       ["150 reals × 3 iterations, full period, 30 % failures", "600", "≈ 10 h + retries", "the realistic case"],
       ["50 reals × 3 iterations (pilot)", "200", "≈ 3.5 hours", "**start here** — proves the interface"],
       ["150 × 3 on a shortened window (e.g. the last 10 years)", "600", "≈ 2.5 hours",
        "the fallback if the full period is too slow; CdL shortened its window for exactly this reason"],
       ["Aggregated stress periods instead of daily", "600", "much faster",
        "**check first** whether the soil-moisture group survives aggregation — it probably does not"]],
      [5.4, 1.6, 4.4, 4.8], font=8)

callout("Measure before you commit:",
        "step 7.4's first job is to time one forward run end to end on the target machine, with the "
        "post-processing included. Every number above is an estimate; the calibration window and the ensemble "
        "size follow from the measurement, not from the plan.")


# =============================================================== 9
h(1, "9  Risk register")
table(["Risk", "Impact", "Mitigation"],
      [["UZF `PET` is not exposed as a writable API array in this MF6 build", "WP2 blocked",
        "probe first (`--probe` equivalent) before coding. Fallback: rewrite the UZF period data per stress period "
        "and reload — slower but certain. Trap 5 applies: validate by size, not by the variable list"],
       ["The residual-PET chain double-counts or starves an ET source", "wrong ET, silently",
        "the PET-closure test (2.8), the GWET guard (2.6) and the runtime PET-balance line (2.5b)"],
       ["**ETg's residual is computed from the PET written to UZF instead of ETuzf ACTUAL**",
        "ETg is starved exactly when the deep unsaturated zone is dry — i.e. when phreatic uptake matters most. "
        "Totals look plausible; the **sourcing** is wrong",
        "step 5 consumes `etuzf_actual_prev` read back from the UZF budget (2.4, 2.5); the runtime PET-balance "
        "line and the closure test both catch a regression"],
       ["CRR's topographic ordering breaks the cell-list contract", "coupler and DISV mapping desynchronise",
        "store a permutation at build time, never re-sort; assert it is a permutation; the existing cell-order "
        "tests must stay green"],
       ["LAK does not converge coupled (it never has here)", "WP4 stalls",
        "carry the CdL settings that fixed it — surfdep smoothing, `maxiter 200`, `stagechg 1e-4`, equilibrium "
        "initial stage. Bisect: LAK alone before LAK + SFR + MVR"],
       ["The PEST forward run is too slow for a useful ensemble", "WP7 unusable",
        "§8.4 — measure, then choose the window and the ensemble size. Pilot at 50 realisations first"],
       ["IES overfits to a checkerboard K field", "a data-fitting parameter set, not a physical one",
        "the full CdL anti-checkerboard package from the first build (WP7.5 callout)"],
       ["ATS makes the observation files non-deterministic", "instruction files misalign, realisations fail",
        "collapse to period-end rows identically in the reference and in every forward run"],
       ["Geospatial dependencies creep onto the model path", "the build stops being portable",
        "the WP1 converter and the repo-hygiene test"],
       ["The Streamlit app becomes a hidden dependency of the model", "the model stops running headless from "
        "Spyder or from a PEST worker",
        "`code/` may never import `code/app/` or `streamlit` — asserted by the hygiene test (WP1b.10), not by "
        "intent. Streamlit lives in a separate optional environment group"],
       ["A model run is launched in-line from the app", "the app blocks for 11 minutes and the run dies on a "
        "browser refresh", "launch detached, poll a log file and a `status.json` (WP1b.6); Streamlit reruns the "
        "whole script on every interaction, so nothing long-running may sit in the page body"],
       ["A stale cache silently overrides the configuration", "a wasted multi-hour run",
        "configuration hashes on every cache (WP0.6). This has already happened once, on CdL's grid design"],
       ["**Voronoi invalidates La Mata's validated baseline**", "the RMSE, the saved `hi_spinup_*` state and the "
        "MODFLOW-NWT comparison no longer apply like-for-like",
        "accepted deliberately: the `structured` producer is kept **permanently** so the regression stays runnable "
        "on the grid it was established on. Regenerate the spin-up state on the new mesh (WP1c)"],
       ["Centre-sampling a raster onto larger Voronoi cells throws data away",
        "silently wrong soil and aquifer parameter fields",
        "WP1c.3: area-weighted for continuous fields, majority-vote for zone rasters; verify resampled means "
        "against the raster means"],
       ["The figure suite is rewritten for vertex grids", "weeks of validated work redone for no scientific gain",
        "WP1c.7: rasterise for display instead, and leave `MARMITESplot_v3.py` untouched. True polygons live in "
        "the Streamlit app"],
       ["WP4 and WP5 each perform their own surgery on the cell list",
        "the Phase-2 cell-list contract breaks and the coupler desynchronises",
        "one shared per-cell surface descriptor built at build time (WP4.6), consumed by both"],
       ["`MM-MF6` → `master` merge conflicts", "merge stalls",
        "7 known modify/delete conflicts on `trunk/pyEARTH1D/*` (as they are on `master`); decide pyEARTH1D's fate before merging (WP8.5)"]],
      [4.8, 3.4, 8.4], font=8)

# =============================================================== 10
h(1, "10  Milestones and acceptance checklist")
table(["WP", "Milestone", "Green when", "Est."],
      [["WP0", "Configuration file drives the model", "the canonical run reproduces byte-identical MF6 inputs from "
               "`code/configs/lamata.toml`; no flag remains except `--config`, `--set`, `--run-tag`", "DONE"],
       ["WP1", "Data tiers separated", "a fresh clone + `DATA_ROOT` runs; the repo-hygiene test passes; the dataset "
               "README lists every Tier A file", "DONE"],
       ["WP1b", "Streamlit interface (stage 1)", "one setting serves the app on this machine and on the server; "
                "every §6 input is inspectable, shapefiles on a map; a run launches, survives a browser refresh "
                "and is reviewed from the browser; the hygiene test forbids `code/` importing streamlit",
        "DONE"],
       ["WP1c", "One grid path, Voronoi by default", "a full coupled La Mata run completes on a Voronoi mesh, "
                "mass-conserving, whole figure suite drawn with `MARMITESplot_v3.py` untouched; "
                "`kind=\"structured\"` still reproduces the pre-WP1c run exactly; the three-rung validation "
                "ladder passes", "DONE (see note)"],
       ["WP2", "Three-source ET", "one-year run splits ET five ways, closure test green, UZF groundwater-ET term "
               "exactly zero", "3–4 d"],
       ["WP3", "SFR coupled and validated", "one-year and full-period runs converge; outlet hydrograph against "
               "observation; mass balance < 0.1 %", "3–4 d"],
       ["WP4", "LAK coupled and validated", "stages coherent, LAK budget closes, MM↔LAK exchange live, LAK "
               "configuration map drawn", "4–5 d"],
       ["WP5", "CRR live", "synthetic mass-balance test exact; La Mata run converges; the change to the outlet "
               "hydrograph is quantified and explained", "4–5 d"],
       ["WP6", "Figures and observation exports", "`--postproc-only --preproc` emits everything with no skips; "
               "the four observation csvs have deterministic row counts", "3–4 d"],
       ["WP7", "PEST-IES calibrated", "prior MC completes with < 30 % failures; IES reduces phi; the parameter "
               "field is smooth (Moran's I ≥ 0, < 30 % of pilot points at bounds)", "1–2 wk"],
       ["WP8", "Tests, docs, merge", "suite green including the six mf6-binary tests; handoff document rewritten; "
               "merge rehearsed", "2–3 d"],
       ["WP9", "CdL as the second case study", "MMsurf runs from a daily CdL forcing record; a coupled CdL run "
               "converges with SFR/LAK/CRR on; nothing case-specific was added to `code/`",
        "gated by data"]],
      [1.2, 4.6, 8.6, 1.4], font=8)

body("The estimates are working days of implementation, excluding the multi-hour and overnight runs, which are "
     "Alain's to launch in Spyder.")

body("**One qualification on WP1c.** WP1c's milestone says \"the three-rung validation ladder passes\". "
     "Rung (a) passes exactly (DIS vs DISV-from-DIS: identical to 1e-9 mm/y on all fourteen terms). Rungs (b) "
     "and (c) are INCONCLUSIVE rather than passed, and deliberately so: they cannot be judged from a cold "
     "start, where runoff and exfiltration reach 30x precipitation because the water table begins above "
     "ground. They need a spin-up state generated on the mesh, then a run of at least one hydrological year "
     "on each grid — which is the modeller's to launch. The tool refuses to print PASS or FAIL on a "
     "spinning-up pair. See `docs/MM-MF6/WP1c8_grid_validation.md`.")

h(2, "Settled on 2026-09-09")
bullet("**Study area** — develop on La Mata, deliver on CdL (§1.4, WP9).")
bullet("**Grid** — Voronoi by default, one path for both case studies; structured kept permanently as the "
       "regression anchor (D5, WP1c).")
bullet("**SFR** — the CdL / SFRmaker method from the stream shapefile, with the parameter-source pattern and "
       "**drainage-scaled width kept as the default width producer** (§3a).")
bullet("**Open water** — LAK and SFR cells bypass the MM soil column (§4a). On the Voronoi mesh this is wholesale "
       "and exact, because each pond gets its own pond-scale cell; only stream cells keep a share.")
bullet("**Drain datums** — the **outlet** drains sit at the **bottom of the layer**, as in CdL, and only the one "
       "cell that becomes the stream outlet is replaced by a reach. **DRN-SEEP is excluded and stays at the land "
       "surface** (§3.3b, trap 15).")
bullet("**CRR** — cascade over the land surface, reinfiltration into the **topmost** soil layer; distinct from "
       "the bottom-up rejected-infiltration return, which is already coded (D4).")
bullet("**PET depletion** — ETg's residual is computed from **ETuzf actual**, read back from the UZF budget, "
       "never from the demand written to UZF (§2b).")

h(2, "Open questions still to settle")
bullet("**Q1 — Soil-moisture probe depths** — which MM soil layer does each `inputObsSM_*` column correspond to? "
       "The `sm` observation group cannot be built without this mapping. Needed before WP7.3.")
bullet("**Q2 — Observed ET** — is the eddy-covariance record usable directly, or do we take the CdL route "
       "(MODIS MOD16A2GF, aggregated to land-cover zones)? This decides §6.3 and the `et` group.")
bullet("**Q3 — Calibration window** — full 1949 stress periods, or a shortened window as CdL adopted? "
       "Answer after the WP7.4 timing measurement, not before it.")
bullet("**Q4 — CdL daily meteorology** — the WP9 blocker. A station record for the catchment, a national-network "
       "series disaggregated, or a downscaled reanalysis? Worth opening the enquiry now: it has a long lead time "
       "and it gates the whole of WP9.")
bullet("**Q5 — `grid.voronoi.stream_refine`** (new, WP1c) — the default here is `true`, but CdL settled on "
       "`false` (the SFRmaker approach: map the network onto the background) after judging a refined mesh "
       "\"too refined near streams\". Every La Mata run so far used `false`. It changes ncpl and therefore "
       "every downstream number, so decide before rung (c).")
bullet("**Q6 — zone-raster resampling** (new, WP1c.3) — majority vote is better per cell but biased in "
       "aggregate: on La Mata it misassigns 7.64 % of the active area against centre sampling's 3.66 %, and "
       "soil zone 1 drops from 13.1 % of the catchment to 9.5 %. `[grid] resample = \"centre\"` is the "
       "alternative. Which matters more here, per-cell fidelity or zone proportions?")
bullet("**Q7 — `seep.cond` on a mesh** (new, WP1c) — it is a per-cell conductance and therefore "
       "grid-dependent. 10 000 m²/d holds the seepage within ~0.2 m on 50 m cells; the first mesh run wanted "
       "~1.5 × 10⁶. Retune against an equilibrated mesh state, not blind.")


pagebreak()
# =============================================================== 11
h(1, "11  Implementation status — what is done, and how to check it")

rich([("Written 2026-09-10. ", "b"),
      ("WP0, WP1, WP1b (stage 1) and WP1c are implemented and committed on branch "
       "`MM-MF6_SFR_LAK_CRR`. This chapter exists so the work can be VERIFIED rather than taken on trust: "
       "every claim below has a command beside it that reproduces it. Where the implementation departed "
       "from the plan in the chapters that follow, §11.4 says so and why.", "")])

h(2, "11.1  At a glance")
table(["Package", "Status", "What proves it"],
      [["WP0  configuration file", "DONE",
        "The ~40 command-line flags are replaced by `code/configs/lamata.toml`. A structured build is "
        "**byte-identical** to the pre-WP0 one across all 13 MF6 input files (bar flopy's timestamp line)."],
       ["WP1  data tiers + converter", "DONE",
        "`code/tools/gis_to_dataset.py` turns the GIS layers into four GRID-INDEPENDENT tables in "
        "`example/LaMata/`. No shapefile is read by the model. `test_repo_hygiene.py` asserts it."],
       ["WP1b  Streamlit interface", "DONE (stage 1)",
        "Five pages, driven in a real browser. A run can be launched from the UI and survives closing the tab. "
        "Stages 2–3 (richer Results, Calibration) follow WP6/WP7."],
       ["WP1c  one grid path", "DONE (1c.1 – 1c.8)",
        "A **coupled La Mata run completes on a Voronoi mesh**, mass-conserving (cumulative discrepancy "
        "0.000 % with SFR). 60 figures draw, against 61 on the structured grid. `MARMITESplot_v3.py` was not "
        "touched."],
       ["WP2 … WP9", "NOT STARTED", "Unchanged from the plan below."]],
      [3.2, 2.2, 11.2], font=8.5)

callout("The one thing that is NOT done:",
        "**Voronoi is not the default.** `[grid] kind` is still `\"structured\"`, and it stays that way until "
        "rung (c) of the validation ladder is judged on an equilibrated multi-year run — which is yours to "
        "launch. See §11.3 step 6 and `docs/MM-MF6/WP1c8_grid_validation.md`.")

h(2, "11.2  Before anything: activate the environment")
rich([("This matters more than it looks. ", "b"),
      ("Calling `envs\\flopy\\python.exe` directly does NOT put the environment's `Library\\bin` on PATH, so "
       "every delay-loaded native DLL fails to resolve and Windows kills the interpreter with "
       "**0xC06D007F** and no Python traceback. It hits numpy's MKL BLAS — `np.dot` alone is enough — and "
       "therefore matplotlib and rasterio too. Two days of \"pre-existing environment faults\" were this.", "")])
code([
    "# PowerShell — always by PREFIX (there are TWO envs called flopy; the other is empty)",
    '& "C:\\miniconda3\\Scripts\\conda.exe" run -p C:\\miniconda3\\envs\\flopy --no-capture-output python -u <script>',
    "",
    "# or activate the shell once, then plain python works",
    "conda activate C:\\miniconda3\\envs\\flopy",
])

h(2, "11.3  Six checks, in order")

rich([("1.  The test suite. ", "b"), ("Expect ", ""), ("353 passed / 6 skipped", "b"),
      (". If you see ~218 with two modules excluded, the environment was not activated (see §11.2).", "")])
code(["cd E:\\00code\\MM-MF6_SFR_LAK_CRR",
      "python -m pytest code/tests -q"])

rich([("2.  WP0 — the structured grid still writes the same model. ", "b"),
      ("This is the regression anchor: everything WP1c added must leave it untouched. Build into a scratch "
       "workspace and compare against a known-good one; all 13 files match once flopy's timestamp comment is "
       "ignored.", "")])
code(["python code/tests/run_lamata_mf6.py --config code/configs/lamata.toml \\",
      "       --set run.build_only=true --set run.nsp=10 --set paths.ws=<TMP>"])

rich([("3.  WP1 — the dataset is grid-independent. ", "b"),
      ("`inputSTREAM.csv`, `inputSTREAM_param.csv`, `inputPONDS.csv` and `inputWATERSHED.csv` carry model-CRS "
       "geometry with a provenance header naming the shapefile and its mtime. Re-generate them and confirm "
       "nothing changes:", "")])
code(["python code/tools/gis_to_dataset.py --case LaMata --dry-run"])

rich([("4.  WP1c — a coupled run on a Voronoi mesh. ", "b"),
      ("Takes a couple of minutes at 10 stress periods. Watch for the mesh line (ncpl, achieved cell size, "
       "resampling and coverage), then the solution check.", "")])
code(["python code/tests/run_lamata_mf6.py --config code/configs/lamata.toml \\",
      "       --set grid.kind=voronoi --set grid.voronoi.stream_refine=false \\",
      "       --set run.nsp=10 --set paths.libmf6=auto \\",
      "       --set spinup.strt_heads= --set spinup.steady_means= \\",
      "       --set postproc.enable=true --set postproc.preproc=true"])
bullet("`mesh: voronoi, ncpl=989` — cell area mean 9845 m², i.e. 99.2 m against the 100 m asked for.")
bullet("`resample=auto, 9.6 source cell(s) per mesh cell, mesh tiles 99.863 % of the grid rectangle`.")
bullet("`solution check: converged, cumulative discrepancy 0.0200%`.")
bullet("`native MARMITESplot: 60 figure(s)` in `<ws-root>/out_<stamp>_<tag>/_output/`.")

rich([("5.  WP1c.8 — rung (a) of the validation ladder. ", "b"),
      ("DIS against DISV-from-DIS: the same geometry re-expressed, so it must be identical. It is — all "
       "fourteen water-balance terms agree to 1e-9 mm/y.", "")])
code(["python code/tools/validate_grid.py --ladder --rungs a --nsp 10"])

rich([("6.  The remaining judgement — rungs (b) and (c). ", "b"),
      ("These cannot be settled from a cold start: at 10 stress periods runoff and exfiltration reach "
       "**30× precipitation**, because the water table begins above ground and recirculates. The tool detects "
       "that and reports INCONCLUSIVE rather than a misleading PASS or FAIL. To settle them you need a "
       "spin-up state generated ON THE MESH (the saved `hi_spinup_*.asc` is structured, and WP0.6 refuses it "
       "on a mesh — correctly), then a run of at least one hydrological year on each grid.", "")])
code(["python code/tools/validate_grid.py --compare <WS_structured> <WS_voronoi> \\",
      "       --labels structured voronoi --tol 15 --report rung_c.md"])

h(2, "11.4  The front end — how to launch it")
body("Streamlit is installed in the **`flopy`** environment (checked first: the conda solve touched nothing in "
     "the model stack). There is no second environment to manage.")
code(["cd E:\\00code\\MM-MF6_SFR_LAK_CRR",
      "streamlit run code/app/Home.py",
      "",
      "# it prints a local URL, normally http://localhost:8501 , and opens a browser",
      "# to serve it from the server instead, run it THERE and reach it over the LAN:",
      "streamlit run code/app/Home.py --server.address 0.0.0.0 --server.port 8501"])
table(["Page", "What it is for"],
      [["Home", "Case selector, machine paths (green = found), recent runs."],
       ["Inputs", "**Tier A** — every file the model reads, with size and a viewer (rasters as images, CSVs as "
                  "tables with their provenance header). **Tier B** — the source cartography on an OpenStreetMap "
                  "backdrop, with the CRS of each layer REPORTED rather than guessed. Also runs the WP1 "
                  "converter, preview first."],
       ["Configuration", "The TOML schema as widgets, every field with its units and meaning — including the "
                         "100 MMsurf parameters enumerated from `__inputMMsurf.ini`, shown pre-filled."],
       ["Run", "Launch a run with one-off `section.key=value` overrides, validated BEFORE launching. The run is "
               "detached: it survives closing the browser. Live log and status."],
       ["Results", "The figures a run wrote, with folder and name filters."],
       ["**Grid**", "**New (WP1c.7).** The TRUE cell polygons of the mesh — not the display raster the figure "
                    "suite uses — coloured by cell area, with the stream network, ponds, catchment boundary and "
                    "observation points overlaid, plus the shared-face topology report. Pick `voronoi` in the "
                    "*Grid kind* selector."]],
      [2.6, 14.0], font=8.5)
callout("If the Grid page says \"no mesh cached\":",
        "the mesh is built the first time a run uses that grid, and cached under "
        "`<ws-root>/MF6_ws_<kind>/_mesh/`. Run check 4 of §11.3 once and it will appear. A structured grid has "
        "no mesh to show, which is why the page opens empty on the default configuration.")

h(2, "11.5  Where the implementation departed from the plan")
table(["Plan said", "What was done", "Why"],
      [["WP1c.1 proves the mesh path with a **quadtree**, then WP1c.2 adds Voronoi.",
        "Voronoi first; quadtree works too but landed second.",
        "flopy's GRIDGEN wrapper needs `pyshp` to read its own output back, and the environment is not "
        "writable without elevation. Voronoi needs only `triangle.exe`, which was already present. `pyshp` "
        "has since been installed and the quadtree producer works (ncpl 7728, refined to 12.5 m along the "
        "streams)."],
       ["Cells `(cid, cid, 0, icell2d)` over `(ncell, 1)` column vectors (the Phase-4 `refined_cell_list`).",
        "`nrow := ncpl, ncol := 1` — arrays are `(ncpl, 1)`, and the ROW INDEX IS the icell2d.",
        "Because `i*ncol+j == i`, `build_cell_list`, `clsMF6._cellid` and `clsMF6._griddata` all work with NO "
        "change. The compacted convention needs an explicit node map and breaks the `_griddata` reshape."],
       ["`code/ppMF6/marmites_voronoi.py`",
        "`code/marmites_meshes.py` (the producers) + `code/marmites_mesh.py` (the projection).",
        "The producers are not Voronoi-specific — one function per `[grid] kind`, all returning the same MF6 "
        "DISV gridprops — and the projection is a separate concern from producing a mesh."],
       ["WP1c.7 rasterises the result vectors for `plotLAYER`.",
        "That, plus three more figure families that failed for unrelated reasons.",
        "The Sankey was reading `<name>.dis.grb`, which a DISV model never writes; the obs time series passed "
        "the wrong cellid arity to flopy; and `delc × delr` gives **1 m² per cell** on a mesh. None of them "
        "was the display raster."],
       ["The MODFLOW-NWT comparison runs on every run.",
        "It is SKIPPED on a mesh, with the reason printed.",
        "The NWT reference is a 65 × 60 structured run. Comparing it to a mesh is not like-for-like, exactly "
        "as §WP1c already anticipated; the new-run figures still draw."]],
      [4.6, 5.4, 6.6], font=8)

h(2, "11.6  Decisions waiting on you")
table(["#", "Question", "Where it bites"],
      [["1", "`grid.voronoi.stream_refine` defaults to **true** (this cookbook's choice), but CdL settled on "
             "**false** — the SFRmaker approach, mapping the network onto the background — after judging a "
             "refined mesh \"too refined near streams\". Every La Mata run so far used `false`.",
        "WP1c / WP3. Changing it changes ncpl and therefore every downstream number."],
       ["2", "Zone rasters: majority vote is better PER CELL but biased in aggregate. On La Mata it "
             "misassigns 7.64 % of the active area against centre sampling's 3.66 %, and soil zone 1 drops "
             "from 13.1 % of the catchment to 9.5 %. `[grid] resample = \"centre\"` is the alternative.",
        "WP1c.3. Affects ETsoil and percolation through the soil parameter set."],
       ["3", "`seep.cond` is per cell and therefore **grid-dependent**. 10 000 m²/d holds the seepage within "
             "~0.2 m on 50 m cells; the first mesh run wanted ~1.5 × 10⁶.",
        "WP1c.8 rungs (b) and (c). Retune it against an equilibrated mesh state, not blind."],
       ["4", "Two observation points, **SM and EC**, fall in the same cell on the 100 m Voronoi mesh, so their "
             "modelled series are identical by construction.",
        "WP7 calibration: two observations constraining one cell."]],
      [0.8, 10.4, 5.4], font=8)


# =============================================================== APPENDIX
h(1, "Appendix A  Flag → configuration key map")
table(["Current flag", "Configuration key", "Current flag", "Configuration key"],
      [["`--ws-root`", "`mm_paths.WS_ROOT`", "`--strt-heads`", "`spinup.strt_heads`"],
       ["`--ws`", "`paths.ws`", "`--steady-means`", "`spinup.steady_means`"],
       ["`--run-tag`", "`meta.name`", "`--save-strt`", "`spinup.save_strt`"],
       ["`--libmf6`", "`mm_paths.LIBMF6`", "`--save-means`", "`spinup.save_means`"],
       ["`--nlay`", "`layers.nlay`", "`--spinup`", "`spinup.cycles`"],
       ["`--aggregate`", "`layers.aggregate`", "`--spinup-tol`", "`spinup.tol`"],
       ["`--grid`", "`grid.kind` (new values)", "`--strt-dem`", "`spinup.strt_dem`"],
       ["`--mode`", "`run.mode`", "`--postproc`", "`postproc.enable`"],
       ["`--relax`", "`run.relax`", "`--preproc`", "`postproc.preproc`"],
       ["`--nsp`", "`run.nsp`", "`--postproc-only`", "`postproc.only`"],
       ["`--aggregated`", "`run.daily = false`", "`--gis-ws`", "`paths.gis_ws`"],
       ["`--no-ats`", "`run.ats = false`", "`--sankey-min-flux`", "`postproc.sankey_min_flux`"],
       ["`--seep`", "`seep.kind`", "`--no-sankey-full`", "`postproc.sankey_full = false`"],
       ["`--seep-cond`", "`seep.cond`", "`--map-days`", "`postproc.map_days`"],
       ["`--uzf-vks-scale`", "`uzf.vks_scale`", "`--sankey-obs-years`", "`postproc.sankey_obs_years`"],
       ["`--sfr`", "`sfr.enable`", "`--build-only`", "`run.build_only`"],
       ["`--sfr-rhk`", "`sfr.rhk`", "`--standalone`", "`run.standalone`"],
       ["`--lak`", "`lak.enable` + `lak.source`", "`--probe`", "(stays a CLI diagnostic)"],
       ["`--lak-bedleak`", "`lak.bedleak`", "`--max-discrepancy`", "`run.max_discrepancy`"],
       ["`--daily`", "`run.daily`", "`--allow-bad-budget`", "`run.allow_bad_budget`"]],
      [3.4, 4.5, 3.9, 4.8], font=8)

body("New keys with no flag today: the whole of `[et]`, `[crr]`, `[ui]`, `[pest]` and `[meta]`; the "
     "`[grid.voronoi]` block and `grid.rebuild` (WP1c); the `[sfr.*]` parameter-source sub-tables — `width` "
     "(with its `drainage` producer), `manning`, `rhk`, `rbth` — plus `sfr.min_slope`, `sfr.source`, "
     "`sfr.monotonic_bed`, `lak.surfdep`, `lak.maxiter`, `lak.stagechg`, `lak.source`. Note `grid.kind` now takes "
     "`voronoi | structured | quadtree` rather than the old `dis | disv`.")

doc.save(OUT)
print("written:", OUT)
