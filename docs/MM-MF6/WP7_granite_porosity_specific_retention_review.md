# Porosity, specific yield and specific retention of weathered and fissured granite

Literature review for the La Mata / Sardón model (2026-10-18). Starting
point: Francés et al. (2014). Purpose: realistic values per aquifer layer
for MODFLOW 6 (Sy) and its UZF package (θs, θr), and bounds for the PEST++
calibration.

## 1. Terms, and how UZF uses them

- **Total porosity n**: the pore volume per unit volume.
- **Specific yield Sy**: the water a saturated material releases by gravity
  drainage, per unit volume (also called drainable or effective porosity).
- **Specific retention Sr = n − Sy**: the water held against gravity,
  which is the field-capacity water content of a drained material.
  Woessner and Poeter (2020) give the relation as n_e = Sy + Sr.
- **Sy depends on the drainage time** (Johnson 1967). A pumping-test,
  water-table-fluctuation or watershed estimate is a different quantity
  from a laboratory column value.
- **UZF in this model** (`uzf.thtr_from = "sy"`): θs is the total
  porosity, and θr = θs − Sy is the specific retention. The kinematic wave
  moves only the water between θr and θs, i.e. the range Sy. UZF's ET
  stops at θr (`extwc_source = "thtr"`). θs on its own changes only the
  absolute water-content levels.

## 2. The weathering profile of granite

The conceptual model (Dewandel et al. 2006; Lachassagne et al. 2021), from
top to bottom:

1. **Saprolite (alterite)**: clayey-sandy in coarse granites, tens of
   metres thick where it hasn't been eroded. A laminated layer at its base.
   It holds most of the storage and has low hydraulic conductivity.
2. **Stratiform fissured layer (SFL)**: dense, mostly subhorizontal
   fracturing, decreasing with depth; about twice the saprolite's
   thickness. It carries most of the transmissivity and stores little.
3. **Fresh basement**: permeable only along deep fractures; its storage is
   negligible at catchment scale.

At Sardón, Francés et al. (2014) found:
- alluvium in the valleys;
- saprolite: about 20 m on the pediment, up to more than 45 m along the
  main fault (F1);
- the fissured layer below, from 20–40 m depth;
- fresh granite above 1000 Ω·m.

## 3. Values from the literature

### 3.1 Saprolite (weathered granite)

| Quantity | Value | Context, method | Source |
|---|---|---|---|
| Total porosity | 5–30 % | coarse granite saprolite, Burkina Faso | Compaoré et al. (1997), via Lachassagne et al. (2021) |
| Total porosity | 11 % (saprolite); 3.6–10.5 % (saprock); 1.7 % (fresh rock) | Viamão granodiorite profile, Brazil | Geoderma (2022), "Weathering of Viamão granodiorite, Part 1"; authors to be added |
| Total porosity | up to 40–50 % near the surface, falling to ~0 at the base of weathering (10–35 m) | seismic refraction and resistivity, Southern Sierra Nevada | Holbrook et al. (2014) |
| Saturated / residual water content | θs 15–37 %, θr 6–24 % | laboratory retention curves on weathered-granite cores | Rouxel et al. (2012), from the abstract |
| Plant-available water (field capacity minus wilting point) | ≥ 10 %; 12 %; 12.4 % | weathered granitic rock, Sierra Nevada | Jones and Graham (1993), via Sternberg et al. (1996); Graham et al. (2010); Hubbert et al. (2001) |
| MRS free water content | ~7 % (3–16 %) | 6 MRS sites, Sardón | Francés et al. (2014) |
| MRS free water content | 1.6–6 % | granite saprolites in France, Madagascar, India and Burkina Faso | Baltassat et al. (2005), Vouillamoz et al. (2005), compiled by Francés et al. (2014) |
| Specific yield | 1.9 % (0.8–4.6 %) | MRS with Cy = 0.28 (Vouillamoz et al. 2007), Sardón | Francés et al. (2014) |
| Specific yield | 1–2 % | granite, India and Burkina Faso | Rushton and Weller (1985); Compaoré et al. (1997), via Durand et al. (2017) |
| Specific yield | 2–5 % | pumping tests in schists and gneisses, Benin | Vouillamoz et al. (2014) |
| Specific yield | 10 % (saprolite) against 2 % (fissured layer) | two-layer model, Plancoët migmatite, Brittany | Durand et al. (2017) |

**Specific retention of saprolite**, Sr = n − Sy:
- With n = 10–30 % and Sy = 1–5 %, Sr is about **8–28 %**.
- The retention measurements agree: θr of 6–24 %, and at least 10 % of
  plant-available water held above the wilting point.
- A fine-grained, low-K saprolite like Sardón's (Francés et al. 2014)
  sits towards the high-Sr, low-Sy end.

### 3.2 Fissured layer (SFL)

| Quantity | Value | Context, method | Source |
|---|---|---|---|
| Drainage porosity (Sy) | 0.5–2 %, about 90 % of it in the rock blocks and clogged fractures | 35 m thick granite SFL, Maheswaram, India | Maréchal et al. (2004) |
| Total Sy, weathered–fractured aquifer | 6.3 × 10⁻³ | Neuman analysis and water-table fluctuation, India | Maréchal et al. (2004) |
| Sy against depth | 1.5 % → 0.2 % over 30–40 m below the base of the saprolite | watershed scale, India | Dewandel et al. (2010, 2017b) |
| Effective porosity, saprolite + SFL | 0.8–2.7 % (mean 1.5 %) | Maheswaram granite, 53 km² | Dewandel et al. (2012) |
| Effective porosity | 0.3–2 % | two watersheds of 700 and 1000 km², India | Dewandel et al. (2017b) |
| Sy, fissured layer | 2 % (another version of the paper: 4 %, 3–5 %) | Plancoët, Brittany | Durand et al. (2017) |
| Sy, fissured zone | starting value 1.5 %, range tested 0.1–10 % | fissured granite and valley mire, Massif Central | Duranel et al. (2021) |
| Sy, borehole NMR | 1.2 ± 0.5 % (calibrated cutoff); 2.1 % with the default cutoff | fractured A-type granite | Phillips et al. (2024) |
| MRS free water content | < 0.5 % (below detection) | Sardón, outside the fault zones | Francés et al. (2014) |
| MRS free water content | 0.2–2.5 % | fissured granite, Burkina Faso | Vouillamoz et al. (2005), via Francés et al. (2014) |
| Elastic storage coefficient | ~10⁻⁶ | fracture zone, granite, India | Boisson et al. (2015) |
| Share of storage | 80–90 % of the stored groundwater, where the saprolite is eroded | Brittany | Wyns et al. (2004) |

Notes:
- **Total porosity of the SFL** is rarely reported. The saprock and
  weathered-bedrock figures (3.6–10.5 % at Viamão) and the fresh-rock
  1–2 % bracket it at **about 2–10 %**.
- **Sr in the SFL** is therefore a few percent at most.
- **At Sardón**, Sy ≤ the MRS free water content, so a fissured layer
  below MRS detection (< 0.5 %) has Sy **below 0.5 %** outside the fault
  zones.

### 3.3 Fresh granite

- Matrix porosity: about 0.02–1.8 % in rock-physics compilations
  (Amadei, lecture notes), and 1.7 % at Viamão.
- Storage: negligible at catchment scale (Maréchal et al. 2004).

## 4. What it means for the La Mata model

**Current base configuration**, both layers alike:
- θs = 0.45 and Sy = 0.01, so θr = 0.44.
- The saprolite's specific retention of 0.44 is at the top of the
  literature range.
- θs = 0.45 is far too high for the fissured layer. Because the mobile
  range stays Sy, it doesn't change UZF's flow, but θr (and UZF's ET
  limit) is not realistic.

**Values consistent with the literature and Sardón's own data:**

| Layer | θs = n | Sy | θr = Sr = n − Sy | Basis |
|---|---|---|---|---|
| L1 saprolite | 0.25 (0.15–0.35) | 0.02 (0.008–0.046) | 0.23 | Sardón MRS Sy; saprolite porosity and retention data |
| L2 fissured layer | 0.04 (0.02–0.10) | 0.004 (0.001–0.02), decreasing with depth | 0.036 | SFL drainage porosity; Sardón MRS below 0.5 % |
| Fault zones (F1, F2/3) | higher in both layers | towards the high end of each range | | Francés et al. (2014), Table 5 |

**For the PEST++ calibration:**
- The current Sy multipliers (×0.03–30 on 0.01, i.e. 0.0003–0.3) allow
  values far outside the literature: 0.3 is a gravel's specific yield.
- Iteration 1 drove layer-1 Sy as low as 0.0007, which empties UZF's room
  for infiltration and feeds the rejected-infiltration flip-flop.
- Literature-based bounds would be about 0.005–0.05 for the saprolite and
  0.001–0.02 for the fissured layer.

**θs per layer** can already be set with per-layer rasters
(`uzf.thts = "thts_%d.asc"`).

## 5. Caveats

- Few studies measure n, Sy and Sr on the same material. Most Sr values
  here come from n − Sy across different studies, or from retention curves.
- MRS gives the free water content, not Sy. The conversion factor
  Cy = 0.28 comes from Burkina Faso granite (Vouillamoz et al. 2007).
  Vouillamoz et al. (2014) refined it for hard rocks in Benin.
- Scale matters: cores give the matrix, pumping tests and water-table
  fluctuation give the aquifer, and watershed budgets integrate the faults.
- The Hassan et al. (2014) GSFLOW model of Sardón is calibrated on 18 years
  of heads and flows. Its Sy values were not retrieved here and should be
  added from the paper.

## References

Checked against the article or a publisher/repository record in this
review: authors, titles and journals of Francés et al. (2014),
Lachassagne et al. (2021), Durand et al. (2017, with volume and pages),
Duranel et al. (2021), Hassan et al. (2014), Rouxel et al. (2012), Johnson
(1967) and the hard-rock titles from Lachassagne et al.'s reference list.
The other volume and page numbers are not all verified: check them before
citing in a publication.

Baltassat JM, Legchenko A, Ambroise B, Mathieu F, Lachassagne P, Wyns R,
Mercier JL, Schott JJ (2005). Magnetic resonance sounding (MRS) and
resistivity characterisation of a mountain hard rock aquifer: the
Ringelbach catchment, Vosges Massif, France. Near Surface Geophysics 3,
267–274.

Boisson A, Guihéneuf N, Perrin J, Bour O, Dewandel B, Dausse A, Viossanges
M, Ahmed S, Maréchal JC (2015). Determining the vertical evolution of
hydrodynamic parameters in weathered and fractured south Indian
crystalline-rock aquifers: insights from a study on an instrumented site.
Hydrogeology Journal 23, 757–773.

Compaoré G, Lachassagne P, Pointet T, Travi Y (1997). Évaluation du stock
d'eau des altérites : expérimentation sur le site granitique de Sanon
(Burkina Faso). IAHS Publication 241 (Hard Rock Hydrosystems).

Dewandel B, Lachassagne P, Wyns R, Maréchal JC, Krishnamurthy NS (2006). A
generalized 3-D geological and hydrogeological conceptual model of granite
aquifers controlled by single or multiphase weathering. Journal of
Hydrology 330, 260–284.

Dewandel B, Perrin J, Ahmed S, Aulong S, Hrkal Z, Lachassagne P, Samad M,
Massuel S (2010). Development of a tool for managing groundwater resources
in semi-arid hard rock regions: application to a rural watershed in South
India. Hydrological Processes 24, 2784–2797.

Dewandel B, Maréchal JC, Bour O, Ladouche B, Ahmed S, Chandra S, Pauwels H
(2012). Upscaling and regionalizing hydraulic conductivity and effective
porosity at watershed scale in deeply weathered crystalline aquifers.
Journal of Hydrology 416–417, 83–97.

Dewandel B, Caballero Y, Perrin J, Boisson A, Dazin F, Ferrant S, Chandra
S, Maréchal JC (2017b). A methodology for regionalizing 3-D effective
porosity at watershed scale in crystalline aquifers. Hydrological
Processes 31.

Duranel A, Thompson JR, Burningham H, Durepaire P, Garambois S, Wyns R,
Cubizolle H (2021). Modelling the hydrological interactions between a
fissured granite aquifer and a valley mire in the Massif Central, France.
Hydrology and Earth System Sciences 25, 291–319.

Durand V, Léonardi V, de Marsily G, Lachassagne P (2017). A new
quality-of-fit criterion to quantify the specific yield in a two-layer
hard-rock aquifer model. Journal of Hydrology 551, 328–339.

Francés AP, Lubczynski MW, Roy J, Santos FAM, Mahmoudzadeh Ardekani MR
(2014). Hydrogeophysics and remote sensing for the design of
hydrogeological conceptual models in hard rocks – Sardón catchment
(Spain). Journal of Applied Geophysics 110, 63–81.

Graham RC, Rossi AM, Hubbert KR (2010). Rock to regolith conversion:
producing hospitable substrates for terrestrial ecosystems. GSA Today
20(2), 4–9.

Hassan SMT, Lubczynski MW, Niswonger RG, Su Z (2014). Surface–groundwater
interactions in hard rocks in Sardon Catchment of western Spain: an
integrated modeling approach. Journal of Hydrology 517, 390–410.

Holbrook WS, Riebe CS, Elwaseif M, Hayes JL, Basler-Reeder K, Harry DL,
Malazian A, Dosseto A, Hartsough PC, Hopmans JW (2014). Geophysical
constraints on deep weathering and water storage potential in the Southern
Sierra Critical Zone Observatory. Earth Surface Processes and Landforms
39, 366–380.

Hubbert KR, Beyers JL, Graham RC (2001). Roles of weathered bedrock and
soil in seasonal water relations of Pinus jeffreyi and Arctostaphylos
patula. Canadian Journal of Forest Research 31, 1947–1957.

Johnson AI (1967). Specific yield — compilation of specific yields for
various materials. USGS Water-Supply Paper 1662-D. (Averages also in USGS
Open-File Report 63-59: clay 2 %, silt 8 %, sandy clay 7 %, fine sand
21 %, medium sand 26 %.)

Jones DP, Graham RC (1993). Water-holding characteristics of weathered
granitic rock in chaparral and forest ecosystems. Soil Science Society of
America Journal 57, 256–261.

Lachassagne P, Dewandel B, Wyns R (2021). Review: Hydrogeology of
weathered crystalline/hard-rock aquifers — guidelines for the operational
survey and management of their groundwater resources. Hydrogeology Journal
29, 2561–2594.

Maréchal JC, Dewandel B, Subrahmanyam K (2004). Use of hydraulic tests at
different scales to characterize fracture network properties in the
weathered-fractured layer of a hard rock aquifer. Water Resources Research
40, W11508.

Phillips S, et al. (2024). Borehole nuclear magnetic resonance estimation
of specific yield in a fractured granite aquifer. Groundwater.

Rouxel M, Ruiz L, Molénat J, Hamon Y, Chirié G, Michot D (2012).
Experimental determination of hydrodynamic properties of weathered
granite. Vadose Zone Journal 11(3).

Rushton KR, Weller J (1985). Response to pumping of a weathered-fractured
granite aquifer. Journal of Hydrology 80, 299–309.

Vouillamoz JM, Descloitres M, Toe G, Legchenko A (2005). Characterization
of crystalline basement aquifers with MRS: comparison with boreholes and
pumping tests data in Burkina Faso. Near Surface Geophysics 3, 205–213.

Vouillamoz JM, Baltassat JM, Girard JF, Plata J, Legchenko A (2007).
Hydrogeological experience in the use of MRS. Boletín Geológico y Minero
118, 531–550.

Vouillamoz JM, Lawson FMA, Yalo N, Descloitres M (2014). The use of
magnetic resonance sounding for quantifying specific yield and
transmissivity in hard rock aquifers: the example of Benin. Journal of
Applied Geophysics 107, 16–24.

Woessner WW, Poeter EP (2020). Hydrogeologic Properties of Earth Materials
and Principles of Groundwater Flow. The Groundwater Project, Guelph.

Wyns R, Baltassat JM, Lachassagne P, Legchenko A, Vairon J, Mathieu F
(2004). Application of proton magnetic resonance soundings to groundwater
reserve mapping in weathered basement rocks (Brittany, France). Bulletin
de la Société Géologique de France 175, 21–34.
