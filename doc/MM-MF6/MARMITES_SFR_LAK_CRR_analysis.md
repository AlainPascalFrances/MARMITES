# SFR / LAK / CRR in MARMITES cells — design analysis

**Date:** 2026-07-22 · analysis before implementation (task 3)

## 1. The central issue

SFR and LAK do not simply *add* physics to MARMITES — they **overlap** it.
MARMITES already owns a surface-water domain. Per Francés & Lubczynski (2023)
§2.1, the land surface is split into a vegetation sub-domain and *"the surface
water domain that handles water storage in land depressions **and stream
channels**. When the maximum capacity of the surface water storage is reached,
then Ro occurs."*

So before writing any package, the question is not "how do I add SFR" but
**"who owns the water in a cell where both MARMITES and SFR/LAK are active"**.

## 2. What the La Mata data already contains

Inspection of the dataset (not assumptions):

| Fact | Value |
|---|---|
| Channel cells (`inputPONDw.asc` > 0) | **244** (13.0% of 1,870 active cells) |
| Channel widths | 1.5 and 3.0 m |
| Channel max depth (`inputPONDhmax.asc`) | 1.0 and 1.5 m |
| Network topology | **single connected component** (8-neighbour), dendritic |
| Channel elevation range | 800.2 m → **735.3 m** (outlet) |
| Channel fraction of a 50 m cell | **3–6 %** |
| `Ssurf_max` in a 3 m-wide channel cell | ≈ 6.3 mm |
| Stream in the groundwater model today | **6 DRN cells at the outlet only** |
| Fate of MARMITES runoff today | **discarded** — Ro leaves the model unrouted |

Two conclusions follow immediately:

1. **The SFR network already exists in the MARMITES inputs.** The 244 cells
   with `PONDw > 0` form a connected, downhill-consistent stream network with
   width and depth per cell. No new delineation is needed.
2. **The real gap is routing, not storage.** Runoff is generated everywhere and
   then vanishes; stream–aquifer exchange exists only at 6 outlet cells. That
   is what SFR + CRR would fix.

## 3. The four genuine overlaps, and how each resolves

### 3.1 Channel storage — *no conflict*
MODFLOW 6 SFR routes flow as **steady within a time step**: the reach balance
is inflow = outflow + leakage + ET ± runoff, with no `dS/dt` (the `storage`
keyword only activates the developmental `dev_storage_weight` option). MF6 SFR
therefore does **not** duplicate MARMITES' `Ssurf`. MARMITES keeps depression
and channel storage; SFR takes the water that spills out of it.

### 3.2 Open-water evaporation — *conflict, resolve by disabling it in SFR*
MARMITES computes `Eow` from `Ssurf` using the channel geometry
(`Eosurf_max = cell_width · w · shapeFactor · Eo`). SFR can also evaporate.
**Recommendation: leave evaporation in MARMITES, set SFR `EVAPORATION = 0`.**
It keeps the ET partitioning/sourcing (the paper's whole purpose) in one place.

### 3.3 Infiltration pathway — *the subtle one*
Channel water in MARMITES infiltrates through `I` into **soil layer 1**
(Eq. 1b), then percolates, then reaches the aquifer through UZF. SFR streambed
leakage instead goes **directly to the aquifer**, bypassing soil and UZF.
Both are defensible; they are not interchangeable, and using both double-counts.

**Recommendation:** water still *stored in the cell* infiltrates via MARMITES
(`I`, unchanged); water *routed in the channel* exchanges via the SFR streambed.
This is consistent because they are physically different water: the first is
what `Ssurf` holds, the second is what has already spilled into the network.

### 3.4 Cell occupancy — *SFR is fractional, LAK is not*
A channel occupies **3–6 %** of a La Mata cell: switching MARMITES off in
channel cells would delete the soil water balance of 94–97 % of that cell's
area. So **SFR must be fractional — MARMITES stays fully active in channel
cells.**

A lake is the opposite: it occupies ~100 % of its cells. There, MARMITES has
nothing meaningful to compute (no soil column under an open water body).
**LAK is binary — MARMITES must be switched off in lake cells**, and LAK owns
precipitation, evaporation and seepage there. Those cells must then be excluded
from `build_cell_list()` so they never enter the soil model, and the catchment
water balance becomes MM(land) + LAK(water).

*La Mata has no lake*, so LAK is infrastructure for other catchments; it should
be implemented but will be inactive here.

## 4. Recommended architecture

```
MARMITES (Python)                         MODFLOW 6
─────────────────                         ─────────
per cell: Pe, Ssurf, Eow, I, Ro   ──┐
                                    │  Ro (m3/d per reach)
CRR cascade (new, in MARMITES):     ├──────────────────────►  SFR
  D8 over the DEM, re-infiltration  │                         routes flow,
  into each downslope soil column   │                         streambed leakage
  residual delivered to channel ────┘                         (EVAPORATION off)
                                                                   │
percolation ──────────────────────────────────────────────►  UZF ──┤
ETg ───────────────────────────────────────────────────────►  WEL  ├─► GWF
                                                              LAK ─┘ (lake cells:
                                                                     MM disabled)
```

### Why CRR belongs in MARMITES, not in MF6
Cascade routing with re-infiltration must interact with the **soil column**:
water arriving at a downslope cell re-infiltrates subject to Eq. 1b,
`I = min[Ssurf ; D_soil,1(φ1 − θ1)]`. MF6 has no cell-to-cell overland cascade
(UZF6 dropped `IRUNFLG`; MVR moves water *between packages*, not across the
land surface). MARMITES already holds the DEM, the soil state and the
infiltration law — CRR is a natural extension of Eq. 1, and the water it
re-infiltrates stays inside the MARMITES budget where it can be partitioned.

**Architectural consequence to flag now:** CRR makes `Ro` non-local. Cells must
be evaluated in **descending topographic order** within a stress period, so an
upslope cell is solved before the cell that receives its runoff. Today
`build_cell_list()` returns raster (i, j) order. This is a real change to the
Phase-2 cell-list contract — the coupler and DISV mapping both depend on cell
order, so the ordering must be introduced deliberately (a permutation stored
once at build time, not a re-sort inside the loop).

### What crosses the API each stress period
| Direction | Quantity | Target |
|---|---|---|
| MM → MF6 | percolation (m/d) | UZF `FINF` *(exists)* |
| MM → MF6 | ETg (m³/d) | WEL `Q` *(exists)* |
| **MM → MF6** | **routed runoff per reach (m³/d)** | **SFR `RUNOFF`** *(new)* |
| **MM → MF6** | lake inflow (m³/d), if any | **LAK `RUNOFF`** *(new)* |
| MF6 → MM | heads, exfiltration, rejected infiltration | *(exists)* |
| MF6 → MM | *(optional)* stream stage | for future channel-storage feedback |

## 5. Interaction with the rejected-infiltration decision

We recently routed UZF's rejected infiltration and groundwater discharge back
into the **bottom soil layer**, per Eq. 1 / §2.3. MF6 offers an alternative:
`GWDTOMVR` / `REJINFTOMVR` send that water straight to SFR through the MVR
package. **These must not both be used.**

Recommendation: **keep the current behaviour**. Refused water saturates the
soil column from below, emerges as `Exf_g1`, becomes Dunnian runoff `Ro`, and
*then* enters CRR/SFR. That chain is the paper's formulation and keeps the
water inside the MARMITES partitioning until it genuinely becomes overland
flow. MVR would short-circuit the soil column and bypass ET partitioning.

## 6. Decisions needed before implementation

1. **SFR reach geometry** — derive from `PONDw`/`PONDhmax` (width, depth) plus
   the DEM for slope, or do you have a separate reach table (lengths,
   streambed K and thickness, Manning n)? Streambed conductance in particular
   is not in the current dataset; the 6 DRN cells have conductance 0.025–0.035.
2. **Should the 6 outlet DRN cells be replaced by the SFR outlet**, or kept
   alongside it? Keeping both risks draining the same water twice.
3. **CRR routing rule** — D8 (single steepest descent) or multiple-flow-direction?
   And re-infiltration limited by Eq. 1b only, or additionally by a
   user-specified infiltration capacity as in Daoud et al. (2022)?
4. **Confirm evaporation stays with MARMITES** (SFR `EVAPORATION = 0`).
5. **Is a lake ever expected in La Mata**, or is LAK purely for other sites?
   (Affects how much of the LAK path we validate now.)

---

# REVISION after reading `trunk/SFR_LAK_CRR/cdl_gwf_model_fable_v2.py`

The reference implementation changes two of my recommendations and answers a
question I had left open. Three findings matter.

## R1. SFR leakage does *not* pass through UZF — and that is an MF6 limitation

An SFR reach is bound to a GWF cell by `cellid`; the `SFR-GWF` term exchanges
**directly with the aquifer**. MF6 has no unsaturated zone beneath streams
(MODFLOW-NWT's SFR2+UZF1 did have one; MF6 did not carry it over).

MVR cannot restore it either: what a package *provides* to MVR is its available
**water** (runoff / outflow / excess), not its streambed leakage, which is an
intrinsic package↔GWF connection. So SFR→UZF via MVR is possible for *delivering
surface water*, but not for streambed seepage.

Consequence for MARMITES: in a stream cell the water table is typically shallow
(La Mata's channel follows the valley floor), so bypassing the soil column is
acceptable there. Where it is *not* acceptable, the honest options are to give
the reach no streambed connection and let its water enter MARMITES' surface
store instead — a modelling choice, not something MF6 can do for us.

## R2. CRR: reuse Daoud's weights, but the receiver logic must change

The reference builds CRR as **static MVR movers** computed once from the grid
topology (`CRR_ENABLE`, `CRR_BETA`): each cell's rejected infiltration (UZF
provider) and seepage (DRN-SEEP provider) is split among downslope neighbours
with MFD weights `α_ij = β·S_ij/ΣS_ij`, receiver priority **LAK > SFR > UZF**,
and topographic sinks evaporate. That is elegant and cheap — the movers are
time-invariant, so the `.mvr` file is written once.

**But its reinfiltration receiver is a UZF cell, and that does not transfer to
MARMITES.** In the CdL model UZF *is* the soil zone. In MARMITES-MODFLOW, UZF is
only the percolation zone **below** the MARMITES soil column. Re-infiltrating
into UZF would inject water beneath the soil, bypassing `Ssoil`, `Esoil`,
`Tsoil` and Eq. 1b — i.e. it would silently defeat the ET partitioning and
sourcing that is the entire purpose of the model.

Since MVR can only move water between *MF6 packages*, and MARMITES is not one,
**the cascade must run in Python** — as originally proposed, but now explicitly
reusing Daoud's Eq. 23 formulation:

* provider = MARMITES `Ro` (which already contains the returned rejected
  infiltration and exfiltration, per the Eq. 1 / §2.3 routing we implemented);
* weights = `α_ij = β·S_ij/ΣS_ij` over downslope neighbours, `CRR_BETA` exposed;
* receiver priority **LAK > SFR > downslope MARMITES soil column**;
* topographic sinks evaporate (Daoud's convention);
* the residual reaching channel/lake cells is delivered to SFR/LAK through the
  API as reach/lake inflow.

Same physics and same parameters as the reference; only the reinfiltration
target differs, because in MARMITES the soil column is not UZF.

## R3. The reference *abandoned* `simulate_gwseep` — and we just enabled it

Verbatim from the script:

> *Surface-seepage drain — replaces UZF simulate_gwseep (deprecated + non-smooth,
> which caused single-cell limit cycles). A DRN at land surface in the top layer
> of every cell, with cubic smoothing over DDRN so groundwater discharge ramps in
> gradually instead of switching on/off.*

with `DRN_SEEP_COND = 10000 m²/d` and `DRN_SEEP_DDRN = 3.5 m` (a smoothing
depth), `mover=True`.

This is directly relevant to the convergence behaviour we are already seeing
(SP1 needed ~100 outer iterations). `SIMULATE_GWSEEP` switches discharge on and
off discontinuously; the smoothed drain ramps it in. **Recommendation: add a
`seep='uzf'|'drn'` option to `clsMF6`.** With `drn`, exfiltration is read from
the DRN-SEEP flows instead of UZF `GWD` and — unlike the reference — is returned
to the **MARMITES soil column**, not moved to SFR, preserving the Eq. 1 / §2.3
pathway. Worth testing against the current run: if it removes the iteration
spikes, it is the better mechanism.

## R4. Two-layer model

The MF ini already declares the intent: `Mnlay = 2` with
`Mlay = [1, 1, 1, 2, 2, 2]` — the six numerical layers map to **two
hydrogeological units**. The reference model follows the same philosophy
("nlay=3 — ONE numerical layer per geologic unit, Daoud-style").

Proposed aggregation, 6 → 2, driven by `Mlay`:

| Property | Rule | Rationale |
|---|---|---|
| `botm` | bottom of the unit's lowest layer | geometry preserved |
| `k` (horizontal) | thickness-weighted **arithmetic** mean | preserves transmissivity |
| `k33` (vertical) | thickness-weighted **harmonic** mean | preserves vertical resistance |
| `ss` | thickness-weighted arithmetic mean | storage preserved |
| `sy` | value of the unit's **uppermost** layer | `sy` acts at the water table |
| `idomain` | active if **any** constituent layer is active | keeps the footprint |
| `icelltype` | 1 (convertible) | unchanged |
| `outcropL` | recomputed on the 2-layer grid | drives the MARMITES cell list |

Benefits: ~3× fewer cells and UZF objects (11,472 → ~3,700), faster solves, and
a model whose numerical layering matches its conceptual layering. Risk to check:
vertical discretisation of the water-table zone becomes coarser, so `sy`
behaviour and the outcrop mapping should be compared against the 6-layer run
before adopting it as the reference configuration.

## 7. Revised decision list

1. **Two layers** — confirm the `Mlay = [1,1,1,2,2,2]` aggregation above, and
   whether it becomes the default or an option (`--nlay 2|6`) so the 6-layer
   run stays available for comparison.
2. **Seepage mechanism** — `simulate_gwseep` (current) or the smoothed
   `DRN-SEEP` of the reference? I recommend implementing both and testing,
   given the convergence spikes.
3. **Streambed conductance** for the 244 reaches — not in the dataset. The
   reference uses a seepage-drain conductance of 10,000 m²/d; La Mata's six
   outlet drains use 0.025–0.035. What value/derivation for the streambed?
4. **The six outlet DRN cells** — replace by the SFR outlet, or keep both?
   (Keeping both drains the same water twice.)
5. **`CRR_BETA`** — reference uses 1.0 (Daoud calibrates 0.8–1.0). Same here?
6. **Evaporation** stays with MARMITES (`SFR EVAPORATION = 0`) — confirm.
7. **LAK** — no lake in La Mata; implement the path but leave it inactive?

---

# 8. IMPLEMENTATION RECORD (decisions taken, code written)

## 8.1 Decisions

| # | Question | Decision | Where |
|---|----------|----------|-------|
| 1 | Two layers | Yes, `Mlay=[1,1,1,2,2,2]` → 2 units, as `--nlay 2` (6-layer run stays available) | `ppMF6/marmites_layers.py` |
| 2 | Seepage mechanism | Implement **both**, selectable `--seep uzf\|drn` | `marmites_mf6.py` |
| 3 | Drain conductance | **10,000 m²/d** per cell — revised from 10 after the first coupled run (see §8.6) | `drn_seep_cond` |
| 4 | Outlet | SFR outlet reaches **replace** the 6 outlet DRN cells; DRN-SEEP on non-SFR cells | `marmites_mf6.build()` |
| 5 | `CRR_BETA` | 1.0 | pending |
| 6 | Evaporation | MARMITES keeps it: SFR `EVAPORATION=0`, LAK has no evaporation | `_add_sfr_package`, `_add_lak_package` |
| 7 | LAK | Built from `GIS/lm_ponds.shp` — see §8.3 | `ppMF6/marmites_lak.py` |

## 8.2 SFR

The stream network was already in the MARMITES inputs: `inputPONDw.asc` gives a
channel width for 244 cells (1.5–3.0 m) and `inputPONDhmax.asc` a channel depth
(1.0–1.5 m). One reach per stream cell.

**Routing had to be reconsidered.** Greedy steepest descent restricted to stream
cells does not work on La Mata: the channel crosses flat stretches where every
neighbour shares the same DEM value, and descent terminates in local sinks —
it stranded **211 of the 244 cells**. The routing is instead a *priority flood
grown upslope from the outlets*: the frontier is always expanded at its lowest
cell and each newly reached cell takes as receiver the cell it was reached
from. All 244 cells then route, flats included, with one dominant trunk
(220 cells) discharging at (7,0).

Reach 0 is deliberately the last cell popped — necessarily a headwater — because
MF6 writes a downstream connection as `-rno` and `-0` has no sign.

Bed tops are the DEM incised by the channel depth, then smoothed to be
downstream-monotonic (SFRmaker rule, Leaf et al. 2021): DEM noise left 7
receivers sitting higher than their own cell.

MARMITES runoff on channel cells is injected as SFR `INFLOW` through the API
each stress period (`MF6Coupler._write_runoff`).

## 8.3 LAK — the ponds are sub-grid

`lm_ponds.shp` holds 12 polygons of **341–2,036 m² against a 2,500 m² cell**:
every pond is smaller than one 50 m cell, and ponds 9 and 13 contain no cell
centre at all. Lakes therefore cannot be made by excavating cells — there are no
cells to excavate.

Resolution (user decision, matching the CdL reference model's own conclusion):
each pond is **one `EMBEDDEDV` lake inside one host cell**, with its real
geometry supplied by a stage-volume-area table rather than inferred from the
cell. The host is the cell containing the polygon centroid — a centroid always
lands somewhere, whereas a "cell centres inside the polygon" rule would drop
ponds 9 and 13. Two ponds sharing a host is rejected loudly, since `EMBEDDEDV`
requires the lake to be the only lake connection in its cell.

The table assumes a **wedge bathymetry** (area growing from ~1 % of the
footprint at the bed to the full polygon area at the rim) rather than a flat
bottom, so the wetted area falls smoothly to zero and a dry pond does not
degenerate numerically. `barea = sarea`: the wetted bed is the exchange area.

Each lake starts at `clip(water table, bed + 0.1, rim)`, so a perched pond
starts nearly empty and one in contact with the water table nearly full.

## 8.4 MVR — the stream runs through the ponds

11 of the 12 ponds sit on the channel. The reach hosting a pond hands its whole
flow to the lake (`sfr → lak`, FACTOR 1.0) and the lake's Manning outlet at the
rim spills into the reach downstream of the pond cell (`lak → sfr`). Because the
lakes are embedded rather than excavated, the reach can stay in the cell and no
reach excision or renumbering is needed — unlike the reference model.

## 8.5 Status

`--nlay 2 --seep drn --sfr --lak` builds and reloads: 2 units, 3,824 UZF objects
(from 11,472), 1,710 seepage drains, 244 reaches, 12 lakes, 22 movers.
Test suite: **174 passing**.

Still open: **CRR** (§5) — the Daoud et al. (2022) cascade with `β = 1.0`, which
needs a descending-topographic cell ordering and therefore changes the Phase-2
cell-list order contract.

## 8.6 Correction: the seepage-drain conductance

The first coupled run (`--nlay 2 --seep drn --strt-dem 0.9995 -2.0`, 1949 SPs)
was diagnostic. Convergence improved markedly — mean outer iterations **5.3 ->
2.23**, and 1.6–1.9 in the last two years — and exfiltration was healthy
(1,954 cells, 1,133 of 1,949 days). But the water table rose to **205 m above
the land surface**, with 1,903 of 1,954 cells above ground at some point. The
overshoot decayed away (228 of the 259 affected days fall in year 1, none after
SP 1023), but the first three years are not physically meaningful.

Cause: `drn_seep_cond = 10 m²/d` was far too small. Peak seepage on La Mata is
823.9 mm/d, i.e. **2,060 m³/d** over a 2,500 m² cell; a drain with C = 10 m²/d
only passes that flux when the head stands 206 m above the drain — which is
precisely what happened.

The underlying misreading is worth recording: **a seepage-face drain conductance
is a numerical device, not a physical property.** It is not a streambed or
aquitard conductance. Its only job is to pin the head at the land surface and
carry off whatever excess arrives, so it must be effectively free-draining. The
CdL reference model's own comment records the same finding from the other
direction: C = 1,000 "let it float", and 10,000 was kept.

| C (m²/d) | head excess at peak seepage |
|---|---|
| 10 (first run) | 206 m |
| 412 | 5 m |
| 2,060 | 1 m |
| 10,000 | 0.2 m |
| 20,600 | 0.1 m |

**Resolution:** default raised to **10,000 m²/d**. Two guards added so this
cannot pass silently again: `MF6Coupler.run()` now reports any water table
standing more than 1 m above the land surface (naming the conductance and the
value needed), and a regression test asserts the default stays free-draining.

## 8.7 The run that looked fine and was not

The second coupled run (C = 10,000 m²/d) was reported as an improvement on the
strength of its head field and iteration counts. It was not a solution at all.
MF6's own listing said so, in two places neither of which was read:

```
mfsim.lst : Solution 1 did not converge for stress period 2 and time step 1
lamatamm.lst:
    IN:   UZF-GWRCH =    1,205,303      <- all recharge, whole 1949-day run
    OUT:  DRN_SEEP  =    9,984,292      <- 8x the recharge ever received
    PERCENT DISCREPANCY = -131.69
```

The single non-converged stress period — the first transient day — discharged
**9.6 × 10⁶ m³**, 26.6× the exfiltration of the other 1,948 days combined. It
emptied the aquifer in one step. Everything downstream followed from that: the
water table fell monotonically from 3 m to 17.8 m below ground and never
recovered, because 5.3 years of recharge is a fraction of what was lost. The
NWT reference, by contrast, oscillates seasonally between 1.2 and 4.2 m below
ground in dynamic equilibrium.

The C = 10 m²/d run had the same defect in milder form: it spread the failure
over three years instead of one day, which is exactly why it read as a
plausible spin-up.

**Cause.** The steady state is solved with `perc_user` = 0.2 mm/d uniform and
leaves the water table 3.16 m below ground, at the surface in places. Day one
then applies real episodic forcing against a free-draining seepage boundary in
a single 1-day step. Newton cannot follow it, hits its iteration limit, and MF6
continues with whatever the solver last held.

**Fixes.**

1. **ATS** (`ModflowUtlats`, on by default, `--no-ats` to disable): MF6 may
   subdivide any transient period down to `dtmin = 1e-4 d`. The steady-state
   period is excluded.
2. **The coupler now follows MF6 across sub-steps.** `_advance()` previously
   assumed one time step per stress period; with ATS that would have left MF6
   behind MARMITES by a growing amount, the two models silently simulating
   different days. It now steps until MF6's own clock reaches the end of the
   period. The soil model still runs once per day and the exchanged fluxes,
   being period rates, are written once and left alone.
3. **`check_solution()`**, called by the runner and **raising by default**:
   parses `mfsim.lst` for non-converged periods and the listing for the
   cumulative `PERCENT DISCREPANCY`, failing above 1 %.

**Method note worth keeping.** Plausible heads and modest iteration counts are
not evidence that MF6 solved anything. The mass balance and the convergence log
are, and they are the first things to read — before the results, not after
someone doubts them.

Test suite: **190 passing**.

## 8.8 Two-layer parameters: read the file, don't aggregate it

A wrong assumption on my part, corrected. I had been deriving the 2-layer model
by collapsing the 6-layer one (harmonic/arithmetic means). But La Mata has an
authoritative 2-layer parameter file, `__inputMF_flopy_v3_2s1L.ini`, and its
values are hand-set, not aggregates:

  * Sy = 0.01 uniform in both layers (the 6-layer file has 0.05 / 0.01);
  * layer-2 Ss ~= 1e-7 (a confining value no thickness-weighted mean of the
    6-layer Ss could produce).

`--nlay 2` now reads that file directly. `--aggregate` keeps the old
collapse-the-6-layer path for comparison only. (Three stale lines in the
2-layer ini -- Mlay/h_plt/h_lbl still had six entries -- were trimmed to two so
the parser accepts it. A separate corruption of the 6-layer ini's `thick` line
was also repaired.)

## 8.9 Pre- and post-processing (marmites_postprocess.py)

Ported the CdL pre/post figure set to La Mata's structured grid -- self
contained, no rasterio/shapely/geopandas/pyproj (none of the CdL Voronoi
machinery is needed on a DIS grid), and validated against the real La Mata MF6
output already on disk.

preproc (<ws>/preproc/): MARMITES input maps (soil/meteo/irrigation zones, soil
thickness, pond width & depth, vegetation areas) AND the aquifer maps (top,
per-layer K / K33 / Ss / Sy / bottom, idomain, stream+pond overlay).

postproc (<ws>/postproc/): observed-vs-computed head time series at the active
piezometers (P0, C1-C3, ... from inputObsHEADS_*), mean-head and
mean-depth-to-water maps per layer, water budget by compartment from the
listing, UZF and SFR internal budgets, per-layer storage change, and the native
MARMITES plotLAYER head map. Each figure has a tidy CSV alongside.

Run with `--preproc --postproc`. The .cbc means are taken over an even
subsample (flopy needs ~30 s just to index a 1583-step budget file; the mean is
unaffected). No other CdL scripts are required.

Deferred: the native Sankey water-balance diagram (plotWBsankey). It reads the
MARMITES water-balance H5 with a bespoke per-flux structure and is a separate
wiring job from the groundwater-side figures done here.

Test suite: **198 passing, 6 skipped** (the skips need an mf6 binary, absent in
CI; they run locally).

## 8.10 Why ATS did not rescue the first transient day (and the fix)

The guard did its job on the `--nlay 2` run: it refused a solution with a
non-converged period and a 130% cumulative discrepancy, rather than reporting
it. Reading the listing showed the real mechanism:

  * ATS *was* active on stress period 2 (`ATS IS OVERRIDING TIME STEPPING`),
    but it set dt = 1.0 and never sub-divided. The residual sat around 14,000
    with head changes of 30-115 m per outer iteration, then the step advanced
    to SP 3 unconverged.

The cause was in the coupler, not MF6. The manual outer-iteration loop capped
at `max_outer = 200`, while the IMS solver's `OUTER_MAXIMUM` was 500. So the
loop stopped calling `solve()` and finalized the step *before* MF6 reached its
own non-convergence limit. MF6 therefore never registered the step as failed,
and ATS's retry-with-smaller-dt is triggered only by a registered failure. The
step was quietly accepted and the clock advanced.

Fix: `clsMF6` now exposes `outer_maximum`, and the coupler defaults
`max_outer` to it (500 here) instead of 200. MF6 now runs its full outer
budget, returns non-convergence, and ATS reduces dt (1.0 -> 0.2 -> 0.04 -> ...
down to dtmin = 1e-4) and retries the period. At small dt the storage term
dominates the diagonal and the step conditions well, so the stiff first day
should converge. A regression test asserts `max_outer` tracks the solver limit.

Note on the stiffness itself: the 2-layer file uses Sy = 0.01, a very small
water-table storage, so a small flux imbalance swings the head metres and the
free-draining seepage face then pulls it back -- a stiff feedback, worst on day
one when the state jumps from the steady-state IC to real forcing. ATS is the
right tool for it; the cap was preventing ATS from engaging. If the first day
still cannot converge even at dtmin after this fix, the next lever is the
steady-state initial condition (seed it with the mean MARMITES percolation
rather than the uniform perc_user), not the solver.

## 8.11 First valid coupled run, and the spin-up it revealed

With the max_outer/ATS fix the `--nlay 2 --seep drn --strt-dem 0.9995 -2.0` run
is numerically valid: 0 non-converged periods, cumulative mass balance 0.04%,
mean 2.3 outer iterations, exfiltration healthy (1,951 cells, 1,092 days).

But the water table declines monotonically from ~0 at day 1 to 18 m below ground
(mean) by year 5.3 and is still falling, reaching 34 m in the uplands. The NWT
reference sits at 1-4 m (mean) in seasonal equilibrium. The budget (post-steady
mean, m3/d) shows why:

| IN | m3/d | OUT | m3/d |
|----|------|-----|------|
| recharge (UZF->GW) | 618 | seepage face (DRN_SEEP) | 600 |
| storage release Sy | 508 | groundwater ET (WEL)    | 405 |
| storage release Ss |  52 | storage uptake          | 148 |
|                    |     | outlet drain            |  23 |

Recharge (618) cannot cover seepage + groundwater ET (~1005). The ~410 m3/d
shortfall drains from storage, and with **Sy = 0.01** that is ~3 m/yr of
water-table fall -- matching the observed ~18 m over 5.3 yr. The decline
decelerates (7.2 -> 3.2 -> 2.6 -> 1.9 m per successive year), so it is asymptoting
toward equilibrium, not diverging.

Diagnosis: this is a **spin-up transient from a non-equilibrium initial
condition**. The DEM regression (0.9995*elev-2) puts the water table ~2 m below
surface everywhere, far too high in the uplands where equilibrium is 15-34 m
deep. It sped per-step convergence (its purpose) but is not a steady state, so
the model spends the record draining toward one.

Fix (user decision): an **iterated spin-up loop**. `run_lamata_mf6.py --spinup N
[--spinup-tol M]` repeats the full forcing up to N times, each cycle starting
from the previous cycle's final head field (`clsMF6.strt_array`), until the mean
between-cycle water-table change drops below the tolerance (default 0.05 m). The
solution guard runs every cycle, so spin-up never iterates from an invalid
state. The equilibrated head field then serves as the IC for the reported run.

## 8.12 Persisting the equilibrated heads (spin-up once, reuse forever)

The spin-up result is now saved so it need not be repeated:

* `clsMF6.save_heads_asc(heads, prefix)` writes the final head field as one
  ESRI-ASCII grid per layer (`<prefix>_l1.asc` ..), inactive cells as nodata;
  `load_heads_asc(prefix)` reads it back.
* The runner **auto-saves** after any `--spinup > 1` run (default prefix
  `hi_spinup`, in `MF_ws/`), and `--save-strt PREFIX` saves after a single run.
* `--strt-heads PREFIX` loads such a field as the IC (via `strt_array`),
  overriding `--strt-dem`, so later runs start already equilibrated and skip
  the spin-up entirely.

Workflow:

```
# once: spin up to equilibrium and save the field
run_lamata_mf6.py --libmf6 ... --nlay 2 --seep drn --strt-dem 0.9995 -2.0 \
    --spinup 10 --spinup-tol 0.05 --postproc        # writes hi_spinup_l1.asc ..

# thereafter: reuse it, no spin-up
run_lamata_mf6.py --libmf6 ... --nlay 2 --seep drn --strt-heads hi_spinup --postproc
```

A round-trip test asserts the saved field reloads unchanged on active cells and
drives the IC when fed back.

## 8.13 Outlet-drainage check, steady-state fix, and native plots

**Outlet drainage is not the cause of the drawdown.** Put in the reference's
mm/yr terms, the current 2-layer DRN-seep run discharges *less* than the
pre-SFR/LAK reference, not more:

| flux (mm/yr) | current | reference |
|--------------|---------|-----------|
| outlet DRN (6 cells) | 1.7 | -- |
| seepage face -> soil (EXFg) | 44.9 | 51.3 |
| groundwater ET | 30.6 | 37.9 |
| **total GW discharge** | **46.6** | **51.3** |

The 6 outlet drains remove 1.7 mm/yr (negligible); the seepage face returns
44.9 mm/yr to the soil (not lost). Total discharge is below the reference.

**The real cause: the steady-state SP0 wipes the initial-head file.** A
steady-state period solves dh/dt = 0, so it ignores STRT and re-solves to its
own equilibrium under the uniform `perc_user` -- a too-wet, near-surface table
(3.2 m). The transient then drains from there every run. This also defeated the
spin-up loop: each cycle re-solved the same SP0, so "equilibrium after 2 runs"
was just the same transient repeated. The water table was still falling at
-1.1 m/yr at the end of each cycle.

Fix (user decision): drive SP0 with the mean actual recharge and groundwater ET
instead of `perc_user`.

* `MF6Coupler.steady_perc` / `steady_etg` -- per-cell means used at the steady
  step; fall back to `perc_user` when unset.
* Runner: `--save-means [PREFIX]` writes `<PREFIX>_perc.asc` / `_etg.asc`
  (auto after a spin-up); `--steady-means PREFIX` loads them so SP0 lands near
  the dynamic equilibrium. Combined with `--strt-heads`, a run starts
  equilibrated and stays there.

**Native MARMITESplot now runs in --postproc** (was only plotLAYER). A new
`native_suite()` drives the native module from the coupled run's in-memory
arrays and writes into `postproc/`:

* `native_wb_catchment.png` -- `plotTIMESERIES_CATCH`, the catchment
  water-balance series (combined MM + soil-layer flux array reassembled from
  `wb_ts` + `wb_ts_soil`).
* `native_native_map_{recharge,exfiltration,ETg,runoff}.png` -- `plotLAYER`
  time-mean flux maps.

Still pending: `plotWBsankey` (Sankey) and `plotCALIBCRIT` (obs-vs-computed
calibration) -- both need the driver's observation/soil data structures and are
the remaining native-suite items.

## 8.14 Correction: the drawdown is a transient recharge deficit, not the IC

My §8.11 diagnosis (IC / spin-up) was incomplete. Fixing the steady state does
not help, because the TRANSIENT itself loses water regardless of where it
starts. The spun-up run's aquifer budget:

  recharge reaching the water table (UZF-GWRCH) : 68 mm/yr
  discharge (seepage 66 + gwET 36 + outlet 2)   : 104 mm/yr
  deficit drained from storage                   : 36 mm/yr

The table falls until the thickening unsaturated zone throttles recharge into
balance (~18 m). Applied percolation is 165 mm/yr (matching the NWT reference's
164), but only 68 reaches the water table.

Cause: UZF6 forbids EPSILON < 3.5, so the NWT model's EPSILON = 2.0 was clamped
to 3.5. A higher Brooks-Corey exponent lowers the unsaturated relative
permeability K(theta) = VKS * Se^EPSILON, so percolation moves more slowly
through the deep unsaturated column and less reaches the water table -- a
drainage feedback the NWT model did not have.

Fix (user decision): raise the UZF unsaturated VKS to restore the recharge rate.

* `clsMF6.uzf_vks_scale` (runner `--uzf-vks-scale F`) multiplies ONLY the UZF
  packagedata vks column, never the aquifer NPF k33 (tested).
* The runner now prints an "aquifer balance" line after every run
  (recharge to WT vs discharge, mm/yr) so the scale can be calibrated to close
  the deficit.

Calibration: raise the scale until recharge ~ discharge (deficit -> 0), then
run the full `--spinup` to equilibrate and save `hi_spinup`. Because
K(theta) = VKS*Se^eps is non-linear (raising VKS lowers the equilibrium Se,
partly self-cancelling), the scale is found iteratively, not analytically; a
first try of ~2-3 is reasonable given the eps 2.0->3.5 change.

## 8.15 THE ACTUAL BUG: UZF infiltration bound to FINF, not SINF

VKS x20 changed nothing -- because VKS was never the throttle. Reading the UZF
budget of a run exposed it: UZF INFILTRATION is EXACTLY 977 m3/d on every single
stress period. 977 m3/d over the catchment is 0.2 mm/d = **perc_user**, the
static build-time value. So UZF has been applying a CONSTANT recharge the whole
time, ignoring MARMITES' daily percolation entirely.

Root cause: the coupler bound the UZF infiltration pointer as **FINF**. In
MODFLOW 6's memory manager the operative infiltration array is **SINF**; FINF is
only the flopy input keyword. Binding FINF returns a valid-but-non-operative
pointer, so every daily write was silently discarded and UZF kept using the
period-0 perc_user. WEL/Q bound correctly (that is why groundwater ET coupled
fine and only the recharge was wrong), which is what made it look plausible.

This invalidates the whole "drawdown" investigation above (§8.11-8.14): the IC,
the steady state, EPSILON and VKS were all red herrings. The water table drained
because the aquifer received a fixed 0.2 mm/d of recharge regardless of rainfall.

Also note: the numbers in §8.13-8.14 (68 vs 104 mm/yr etc.) were partly read
from mismatched files as the diagnosis iterated -- the on-disk .lst/.uzf.cbc
were from a 1-year `--nsp 365` test while the h5 was from a full run. Treat them
as unreliable; they are superseded by this section.

Fix:
* `MF6Coupler._bind` now tries `SINF` first, then `FINF` as a fallback.
* `check_solution` gains an infiltration-fidelity guard: if the written
  percolation varies (CV > 0.2) but UZF INFILTRATION is nearly constant
  (CV < 0.02), the run FAILS -- this exact silent failure can never recur.

Verification pending (workspace shell was unavailable when the fix landed):
1. `--probe` to confirm MF6 6.7.0 exposes `<MODEL>/UZF/SINF`.
2. a coupled run to confirm UZF INFILTRATION now tracks the daily percolation
   and the water table stabilises near the NWT 1-4 m.
3. once recharge couples, revert `--uzf-vks-scale` to 1.0 and re-spin-up; the
   EPSILON=3.5 clamp may still matter, but only as a second-order effect now.
Full test suite not re-run since the edit (shell jammed).

### 8.15.1 SINF was necessary but NOT sufficient -- the write ORDER

Probe confirmed `LAMATAMM/UZF/SINF` exists, and the coupler now binds it. But
the run STILL failed the new guard: UZF INFILTRATION stayed at 977 (perc_user).
So writing to SINF was not enough -- the second half of the bug is the write
ORDER.

In lagged mode the coupler wrote the fluxes and THEN called `_advance`, whose
`prepare_time_step` runs each package's read-prepare (`rp`). UZF's `rp`
re-loads SINF from the input file's period data (the build-time perc_user),
overwriting the value just written. So every daily write was reverted before the
solve. (WEL/Q looked fine earlier only because those numbers came from
mismatched files -- Q was being clobbered too.)

Fix: the exchanged fluxes are now written AFTER `prepare_time_step`, via a
`write_cb` threaded through `_advance` -> `_one_step`, re-applied on every ATS
sub-step (idempotent). Applied to the steady SP0, the lagged march, and the
iterative ATS continuation. This is the canonical modflowapi ordering
(prepare_time_step -> set inputs -> prepare_solve -> solve).

So the recharge decoupling had TWO causes, both now fixed: (1) wrong variable
(FINF vs SINF), (2) wrong write order (before vs after rp). The guard added in
8.15 catches either recurrence.

### 8.15.2 Pinned exactly: write after prepare_solve, not prepare_time_step

Writing after prepare_time_step (8.15.1) STILL failed -- UZF INFILTRATION stayed
at 977 = perc_user. `tests/diag_sinf.py` settled it by writing a known value at
three positions and checking whether it survived the solve:

    position                     stuck?
    A after prepare_time_step     no    (reverts to perc_user)
    B after prepare_solve         YES
    C before every solve()        YES

So `get_value_ptr` returns a live view (not a copy: the write is visible, MF6
overwrites it), and MF6 re-derives UZF SINF from its stored period data during
BOTH prepare_time_step and prepare_solve. Only a value written AFTER
prepare_solve is the one the solve uses.

Final fix: `_one_step` calls the flux `write_cb` after `api.prepare_solve(1)`
(position B). This covers the steady SP0, the lagged march and the iterative ATS
continuation. SFR runoff was folded into the same callback (it was previously
written after the advance -- position A -- so it would have failed identically
once SFR was used). The iterative outer loop already wrote after prepare_solve,
which is why that path was never the one under test.

Root-cause chain, complete: constant recharge <- SINF write not applied <-
written at the wrong point in the BMI step sequence (and, separately, had been
bound as FINF). All fixed; diag_sinf.py + the 8.15 guard prevent regression.

### 8.15.3 Why the fix still did not take: the bind gate rejected SINF

After 8.15.1/8.15.2 the run STILL gave a constant 977, yet diag_sinf.py proved a
post-prepare_solve write to SINF sticks. The difference was the BIND, not the
write. diag_sinf.py binds SINF with get_var_address + get_value_ptr directly;
the coupler's `_bind_first` additionally required the address to appear in
`get_input_var_names()`. UZF SINF is an advanced-package variable: reachable via
get_value_ptr but ABSENT from that input list. So `_bind_first` silently
rejected SINF and fell back to FINF (which IS listed) -- the coupler was writing
the non-operative array the whole time, correct position or not.

Fix: `_bind_first` no longer requires input-list membership. It validates a
candidate by SIZE (a wrong address returns a wrong-sized pointer -> heap
corruption, the original reason for the gate) via a `min_size` argument, and
accepts a reachable, correctly-sized pointer even if MF6 omits it from the input
list. The infiltration bind passes `min_size=ncell`, and the coupler now prints
`coupler: UZF infiltration bound to <addr>` so the bound variable is visible.

Confirmation to look for in the next run: that line must read
`.../UZF/SINF` (not FINF), and the aquifer-balance recharge must vary with
rainfall instead of sitting at 977 m3/d.

### 8.15.4 CONFIRMED FIXED, and the residual is a calibration matter (not code)

The run after 8.15.3 printed `coupler: UZF infiltration bound to
LAMATAMM/UZF/SINF`, passed the infiltration guard, and gave a recharge that
VARIES with rainfall (67.1 mm/yr). The water-table map now follows topography
(2-5 m in the valley network, 15-25 m in the uplands) instead of a uniform
drain to 18 m. obs_heads.png shows the computed heads START on the observations
at all four piezometers (P0, C1, C2, C3). The recharge coupling is fixed.

Residual: the computed heads then drift down ~2-3 m/yr away from the (stable)
observations -- recharge 67 vs discharge 97 mm/yr, a 30 mm/yr deficit drained
from storage. Two findings settle what this is:

1. `--uzf-vks-scale` 3 and 20 do NOT change it. Now that SINF is genuinely
   coupled, that means recharge-to-water-table is NOT UZF-throttled: the ~67
   mm/yr is what the calibrated MARMITES soil delivers net, and UZF passes it
   through ~1:1. VKS is a dead lever here.
2. The modeller reports the MF6 result is "pretty similar to the previous
   version" (MODFLOW-NWT). The NWT reference ALSO drains over the record
   (06_heads.png of figures_NWT2MF6: catchment head 778 -> 766 m). So the
   drawdown is a property SHARED by both models, present before the MF6
   conversion.

Conclusion: the Python-3 / MODFLOW-6 / UZF6 conversion is complete and behaves
like the NWT model it replaced -- which was the goal. The water table declining
below the observed levels is a pre-existing LA MATA CALIBRATION issue (recharge
vs discharge / boundary conditions / soil parameters), shared by the NWT and MF6
models alike, and is a hydrological-calibration task for the modeller, not a
coupling defect. `plot_water_budget.py` (now wired into --postproc) quantifies
the MF6-vs-NWT agreement flux by flux in 00_summary.txt / 03_wb_totals.png.
