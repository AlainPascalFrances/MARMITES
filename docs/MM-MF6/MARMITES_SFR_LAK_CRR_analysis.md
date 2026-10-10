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

> **Superseded by §8.16 (WP3, 2026-09-26)** for the network source, the
> outlet and the reach geometry. The routing principle below still holds.

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

_input (<out-dir>/_input/): every parameter field as a native `plotLAYER`
page, `IN_<nnn>_<name>` -- geometry (elev/top/botm/thick/strt), aquifer
properties (hk/T/Ss/Sy/vka), the DRN and GHB arrays, the UZF soil parameters,
the model footprint (ibound) and the MARMITES zoning (soil / meteo /
irrigation zones, vegetation areas) -- plus `IN_000_general_map.png`, the
site's general map rebuilt from the ArcMap GIS layers (soil types, irrigation
plots, ponds, hydrography, catchment boundary, monitoring and observation
points over shaded relief and elevation contours), after
`GIS/LaMata_MM_MF_202109.png`. That one needs geopandas and the GIS
workspace, and is skipped without them.

A second, plainer set of `aq_*` / `mm_*` imshow maps used to be drawn
alongside. It duplicated the native set field for field, in a style borrowed
from another project and with no coordinate frame, and has been removed.

_output (<out-dir>/_output/): observed-vs-computed head time series at every
monitoring point, water budget by compartment from the listing, UZF and SFR
internal budgets, and the native figure set (per-point time series, Sankeys,
calibration criteria, the maps). `MMmap_*` (the MM soil-water-balance fluxes)
and `GWmap_*` (the aquifer terms) come from ONE function, `_native_result_maps`,
so the two sets share a layout. The mean head, mean
depth-to-water and per-layer storage change are written as CSV only: the
native GWmap_head / MMmap_dgwt draw the same fields with the full axes.

Every map, in both folders, carries the same two frames: MODFLOW row/column
indices on the top and right, projected coordinates in km on the bottom and
left (`MARMITESplot_v3.add_real_coord_axes`, shared by the plotLAYER pages and
the imshow overlay).

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
  `wb_ts` + `wb_ts_soil`). **Removed 2026-09-09**: the reassembled array never
  filled the curves and the per-point series carry the same fluxes.
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

## 8.16 WP3: the network rebuilt the SFRmaker way (2026-09-26)

The network is no longer the 244-cell raster of §8.2. It is the mapped
hydrography (`inputSTREAM.csv`, 97 segments) burned onto the model grid.
On La Mata's Voronoi mesh that gives 2663 reaches and 14,238 m, equal to the
mapped length. The attributes follow SFRmaker (Leaf et al. 2021):

- **One outlet**, where the network leaves the catchment: the lowest stream
  cell with a face on the catchment boundary (mesh cell 4340). The old rule
  (every stream cell that was also an outlet DRN cell) gave 29 exits once
  the outlet drain became a line.
- **Reach length** = the channel length mapped inside the cell (0.1-10.5 m
  on the mesh). The old rule took index differences times delr/delc, which
  on the mesh's (ncpl, 1) proxy grid meant differences of *cell numbers*:
  1,289 km of reaches for 14.2 km of channel.
- **Slope** over the centroid spacing to the next reach. **Bed**, per the
  user rule of 2026-09-27: the stream's total depth below the land surface
  is the soil depth of its cell + the channel depth + the streambed
  thickness. So the channel is cut through the MMsoil column: its bed top
  sits one channel depth below the aquifer top (land surface minus soil
  thickness), then made downstream-monotonic (518 tops lowered, 578 reaches
  on the 1e-4 floor on the aquifer-top datum).
  *Why:* measured from the land surface (the first WP3 version), La Mata's
  1.5 m of soil equalled the 1.0 m channel + 0.5 m streambed. That put 2024
  of 2663 streambed bottoms exactly on the aquifer top, where the seepage
  drains hold the water table. MF6's stream-aquifer exchange switches at
  the streambed bottom, and the first SFR year took 40-500 outer iterations
  a period (5.5 without the stream), about 10 h a spin-up cycle. With the
  rule, every streambed bottom sits 1.5-2.96 m below the aquifer top.
- Manning's n, streambed K and thickness per segment from
  `inputSTREAM_param.csv`, with the panel value where the table has none.
- **Corner-only pieces.** A mapped line that crosses the mesh exactly at a
  vertex leaves two cells touching only at a corner. Two pieces of segment 5
  were joined this way to the nearest reach that is not higher (cell 6101 of
  segment 4, 68.6 and 31.5 m away). Segment 5 is the channel just below the
  delineated outlet, 36 of its 158 m outside the domain, so draining it
  through segment 4 is consistent with the DEM.
- **Outlet drain (decision 4 revised, cookbook 3.3).** Only the outlet
  reach's own DRN record (its cell and layer) is removed: 1 of the 164
  line-drain records. The layer-2 drain beneath it and the rest of the line
  are legitimate boundary drainage and stay.
- **Burning onto the mesh.** `TargetGrid.from_cMF` now reads
  `mesh_gridprops`. Before that fix it burned the 97 streams onto 4 cells
  of the proxy grid.

**Observations** (`<name>.obs.sfr.csv`). The outlet reach's ext-outflow,
stage and leakage are recorded, plus network totals that MF6 sums over the
boundname `network`: ext-inflow (MM runoff), evaporation, leakage,
ext-outflow, and to/from-mvr when the ponds are on. **Sign:** the `sfr`
observation is positive when the stream *loses* to the aquifer. The reach
solver negates the leakage and routes `qd = qsrc - qgwf`
(`gwf-sfr-steady.f90`); a five-reach test model closes as
500 in - 0.1 evaporated - 1.633 `sfr` = 498.267 out.

**Post-processing.**
- `sfr_observations()` collapses the csv to one row per stress period, as
  the time-weighted mean of the ATS steps, with the steady period dropped.
- `outlet_streamflow.png`/`.csv` compares the outlet discharge (m3/d, with a
  mm/yr axis) to `inputObsRo_catchment.txt` (mm/d), with NSE, r and volume
  bias.
- `budget_sfr_ts.png` shows the network budget month by month.

**Expect storm peaks below the gauge until CRR (WP5).** Only the runoff
generated on channel cells is delivered to the reaches today; off-channel
runoff still leaves the model as in the legacy MARMITES.

**Known limitation (R1).** MF6 has no unsaturated zone beneath a reach, so
streambed seepage reaches the water table directly. This is now stated in
the SFR panel help and in `docs/MARMITES_overview.md`.

## 8.17 WP4.1: the ponds on the mesh (2026-09-26)

**The bug.** The model map (`IN_000_model_map.png`) showed the LAK host
cells away from the ponds. `_build_ponds` placed the ponds with
`cMF.delr/delc/nrow/ncol`, which on a projected mesh describe the
(ncpl, 1) proxy grid of 1 m squares ("cell 1 m2" in the log). Every host
cell was an unrelated mesh cell, and only 3 ponds counted as on-channel.
This is the same class of bug as the stream burn fixed in §8.16.

**The mesh does not give one pond-scale cell per pond.** Each La Mata pond
(341-2036 m2) is covered by 15-83 cells of about 25 m2. The single seeded
cell (22-34 m2) is a stream cell for only one pond, although the mapped
stream crosses all 11 ponds inside the catchment. So a host-cell-only rule
leaves the stream bypassing the ponds.

**The fix follows the CdL design** (`cdl_gwf_model_fable_v2` §5b/§6,
converged over 45 years), on the grid the model actually uses (a
`TargetGrid`, structured or mesh):

- **Footprint:** the active cells whose centre lies inside the pond. On the
  50 m grid this is the host cell.
- **Host:** the cell holding the pond centroid (the seeded cell on the
  mesh). This is the one EMBEDDEDV connection.
- **Rim:** the mean model top over the footprint; bottom = rim - depth.
- **Stream through the pond.** The reaches in the footprint are excised
  after routing (`marmites_sfr.excise_reaches`). Every reach that drained
  into the footprint hands its flow to the lake through MVR, and the lake's
  Manning outlet spills into the reach(es) leaving it. The passage runs
  from the first entry to the **last** exit: one mapped line on the mesh
  leaves a footprint for a single cell and re-enters it. Without cutting
  that detour, the lake would have spilled half its outflow into a reach
  feeding it again, an MVR loop.
- **Ponds outside the domain.** A pond wholly outside the active domain is
  not a lake of the model. Pond 8 is dropped on the mesh and on the 50 m
  grid, so there are 11 lakes.

On La Mata's mesh (build only, not run):

- 96 reaches excised, leaving 2567 reaches;
- 15 reaches hand their flow to a pond and 11 spills go back to the stream;
- every spill path ends at the catchment outlet, through the chain
  0 -> 1 -> 2 -> 10 -> 7 -> outlet, 3 -> 4 -> 5 -> 10, 8 -> 9 -> 5 and
  6 -> 7.

**Consequences to carry into WP4:**
- Runoff on excised footprint cells is no longer delivered to SFR; it waits
  for 4.4 (MM runoff to LAK RUNOFF).
- In 4.6, f_lake = 1 over the whole footprint, not only the host cell.
- **Rim datum (user decision 2026-09-26): the LAND SURFACE**, the mean of
  `cMF.elev` over the footprint. **Bed (user rule 2026-09-27, the streams'
  rule):** the pond's total depth below the land surface is the soil depth
  of its footprint + the pond depth. So the bed sits one pond depth below
  the aquifer top, the mean of the model top over the footprint.
  Measured from the land surface, La Mata's 1.5 m of soil equalled the
  1.5 m pond, and every bed sat on the aquifer top (-0.44..+0.26 m at the
  host cell), where the seepage drains hold the water table. That is the
  switch that made the streams crawl. Build-only check on the mesh: the beds
  now sit 0.93-1.52 m below the aquifer top at the host cell, and each pond
  is 2.7-3.0 m deep from rim to bed.

## 8.18 WP4.6 + WP4.4: open water bypasses the soil column (2026-09-26)

This implements the user decision of 2026-09-09 (cookbook §4a), with the
share-based rule for streams (user, 2026-09-26).

**One per-cell surface descriptor**, built once from the network and the
ponds as written, on the MM cell list as ordered:
`f_lake + f_stream + f_soil = 1` (`clsMF6.surface_fractions`, attached to
`ctx.f_open` by the coupler).

- `f_stream = rwid * rlen / A`, the channel's own surface.
- `f_lake` spreads the pond's polygon area over its footprint: about 1 on
  the mesh, where the footprint *is* the pond, and the pond's share of its
  one host cell on the 50 m grid (0.14-0.81).
- Both are capped so no cell is more than all open water.

**MMsoil** (`_cell_step`):
- The column is computed exactly as before, per unit of soil area; its
  carried state is its own.
- Every output flux is multiplied by `f_soil`.
- Over the open fraction there is no interception, soil ET, groundwater ET,
  percolation or UZF demand. Its rain, and any groundwater seeping up from
  a seep drain, is handed over as runoff: `Ro = f_soil*Ro_col + f_open*(P + exf)`.
- The column's inputs from below (UZF rejected infiltration, the previous
  period's UZF ET) came only from under the soil fraction, so they are
  divided by `f_soil`.
- Checks: `P = Ei + Pe` still holds per cell, and the ET <= PET check holds
  as it did, now over the soil fraction's demand.

**Coupler:**
- Runoff goes to SFR INFLOW (channel cells) and, new (WP4.4), to LAK
  RUNOFF (the sum over each pond's footprint), in BOTH coupling modes. The
  iterative mode had delivered no runoff to the stream at all.
- A pond's evaporation is booked over its footprint by the pond area each
  cell holds. It used to go on the host cell alone: a 2036 m2 pond's
  evaporation on 25 m2.

**La Mata's mesh (build only, not run):**
- the soil column runs on 99.28 % of the 4.80 km2;
- 463 pond cells (`f_lake` 0.887-1.000, 11,769 m2 of the ponds' 11,812 m2);
- 2567 channel cells (`f_stream` median 0.31, 22,984 m2 of channel).

**What changes in a run:**
- MM runoff grows by the rain on 0.72 % of the catchment, and it now
  reaches the ponds.
- The PET demand line (PT + PE) covers the soil fraction only. The open
  fraction's evaporation is Eow, from MF6.
- In a pond cell, `iSsoil_pc` and `idgwt` still report the (uncounted)
  column's state.

## 8.19 The botm_l0 fix on a real EVT run (2026-10-07)

Run `20261006164429_2lay_evt_fix` (EVT route, code at 84d487f, so with
6b9defa) against `20261005220110` (EVT, before the fix) and
`20261005081712` (WEL route). All three start from the saved state
`hi_voronoi_lamata_rbth02`, and each converged after 2 spin-up cycles. The
old runs' binary outputs were overwritten, so they are compared through
their logs, CSVs, `mfsim.lst` and saved states.

**The fix took effect.** In the 181 layer-2-outcrop cells (4.15 % of the
area), MMsoil's head now equals MF6's start-of-day layer-2 head (+0.004 m on
average). Before the fix it read the soil base, 3.81 m higher on average:
on 99 % of cell-days the layer-2 head lies below the soil base.

**Under EVT the fix changes almost nothing, and that is expected.**
- The buggy curve was anchored at the soil base with the full-potential
  rate. Below that it fell along Shah's curve (Eg) and the root ramps (Tg).
- MF6 evaluated it at the real head, which gives the same rate as a curve
  anchored at the real head.
- So the bug inflated only the reserve (RATE at SURFACE) that UZF's demand
  had to leave room for.
- UZF in those cells is limited by water, not by demand, so freeing the
  reserve changed little: ETuzf +0.3 mm/yr at G2 (the one observation point
  in an L2 cell) and +0.05 mm/yr over the catchment.

**Catchment, EVT before → after the fix (mm/yr):**
- ETg 24.4 → 24.4 (Eg 15.4 → 15.3, Tg 9.1), ETuzf 13.70 → 13.75.
- Rp 165.7 and EXFg 13.0 in both.
- Every GWF listing-budget term is within 0.05 mm/yr.
- Outlet 91 mm/yr, NSE −0.59, bias −28 % in both.
- Per-cell mean heads: layer-1 cells max |Δ| 9 mm. L2-outcrop cells −3 mm
  (max 16 mm), slightly lower because UZF now takes a little more.

**L2 cells against L1 cells at the same water-table depth.** The comparison
has to be split by soil zone: 132 of the 181 are zone 3 (outcrop, Shah
'sand', extinction 0.5 m), against 107 of the 4343 L1 cells. Annual means,
area-weighted:

| group          | cells | ETg  | Eg   | Tg   | WT below land (m) | soil (m) |
|----------------|------:|-----:|-----:|-----:|------------------:|---------:|
| zones 1-2, L1  | 4236  | 25.2 | 16.0 |  9.2 | 3.37              | 0.63     |
| zones 1-2, L2  |   49  | 30.5 | 23.2 |  7.3 | 3.64              | 0.23     |
| zone 3, L1     |  107  | 10.7 |  0.0 | 10.7 | 4.40              | 0.93     |
| zone 3, L2     |  132  |  3.9 |  0.0 |  3.9 | 5.06              | 0.14     |

- **Zones 1-2:** at equal daily depth (2-5 m), L2/L1 = 0.85-0.93.
  Shallower than 2 m, the L2 cells take more, because their thin soil
  leaves more PE below it (PETuzf 278 vs 218 mm/yr).
- **Zone 3:** Eg is zero in both groups (the water table is deeper than the
  extinction depth). Tg per % of tree cover is 0.59 vs 1.26 mm/yr, with a
  0.14 m soil against 0.93 m and Q. pyrenaica cover of 0.14 % against
  0.78 %. These are surface differences, not the bug: MMsoil reads the
  real head there.

**The WEL-EVT gap of §6, cell by cell** (steady means per cell, WEL vs
EVT-fix):
- L1 cells: 24.9 vs 25.0 mm/yr. L2 cells: 195.7 vs 9.9 mm/yr (WEL median
  214).
- Contribution to the catchment mean: L1 23.9 / 24.0, L2 8.1 / 0.4. The
  whole 32.0 vs 24.4 gap is the bug under WEL.
- End-of-cycle heads, EVT-fix − WEL: L2-outcrop cells +5 cm (max 0.28 m),
  L1 cells −3.6 cm. Near the streams, the extra WEL pumping came out of the
  baseflow (SFR_OUT 59.2 vs 73.5 mm/yr) rather than out of the head.

**Convergence.**
- Extra ATS sub-steps per cycle: WEL 986/989, EVT 530/615, EVT-fix 498/667.
- In cycle 2 of both EVT runs, the same 18 SPs fail (257 vs 277 failed
  steps). 8 SPs changed their step count, net +52; SP 360 alone went from
  3 to 27.
- At a failed step, the largest change is in a layer-1 stream-reach cell:
  the 10 most frequent cells are all SFR cells. It sits in an L2-outcrop
  cell only 2-4 times.
- So the convergence trouble is in the stream-aquifer exchange on storm
  days, not in groundwater ET.

**Post-processing trap met on the way.**
- `budget_uzf/sfr/lak.csv` come from `package_budget(max_samples=120)`, a
  time-weighted mean over an even subsample of 120 of about 1000 ATS
  records. Two runs with different step lists therefore sample different
  days.
- On the EVT-fix run, subsample vs all records (m3/d):

  | term              | 120 records | all records |
  |-------------------|------------:|------------:|
  | UZF GWF           | −1316.5     | −1406.2 (−6.4 %) |
  | UZF STORAGE       | −18.1       | 0.0         |
  | LAK STORAGE       | −5.6        | 0.0         |
  | SFR EXT-OUTFLOW   | −1184       | −1201 (the log's value) |

- Reading every record costs about 40 s per file.
- `layer_storage_change` also subsamples, and it weights each record
  equally.
- `budget_terms.csv` (from the listing) and the observation exports are
  exact.

## 8.20 The WEL route and iterative mode removed, on a real run (2026-10-08)

Run `20261007132734_2lay_evt_fix` (code with dbbde15, before e8673a3: it
still read the parameter file) against the EVT-fix run `20261006164429`
(§8.19). Same start (`hi_voronoi_lamata_rbth02`). Same configuration but for
the three retired keys: the two `resolved_config.toml` differ only in
`run.mode`, `run.relax` and `et.gw_route`. Same dataset: the converter
rewrote it at 13:24 with the new GIS path, but only the provenance headers
changed, and the vegetation overlay it recomputed for that reason is
bit-identical to the cached one.

**What changed in the model.** Only the ETg wells: 4524 WEL records at
q = 0 with AUTO_FLOW_REDUCE. In the EVT-fix run their flow was exactly 0.0,
and the coupler wrote zeros to them, so removing them changes no term. The
lagged path is otherwise unchanged.

**Results.**
- Catchment budget (`budget_terms.csv`, exact): every term within
  0.003 mm/yr. The NWT comparison table is identical to 0.1 mm/yr (ETg 24.4,
  EXFg 13.0, Ro 31.4).
- Mean heads per cell: mean |Δ| 0.05 mm, max 2.7 mm (layer 1) and 0.8 mm
  (layer 2). Observation heads within 3 mm; h RMSE 1.49 → 1.48 m.
- Exact per-period SFR and LAK budgets: run means within 0.003 % (LAK
  EXT-OUTFLOW 0.3 %, i.e. 0.03 m3/d).
- Outlet flow: mean 1201.08 vs 1201.11 m3/d. 10 days differ by more than
  1 %, the worst 2008-10-29 (275 vs 241 m3/d), all storm days.
- Convergence: extra ATS sub-steps 496/642 (cycles 1/2) vs 498/667. In
  cycle 2 the same 18 SPs fail (268 vs 277 failed steps), SP 360 again takes
  27 steps, and the largest change at a failure is again in the same
  layer-1 SFR cells.

**So the runs agree, but not to the digit.** The differences are round-off
sized. They grow only where ATS takes a different discrete decision (a step
that barely fails in one run and barely converges in the other), which is
on storm days. Which round-off differs was not found. A repeat of one run
would tell whether a run is bit-reproducible at all.

**Expected differences in the outputs.**
- `_sp_plt_GWmap_ETg` is no longer drawn: it plotted the WEL term, all
  zeros on the EVT route. Groundwater ET is still mapped as `GWmap_Eg`,
  `GWmap_Tg` and `MMmap_ETg`.
- `budget_uzf/sfr/lak.csv` and `storage_change_L*.csv` differ by up to
  15 % (SFR TO-MVR 2878 vs 2491 m3/d). This is the 120-record subsample of
  §8.19: it picks other records when the step list changes. The exact
  per-period files agree.
- A SyntaxWarning from `MARMITESplot_v3.py:1849` (`'%s\%s_%s'`, the ffmpeg
  batch line) appears once when the module is recompiled after a copy. It
  is harmless.

The parameter-file removal (e830713) has not been run yet.

## 8.21 The pond outlets leak past the mover: the solver tolerance (2026-10-09)

Run `20261008114055`: LAK EXT-OUTFLOW at pond10 -2,748 m3/yr (on 339 of
365 days) and pond13 -535, although every outlet is moved to the stream at
FACTOR 1; together ~0.75 % of the outlet streamflow.

**Mechanism (MF6 6.7.0 source).**
- EXT-OUTFLOW = the outlet's discharge minus what the mover took from it
  (`gwf-lak.f90` `lak_get_external_outlet` + `lak_get_external_mover`).
- The mover moves the PREVIOUS outer iteration's provider flows: `mvr_fc`
  runs before the packages' fc (`gwf.f90`), and LAK fills `qformvr` at the
  end of `lak_solve`.
- Those flows are zeroed at the start of every time step
  (`PackageMover%ad`). The day's water reaches a pond N links down a chain
  only after N outer iterations, always from below -- hence a loss on
  almost every day, and none at a pond with nothing upstream.
- LAK's convergence check turns the change of outlet discharge into a depth
  over the pond's area per step (`lak_cc`: dqout x delt / area); the solver
  accepts it below `outer_dvclose` (`sln_package_convergence`). The Run
  panel has 0.025 m.

**On La Mata.** Each pond's daily |EXT-OUTFLOW| against that bound
(0.025 m x surface area / 1 d): no pond ever exceeds it; pond10 reaches
96 % (37.2 of 38.7 m3/d) and pond13 92 % (21.8 of 23.8). They are the two
at the end of the longest chains (4 and 3 ponds upstream) and the two with
the largest spill (877 and 769 m3/d). The other nine stay below 20 %.

**Toy model** (`code/tests/diag_lak_mover_leak.py`: 5 ponds in series on a
stream over an aquifer with La Mata's K 0.05 and Sy 0.01, 60 days, 600 to
1,000 m3/d plus storms). EXT-OUTFLOW of the last pond:

| outer_dvclose | last pond, 60 d | worst day | outer iterations/step |
|---------------|-----------------|-----------|-----------------------|
| 0.025         | -55 m3          | -45 m3/d  | 6.7                   |
| 0.001         | -12 m3          | -1.5 m3/d | 8.0                   |
| 1e-5          | -0.2 m3         | 0.0       | 10.8                  |

The loss grows down the chain (-0.1, -0.6, -3.5, -33, -55 m3 at 0.025).

**What to do.** The lever is `solver.outer_dvclose` (Run panel); 0.001 m
(the approved default) cuts pond10's bound from 38.7 to 1.5 m3/d. The cost
is more outer iterations per day (+20 % on the toy). Nothing in MF6 tightens
LAK alone: `MAXIMUM_STAGE_CHANGE` governs LAK's internal stage loop, not
this check.

## 8.22 No groundwater into the ponds: their stage, not their bed (2026-10-09)

The pond budgets per year show every pond losing to the aquifer and never
gaining. MF6 LAK, EMBEDDEDV: flow into the pond = cond x (max(head, bed) -
max(stage, bed)), cond = bedleak x wetted area (`lak_calculate_conn_exchange`;
for EMBEDDEDV `belev` is the bottom of the stage table, the bed). So the
direction is set by the aquifer head in the host cell against the pond's
stage; the conductance (bedleak 0.001 1/d) only scales it.

Run `20261008114055`, per pond and day:
- the aquifer head is BELOW the stage on every day at every pond: by
  0.8-1.6 m on average (0.1 to 2.4 m on single days);
- it is ABOVE the pond bed, by 1.2-2.1 m: the ponds are connected to the
  aquifer, not perched over it;
- the stage stays at the outlet sill (the pond rim from the DEM) all year:
  the stream passes through every pond (on-channel, FROM-MVR ~100x every
  other term), so its level cannot fall to the water table. The bed is
  3.0 m below the sill (2.5-3.0), 1.5 m below the cell top.

The pond piezometers C1-C3 (in the cells of pond1, pond3, pond10) are no
help: their "observed" series repeat the same four values every year,
2007-2012, = h0 + (-2.0, -0.5, 0.0, -1.0) m for all three (h0 from
inputObs.txt). They look like placeholders, not measurements. Taken at face
value they too put the water table 0.1-2.4 m below the pond rims, and the
simulated heads are within 0.1-0.5 m of them.

So it is conceptual before it is calibration: the model holds the ponds
full to their rims with stream water, above the water table. For a
groundwater-dependent pond the level follows the water table. Options (the
user's decision):
1. measured pond levels or depths, to check the rims and beds against;
2. if the real level is below the rim most of the year, the spill level
   (outlet sill) at that level, not at the DEM rim;
3. if the stream does not flow through the ponds all year, ponds off the
   channel -- fed by runoff, rain and groundwater, spilling only when full;
4. only then calibration (K, Sy, recharge) to raise the valley heads.
A larger bed conductance would make the losses larger, not reverse them.

## 8.23 Why 0.001 does not converge: SFR evaporation flip-flops on nearly dry reaches (2026-10-09)

Run `20261009165953` (outer_dvclose 0.001, inner_dvclose 1e-4,
outer_maximum 100; stopped by the user after 44 SPs in 5 h):
- 703 time steps, plus 266 attempts that ran all 100 outer iterations and
  failed; 7 SPs failed (8, 9, 11, 13, 25, 42, 43), the other 37 converged
  in one step;
- at the last iteration of a failed attempt, the largest head change is
  always in a layer-1 stream-reach cell: 3443 (reach 46), 3762 (4), 3269
  (16), 2932 (449), 3243 (539), 3185 (526), 1720 (410), 3417 (11), 3170
  (275), 2948 (84). The median of that change at the last iteration is
  0.32 m (the largest 41 m);
- these cells' reaches are nearly dry: depth 0 to 9 mm, often written as
  dry. Their heads sit between the streambed bottom and top, a few mm below
  the top at 3762, 3417 and 1720;
- the failures are not storm-only: SP 42 is a quiet day at 3762 (no
  infiltration, no UZF recharge, the head 1-2 mm below the streambed top,
  the reach losing 0.06 m3/d, EVT 0.14 m3/d). Every step longer than about
  0.005 d failed, so ATS ran the day in 141 steps of 0.001-0.005 d.

**Two iteration signatures, one cause.**
- A creep (SP 43 step 111, dt 0.009 d): the change at 3762 grows by about
  8 % per outer iteration, always upward (0.05 -> 1.38 m at iteration 100),
  and the residual grows with it (5.8 -> ~130). Backtracking cuts every
  step and gets the residual back to about 7-18, but the next Newton step
  raises it again: the Newton direction is not a descent direction.
- A stall (SP 42 step 3, dt 0.4 d): from iteration 5 to 100 the proposed
  change at 2948 (reach 84) stays at 6.2 mm, with the residual slowly
  rising.
At 0.025 m the stall is accepted at iteration 3, which is why 0.025
"converged".

**Mechanism (MF6 6.7.0 source; unchanged in develop on 2026-10-09).**
- `sfr_solve` takes the reach evaporation from the depth STORED by the
  previous solve: `qe = evap * calc_surface_area_wet(n, this%depth(n))`.
  That wetted area ramps from 0 to w x L over the first 1e-5 m of depth
  (`sCubicSaturation`).
- Take a reach whose inflow (upstream + runoff + mover) is below its
  potential open-water evaporation E0 x w x L (about 0.1 m3/d on La Mata),
  with the head below the streambed top. Two states alternate:
  - stored depth > 1e-5 m: the reach evaporates all its inflow, leaks
    nothing, and its new depth is 0;
  - stored depth 0: it evaporates nothing, leaks all its inflow
    (`gwf-sfr-steady.f90` shortcut, `isolve = 0`), and its new depth is
    about 1.2e-5 m.
- `sfr_fc`'s Picard loop (MAXSFRPICARD 100) never settles. After its 100
  passes it ends in the state it started from.
- `sfr_fn` then perturbs the head by DEM4 = 1e-4 m and solves once
  (`update=.false.`). It lands in the OTHER state, and hands GWF a
  derivative of +-inflow / 1e-4 m (500 m2/d for 0.05 m3/d; 775 m2/d at
  3762) where the true one is about 0.
- The corrupted Jacobian makes the Newton step at that cell crawl or
  creep. It converges only once the cell's storage Sy x A / dt outweighs
  the corrupted derivative. At 3762 that is 3.62 m2 / dt against 775 m2/d,
  i.e. dt below about 0.005 d -- exactly where ATS's steps stopped
  failing.
- The flip-flop was checked with a Python replica of `sfr_calc_steady` /
  `sfr_solve` / `sfr_fn` (`X:\tmp_claude\helpers\sfr_replica.py`), at the
  failing toy's state: stored depth 1.17e-5 m, inflow 0.05 m3/d, all of
  it evaporated, hcof = rhs = 0.

**On La Mata (run output, read-only).** Reaches in that state (inflow > 0
but below E0 x w x L, head below the streambed top) at least once in the
SP's saved steps:

| SP | steps | reaches flipping | of which hold a failing cell |
|----|------:|-----------------:|-----------------------------|
| 6, 7, 12, 14, 41 | 1-3 | 1-2 | 0-1 |
| 8  | 101 | 45 | 6 (1, 46, 159, 490, 501, 505) |
| 13 | 222 | 53 | 9 |
| 42 | 141 | 2 | both (84, 275) |
| 43 | 118 | 4 | 3 (4, 11, 275) |

Both states appear in the saved steps:
- SP 43, reaches 4, 11 and 275: 0.02-0.09 m3/d in against 0.09-0.11 of
  potential; all of it evaporated, nothing to the aquifer.
- SP 42, reach 84: nothing evaporated, all of its 0.011 m3/d to the
  aquifer.

After storms many reaches recede to a trickle, hence 45-53 reaches on
SPs 8 and 13.

**Toy model** (`code/tests/diag_sfr_evap_flipflop.py`): a valley with La
Mata's stream-cell setup (Sy 0.01, streambed 0.5 m below land, rhk 0.1,
rbth 0.2, w 1.5 m) and headwater reaches fed 0.05 m3/d against E0 0.004
m/d x 30 m2. It is driven through the API like the coupler (ATS, EVAP
written after prepare_time_step).

| EVAP | outer_dvclose | failed attempts | outer iterations |
|------|---------------|----------------:|-----------------:|
| E0 | 0.001 | 3 | 353 |
| E0 | 0.025 | 0 | 41 |
| min(E0, 0.5 x inflow / (w L)) | 0.001 | 0 | 43 |

In the fuller toy (with UZF, EVT and a 20-day run,
`X:\tmp_claude\helpers\toy_stream_conv.py`):
- **Still fails at 0.001:** without under-relaxation, without
  backtracking, at MODERATE, with a seepage drain on the stream cells, and
  without EVT. So DBD, UZF and EVT are not the cause.
- **Converges at 0.001:** without the reaches' evaporation, or without the
  trickle of runoff.

Trap for toy builders: flopy needs a second word after inner_rclose. With
STRICT, an outer iteration is accepted only when the linear solve
converges at its FIRST inner iteration (`ImsLinearBase` `testcnvg`). La
Mata's IMS has none.

**What to do (the user's decision).**
1. Cap the stream evaporation in the coupler. Write
   EVAP = min(E0, f x qin / (w L)), where qin = the step's INFLOW (MM
   runoff) + USFLOW + QFROMMVR left by the last solve (all at the SFR
   memory path), with f = 0.5.
   - A reach can then never evaporate more than half of what reaches it,
     so it stays away from the switch.
   - Over the run's 43 saved days it would remove 0.3 % of the stream
     evaporation (61 m3/d; 0.1 % at f = 0.8).
   - It keeps outer_dvclose 0.001 and with it the pond-leak fix of §8.21.
   - It needs a config field (with front-end help), the coupler code and
     tests.
2. Back to 0.025, accepting the pond leak (§8.21), or an intermediate
   value. On La Mata the stall reached 6 mm, so 0.005 is not safe.
3. Report it to MODFLOW 6: compute the evaporation from the depth being
   solved, or do sfr_fn's perturbation from the converged state. The toy
   reproduces it in seconds.

### 8.23.1 Implemented: the cap on the stream evaporation (2026-10-09)

The user chose option 1, with f = 0.5.

**Configuration.** `sfr.evap_inflow_fraction`, default 0.5, on panel 4 in
the SFR form, with help text in `app/lib/schema.py`.
- Valid strictly between 0 and 1. At 1 a reach can still evaporate all it
  receives, which is the flip-flop again.
- `tests/run_lamata_mf6.py` hands it to the builder as `b.sfr_evap_frac`.
  The coupler reads it there; the builder's default is 0.5.

**Coupler** (`marmites_coupler.py`).
- `_bind` also binds SFR `USFLOW`, `QFROMMVR` (only with MOVER), `LENGTH`
  and `WIDTH`, and keeps w x L per reach. If one of USFLOW, LENGTH or
  WIDTH is missing it warns, and EVAP is written uncapped.
- `_keep_sfr_inflow`, after every step of `_advance` (a period's first step
  and its ATS sub-steps), copies USFLOW + QFROMMVR. It is a copy because
  MF6 zeroes USFLOW when it advances the next step (`sfr_ad`, inside
  prepare_solve), which the older prepare_solve route runs before the
  coupler writes.
- `_cap_sfr_evap`, inside `_write_openwater_evap`, writes
  EVAP = min(Eo, f x qin / (w L)), with qin = this step's INFLOW (the
  runoff `_write_runoff` has just written) + the copy.
  - The reach's own groundwater discharge is not counted: a reach the
    aquifer feeds does not flip-flop, and one it stops feeding must not be
    left evaporating it.
  - The first step of a run counts its runoff only.
- **Log lines:**
  - at bind: `coupler: stream evaporation capped at 0.5 x each reach's
    inflow (sfr.evap_inflow_fraction)`;
  - at the end: `stream evaporation: capped at 0.5 x the inflow on N
    reach-step(s), P % of those written`.

**Checked.**
- A scratch La Mata build (initialize only): all six SFR arrays bind at the
  coupler's addresses, 660 reaches each, with ponds and mover. w x L from
  MF6 equals the package file exactly (median 32.9 m2).
- `tests/test_sfr_evap_cap.py` drives the toy through libmf6 with the
  coupler's own methods at outer_dvclose 0.001:
  - capped: 0 failed attempts, and no reach evaporates more than half its
    inflow;
  - uncapped: it still fails. This pins the MF6 behaviour; if that test
    starts failing, MF6 has changed and the cap may no longer be needed.
- `tests/test_openwater_evap.py`: the cap, the first step, the missing
  areas, the copy after each step. `tests/test_config.py`: the default and
  the range.

**On the next run (outer_dvclose 0.001).**
- Expect few or no failed attempts at the stream cells, and far fewer ATS
  sub-steps.
- Expect the stream evaporation within about 0.3 % of before.
- A failure left at a stream cell would mean its inflow fell by more than
  half within one step, because the cap uses the previous step's upstream
  flow. The end-of-run line says how often the cap acted.

### 8.23.2 First run with the cap: 20261009234407 (2026-10-10)

**Setup.** 500 SPs (2008-05-31 to 2009-10-12), outer_dvclose 0.001,
inner_dvclose 1e-4, sfr.evap_inflow_fraction 0.5. Three spin-up cycles
(1h58, 1h32, 1h40), converged at cycle 3 (mean |dWT| 0.094, then 0.019 m).
Finished at 05:02, 5h18 with post-processing.

**Convergence.** From the solution check and the last cycle's `mfsim.lst`:
- extra ATS sub-steps per cycle: 685, 406, 451 over 500 days. Run
  `20261008114055` at 0.025 without the cap had 526 and 623 over 365 days,
  so the tight run takes fewer sub-steps per day;
- cumulative discrepancy 0.01 %;
- last cycle: 951 time steps, 211 failed attempts, 51 SPs needing more
  than one step. The uncapped run at 0.001 had 266 failed attempts in its
  first 44 SPs.

Every failure still ends at a layer-1 stream cell, by two routes:
1. **The cap does not hold: inflow collapsing within a step.** SPs
   496-499, the October 2009 storm recession:
   - reach 441's inflow falls from 284 to 0.007 m3/d within SP 498, and
     reach 140's from 0.36 to 0.006;
   - reaches 376, 644 and 323 evaporate all of it (evap/qin = 1.00);
   - this is the flip-flop again, with the creep signature (60 m at cell
     1711);
   - it is the case §8.23.1 named: the cap uses the previous step's
     upstream flow.
2. **The cap holds and the reach sits on the switch.** SP 477 (48 failed
   attempts, reach 314), SPs 458-459 (reaches 83, 86), and the summer
   days 404-447:
   - evaporation is exactly 0.5 x the inflow, and the rest leaks;
   - the depth stays at 0.010 mm, the edge of MF6's 1e-5 m wetted-area
     ramp;
   - signature: oscillation at 1-17 mm, 50-90 sign changes in 100
     iterations, around the tolerance.

**Open.** At those reaches the saved head is 0.12-0.27 m ABOVE the
streambed top while SFR reports them losing (checked in both budget files,
user node numbers in both). The replica of `sfr_calc_steady` reproduces
the loss only with the head below the bed. To pin down.

**Results.**
- **Pond leak:** LAK EXT-OUTFLOW went from -9.1 m3/d (0.025) to -0.1 m3/d.
- **Outlet** over the 329 gauged days the two runs share: NSE -0.54
  against -0.59, bias -37 % against -28 %, 162 m3/d less simulated flow.
  The 500-day NSE of -3.98 comes from the storm of 7-9 October 2009, first
  reached by this run: simulated 104,048 and 36,343 m3/d against 9,493 and
  5,539 observed. That is runoff generation, not convergence.
- **Stream evaporation:** 2.40 mm/yr (2.55 at 0.025 over a different
  window). UZF rejected 21.4 % of the applied percolation.
