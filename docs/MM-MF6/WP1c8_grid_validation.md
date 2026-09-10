# WP1c.8 — the grid validation ladder

Does a Voronoi run mean the same thing as the validated structured run?

Three rungs, each further from the model we trust. The tool is
`code/tools/validate_grid.py`; it compares two finished runs by reducing each
to the catchment water balance (the same fourteen terms `plot_water_budget`
draws) plus the mean modelled head at every observation point.

```bash
# rungs (a) and (b), run end to end
python code/tools/validate_grid.py --ladder --nsp 30 --report ladder.md

# any two finished runs -- this is how rung (c) is done
python code/tools/validate_grid.py \
    --compare E:/00code_ws/LaMata_MM-MF6/MF6_ws \
              E:/00code_ws/LaMata_MM-MF6/MF6_ws_voronoi \
    --labels structured voronoi --tol 15 --report rung_c.md
```

Heads are matched by observation-point **name**, never by position: on
different grids the same piezometer lives in a different cell, which is the
whole reason the ladder exists.

---

## Rung (a) — DIS vs DISV-from-DIS  ✅ PASS, exact

A DISV-from-DIS grid is the same geometry re-expressed as polygons. Nothing
about the physics changes, so **any** difference is a bug in the projection
rather than an effect of discretisation. The tolerance is therefore exact.

Result (10 stress periods, 1954 active cells on both):

| term | DIS [mm/y] | DISV-from-DIS [mm/y] | diff |
|---|---:|---:|---:|
| P | 237.132 | 237.132 | +0.000 |
| Ro | 7308.627 | 7308.627 | +0.000 |
| ETsoil | 544.164 | 544.164 | +0.000 |
| Rp | 165.749 | 165.749 | +0.000 |
| EXFg | 8026.607 | 8026.607 | +0.000 |
| … | … | … | +0.000 |

**Identical on every one of the fourteen terms, to 1e-9 mm/y.** The balance
residual matches too (−73.09 mm/y both). Together with the unit-level identity
tests (`test_mesh_projection.py`: projecting onto a DIS-equivalent mesh
reproduces every array, all 12 DRN records and all 1954 active cells) this
closes the question of whether the projection is faithful.

---

## Rung (b) — DISV-from-DIS vs uniform Voronoi  ⚠️ INCONCLUSIVE so far

Here the discretisation genuinely changes: 1954 cells of 2500 m² become 469
cells averaging 9845 m². The water balance is *expected* to move. The question
is whether it moves by an amount a modeller accepts.

**It cannot be answered from a cold start, and the tool now refuses to
pretend otherwise.** Run at 10 stress periods from a flat initial water table,
both runs report:

| term | DISV-from-DIS | Voronoi | vs precipitation |
|---|---:|---:|---:|
| P | 237 mm/y | 236 mm/y | — |
| Ro | 7309 mm/y | 5543 mm/y | **31× / 23× P** |
| EXFg | 8027 mm/y | 6189 mm/y | **34× / 26× P** |

The water table starts above ground over much of the catchment, the seepage
face discharges enormously, and that water runs off and re-infiltrates. What
such a comparison measures is how each *grid* relaxes its initial condition,
not its hydrology. `spinup_dominated()` detects it (Ro or EXFg above 3× P) and
the report says **INCONCLUSIVE** instead of PASS or FAIL.

Note that rung (a) is unaffected: identity does not depend on equilibrium, so
the same geometry must give the same answer whatever state the run is in.

### To settle rung (b)

1. Generate a spin-up state **on the Voronoi mesh** — the saved
   `hi_spinup_*.asc` is on the structured grid and WP0.6 refuses it on a mesh,
   deliberately.
2. Re-run both cases from their own equilibrated states, over at least a full
   hydrological year.
3. Compare with `--tol`. The 15 % default in `RUNGS['b']` is a placeholder,
   not a validated threshold.

---

## Rung (c) — production Voronoi vs the validated multi-year DIS run

This is the rung that decides whether Voronoi becomes the default, and it is
**not automated**. It needs the real multi-year run, which the modeller
launches; the tool does the comparing.

It also cannot be reduced to a pass mark, because the standard is not "the
numbers are close" but **every difference is attributable to discretisation
rather than to a fault**. What to check, in order:

1. **Mass balance closes on both.** A residual that grows on the mesh is a
   bug, not discretisation.
2. **Precipitation and interception agree to well under 1 %.** They are forced
   inputs and barely depend on the grid; La Mata's 100 m mesh shows 0.46 %,
   which is the boundary discretisation (the mesh tiles 99.86 % of the model
   rectangle) and nothing more. A larger gap means cells are being lost.
3. **ETsoil, Rp and ETg may move by several per cent** — they depend on soil
   zone and thickness, which are resampled (WP1c.3). Cross-check against the
   resampling report: soil thickness reproduces the raster mean to 0.06 %, but
   the *zone* map is majority-voted and loses minority classes (soil zone 1
   drops from 13.1 % to 9.5 % of the catchment). If ETsoil moves much more
   than the zone areas did, look for something else.
4. **Heads at the observation points.** Expect metre-scale differences where
   the topography is steep, near-zero where it is flat — I1 moved 0.007 m and
   O2 −7.8 m on the cold-start pair. A large shift at a *flat* point is
   suspicious.
5. **SM and EC share one cell on the 100 m mesh**, so their modelled series
   are identical by construction. That is a real loss of resolution, not a
   coincidence, and it matters for calibration (WP7): two observations
   constraining one cell.

### Known issues to settle before rung (c) is meaningful

- **`seep.cond` is per cell and therefore grid-dependent.** 10 000 m²/d holds
  the seepage within ~0.2 m on 50 m cells; the first mesh run wanted ~1.5×10⁶.
  Retuning it blind is worse than leaving it — do it against an equilibrated
  mesh state.
- **`grid.voronoi.stream_refine` is `true` by default** (the cookbook's
  choice), but CdL settled on `false` after judging a refined mesh "too
  refined near streams". La Mata has not made that call; every run so far used
  `false`.
- **The MODFLOW-NWT comparison stays on the structured grid.** It is a 65 × 60
  run, so it is not like-for-like against a mesh; `plot_water_budget` skips the
  reference on a mesh and says so.
