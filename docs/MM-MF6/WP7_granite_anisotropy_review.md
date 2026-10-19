# Vertical against horizontal hydraulic conductivity in saprolite and weathered granite

Literature review for the La Mata model (WP7), 2026-10-19, at the user's
request. Companion of `WP7_granite_porosity_specific_retention_review.md`.

## What the model uses today

- `[layers] k33_as_ratio = true`, `[layers.k33] value = 2.0`: **Kv = Kh / 2
  in both layers** (MF6 NPF k33 = k / 2).
- It comes from the legacy MODFLOW-NWT model of La Mata: UPW `layvka = 1`,
  `vka = 2` in both layers. The LPF block of the same parameter file (read
  only by non-NWT versions) had `vka = 10` (Kv = Kh / 10).
- UZF's vertical saturated conductivity is the layer's k33 times
  `uzf.vks_scale` (`uzf.vks_from = layer`).
- The calibration adjusts K; Kv follows it at the fixed ratio.

## The fissured (stratiform fractured) layer: Kh about 10 Kv

The anisotropy of the weathering profile of granites is that of its
fissured layer, and it has a geological origin: the weathering-induced
fractures are subhorizontal, parallel to the paleo-weathering surface.

- Lachassagne, Dewandel & Wyns (2021), the reference review: in granites the
  subhorizontal fractures are more permeable than the subvertical ones,
  "Kxy ≈ 10 Kz" (citing Maréchal et al. 2004); the horizontal fracture set
  is on average 10 times more permeable, and more numerous, than the
  subvertical joints.
- Maréchal, Wyns, Lachassagne & Subrahmanyam (2004), Maheshwaram granite
  watershed (India), pumping tests with observation wells: the horizontal
  permeability is 2 to 30 times the vertical one, systematically.
- Maréchal, Dewandel & Subrahmanyam (2004), same site, hydraulic tests at
  several scales: horizontal/vertical ratio close to 10 for the
  weathered-fractured layer (Neuman's anisotropic unconfined solution).

## The saprolite: weakly anisotropic, little measured

No field measurement of Kh/Kv in granite saprolite was found. What there is:

- In the conceptual model (Dewandel et al. 2006; Lachassagne et al. 2021)
  the saprolite is the capacitive, low-permeability layer of the composite
  aquifer: K geometric mean 2·10⁻⁶ m/s (0.17 m/d), range 10⁻⁷ to
  5·10⁻⁶ m/s (0.009-0.43 m/d) in Indian granite; lower in the basal
  LAMINATED layer, the one part of the saprolite that keeps subhorizontal
  fissures -- and so the one likely to be anisotropic like the fissured
  layer below.
- In the isalterite the rock structure is preserved (little volume change),
  so the parent rock's fabric carries into the saprolite: in Piedmont
  saprolite (Virginia) a preserved steep foliation inhibited vertical flow
  in the laboratory; in granite, without foliation, there is no such fabric
  except the subhorizontal fissuring of the laminated base.
- Laboratory limits for unfractured earth materials: homogeneous granular
  materials hardly exceed Kh/Kv ≈ 2.5; cohesive soils 0.7-4, central
  value ≈ 1.6 (Guo 2015). Compacted completely decomposed granite: kh/kv
  above 3 at a void ratio of 0.70 (cube permeameter).
- Field values are controlled by macropores, root channels and relict
  structures that cores miss (Vepraskas et al. 1991: channels carried ~93 %
  of the flow in a saprolite).

## For La Mata

| Layer | Literature Kh/Kv | Model now |
|---|---|---|
| L1 saprolite | ~1 to 3 (no fabric); higher in the laminated base | 2 |
| L2 fissured | ~10 (2 to 30) | 2 |

The saprolite's ratio of 2 is in the range; the fissured layer's is low:
the literature points to about 10 for a granite fissured layer.

## Sources

- Lachassagne P., Dewandel B., Wyns R. (2021) Review: Hydrogeology of
  weathered crystalline/hard-rock aquifers -- guidelines for the operational
  survey and management of their groundwater resources. Hydrogeology
  Journal 29:2561-2594. https://link.springer.com/article/10.1007/s10040-021-02339-7
- Maréchal J.C., Wyns R., Lachassagne P., Subrahmanyam K. (2004) Vertical
  anisotropy of hydraulic conductivity in the fissured layer of hard-rock
  aquifers due to the geological structure of weathering profiles. Journal
  of the Geological Society of India 63:545-550.
  https://www.geosocindia.org/index.php/jgsi/article/download/82762/63884/142380
- Maréchal J.C., Dewandel B., Subrahmanyam K. (2004) Use of hydraulic tests
  at different scales to characterize fracture network properties in the
  weathered-fractured layer of a hard rock aquifer. Water Resources
  Research 40, W11508. https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/2004WR003137
- Dewandel B., Lachassagne P., Wyns R., Maréchal J.C., Krishnamurthy N.S.
  (2006) A generalized 3-D geological and hydrogeological conceptual model
  of granite aquifers controlled by single or multiphase weathering.
  Journal of Hydrology 330:260-284.
  https://www.sciencedirect.com/science/article/abs/pii/S0022169406001673
- Guo (2015) Limit of anisotropic hydraulic conductivity ratio of
  homogeneous granular materials. Vadose Zone Journal.
  https://acsess.onlinelibrary.wiley.com/doi/10.2136/vzj2015.01.0015
- Vepraskas M.J. et al. (1991) Hydraulic conductivity of saprolite as
  determined by channels and porous groundmass. SSSA Journal 55.
  https://acsess.onlinelibrary.wiley.com/doi/pdfdirect/10.2136/sssaj1991.03615995005500040006x
- The influence of saprolite composition, texture and structure on
  groundwater flow in the Piedmont Province, Henry County, Virginia (Old
  Dominion University thesis, 1994). https://digitalcommons.odu.edu/oeas_etds/266/
- Kuna Raj J. (2021) Saturated hydraulic conductivity of earth materials in
  the weathering profile over a porphyritic biotite granite, Kuala
  Lumpur-Karak Highway. Bulletin of the Geological Society of Malaysia 71.
  https://gsm.org.my/articles/702001-101903/
