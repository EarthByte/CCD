# Step 2 — Reconstructed ocean-basin areas

Produces the time-dependent Atlantic / Pacific / Indian area fractions that step 3
uses to weight the regional CCD curves into a global mean.

The whole step is one notebook, `ocean_basin_areas.ipynb`.

**This step is not run by `run_ccd_core.py`.** It is a notebook, and its result
reaches step 3 through `data/ccds_and_oceanbasin_fractions.xlsx` by a manual copy
(see *Getting the weights into step 3* below). Re-running the core workflow does
not regenerate the weights.

## Plate model

> **Müller et al. (2019)**, obtained through GPlately's `PlateModelManager`
> (`MODEL_NAME = "Muller2019"`), with continental polygons
> `muller2019/ContPolygons/Global_PresentDay_ContPolygons_2019_v1.shp`.

The rest of the workflow (steps 8–10, Figure 3) uses Alfonso et al. (2025, mod.
Clennett–Müller). That difference does not matter here: the two models are identical
in the shape of the large ocean basins, and the Alfonso modifications concern island
arcs in the NE Pacific, which have no bearing on basin-scale areas. The model is
recorded for provenance, not as a caveat.

## Method

Everything is computed on the unit sphere (radius 6371 km), never on a projection,
so there are no antimeridian or area-distortion artefacts.

1. **Grid.** A 1° × 1° longitude/latitude grid of nodes; each node carries its own
   spherical area element `dA = R² cos(lat) dlon dlat`.
2. **Dynamic land mask.** The Müller 2019 continental polygons are reconstructed to
   each time step and rasterised onto the grid. Ocean = every node not inside a
   reconstructed continental polygon. Antarctica therefore forms the southern
   boundary naturally; no separate Southern Ocean cut is imposed.
3. **Basin boundaries.** Ocean nodes are split between basins by three dividing
   "curtains", each a mask of half-width 1.5° (2.0° for Tasmania) about a meridian:
   - **Africa** (Atlantic | Indian), near 20°E;
   - **West Pacific** (Indian | Pacific), near 136°E, from 30°S up to a
     time-dependent northern cap (12° at ≤10 Ma, 20° at ≤30 Ma, 35° earlier);
   - **Tasmania**, closing the Indian–Pacific boundary in the south, with an
     Antarctic Peninsula anchor at (58.5°W, 63.5°S).

   The present-day meridians are inferred once from Natural Earth ocean polygons,
   then each is pinned to an anchor point whose plate ID is resolved from the static
   polygons, so the boundaries **move with the plates** rather than staying fixed in
   longitude.
4. **Integration.** Basin area at each time = sum of the node areas classified into
   that basin.
5. **Times.** 0–70 Ma at 5 Myr, then spline-interpolated to 0.5 Myr.

## Outputs (`ocean_basin_area_outputs/`)

| file | contents |
|---|---|
| `ocean_basin_areas_0_70Ma_million_km2.csv` | the 5 Myr areas — **this is the published set** (Dataset S1, sheet "Ocean basin areas") |
| `ocean_basin_areas_0_70Ma.csv` | same, in km² |
| `ocean_basin_areas_0_70Ma_0.5my_million_km2.csv` | spline interpolation to 0.5 Myr |
| `basins_*Ma.png`, `areas_timeseries*.png` | per-time basin maps and the area time series |

Older generations (`*_0_65Ma*`, `*_v6h_*`, and the `basin_areas/` folder) are
superseded runs, kept only for reference.

## Getting the weights into step 3

The published fractions in `data/ccds_and_oceanbasin_fractions.xlsx` (columns 7–9)
and in Dataset S1 were derived from the 5 Myr areas above, linearly interpolated to
1 Myr, with the Indian Ocean handled according to how many basins still have CCD
data at that age:

- **0–23 Ma, all three basins:** ordinary area fractions,
  `f_i = A_i / (A_Atl + A_Pac + A_Ind)`.
- **24–52 Ma, Indian CCD record ended:** the Indian Ocean's *area* is split evenly
  between the two remaining basins rather than discarded, i.e.
  `f_Atl = (A_Atl + A_Ind/2) / (A_Atl + A_Pac + A_Ind)` and likewise for the Pacific.
  Verified to reproduce the published fractions to four decimal places.

  **What this assumption actually is.** Splitting the Indian area 50:50 is
  algebraically identical to imputing the unsampled Indian CCD as the midpoint of
  the Atlantic and Pacific curves:

  ```
  (A_Atl + A_Ind/2)·CCD_Atl + (A_Pac + A_Ind/2)·CCD_Pac
  ------------------------------------------------------  ==
                A_Atl + A_Pac + A_Ind

      A_Atl·CCD_Atl + A_Pac·CCD_Pac + A_Ind·(CCD_Atl + CCD_Pac)/2
      -----------------------------------------------------------
                        A_Atl + A_Pac + A_Ind
  ```

  (verified numerically: the two agree to 0.000000 m at every timestep). Stating it
  as an imputation rather than as an area manipulation makes the commitment
  explicit — the unsampled basin is assigned the mean behaviour of its sampled
  neighbours, which is the least-committed choice available: it assumes the Indian
  resembled neither the Atlantic nor the Pacific preferentially.

  **Why not simply drop it.** The Indian Ocean does not stop existing when its CCD
  record ends; it becomes unsampled. Discarding its area would assert that the world
  ocean is only as large as the part we can measure. The choice is not cosmetic:
  against dropping the area, the Atlantic weight over this interval is 0.31 rather
  than 0.26, and because the Atlantic and Pacific CCDs differ by ~600 m here, the
  global curve differs by a mean of 27 m and up to 123 m at individual timesteps.
- **53 Ma and older, Atlantic only:** weight 1.0.

Step 3 renormalises whatever weights it is given over the basins present at each
timestep, so the renormalisation is idempotent and applying it twice is harmless.

**Beware:** because the join is manual, the fractions in the xlsx do *not*
automatically follow a re-run of this notebook. If the areas are recomputed, the
xlsx columns must be rebuilt by the rule above, and the 5 Myr areas sheet of
Dataset S1 updated to match.
