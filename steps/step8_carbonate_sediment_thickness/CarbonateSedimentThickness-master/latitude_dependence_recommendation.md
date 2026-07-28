# Latitude-dependent reduction of predicted carbonate thickness

## 1. Why the current workflow over-predicts high-latitude carbonate

In `carbonate_sediment_thickness.py` the decompacted carbonate sediment rate at every reconstruction step is computed as:

```
max_carbonate_decompacted_sediment_rate = 10 * max_carbonate_decomp_sed_rate_cm_per_ky_curve(reconstruction_time)
carbonate_decompacted_sediment_rate = (
    max_carbonate_decompacted_sediment_rate *
    min(bathymetry_at_reconstruction_time - ccd_depth_at_reconstruction_time, CCD_DISSOLUTION_DISTANCE) /
    CCD_DISSOLUTION_DISTANCE
)
```

The only spatial control is bathymetry relative to the (globally uniform) CCD. There is no dependence on latitude, surface temperature, carbonate ion saturation state, or productivity. Because the Arctic crust is young enough (and therefore shallow enough) to sit above the prescribed CCD for much of the past 170 Myr, the model deposits thick pelagic carbonate there — which is geologically wrong. In reality, pelagic CaCO₃ accumulation drops sharply poleward of roughly 50–60° because:

- Surface waters are cold and have a low calcite/aragonite saturation state — calcifying plankton (coccolithophores, planktonic foraminifera, pteropods) are scarce there and what they produce dissolves more readily.
- The CCD/lysocline shoals strongly at high latitudes (the global average curve underestimates this shoaling).
- The Arctic is also poorly ventilated, biologically siliceous- and clay-dominated, and receives large IRD/terrigenous input that suppresses the carbonate fraction.

## 2. Recommended fix: a latitude-dependent multiplier

The cleanest, lowest-risk change is to multiply the per-step `carbonate_decompacted_sediment_rate` by a factor `f(|φ|)` ∈ [0, 1] where φ is the **paleo-latitude of the reconstructed point at `reconstruction_time`** (not the present-day latitude — this matters because Arctic crust that originated at lower latitudes should still receive carbonate from its low-latitude history). The paleo-latitude is already available in the loop: `reconstructed_points[index].to_lat_lon()[0]`.

Because you don't want any reduction at low/mid latitudes, the function should be exactly 1 below some threshold and roll off smoothly above it. I recommend a **half-cosine taper** with two parameters:

```
def latitude_factor(lat_deg, lat_full=60.0, lat_zero=80.0):
    """
    Latitude weighting for pelagic carbonate accumulation.
        |lat| <= lat_full   -> 1   (no reduction; tropics through mid-latitudes)
        lat_full < |lat| < lat_zero -> smooth cosine roll-off
        |lat| >= lat_zero   -> 0   (polar; no pelagic CaCO3)
    """
    abs_lat = abs(lat_deg)
    if abs_lat <= lat_full:
        return 1.0
    if abs_lat >= lat_zero:
        return 0.0
    x = (abs_lat - lat_full) / (lat_zero - lat_full)
    return 0.5 * (1.0 + math.cos(math.pi * x))
```

With `lat_full=60°`, `lat_zero=80°` this gives 1.0 at 60°, 0.85 at 65°, 0.50 at 70°, 0.15 at 75°, 0.0 at 80° — a smooth taper that kills Arctic carbonate but leaves the productive belt untouched.

Why a half-cosine rather than a sigmoid or piecewise-linear:

- It is exactly 1 below `lat_full` and exactly 0 above `lat_zero`, satisfying your "no effect at low/mid latitudes" requirement (a logistic always shaves a bit off everywhere).
- It is C¹ continuous, so no kinks propagate into the gridded output as artefacts.
- It has only two physically meaningful knobs (`lat_full`, `lat_zero`) that are easy to justify and to vary in sensitivity tests.

Defensible parameter choices from modern observations:

- `lat_full = 55–60°`: poleward edge of the subpolar productive belt; coccolith and foram fluxes from sediment-trap compilations (e.g. Honjo et al., 2008) are still close to their mid-latitude values here.
- `lat_zero = 75–80°`: by 80° both calcite and aragonite are seasonally undersaturated in surface waters and pelagic carbonate burial in the Arctic Ocean (e.g. Stein 2008, *Arctic Ocean Sediments*) is essentially nil.

A useful sensitivity sweep would be (lat_full, lat_zero) = (50, 70), (60, 80), (65, 85).

## 3. Where to insert it

In `calc_carbonate_decompacted_sediment_thickness` (around lines 519–525):

```python
# paleo-latitude of this point at reconstruction_time (carbonate ref. frame is mantle/abs.)
paleo_lat, _ = reconstructed_points[index].to_lat_lon()
lat_weight = latitude_factor(paleo_lat)

max_carbonate_decompacted_sediment_rate = 10 * max_carbonate_decomp_sed_rate_cm_per_ky_curve(reconstruction_time)
carbonate_decompacted_sediment_rate = (
    max_carbonate_decompacted_sediment_rate *
    min(bathymetry_at_reconstruction_time - ccd_depth_at_reconstruction_time, CCD_DISSOLUTION_DISTANCE) /
    CCD_DISSOLUTION_DISTANCE *
    lat_weight
)
```

`reconstructed_points` are returned by `topological_snapshot.reconstruct_points(...)` (or the older `reconstruct_geometry` path) and live in the carbonate reference frame, which — since `carbonate_anchor_plate_id = 701` and the rotation model is the Müller 2019 mantle-reference frame — corresponds to paleo-positions in absolute coordinates. So `reconstructed_points[index].to_lat_lon()[0]` is the true paleo-latitude.

Also expose `lat_full` and `lat_zero` as module-level constants (next to `CCD_DISSOLUTION_DISTANCE`) so they can be tuned without editing the inner loop.

## 4. Caveats and possible extensions

- This is a *kinematic* fix, not a biogeochemical one. It implicitly bundles together the temperature/productivity/saturation-state controls into one latitude proxy. If you later want a more physical model, the natural next step is to make the multiplier a function of paleo-SST (from a paleoclimate model) or of latitude-dependent CCD shoaling, but the half-cosine in latitude captures ~80% of the signal at trivial computational cost.
- The threshold should probably *not* be made time-dependent in this paper unless you also vary the CCD curve regionally; coupling two poorly-constrained latitude dependences invites over-fitting.
- Equatorial upwelling zones (carbonate-poor due to high silica productivity and shoaled lysocline along the equator) are a separate, smaller artefact that is not addressed by a `|lat|` taper — flag this in the methods if reviewers raise it.
