# Evolution and drivers of the global carbonate compensation depth over the past 170 million years — model grids and animations

Dutkiewicz, A., and Müller, R.D.
EarthByte Group, School of Geosciences, The University of Sydney

**DOI:** https://doi.org/10.5281/zenodo.22785465

This archive holds the reconstructed grids and animations behind the paper. The
code that produces them, the derived curves, the figures and the supplementary
datasets are in the companion repository:

> **Code:** https://github.com/EarthByte/CCD

The two are meant to be used together: clone the repository, download this
archive, and point `pipeline_carbon/config.sh` at the unpacked folders. Everything
in the repository that does not need a grid runs without this archive.

---

## Contents

| Archive | Contents | Files | Unpacked |
|---------|----------|------:|---------:|
| `carbonate_sediment_thickness_mean.tar.gz` | Carbonate sediment thickness, central sedimentation-rate scenario | 513 | 451 MB |
| `carbonate_sediment_thickness_min.tar.gz` | Carbonate sediment thickness, minimum sedimentation-rate scenario | 513 | 451 MB |
| `carbonate_sediment_thickness_max.tar.gz` | Carbonate sediment thickness, maximum sedimentation-rate scenario | 513 | 451 MB |
| `paleobathymetry.tar.gz` | Merged paleobathymetry driving the thickness model | 171 | 396 MB |
| `continental_masks.tar.gz` | Reconstructed continental masks, used to exclude continental crust | 171 | 11 MB |
| `plate_model.tar.gz` | The plate model the reconstructions run on | — | 149 MB |
| `present_day_validation_grids.tar.gz` | Observed present-day carbonate thickness, min / mean / max | 3 | 5 MB |
| `animations.tar.gz` | The two supplementary animations, each rendered forward and backward in time | 4 | 44 MB |
| `MANIFEST.txt` | Every file in the archive with its size and SHA-256 checksum | — | — |

`MANIFEST.txt` is the authority on what a download should contain. Verify a
download against it before using the grids:

```bash
sha256sum -c MANIFEST.txt          # Linux
shasum -a 256 -c MANIFEST.txt      # macOS
```

---

## Carbonate sediment thickness grids

These are the archive's principal product: the compacted and decompacted thickness
of pelagic carbonate on the ocean floor, reconstructed at every million years from
170 Ma to the present, under three sedimentation-rate scenarios.

Each of the three scenario archives contains 171 time steps of three grids:

| File pattern | Quantity | Units |
|--------------|----------|-------|
| `compacted_sediment_thickness_0.25_<age>.nc` | Carbonate thickness after compaction — the thickness a drill core would measure | m |
| `decompacted_sediment_thickness_0.25_<age>.nc` | The same carbonate at its original, uncompacted thickness | m |
| `deposition_mask_0.25_<age>.nc` | Where carbonate is deposited at that reconstruction time: seafloor above the CCD and younger than the onset of pelagic carbonate deposition | 1 = deposition, 0 = none |

`<age>` is the reconstruction age in Ma, an integer from 0 to 170. All grids are
global NetCDF (CF-1.7), 0.25° by 0.25°, 1441 × 721 nodes, longitude −180° to 180°
and latitude −90° to 90°, with the data in the variable `z` as 32-bit floats and
land and undeposited ocean floor as NaN. They are in **reconstructed (paleo)
coordinates**: each grid is already in the frame of its own reconstruction time and
must not be rotated again.

### The three scenarios

All three use the same CCD — the hybrid curve of this study — and the same
paleobathymetry. They differ only in the maximum carbonate sedimentation-rate
curve applied above the CCD, which is the dominant uncertainty in converting a CCD
into a thickness:

| Scenario | Sedimentation-rate curve | Directory in the repository |
|----------|--------------------------|------------------------------|
| minimum | `sed_rate_min.txt` | `carbonate_sed_thickness_min_DM2026` |
| central | `sed_rate_best.txt` | `carbonate_sed_thickness_DM2026` |
| maximum | `sed_rate_max.txt` | `carbonate_sed_thickness_max_DM2026` |

The central scenario is the one the paper's figures, volumes and carbon budget are
built from. The minimum and maximum scenarios bound it, and are what the
ground-truth comparison against 38 DSDP/ODP basement sites (Figs S2–S4) and the
end-member carbon-flux scenarios (Fig. S5) use. **Use all three if you are
propagating uncertainty; use the central one if you want a single best estimate.**

---

## Paleobathymetry grids

`paleobathymetry_<age>Ma.nc`, 171 grids from 0 to 170 Ma at 1 Myr, global, in
metres below sea level and in reconstructed coordinates. These set the depth that
the CCD is compared against at every node, so the thickness grids cannot be
reproduced without them.

## Continental masks

`Alfonso2024_continent_mask_<age>.00Ma.nc`, 171 grids matching the reconstruction
times. Non-zero marks continental crust, which carries no pelagic carbonate and is
excluded from every volume and area statistic. These also supply the grey
continental background of the maps in Figure 3, where the rigid continental
polygons alone would leave the deforming regions blank.

## Plate model

The plate model the reconstructions run on — rotations, topological plate
boundaries, deforming meshes, static polygons, coastlines, continental polygons,
isochrons and terranes. It is included so the archive is self-contained and the
reconstructions can be repeated against exactly the model that produced them,
rather than against whichever version is current.

## Present-day validation grids

Observed present-day compacted carbonate thickness on the same 0.25° grid, as
minimum, mean and maximum estimates. These are what the modelled present-day
thickness is tested against.

## Animations

| File | Shows |
|------|-------|
| `carbonate_thickness_back_in_time.mp4` | Compacted carbonate sediment thickness, present day back to 170 Ma |
| `carbonate_thickness_forward_in_time.mp4` | The same, 170 Ma forward to the present |
| `paleobathymetry_back_in_time.mp4` | Global paleobathymetry, present day back to 170 Ma |
| `paleobathymetry_forward_in_time.mp4` | The same, 170 Ma forward to the present |

171 frames each at 8 frames per second. The carbonate-thickness animation uses the
central scenario and the same non-linear colour scale as Figure 3, so a colour in
the movie means the same thickness as that colour in the figure.

---

## Using the archive with the code

```bash
git clone https://github.com/EarthByte/CCD.git
cd CCD
pip install -r requirements.txt

python run_ccd_core.py            # steps 1-7; needs nothing from this archive
```

To repeat the grid-producing stages, unpack the archives and give
`pipeline_carbon/config.sh` their locations:

```bash
mkdir -p ccd_grids && cd ccd_grids
find .. -maxdepth 1 -name '*.tar.gz' -exec tar xzf {} \;
```

The thickness scenarios belong under
`steps/step8_carbonate_sediment_thickness/`, the paleobathymetry under its
`input_grids/`, the plate model under
`steps/step9_carbonate_volume_analysis/input/`, and the continental masks under
`steps/step10_carbon_cycle_degassing/.../Grids/InputGrids/`. `config.sh` names
every one of these paths and is the single place to change them.

---

## Not in this archive

The carbon-cycle stage (step 10) also reads seafloor-age, spreading-rate, total
sediment thickness, crustal-carbon and reservoir grids. Those are third-party
products with their own citations and download locations, given in the repository
at `steps/step10_carbon_cycle_degassing/`, and they are not redistributed here.

---

## Licence and citation

The grids and animations are released under CC BY 4.0. If you use them, cite the
paper, this archive, and the plate model:

> Dutkiewicz, A., and Müller, R.D. *Evolution and drivers of the global carbonate
> compensation depth over the past 170 million years.* Geology.

> Dutkiewicz, A., and Müller, R.D. *Evolution and drivers of the global carbonate
> compensation depth over the past 170 million years — model grids and animations.*
> Zenodo, https://doi.org/10.5281/zenodo.22785465

The carbonate sediment thickness model follows the method of Dutkiewicz, A.,
Müller, R.D., Cannon, J., Vaughan, S. and Zahirovic, S., 2019, Sequestration and
subduction of deep-sea carbonate in the global ocean since the Early Cretaceous:
Geology, v. 47, p. 91–94, https://doi.org/10.1130/G45424.1.
