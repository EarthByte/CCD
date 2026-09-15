#!/usr/bin/env python3
"""Figure 3 (pyGMT), full-page width for Geology (18.5 cm), as a 2 x 2 block.

(A) 0 Ma, (B) 34 Ma, (C) 115 Ma reconstructed compacted carbonate sediment
thickness on the Alfonso et al. (2025) plate model, and (D) the carbonate budget.

The colour scale is inverted roma stepped one colour per class, on the class
boundaries the thickness scale has always used, so the ten thin classes below
50 m keep their contrast instead of being squeezed into one end of a ramp.

The reconstructed vector layers - coastlines, continents, plate boundaries and
subduction polarity - are exported once to Figures/_fig3_layers/ and re-read on
later runs, so redrawing the figure does not re-run the plate reconstruction.
Pass --refresh to rebuild them.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pygmt

import paper_gmt as S

CW, OUT = S.CW, S.FIGDIR
TIMES = [0, 34, 115]                     # 34 Ma is the Eocene-Oligocene transition
CACHE = OUT / "_fig3_layers"
STEP8 = CW / "steps/step8_carbonate_sediment_thickness/carbonate_sed_thickness_DM2026"
STAGE = CW / "steps/step9_carbonate_volume_analysis/_cloud_stage"
MODEL = CW / "steps/step9_carbonate_volume_analysis/input/Alfonso_etal_2024_modClennettMuller"
CMASK = CW / "steps/step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26/Grids/InputGrids/ContinentalMasks"
BOUNDS = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 160, 210, 270, 320]
GREY = "189/189/189"                     # continental crust, and any cell with no data

# ---- geometry ---------------------------------------------------------------
W = S.W_PAGE
MAPW = 8.90                              # cm, Mollweide width; its height is half that
COLGAP = W - 2 * MAPW                    # 0.7 cm between the columns
MAPH = MAPW / 2.0
BUDW, BUDH = 7.00, 3.85                  # cm, panel (D) frame
ROWGAP = 0.80                            # cm between the rows: enough for the part letter
BOTTOM = 1.45                            # cm of clear page under the bottom row, which
                                         # carries the colour bar and panel (D)'s age axis
CBAR_H = 0.30


def thickness_grid(t):
    live = STEP8 / f"compacted_sediment_thickness_0.25_{t}.nc"
    if live.exists():
        return live
    snap = STAGE / f"carb_thick_{t}Ma.nc"
    if snap.exists():
        print(f"  [fig3] {t} Ma: using the _cloud_stage snapshot, not the current step-8 grid")
        return snap
    raise SystemExit(f"no carbonate-thickness grid for {t} Ma (looked in {STEP8} and {STAGE})")


# ---- reconstructed vector layers -------------------------------------------
LAYERS = ("coastlines", "continents", "boundaries", "subduction_left", "subduction_right")


def layer(t, name):
    return CACHE / f"{name}_{t}Ma.gmt"


def export_layers(refresh=False):
    """Reconstruct the vector layers once and write them as OGR_GMT files."""
    need = [t for t in TIMES if refresh or not all(layer(t, n).exists() for n in LAYERS)]
    if not need:
        return
    import gplately
    import pygplates

    rot = [str(p) for p in (MODEL / "Rotations").glob("*.rot")]
    topo = (sorted(str(p) for p in (MODEL / "DeformingMeshes").glob("*.gpml"))
            + sorted(str(p) for p in (MODEL / "PlateBoundaries").glob("*.gpml")))
    static = str(MODEL / "StaticPolygons/Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons.gpml")
    coast = str(MODEL / "Coastlines/Clennett__etal_2020_Coastlines.gpml")
    cont = str(MODEL / "ContinentalPolygons/Global_PresentDay_ContPolygons_2019_v2.gpml")
    model = gplately.PlateReconstruction(pygplates.RotationModel(rot), topology_features=topo,
                                         static_polygons=static, anchor_plate_id=0)
    CACHE.mkdir(parents=True, exist_ok=True)
    for t in need:
        gp = gplately.PlotTopologies(model, coastlines=coast, continents=cont, time=t)
        left, right = gp.get_subduction_direction()
        for name, gdf in (("coastlines", gp.get_coastlines()),
                          ("continents", gp.get_continents()),
                          ("boundaries", gp.get_all_topological_sections()),
                          ("subduction_left", left),
                          ("subduction_right", right)):
            path = layer(t, name)
            if gdf is None or len(gdf) == 0:
                path.write_text("# empty\n")
                continue
            gdf.to_file(path, driver="OGR_GMT")
        print(f"  [fig3] exported layers for {t} Ma")


# Only the boundaries that are plate boundaries in the ordinary sense. The exported
# file also holds the edges of the deforming networks and the various crust-type and
# terrane lines the model carries (unclassified, extended continental crust, slab edges,
# inferred palaeo-boundaries), which crowd the maps without saying anything about
# carbonate. Subduction zones stay in so a trench still draws where the polarity layers
# have no segment for it; where they do, the teeth are drawn on top of the same line.
BOUNDARY_TYPES = ("gpml:MidOceanRidge", "gpml:Transform", "gpml:SubductionZone")


def filtered_boundaries(t):
    """Copy of the boundary file holding only BOUNDARY_TYPES, written beside it.

    The OGR_GMT export carries feature_type in each segment header, so the selection is
    made here rather than by re-running the reconstruction.
    """
    src = layer(t, "boundaries")
    if not src.exists():
        return None
    out = CACHE / f"boundaries_plate_{t}Ma.gmt"
    lines = src.read_text(errors="replace").splitlines()
    kept, keep, head = [], False, []
    for line in lines:
        if line.startswith("#") and not line.startswith("# @D"):
            if not kept and not head:
                head.append(line)
            elif line == "# FEATURE_DATA":
                head.append(line)
            continue
        if line.startswith(">"):
            keep = False
            pending = [line]
            continue
        if line.startswith("# @D"):
            keep = any(f"|{ft}|" in line for ft in BOUNDARY_TYPES)
            if keep:
                kept.extend(pending + [line])
            continue
        if keep:
            kept.append(line)
    out.write_text("\n".join(["# @VGMT1.0", "# @GLINESTRING", "# FEATURE_DATA"] + kept) + "\n")
    return out


# ---- colour -----------------------------------------------------------------
def continent_mask(t):
    """The mask grid holds 0 over ocean, and GMT paints a below-range value with the
    CPT's background colour - which covered the thickness grid with white. Clip the
    zeros to NaN so only the continental cells are painted."""
    src = CMASK / f"Alfonso2024_continent_mask_{float(t):.2f}Ma.nc"
    if not src.exists():
        return None
    out = CACHE / f"continent_mask_{t}Ma.nc"
    if not out.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        pygmt.grdclip(grid=str(src), outgrid=str(out), below=[0.5, "NaN"])
    return out


def build_cpts():
    # Washed out to the same degree as the scale this replaces, so the plate boundaries
    # and coastlines drawn over it still read: chroma cut to 0.45 and lightness lifted.
    cpt = S.stepped_cpt("roma", BOUNDS, OUT / "carbonate_thickness_roma_stepped.cpt",
                        reverse=True, nan=GREY, pale=0.45,
                        header="Compacted carbonate sediment thickness (m), pale")
    grey = OUT / "_continent_grey.cpt"
    # The NaN colour has to differ from the colour the mask itself is drawn in. GMT
    # makes NaN cells transparent by colour masking, and when the two are the same it
    # cannot tell them apart - which painted every ocean cell of the mask solid red on
    # top of the thickness grid. Magenta appears nowhere else in the figure.
    grey.write_text(f"# continental mask: anything above zero is continental crust\n"
                    f"0.5\t{GREY}\t1e6\t{GREY}\nB\t{GREY}\nF\t{GREY}\nN\t255/0/255\n")
    return cpt, grey


def main():
    export_layers(refresh="--refresh" in sys.argv)
    cpt, grey_cpt = build_cpts()
    budget = pd.read_csv(OUT / "Fig3_carbonate_budget.csv")

    fig = pygmt.Figure()
    pygmt.config(**S.defaults())
    S.begin(fig)
    checks = []

    # -------- the three maps and the budget, filled row by row from the bottom
    cells = {"A": (0, BOTTOM + MAPH + ROWGAP), "B": (MAPW + COLGAP, BOTTOM + MAPH + ROWGAP),
             "C": (0, BOTTOM), "D": (MAPW + COLGAP, BOTTOM)}

    for letter, t in zip("ABC", TIMES):
        x0, y0 = cells[letter]
        fig.shift_origin(xshift=f"{x0}c", yshift=f"{y0}c")
        proj = f"W{MAPW}c"
        fig.basemap(region="d", projection=proj, frame=["+g" + GREY])
        fig.grdimage(grid=str(thickness_grid(t)), cmap=str(cpt), projection=proj,
                     region="d", nan_transparent=True)
        mask = continent_mask(t)
        if mask is not None:
            # The deforming networks are not covered by the rigid continental polygons;
            # without the mask they render as holes in the continents.
            fig.grdimage(grid=str(mask), cmap=str(grey_cpt), projection=proj, region="d",
                         nan_transparent=True)
        else:
            print(f"  warning: no continental mask for {t} Ma")
        if layer(t, "continents").exists():
            fig.plot(data=str(layer(t, "continents")), fill=GREY, projection=proj, region="d")
        if layer(t, "coastlines").exists():
            fig.plot(data=str(layer(t, "coastlines")), pen="0.15p,gray40", projection=proj, region="d")
        bnd = filtered_boundaries(t)
        if bnd is not None:
            for pen in ("1.1p,white", "0.45p,black"):
                fig.plot(data=str(bnd), pen=pen, projection=proj, region="d")
        for side, flag in (("subduction_left", "+l"), ("subduction_right", "+r")):
            f = layer(t, side)
            if f.exists() and f.stat().st_size > 20:
                fig.plot(data=str(f), pen="0.45p,black", fill="black",
                         style=f"f0.45c/0.09c{flag}+t", projection=proj, region="d")
        fig.basemap(region="d", projection=proj, frame=["xg60", "yg30"])
        # Age label inside the top left of the map box, where the Mollweide outline
        # leaves the corner empty.
        fig.text(x=-168, y=76, text=f"{t} Ma", justify="LT", font=f"{S.PT_LABEL}p,{S.FONT},black",
                 fill="white@30", clearance="0.04c/0.02c", projection=proj, region="d", no_clip=True)
        fig.text(x=0.0, y=MAPH + 0.12, text=letter, justify="LB",
                 font=f"{S.PT_TAG}p,{S.FONT}-Bold,black", no_clip=True,
                 region=[0, MAPW, 0, MAPH + 0.9], projection=f"X{MAPW}c/{MAPH + 0.9}c")
        fig.shift_origin(xshift=f"-{x0}c", yshift=f"-{y0}c")

    # -------- shared colour bar, under the 115 Ma map ------------------------
    x0, y0 = cells["C"]
    CB_BLOCK = 1.40                       # cm of page under the 115 Ma map
    fig.shift_origin(xshift=f"{x0}c", yshift=f"{y0 - CB_BLOCK}c")
    # Equal-width classes: the ten classes below 50 m are what the reader needs to tell
    # apart, and on a bar proportional to thickness they take a sixth of its length.
    # Placed in paper coordinates below the map, with its own label drawn underneath -
    # GMT's -L (equal classes) refuses a -B that sets increments, so the label cannot
    # ride on the bar's own axis.
    CBW = MAPW - 1.9
    fig.colorbar(cmap=str(cpt), equalsize=0.0,
                 position=f"x{MAPW / 2}c/{CB_BLOCK - 0.30}c+w{CBW}c/{CBAR_H}c+h+jTC")
    fig.text(x=MAPW / 2, y=0.12, text="Compacted carbonate sediment thickness (m)",
             justify="CB", font=f"{S.PT_LABEL}p,{S.FONT},black", no_clip=True,
             region=[0, MAPW, 0, CB_BLOCK], projection=f"X{MAPW}c/{CB_BLOCK}c")
    fig.shift_origin(xshift=f"-{x0}c", yshift=f"-{y0 - CB_BLOCK}c")

    # -------- (D) the carbonate budget ---------------------------------------
    x0, y0 = cells["D"]
    # The budget panel is narrower than a map because it carries a two-line label on
    # each side; it is centred in the map's column so the two rows still read as a grid.
    dx = x0 + (MAPW - BUDW) / 2 + 0.35
    fig.shift_origin(xshift=f"{dx}c", yshift=f"{y0}c")
    VOL, AREA = S.rgb("#1f6f8b"), S.rgb("#b3202c")
    g = budget.dropna(subset=["net_gain_1e6km3_per_myr_smoothed"])
    gmax = float(np.nanmax(np.abs(g["net_gain_1e6km3_per_myr_smoothed"]))) * 1.15
    pD = S.Panel(region=(0, 170, -gmax, gmax), width=BUDW, height=BUDH, x_reversed=True)
    RD, PD = list(pD.region), pD.projection
    fig.basemap(region=RD, projection=PD, frame=["WSn", "xa20f10+lAge (Ma)", "ya0.5f0.25"])
    fig.plot(x=[0, 170], y=[0, 0], pen="0.4p,gray70", region=RD, projection=PD)
    import importlib.util as ilu
    _spec = ilu.spec_from_file_location("paper_events", S._HERE / "paper_events.py")
    _ev = ilu.module_from_spec(_spec)
    _spec.loader.exec_module(_ev)
    _ev.draw_events_gmt(fig, pD, pt=S.PT_ANNOT)
    fig.plot(x=g["Age_Ma"], y=g["net_gain_1e6km3_per_myr_smoothed"], pen=f"1.0p,{VOL}",
             region=RD, projection=PD)

    amax = float(np.nanmax(budget["area_above_CCD_1e6km2"])) * 1.10
    pA2 = S.Panel(region=(0, 170, 0, amax), width=BUDW, height=BUDH, x_reversed=True)
    RA2 = list(pA2.region)
    fig.basemap(region=RA2, projection=PD, frame=["E", "ya50f25"])
    fig.plot(x=budget["Age_Ma"], y=budget["area_above_CCD_1e6km2"], pen=f"0.9p,{AREA},4_2:0",
             region=RA2, projection=PD)
    fig.text(x=0.0, y=BUDH + 0.12, text="D", justify="LB",
             font=f"{S.PT_TAG}p,{S.FONT}-Bold,black", no_clip=True,
             region=[0, BUDW, 0, BUDH + 0.9], projection=f"X{BUDW}c/{BUDH + 0.9}c")
    # Two lines each, because the quantity and its units together run longer than the
    # panel is tall; coloured like the curve each one belongs to.
    lb = S.rotated_label(fig, pD, ["Net carbonate volume gain", "(10@+6@+ km@+3@+ Myr@+-1@+)"],
                         side="left", colour=VOL, gap=0.62)
    rb = S.rotated_label(fig, pD, ["Seafloor area above the CCD", "(10@+6@+ km@+2@+)"],
                         side="right", colour=AREA, gap=0.62)
    for name, text in (("left label", "Net carbonate volume gain"),
                       ("right label", "Seafloor area above the CCD")):
        if S.text_width_cm(text, S.PT_LABEL) > BUDH:
            checks.append(f"(D) {name} is longer than the {BUDH:.1f} cm panel height")
    fig.shift_origin(xshift=f"-{dx}c", yshift=f"-{y0}c")

    S.report(checks, "layout")
    small = S.too_small(S.PT_ANNOT, S.PT_LABEL, S.PT_LEG)
    print("  type below the 7 pt journal floor: " + (str(small) if small else "none"))
    print(f"  page {W:.2f} cm wide: maps {MAPW:.2f} x {MAPH:.2f} cm, budget {BUDW:.2f} x {BUDH:.2f} cm")

    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"Fig3_carbonate_thickness_maps_gmt.{ext}", dpi=600, crop=True)
    print("wrote Fig3_carbonate_thickness_maps_gmt")


if __name__ == "__main__":
    main()
