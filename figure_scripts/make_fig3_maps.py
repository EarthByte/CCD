import glob, time, warnings, numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")

# ---- paths (resolved from this script's location; runnable from any directory) ----
# ---- path resolution -------------------------------------------------------
# Works both in the working tree (this script in Paper/, the workflow in a
# sibling CCD_workflow_clean/) and in the public repo (this script in
# figure_scripts/, the workflow at the repo root). CW = workflow root,
# FIGDIR = where figures and the small figure-input curves live.
_HERE = Path(__file__).resolve().parent
def _workflow_root():
    for _p in (_HERE.parent / "CCD_workflow_clean", _HERE.parent, _HERE):
        if (_p / "steps").is_dir() and (_p / "ccdworkflow").is_dir():
            return _p
    raise SystemExit("cannot locate the CCD workflow root (expected steps/ + ccdworkflow/)")
CW = _workflow_root()
FIGDIR = _HERE / "Figures" if (_HERE / "Figures").is_dir() else CW / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)
STEP9 = CW/"steps/step9_carbonate_volume_analysis"
# Carbonate-thickness grids. The live source is the carbonate-thickness step, so the
# figure always reflects the current run; the hand-staged _cloud_stage snapshot is
# kept only as a fallback for a checkout without the grids. Reading the snapshot by
# default silently froze this figure at the grids of whichever run staged it.
STEP8 = CW/"steps/step8_carbonate_sediment_thickness/carbonate_sed_thickness_DM2026"
STAGE = STEP9/"_cloud_stage"                       # fallback: carb_thick_*Ma.nc


def _thickness_grid(T):
    live = STEP8/f"compacted_sediment_thickness_0.25_{T}.nc"
    if live.exists():
        return live
    snap = STAGE/f"carb_thick_{T}Ma.nc"
    if snap.exists():
        print(f"  [fig3] {T} Ma: using the _cloud_stage snapshot, not the current step-8 grid")
        return snap
    raise SystemExit(f"no carbonate-thickness grid for {T} Ma (looked in {STEP8} and {STAGE})")
MODEL = STEP9/"input/Alfonso_etal_2024_modClennettMuller"
CPT   = CW/"data/carbonate_thickness_roma_pale.cpt"
CMASK = CW/"steps/step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26/Grids/InputGrids/ContinentalMasks"
OUT   = FIGDIR
for _p in (MODEL, CPT):
    if not Path(_p).exists(): raise SystemExit(f"missing input: {_p}")
if not (STEP8.is_dir() or STAGE.is_dir()):
    raise SystemExit(f"missing input: no carbonate-thickness grids at {STEP8} or {STAGE}")
OUT.mkdir(exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
# Geology asks for Helvetica or Arial and for lettering of 7-12 pt, with part labels of
# 13-16 pt, at the size the figure is printed. This figure is built at the journal's
# full page width, so the sizes below are the sizes the reader sees.
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "Nimbus Sans",
                                          "Liberation Sans", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42      # embed as TrueType, not as outlines
matplotlib.rcParams["axes.unicode_minus"] = False
# Superscripts and subscripts - the 10^6 km^3 Myr^-1 of panel (D)'s axis labels - are set
# by matplotlib's maths typesetter, which ignores the family above and has its own font.
# Left alone it renders them in DejaVu Sans and embeds a second face in the PDF, so a
# label reads in two typefaces and Illustrator inherits both. Pointing every maths face
# at the sans-serif chain puts the whole figure in one font.
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.default"] = "regular"
for _k, _v in (("mathtext.rm", "sans"), ("mathtext.it", "sans:italic"),
               ("mathtext.bf", "sans:bold"), ("mathtext.sf", "sans"),
               ("mathtext.tt", "sans"), ("mathtext.cal", "sans:italic")):
    matplotlib.rcParams[_k] = _v

# Which face matplotlib actually resolved. font.sans-serif is a fallback chain, so a
# Helvetica it cannot find is not an error but a silent substitution, and the first the
# figure would know of it is the journal asking why the PDF carries DejaVu Sans.
from matplotlib import font_manager as _fm
try:
    _font_file = _fm.findfont(_fm.FontProperties(family=matplotlib.rcParams["font.sans-serif"]))
    _font_name = _fm.get_font(_font_file).family_name
except Exception as _e:                       # never let a font query stop the figure
    _font_file, _font_name = "", f"unresolved ({_e})"
print(f"  text font: {_font_name}" + (f" [{Path(_font_file).name}]" if _font_file else ""))
if _font_name.split(" ")[0] not in ("Helvetica", "Arial"):
    print("  WARNING: Geology asks for Helvetica or Arial, and matplotlib resolved "
          f"{_font_name}. Put Helvetica or Arial where matplotlib can see it, delete its "
          "font cache (~/.matplotlib/fontlist-*.json), and run this again.")
PT_LABEL, PT_TICK, PT_EVENT, PT_AGE, PT_LETTER = 8.0, 7.0, 7.0, 8.0, 13.0
W_PAGE_MM = 185.0                             # Geology, full page

# ---- line weights ----------------------------------------------------------
# Every rule on the page is drawn through thin(), so the figure carries one line-weight
# scale instead of a dozen independent numbers. The PDF hands Illustrator whatever
# weights it was drawn with, and at 185 mm the earlier ones printed heavier than the
# type they sit beside. Scaling them together leaves the relative emphasis between a
# plate boundary, a coastline and a graticule exactly as it was.
LW_SCALE = 0.75
LW_FLOOR = 0.20                               # pt, the weight below which a rule breaks up in print


def thin(w: float) -> float:
    """A drawn line weight in points, on the figure's common scale."""
    return max(w * LW_SCALE, LW_FLOOR)


# The rules matplotlib draws on its own - the spines of panel (D), the box around the
# colour bar and every tick - take the 0.8 pt default unless they are told otherwise,
# which is why the colour-bar box printed heavier than the maps beside it.
matplotlib.rcParams["axes.linewidth"] = thin(0.8)
matplotlib.rcParams["patch.linewidth"] = thin(1.0)
matplotlib.rcParams["grid.linewidth"] = thin(0.8)
for _k in ("xtick.major.width", "ytick.major.width"):
    matplotlib.rcParams[_k] = thin(0.8)
for _k in ("xtick.minor.width", "ytick.minor.width"):
    matplotlib.rcParams[_k] = thin(0.6)
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
import cartopy.crs as ccrs
import gplately, pygplates, xarray as xr
from matplotlib.colors import ListedColormap as _LCM
import matplotlib.patheffects as _pe

# ---- parse GMT cpt ----
bounds=[]; colors=[]; over=under=None
for line in open(CPT):
    s=line.split()
    if not s or s[0] in ("#",): continue
    if s[0]=="B": under=tuple(int(x)/255 for x in s[1].split("/")); continue
    if s[0]=="F": over=tuple(int(x)/255 for x in s[1].split("/")); continue
    if s[0]=="N": continue
    try: z0=float(s[0]); z1=float(s[2])
    except: continue
    c=tuple(int(x)/255 for x in s[1].split("/"))
    if not bounds: bounds.append(z0)
    bounds.append(z1); colors.append(c)
cmap=ListedColormap(colors); cmap.set_over(over); cmap.set_under(under); cmap.set_bad("none")   # NaN transparent, so the grey continental background shows through
norm=BoundaryNorm(bounds, cmap.N)
print("cpt bounds", bounds)

# ---- model ----
rot = sorted(glob.glob(str(MODEL/"Rotations/*.rot"))) + sorted(glob.glob(str(MODEL/"DeformingMeshes/*.rot")))
topo = sorted(glob.glob(str(MODEL/"DeformingMeshes/*.gpml"))) + sorted(glob.glob(str(MODEL/"PlateBoundaries/*.gpml")))
static = str(MODEL/"StaticPolygons/Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons.gpml")
coast = str(MODEL/"Coastlines/Clennett__etal_2020_Coastlines.gpml")
cont  = str(MODEL/"ContinentalPolygons/Global_PresentDay_ContPolygons_2019_v2.gpml")
rotm = pygplates.RotationModel(rot)
model = gplately.PlateReconstruction(rotm, topology_features=topo, static_polygons=static, anchor_plate_id=0)

# The carbonate-thickness grids are stored in reconstructed (paleo)
# coordinates, so they must NOT be rotated again. Set False only if you swap in
# grids that are still in present-day coordinates.
GRIDS_ARE_RECONSTRUCTED = True

# Two by two: three Mollweide maps and the carbonate budget. A Mollweide map is twice
# as wide as it is tall, so one column of three plus the budget ran 225 mm down the
# page; side by side they fit the page width with room for the budget in the fourth
# cell. Axes are placed in figure fractions computed from millimetres, so the panel
# sizes are the printed ones rather than whatever a gridspec leaves over.
OUT_DPI = 600      # the raster the journal wants for the assembled artwork
TIMES=[0,34,115]   # 34 Ma is the Eocene-Oligocene transition
proj=ccrs.Mollweide(central_longitude=0)

MAP_W, COL_GAP, LEFT_PAD = 89.0, 5.5, 1.5     # mm
MAP_H = MAP_W/2.0
ROW_GAP = 5.0
CB_BLOCK = 13.0                               # mm under the 115 Ma map for the scale
# The budget panel carries a two-line axis label and its tick labels on BOTH sides,
# about 13.5 mm each, so its frame is that much narrower than the map beside it.
BUD_W, BUD_H = 58.0, 35.0                     # mm, panel (D) frame
BUD_LEFT, BUD_BOTTOM = 15.5, 9.5              # mm of its cell taken by labels
FIG_W = LEFT_PAD + 2*MAP_W + COL_GAP
FIG_H = CB_BLOCK + MAP_H + ROW_GAP + MAP_H + 6.0   # the last term is the part letters' row
fig=plt.figure(figsize=(FIG_W/25.4, FIG_H/25.4))

def _cell(x_mm, y_mm, w_mm, h_mm):
    """Axes rectangle in figure fractions, from millimetres on the page."""
    return [x_mm/FIG_W, y_mm/FIG_H, w_mm/FIG_W, h_mm/FIG_H]

_row_top = CB_BLOCK + MAP_H + ROW_GAP
_col2 = LEFT_PAD + MAP_W + COL_GAP
_cells = {"A": (LEFT_PAD, _row_top), "B": (_col2, _row_top),
          "C": (LEFT_PAD, CB_BLOCK), "D": (_col2, CB_BLOCK)}
axes=[fig.add_axes(_cell(*_cells[l], MAP_W, MAP_H), projection=proj)
      for l in "ABC"]
axd=fig.add_axes(_cell(_cells["D"][0]+BUD_LEFT, _cells["D"][1]+BUD_BOTTOM, BUD_W, BUD_H))

for ax,T in zip(axes,TIMES):
    t0=time.time()
    ax.set_global(); ax.spines["geo"].set_linewidth(thin(0.7))
    ax.set_facecolor("0.74")   # any cell with no grid data reads as continental crust, not white
    d=xr.open_dataset(_thickness_grid(T)); z=d["z"].values
    if GRIDS_ARE_RECONSTRUCTED or T==0:
        data=z                      # grid is already in paleo-coordinates: plot as is
    else:                           # only for grids supplied in present-day coordinates
        r=gplately.Raster(data=z, plate_reconstruction=model, extent=[-180,180,-90,90], origin="lower")
        data=r.reconstruct(time=T, partitioning_features=static).data
    ax.imshow(data, origin="lower", extent=[-180,180,-90,90], transform=ccrs.PlateCarree(),
              cmap=cmap, norm=norm, zorder=1)
    # Alfonso2024 continental mask, same grey as the continents. Without this the
    # deforming-network areas (not covered by the rigid continental polygons) render
    # white and punch holes in the map.
    _cm = CMASK/f"Alfonso2024_continent_mask_{float(T):.2f}Ma.nc"
    if Path(_cm).exists():
        _m = xr.open_dataset(_cm)["z"].values
        ax.imshow(np.where(_m>0, 1.0, np.nan), origin="lower", extent=[-180,180,-90,90],
                  transform=ccrs.PlateCarree(), cmap=_LCM(["0.74"]), vmin=0, vmax=1, zorder=1.6)
    else:
        print(f"  warning: no continental mask for {T} Ma ({_cm.name}) - deforming areas will be white")
    gp=gplately.PlotTopologies(model, coastlines=coast, continents=cont, time=T)
    gp.plot_continents(ax, facecolor="0.74", edgecolor="none", zorder=2)
    gp.plot_coastlines(ax, color="0.4", linewidth=thin(0.22), zorder=3)
    for w,c,zz in [(thin(1.5),"white",4),(thin(0.65),"black",5)]:
        gp.plot_topological_plate_boundaries(ax, color=c, linewidth=w, zorder=zz)
    gp.plot_subduction_teeth(ax, color="black", zorder=6)
    # graticule: parallels every 30 deg and meridians every 60 deg. The +-180
    # meridians are omitted because the Mollweide outline already draws them.
    ax.gridlines(ylocs=[-60, -30, 0, 30, 60], xlocs=[-120, -60, 0, 60, 120],
                 color="0.30", linewidth=thin(0.35), alpha=0.55, zorder=8)
    for _lat in (-60, -30, 0, 30, 60):
        ax.plot([-180, -174], [_lat, _lat], transform=ccrs.PlateCarree(),
                color="black", lw=thin(0.9), solid_capstyle="butt", zorder=9)
        _t = "0\u00b0" if _lat == 0 else f"{abs(_lat)}\u00b0{'N' if _lat > 0 else 'S'}"
        ax.text(-171, _lat, _t, transform=ccrs.PlateCarree(), fontsize=6.0,
                va="center", ha="left", zorder=10,
                path_effects=[_pe.withStroke(linewidth=thin(1.8), foreground="white")])
    # Panel letters are placed later, in figure coordinates, so that all four line up.
    # age label in the top left of the map box, where the Mollweide outline leaves the
    # corner empty
    ax.text(0.085, 0.975, f"{T} Ma", transform=ax.transAxes, fontsize=PT_AGE,
            va="top", ha="left", zorder=10,
            path_effects=[_pe.withStroke(linewidth=thin(2.0), foreground="white")])
    print(f"{T} Ma done {round(time.time()-t0,1)}s")

# Shared colour bar, in the band under the 115 Ma map. Uniform spacing draws every
# class at the same width: the ten classes below 50 m are the ones a reader has to tell
# apart, and on a bar proportional to thickness they share a sixth of its length. The
# boundaries are all labelled, so the scale stays readable as a non-linear one.
cax=fig.add_axes(_cell(LEFT_PAD+9.0, CB_BLOCK-5.0, MAP_W-18.0, 3.0))
sm=ScalarMappable(cmap=cmap,norm=norm)
cb=fig.colorbar(sm, cax=cax, orientation="horizontal", extend="max", spacing="uniform")
cb.set_label("Compacted carbonate sediment thickness (m)", fontsize=PT_LABEL, labelpad=2)
cb.set_ticks(bounds[:-1]); cb.ax.tick_params(labelsize=PT_TICK, pad=1.5)

# (d) the carbonate budget, drawn by carbonate_budget.py so the panel and the data
# file cannot drift apart
import importlib.util as _ilu
_spec=_ilu.spec_from_file_location("carbonate_budget", _HERE/"carbonate_budget.py")
_cb=_ilu.module_from_spec(_spec); _spec.loader.exec_module(_cb)
_series=_cb.budget_series()
_bx=_cb.draw_budget(axd, _series, label_size=PT_LABEL, tick_size=PT_TICK,
                    event_size=PT_EVENT, lw_scale=LW_SCALE)
_cb.write_csv(_series)

# Part letters at the top left of each cell, in figure fractions: the maps and the
# budget panel have different drawn widths, so axes coordinates would not line them up.
for _l in "ABCD":
    _x, _y = _cells[_l]
    fig.text(_x/FIG_W, (_y+MAP_H)/FIG_H + 0.004, _l, fontsize=PT_LETTER,
             fontweight="bold", va="bottom", ha="left", zorder=10)
print(f"  budget panel: r = {_series['r']:.2f} between net gain and area above the CCD")
_bad=_cb.text_collisions(fig,[(axd,True),(_bx,False),(cax,False)])
print("  text collisions (panel d): "+(", ".join(_bad) if _bad else "none"))

# The weights the page actually carries, reported rather than eyeballed: a change of
# LW_SCALE cannot quietly take a rule under what a press can hold, and any weight the
# floor had to catch is named.
_WEIGHTS = {"plate boundaries": 0.65, "plate-boundary halo": 1.5, "coastlines": 0.22,
            "map outline": 0.7, "graticule": 0.35, "latitude ticks": 0.9,
            "colour-bar box, panel (D) frame and ticks": 0.8,
            "budget curves": 1.6}
print("  line weights at printed size (pt): "
      + ", ".join(f"{k} {thin(w):.2f}" for k, w in _WEIGHTS.items()))
_held = [k for k, w in _WEIGHTS.items() if w * LW_SCALE < LW_FLOOR]
if _held:
    print(f"  held at the {LW_FLOOR:.2f} pt floor rather than scaled: " + ", ".join(_held))

for ext in ("png","pdf"):
    # No tight bounding box: the canvas IS the printed page, 185 mm wide, so the type
    # sizes above survive placement unscaled. A tight box trims to the ink and the
    # journal then scales whatever it gets.
    fig.savefig(OUT/f"Fig3_carbonate_thickness_maps.{ext}", dpi=OUT_DPI)

# The PDF is what Illustrator opens, so the font is checked on the artwork rather than
# on the settings that were meant to produce it: every face the file embeds is named,
# and a second family means a run of text got past the sans-serif chain.
import re as _re
_pdf = (OUT/"Fig3_carbonate_thickness_maps.pdf").read_bytes()
_faces = sorted({m.decode().split("+")[-1]
                 for m in _re.findall(rb"/BaseFont\s*/([A-Za-z0-9+#.\-]+)", _pdf)})
print("  fonts embedded in the PDF: " + (", ".join(_faces) if _faces else "none found"))
_fams = {f.split("-")[0] for f in _faces}
if len(_fams) > 1:
    print("  WARNING: the PDF carries more than one typeface: " + ", ".join(sorted(_fams)))
_right_edge = _cells["D"][0] + BUD_LEFT + BUD_W
if _right_edge + 15.0 > FIG_W + 0.5:
    print(f"  WARNING: panel (D) leaves {FIG_W-_right_edge:.1f} mm for its right-hand "
          f"label and tick labels, which need about 15 mm")
print(f"  page {FIG_W:.0f} x {FIG_H:.0f} mm: maps {MAP_W:.0f} x {MAP_H:.0f} mm, "
      f"budget {BUD_W:.0f} x {BUD_H:.0f} mm; type 7-13 pt at printed size")
print("saved Fig3 ->", OUT/"Fig3_carbonate_thickness_maps.png")
