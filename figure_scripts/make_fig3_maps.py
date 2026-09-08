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
CPT   = CW/"data/carbonate_thickness_blue_orange_red_pale.cpt"
CMASK = CW/"steps/step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26/Grids/InputGrids/ContinentalMasks"
OUT   = FIGDIR
for _p in (MODEL, CPT):
    if not Path(_p).exists(): raise SystemExit(f"missing input: {_p}")
if not (STEP8.is_dir() or STAGE.is_dir()):
    raise SystemExit(f"missing input: no carbonate-thickness grids at {STEP8} or {STAGE}")
OUT.mkdir(exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

# Three maps stacked vertically, with the carbonate budget below them at the same
# width. A Mollweide map is twice as wide as it is tall, so at the journal's full
# 190 mm column three of them plus the budget would run past the page; the figure is
# sized to fit a page instead, 117 mm wide by about 225 mm tall.
OUT_DPI = 300      # the fit measures at this dpi too, so the two cannot drift
TIMES=[0,25,115]; LABELS=["a","b","c"]
proj=ccrs.Mollweide(central_longitude=0)
# The maps get the full column. Panel (d) starts from the same cell and is then fitted
# to them further down, so that its axis labels finish flush with the map edges rather
# than overhanging. Row 4 is an empty spacer that the shared colour bar sits in, so it
# cannot land on top of the last map.
fig=plt.figure(figsize=(5.4,8.85))
gs=fig.add_gridspec(5,1,height_ratios=[1.0,1.0,1.0,0.34,0.72],
                    left=0.165,right=0.855,top=0.995,bottom=0.070,hspace=0.05)
axes=[fig.add_subplot(gs[i,0],projection=proj) for i in range(3)]
axd=fig.add_subplot(gs[4,0])
_spacer=fig.add_subplot(gs[3,0]); _spacer.axis("off")

for ax,T,lab in zip(axes,TIMES,LABELS):
    t0=time.time()
    ax.set_global(); ax.spines["geo"].set_linewidth(0.7)
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
    gp.plot_coastlines(ax, color="0.4", linewidth=0.22, zorder=3)
    for w,c,zz in [(1.5,"white",4),(0.65,"black",5)]:
        gp.plot_topological_plate_boundaries(ax, color=c, linewidth=w, zorder=zz)
    gp.plot_subduction_teeth(ax, color="black", zorder=6)
    # graticule: parallels every 30 deg and meridians every 60 deg. The +-180
    # meridians are omitted because the Mollweide outline already draws them.
    ax.gridlines(ylocs=[-60, -30, 0, 30, 60], xlocs=[-120, -60, 0, 60, 120],
                 color="0.30", linewidth=0.35, alpha=0.55, zorder=8)
    for _lat in (-60, -30, 0, 30, 60):
        ax.plot([-180, -174], [_lat, _lat], transform=ccrs.PlateCarree(),
                color="black", lw=0.9, solid_capstyle="butt", zorder=9)
        _t = "0\u00b0" if _lat == 0 else f"{abs(_lat)}\u00b0{'N' if _lat > 0 else 'S'}"
        ax.text(-171, _lat, _t, transform=ccrs.PlateCarree(), fontsize=6.0,
                va="center", ha="left", zorder=10,
                path_effects=[_pe.withStroke(linewidth=1.8, foreground="white")])
    # Panel letters are placed later, in figure coordinates, so that all four line up.
    # age label shifted ~4 mm left of its former position so it clears the Mollweide outline
    _panel_mm = ax.get_position().width*fig.get_size_inches()[0]*25.4
    ax.text(0.145-4.0/_panel_mm, 0.972, f"{T} Ma", transform=ax.transAxes, fontsize=10,
            va="top", ha="left", zorder=10)
    print(f"{T} Ma done {round(time.time()-t0,1)}s")

# shared colorbar for the three maps, inside the spacer row
_sp=_spacer.get_position()
# High in the spacer row: its lower part holds panel (d) event labels.
cax=fig.add_axes([_sp.x0+0.17*_sp.width, _sp.y0+0.80*_sp.height, 0.66*_sp.width, 0.009])
sm=ScalarMappable(cmap=cmap,norm=norm)
cb=fig.colorbar(sm, cax=cax, orientation="horizontal", extend="max", spacing="proportional")
cb.set_label("Compacted carbonate sediment thickness (m)", fontsize=7.5, labelpad=2)
# The 0-50 m bands are narrow on a proportional bar, so labelling every one of them
# ran the numbers together; label every second.
cb.set_ticks([0,20,40,100,160,210,270,320]); cb.ax.tick_params(labelsize=6.5, pad=1.5)

# (d) the carbonate budget, drawn by carbonate_budget.py so the panel and the data
# file cannot drift apart
import importlib.util as _ilu
_spec=_ilu.spec_from_file_location("carbonate_budget", _HERE/"carbonate_budget.py")
_cb=_ilu.module_from_spec(_spec); _spec.loader.exec_module(_cb)
_series=_cb.budget_series()
_bx=_cb.draw_budget(axd, _series, label_size=7.5, tick_size=6.5)
_cb.write_csv(_series)

# The maps are ellipses that reach the edge of their axes box; panel (d) is a
# rectangle carrying a two-line y-label on each side. Given the same box, those
# labels overhang the maps by about 13 mm a side and the panel reads as wider than
# them. Shrink its box until its own drawn ink - tick labels and axis labels
# included - spans exactly what the maps span. Font sizes are in points and are not
# touched by this, so nothing gets smaller to read.
def _ink_span(_axs, _pad=2):
    """Left and right edge of these axes' drawn ink, as fractions of figure width.

    Measured off the rendered canvas, not from get_tightbbox: that pads around glyph
    boxes, which left the panel a few per cent wider than the maps. The canvas is
    rendered at the dpi the figure is saved at, because glyph advances are hinted to
    whole pixels - measured at the default 100 dpi, the two-line y-label came out
    narrower than it renders at 300 and the panel was fitted too far left.
    """
    _old_dpi = fig.get_dpi()
    fig.set_dpi(OUT_DPI)
    try:
        fig.canvas.draw()
        _buf = np.asarray(fig.canvas.buffer_rgba())[..., :3].mean(axis=2)
        _H, _W = _buf.shape
        _r = fig.canvas.get_renderer()
        _bb = [a.get_tightbbox(_r) for a in _axs]
        _y0 = max(0, int(_H - max(b.y1 for b in _bb)) - _pad)
        _y1 = min(_H, int(_H - min(b.y0 for b in _bb)) + _pad)
        _cols = np.where((_buf[_y0:_y1, :] < 245).any(axis=0))[0]
        return _cols.min()/_W, (_cols.max()+1)/_W
    finally:
        fig.set_dpi(_old_dpi)

_t0, _t1 = _ink_span(axes[:1])
for _ in range(5):
    _d0, _d1 = _ink_span([axd, _bx])
    if abs(_d0-_t0) < 3e-4 and abs(_d1-_t1) < 3e-4:
        break
    _p = axd.get_position()
    _nx0 = _p.x0 + (_t0-_d0); _nx1 = _p.x1 - (_d1-_t1)
    for _a in (axd, _bx):
        _a.set_position([_nx0, _p.y0, _nx1-_nx0, _p.height])
_mm = fig.get_size_inches()[0]*25.4
_d0, _d1 = _ink_span([axd, _bx])
print(f"  panel (d) fitted to the maps: {axd.get_position().width*_mm:.0f} mm frame, "
      f"ink {(_d1-_d0)*_mm:.1f} mm against the maps' {(_t1-_t0)*_mm:.1f} mm")

# Panel letters in FIGURE coordinates at one x. Axes coordinates will not do it: the
# Mollweide panels are aspect-constrained, so their drawn box is narrower than the
# gridspec cell and the same transAxes offset lands in a different place on each.
# Left of panel (d) two-line y-label, which starts at x = 0.058 of the figure width;
# the gridspec left margin was widened to 0.165 to make that room.
_LETTER_X = 0.018
for _a, _lab in zip(axes + [axd], "abcd"):
    fig.text(_LETTER_X, _a.get_position().y1, _lab, fontsize=13, fontweight="bold",
             va="top", ha="left", zorder=10)
print(f"  budget panel: r = {_series['r']:.2f} between net gain and area above the CCD")
_bad=_cb.text_collisions(fig,[(axd,True),(_bx,False),(cax,False)])
print("  text collisions (panel d): "+(", ".join(_bad) if _bad else "none"))

for ext in ("png","pdf"): fig.savefig(OUT/f"Fig3_carbonate_thickness_maps.{ext}", dpi=OUT_DPI)
print("saved Fig3 ->", OUT/"Fig3_carbonate_thickness_maps.png")
