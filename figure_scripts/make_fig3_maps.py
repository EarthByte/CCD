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
STAGE = STEP9/"_cloud_stage"                       # carb_thick_*Ma.nc
MODEL = STEP9/"input/Alfonso_etal_2024_modClennettMuller"
CPT   = CW/"data/carbonate_thickness_blue_orange_red_pale.cpt"
CMASK = CW/"steps/step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26/Grids/InputGrids/ContinentalMasks"
OUT   = FIGDIR
for _p in (STAGE, MODEL, CPT):
    if not Path(_p).exists(): raise SystemExit(f"missing input: {_p}")
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

# The carb_thick_*Ma.nc grids in _cloud_stage are stored in reconstructed (paleo)
# coordinates, so they must NOT be rotated again. Set False only if you swap in
# grids that are still in present-day coordinates.
GRIDS_ARE_RECONSTRUCTED = True

TIMES=[0,35,50,80]; LABELS=["a","b","c","d"]
proj=ccrs.Mollweide(central_longitude=0)
fig,axes=plt.subplots(2,2, figsize=(10.4,6.1), subplot_kw={"projection":proj})
fig.subplots_adjust(left=0.01,right=0.99,top=0.985,bottom=0.11,wspace=0.04,hspace=0.02)

for ax,T,lab in zip(axes.ravel(),TIMES,LABELS):
    t0=time.time()
    ax.set_global(); ax.spines["geo"].set_linewidth(0.7)
    ax.set_facecolor("0.74")   # any cell with no grid data reads as continental crust, not white
    d=xr.open_dataset(STAGE/f"carb_thick_{T}Ma.nc"); z=d["z"].values
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
        ax.text(-171, _lat, _t, transform=ccrs.PlateCarree(), fontsize=7.5,
                va="center", ha="left", zorder=10,
                path_effects=[_pe.withStroke(linewidth=1.8, foreground="white")])
    ax.text(0.015,0.985, f"{lab}", transform=ax.transAxes, fontsize=13, fontweight="bold",
            va="top", ha="left", zorder=10)
    # age label shifted ~4 mm left of its former position so it clears the Mollweide outline
    _panel_mm = ax.get_position().width*fig.get_size_inches()[0]*25.4
    ax.text(0.145-4.0/_panel_mm, 0.972, f"{T} Ma", transform=ax.transAxes, fontsize=12,
            va="top", ha="left", zorder=10)
    print(f"{T} Ma done {round(time.time()-t0,1)}s")

# shared colorbar
cax=fig.add_axes([0.24,0.065,0.52,0.028])
sm=ScalarMappable(cmap=cmap,norm=norm)
cb=fig.colorbar(sm, cax=cax, orientation="horizontal", extend="max", spacing="proportional")
cb.set_label("Compacted carbonate sediment thickness (m)", fontsize=10)
cb.set_ticks([0,10,20,30,40,50,100,160,210,270,320]); cb.ax.tick_params(labelsize=8)
for ext in ("png","pdf"): fig.savefig(OUT/f"Fig3_carbonate_thickness_maps.{ext}", dpi=300, bbox_inches="tight")
print("saved Fig3 ->", OUT/"Fig3_carbonate_thickness_maps.png")
