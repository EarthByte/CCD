#!/usr/bin/env python3
"""
Shared helpers for the two CCD-project supplementary videos
(compacted carbonate thickness + paleobathymetry), both on a Mollweide
projection with GPlately-reconstructed plate boundaries, subduction
teeth, gray continents and thin coastlines.

The netCDF grids are ALREADY reconstructed into Alfonso et al. (2025)
palaeo-coordinates at each Ma, so they are plotted as-is (no raster
reconstruction) with the vector overlays reconstructed to the same age.
Longitude convention differs between products, though: carbonate thickness
is -180..180, paleobathymetry is 0..360. _grid_south_up() normalises both
to -180..180 to match the plotting extent.

Requires: gplately, pygplates, cartopy, xarray, matplotlib, ffmpeg on PATH.
"""
from __future__ import annotations
import glob, os, subprocess, shutil, time as _t
import numpy as np, xarray as xr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
import cartopy.crs as ccrs
import gplately, pygplates

CONTINENT_GRAY = "0.74"          # fill for continental polygons AND grid NaNs
CBAR_GAP_MM    = 4.0             # clear space between the globe and the colour bar
_CBAR_H        = 0.028           # colour-bar height, figure fraction

def build_model(MODEL_DIR):
    """Return (rotation_model, PlateReconstruction, coastlines_path, continents_path)."""
    rot  = sorted(glob.glob(os.path.join(MODEL_DIR, "Rotations", "*.rot")) +
                  glob.glob(os.path.join(MODEL_DIR, "DeformingMeshes", "*.rot")))
    topo = sorted(glob.glob(os.path.join(MODEL_DIR, "DeformingMeshes", "*.gpml")) +
                  glob.glob(os.path.join(MODEL_DIR, "PlateBoundaries", "*.gpml")))
    static = os.path.join(MODEL_DIR, "StaticPolygons",
                          "Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons.gpml")
    coast  = os.path.join(MODEL_DIR, "Coastlines", "Clennett__etal_2020_Coastlines.gpml")
    cont   = os.path.join(MODEL_DIR, "ContinentalPolygons",
                          "Global_PresentDay_ContPolygons_2019_v2.gpml")
    rotm = pygplates.RotationModel(rot)
    model = gplately.PlateReconstruction(rotm, topology_features=topo,
                                         static_polygons=static, anchor_plate_id=0)
    return rotm, model, coast, cont, static

def gmt_cpt(path):
    """Parse a GMT .cpt into (ListedColormap, BoundaryNorm) with NaN->gray."""
    bounds, colors, over, under = [], [], None, None
    for line in open(path):
        s = line.split()
        if not s or s[0].startswith("#"):
            continue
        if s[0] == "B": under = tuple(int(x)/255 for x in s[1].split("/")); continue
        if s[0] == "F": over  = tuple(int(x)/255 for x in s[1].split("/")); continue
        if s[0] == "N": continue
        try: z0, z1 = float(s[0]), float(s[2])
        except Exception: continue
        c = tuple(int(x)/255 for x in s[1].split("/"))
        if not bounds: bounds.append(z0)
        bounds.append(z1); colors.append(c)
    cmap = ListedColormap(colors); cmap.set_over(over); cmap.set_under(under)
    cmap.set_bad(CONTINENT_GRAY)
    return cmap, BoundaryNorm(bounds, cmap.N)

def _grid_south_up(path):
    """Load a 2-D grid as south-up on a -180..180 longitude axis.

    The carbonate-thickness grids are already -180..180, but the Alfonso 2024
    pyBacktrack paleobathymetry grids are stored on a 0..360 axis. Both are
    plotted with extent=[-180,180,-90,90], so a 0..360 grid must be rolled by
    half a globe first - otherwise the raster lands 180 degrees away from the
    reconstructed overlays and the frame shows two sets of continents.
    """
    d = xr.open_dataset(path)
    v = "z" if "z" in d.data_vars else list(d.data_vars)[0]
    da = d[v]
    latname = "lat" if "lat" in da.coords else [c for c in da.coords if da[c].ndim == 1 and "lat" in c.lower()][0]
    lonname = "lon" if "lon" in da.coords else [c for c in da.coords if da[c].ndim == 1 and "lon" in c.lower()][0]
    lat = da[latname].values
    lon = da[lonname].values
    z = da.values
    if lat[0] > lat[-1]:               # north-first -> flip to south-up
        z = z[::-1, :]
    if float(np.nanmax(lon)) > 180.0:  # 0..360 -> -180..180
        if np.isclose(lon[-1] - lon[0], 360.0):   # drop the duplicated wrap column
            z = z[:, :-1]; lon = lon[:-1]
        k = int(np.searchsorted(lon, 180.0))
        z = np.concatenate([z[:, k:], z[:, :k]], axis=1)
    return z

def render_video(times, grid_fmt, out_dir, out_prefix, cbar_label,
                 MODEL_DIR, cmap, norm=None, vmin=None, vmax=None,
                 central_lon=15.0, framerate=8, dpi=150,
                 cbar_ticks=None, cbar_extend="max", figsize=(10.0, 5.8),
                 force=False):
    rotm, model, coast, cont, static = build_model(MODEL_DIR)
    frame_dir = os.path.join(out_dir, f"frames_{out_prefix}")
    os.makedirs(frame_dir, exist_ok=True)
    proj = ccrs.Mollweide(central_longitude=central_lon)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    if norm is None: sm.set_clim(vmin, vmax)

    n_done = 0
    for T in times:
        out = os.path.join(frame_dir, f"frame_{int(T):04d}.png")
        if os.path.exists(out) and not force:
            n_done += 1; continue
        gp_path = grid_fmt.format(t=int(T))
        if not os.path.exists(gp_path):
            print(f"  ! {T} Ma grid missing: {gp_path}"); continue
        z = _grid_south_up(gp_path)
        fig = plt.figure(figsize=figsize)
        # Layout. The Mollweide globe is width-limited, so its height is half the map
        # box width; size the box to the globe and hang the colour bar CBAR_GAP_MM
        # below it. The old fixed boxes (map bottom 0.10, bar top 0.105) made the bar
        # overlap the bottom of the globe.
        _fig_mm = figsize[1] * 25.4
        _map_h  = (0.98 * figsize[0] / 2.0) / figsize[1]
        _map_b  = 0.985 - _map_h
        ax = fig.add_axes([0.01, _map_b, 0.98, _map_h], projection=proj); ax.set_global()
        ax.spines["geo"].set_linewidth(0.7)
        imkw = dict(origin="lower", extent=[-180, 180, -90, 90],
                    transform=ccrs.PlateCarree(), zorder=1)
        if norm is not None: ax.imshow(z, cmap=cmap, norm=norm, **imkw)
        else:               ax.imshow(z, cmap=cmap, vmin=vmin, vmax=vmax, **imkw)
        gp = gplately.PlotTopologies(model, coastlines=coast, continents=cont, time=T)
        gp.plot_continents(ax, facecolor=CONTINENT_GRAY, edgecolor="none", zorder=2)
        gp.plot_coastlines(ax, color="0.4", linewidth=0.22, zorder=3)
        for w, c, zz in [(1.5, "white", 4), (0.65, "black", 5)]:
            gp.plot_topological_plate_boundaries(ax, color=c, linewidth=w, zorder=zz)
        gp.plot_subduction_teeth(ax, color="black", zorder=6)
        ax.text(0.012, 0.97, f"{int(T)} Ma", transform=ax.transAxes, fontsize=15,
                fontweight="bold", va="top", ha="left", zorder=10,
                bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=1.0))
        cax = fig.add_axes([0.25, _map_b - CBAR_GAP_MM/_fig_mm - _CBAR_H, 0.5, _CBAR_H])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal", extend=cbar_extend,
                          spacing="proportional")
        cb.set_label(cbar_label, fontsize=11)
        if cbar_ticks is not None: cb.set_ticks(cbar_ticks)
        cb.ax.tick_params(labelsize=9)
        fig.savefig(out, dpi=dpi); plt.close(fig)
        n_done += 1
        if int(T) % 10 == 0: print(f"  {int(T):3d} Ma -> {os.path.basename(out)}")
    print(f"  frames ready: {n_done}")

    if shutil.which("ffmpeg") is None:
        print("  ffmpeg not on PATH; frames left in", frame_dir); return
    back = os.path.join(out_dir, f"{out_prefix}_back_in_time.mp4")
    fwd  = os.path.join(out_dir, f"{out_prefix}_forward_in_time.mp4")
    subprocess.run(["ffmpeg","-y","-framerate",str(framerate),
        "-i",os.path.join(frame_dir,"frame_%04d.png"),"-pix_fmt","yuv420p",
        "-vcodec","libx264","-crf","20","-vf","scale=trunc(iw/2)*2:trunc(ih/2)*2",back],check=True)
    subprocess.run(["ffmpeg","-y","-i",back,"-vf","reverse","-vcodec","libx264",
        "-crf","20","-pix_fmt","yuv420p",fwd],check=True)
    print("  wrote", back, "and", fwd)
