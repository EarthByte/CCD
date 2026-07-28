"""
Reconstructed coastlines helper.

Uses gplately + plate_model_manager to fetch the Alfonso 2024 plate model
(based on Clennett 2020 extending Müller 2019, 0-170 Ma, mantle reference
frame so anchor plate id = 0) and reconstruct present-day coastlines to a
target age t.  The result is returned as a numpy array suitable for
pyGMT's ``fig.plot(x=..., y=...)`` interface — multi-polygon segments are
separated by NaN rows so a single ``fig.plot`` call draws every polygon.

Cached so repeated calls at the same age don't re-trigger the gplately
reconstruction step.

Used by 01_render_videos.py (per-frame coastlines on top of the carbonate-
thickness rasters) and 06_delta_panel.py (coastlines on the 2x2 picked-
time difference panel).
"""
from __future__ import annotations

import functools
import os
import warnings
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Plate model and cache location.  The Alfonso 2024 model is fetched by
# plate-model-manager into a per-script cache so repeated runs share it.
# ---------------------------------------------------------------------------
PLATE_MODEL_NAME = "Alfonso2024"
PLATE_MODEL_ANCHOR = 0                  # mantle reference frame

HERE = Path(__file__).resolve().parent
PLATE_MODEL_CACHE = HERE / "output" / "plate-model-cache"


# ---------------------------------------------------------------------------
# Lazy initialisation -- the first call to `get_coastlines_xy(t)` builds
# the reconstruction driver once.  Subsequent calls reuse it.
# ---------------------------------------------------------------------------
_DRIVER = None
_AVAILABLE = None


def _init_driver():
    """Return a (rotation_model, coastline_features) tuple, or (None, None)
    if gplately / plate_model_manager / pygplates aren't installed.
    """
    global _DRIVER, _AVAILABLE
    if _AVAILABLE is False:
        return None, None
    if _DRIVER is not None:
        return _DRIVER
    try:
        import pygplates                                       # noqa: F401
        from plate_model_manager import PlateModelManager
    except ImportError as exc:
        warnings.warn(
            f"reconstructed-coastline overlay disabled "
            f"({exc.__class__.__name__}: {exc}); install gplately + "
            "plate-model-manager to enable.",
            stacklevel=2,
        )
        _AVAILABLE = False
        return None, None

    PLATE_MODEL_CACHE.mkdir(parents=True, exist_ok=True)
    pm = PlateModelManager()
    plate_model = pm.get_model(PLATE_MODEL_NAME, data_dir=str(PLATE_MODEL_CACHE))
    rotation_model = plate_model.get_rotation_model()
    coastlines = plate_model.get_coastlines()
    _DRIVER = (rotation_model, coastlines)
    _AVAILABLE = True
    return _DRIVER


@functools.lru_cache(maxsize=200)
def get_coastlines_xy(t: int) -> np.ndarray | None:
    """Reconstruct the Alfonso 2024 coastlines to age t and return them
    as a single ``(N, 2)`` ``[lon, lat]`` array with NaN-row separators
    between polygons.  Returns None if the plate model isn't available.

    Pass straight to ``fig.plot(x=arr[:, 0], y=arr[:, 1], pen=...)`` -- the
    NaN rows make pyGMT lift the pen between polygons.
    """
    rotation_model, coastlines = _init_driver()
    if rotation_model is None:
        return None

    import pygplates
    reconstructed = []
    pygplates.reconstruct(
        coastlines, rotation_model, reconstructed, float(t),
        anchor_plate_id=PLATE_MODEL_ANCHOR,
    )

    segments = []
    for feat in reconstructed:
        geom = feat.get_reconstructed_geometry()
        if geom is None:
            continue
        # to_lat_lon_array() returns (N, 2) with cols = (lat, lon)
        ll = geom.to_lat_lon_array()
        if ll.size == 0:
            continue
        # swap to (lon, lat) for pyGMT's x/y convention
        lonlat = np.column_stack([ll[:, 1], ll[:, 0]])
        # close polygons so the last segment draws
        if not np.allclose(lonlat[0], lonlat[-1]):
            lonlat = np.vstack([lonlat, lonlat[0]])
        segments.append(lonlat)
        segments.append(np.array([[np.nan, np.nan]]))      # pen-up separator

    if not segments:
        return np.empty((0, 2))
    return np.vstack(segments)


def plot_on(fig, t: int, *, projection: str, region: str | list,
            pen: str = "0.4p,gray30"):
    """Convenience wrapper -- reconstruct + plot in one call.  Silently
    no-ops if the plate model isn't available.
    """
    xy = get_coastlines_xy(t)
    if xy is None or len(xy) == 0:
        return
    fig.plot(x=xy[:, 0], y=xy[:, 1],
             pen=pen, projection=projection, region=region)
