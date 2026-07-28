"""
Reconstructed continental polygon AND coastline outlines helper
(local-files-only, no plate-model-manager fetch).

Loads the local Alfonso et al. (2024) modified Clennett-Müller plate model
directly via pygplates from the files shipped under
    input/Alfonso_etal_2024_modClennettMuller/Rotations/*.rot
    input/Alfonso_etal_2024_modClennettMuller/ContinentalPolygons/*.gpml
    input/Alfonso_etal_2024_modClennettMuller/Coastlines/*.gpml

Both rotation files are loaded into a single ``pygplates.RotationModel`` as
required by the model documentation (the Clennett file holds the western
North America refinements, the global 250-0 Ma file holds everything else).
The same RotationModel is shared between the continental-polygon and the
coastline reconstructions, so the two overlays line up exactly and the
rotations are read off disk only once per process.

Two functions are exposed:
    get_continents_xy(t), plot_continents_on(fig, t, ...)  - continental crust outlines
    get_coastlines_xy(t), plot_coastlines_on(fig, t, ...)  - modern shoreline outlines

Both return a single ``(N, 2)`` ``[lon, lat]`` array with NaN-row separators
between polygons, suitable for a single ``fig.plot(x=arr[:, 0], y=arr[:, 1])``
call (the NaN rows make pyGMT lift the pen between polygons).

For backwards compatibility ``plot_on`` is kept as an alias of
``plot_continents_on``.

This module is intentionally separate from ``paleo_coastlines.py``. The
latter is the original plate-model-manager-backed helper used by
``01_render_videos.py`` and ``06_delta_panel.py``. The script
``07_render_paleobathymetry_video.py`` uses this local-files-only
implementation instead so it runs entirely offline against the files the
user supplied.
"""
from __future__ import annotations

import functools
import warnings
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Local plate model files. Both rotation files MUST be loaded into the same
# RotationModel - the Clennett file contains the western-North-America
# refinements, the global file contains everything else, and pygplates
# resolves cross-references between them only when both are present in the
# same model (see _init() below where the merge happens).
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
PLATE_MODEL_DIR = HERE / "input" / "Alfonso_etal_2024_modClennettMuller"

ROTATION_FILES = [
    PLATE_MODEL_DIR / "Rotations" / "Clennett_etal_2020_Rotations.rot",
    PLATE_MODEL_DIR / "Rotations" / "Global_250-0Ma_Rotations.rot",
]
CONTPOLYGONS_FILE = (
    PLATE_MODEL_DIR
    / "ContinentalPolygons"
    # Note: the v2 file shipped here has been edited in place - the South
    # Cuba polygon's plate ID was changed from 2022 -> 254 (the original
    # upstream v2, which had 2022, was deleted). Don't replace this file
    # with a fresh upstream v2 without re-applying that edit.
    / "Global_PresentDay_ContPolygons_2019_v2.gpml"
)
COASTLINES_FILE = (
    PLATE_MODEL_DIR
    / "Coastlines"
    / "Clennett__etal_2020_Coastlines.gpml"
)

# Mantle reference frame, matching the carbonate paleobathymetry input grids.
PLATE_MODEL_ANCHOR = 0


# ---------------------------------------------------------------------------
# Lazy initialisation - the first call to either get_*_xy(t) function builds
# the shared RotationModel and both FeatureCollections once. Subsequent
# calls reuse the cached state.
# ---------------------------------------------------------------------------
_ROTATION_MODEL = None
_CONT_FEATURES = None
_COAST_FEATURES = None
_AVAILABLE: bool | None = None


def _init():
    """Return True if pygplates is importable and all input files exist.

    Builds a single ``pygplates.RotationModel`` from BOTH rotation files
    listed in ``ROTATION_FILES``. This is the connection step the user
    asked for: passing the list (rather than a single file) makes
    pygplates merge the two rotation hierarchies so that any plate ID
    defined in one file can be resolved against parents defined in the
    other. Reconstructing with only one of the two would silently drop
    plate IDs and produce wrong geometries.
    """
    global _ROTATION_MODEL, _CONT_FEATURES, _COAST_FEATURES, _AVAILABLE
    if _AVAILABLE is False:
        return False
    if _ROTATION_MODEL is not None:
        return True

    try:
        import pygplates
    except ImportError as exc:
        warnings.warn(
            f"local-plate-model overlay disabled "
            f"({exc.__class__.__name__}: {exc}); install pygplates to enable.",
            stacklevel=2,
        )
        _AVAILABLE = False
        return False

    required = (*ROTATION_FILES, CONTPOLYGONS_FILE, COASTLINES_FILE)
    missing = [p for p in required if not p.exists()]
    if missing:
        warnings.warn(
            "local-plate-model overlay disabled; missing input file(s): "
            + ", ".join(str(p) for p in missing),
            stacklevel=2,
        )
        _AVAILABLE = False
        return False

    # Merge BOTH .rot files into a single RotationModel so the plate-ID
    # hierarchy is fully connected. This is the key step: pygplates only
    # cross-resolves rotations defined across multiple files when they are
    # loaded together as a single model.
    _ROTATION_MODEL = pygplates.RotationModel([str(p) for p in ROTATION_FILES])
    _CONT_FEATURES = pygplates.FeatureCollection(str(CONTPOLYGONS_FILE))
    _COAST_FEATURES = pygplates.FeatureCollection(str(COASTLINES_FILE))
    _AVAILABLE = True
    return True


def _reconstruct_to_segments(features, t: int) -> np.ndarray:
    """Reconstruct ``features`` to age ``t`` using the merged local
    RotationModel and return an (N, 2) [lon, lat] array with NaN-row
    separators between polygons. Internal helper shared by the
    continents and coastlines paths.
    """
    import pygplates

    reconstructed = []
    pygplates.reconstruct(
        features, _ROTATION_MODEL, reconstructed, float(t),
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


@functools.lru_cache(maxsize=200)
def get_continents_xy(t: int) -> np.ndarray | None:
    """Reconstruct the continental polygons to age ``t`` Ma and return
    them as an (N, 2) [lon, lat] array with NaN-row separators between
    polygons. Returns None if the plate model isn't available.
    """
    if not _init():
        return None
    return _reconstruct_to_segments(_CONT_FEATURES, t)


@functools.lru_cache(maxsize=200)
def get_coastlines_xy(t: int) -> np.ndarray | None:
    """Reconstruct the (local) coastlines to age ``t`` Ma and return them
    as an (N, 2) [lon, lat] array with NaN-row separators between
    polygons. Returns None if the plate model isn't available.
    """
    if not _init():
        return None
    return _reconstruct_to_segments(_COAST_FEATURES, t)


def plot_continents_on(fig, t: int, *, projection: str, region: str | list,
                       pen: str = "0.6p,red"):
    """Convenience wrapper - reconstruct continental polygons + plot in
    one call. Silently no-ops if the plate model isn't available.
    """
    xy = get_continents_xy(t)
    if xy is None or len(xy) == 0:
        return
    fig.plot(x=xy[:, 0], y=xy[:, 1],
             pen=pen, projection=projection, region=region)


def plot_coastlines_on(fig, t: int, *, projection: str, region: str | list,
                       pen: str = "0.4p,green"):
    """Convenience wrapper - reconstruct (local) coastlines + plot in one
    call. Silently no-ops if the plate model isn't available.
    """
    xy = get_coastlines_xy(t)
    if xy is None or len(xy) == 0:
        return
    fig.plot(x=xy[:, 0], y=xy[:, 1],
             pen=pen, projection=projection, region=region)


# Backwards-compatible alias - earlier revisions exposed only ``plot_on``
# for the continental-polygon case.
plot_on = plot_continents_on
