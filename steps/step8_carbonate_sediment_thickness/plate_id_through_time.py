"""
plate_id_through_time.py
========================

For a given present-day (lon, lat) point, report the plate ID assigned by a
topological plate model at each 1-My time step between a user-supplied start
and end time.

Approach
--------
1.  Load the same rotation + topology files used by
    ``carbonate_sediment_thickness_2026.ipynb`` (default: ``Alfonso2024`` via
    ``plate-model-manager``; or supply local files via ``--rotation-files`` /
    ``--topology-files``).
2.  Resolve topologies at present day (t = 0) and find the resolved-topology
    polygon that contains the input lon/lat.  The plate ID of that polygon is
    the "present-day" / carrier plate ID (P0).
3.  For each time ``t`` from ``start_time`` to ``end_time`` in 1-My steps:
        a.  Reconstruct the input lon/lat back to time ``t`` using plate ID
            P0 (the point "moves with" its present-day plate).
        b.  Resolve topologies at time ``t`` and look up the plate ID (Pt) of
            the polygon containing the reconstructed paleo location.
4.  Print one row per time step to stdout and also write a CSV.

Notes
-----
*   The reconstruction uses ``carbonate_anchor_plate_id = 0`` (the same anchor
    plate used in the notebook for the carbonate reference frame).
*   If the point at the given paleo time falls in a "sliver" outside all
    resolved topologies, ``Pt`` is reported as 0 (matching the notebook's
    fallback in ``carbonate_sediment_thickness.py``).
*   ``end_time`` may be either greater or less than ``start_time``; the
    script always iterates inclusively in 1-My steps.

Usage
-----
    python plate_id_through_time.py LON LAT START_TIME END_TIME [options]

Examples
--------
    # Mid-Atlantic point, 0 to 100 Ma, default Alfonso2024 model:
    python plate_id_through_time.py -30 0 0 100

    # Use a local model:
    python plate_id_through_time.py -30 0 0 170 \\
        --rotation-files input_data/topology_model/2019_v2/*.rot \\
        --topology-files input_data/topology_model/2019_v2/*.gpmlz

    # Custom output CSV:
    python plate_id_through_time.py 150 -30 0 80 -o my_point.csv
"""

import argparse
import csv
import glob
import os
import sys

import pygplates


# Default carbonate-reference-frame anchor plate ID, matching the notebook.
DEFAULT_ANCHOR_PLATE_ID = 0


# --------------------------------------------------------------------------- #
# Topological-model loading                                                   #
# --------------------------------------------------------------------------- #
def load_rotation_and_topology_files(args):
    """
    Returns (rotation_filenames, topology_filenames) lists.

    Resolution order:
        1.  --rotation-files / --topology-files (explicit local files)
        2.  --plate-model (downloaded via plate-model-manager)
    """
    if args.rotation_files or args.topology_files:
        if not (args.rotation_files and args.topology_files):
            raise SystemExit(
                "Must supply BOTH --rotation-files and --topology-files when "
                "using local files."
            )

        # Expand any glob patterns the shell did not.
        rotation_filenames = []
        for pattern in args.rotation_files:
            matched = glob.glob(pattern)
            rotation_filenames.extend(matched if matched else [pattern])
        topology_filenames = []
        for pattern in args.topology_files:
            matched = glob.glob(pattern)
            topology_filenames.extend(matched if matched else [pattern])

        if not rotation_filenames:
            raise SystemExit("No rotation files found.")
        if not topology_filenames:
            raise SystemExit("No topology files found.")

        return rotation_filenames, topology_filenames

    # Fall back to plate-model-manager (same as the notebook).
    try:
        from plate_model_manager import PlateModelManager
    except ImportError as e:
        raise SystemExit(
            "plate-model-manager is required when no local files are supplied. "
            "Install it with `pip install plate-model-manager`, or pass "
            "--rotation-files / --topology-files."
        ) from e

    pmm = PlateModelManager()
    plate_model = pmm.get_model(args.plate_model, data_dir=args.plate_model_dir)
    return plate_model.get_rotation_model(), plate_model.get_topologies()


# --------------------------------------------------------------------------- #
# Plate-ID lookup                                                             #
# --------------------------------------------------------------------------- #
def plate_id_at_point(point, topology_filenames, rotation_model, time,
                      anchor_plate_id=DEFAULT_ANCHOR_PLATE_ID):
    """
    Resolve topologies at ``time`` (with ``anchor_plate_id`` as the reference
    frame) and return the plate ID of the resolved-topology polygon that
    contains ``point``.  Returns 0 if the point falls outside all topologies.
    """
    resolved_topologies = []
    pygplates.resolve_topologies(
        topology_filenames,
        rotation_model,
        resolved_topologies,
        time,
        anchor_plate_id=anchor_plate_id,
    )

    for resolved_topology in resolved_topologies:
        boundary = resolved_topology.get_resolved_boundary()
        if boundary.is_point_in_polygon(point):
            return resolved_topology.get_feature().get_reconstruction_plate_id()

    return 0  # Fallback used in carbonate_sediment_thickness.py


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Report plate ID through time for a present-day (lon, lat) point, "
            "using a topological plate model."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("lon", type=float, help="Present-day longitude (-180, 180]")
    parser.add_argument("lat", type=float, help="Present-day latitude  [-90, 90]")
    parser.add_argument("start_time", type=float, help="Start time (Ma, >= 0)")
    parser.add_argument("end_time", type=float, help="End time   (Ma, >= 0)")
    parser.add_argument(
        "-o", "--output-csv",
        default=None,
        help=("CSV output file (default: "
              "plate_id_through_time_<lon>_<lat>_<start>_<end>.csv "
              "in the current directory)"),
    )
    parser.add_argument(
        "--plate-model",
        default="Alfonso2024",
        help="plate-model-manager model name (default: Alfonso2024)",
    )
    parser.add_argument(
        "--plate-model-dir",
        default="plate-model-repo",
        help="Directory to cache the downloaded plate model (default: plate-model-repo)",
    )
    parser.add_argument(
        "--rotation-files",
        nargs="+",
        default=None,
        help="Explicit list of local rotation files (.rot). "
             "If supplied you must also pass --topology-files.",
    )
    parser.add_argument(
        "--topology-files",
        nargs="+",
        default=None,
        help="Explicit list of local topology files (.gpml/.gpmlz). "
             "If supplied you must also pass --rotation-files.",
    )
    parser.add_argument(
        "--anchor-plate-id",
        type=int,
        default=DEFAULT_ANCHOR_PLATE_ID,
        help=(f"Anchor plate ID for resolving topologies and reconstructing "
              f"(default: {DEFAULT_ANCHOR_PLATE_ID}, matching the notebook's "
              f"carbonate reference frame)."),
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Don't print per-time-step rows to stdout (still writes CSV).",
    )

    args = parser.parse_args()

    # ---------------------- input validation ------------------------------- #
    if not (-180.0 <= args.lon <= 360.0):
        parser.error("lon must be in [-180, 360]")
    if not (-90.0 <= args.lat <= 90.0):
        parser.error("lat must be in [-90, 90]")
    if args.start_time < 0 or args.end_time < 0:
        parser.error("start_time and end_time must be >= 0 Ma")

    # Build inclusive 1-My step list (works whether start <= end or start > end).
    if args.end_time >= args.start_time:
        times = [args.start_time + i for i in range(int(args.end_time - args.start_time) + 1)]
    else:
        times = [args.start_time - i for i in range(int(args.start_time - args.end_time) + 1)]

    # ---------------------- load model ------------------------------------- #
    print("Loading rotation and topology files...", file=sys.stderr)
    rotation_filenames, topology_filenames = load_rotation_and_topology_files(args)
    print(f"  {len(rotation_filenames)} rotation file(s), "
          f"{len(topology_filenames)} topology file(s)", file=sys.stderr)
    rotation_model = pygplates.RotationModel(rotation_filenames)

    # ---------------------- present-day plate ID --------------------------- #
    # PointOnSphere takes (lat, lon) — same convention as the notebook.
    present_point = pygplates.PointOnSphere(args.lat, args.lon)

    print(f"Determining present-day plate ID at "
          f"(lon={args.lon}, lat={args.lat}) ...", file=sys.stderr)
    present_plate_id = plate_id_at_point(
        present_point,
        topology_filenames,
        rotation_model,
        time=0,
        anchor_plate_id=args.anchor_plate_id,
    )
    print(f"  present-day (carrier) plate ID = {present_plate_id}", file=sys.stderr)

    # ---------------------- output setup ----------------------------------- #
    if args.output_csv is None:
        args.output_csv = (f"plate_id_through_time_"
                           f"{args.lon}_{args.lat}_"
                           f"{int(args.start_time)}_{int(args.end_time)}.csv")

    out_dir = os.path.dirname(args.output_csv)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    header = [
        "time_Ma",
        "input_lon", "input_lat",
        "reconstructed_lon", "reconstructed_lat",
        "carrier_plate_id",          # P0 (present-day plate ID, constant)
        "plate_id_at_time",          # Pt (topology plate ID at paleo location)
    ]

    if not args.quiet:
        print("\t".join(header))

    rows = []

    # ---------------------- iterate over times ----------------------------- #
    for t in times:
        if t == 0:
            recon_lat, recon_lon = args.lat, args.lon
        else:
            # Rotation that takes the present-day point on plate `present_plate_id`
            # back to its paleo location at time `t` (in the same anchor plate frame).
            stage_rotation = rotation_model.get_rotation(
                t, present_plate_id, 0, anchor_plate_id=args.anchor_plate_id,
            )
            paleo_point = stage_rotation * present_point
            recon_lat, recon_lon = paleo_point.to_lat_lon()

        # Plate ID of the resolved topology polygon containing the paleo location at time t.
        if t == 0:
            plate_id_at_t = present_plate_id
            paleo_point_for_lookup = present_point
        else:
            paleo_point_for_lookup = pygplates.PointOnSphere(recon_lat, recon_lon)
            plate_id_at_t = plate_id_at_point(
                paleo_point_for_lookup,
                topology_filenames,
                rotation_model,
                time=t,
                anchor_plate_id=args.anchor_plate_id,
            )

        row = [
            f"{t:g}",
            f"{args.lon:g}", f"{args.lat:g}",
            f"{recon_lon:.6f}", f"{recon_lat:.6f}",
            present_plate_id,
            plate_id_at_t,
        ]
        rows.append(row)

        if not args.quiet:
            print("\t".join(str(c) for c in row))

    # ---------------------- write CSV -------------------------------------- #
    with open(args.output_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} rows to: {args.output_csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
