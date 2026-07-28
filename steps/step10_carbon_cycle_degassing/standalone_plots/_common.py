"""
Shared configuration + CSV loaders for the standalone sedimentary-carbon
plotting scripts.

Each plotting script in this folder uses these helpers to load the same set
of CSV files for both the DM2026 and BW1991 models, so the plots stay
strictly comparable across the two scenarios.

Folder layout (input):
    ../output/<MODEL>_carbon/Notebook02/csv/02_plate_influx.csv
    ../output/<MODEL>_carbon/Notebook02/csv/02_subducted_carbon.csv
    ../output/<MODEL>_carbon/Notebook05/csv/05_atmospheric_influx_all_sources.csv
    ../output/<MODEL>_carbon/Notebook05/csv/05_carbon_flux.csv
    ../output/<MODEL>_carbon/Notebook05/csv/05_global_temperature.csv
    ../output/<MODEL>_carbon/Notebook05/figures/05_overriding_plate_storage.csv
    ../utils/GAST_800Ma_Scotese.xlsx                       (Scotese GAST curve)

Folder layout (output):
    ./output/<plot_name>_<MODEL>.png/pdf
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent                              # 9_CO2_Alfonso2024/
# Read the CURRENT run directly (the old per-model output/ duplicate is retired).
DM26_OUTPUTS = PROJECT_ROOT / "Alfonso_etal_2024_DM26" / "Outputs"
UTILS = PROJECT_ROOT / "utils"
SCOTESE_XLSX = UTILS / "GAST_800Ma_Scotese.xlsx"

PLOTS_OUT = HERE / "output"
PLOTS_OUT.mkdir(parents=True, exist_ok=True)

MODELS = ("DM2026",)   # BW1991 (superseded old-CCD comparison) retired: only the current run


# ---------------------------------------------------------------------------
# Per-model CSV path helpers
# ---------------------------------------------------------------------------
def _model_dir(model: str) -> Path:
    if model not in MODELS:
        raise ValueError(f"unknown model {model!r}; expected one of {MODELS}")
    return DM26_OUTPUTS   # current run (Alfonso_etal_2024_DM26/Outputs)


def csv_path(model: str, notebook: str, name: str) -> Path:
    """e.g. csv_path('DM2026', 'Notebook02', '02_plate_influx.csv')"""
    return _model_dir(model) / notebook / "csv" / name


def fig_csv_path(model: str, notebook: str, name: str) -> Path:
    """A handful of CSVs live in the figures/ folder rather than csv/
    (notably 05_overriding_plate_storage.csv)."""
    return _model_dir(model) / notebook / "figures" / name


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def load_plate_influx(model: str) -> pd.DataFrame:
    """02_plate_influx.csv — multi-index columns
    (component, {min,mean,max}); index = age (Ma).
    Components: sediments, crust, lithosphere, serpentinite_total,
                serpentinite_bending, serpentinite_mor, organic_sediments,
                total_influx
    Units: Mt C / yr (already scaled inside the source notebook).
    """
    df = pd.read_csv(csv_path(model, "Notebook02", "02_plate_influx.csv"),
                     header=[0, 1], index_col=0)
    df.index = df.index.astype(float)
    df.index.name = "age"
    return df


def load_subducted_carbon(model: str) -> pd.DataFrame:
    """02_subducted_carbon.csv — same multi-index column shape.
    Components: sediments, crust, lithosphere, serpentinite,
                serpentinite_bending, serpentinite_mor, organic_sediments,
                total, cumulative_subducted
    Units: Mt C / yr.
    """
    df = pd.read_csv(csv_path(model, "Notebook02", "02_subducted_carbon.csv"),
                     header=[0, 1], index_col=0)
    df.index = df.index.astype(float)
    df.index.name = "age"
    return df


def load_atmospheric_influx(model: str) -> pd.DataFrame:
    """05_atmospheric_influx_all_sources.csv — flat columns indexed by age.
    Columns include ridge/subduction/carbonate_platform/rift_*/intraplate_*
    outflux components (each as _min/_mean/_max), plus the various
    gross_atmospheric_outflux_* aggregates and the gross_upper_plate_influx_*
    that the net-atmospheric-influx plot needs.
    """
    df = pd.read_csv(csv_path(model, "Notebook05",
                              "05_atmospheric_influx_all_sources.csv"))
    df = df.set_index("Age (Ma)")
    return df


def load_net_carbon_outflux(model: str) -> pd.DataFrame:
    """05_net_carbon_outflux.csv - multi-index (variable, {min,mean,max}).
    ``net_carbon_outflux`` = gross atmospheric outflux minus the upper-plate
    storage computed WITHOUT the carbonate-sediment reservoir (serpentinite +
    crust only); i.e. the growing pelagic carbonate reservoir is treated as a
    shallow-to-deep redistribution, not a new sink."""
    df = pd.read_csv(csv_path(model, "Notebook05", "05_net_carbon_outflux.csv"),
                     header=[0, 1], index_col=0)
    df.index = df.index.astype(float); df.index.name = "age"
    return df


def load_carbon_flux(model: str) -> pd.DataFrame:
    """05_carbon_flux.csv — multi-index columns (variable, {min,mean,max}).
    Carries both the standard and Ratschbacher variants of the gross/net
    atmospheric influx series plus GAST. This is the source the
    170-0_diff_atmospheric_flux_* panels pull from.
    """
    df = pd.read_csv(csv_path(model, "Notebook05", "05_carbon_flux.csv"),
                     header=[0, 1], index_col=0)
    df.index = df.index.astype(float)
    df.index.name = "age"
    return df


def load_global_temperature(model: str) -> pd.DataFrame:
    """05_global_temperature.csv — flat columns (GAST_min/mean/max +
    average_min/mean/max), index = age (Ma). 'GAST' = Mills curve."""
    df = pd.read_csv(csv_path(model, "Notebook05", "05_global_temperature.csv"))
    df = df.set_index("Age (Ma)")
    return df


def load_overriding_plate_storage(model: str) -> pd.DataFrame:
    """05_overriding_plate_storage.csv (lives in figures/, not csv/).
    Provides slab_storage_*, overriding_plate_storage_* and
    subduction_atmospheric_influx_* — the three series the reservoir-flux
    plot draws."""
    df = pd.read_csv(fig_csv_path(model, "Notebook05",
                                  "05_overriding_plate_storage.csv"))
    df = df.set_index("Age (Ma)")
    return df


def load_scotese_gast() -> pd.DataFrame:
    """GAST_800Ma_Scotese.xlsx — columns: 'Age (Ma)', 'GAST (deg C) Scotese'."""
    df = pd.read_excel(SCOTESE_XLSX)
    df = df.rename(columns={"Age (Ma)": "age",
                            "GAST (deg C) Scotese": "scotese_gast"})
    return df


# ---------------------------------------------------------------------------
# Save helper - emits both PNG and PDF to ./output/
# ---------------------------------------------------------------------------
def save_figure(fig, stem: str, *, dpi: int = 300, also_pdf: bool = True):
    out_png = PLOTS_OUT / f"{stem}.png"
    fig.savefig(out_png, bbox_inches="tight", dpi=dpi)
    print(f"  wrote {out_png}")
    if also_pdf:
        out_pdf = PLOTS_OUT / f"{stem}.pdf"
        fig.savefig(out_pdf, bbox_inches="tight", dpi=dpi)
        print(f"  wrote {out_pdf}")
