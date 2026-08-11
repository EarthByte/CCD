"""
Standard I/O helpers.

Every tabular product in this workflow is a whitespace/tab-delimited text file
with an optional leading ``#``-commented header. These helpers give every step a
single, consistent way to read and write such files so that output file endings
and formats no longer drift between ``.tsv`` / ``.txt`` / ``.xy`` / ``.csv``.

Convention
----------
* On disk: one ``.txt`` file per product, tab-separated, with a single
  ``# col1<TAB>col2 ...`` header line.
* In memory: a :class:`pandas.DataFrame` with explicit column names.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


def read_table(path: str | Path, names: Sequence[str] | None = None) -> pd.DataFrame:
    """Read a commented, whitespace/comma/tab-delimited numeric table.

    Comment lines (starting with ``#``) and any header row are skipped; columns
    are assigned from ``names`` (or ``col0, col1, ...`` if not given). Rows that
    cannot be parsed as numbers are dropped.
    """
    df = pd.read_csv(
        path,
        sep=r"\s+|\t+|,",
        engine="python",
        comment="#",
        header=None,
    )
    if names is not None:
        df = df.iloc[:, : len(names)]
        df.columns = list(names)
    else:
        df.columns = [f"col{i}" for i in range(df.shape[1])]
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(how="all").reset_index(drop=True)


def read_xy(path: str | Path, names: Sequence[str] = ("Age_Ma", "Value")) -> pd.DataFrame:
    """Read a two-column (age, value) table, sorted by the first column."""
    df = read_table(path, names=names).dropna().reset_index(drop=True)
    return df.sort_values(list(names)[0]).reset_index(drop=True)


def write_table(
    path: str | Path,
    df: pd.DataFrame,
    header: Sequence[str] | None = None,
    float_format: str = "%.3f",
) -> Path:
    """Write ``df`` as a tab-separated ``.txt`` file with a ``# ``-commented header.

    Parameters
    ----------
    path : output path (``.txt`` enforced).
    df : data to write.
    header : column labels for the ``# ``-commented header line. Defaults to
        the DataFrame column names.
    float_format : printf-style format for floating point values.
    """
    path = Path(path).with_suffix(".txt")
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(header) if header is not None else list(df.columns)
    with open(path, "w") as f:
        f.write("# " + "\t".join(str(c) for c in cols) + "\n")
        for _, row in df.iterrows():
            vals = []
            for v in row:
                if isinstance(v, (int, np.integer)):
                    vals.append(str(int(v)))
                elif v is None or (isinstance(v, float) and not np.isfinite(v)):
                    vals.append("NaN")
                else:
                    vals.append(float_format % float(v))
            f.write("\t".join(vals) + "\n")
    return path


def _figure_targets(name: str, step: str | None, paper_name: str | None = None):
    """Return (png, pdf) path pairs a figure should be written to: always the
    shared ``Figures/`` folder, plus the per-step ``outputs/<step>/figures`` folder
    when ``step`` is given, plus the paper's figure folder under ``paper_name``
    when that is given. Directories are created as needed.
    """
    from . import config
    targets = [(config.FIGURES, name)]
    if step:
        targets.append((config.step_figures_dir(step), name))
    if paper_name:
        if config.PAPER_FIGURES is None:
            print(f"  [figures] paper figure folder not found; "
                  f"{paper_name} not propagated (set CCD_PAPER_FIGURES)")
        else:
            targets.append((config.PAPER_FIGURES, paper_name))
    pairs = []
    for d, stem in targets:
        d.mkdir(parents=True, exist_ok=True)
        pairs.append((d / f"{stem}.png", d / f"{stem}.pdf"))
    return pairs


def save_matplotlib_figure(fig, name: str, dpi: int = 300, step: str | None = None,
                           paper_name: str | None = None):
    """Save a matplotlib figure as both PNG and PDF at ``dpi``. Always writes to the
    shared ``Figures/`` folder; if ``step`` is given (e.g. ``"step5"``) it is also
    written to that step's ``outputs/<step>/figures`` folder. ``name`` is a bare
    stem (no extension). Returns the shared-folder (png, pdf) paths.
    """
    pairs = _figure_targets(name, step, paper_name)
    for png, pdf in pairs:
        fig.savefig(png, dpi=dpi, bbox_inches="tight")
        fig.savefig(pdf, dpi=dpi, bbox_inches="tight")
    return pairs[0]


def save_pygmt_figure(fig, name: str, dpi: int = 300, step: str | None = None):
    """Save a pyGMT figure as both PNG and PDF at ``dpi``. Always writes to the
    shared ``Figures/`` folder; if ``step`` is given it is also written to that
    step's ``outputs/<step>/figures`` folder. ``name`` is a bare stem (no
    extension). Returns the shared-folder (png, pdf) paths.
    """
    pairs = _figure_targets(name, step)
    for png, pdf in pairs:
        fig.savefig(str(pdf))
        fig.savefig(str(png), dpi=dpi)
    return pairs[0]


def write_polygon(path: str | Path, age: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> Path:
    """Write a closed (age, y) uncertainty-envelope polygon as a ``.txt`` file.

    Lower bound is written forward, upper bound backward, then the ring is
    closed — the format expected by GMT/pyGMT ``plot`` with a fill.
    """
    path = Path(path).with_suffix(".txt")
    path.parent.mkdir(parents=True, exist_ok=True)
    age = np.asarray(age, float)
    lo = np.asarray(lo, float)
    hi = np.asarray(hi, float)
    with open(path, "w") as f:
        f.write("# Age_Ma\tCCD_m  (closed uncertainty-envelope polygon)\n")
        for xi, yi in zip(age, lo):
            f.write(f"{xi:.6f}\t{yi:.6f}\n")
        for xi, yi in zip(age[::-1], hi[::-1]):
            f.write(f"{xi:.6f}\t{yi:.6f}\n")
        f.write(f"{age[0]:.6f}\t{lo[0]:.6f}\n")
    return path
