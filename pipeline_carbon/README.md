# Carbonate + CO₂ batch pipeline (steps 0→7→8→9)

Part of the consolidated **CCD_workflow_clean** repository. One command runs the
whole downstream chain, non-interactively and DM2026-only, always starting from
the up-to-date **area-weighted** hybrid CCD so nothing stale propagates:

```bash
./run_all.sh                 # regenerate CCD -> grids -> stats -> CO2 -> figures
./run_all.sh --skip-ccd      # keep the current CCD; start at step 7
./run_all.sh --skip-step7    # reuse existing grids; stats + CO2 only
./run_all.sh --only-step9    # CO2 analysis + figures only
```

**Resume the CO₂ notebooks after a fixed error** (skip the notebooks that already
succeeded, and don't repeat steps 0/7/8):

```bash
CO2_START_FROM=02 ./run_all.sh --only-step9        # start at 02-Subducted-Carbon
# or run just the notebooks directly:
CO2_START_FROM=02-Subducted-Carbon.ipynb bash 03_run_co2_notebooks.sh
```

`CO2_START_FROM` accepts a notebook name or a prefix (`02`, `05`, …); every
notebook before it is skipped. The notebooks write their outputs to disk as they
run, so the ones already completed do not need to be re-executed.

> **macOS overnight runs — sleep-safe automatically.** Step 7 takes hours; if the
> Mac idle-sleeps the process is suspended and the run stalls. On macOS the script
> now **re-execs itself under `caffeinate -i`**, so a bare `./run_all.sh` never
> stalls on sleep — no wrapper to remember. The no-sleep assertion is released as
> soon as the run exits. Opt out with `./run_all.sh --no-caffeinate` (or
> `NO_CAFFEINATE=1`) if you're wrapping it yourself. Note: on battery,
> clamshell/lid-close sleep is *not* overridden by caffeinate — keep the Mac on AC
> power for overnight runs.

> **Runs in your local geoscience environment, not the cloud sandbox.** Plus the
> input grids (see *Data availability* below). It was wired but **not executed**
> here — test it and confirm the `>>> CONFIRM` paths in `config.sh`.

## Requirements (steps 7–9)

Steps 7–9 need a local geoscience Python environment (conda recommended). Beyond
the standard scientific stack (`numpy`, `pandas`, `scipy`, `matplotlib`), the
step-specific dependencies are:

| Package | Used by | Install |
|---|---|---|
| `pygplates`, `gplately` | step 7 grids; step 9 plate reconstructions | conda-forge |
| `pygmt` + GMT | step 7 gridding; pyGMT figures | `conda install -c conda-forge pygmt` |
| `xarray`, `netCDF4` | reading/writing the `.nc` grids (steps 7–9) | conda-forge |
| `jupyter` / `nbconvert` (or `papermill`) | headless step-9 notebook execution | conda-forge |
| **`pytables`** (imports as `tables`) | step-9 notebooks read HDF5 via `pd.read_hdf` | `conda install -c conda-forge pytables` |
| `joblib` | step-9 parallel loops (`Parallel`) | conda-forge |
| `plate-model-manager` | downloads the Alfonso2024 plate model (step 7) | `pip install plate-model-manager` |
| `openpyxl` | reading the `.xlsx` inputs | conda-forge |
| `moviepy` **≥ 2.0** | step-9 panel-plot / animation notebooks (video output) | `pip install -U "moviepy>=2.0"` |
| **`melt`** (Ben Mather) | step-9 mantle-melt / degassing calculations | not on PyPI — install from GitHub (below) |

**`moviepy` version.** The animation notebooks require **moviepy ≥ 2.0** (they use
the class-based effects API, e.g. `from moviepy.video.fx.FadeIn import FadeIn`).
Install/upgrade with `pip install -U "moviepy>=2.0"`, then restart the kernel.
Two API changes between 1.x and 2.x show up here and are already fixed in the
`utils/0D-Carbon-Panel-Plot-Videos*.ipynb` notebooks:

- `No module named 'moviepy.video.fx.FadeIn'` → the effects API is class-based in
  2.x (1.x used lowercase function modules).
- `No module named 'moviepy.editor'` → `moviepy.editor` was removed in 2.x; use
  `import moviepy as mpy` (all clip classes are now top-level).

moviepy renders via FFmpeg — recent versions bundle it through `imageio-ffmpeg`,
but if export complains FFmpeg is missing, add it with
`conda install -c conda-forge ffmpeg`.

**`melt` (Ben Mather).** The step-9 notebooks import the `melt` package from Ben
Mather's GitHub (`brmather`). It is not on PyPI/conda, so install it from source
into the same environment, e.g.:

```bash
pip install git+https://github.com/brmather/melt.git
# or clone and install editable:
git clone https://github.com/brmather/melt.git && pip install -e melt
```

(Confirm the exact repository URL — it is a research repo, not a package index
entry. If import still fails, add the cloned repo to `PYTHONPATH`.)

If a step-9 notebook aborts with `Missing optional dependency 'pytables'` or
`No module named 'tables'`, install `pytables` and resume with
`CO2_START_FROM=<notebook>` (above) rather than re-running the whole chain.

## Pipeline

| Step | File | Role |
|------|------|------|
| 0 | `00_update_hybrid_ccd.py` | Re-runs the clean CCD workflow (`../run_all.py`, steps 0–6) and writes the **area-weighted** hybrid CCD into step 7 input (`CCD_sl_hybrid_2026.txt`) and `Paper/Figures/`. Guarantees the current CCD flows through everything. |
| 7 | `01_carbonate_thickness_min_mean_max.py` | De-notebooked `carbonate_sediment_thickness_2026.ipynb`; runs it 3× (min/mean/max = `sed_rate_min/best/max.txt`; CCD fixed at the hybrid) → three grid sets. |
| 8 | `02_carbonate_volume_stats.py` | DM2026-only stats: global carbonate volume & area through time (mean + min–max envelope). Plus the well-data ground-truth (`08_...py`). |
| 9 | `03_run_co2_notebooks.sh` | Headless-executes CO₂ notebooks 01–06 in place (`jupyter nbconvert --execute` / papermill), then `standalone_plots/run_all.py`. |
| — | `run_all.sh` | Orchestrates 0→7→8→9. |
| — | `config.sh` | Paths (incl. `STEPS_ROOT`), notebook order, toggles. |

## Where the step 7/8/9 code lives (`STEPS_ROOT`)

The step-7/8/9 code and grids have been consolidated into this repo. `STEPS_ROOT`
(in `config.sh`) now defaults to **`$CLEAN_REPO/steps_carbon`** — everything
(`step7_carbonate`, `step8_analysis`, `step9_co2/Alfonso_etal_2024_DM26`) resolves
inside the repo, so the old `../scripts/` tree is no longer referenced and can be
retired. Override `STEPS_ROOT` only if you relocate the heavy grids elsewhere.

## min / mean / max convention

Scenarios are the **maximum carbonate sedimentation-rate curve**
(`sed_rate_min/best/max.txt`), NOT the CCD; the CCD is fixed at the area-weighted
DM2026 hybrid for all three. (Model choice DM2026-vs-BW1991 would be the CCD
curve — BW1991 is dropped.)

## Data availability (large grids)

The input grids (Alfonso 2024 seafloor-age, spreading-rate, and paleobathymetry
grids; the plate model; and the reservoir/subducted-carbon grids under
`Alfonso_etal_2024_DM26/`) are too large for version control. They are **git-
ignored** (see the repo `.gitignore`) and will be archived on **Zenodo**; download
them and place them at the paths in `config.sh` / the step scripts before running.
Small text inputs (CCD curves, sed-rate curves, reference CCDs, the well table)
are tracked in the repo.

## Confirm before first run

- `STEPS_ROOT` (and, after consolidation, the in-repo grid locations).
- `jupyter`/`papermill` + a `python3` kernel available.
- The Alfonso2024 plate model downloads via `plate-model-manager` (needs network)
  unless you set `USE_LOCAL_MODEL = True` in `01_...py`.
- CO₂ notebook order and any prep notebooks (`CO2_PREP_NOTEBOOKS` in `config.sh`).
