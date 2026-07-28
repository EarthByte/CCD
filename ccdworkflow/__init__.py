"""
ccdworkflow
===========

Shared configuration and I/O helpers for the Dutkiewicz & Müller CCD–sea-level
workflow (steps 0–6). Import :mod:`ccdworkflow.config` for canonical paths and
parameters, and :mod:`ccdworkflow.io` for the standard text-table readers and
writers used throughout the pipeline.

Design principles for this clean re-implementation
--------------------------------------------------
* Every step reads its inputs from the *declared output location* of the
  previous step (see ``config.py``). No file is ever hand-copied between steps.
* Every tabular (x-y / x-y-z) product is written as a single ``.txt`` file with
  a leading ``# ``-commented, tab-separated header. No parallel ``.tsv`` / ``.xy``
  / ``.csv`` duplicates of the same data.
* Numerical parameters that were previously tuned interactively (ipywidgets) are
  fixed here as named constants in ``config.Params`` and documented.
"""

from . import config, io  # noqa: F401

__all__ = ["config", "io"]
