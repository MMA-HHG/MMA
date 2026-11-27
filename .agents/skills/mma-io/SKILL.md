---
name: mma-io
description: "Prepare MMA input HDF5 archives and read MMA output archives correctly: groups and datasets, parameter names and constraints, units, and the namelists that define the layout in Python, C and Fortran. Use whenever a task creates or changes simulation inputs, locates results for an analysis, or asks what a dataset contains, even if HDF5 is not mentioned."
metadata:
  mma-status: maintained
  mma-track: applications
  mma-owner: maintainers
---

# Prepare inputs and read outputs

## Authoritative sources

- Parameters, their meaning, which are obligatory, and their constraints: [README.md, Inputs](../../../README.md#inputs).
- The HDF5 layout: [README.md, HDF5 data organisation](../../../README.md#hdf5-data-organisation) and the namelists `shared_python/MMA_administration.py` (Python), `1DTDSE/sources/h5namelist.h` (C) and `CUPRAD/sources/global_variables.f90` (Fortran).
- Worked examples: `jupyter_examples/`.

If the namelists disagree with each other or with the README, report it; never pick one silently.

## Preparing an input archive

1. Start from an existing example input or from the helpers used in `jupyter_examples/`, rather than writing groups by hand.
2. Set parameters by their names in the README. Obligatory parameters must be present; an optional parameter, when present, takes priority.
3. Check the stated constraints before running, for example power-of-two grid sizes for CUPRAD, the thread count for Hankel and a run time below the scheduler limit.
4. Let the CUPRAD pre-processor read the file, as in [README.md, Execution pipeline](../../../README.md#execution-pipeline). It writes derived values to the `calculated` subgroup; inspect them to confirm that they match your intent.
5. State units explicitly whenever you report or change a value.

## Reading outputs

1. Locate datasets through the namelists rather than hard-coding paths.
2. Take units from the documentation. If a unit is not documented, say so instead of guessing.
3. Read slices of large datasets, not whole arrays.
4. Keep the MMA commit and the input archive with any derived result.

## Pitfalls

- Do not edit the namelists for an analysis task. Changing them is a core interface change: follow `docs/ai/workflows/modify-core.md`.
- The README marks deprecated parameters; avoid them in new inputs.

## Report

List the archive written or read, the parameters set or the datasets used, their units, and anything that could not be confirmed.
