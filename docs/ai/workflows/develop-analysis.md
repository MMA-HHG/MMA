---
description: "Route for developing a scientific analysis, parameter study or application on top of MMA's existing interfaces (HDF5 archives, shared Python code, the interactive CTDSE library) without changing core modules. Use when the task is post-processing, plotting, parameter scans or new scripts that consume MMA results."
mma-status: maintained
mma-track: applications
mma-owner: maintainers
---

# Workflow: develop-analysis

## Purpose

Build analysis code that uses MMA as it is; core modules stay unchanged. Not for: running new simulations ([run-simulation](run-simulation.md)) or extending a module ([modify-core](modify-core.md)).

## Entry conditions

- The scientific question and the quantities needed are stated.
- Existing outputs are available, or simulations can be run with [run-simulation](run-simulation.md).

## Route

1. Find the interfaces: the HDF5 layout through the skill `mma-io`, the shared Python code in `shared_python/`, and the examples in `jupyter_examples/`.
2. Take units and conventions from the authoritative sources in the ownership map of `docs/ai/README.md`; never assume them.
3. Write the analysis outside the core-module directories.
4. Validate it on a small or known case before scaling up.
5. Document how to reproduce it: MMA commit, inputs and scripts.

## Decision points and transitions

- An interface is missing, such as a quantity that is not written or an option that does not exist: stop, state the transition to core development and the reason, then follow [modify-core](modify-core.md).
- New simulations are needed: [run-simulation](run-simulation.md).

## Exit criteria

- The analysis runs from documented inputs and states its units.
- No core module was edited.
- Limitations and unverified assumptions are listed.
