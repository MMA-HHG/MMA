---
description: "Route for configuring, executing and checking an MMA simulation with the existing pipeline: CUPRAD pre-processor and propagation, TDSE preparation, TDSE, merge and Hankel transform. Use when the task is to set up inputs, run the pipeline locally or on a cluster, or check its outputs."
mma-status: maintained
mma-track: operation, applications
mma-owner: maintainers
---

# Workflow: run-simulation

## Purpose

Produce a simulation result with the existing pipeline and check that it is sound; module code stays unchanged. Not for: installation problems ([deploy-hpc](deploy-hpc.md)) or new analyses of existing results ([develop-analysis](develop-analysis.md)).

## Entry conditions

- A working build on the target machine; otherwise start with [deploy-hpc](deploy-hpc.md).
- The physical setup is stated: gas, laser, geometry and the quantities of interest.

## Route

1. Prepare the input HDF5 archive with the skill `mma-io` (`.agents/skills/mma-io/SKILL.md`).
2. Before running, check the constraints in [README.md, Inputs](../../../README.md#inputs), in particular on MPI ranks, grid sizes, thread counts and run time.
3. Run the steps in the order of [README.md, Execution pipeline](../../../README.md#execution-pipeline). On clusters, follow `multiscale/scripts/README.md`.
4. On a cluster, start with a reduced case in a debug allocation. Larger submissions need the user's authorisation.
5. Check the outputs: the expected groups exist (see `mma-io`), every step finished within its time limit, and the results are physically plausible. Report what was checked and what was not.
6. Keep the MMA commit, the input archive and a description of the environment together with the results.

## Decision points and transitions

- Environment or build problems: [deploy-hpc](deploy-hpc.md).
- Post-processing or new analyses: [develop-analysis](develop-analysis.md).
- The study needs a quantity MMA does not produce, or any code change: stop, state the transition, then follow [modify-core](modify-core.md).

## Exit criteria

- Outputs exist and pass the stated checks.
- Inputs, commit and environment are recorded with the results.
- Open doubts are listed.
