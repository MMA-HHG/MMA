---
description: "Route for preparing, installing and diagnosing MMA on a computing system, from the reference Docker image to an HPC cluster. Use when the task is to build MMA, set up its environment or make it run on a particular machine."
mma-status: maintained
mma-track: operation
mma-owner: maintainers
---

# Workflow: deploy-hpc

## Purpose

Obtain a working MMA build and environment on a given machine, and leave a dated record of what worked. Not for: designing simulations ([run-simulation](run-simulation.md)) or changing module code ([modify-core](modify-core.md)).

## Entry conditions

- The target is known: the reference Docker image, a workstation or a named cluster.
- On a cluster: access to a login node, and the user has said which allocation may be used.

## Route

1. Look for a record of the machine in `docs/ai/records/` and for its module set in `Modules/load_modules.sh`, selected by `HPC`. A `verified` record is the starting point; an `unverified` one is a hint.
2. Without machine-specific knowledge, start from the reference environment: the `Dockerfile` and [README.md, Installations](../../../README.md#installations).
3. Build as described in [README.md, Custom installations](../../../README.md#custom-installations). Re-create the environment inside every command, as `AGENTS.md` explains.
4. On any failure, follow the skill `mma-diagnose-hpc` (`.agents/skills/mma-diagnose-hpc/SKILL.md`).
5. Run a small case through the pipeline ([README.md, Execution pipeline](../../../README.md#execution-pipeline)). On a cluster, use an interactive or debug allocation, never the login node.
6. Write or update the record from `docs/ai/templates/record-hpc.md`, with record status `unverified` unless a maintainer confirms it.

## Decision points and transitions

- The failure lies in module code rather than in the environment: [modify-core](modify-core.md).
- The build works and a study should follow: [run-simulation](run-simulation.md).
- A job larger than a debug-size allocation needs the user's authorisation (`AGENTS.md`, Collaboration and autonomy).

## Exit criteria

- MMA builds, a small case completes, and the commands used are reproducible in a fresh shell.
- A record exists or was updated, with date, MMA commit, module set and record status.
- Everything not verified is stated.
