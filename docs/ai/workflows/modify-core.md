---
description: "Route for implementing and verifying a change to an MMA numerical-core module (CUPRAD, CTDSE, Hankel) or to the HDF5 namelists. Use whenever a task edits module code or anything that can change numerical results or the structure of output files."
mma-status: maintained
mma-track: core
mma-owner: maintainers
---

# Workflow: modify-core

## Purpose

Change a core module so that the change is correct, reviewed and traceable. Not for: analyses that the existing interfaces already support ([develop-analysis](develop-analysis.md)).

## Entry conditions

- The change and its reason are stated: bug fix, feature, optimisation or refactoring.
- The affected module is in the core block of `.github/CODEOWNERS`.

## Route

1. State the module, the change and the expected numerical-impact class (`none`, `roundoff` or `intended`; see [CONTRIBUTING.md](../../../CONTRIBUTING.md)).
2. Work on a branch. Keep refactoring, numerical changes, reference-data updates and documentation in separate commits.
3. Implement. If the HDF5 layout changes, inspect the affected producers and consumers; update the relevant namelists, interface documentation and [README.md, Inputs](../../../README.md#inputs) as needed.
4. Verify with the existing tests ([README.md, Running unit tests](../../../README.md#running-unit-tests)) and examples for the module. Until reduced reference cases exist, state exactly which checks ran and which could not.
5. If the change departs from the model or approach described in the CPC paper, add an entry to `docs/deviations-from-paper.md`. Document the detailed analysis in the module's README or the source, not in AI-support files.
6. Check the documentation with the skill `mma-documentation-consistency`.
7. Prepare a reviewable diff and a filled-in pull-request description. Create commits with the trailers from [CONTRIBUTING.md](../../../CONTRIBUTING.md), push or open a pull request only within the user-authorised scope.

## Decision points and transitions

- Environment problems: [deploy-hpc](deploy-hpc.md) and the skill `mma-diagnose-hpc`.
- The change turns out to be unnecessary because existing interfaces suffice: [develop-analysis](develop-analysis.md).
- A test fails after the change: fix the code, or explain the change of results. Never alter reference data or tolerances to make a test pass, except in a separate, justified commit.

## Exit criteria

- The reviewable diff and pull-request description state the impact, the evidence, what was not run and the AI assistance. Any created commits carry the required trailers.
- Namelists, README and module documentation agree.
- The deviations log is updated when the model departs from the paper.
