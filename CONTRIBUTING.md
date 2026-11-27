# Contributing to MMA

Contributions are welcome through pull requests. This file defines the conventions that keep MMA's history traceable. They apply to human and AI-assisted work alike, whatever degree of autonomy a contributor gives an AI agent.

## Core modules

A core module is any module with defined inputs that performs calculations and provides outputs. Today these are CUPRAD (`CUPRAD/`), CTDSE (`1DTDSE/`) and Hankel (`Hankel/`). The shared HDF5 namelists define their inputs and outputs, so changing them is a core interface change.

The authoritative list is the core block in `.github/CODEOWNERS`. A new module joins it when it publishes a module contract, written from `docs/ai/templates/module-contract.md`.

## Commits

- Keep refactoring, numerical changes, reference-data updates and documentation in separate commits.
- End every commit that changes core-module code with a numerical-impact trailer:
  - `Numerical-impact: none`: results are identical;
  - `Numerical-impact: roundoff`: results differ only at floating-point level;
  - `Numerical-impact: intended`: results change on purpose; the pull request explains why and shows before-and-after evidence.
- Shared build files and helpers outside the modules carry the same trailer whenever they can change results. Commits that change only documentation need no trailer, unless they are AI-assisted (see below).
- Change reference data and test tolerances only in their own commit, with a justification in the message.

## AI-assisted contributions

AI assistance is welcome. It must be traced:

- End every commit that contains AI-generated or AI-modified content with one trailer per agent and model: `Assisted-by: <agent>:<model>`, for example `Assisted-by: claude.ai:claude-opus-5-5`. Trailers that a tool adds by itself may stay alongside.
- Every AI-assisted commit also carries a `Numerical-impact:` trailer, even when it changes no code; `none` then states explicitly that results cannot change.
- In the pull request, state the degree of assistance (none, assisted, substantially generated), the agents and models, a summary of the prompts or session, the affected parts and how they were checked.
- A named human author is responsible for every change. AI tools never add `Signed-off-by` lines.

## Pull requests

- Use the pull-request template; it asks for the numerical impact, the evidence and the AI assistance.
- Changes to core modules are reviewed by their owners, as listed in `.github/CODEOWNERS`.
- Departures from the model or approach described in the published CPC paper are recorded in `docs/deviations-from-paper.md`; detailed analyses go into the module READMEs or the source.
- Before opening a pull request that touches documentation or AI-support files, run `python3 ai/check_ai_layer.py`.
