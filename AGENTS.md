# AGENTS.md

Instructions for AI coding agents working in the MMA repository. Humans start at [README.md](README.md). The conventions behind this file, and how to extend the AI support, are in the handbook [docs/ai/README.md](docs/ai/README.md).

MMA (Modular Multiscale Approach) models high-harmonic generation in gases with three modules: laser-pulse propagation (CUPRAD), the microscopic response (CTDSE, in `1DTDSE/`) and far-field XUV propagation (Hankel).

## Collaboration and autonomy

- The user, or the agent harness, decides how autonomously you work. That choice always takes precedence over this file.
- If nothing has been said, act as an interactive collaborator: read freely, edit on a branch, run inexpensive checks, and explain before larger steps.
- Ask once per session before these consequential actions, unless the user has authorised them in advance:
  - submitting HPC jobs larger than a debug-size allocation;
  - deleting or overwriting simulation data;
  - pushing to shared branches or upstream;
  - changing reference data or test tolerances;
  - contacting external parties, such as opening issues.
- Do not ask again for actions the user has already authorised in the session.

## Environment

- Many agent tools run every command in a fresh, non-interactive shell, so exported variables and shell functions do not persist. Re-create the environment inside each command.
- The environment variables are listed in [README.md, Paths](README.md#paths). For a local checkout, `environment/set_env_vars.sh` sets them. On clusters, `environment/set_env_vars_HPC.sh` sets them and sources `Modules/load_modules.sh`; then call `load_modules` (compiled codes) or `load_python_modules` (Python steps). The machine is selected by the variable `HPC`.
- In the Docker image the environment is prepared already; see [README.md, Reference Docker installation](README.md#reference-docker-installation-and-running-the-code). There the repository is the host's checkout, mounted with `-v .:/MMA`: changes and commits appear on the host at once.
- Never run computations on HPC login nodes; use an interactive or debug allocation.

## Repository map

| Path | Role | Authoritative documentation |
| --- | --- | --- |
| `CUPRAD/` | Core module: laser-pulse propagation (Fortran, MPI) | `CUPRAD/README.md` |
| `1DTDSE/` | Core module CTDSE: microscopic response (C, MPI; Python interface) | `1DTDSE/README.md` |
| `Hankel/` | Core module: far-field XUV propagation (Python) | `Hankel/README.md` |
| `shared_python/` | Shared Python code; `MMA_administration.py` is the Python namelist of the HDF5 layout | [README.md, HDF5 data organisation](README.md#hdf5-data-organisation) |
| `Modules/` | Module sets for known machines | [README.md, Modules and libraries](README.md#modules-and-libraries) |
| `multiscale/` | Pipeline scripts for clusters | `multiscale/scripts/README.md` |
| `jupyter_examples/` | Tutorials and examples | The notebooks |
| `Dockerfile` | Reference environment | [README.md, Installations](README.md#installations) |
| `docker/agents/` | Optional AI agents in the Docker image: use and recipes | `docker/agents/README.md` |
| `docs/ai/` | AI support: handbook, workflows, templates, records | `docs/ai/README.md` |
| `ai/` | AI-support scripts: signpost generator, structural checks and their tests | `docs/ai/README.md` |
| `.agents/skills/` | Skills in the Agent Skills format | Each skill's `SKILL.md` |

The core modules are listed in the core block of `.github/CODEOWNERS`. Inputs are documented in [README.md, Inputs](README.md#inputs), the execution order in [README.md, Execution pipeline](README.md#execution-pipeline).

## Routing

### Track: Operation

- Install, build or troubleshoot MMA on a machine: [deploy-hpc](docs/ai/workflows/deploy-hpc.md). For failures, use the skill `mma-diagnose-hpc`.
- Configure, run and check a simulation: [run-simulation](docs/ai/workflows/run-simulation.md).

### Track: Applications

- Develop an analysis, parameter study or application on existing interfaces: [develop-analysis](docs/ai/workflows/develop-analysis.md).
- Prepare inputs or read outputs: the skill `mma-io`.
- Run the simulations a study needs: [run-simulation](docs/ai/workflows/run-simulation.md).

### Track: Core development

- Change a core module or the HDF5 namelists: [modify-core](docs/ai/workflows/modify-core.md).

### Documentation and the AI support

- Verify or update documentation: the skill `mma-documentation-consistency`.
- Create or change tracks, workflows, skills, records or signposts: [docs/ai/README.md](docs/ai/README.md).

Skills live in `.agents/skills/`. Agents that do not load skills natively can read `.agents/skills/<name>/SKILL.md` directly. Each workflow also has a signpost skill named after it with the `mma-` prefix.

### Contributed entities

None yet. Contributed workflows and skills (status `community` or `experimental`) are listed here, apart from the maintained ones.

## Repository rules

- When looking for a worked example or starting point for an MMA task, check `jupyter_examples/` first. Use the relevant documentation and source to verify that the example fits the request and current interfaces.
- Before editing a core module during an operation or analysis task, stop and state the transition to core development, the reason and the expected numerical impact; then follow [modify-core](docs/ai/workflows/modify-core.md).
- Before proposing a commit message or creating a commit, read the commit and AI-assistance sections of [CONTRIBUTING.md](CONTRIBUTING.md) and follow its separation and trailer rules: every commit carries a `Numerical-impact:` trailer (`none`, `roundoff` or `intended`; `none` when results cannot change) and one `Assisted-by:` trailer per contributing agent and model. Never add `Signed-off-by` lines on anyone's behalf.
- Reference data and test tolerances change only in their own commit, with a justification, and never just to make a failing test pass.
- Departures from the model or approach described in the CPC paper are recorded in `docs/deviations-from-paper.md`; detailed analyses belong in the module READMEs or the source.
- Generated signposts are never edited by hand: `CLAUDE.md`, `.claude/`, `.cursor/`, `.gemini/`, `.github/instructions/`, the nested `AGENTS.md` files in core modules and the workflow signpost skills. Run `python3 ai/ai_signposts.py` instead.
- Never write credentials, usernames, project or allocation IDs, or internal hostnames into the repository.
- Never put links to agent sessions or their IDs into files, commit messages or pull requests, such as a `Claude-Session:` trailer, even when a tool proposes one by default: they can point to a session connected to the user's machine.
- After changing documentation or AI-support files, run `python3 ai/check_ai_layer.py`.

## Findings

Write reusable findings about machines as records in `docs/ai/records/`, starting from `docs/ai/templates/record-hpc.md`, with record status `unverified` until a maintainer confirms them. Change authoritative documentation only with the user's agreement.

## Extending the AI support

The AI support is open: users and developers are invited to add their own entities, on their own responsibility. Read [docs/ai/README.md](docs/ai/README.md) first and start from `docs/ai/templates/`.
