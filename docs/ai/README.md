# MMA AI support: handbook

This handbook defines how MMA supports AI agents: which kinds of files exist, where they live, how they are named and checked, and how anyone can extend them. It is the authoritative source for these conventions. Agents enter through [AGENTS.md](../../AGENTS.md); humans start here.

Conventions version: 2 (see the [changelog](#14-changelog)).

## 1. Principles

- **Provider independence.** Everything is plain Markdown in two open formats: [AGENTS.md](https://agents.md) for instructions and [Agent Skills](https://agentskills.io/specification) for skills. Agent-specific files are generated signposts that point to these.
- **One authoritative source per fact.** Other files summarise in at most a sentence and link to it; the [ownership map](#12-ownership-map) says which file owns which fact.
- **Layers of knowledge.** The CPC paper introduces MMA's core philosophy and physical model, as of version v1.1.0. Decisive analyses and details live in the lower layers, the module READMEs and the source, which are authoritative for them.
- **Links run one way.** AI-support files link into authoritative documentation, never the reverse. The one exception is a short section in `README.md` that introduces the AI support: it says that autonomy is the user's choice, that AI-assisted work is traced and that anyone can extend the AI support, and it points humans here, agents to `AGENTS.md` and Docker users to `docker/agents/README.md`.
- **Autonomy is the user's choice.** Traceability is a property of commits and pull requests, defined in [CONTRIBUTING.md](../../CONTRIBUTING.md), and holds at every level of autonomy.
- **Open by design.** The shipped entities are templates and intended uses, not a closed set. Users and developers are invited to add entities of every kind, on their own responsibility.
- **Deterministic checks gate; judgement advises.** CI runs only reproducible checks. AI-assisted reviews produce findings for humans.

## 2. Entity model

| Entity | Purpose | Lives in | Contains | Must not contain |
| --- | --- | --- | --- | --- |
| Router | The one file every agent loads | `AGENTS.md` | Map, routing, repository rules, autonomy defaults | Technical content owned elsewhere |
| Track | Entry point for a broad use case | A heading under Routing in `AGENTS.md` | Routes to workflows and skills | Procedures or facts |
| Workflow | Route to an objective | `docs/ai/workflows/<name>.md` | Entry conditions, route, decisions, exit criteria, links | Commands or technical facts |
| Skill | Reusable procedure | `.agents/skills/<name>/SKILL.md` | Steps, commands, interpretation of results, pitfalls | Facts owned elsewhere |
| Knowledge | What is true about MMA | `README.md`, module READMEs, namelists, source, the CPC paper | Physics, numerics, interfaces, formats, build facts | Agent-specific instructions |
| Record | Dated observation, such as experience with a machine | `docs/ai/records/` | Machine, date, commit, status, findings | Copies of module lists or install steps; personal or sensitive data |
| Module contract | Interface of a core module | A "Module contract" section in the module's README | Inputs, outputs, invariants, tests, owners | Implementation detail |
| Signpost | Lets a given agent find an entity | Each agent's standard paths, generated | An import line, a link or a one-line pointer | Any content of its own |

Glossary:

- **Core module**: a module with defined inputs that performs calculations and provides outputs; today CUPRAD, CTDSE and Hankel, listed in the core block of `.github/CODEOWNERS`.
- **Numerical-impact class**: `none`, `roundoff` or `intended`, as defined in [CONTRIBUTING.md](../../CONTRIBUTING.md).
- **Consequential action**: an action that needs explicit authorisation unless granted in advance; the list is in [AGENTS.md](../../AGENTS.md).
- **Tier**: the level of support for an agent; see [section 10](#10-signposts-and-supported-agents).

## 3. Where things live

```text
AGENTS.md                        router (hand-written)
CONTRIBUTING.md                  commit and pull-request conventions
.agents/skills/<name>/           skills; workflow signpost skills are generated here too
docs/ai/README.md                this handbook
docs/ai/workflows/               workflows
docs/ai/templates/               one template per entity type
docs/ai/records/                 records
docs/ai/signposts.yml            signpost manifest
docs/ai/local/                   personal entities, ignored by Git
docs/deviations-from-paper.md    departures from the published version
docker/agents/                   optional AI agents in the Docker image: use and recipes
ai/ai_signposts.py               generates and verifies signposts
ai/check_ai_layer.py             structural checks
ai/tests/                        tests of the AI support, such as the offline installer tests
```

Generated signposts live at each agent's standard paths, listed in [section 10](#10-signposts-and-supported-agents).

## 4. Naming

- Use lowercase letters, digits and single hyphens. A skill's name equals its folder name, as the Agent Skills specification requires; a workflow's name equals its file name.
- Every repository skill carries the `mma-` prefix, and personal skills must not use it. In Claude Code a personal skill overrides a project skill of the same name, and a project skill overrides a bundled one; the prefix prevents both silent replacements.
- The prefix is not exclusive to skills: other files of the package may use it too, such as the container commands `mma-build` and `mma-jupyter`.
- Skill names are unique. Workflows are named after their objective (`modify-core`) and their signpost skills add the prefix (`mma-modify-core`), so no hand-written skill may take that name.
- Name entities after what they do, not after a tool or a person.

## 5. Metadata and status

Workflows start with YAML front matter; skills put the same keys in the specification's `metadata` map. Values are strings; quote them when they contain a colon.

| Key | Meaning | Values |
| --- | --- | --- |
| `description` | When to use the entity. For workflows it also becomes the trigger of their signpost skill; skills use the standard `description` field. | One or two sentences, at most 1024 characters |
| `mma-status` | Support level | `maintained`, `community`, `experimental`, `generated` |
| `mma-track` | Track or tracks served | `operation`, `applications`, `core`, `documentation`, comma-separated |
| `mma-owner` | Who maintains it | A GitHub handle, a team or `maintainers` |

Records use their own keys, including `record-status` (`unverified` or `verified`); see `docs/ai/templates/record-hpc.md`.

## 6. Scopes and precedence

| Scope | Where | Status label | What the project guarantees |
| --- | --- | --- | --- |
| Maintained | Repository, owned by a maintainer | `maintained` | Structure checks and a named owner who keeps it current |
| Contributed | Repository, merged through a pull request | `community` or `experimental` | Structure checks only; listed apart in `AGENTS.md` |
| Personal | The agent's user directory, such as `~/.agents/skills/`, or `docs/ai/local/` | None | Nothing; the user's own responsibility |

Precedence between scopes differs from agent to agent. Do not rely on it: use distinct names, and keep the `mma-` prefix for repository skills.

## 7. Templates

| Entity | Template |
| --- | --- |
| Track | `docs/ai/templates/track.md` |
| Workflow | `docs/ai/templates/workflow.md` |
| Skill | `docs/ai/templates/skill.md`, copied to `.agents/skills/<name>/SKILL.md` |
| Record | `docs/ai/templates/record-hpc.md` |
| Module contract | `docs/ai/templates/module-contract.md` |
| Signpost entry | `docs/ai/templates/signpost-entry.yml` |

## 8. Adding or changing an entity

1. Pick the entity type from the table in [section 2](#2-entity-model). If none fits, see [section 9](#9-proposing-a-new-entity-type).
2. Copy the template, choose a name ([section 4](#4-naming)) and fill in the metadata ([section 5](#5-metadata-and-status)).
3. Link to authoritative documents instead of restating them.
4. For a workflow or skill meant for others, add one line to the Routing section of `AGENTS.md`; contributed entities go under "Contributed entities".
5. Run `python3 ai/ai_signposts.py`, then `python3 ai/check_ai_layer.py`, and fix every error.
6. Open a pull request following [CONTRIBUTING.md](../../CONTRIBUTING.md), or keep the entity personal.

## 9. Proposing a new entity type

Open a pull request that adds all of the following, so that the new type is as well defined as the existing ones:

1. a row in the entity model ([section 2](#2-entity-model));
2. a template in `docs/ai/templates/`;
3. its location and naming rule;
4. how agents discover it, including any new signpost kind in `ai/ai_signposts.py`;
5. its structural check in `ai/check_ai_layer.py`;
6. an entry in the [changelog](#14-changelog).

## 10. Signposts and supported agents

Signposts let each agent find MMA's entities at the paths it already knows. They are generated from `docs/ai/signposts.yml` by `ai/ai_signposts.py` and contain no content of their own: an import line, a link or a one-line pointer. Never edit them by hand; CI fails when they drift from the manifest.

Tier 1 agents get generated signposts from the project. Consistent files do not prove that an agent loads them; the verification table below records what was actually observed.

| Entity | Canonical | Claude Code | Codex | GitHub Copilot | Gemini CLI | Cursor |
| --- | --- | --- | --- | --- | --- | --- |
| Router | `AGENTS.md` | `CLAUDE.md` imports it | Native | Native | `.gemini/settings.json` lists it | Native |
| Skills | `.agents/skills/` | Pointer skills in `.claude/skills/` | Native | Native | Native | Native |
| Workflows | `docs/ai/workflows/` | Signpost skills | Signpost skills | Signpost skills | Signpost skills | Signpost skills |
| Core-module rule | `AGENTS.md` | `.claude/rules/mma-numerical-core.md` | Nested `AGENTS.md` in each module | `.github/instructions/mma-numerical-core.instructions.md` | Nested `AGENTS.md` in each module | `.cursor/rules/mma-numerical-core.mdc` |

For Claude Code, the generator also writes `.claude/.gitignore`, so that the files that Claude Code may create there for a single user, its local settings and worktrees, stay out of Git.

Verified discovery, recorded after testing an agent against this repository:

| Agent | Version and date | Router loaded | Skills listed | Core-module rule |
| --- | --- | --- | --- | --- |
| Claude Code | not yet tested | | | |
| Codex | not yet tested | | | |
| GitHub Copilot | not yet tested | | | |
| Gemini CLI | not yet tested | | | |
| Cursor | not yet tested | | | |

Tier 2 is open to everyone: add an entry of kind `pointer-file` to the manifest (template: `docs/ai/templates/signpost-entry.yml`), mark it `community` and say how you tested it. The project does not test Tier 2 signposts.

Signposts only let an installed agent find MMA's entities; they do not install it. Installing agents in the Docker image is separate: one recipe per agent in `docker/agents/`, with the contract in `docker/agents/README.md`. A new agent may need a signpost, a recipe or both.

Known caveats:

- GitHub Copilot reads both `.claude/skills/` and `.agents/skills/` and may list a skill twice.
- Claude Code reads `CLAUDE.md` instead of `AGENTS.md` when both exist, so the generated `CLAUDE.md` holds only the import line. GitHub Copilot also reads `CLAUDE.md`; the import adds no content of its own.
- Codex builds its initial instruction chain from the repository root to its working directory. Keep cross-module rules in the root `AGENTS.md`; nested files provide local reminders. Launch Codex in the repository root for work spanning modules.

## 11. Verification

- **Structural checks.** `python3 ai/check_ai_layer.py` checks skills, workflows, `AGENTS.md`, links and anchors, skill references, signposts, the core block and records. CI runs it on every pull request. With `--with-readme` it also checks the links, anchors and skill references in `README.md`, reporting findings there as warnings.
- **Semantic checks.** The skill `mma-documentation-consistency` and its verification summary, `.agents/skills/mma-documentation-consistency/references/verification-summary.md`. Their findings are advisory.

## 12. Ownership map

| Knowledge | Authoritative source |
| --- | --- |
| Core philosophy, physical model and modelling assumptions | The CPC paper on version v1.1.0 ([doi:10.1016/j.cpc.2026.110417](https://doi.org/10.1016/j.cpc.2026.110417), open access) |
| Departures from the paper | `docs/deviations-from-paper.md` |
| Module concepts, numerics and decisive analyses | Module READMEs |
| Module interfaces | The "Module contract" section of each core module's README |
| Input parameters and constraints | [README.md, Inputs](../../README.md#inputs) |
| HDF5 layout | The namelists `shared_python/MMA_administration.py`, `1DTDSE/sources/h5namelist.h` and `CUPRAD/sources/global_variables.f90` |
| Installation and environment | [README.md, Installations](../../README.md#installations) and the `Dockerfile` |
| AI agents in the Docker image | `docker/agents/README.md`: use, recipes and their contract |
| Execution order | [README.md, Execution pipeline](../../README.md#execution-pipeline) |
| Machine configurations | `Modules/load_modules.sh` |
| Experience with machines | Records in `docs/ai/records/` |
| Commit and pull-request conventions | [CONTRIBUTING.md](../../CONTRIBUTING.md) |
| Core modules and their owners | The core block of `.github/CODEOWNERS` |
| Autonomy defaults and repository rules for agents | [AGENTS.md](../../AGENTS.md) |
| AI-support conventions | This handbook |
| Agent entry points | Signposts generated from `docs/ai/signposts.yml` |

## 13. Planned entities

These entities are planned but do not exist yet; their names are reserved.

| Name | Kind | Waits for |
| --- | --- | --- |
| `mma-numerical-validation` | Skill | Reduced reference cases and a comparison script |
| `mma-extend-ai` | Skill | Templates that have been used in practice |
| `mma-doxygen-compliance` | Skill | Doxygen conventions and a pilot module |
| `environment/mma-env.sh` | Tool | A wrapper that re-creates the environment in non-interactive shells |

## 14. Changelog

- **Version 2** (2026-11-27):
  - The AI-support scripts and their tests moved from the former tools directory to `ai/`. The planned environment wrapper is reserved in `environment/` ([section 13](#13-planned-entities)).
- **Version 1** (2026-11-15): initial conventions.
