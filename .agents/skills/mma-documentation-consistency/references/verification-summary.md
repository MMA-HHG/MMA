# Verification summary

Use this summary to verify MMA's documentation and AI-support files. It condenses the conventions of the handbook, `docs/ai/README.md`, into checks. If the two disagree, the handbook wins and this summary must be corrected.

## A. The design in brief

1. **One router.** `AGENTS.md` is the only file every agent loads. The tracks (Operation, Applications, Core development) are headings in it. It stays under about 150 lines.
2. **Entities.** Workflows are routes, skills are procedures, knowledge is authoritative human documentation, records are dated observations, module contracts describe the interfaces of core modules, and signposts are generated pointers. Definitions: handbook, section 2.
3. **One source per fact.** Every fact has one authoritative home, given by the ownership map (handbook, section 12). Other files summarise in at most a sentence and link to it.
4. **Links run one way**, from AI-support files into authoritative documents. The only exception is the README section on the AI support, which points to `AGENTS.md`, the handbook and `docker/agents/README.md`.
5. **Numerical core.** Every module with defined inputs, calculations and outputs: today CUPRAD, CTDSE (`1DTDSE/`) and Hankel, plus the shared namelists that define their inputs and outputs. The core block of `.github/CODEOWNERS` is the authoritative list.
6. **Autonomy and traceability.** The user chooses how autonomously an agent works. Traceability is required at every level of autonomy: the `Numerical-impact:` and `Assisted-by:` trailers and the pull-request fields defined in `CONTRIBUTING.md`.
7. **Openness.** Users and developers may add entities of every kind, on their own responsibility. The scopes are maintained, contributed (`community`, `experimental`) and personal. Every repository skill carries the `mma-` prefix and personal skills must not use it; other files of the package, such as container commands, may use the prefix too.
8. **Signposts** are generated from `docs/ai/signposts.yml`, hold no content and are never edited by hand.
9. **The paper.** The CPC paper introduces MMA's core philosophy and physical model (version v1.1.0); decisive analyses belong in the module READMEs and the source. A departure from what the paper describes is justified only when it is recorded in `docs/deviations-from-paper.md`.

## B. Structural checks

Run `python3 tools/check_ai_layer.py`, adding `--with-readme` when `README.md` changed, and report its output verbatim.

| ID | What it checks |
| --- | --- |
| S1 | Skills: front matter, name equal to the folder, name format, description length, status label, `mma-` prefix, body length, no clash with a workflow's signpost skill |
| S2 | Workflows: front matter keys, required headings, length, routed from `AGENTS.md` |
| S3 | `AGENTS.md`: required sections, length, every maintained workflow and skill routed |
| S4 | Links, anchors and paths: relative Markdown links, their anchors, and relative code-span paths with a directory part exist; absolute paths, such as container paths, are skipped; in `README.md` only links and anchors are checked |
| S5 | Skill references: every `mma-` name introduced as a skill, as in "the skill `mma-io`", is an existing skill, a workflow's signpost skill or a planned entity |
| S6 | Signposts: every generated file exists and matches the manifest; no stale generated files remain |
| S7 | Core block in `.github/CODEOWNERS`: present, non-empty and pointing to existing paths; placeholder owners reported |
| S8 | Records: required keys, record status and staleness |

## C. Semantic checks

For each changed document or statement, answer:

| ID | Question |
| --- | --- |
| M1 | Does it say the same as its authoritative source in the ownership map? |
| M2 | Does it restate facts owned elsewhere beyond a one-sentence summary? If so, replace the text with a link. |
| M3 | Are entity names, paths, trailers, status labels and track names exactly those of the handbook and `CONTRIBUTING.md`? |
| M4 | Do statements about core modules, trailers, reference data and AI assistance agree with `CONTRIBUTING.md` and `.github/CODEOWNERS`? |
| M5 | Does it neither restrict the user's choice of autonomy beyond the consequential actions in `AGENTS.md` nor grant permissions beyond them? |
| M6 | Is every departure from the paper's model or approach backed by an entry in `docs/deviations-from-paper.md`, and is each detailed analysis placed in a module README or the source rather than in an AI-support file? |
| M7 | Is extension described as open, with the responsibility for contributed and personal entities stated? |
| M8 | Does the text address its audience: humans in READMEs; agents in `AGENTS.md`, workflows and skills? |
| M9 | Do statements about agents and their file names match `docs/ai/signposts.yml` and the handbook, section 10, and does any claim that an agent finds MMA's files rest on the verification table there? |
| M10 | If the README section on the AI support changed: is it short, does it say that autonomy is the user's choice, that AI-assisted work is traced and that anyone can extend the AI support, and does it point humans to the handbook, agents to `AGENTS.md` and Docker users to `docker/agents/README.md`? |

## D. Classification and report

Classify each finding as an inconsistency, a justified difference (with its reason), a duplication or unclear, and report it in the table defined by the skill `mma-documentation-consistency`. The report is advisory; only the structural checks gate CI.
