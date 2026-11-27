---
name: mma-documentation-consistency
description: "Check that MMA's documentation, AI-support files and code agree: links, paths and anchors, entity structure and names, signposts, namelists, the input list, the deviations log and the meaning of statements. Use after editing README.md, a module README, AGENTS.md, anything under docs/ai or .agents/skills, the namelists or a core module, and whenever someone asks to verify, review or proofread MMA documentation."
metadata:
  mma-status: maintained
  mma-track: documentation
  mma-owner: maintainers
---

# Check documentation consistency

Two kinds of checks complement each other. Structural checks are deterministic and gate CI. Semantic checks need judgement; their findings are advisory and go to a human.

## Procedure

1. **Scope.** List the changed files, for example with `git diff --name-only` against the base branch, and the documents that describe them. The ownership map in `docs/ai/README.md` says which document owns which fact.
2. **Structural checks.** Run `python3 ai/check_ai_layer.py`, adding `--with-readme` when `README.md` changed. Report errors and warnings verbatim. Never fix generated signposts by hand; run `python3 ai/ai_signposts.py`.
3. **Semantic checks.** Work through section C of `references/verification-summary.md` for each changed statement. Compare the statement with its authoritative source, not with another summary.
4. **Classify** every finding as one of:
   - *inconsistency*: the statement contradicts its source; propose a fix;
   - *justified difference*: cite the reason; a departure from the CPC paper's model or approach needs an entry in `docs/deviations-from-paper.md`;
   - *duplication*: the statement restates a fact owned elsewhere; propose a link;
   - *unclear*: the source is ambiguous or missing; ask.
5. **Fix or propose.** Apply fixes only within the autonomy the user has given. Changes to authoritative documents need the user's agreement unless the user has authorised them.

## Report

Use this table, followed by the output of the structural checks:

| File and line | Statement | Authoritative source | Classification | Proposed fix |
| --- | --- | --- | --- | --- |

State which checks could not be performed, for example because a referenced file is missing.
