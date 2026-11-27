## Summary

<!-- What changes, and why. -->

## Scope

- [ ] Changes a core module (CUPRAD, CTDSE in `1DTDSE/`, Hankel) or the shared namelist
- [ ] Documentation only
- [ ] AI-support files (`AGENTS.md`, `docs/ai/`, `.agents/skills/`, signposts)

## Numerical impact

<!-- Core changes only. Use none, roundoff or intended, matching the Numerical-impact trailers of the commits (see CONTRIBUTING.md). -->

- Impact:
- Evidence (tests and examples run, machine, compiler, MPI, results):
- Not run, and why:
- Deviations log entry, if a modelling assumption changed:

## AI assistance

<!-- Required for every pull request. Keep one of: none, assisted, substantially generated. -->

- AI assistance:
- Agents and models (as in the `Assisted-by` trailers):
- Summary of the prompts or session:
- Parts affected by AI assistance:
- How those parts were checked:

## Checklist

- [ ] Commits carry the trailers required by `CONTRIBUTING.md`
- [ ] Reference data and tolerances are unchanged, or changed in a separate, justified commit
- [ ] `python3 ai/check_ai_layer.py` passes when documentation or AI-support files changed
- [ ] A human author takes responsibility for every change
