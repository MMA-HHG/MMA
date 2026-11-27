# Deviations from the published version

The CPC paper ([doi:10.1016/j.cpc.2026.110417](https://doi.org/10.1016/j.cpc.2026.110417), open access) describes MMA version v1.1.0, whose program files are archived in the CPC Library ([doi:10.17632/g422fdhx7z.1](https://doi.org/10.17632/g422fdhx7z.1)). It introduces MMA's core philosophy: the multiscale approach, the physical model and its modelling assumptions. The decisive analyses, such as numerical details, validation and implementation decisions, are documented in the lower layers: the module READMEs and the source.

This log is the authoritative record of later departures from what the paper describes. Such a departure is justified only when it is recorded here; changes below that level are documented in the lower layers as usual.

## Entries

| Date | Change | Paper section | Effect on results | Pull request |
| --- | --- | --- | --- | --- |
| — | No deviations recorded yet. | — | — | — |

## Adding an entry

Add one row per change, in the same pull request as the change itself:

- **Change**: what differs from the paper, in one or two sentences.
- **Paper section**: the section, equation or table concerned.
- **Effect on results**: none, roundoff or intended, matching the commit's numerical-impact trailer, with a short note for intended changes.
- **Pull request**: a link to the pull request that holds the evidence.
