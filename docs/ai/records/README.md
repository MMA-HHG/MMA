# Records

Records are dated observations, such as experience with installing and running MMA on a particular machine. Unlike authoritative documentation, anyone can write them, agents included, and they decay as machines change.

- One file per machine: `docs/ai/records/hpc/<machine>.md`, written from `docs/ai/templates/record-hpc.md`.
- A new record starts with record status `unverified`; a maintainer changes it to `verified` after reproducing the result.
- A record points to the machine's entry in `Modules/load_modules.sh` instead of copying module lists.
- Never record usernames, project or allocation IDs, internal hostnames or licence-server addresses.

The conventions are in the handbook, [docs/ai/README.md](../README.md).
