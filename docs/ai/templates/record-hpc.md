---
machine: "<public name of the cluster, as used by its computing centre>"
hpc-entry: "<value of HPC that selects the module set in Modules/load_modules.sh, or none>"
record-status: unverified
last-verified: "<YYYY-MM-DD>"
mma-commit: "<commit hash of the MMA version used>"
verified-by: "<agent and model, or the maintainer's handle>"
---

# Record: <machine>

<!-- Copy to docs/ai/records/hpc/<machine>.md.
     Never record usernames, project or allocation IDs, internal hostnames or licence-server addresses. -->

## Environment

<Compilers, MPI, HDF5 and FFTW versions observed. Reference the module set in Modules/load_modules.sh instead of copying it.>

## Build

<What worked, which README steps were followed, and any deviation from them.>

## Running

<Scheduler, partitions or queues suitable for tests, and limits that matter for MMA: wall-clock time, memory, MPI ranks.>

## Known issues

- <Date: symptom, cause, fix.>

## Not verified

- <What was not tested.>
