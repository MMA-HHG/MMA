#!/usr/bin/env bash
set -euo pipefail

project_root="${MMA_PATH:-/MMA}"

die() {
    echo "[entrypoint:error] $*" >&2
    exit 1
}

set_num_proc_default() {
    local mpi_slots
    local n

    # Ask mpirun directly how many MPI slots it sees.
    mpi_slots=$(mpirun hostname | wc -l)

    export NUM_PROC_DOCKER="${mpi_slots}"

    if [[ -z "${NUM_PROC_DEFAULT_TDSE_1D:-}" ]]; then
        export NUM_PROC_DEFAULT_TDSE_1D="${mpi_slots}"
    fi

    if [[ -z "${NUM_PROC_DEFAULT_CUPRAD:-}" ]]; then
        n=1
        while (( 2 * n <= mpi_slots )); do
            n=$((2 * n))
        done
        export NUM_PROC_DEFAULT_CUPRAD="${n}"
    fi

    echo "[entrypoint] MPI slots detected by mpirun: ${mpi_slots}"
    echo "[entrypoint] NUM_PROC_DEFAULT_CUPRAD=${NUM_PROC_DEFAULT_CUPRAD}"
    echo "[entrypoint] NUM_PROC_DEFAULT_TDSE_1D=${NUM_PROC_DEFAULT_TDSE_1D}"

    if (( NUM_PROC_DEFAULT_CUPRAD < 2 )); then
        echo "[entrypoint:warning] CUPRAD requires at least 2 MPI processes." >&2
    fi
}

check_project_root() {
    [[ -d "${project_root}" ]] \
        || die "Project directory does not exist: ${project_root}"

    [[ -f "${project_root}/CMakeLists.txt" ]] \
        || die "${project_root} does not look like the MMA project root. Did you forget -v .:/MMA ?"
}

check_project_writable() {
    local test_file

    test_file="$(mktemp "${project_root}/.mma-write-test.XXXXXX")" || {
        cat >&2 <<ERROR
[entrypoint:error] The mounted MMA project directory is not writable.

Directory: ${project_root}
User: $(id -un)
UID: $(id -u)
GID: $(id -g)

This usually means that the container user does not match the owner of the
mounted repository on the host system.

On Linux, rebuild the image with matching UID/GID:

  docker build \\
    --build-arg USER_UID=\$(id -u) \\
    --build-arg USER_GID=\$(id -g) \\
    -t mma .

See the Docker documentation for build arguments and --build-arg:
https://docs.docker.com/build/building/variables/#build-arguments
ERROR
        exit 1
    }

    rm -f "${test_file}"
}

expected_executables_present() {
    [[ -x "${CUPRAD_BUILD}/make_start.e" ]] \
        && [[ -x "${CUPRAD_BUILD}/cuprad.e" ]] \
        && [[ -x "${TDSE_1D_BUILD}/TDSE.e" ]]
}

run_auto_build_if_needed() {
    if [[ "${MMA_AUTO_BUILD:-1}" == "0" ]]; then
        echo "[entrypoint] Automatic MMA build disabled by MMA_AUTO_BUILD=0."
        return
    fi

    check_project_root
    check_project_writable

    if [[ "${MMA_FORCE_REBUILD:-0}" != "1" ]] && expected_executables_present; then
        echo "[entrypoint] MMA appears to be already compiled; skipping automatic build."
        return
    fi

    if [[ "${MMA_FORCE_REBUILD:-0}" == "1" ]]; then
        echo "[entrypoint] MMA_FORCE_REBUILD=1; running build."
    else
        echo "[entrypoint] MMA executables were not found; running automatic build."
    fi

    mma-build

    expected_executables_present \
        || die "Build finished, but one or more expected executables are missing."
}

set_num_proc_default
run_auto_build_if_needed

# Default command if the user did not provide one.
if [[ "$#" -eq 0 ]]; then
    set -- bash
fi

exec "$@"

# Ensure the working directory used by the MMA basics tutorial exists.
mkdir -p "${MULTISCALE_WORK_DIR:-/MMA/work_dir}/mma_basics"
