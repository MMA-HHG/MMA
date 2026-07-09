#!/usr/bin/env bash
set -euo pipefail

project_root="${MSM_PATH:-/MMA}"

is_disabled() {
    case "${1:-}" in
        0|false|False|FALSE|no|No|NO) return 0 ;;
        *) return 1 ;;
    esac
}

die() {
    echo "[entrypoint:error] $*" >&2
    exit 1
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
    if is_disabled "${MMA_AUTO_BUILD:-1}"; then
        echo "[entrypoint] Automatic MMA build disabled by MMA_AUTO_BUILD=${MMA_AUTO_BUILD}."
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

run_auto_build_if_needed

# Default command if the user did not provide one.
if [[ "$#" -eq 0 ]]; then
    set -- bash
fi

exec "$@"
