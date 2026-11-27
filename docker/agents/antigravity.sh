#!/usr/bin/env bash
# Recipe: Antigravity CLI, command agy. Contract: docker/agents/README.md.
# Installation method: Google's installer, https://antigravity.google/download/
#   curl -fsSL https://antigravity.google/cli/install.sh | bash
#
# Usage (through mma-agents): antigravity.sh [VERSION]
#   (none), latest  the current release, which then updates itself in the background
#   X.Y.Z           keep that version, with background updates off. The installer offers
#                   only the current release, so X.Y.Z must be installed or be current.
set -euo pipefail

request="${1:-}"
state_dir="${MMA_AGENTS_STATE:-${HOME}/.mma-agents}/antigravity"
cli_state="${HOME}/.gemini/antigravity-cli"
command="${HOME}/.local/bin/agy"
keep_dir="${HOME}/.local/share/mma-agents/antigravity"
kept="${keep_dir}/agy"
# An existing installation waits here, off PATH, while the installer runs.
previous="${HOME}/.local/share/mma-agents/antigravity.previous"
work=""
marker="# mma-agents:antigravity launcher"
installer_url="https://antigravity.google/cli/install.sh"

say() { echo "[mma-agents:antigravity] $*"; }
fail() { echo "[mma-agents:antigravity:error] $*" >&2; exit 1; }

case "${request}" in
    ""|latest) pinned=0 ;;
    *) [[ "${request}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] \
           || fail "Version must be X.Y.Z or 'latest', not '${request}'."
       pinned=1 ;;
esac

# The CLI keeps settings, logs and conversations in ~/.gemini/antigravity-cli and offers
# no setting to move them; a link puts them into the agents' state volume.
link_state() {
    mkdir -p "${state_dir}" "$(dirname "${cli_state}")"
    local aside
    if [[ -L "${cli_state}" ]]; then
        [[ "$(readlink "${cli_state}")" == "${state_dir}" ]] && return 0
        say "${cli_state} pointed to $(readlink "${cli_state}"); it now points to ${state_dir}."
        rm -f "${cli_state}"
    elif [[ -d "${cli_state}" ]]; then
        if [[ -z "$(ls -A "${state_dir}")" ]]; then
            rmdir "${state_dir}"
            mv "${cli_state}" "${state_dir}"
            say "Moved the existing ${cli_state} into ${state_dir}."
        else
            aside="${cli_state}.before-mma-agents.$(date +%Y%m%d-%H%M%S)"
            if [[ -e "${aside}" ]]; then
                aside="${aside}.$$"
            fi
            mv "${cli_state}" "${aside}"
            say "Kept the state in ${state_dir}; the existing ${cli_state} is now ${aside}."
        fi
    elif [[ -e "${cli_state}" ]]; then
        fail "${cli_state} is neither a directory nor a link; move it away and retry."
    fi
    ln -s "${state_dir}" "${cli_state}"
}

is_launcher() {
    [[ -f "$1" && "$(head -c 512 "$1" 2>/dev/null | tr -d '\0')" == *"${marker}"* ]]
}

# The executable: the installer's own file, or the kept one behind the launcher.
executable() {
    if is_launcher "${command}"; then
        echo "${kept}"
    else
        echo "${command}"
    fi
}

# Read the version without starting an update in the background.
version_of() {
    [[ -x "$1" ]] || return 0
    AGY_CLI_DISABLE_AUTO_UPDATE=true timeout 60 "$1" --version 2>/dev/null \
        | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' | head -n 1 || true
}

# Keep the installed version: the executable moves aside, and a launcher sets the CLI's
# documented AGY_CLI_DISABLE_AUTO_UPDATE, so that it does not replace itself.
keep_version() {
    mkdir -p "${keep_dir}"
    if ! is_launcher "${command}"; then
        mv -f "${command}" "${kept}"
    fi
    {
        printf '#!/usr/bin/env bash\n%s\n' "${marker}"
        printf '# Keeps Antigravity CLI %s; see docker/agents/README.md.\n' "$1"
        printf 'export AGY_CLI_DISABLE_AUTO_UPDATE=true\n'
        printf 'exec %q "$@"\n' "${kept}"
    } > "${command}.tmp"
    chmod 0755 "${command}.tmp"
    mv -f "${command}.tmp" "${command}"
}

# Return to the installer's layout, in which the CLI updates itself.
release_version() {
    if is_launcher "${command}"; then
        mv -f "${kept}" "${command}"
    fi
    rm -rf "${keep_dir}"
}

exists() { [[ -e "$1" || -L "$1" ]]; }

# Put back the installation that was set aside: after a failure or an interruption
# (exit trap), or at the next run if the previous one was killed outright.
restore_previous() {
    rm -f "${command}"
    if exists "${previous}"; then
        mv -f "${previous}" "${command}"
    fi
    if [[ -n "${work}" ]]; then
        rm -rf "${work}"
    fi
}

# True if version $1 is newer than version $2.
newer() {
    [[ "$1" != "$2" && "$(printf '%s\n%s\n' "$1" "$2" | sort -V | tail -n 1)" == "$1" ]]
}

# The installer offers only the current release and does nothing while agy exists, so
# an existing installation is set aside first and restored if anything fails.
install_current() {
    local version
    command -v curl >/dev/null \
        || fail "curl is missing from the image; it belongs in the Dockerfile's package list."
    mkdir -p "$(dirname "${command}")" "$(dirname "${previous}")"
    work="$(mktemp -d)"
    if exists "${command}"; then
        mv -f "${command}" "${previous}"
    fi
    trap restore_previous EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
    trap 'exit 129' HUP

    say "Installing the current Antigravity CLI release with Google's installer."
    curl -fsSL --retry 3 --connect-timeout 20 "${installer_url}" -o "${work}/install.sh" \
        || fail "Downloading the installer failed; check network access."
    # As with the documented "curl | bash", the installer gets no terminal input.
    bash "${work}/install.sh" < /dev/null \
        || fail "The installer failed; see its messages above."
    version="$(version_of "${command}")"
    [[ -n "${version}" ]] \
        || fail "The installer finished, but ${command} reports no version."
    if (( pinned )) && [[ "${version}" != "${request}" ]]; then
        fail "Antigravity CLI ${request} is not available: the installer offers only the" \
             "current release, ${version}. Request 'antigravity' or 'antigravity@${version}'."
    fi

    trap - EXIT INT TERM HUP
    rm -f "${previous}"
    rm -rf "${work}"
    work=""
    installed="${version}"
}

link_state
if exists "${previous}"; then
    restore_previous
    say "Restored the installation that an interrupted run had set aside."
fi
installed="$(version_of "$(executable)")"
# The CLI only updates forwards, so a version older than the installed one is no longer
# the current release and cannot be installed; this needs no download to tell.
if (( pinned )) && [[ -n "${installed}" ]] && newer "${installed}" "${request}"; then
    fail "Antigravity CLI ${request} is not available: ${installed} is installed, and the" \
         "installer offers only the current release. Request 'antigravity' or 'antigravity@${installed}'."
fi
action="is already installed"
if [[ -z "${installed}" ]] || { (( pinned )) && [[ "${installed}" != "${request}" ]]; }; then
    install_current
    action="was installed"
fi

if (( pinned )); then
    keep_version "${installed}"
else
    release_version
fi
current="$(version_of "${command}")"
[[ "${current}" == "${installed}" ]] \
    || fail "${command} reports '${current:-no version}' instead of ${installed}."

if (( pinned )); then
    say "Antigravity CLI ${installed} ${action} and is kept; background updates are off."
else
    say "Antigravity CLI ${installed} ${action}; it updates itself in the background."
fi
