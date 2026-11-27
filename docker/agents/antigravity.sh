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
    if [[ -L "${cli_state}" ]]; then
        [[ "$(readlink "${cli_state}")" == "${state_dir}" ]] && return 0
        rm -f "${cli_state}"
    elif [[ -d "${cli_state}" ]]; then
        if [[ -z "$(ls -A "${state_dir}")" ]]; then
            rmdir "${state_dir}"
            mv "${cli_state}" "${state_dir}"
            say "Moved the existing ${cli_state} into ${state_dir}."
        else
            [[ ! -e "${cli_state}.before-mma-agents" ]] \
                || fail "${cli_state} and ${state_dir} both hold state, and" \
                        "${cli_state}.before-mma-agents exists; resolve this by hand."
            mv "${cli_state}" "${cli_state}.before-mma-agents"
            say "Kept the state in ${state_dir}; the existing ${cli_state}" \
                "is now ${cli_state}.before-mma-agents."
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

# The installer offers only the current release and does nothing while agy exists, so
# an existing installation is set aside first and restored if anything fails.
install_current() {
    command -v curl >/dev/null \
        || fail "curl is missing from the image; it belongs in the Dockerfile's package list."
    mkdir -p "$(dirname "${command}")"
    local work backup="" version
    work="$(mktemp -d)"
    if [[ -e "${command}" || -L "${command}" ]]; then
        backup="$(mktemp "${command}.previous.XXXXXX")"
        mv -f "${command}" "${backup}"
    fi
    restore() {
        rm -f "${command}"
        if [[ -n "${backup}" ]]; then
            mv -f "${backup}" "${command}"
        fi
        rm -rf "${work}"
    }

    say "Installing the current Antigravity CLI release with Google's installer."
    if ! curl -fsSL --retry 3 --connect-timeout 20 "${installer_url}" -o "${work}/install.sh"; then
        restore
        fail "Downloading the installer failed; check network access."
    fi
    if ! bash "${work}/install.sh"; then
        restore
        fail "The installer failed; see its messages above."
    fi
    version="$(version_of "${command}")"
    if [[ -z "${version}" ]]; then
        restore
        fail "The installer finished, but ${command} reports no version."
    fi
    if (( pinned )) && [[ "${version}" != "${request}" ]]; then
        restore
        fail "Antigravity CLI ${request} is not available: the installer offers only the" \
             "current release, ${version}. Request 'antigravity' or 'antigravity@${version}'."
    fi
    if [[ -n "${backup}" ]]; then
        rm -f "${backup}"
    fi
    rm -rf "${work}"
    installed="${version}"
}

link_state
installed="$(version_of "$(executable)")"
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
