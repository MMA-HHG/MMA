#!/usr/bin/env bash
# Recipe: Codex CLI. Contract: docker/agents/README.md.
# OpenAI documents standalone release binaries at https://github.com/openai/codex#quickstart.
# No Node.js, npm, root privileges or shell start-up changes are needed.
# Runs without the shared background server; requires support for --no-daemon.
set -euo pipefail

request="${1:-}"
export CODEX_HOME="${CODEX_HOME:-${MMA_AGENTS_STATE:-${HOME}/.mma-agents}/codex}"
install_dir="${HOME}/.local/share/mma-agents/codex"
binary="${install_dir}/codex"
launcher="${HOME}/.local/bin/codex"
request_file="${install_dir}/request"

say() { echo "[mma-agents:codex] $*"; }
fail() { echo "[mma-agents:codex:error] $*" >&2; exit 1; }

case "${request}" in
    ""|latest) wanted=current; pinned=0 ;;
    *) [[ "${request}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] \
           || fail "Version must be X.Y.Z or 'latest', not '${request}'."
       wanted="${request}"; pinned=1 ;;
esac

case "$(uname -s)" in
    Linux) ;;
    *) fail "This container recipe supports Linux only." ;;
esac
case "$(uname -m)" in
    x86_64|amd64) target=x86_64-unknown-linux-musl ;;
    aarch64|arm64) target=aarch64-unknown-linux-musl ;;
    *) fail "Unsupported architecture: $(uname -m)." ;;
esac

installed_version() {
    [[ -x "$1" ]] || return 0
    "$1" --version 2>/dev/null | awk 'NR == 1 && $1 == "codex-cli" { print $2 }' || true
}

# Keep version policy beside the executable, outside the mounted credentials volume.
# The launcher leaves config.toml, credentials, model and permission settings untouched.
write_launcher() {
    {
        printf '#!/usr/bin/env bash\nset -euo pipefail\n'
        printf 'export CODEX_HOME="${CODEX_HOME:-${MMA_AGENTS_STATE:-${HOME}/.mma-agents}/codex}"\n'
        # Run the single executable without the shared background server.
        printf 'exec %q --no-daemon ' "${binary}"
        if (( pinned )); then
            printf '%s ' '-c' 'check_for_update_on_startup=false'
        fi
        printf '"$@"\n'
    } > "${launcher}.tmp"
    chmod 0755 "${launcher}.tmp"
    mv -f "${launcher}.tmp" "${launcher}"
}

mkdir -p "${CODEX_HOME}" "${install_dir}" "$(dirname "${launcher}")"
current="$(installed_version "${binary}")"
if [[ -n "${current}" && "$(cat "${request_file}" 2>/dev/null || true)" == "${wanted}" ]] \
        && { (( ! pinned )) || [[ "${current}" == "${request}" ]]; }; then
    write_launcher
    say "Codex CLI ${current} is already installed (${wanted})."
    exit 0
fi

command -v curl >/dev/null || fail "curl is missing from the image."
command -v tar >/dev/null || fail "tar is missing from the image."
asset="codex-${target}"
release_url=https://github.com/openai/codex/releases
if (( pinned )); then
    release_url+="/download/rust-v${request}"
else
    release_url+="/latest/download"
fi

# Download and validate before replacing a working installation.
stage="$(mktemp -d "${install_dir}/.install.XXXXXX")"
trap 'rm -rf "${stage}"' EXIT
say "Installing Codex CLI (${wanted}) from OpenAI's standalone release."
curl -fsSL --retry 3 --connect-timeout 20 --max-time 600 \
    "${release_url}/${asset}.tar.gz" -o "${stage}/codex.tar.gz" \
    || fail "Download failed; check network access and the requested release."
tar -xzf "${stage}/codex.tar.gz" -C "${stage}" -- "${asset}" \
    || fail "The release archive does not contain ${asset}."
current="$(installed_version "${stage}/${asset}")"
[[ "${current}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] \
    || fail "The downloaded executable reports no valid stable Codex CLI version."
(( ! pinned )) || [[ "${current}" == "${request}" ]] \
    || fail "Requested ${request}, but the downloaded executable reports ${current}."
chmod 0755 "${stage}/${asset}"
mv -f "${stage}/${asset}" "${binary}"
write_launcher
printf '%s\n' "${wanted}" > "${request_file}"
say "Installed Codex CLI ${current} (${wanted}); no login or permission settings were changed."
