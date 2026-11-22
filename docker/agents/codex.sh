#!/usr/bin/env bash
# Recipe: Codex CLI. Contract: docker/agents/README.md.
# Install OpenAI's complete codex-package release, not the single-binary asset.
# No Node.js, npm, root privileges or shell start-up changes are needed.
# Preserve the package layout and run without the shared background server.
set -euo pipefail

request="${1:-}"
codex_state="${CODEX_HOME:-${MMA_AGENTS_STATE:-${HOME}/.mma-agents}/codex}"
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

# Helpers are resolved beside the real executable, not the public symlink.
# These checks do not sign in, call a model, or change user configuration.
package_works() {
    local executable directory helper
    executable="$(readlink -f "$1")" || return 1
    directory="$(dirname "${executable}")"
    [[ -x "${executable}" && -x "${directory}/codex-code-mode-host" ]] || return 1
    "${executable}" --no-daemon --help >/dev/null 2>&1 || return 1
    "${directory}/codex-code-mode-host" --help >/dev/null 2>&1 || return 1
    for helper in codex-exec-server codex-app-server; do
        if [[ -e "${directory}/${helper}" ]]; then
            [[ -x "${directory}/${helper}" ]] || return 1
            "${directory}/${helper}" --help >/dev/null 2>&1 || return 1
        fi
    done
}

# Keep version policy outside the credentials volume. Retain --no-daemon
# and forward user arguments without editing config.toml or credentials.
write_launcher() {
    {
        printf '#!/usr/bin/env bash\nset -euo pipefail\n'
        printf 'export CODEX_HOME="${CODEX_HOME:-${MMA_AGENTS_STATE:-${HOME}/.mma-agents}/codex}"\n'
        printf 'exec "$(readlink -f %q)" --no-daemon ' "${binary}"
        if (( pinned )); then
            printf '%s ' '-c' 'check_for_update_on_startup=false'
        fi
        printf '"$@"\n'
    } > "${launcher}.tmp"
    chmod 0755 "${launcher}.tmp"
    mv -f "${launcher}.tmp" "${launcher}"
}

mkdir -p "${codex_state}" "${install_dir}/packages" "$(dirname "${launcher}")"
current="$(installed_version "${binary}")"
if [[ "${current}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ \
        && -L "${binary}" \
        && "$(cat "${request_file}" 2>/dev/null || true)" == "${wanted}" ]] \
        && { (( ! pinned )) || [[ "${current}" == "${request}" ]]; } \
        && package_works "${binary}"; then
    write_launcher
    say "Codex CLI ${current} and its helpers are already installed (${wanted})."
    exit 0
fi

command -v curl >/dev/null || fail "curl is missing from the image."
command -v tar >/dev/null || fail "tar is missing from the image."
asset="codex-package-${target}.tar.gz"
release_url=https://github.com/openai/codex/releases

# Repair an incomplete unchanged request at its existing version.
# A repair must not silently become a CLI upgrade.
if (( pinned )); then
    release_version="${request}"
elif [[ "${current}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ \
        && "$(cat "${request_file}" 2>/dev/null || true)" == "${wanted}" ]]; then
    release_version="${current}"
else
    release_version=""
fi
if [[ -n "${release_version}" ]]; then
    release_url+="/download/rust-v${release_version}"
else
    release_url+="/latest/download"
fi

stage="$(mktemp -d "${install_dir}/.install.XXXXXX")"
release_dir=""
activated=0
cleanup() {
    rm -rf "${stage}"
    if (( ! activated )) && [[ -n "${release_dir}" ]]; then
        rm -rf "${release_dir}"
    fi
}
trap cleanup EXIT
mkdir "${stage}/package"
say "Installing the complete Codex CLI package (${release_version:-latest stable})."
curl -fsSL --retry 3 --connect-timeout 20 --max-time 600 \
    "${release_url}/${asset}" -o "${stage}/codex.tar.gz" \
    || fail "Package download failed; check network access and release availability."
tar -xzf "${stage}/codex.tar.gz" -C "${stage}/package" --no-same-owner \
    || fail "Cannot extract the Codex package."

# Accept a flat package, or an enclosing directory and/or bin/.
# Preserve the entire tree, including resources and executable permissions.
mapfile -t executables < <(find "${stage}/package" -type f -name codex -print)
(( ${#executables[@]} == 1 )) \
    || fail "Expected exactly one executable named codex in ${asset}."
entry="${executables[0]}"
current="$(installed_version "${entry}")"
[[ "${current}" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] \
    || fail "The package reports no valid stable Codex CLI version."
[[ -z "${release_version}" || "${current}" == "${release_version}" ]] \
    || fail "Requested ${release_version}, but the package reports ${current}."
package_works "${entry}" \
    || fail "Package check failed: codex must support --no-daemon and its code-mode host must run."
entry_relative="${entry#"${stage}/package/"}"

# Validate before switching. Keep old packages so running sessions retain helpers.
release_dir="${install_dir}/packages/${current}-${target}-${stage##*.}"
mv "${stage}/package" "${release_dir}"
ln -s "${release_dir}/${entry_relative}" "${stage}/codex"
mv -Tf "${stage}/codex" "${binary}"
activated=1
write_launcher
printf '%s\n' "${wanted}" > "${request_file}"
say "Installed Codex CLI ${current} with its complete package (${wanted})."
